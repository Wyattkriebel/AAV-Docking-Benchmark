"""Verified temporary source fixtures; no report, native or prediction calls."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

from .replay_workspace import hydrate_verified_archives, relative_path

ROOT = Path(__file__).resolve().parents[1]
DISTRIBUTION = Path("benchmarks/method_protocol_public_distribution.json")
DISTRIBUTION_SHA256 = "bed3115a6ab2685aa3c9cd790ff5080e931ee813c37c720188b4e44901570beb"
RECEIPT = ".method-source-fixture.json"
GROUPS = {"prepared", "haddock", "development"}
AMENDMENT = Path("results/method_protocol_source_distribution/publication_amendment.json")
AMENDMENT_SHA256 = "fe1d02b27fb05400b3e4fd214642114129c37611afee88c1c51c8cdfc61c1e6c"


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _checked(root, record):
    relative_path(record["path"])
    root = Path(root)
    path = root / record["path"]
    parts = Path(record["path"]).parts
    if any(root.joinpath(*parts[:i]).is_symlink() for i in range(1, len(parts) + 1)) or not path.is_file():
        raise ValueError("Missing regular public source: " + record["path"])
    if path.stat().st_size != record["bytes"] or _digest(path) != record["sha256"]:
        raise ValueError("Public source pin changed: " + record["path"])
    return path


def _distribution(root):
    root = Path(root)
    path = root / DISTRIBUTION
    if (any(root.joinpath(*DISTRIBUTION.parts[:i]).is_symlink()
            for i in range(1, len(DISTRIBUTION.parts) + 1))
            or _digest(path) != DISTRIBUTION_SHA256):
        raise ValueError("Method source distribution changed")
    value = json.loads(path.read_text())
    if (value.get("schema_version") != 1
            or value.get("frozen_implementations_unchanged") is not True
            or value.get("no_new_predictions_native_scoring_or_accuracy_calls") is not True):
        raise ValueError("Unexpected method source distribution")
    names = [record["path"] for record in value["root_files"]]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate public source identity")
    return value


def verify_source_distribution(source_root=ROOT):
    """Check all declared files and the unchanged 14/59 code bindings."""
    root = Path(source_root).resolve()
    value = _distribution(root)
    for record in value["root_files"]:
        _checked(root, record)
    if len(value["frozen_implementation_sources"]) != 14:
        raise ValueError("All fourteen frozen implementations are required")
    for record in value["frozen_implementation_sources"]:
        _checked(root, record)
    entries = 0
    for record in value["current_frozen_source_maps"]:
        plan = json.loads(_checked(root, record["plan"]).read_text())
        if plan[record["field"]] != record["sources"]:
            raise ValueError("Frozen method source map changed")
        for name, digest in record["sources"].items():
            relative_path(name)
            if _digest(root / name) != digest:
                raise ValueError("Frozen method implementation changed: " + name)
            entries += 1
    if entries != 59:
        raise ValueError("All fifty-nine source-map entries are required")
    return value


def verify_publication_amendment(source_root=ROOT):
    """Keep current public tests distinct from exact historical archive copies."""
    root = Path(source_root).resolve()
    if _digest(root / AMENDMENT) != AMENDMENT_SHA256:
        raise ValueError("Public test adaptation amendment changed")
    value = json.loads((root / AMENDMENT).read_text())
    distribution = _distribution(root)
    if (value["frozen_implementation_sources"] != distribution["frozen_implementation_sources"]
            or value.get("current_source_map_entries_unchanged") != 59
            or value.get("original_scientific_freezes_and_receipts_unchanged") is not True):
        raise ValueError("Public test adaptation changed the frozen method contract")
    for record in value["test_adaptations"]:
        _checked(root, record["active_public_test"])
        with tarfile.open(_checked(root, record["frozen_historical_archive"]), "r:gz") as archive:
            members = [m for m in archive.getmembers() if m.name == record["historical_member"]]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError("Historical test has no unique ordinary archive member")
            raw = archive.extractfile(members[0]).read()
        if (len(raw) != record["historical_test_bytes"]
                or hashlib.sha256(raw).hexdigest() != record["historical_test_sha256"]):
            raise ValueError("Historical public-test provenance changed")
    if len(value["test_adaptations"]) != 4:
        raise ValueError("All four public test adaptations must retain their historical pins")
    return {"adapted_public_tests_verified": 4, "historical_archived_test_pins_verified": 4}


def _private_root(path):
    path = Path(path).resolve()
    bases = {Path(tempfile.gettempdir()).resolve(), Path("/tmp").resolve()}
    if not any(path.is_relative_to(base) and path != base for base in bases):
        raise ValueError("Source fixtures require a private temporary directory")
    if any((parent / ".git").exists()
           or ((parent / "HEAD").is_file() and (parent / "config").is_file()
               and (parent / "objects").is_dir()) for parent in (path, *path.parents)):
        raise ValueError("Source fixtures cannot be inside a Git checkout")
    return path


def prepare_source_fixture(destination, *, source_root=ROOT, groups=("prepared",)):
    """Copy declared bytes and hydrate selected archives into empty TEMP only."""
    if Path(destination).is_symlink():
        raise ValueError("Source fixture cannot be a symlink")
    destination, source_root = _private_root(destination), Path(source_root).resolve()
    if destination.is_relative_to(source_root):
        raise ValueError("Source fixture cannot modify the source repository")
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise ValueError("Use an empty private temporary source fixture")
    selected = set(groups)
    if not selected or not selected <= GROUPS:
        raise ValueError("Unknown or empty source fixture group")
    if "haddock" in selected:
        selected.add("prepared")
    value = verify_source_distribution(source_root)
    destination.mkdir(parents=True, exist_ok=True)
    for record in value["root_files"]:
        source = _checked(source_root, record)
        target = destination / record["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        with source.open("rb") as stream, target.open("xb") as output:
            shutil.copyfileobj(stream, output)
        _checked(destination, record)
    with (source_root / DISTRIBUTION).open("rb") as stream, (destination / DISTRIBUTION).open("xb") as output:
        shutil.copyfileobj(stream, output)
    entries = [{key: entry[key] for key in ("archive", "manifest", "target_prefix")}
               for entry in value["archives"] if entry["group"] in selected]
    hydration = hydrate_verified_archives(destination, entries)
    receipt = {"schema_version": 1, "distribution_sha256": DISTRIBUTION_SHA256,
               "groups": sorted(selected), "copied_public_files": len(value["root_files"]),
               "hydration": hydration, "accuracy_computed": False,
               "native_executables_used": False, "new_predictions_or_scoring": False}
    with (destination / RECEIPT).open("x") as output:
        json.dump(receipt, output, indent=2, sort_keys=True, allow_nan=False)
        output.write("\n")
    return receipt


def restore_preflight_receptor_aliases(root):
    """Restore exactly twelve historical receptor aliases inside TEMP only.

    Bytes come from public source archives, never the ignored processed tree.
    The generic hydrator continues to prohibit processed-tree destinations.
    """
    root = _private_root(root)
    receipt = json.loads((root / RECEIPT).read_text())
    if (receipt.get("distribution_sha256") != DISTRIBUTION_SHA256
            or "haddock" not in receipt.get("groups", [])):
        raise ValueError("Aliases require a verified private HADDOCK fixture")
    aliases = _distribution(root)["historical_preflight_receptor_aliases"]
    expected = {f"data/processed/heldout_pose_v3/{profile}/{case}/receptor.pdb"
                for profile in ("strict", "repair")
                for case in ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")}
    if len(aliases) != 12 or {record["destination"] for record in aliases} != expected:
        raise ValueError("Only twelve exact historical receptor aliases are permitted")
    for record in aliases:
        manifest = json.loads(_checked(root, record["source_manifest"]).read_text())
        archive_path = _checked(root, record["source_archive"])
        rows = [row for row in manifest["files"] if row["path"] == record["member"]]
        if (len(rows) != 1 or rows[0]["sha256"] != record["sha256"]
                or rows[0]["bytes"] != record["bytes"]):
            raise ValueError("Historical receptor alias differs from public source")
        with tarfile.open(archive_path, "r:gz") as archive:
            members = [m for m in archive.getmembers() if m.name == record["member"]]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError("Historical receptor alias has no unique ordinary source")
            raw = archive.extractfile(members[0]).read()
        if len(raw) != record["bytes"] or hashlib.sha256(raw).hexdigest() != record["sha256"]:
            raise ValueError("Historical receptor archive bytes changed")
        target = root / record["destination"]
        if any(parent.is_symlink() for parent in (target, *target.parents)):
            raise ValueError("Historical receptor alias traverses a symlink")
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as output:
            output.write(raw)
    return {"restored_historical_receptor_aliases": 12, "copied_from_ignored_live_tree": False,
            "native_executables_used": False, "accuracy_computed": False}
