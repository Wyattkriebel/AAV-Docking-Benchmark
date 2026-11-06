"""Strict read-only archive/provenance replay; no native or accuracy execution."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tarfile

import gemmi
import pytest

import aav_benchmark.heldout_serialization as serialization


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("serialization_verifier", ROOT / "scripts/verify_serialization_archives.py")
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def archive_bytes(entries):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as archive:
        for name, content, kind in entries:
            item = tarfile.TarInfo(name)
            item.type = kind
            item.size = len(content)
            if kind in (tarfile.SYMTYPE, tarfile.LNKTYPE):
                item.linkname = "target"
            archive.addfile(item, io.BytesIO(content))
    raw = stream.getvalue()
    manifest = {"archive_sha256": verifier.digest(raw), "files": [
        {"path": name, "bytes": len(content), "sha256": verifier.digest(content)} for name, content, _ in entries]}
    return raw, manifest


def test_valid_archive_is_read_in_memory_without_extracting(monkeypatch):
    raw, manifest = archive_bytes([("case/input.cif", b"source", tarfile.REGTYPE)])
    def forbidden(*args, **kwargs):
        pytest.fail("Archive verification must never extract onto the filesystem")
    monkeypatch.setattr(tarfile.TarFile, "extract", forbidden)
    monkeypatch.setattr(tarfile.TarFile, "extractall", forbidden)
    assert verifier.read_archive(raw, manifest, expected_count=1) == {"case/input.cif": b"source"}


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a/../escape", "./case", "a//b", "a\\b", "."])
def test_unsafe_or_noncanonical_member_paths_are_rejected(name):
    raw, manifest = archive_bytes([(name, b"source", tarfile.REGTYPE)])
    with pytest.raises(ValueError, match="member"):
        verifier.read_archive(raw, manifest, expected_count=1)


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.DIRTYPE, tarfile.FIFOTYPE])
def test_nonregular_archive_members_are_rejected(kind):
    raw, manifest = archive_bytes([("case", b"", kind)])
    with pytest.raises(ValueError, match="regular"):
        verifier.read_archive(raw, manifest, expected_count=1)


def test_duplicate_members_are_rejected_even_with_identical_bytes():
    raw, manifest = archive_bytes([("case", b"source", tarfile.REGTYPE)] * 2)
    manifest["files"] = manifest["files"][:1]
    with pytest.raises(ValueError, match="Duplicate"):
        verifier.read_archive(raw, manifest, expected_count=1)
    manifest["files"] *= 2
    with pytest.raises(ValueError, match="duplicate archive manifest"):
        verifier.read_archive(raw, manifest, expected_count=1)


@pytest.mark.parametrize("change", ["archive_hash", "member_hash", "member_size", "missing", "extra", "denominator"])
def test_archive_checksum_exact_member_set_and_sizes_are_enforced(change):
    raw, manifest = archive_bytes([("case", b"source", tarfile.REGTYPE)])
    count = 1
    if change == "archive_hash":
        manifest["archive_sha256"] = "0" * 64
    elif change == "member_hash":
        manifest["files"][0]["sha256"] = "0" * 64
    elif change == "member_size":
        manifest["files"][0]["bytes"] += 1
    elif change == "missing":
        manifest["files"].append({"path": "missing", "bytes": 0, "sha256": verifier.digest(b"")}); count = 2
    elif change == "extra":
        manifest["files"][0]["path"] = "other"
    else:
        count = 2
    with pytest.raises(ValueError):
        verifier.read_archive(raw, manifest, expected_count=count)


@pytest.fixture(scope="module")
def fixed_pair():
    archives, receipts = [], []
    for name, count in (("structural_test", 99), ("structural_test_v2", 107)):
        parent = ROOT / "results" / name
        manifest = json.loads((parent / "source_archive_manifest.json").read_bytes())
        archives.append(verifier.read_archive((parent / "source_preparation.tar.gz").read_bytes(), manifest, expected_count=count))
        receipts.append(json.loads((parent / "source_preparation.json").read_bytes()))
    amendment = json.loads((ROOT / "results/structural_test_v2/serialization_amendment.json").read_bytes())
    return (*archives, *receipts, amendment)


def reseal_test_provenance(original, amended, receipt, amendment):
    """Make a mutation internally hash-consistent to test semantic rejection."""
    for case in receipt["cases"]:
        prefix = case["case_id"] + "/"
        case["artifact_sha256"] = {p.removeprefix(prefix): verifier.digest(raw) for p, raw in amended.items()
            if p.startswith(prefix) and p != prefix + "source_preparation.json"}
        amended[prefix + "source_preparation.json"] = json_bytes(case)
    amended["source_preparation.json"] = json_bytes(receipt)
    amendment["source_preparation_v2_sha256"] = verifier.digest(amended["source_preparation.json"])
    amendment["changed_original_files"] = [{"path": p, "prior_sha256": verifier.digest(raw), "sha256": verifier.digest(amended[p])}
        for p, raw in original.items() if raw != amended[p]]
    amendment["unchanged_original_files"] = [{"path": p, "sha256": verifier.digest(raw)}
        for p, raw in original.items() if raw == amended[p]]
    amended["serialization_amendment.json"] = json_bytes(amendment)


@pytest.mark.parametrize("mutation, message", [
    ("coordinate", "protected atom field Cartn_x"),
    ("occupancy", "protected atom field occupancy"),
    ("formatting", "reproduce archived V2 CIF bytes"),
    ("graph", "predeclared metadata files"),
    ("source_decision", "source decision"),
    ("protein_record", "Protein preparation receipts disagree"),
    ("protein_decision", "protein source decision"),
    ("aggregate", "Aggregate source decisions"),
])
def test_internally_resealed_tampering_cannot_change_physical_atoms_graphs_or_decisions(fixed_pair, mutation, message):
    original, amended, original_receipt, amended_receipt, amendment = copy.deepcopy(fixed_pair)
    if mutation in ("coordinate", "occupancy"):
        document = gemmi.cif.read_string(amended["2RDK/receptor.cif"].decode())
        atoms = document.sole_block().get_mmcif_category("_atom_site.")
        field = "Cartn_x" if mutation == "coordinate" else "occupancy"
        atoms[field][0] = str(float(atoms[field][0]) + .25)
        document.sole_block().set_mmcif_category("_atom_site.", atoms)
        amended["2RDK/receptor.cif"] = (document.as_string().rstrip() + "\n#\n").encode()
    elif mutation == "formatting":
        amended["2RDK/receptor.cif"] += b"\n"
    elif mutation == "graph":
        amended["2RDK/chemistry.json"] += b"\n"
    elif mutation == "source_decision":
        amended_receipt["cases"][0]["qualification"] = "changed"
    elif mutation in ("protein_record", "protein_decision"):
        amended_receipt["cases"][0]["protein"]["distance_based_selection"] = True
        if mutation == "protein_decision":
            for name in ("protein_preparation.json", "receptor.cif.manifest.json"):
                amended["2RDK/" + name] = json_bytes(amended_receipt["cases"][0]["protein"])
    else:
        amended_receipt["model_family_heldout_claim"] = True
    if mutation in ("coordinate", "occupancy", "formatting"):
        protein = amended_receipt["cases"][0]["protein"]
        protein["receptor_sha256"] = verifier.digest(amended["2RDK/receptor.cif"])
        for name in ("protein_preparation.json", "receptor.cif.manifest.json"):
            amended["2RDK/" + name] = json_bytes(protein)
    reseal_test_provenance(original, amended, amended_receipt, amendment)
    with pytest.raises(ValueError, match=message):
        verifier.verify_metadata_pair(original, amended, original_receipt, amended_receipt, amendment,
            module_sha256=verifier.PINNED_FILES["aav_benchmark/heldout_serialization.py"])


def test_amendment_must_enumerate_every_changed_and_unchanged_original(fixed_pair):
    original, amended, original_receipt, amended_receipt, amendment = copy.deepcopy(fixed_pair)
    amendment["unchanged_original_files"].pop()
    with pytest.raises(ValueError, match="enumerate the actual"):
        verifier.verify_metadata_pair(original, amended, original_receipt, amended_receipt, amendment,
            module_sha256=verifier.PINNED_FILES["aav_benchmark/heldout_serialization.py"])


def test_complete_public_bundle_replays_read_only_without_tools_processed_or_network(tmp_path, monkeypatch):
    paths = [*verifier.PINNED_FILES, "aav_benchmark/__init__.py", "scripts/verify_serialization_archives.py"]
    for name in paths:
        destination = tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, destination)
    before = {name: verifier.digest((tmp_path / name).read_bytes()) for name in paths}
    def forbidden(*args, **kwargs):
        pytest.fail("Read-only verification must not call network or build_source_v2")
    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setattr(serialization, "build_source_v2", forbidden)
    result = verifier.verify(tmp_path)
    assert result["original_members"] == 99 and result["v2_members"] == 107
    assert result["changed_metadata_members"] == 25 and result["unchanged_original_members"] == 74
    assert result["fresh_CIF_byte_matches"] == result["cases"] == 6
    assert result["unchanged_PDB_PDBQT_SDF_members"] == 34 and result["unchanged_source_graph_members"] == 12
    assert result["atom_counts"] == dict(zip(serialization.CASE_IDS, (760, 1983, 3076, 4769, 1042, 6221)))
    assert all(type(n) is int for n in result["atom_counts"].values())
    assert not (tmp_path / "tools").exists() and not (tmp_path / "data/processed").exists()
    assert before == {name: verifier.digest((tmp_path / name).read_bytes()) for name in paths}
    completed = subprocess.run([sys.executable, "scripts/verify_serialization_archives.py"], cwd=tmp_path,
                               text=True, capture_output=True, check=True)
    assert json.loads(completed.stdout) == result and not completed.stderr
    altered = tmp_path / "results/structural_test_v2/source_preparation.json"
    altered.write_bytes(altered.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Pinned source SHA256"):
        verifier.verify(tmp_path)
