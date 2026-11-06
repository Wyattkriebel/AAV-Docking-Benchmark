"""Reparse archived site predictions without installing or running predictors.

This verifies execution evidence and normalized candidates. It does not use
native ligands, repeat site scoring, or infer biology from workflow failures.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import signal
import stat
import tempfile

from . import p2rank, sites
from .provenance import sha256, write_json
from .site_archive import restore_run


CENTER_TOLERANCE_A = 1e-10
INPUT_REVISIONS = frozenset({"original", "fpocket_compatibility"})


def _filename(value):
    if (not isinstance(value, str) or not value or value in {".", ".."}
            or any(char in value for char in ("/", "\\", ":", "\x00"))):
        raise ValueError("site input archive requires portable filenames without traversal")
    return value


def _digest(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("site input archive requires lowercase SHA256 digests")
    return value


def _regular(path):
    if not stat.S_ISREG(Path(path).lstat().st_mode):
        raise ValueError("site input archive must be a regular file, not a symlink")
    return Path(path)


def restore_site_inputs(root, condition_id, revision, destination):
    """Restore one frozen input revision to a fresh directory and verify bytes.

    Only files declared in ``results/site_discovery/<condition>/inputs`` are
    read. Historical absolute source/map paths inside manifests are metadata;
    they are never opened. Returns a portable receipt with verified file hashes.
    """
    root, destination = Path(root), Path(destination)
    condition_id = _filename(condition_id)
    if revision not in INPUT_REVISIONS:
        raise ValueError("unknown site input revision")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("site input restoration requires a fresh destination")
    directory = root / "results/site_discovery" / condition_id / "inputs"
    manifest_path = _regular(directory / "manifest.json")
    manifest = json.loads(manifest_path.read_text())
    if manifest["condition_id"] != condition_id or not isinstance(manifest["files"], list):
        raise ValueError("site input archive condition or file inventory differs")
    selected, seen = [], set()
    for entry in manifest["files"]:
        if entry["revision"] not in INPUT_REVISIONS:
            raise ValueError("input manifest contains an unknown revision")
        filename, archive_name = _filename(entry["filename"]), _filename(entry["archive"])
        identity = (entry["revision"], filename)
        if identity in seen:
            raise ValueError("duplicate file in site input archive revision")
        seen.add(identity)
        _digest(entry["sha256"])
        _digest(entry["archive_sha256"])
        if entry["revision"] == revision:
            selected.append(dict(entry))
    if not selected or "receptor.cif" not in {entry["filename"] for entry in selected}:
        raise ValueError("site input revision has no receptor.cif")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".site-input-restore.", dir=destination.parent))
    try:
        for entry in selected:
            compressed = _regular(directory / entry["archive"])
            if sha256(compressed) != entry["archive_sha256"]:
                raise ValueError("site input compressed archive checksum mismatch")
            digest = hashlib.sha256()
            with gzip.open(compressed, "rb") as source, (temporary / entry["filename"]).open("xb") as output:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    output.write(chunk)
                    digest.update(chunk)
            if digest.hexdigest() != entry["sha256"]:
                raise ValueError("restored site input checksum mismatch")
        if destination.exists() or destination.is_symlink():
            raise FileExistsError("site input destination appeared during restoration")
        temporary.rename(destination)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return {"condition_id": condition_id, "revision": revision,
            "source_manifest_sha256": sha256(manifest_path), "files": selected}


def _same_predictions(expected, observed, path="predictions"):
    """Only centers get a numerical tolerance; all other values match exactly."""
    if type(expected) is not type(observed):
        raise ValueError(f"raw prediction type differs at {path}")
    if isinstance(expected, dict):
        if set(expected) != set(observed):
            raise ValueError(f"raw prediction fields differ at {path}")
        for key in expected:
            if key in {"center", "center_A"}:
                a, b = expected[key], observed[key]
                if (not isinstance(a, list) or not isinstance(b, list) or len(a) != 3 or len(b) != 3
                        or any(isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v) for v in a + b)
                        or any(abs(x - y) > CENTER_TOLERANCE_A for x, y in zip(a, b))):
                    raise ValueError(f"raw prediction center differs at {path}.{key}")
            else:
                _same_predictions(expected[key], observed[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(observed):
            raise ValueError(f"raw prediction count differs at {path}")
        for index, (left, right) in enumerate(zip(expected, observed)):
            _same_predictions(left, right, f"{path}[{index}]")
    elif expected != observed:
        raise ValueError(f"raw prediction value differs at {path}")


def _validate_command(record, method, lock, lock_sha):
    command = record.get("command")
    if not isinstance(command, list) or not command or any(not isinstance(value, str) for value in command):
        raise ValueError("site command is missing or malformed")
    if method == "fpocket":
        config = lock["methods"]["fpocket"]
        if Path(command[0]).name != "fpocket" or command[1:] != ["-f", "receptor.cif", *config["parameters"]]:
            raise ValueError("fpocket command differs from the frozen parameters")
        for key in ("release", "commit", "source_sha256"):
            if record.get(key) != config[key]:
                raise ValueError("fpocket release/source provenance differs from the lock")
        if record.get("method_lock_sha256") != lock_sha:
            raise ValueError("fpocket method-lock digest differs")
        installation = record.get("installation", {})
        if (installation.get("executable_sha256") != record.get("executable_sha256")
                or any(installation.get(key) != config[key] for key in ("release", "commit", "source_sha256"))):
            raise ValueError("fpocket recorded installation provenance differs")
    else:
        config = lock["parameters"]
        expected = ["-threads", str(config["threads"]), "-visualizations", "1" if config["visualizations"] else "0",
                    "-seed", str(config["seed"]), "-model", config["model"]]
        if (len(command) != 6 + len(expected) or Path(command[0]).name != "prank"
                or command[1:3] != [config["command"], "-f"] or command[4] != "-o"
                or Path(command[3]).name != "receptor.cif" or Path(command[5]).name != "raw"
                or command[6:] != expected or command[3] != record.get("input_path")):
            raise ValueError("P2Rank command differs from the frozen model/thread/seed parameters")
        if (record.get("parameters") != config or record.get("lock_sha256") != lock_sha
                or record.get("version") != lock["version"]
                or record.get("executable_sha256") != lock["resource_sha256"]["prank"]
                or record.get("timeout_seconds") != config["timeout_seconds"]):
            raise ValueError("P2Rank recorded parameters/version/provenance differ from the lock")


def _input_evidence(case, original, corrected):
    old_manifest = json.loads((original / "receptor.cif.manifest.json").read_text())
    new_manifest = json.loads((corrected / "receptor.cif.manifest.json").read_text())
    old_sha, new_sha = sha256(original / "receptor.cif"), sha256(corrected / "receptor.cif")
    for directory, manifest, digest in ((original, old_manifest, old_sha), (corrected, new_manifest, new_sha)):
        if (manifest.get("source_sha256") != case["source_sha256"] or manifest.get("receptor_sha256") != digest
                or manifest.get("native_ligand_access") is not False or manifest.get("distance_based_selection") is not False
                or manifest.get("atom_map_sha256") != sha256(directory / "receptor.cif.atoms.tsv.gz")):
            raise ValueError("archived receptor preparation/source/atom-map evidence differs")
    old_bytes, new_bytes = (original / "receptor.cif").read_bytes(), (corrected / "receptor.cif").read_bytes()
    if new_bytes != old_bytes.rstrip() + b"\n#\n":
        raise ValueError("fpocket compatibility revision changed more than the comment delimiter")
    revision = json.loads((corrected / "receptor.cif.revision.json").read_text())
    if (revision.get("prior_receptor_sha256") != old_sha or revision.get("receptor_sha256") != new_sha
            or revision.get("prior_preparation_manifest_sha256") != sha256(original / "receptor.cif.manifest.json")
            or revision.get("coordinates_exactly_unchanged") is not True
            or revision.get("author_atom_identities_unchanged") is not True
            or revision.get("method_parameters_changed") is not False
            or revision.get("atom_map_sha256") != old_manifest["atom_map_sha256"]
            or new_manifest["atom_map_sha256"] != old_manifest["atom_map_sha256"]
            or old_manifest["protein_heavy_atoms"] != new_manifest["protein_heavy_atoms"]):
        raise ValueError("fpocket compatibility revision receipt differs from the archived inputs")
    return {"original_sha256": old_sha, "corrected_sha256": new_sha,
            "source_sha256": case["source_sha256"], "protein_heavy_atoms": old_manifest["protein_heavy_atoms"]}


def _verify_attempt(run_dir, restored, input_dir, *, method, lock, lock_sha, initial=False, original_sha=None,
                    installation=None):
    record_name = "execution.json" if method == "fpocket" else "result.json"
    archive = restore_run(run_dir, restored)
    if sha256(_regular(run_dir / record_name)) != archive["record_sha256"]:
        raise ValueError("live site record differs from the verbatim archived record")
    record = json.loads((restored / record_name).read_text())
    if record.get("method") != method:
        raise ValueError("archived site method identity differs")
    _validate_command(record, method, lock, lock_sha)
    receptor = input_dir / "receptor.cif"
    if record.get("input_sha256") != sha256(receptor):
        raise ValueError("site attempt used a different frozen receptor revision")
    if method == "p2rank":
        residues, atoms = p2rank._input_residues(receptor)
        if (record.get("input_atom_count") != atoms
                or record.get("input_amino_acid_residue_count") != len(residues)
                or record.get("reference_independent") is not True):
            raise ValueError("P2Rank recorded full-receptor input audit differs")
    elif record.get("information_regime") != "full_surface_site_agnostic":
        raise ValueError("fpocket recorded information regime differs")
    if "predictions" not in record or not isinstance(record["predictions"], list):
        raise ValueError("site attempt has no normalized candidate list")
    if not isinstance(record.get("returncode"), int) or isinstance(record["returncode"], bool):
        raise ValueError("site attempt has no actual process return code")
    logs = {name: {"sha256": sha256(_regular(restored / name)), "bytes": (restored / name).stat().st_size}
            for name in ("stdout.log", "stderr.log")}
    if method == "fpocket":
        if installation is not None and record.get("installation") != installation:
            raise ValueError("fpocket attempt differs from the archived installation/build evidence")
        if sha256(restored / "receptor.cif") != sha256(receptor):
            raise ValueError("fpocket embedded input differs from its separate input archive")
        if sha256(restored / "receptor_preparation.json") != sha256(input_dir / "receptor.cif.manifest.json"):
            raise ValueError("fpocket embedded preparation receipt differs")
        if not initial:
            if (record.get("input_serialization_revision") != 2 or record.get("superseded_input_sha256") != original_sha
                    or sha256(restored / "input_revision.json") != sha256(input_dir / "receptor.cif.revision.json")):
                raise ValueError("corrected fpocket attempt lacks the frozen serialization revision")
    successful = record["status"] == ("completed" if method == "fpocket" else "ok")
    if successful:
        if initial or record["returncode"] != 0:
            raise ValueError("site success/initial-failure status contradicts the return code")
        if method == "fpocket":
            raw = restored / "receptor_out"
            predictions = sites.parse_fpocket(raw, receptor)
            if predictions:
                audit = sites.audit_fpocket_receptor(raw / "receptor_out.cif", receptor)
                if audit != record.get("receptor_output_audit"):
                    raise ValueError("fpocket full-receptor output audit differs")
            elif not record.get("no_pockets_reported"):
                raise ValueError("zero fpocket candidates lack an explicit no-pockets result")
            _same_predictions(record["predictions"], json.loads((restored / "predictions.json").read_text()))
            count_key = "pocket_count"
        else:
            raw = restored / "raw"
            predictions = p2rank.parse_p2rank_predictions(raw / "receptor.cif_predictions.csv",
                                                         raw / "receptor.cif_residues.csv", receptor)
            params = (raw / "params.txt").read_text()
            required = {"model": lock["parameters"]["model"], "threads": str(lock["parameters"]["threads"]),
                        "seed": str(lock["parameters"]["seed"]),
                        "visualizations": str(lock["parameters"]["visualizations"]).lower()}
            for key, value in required.items():
                if re.findall(r"(?m)^\s*" + re.escape(key) + r"\s*=\s*(\S+)\s*$", params) != [value]:
                    raise ValueError(f"P2Rank raw runtime parameter differs: {key}")
            count_key = "prediction_count"
        _same_predictions(record["predictions"], predictions)
        if record.get(count_key) != len(predictions):
            raise ValueError("recorded site candidate count differs from raw output")
    else:
        if record["status"] not in {"failed", "timeout", "cancelled", "canceled", "error"} or record["predictions"] != []:
            raise ValueError("failed site attempt must preserve an empty normalized candidate list")
        if initial and (record["returncode"] != 1 or "Could not get atom number" not in (restored / "stdout.log").read_text()
                        or "Structure reading failed" not in (restored / "stderr.log").read_text()):
            raise ValueError("initial fpocket failure is not the preserved CIF reader failure")
        predictions = []
    negative_code = record["returncode"] < 0
    return {"method": method, "execution_status": record["status"], "returncode": record["returncode"],
            "termination_signal": signal.Signals(-record["returncode"]).name if negative_code else None,
            "prediction_count": len(predictions), "raw_predictions_reparsed": successful,
            "record_sha256": archive["record_sha256"], "archive_sha256": archive["archive_sha256"],
            "artifact_count": len(archive["entries"]) - 1, "input_sha256": record["input_sha256"],
            "logs": logs, "error_type": record.get("error_type"), "error": record.get("error", record.get("reason")),
            "interpretation": "normalized candidates verified against raw output" if successful else "preserved workflow failure; no negative-site biological inference"}


def verify_site_predictions(root, *, output_path=None):
    """Verify all frozen predictor conditions plus the original CIF failures.

    Native-site truth and installed predictor binaries are not accessed. Any
    missing, corrupted, re-ranked, reparameterized or mismapped result raises
    before a successful verification record can be written.
    """
    root = Path(root).resolve()
    config_path = root / "benchmarks/site_experiment.json"
    config = json.loads(config_path.read_text())
    base = root / "results/site_discovery"
    freeze_path = base / "freeze_receipt.json"
    freeze = json.loads(freeze_path.read_text())
    if sha256(config_path) != freeze["configuration_sha256"]:
        raise ValueError("site experiment differs from its frozen configuration")
    sources = {}

    def track(path):
        path = Path(path)
        sources[path.relative_to(root).as_posix()] = sha256(path)

    track(config_path)
    track(freeze_path)
    for path, digest in config["method_locks"].items():
        if sha256(root / path) != digest or freeze["method_lock_sha256"].get(path) != digest:
            raise ValueError("site method lock differs from the frozen digest")
        track(root / path)
    fpocket_lock_path, p2rank_lock_path = root / "benchmarks/site_methods.lock.json", root / "benchmarks/p2rank.lock.json"
    fpocket_lock, p2rank_lock = json.loads(fpocket_lock_path.read_text()), json.loads(p2rank_lock_path.read_text())
    installation_dir = root / "results/site_methods/software/fpocket"
    installation_path = _regular(installation_dir / "installation.json")
    installation = json.loads(installation_path.read_text())
    fpocket_config = fpocket_lock["methods"]["fpocket"]
    if (installation.get("method") != "fpocket" or installation.get("algorithm_source_modifications") != []
            or installation.get("method_lock_sha256") != sha256(fpocket_lock_path)
            or any(installation.get(key) != fpocket_config[key] for key in ("release", "commit", "source_sha256"))):
        raise ValueError("fpocket archived installation differs from the frozen source/build policy")
    _digest(installation["executable_sha256"])
    track(installation_path)
    for filename, digest in installation["build_evidence"].items():
        evidence = _regular(installation_dir / _filename(filename))
        if sha256(evidence) != _digest(digest):
            raise ValueError("fpocket archived build evidence checksum differs")
        track(evidence)
    if any(installation["build_evidence"].get(filename) != digest
           for filename, digest in fpocket_config["license_files"].items()):
        raise ValueError("fpocket archived license differs from its frozen source")
    revision_path = base / "fpocket_serialization_revision.json"
    revision = json.loads(revision_path.read_text())
    if (revision.get("scientific_parameter_changes") is not False
            or revision.get("prediction_inputs_have_identical_protein_atoms_coordinates_and_author_ids") is not True
            or revision.get("method_lock_sha256") != sha256(fpocket_lock_path)):
        raise ValueError("global fpocket compatibility amendment changes scientific parameters")
    track(revision_path)
    conditions = [case["condition_id"] for case in config["cases"]]
    if len(conditions) != len(set(conditions)):
        raise ValueError("duplicate frozen site conditions")
    revised = {entry["condition_id"]: entry for entry in revision["inputs"]}
    if len(revised) != len(revision["inputs"]) or set(revised) != set(conditions):
        raise ValueError("fpocket serialization revision denominator differs")
    corrected_rows, initial_rows, input_rows = [], [], []
    with tempfile.TemporaryDirectory(prefix="aav-site-prediction-verification-") as temporary:
        temporary = Path(temporary)
        for case in config["cases"]:
            name = _filename(case["condition_id"])
            if freeze["input_source_sha256"].get(name) != case["source_sha256"]:
                raise ValueError("case source differs from its original freeze receipt")
            original, corrected = temporary / name / "original", temporary / name / "corrected"
            restore_site_inputs(root, name, "original", original)
            restore_site_inputs(root, name, "fpocket_compatibility", corrected)
            evidence = _input_evidence(case, original, corrected)
            if (revised[name]["original_input_sha256"] != evidence["original_sha256"]
                    or revised[name]["corrected_input_sha256"] != evidence["corrected_sha256"]):
                raise ValueError("case input bytes differ from the serialization amendment")
            input_rows.append({"condition_id": name, **evidence})
            track(base / name / "inputs/manifest.json")
            for method, directory, input_dir, initial in (("fpocket", "fpocket_v2", corrected, False),
                                                         ("p2rank", "p2rank", original, False),
                                                         ("fpocket", "fpocket", original, True)):
                run_dir = base / name / directory
                row = _verify_attempt(run_dir, temporary / name / (directory + "_raw"), input_dir,
                                      method=method, lock=fpocket_lock if method == "fpocket" else p2rank_lock,
                                      lock_sha=sha256(fpocket_lock_path if method == "fpocket" else p2rank_lock_path),
                                      initial=initial, original_sha=evidence["original_sha256"],
                                      installation=installation if method == "fpocket" else None)
                row.update(condition_id=name, attempt_role="initial_CIF_reader_failure" if initial else "corrected_method_condition",
                           source_record=(run_dir / ("execution.json" if method == "fpocket" else "result.json")).relative_to(root).as_posix())
                (initial_rows if initial else corrected_rows).append(row)
                for filename in ("raw_archive.json", "execution.json" if method == "fpocket" else "result.json"):
                    track(run_dir / filename)
    expected = len(conditions) * 2
    if len(corrected_rows) != expected or len(initial_rows) != len(conditions):
        raise ValueError("site prediction verification denominator is incomplete")
    result = {"schema_version": 1, "status": "verified", "scientific_tools_executed": False,
              "native_ligand_or_site_evaluation_access": False, "center_comparison_tolerance_A": CENTER_TOLERANCE_A,
              "all_other_prediction_fields": "exact equality, including scores, rank, IDs, footprints, counts and raw references",
              "eligible_conditions": len(conditions), "verified_corrected_method_conditions": len(corrected_rows),
              "verified_initial_CIF_reader_failures": len(initial_rows),
              "corrected_completed_conditions": sum(row["raw_predictions_reparsed"] for row in corrected_rows),
              "corrected_workflow_failures": sum(not row["raw_predictions_reparsed"] for row in corrected_rows),
              "inputs": input_rows, "corrected_method_conditions": corrected_rows, "initial_CIF_reader_failures": initial_rows,
              "source_sha256": sources, "verifier_sha256": sha256(Path(__file__)),
              "parser_sha256": {"fpocket": sha256(Path(sites.__file__)), "p2rank": sha256(Path(p2rank.__file__))},
              "scope": "archive integrity and raw-to-normalized prediction replay; workflow failures are retained, not biological negatives"}
    if output_path is not None:
        write_json(output_path, result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true", help="compare with the existing output without rewriting it")
    args = parser.parse_args(argv)
    if args.check and args.output is None:
        parser.error("--check requires --output")
    output = args.output if args.output is None or args.output.is_absolute() else args.root / args.output
    result = verify_site_predictions(args.root, output_path=None if args.check else output)
    if args.check and result != json.loads(output.read_text()):
        raise ValueError("site prediction verification record is stale")
    if output is None:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    else:
        print(f"Verified {result['verified_corrected_method_conditions']} corrected method conditions and "
              f"{result['verified_initial_CIF_reader_failures']} initial CIF failures; "
              f"{'checked' if args.check else 'wrote'} {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
