"""Frozen, incremental scoring of the unchanged retained Vina pose pool.

This module never forwards experimental atom coordinates or computes accuracy.
It reuses the source workflow's execution verifier, extracts only names/source
keys from archived chemical maps, and retains every source pose identity.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from .additional_method import (
    export_vina_carb_ligand, export_vina_carb_pose, validate_installation,
    expected_chi_links, verify_chi_recognition,
)
from .provenance import sha256, write_json, write_text_atomic


ROOT = Path(__file__).resolve().parents[1]
CONFIG = "benchmarks/common_pose_rescoring.json"
RECEIPT = "results/common_pose_rescoring/freeze_receipt.json"
ELIGIBLE = ["1GNY", "1KJL", "1QFO", "2J1T", "1PWB", "2VUZ"]
SCORE_FIELDS = {"glycosidic_links", "fully_supported_links", "intramolecular_reference", "base_shared0",
                "chi_coefficient0", "chi_coefficient1", "combined_shared0", "original_score_only_base",
                "pose_intramolecular_energy", "torsion_normalization_scale", "nonbonded_cutoff_A"}


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def validate_plan(root=ROOT):
    root = Path(root).resolve()
    config_path = root / CONFIG
    config = json.loads(config_path.read_text())
    freeze = json.loads((root / RECEIPT).read_text())
    if (sha256(config_path) != freeze["configuration_sha256"]
            or config.get("frozen_before_benchmark_rescoring") is not True
            or freeze.get("frozen_before_benchmark_rescoring") is not True
            or config.get("benchmark_rescoring_outcomes_consulted") is not False
            or freeze.get("new_rescoring_outcomes_consulted") is not False):
        raise ValueError("Common-pool plan differs from its pre-rescoring freeze")
    for key in ("method_lock", "atomic_reference_ledger", "source_experiment"):
        path = Path(config[key])
        if path.is_absolute() or ".." in path.parts or sha256(root / path) != config[key + "_sha256"]:
            raise ValueError("A frozen method, source experiment or reference scope changed")
        if key + "_sha256" in freeze and freeze[key + "_sha256"] != config[key + "_sha256"]:
            raise ValueError("Common-pool receipt and plan disagree")
    if (config["eligible_complexes"] != ELIGIBLE or config["intramolecular_reference"] != 0
            or config["chi_coefficient"] != 1 or config["paired_control_chi_coefficient"] != 0
            or config["chi_cutoff_kcal_mol"] != 2):
        raise ValueError("Frozen cohort or paired scoring convention changed")
    return config


def model_blocks(payload):
    """Hash complete raw MODEL blocks, including every atom and annotation."""
    rows, current, rank = {}, None, None
    for line in payload.splitlines(keepends=True):
        if line.startswith(b"MODEL"):
            if current is not None:
                raise ValueError("Nested source pose MODEL")
            rank = int(line.split()[1])
            if rank != len(rows) + 1:
                raise ValueError("Source pose ranks are not contiguous")
            current = [line]
        elif current is not None:
            current.append(line)
            if line.startswith(b"ENDMDL"):
                rows[rank] = b"".join(current)
                current = None
        elif line.startswith((b"ATOM  ", b"HETATM", b"ENDMDL")):
            raise ValueError("Source atom/model terminator outside MODEL")
    if current is not None or not rows:
        raise ValueError("Missing or truncated source pose models")
    return rows


def _plain_file(value):
    if not isinstance(value, str) or Path(value).name != value or value in {".", ".."}:
        raise ValueError("Source archive requires plain filenames")
    return value


def source_atom_map(root, case_id):
    """Return source identities only; native-coordinate fields are discarded."""
    root = Path(root)
    directory = root / "results/development" / case_id / "inputs"
    if not directory.exists():
        return {"case_id": case_id, "status": "unsupported_missing_parent_preparation", "atom_map": [],
                "native_coordinate_fields_forwarded": False}
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    matches = [entry for entry in manifest["files"] if entry["filename"] == "reference.json"]
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous source atom-map archive")
    entry = matches[0]
    archive = directory / _plain_file(entry["archive"])
    if archive.is_symlink() or sha256(archive) != entry["archive_sha256"]:
        raise ValueError("Source atom-map compressed bytes changed")
    raw = gzip.decompress(archive.read_bytes())
    if hashlib.sha256(raw).hexdigest() != entry["sha256"]:
        raise ValueError("Source atom-map decompressed bytes changed")
    # No coordinate field from this old document is returned or scored.
    identities = json.loads(raw)["atom_map"]
    atoms = [{"name": entry["name"], "source_key": entry["source_key"]} for entry in identities]
    if (not atoms or len({a["name"] for a in atoms}) != len(atoms)
            or len({tuple(a["source_key"]) for a in atoms}) != len(atoms)
            or any(len(a["source_key"]) != 5 for a in atoms)):
        raise ValueError("Source identity map is incomplete or duplicated")
    return {"case_id": case_id, "status": "mapped", "atom_map": atoms,
            "native_coordinate_fields_forwarded": False, "archive_manifest_sha256": sha256(manifest_path),
            "source_archive": str(archive.relative_to(root)), "source_archive_sha256": entry["archive_sha256"],
            "source_reference_sha256": entry["sha256"]}


def _source_context(root):
    # Execution verification is separate from qualified_pose's native-reference
    # and RMSD routines. This call neither loads its atomic ledger nor scores it.
    from .qualified_pose import Evidence, _all_execution
    evidence = Evidence(root)
    config = evidence.json("benchmarks/blind_pose_experiment.json")
    freeze = evidence.json("results/blind_pose/freeze_receipt.json")
    if evidence.sources["benchmarks/blind_pose_experiment.json"]["sha256"] != freeze["configuration_sha256"]:
        raise ValueError("Source execution plan changed")
    inventory = evidence.json("results/blind_pose/planning_inventory.json", freeze["planning_inventory_sha256"])
    oracle = evidence.json("results/blind_pose/oracle_controls.json", freeze["oracle_controls_sha256"])
    conditions, pending, denominator = _all_execution(evidence, config, freeze, inventory, oracle)
    return evidence, config, conditions, pending, denominator


def _pose_id(condition, box, seed, rank):
    return f"{condition['condition_id']}/{condition['method']}/{box['id']}/s{seed}/p{rank:02d}"


def merge_pose_index(previous, current):
    """Append newly terminal source poses without rewriting prior identities."""
    old = {row["pose_id"]: row for row in previous}
    new = {row["pose_id"]: row for row in current}
    if len(old) != len(previous) or len(new) != len(current):
        raise ValueError("Duplicate common-pool pose identity")
    if old.keys() - new.keys():
        raise ValueError("A previously indexed source pose disappeared")
    for identity, row in old.items():
        if row != new[identity]:
            raise ValueError(f"Previously indexed source pose changed: {identity}")
    return sorted(new.values(), key=lambda row: row["pose_id"])


def index_source_pool(root=ROOT, *, output_root=None):
    root = Path(root).resolve()
    plan = validate_plan(root)
    output = Path(output_root).resolve() if output_root else root / "results/common_pose_rescoring"
    output.mkdir(parents=True, exist_ok=True)
    evidence, source_config, conditions, pending, denominator = _source_context(root)
    maps = {}
    for case in plan["eligible_complexes"]:
        mapping = source_atom_map(root, case)
        relative = Path("source_atom_maps") / (case + ".json")
        target = output / relative
        if target.exists() and json.loads(target.read_text()) != mapping:
            raise ValueError("Previously frozen source identity map changed")
        if not target.exists():
            write_json(target, mapping)
        maps[case] = {"file": relative.as_posix(), "sha256": sha256(target), "status": mapping["status"]}
    poses, condition_rows, terminal_subset = [], [], 0
    for condition in conditions:
        if condition["complex_id"] not in plan["eligible_complexes"]:
            continue
        directory = Path("results/blind_pose") / condition["condition_id"] / condition["method"]
        c_row = {"condition_id": condition["condition_id"], "complex_id": condition["complex_id"],
                 "method": condition["method"], "preparation_status": condition["preparation_status"], "boxes": []}
        for box in condition["boxes"]:
            relative = directory / box["id"]
            box_row = {"box_id": box["id"], "candidate_class_rank": box["candidate_class_rank"],
                       "candidate_original_rank": box["candidate_original_rank"], "status": box.get("status", "prepared"),
                       "runs": []}
            if box.get("status") == "preparation_failed":
                c_row["boxes"].append(box_row)
                continue
            files = evidence.bundle(relative)
            for run in box["runs"]:
                box_row["runs"].append({"seed": run["seed"], "status": run["status"]})
                if run["status"] not in {"ok", "failed"}:
                    continue
                terminal_subset += 1
                if run["status"] != "ok":
                    continue
                raw = evidence.read(relative / run["pose_file"], run["pose_sha256"])
                blocks = model_blocks(raw)
                models = box["poses_by_seed"][run["seed"]]
                if len(blocks) != len(models) or len(models) != run["returned_poses"]:
                    raise ValueError("Source retained-model denominator differs")
                for model in models:
                    row = {"pose_id": _pose_id(condition, box, run["seed"], model["rank"]),
                           "condition_id": condition["condition_id"], "complex_id": condition["complex_id"],
                           "method": condition["method"], "box_id": box["id"], "seed": run["seed"],
                           "candidate_class_rank": box["candidate_class_rank"], "candidate_original_rank": box["candidate_original_rank"],
                           "model_rank": model["rank"], "original_vina_score": model["vina_score"],
                           "source_pose_file": (relative / run["pose_file"]).as_posix(),
                           "source_pose_file_sha256": run["pose_sha256"],
                           "source_model_sha256": hashlib.sha256(blocks[model["rank"]]).hexdigest(),
                           "source_run_record_sha256": canonical_sha256(run),
                           "source_box": relative.as_posix(), "source_atom_map_sha256": maps[condition["complex_id"]]["sha256"],
                           "source_ligand_sha256": hashlib.sha256(files["ligand.pdbqt"]).hexdigest(),
                           "source_receptor_sha256": hashlib.sha256(files["neighborhood.pdbqt"]).hexdigest(),
                           "source_preparation_sha256": hashlib.sha256(files["preparation.json"]).hexdigest()}
                    row["pose_identity_sha256"] = canonical_sha256(row)
                    poses.append(row)
            c_row["boxes"].append(box_row)
        condition_rows.append(c_row)
    if len(condition_rows) != 18:
        raise ValueError("Eligible common-pool condition denominator changed")
    target = output / "source_pool.json"
    if target.exists():
        previous = json.loads(target.read_text())
        if previous["configuration_sha256"] != sha256(root / CONFIG):
            raise ValueError("An existing output uses a different frozen common-pool configuration")
        poses = merge_pose_index(previous["poses"], poses)
    else:
        poses = merge_pose_index([], poses)
    index = {"schema_version": 1, "configuration_sha256": sha256(root / CONFIG),
             "method_lock_sha256": plan["method_lock_sha256"], "atomic_reference_ledger_sha256": plan["atomic_reference_ledger_sha256"],
             "source_experiment_sha256": plan["source_experiment_sha256"], "conditions": condition_rows,
             "eligible_complexes": plan["eligible_complexes"], "eligible_conditions": 18,
             "source_global_denominator": denominator, "source_global_pending": pending,
             "source_global_complete_and_verified": not pending and denominator["terminal_seed_runs"] == 156,
             "atomic_subset_terminal_runs": terminal_subset, "source_atom_maps": maps,
             "retained_pose_count": len(poses), "poses": poses,
             "native_coordinates_used_for_scoring": False, "atomic_accuracy_computed": False}
    write_json(target, index)
    preflight_path = output / "case_preflight.json"
    preflight = json.loads(preflight_path.read_text()) if preflight_path.exists() else {}
    for case in index["eligible_complexes"]:
        mapping = json.loads((output / maps[case]["file"]).read_text())
        candidates = [row for row in poses if row["complex_id"] == case]
        if case in preflight and preflight[case].get("source_box"):
            first = next((row for row in candidates if row["source_box"] == preflight[case]["source_box"]), None)
            if first is None:
                raise ValueError("The preflight source input disappeared")
        else:
            first = next(iter(candidates), None)
        if first is not None or mapping["status"] != "mapped":
            decision = _case_preflight(root, evidence, first, mapping)
            if case in preflight and preflight[case] != decision:
                raise ValueError("Static chemical support changed during indexing")
            preflight[case] = decision
    write_json(preflight_path, preflight)
    return index


def _record_path(output, pose_id):
    if Path(pose_id).is_absolute() or ".." in Path(pose_id).parts or not re.fullmatch(r"[A-Za-z0-9_/-]+", pose_id):
        raise ValueError("Unsafe common-pool pose identity")
    return Path(output) / "scores" / (pose_id + ".json")


def verify_score_record(record, source, method_lock_sha256):
    if (record["source"] != source or record["source"]["pose_identity_sha256"] != source["pose_identity_sha256"]
            or record["method_lock_sha256"] != method_lock_sha256
            or record.get("coordinate_optimization") is not False or record.get("native_coordinates_used") is not False):
        raise ValueError("Rescore record source, method or coordinate policy changed")
    if record["status"] != "ok":
        if record["status"] not in {"unsupported_case", "pose_scoring_failed"} or not record.get("error"):
            raise ValueError("Invalid common-pool failure record")
        if "native_attempt" in record:
            attempt = record["native_attempt"]
            _verify_raw_streams(attempt)
            if (attempt["ligand_sha256"] != record["exported_ligand_sha256"]
                    or attempt["receptor_sha256"] != source["source_receptor_sha256"]):
                raise ValueError("Failed native scoring attempt changed input identities")
        return record
    score = record["score"]
    _verify_raw_streams(score)
    parsed = json.loads(score["raw_stdout"].strip().splitlines()[-1])
    if set(parsed) != SCORE_FIELDS or any(parsed[key] != score.get(key) for key in parsed):
        raise ValueError("Parsed fixed-pose score differs from raw output")
    if (score["ligand_sha256"] != record["exported_ligand_sha256"]
            or score["receptor_sha256"] != source["source_receptor_sha256"]
            or score["intramolecular_reference"] != 0 or score["chi_coefficient0"] != 0
            or score["nonbonded_cutoff_A"] != 8 or not 0 < score["torsion_normalization_scale"] <= 1
            or score["glycosidic_links"] != record["expected_glycosidic_links"]
            or score["fully_supported_links"] != record["expected_glycosidic_links"]
            or score["coordinate_optimization"] is not False or score["native_coordinates_accessed"] is not False
            or abs(score["combined_shared0"] - score["base_shared0"] - score["chi_coefficient1"]) > 1e-9
            or not math.isclose(score["base_shared0"] - score["original_score_only_base"],
                                score["torsion_normalization_scale"] * score["pose_intramolecular_energy"],
                                rel_tol=1e-9, abs_tol=1e-8)
            or any(isinstance(v, float) and not math.isfinite(v) for v in score.values())):
        raise ValueError("Fixed-pose coefficient, linkage, normalization or source control failed")
    return record


def _verify_raw_streams(record):
    for field in ("raw_stdout", "raw_stderr", "recognition_log"):
        if hashlib.sha256(record[field].encode()).hexdigest() != record[field + "_sha256"]:
            raise ValueError("Raw fixed-pose scoring stream changed: " + field)


class NativeScoringFailure(RuntimeError):
    def __init__(self, message, attempt):
        super().__init__(message)
        self.attempt = attempt


def _invoke_fixed_pose(receptor, ligand, export_receipt, installation, *, root):
    """Use the frozen helper while preserving streams even when validation fails.

    The invocation and controls match additional_method.score_fixed_pose; this
    wrapper adds failure capture without changing the frozen module or helper.
    """
    installation = Path(installation).resolve()
    validate_installation(installation, root=root)
    scorer = json.loads((installation / "fixed_pose_scorer.json").read_text())
    executable = installation / scorer["executable"]
    if (sha256(executable) != scorer["executable_sha256"]
            or sha256(ligand) != export_receipt["ligand_sha256"]
            or export_receipt["native_coordinates_accessed"] is not False
            or scorer["method_lock_sha256"] != sha256(root / "benchmarks/additional_method.lock.json")):
        raise ValueError("Frozen fixed-pose scorer, ligand or method changed")
    with tempfile.TemporaryDirectory(prefix="vc-fixed-pose-") as temporary:
        command = [str(executable), str(Path(receptor).resolve()), str(Path(ligand).resolve())]
        try:
            result = subprocess.run(command, cwd=temporary, env={**os.environ, "VINA_CARB": str(installation / "VC_1_0")},
                                    capture_output=True, text=True, timeout=120)
            stdout, stderr, returncode, timed_out = result.stdout, result.stderr, result.returncode, False
        except subprocess.TimeoutExpired as error:
            stdout, stderr, returncode, timed_out = error.stdout or "", error.stderr or "", None, True
            if isinstance(stdout, bytes):
                stdout = stdout.decode(errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode(errors="replace")
        recognition_path = Path(temporary) / "VC_log.txt"
        recognition_text = recognition_path.read_text() if recognition_path.exists() else ""
        attempt = {"command": command, "executable_sha256": scorer["executable_sha256"],
                   "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
                   "returncode": returncode, "timed_out": timed_out, "raw_stdout": stdout,
                   "raw_stderr": stderr, "recognition_log": recognition_text}
        attempt.update({key + "_sha256": hashlib.sha256(attempt[key].encode()).hexdigest()
                        for key in ("raw_stdout", "raw_stderr", "recognition_log")})
        if returncode != 0 or timed_out:
            raise NativeScoringFailure("Frozen fixed-pose scorer failed or timed out", attempt)
        try:
            parsed = json.loads(stdout.strip().splitlines()[-1])
            if set(parsed) != SCORE_FIELDS:
                raise ValueError("Unexpected frozen helper output fields")
            recognition = verify_chi_recognition(recognition_path, expected_chi_links(export_receipt))
            score = {**parsed, **attempt, "recognition": recognition,
                     "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
                     "coordinate_optimization": False, "native_coordinates_accessed": False,
                     "score_role": "VC1.0 original scoring functions with common zero intramolecular reference; rank only within one fixed ligand/torsion-tree pool"}
            # The enclosing record validator checks every numerical control too.
            return score
        except (ValueError, IndexError, KeyError, OSError) as error:
            raise NativeScoringFailure(str(error), attempt) from error


def _old_records(output, index, *, recover_orphans=False):
    manifest_path = Path(output) / "scores_manifest.json"
    if not manifest_path.exists():
        if not recover_orphans and (Path(output) / "scores").exists() and any((Path(output) / "scores").rglob("*.json")):
            raise ValueError("Score artifacts exist without their checked manifest")
        manifest = {"schema_version": 1, "configuration_sha256": index["configuration_sha256"], "records": {}}
    else:
        manifest = json.loads(manifest_path.read_text())
    if manifest["configuration_sha256"] != index["configuration_sha256"]:
        raise ValueError("Rescoring configuration changed during resume")
    sources = {row["pose_id"]: row for row in index["poses"]}
    records = {}
    for identity, digest in manifest["records"].items():
        path = _record_path(output, identity)
        if identity not in sources or sha256(path) != digest:
            raise ValueError("A completed score artifact disappeared or changed")
        record = json.loads(path.read_text())
        records[identity] = verify_score_record(record, sources[identity], index["method_lock_sha256"])
    existing = {p.relative_to(Path(output) / "scores").as_posix()[:-5] for p in (Path(output) / "scores").rglob("*.json")}
    if existing != records.keys() and not recover_orphans:
        raise ValueError("Unexpected score artifacts outside the manifest")
    for identity in sorted(existing - records.keys()):
        # A process can stop after the atomic record write and before the
        # manifest update. Recover only a complete, source-bound valid record.
        if identity not in sources:
            raise ValueError("Orphan score is outside the frozen source pool")
        path = _record_path(output, identity)
        records[identity] = verify_score_record(json.loads(path.read_text()), sources[identity], index["method_lock_sha256"])
        manifest["records"][identity] = sha256(path)
    return records, manifest


def _scorer_provenance(root, output, installation, method_lock_sha256):
    installation = Path(installation).resolve()
    scorer = json.loads((installation / "fixed_pose_scorer.json").read_text())
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    if (scorer["method_lock_sha256"] != method_lock_sha256
            or scorer["caller_source_sha256"] != lock["fixed_pose_rescorer"]["source_sha256"]
            or sha256(installation / scorer["executable"]) != scorer["executable_sha256"]):
        raise ValueError("Fixed-pose executable does not match the frozen source scoring caller")
    record = {"schema_version": 1, "method_lock_sha256": method_lock_sha256,
              "caller_source_sha256": scorer["caller_source_sha256"], "executable_sha256": scorer["executable_sha256"],
              "installation_receipt_sha256": sha256(installation / "installation.json"),
              "underlying_scoring_source_changed": False, "driver_source_sha256": sha256(Path(__file__))}
    archived = {"driver_source.py": Path(__file__).read_bytes(),
                "native_installation.json": (installation / "installation.json").read_bytes(),
                "fixed_pose_scorer.json": (installation / "fixed_pose_scorer.json").read_bytes()}
    record["archived_artifacts"] = {}
    for name, payload in archived.items():
        target = output / "scorer_source" / name
        if target.exists() and target.read_bytes() != payload:
            raise ValueError("Archived scoring driver or native receipt changed during resume")
        if not target.exists():
            write_text_atomic(target, payload.decode("utf-8"))
        record["archived_artifacts"][str(target.relative_to(output))] = hashlib.sha256(payload).hexdigest()
    path = output / "scorer_provenance.json"
    if path.exists() and json.loads(path.read_text()) != record:
        raise ValueError("Scoring executable/provenance changed during resume")
    if not path.exists():
        write_json(path, record)
    return record


def _case_preflight(root, evidence, row, mapping):
    if mapping["status"] != "mapped":
        return {"status": "unsupported_case", "error": "No valid full-parent chemical preparation is archived"}
    files = evidence.bundle(row["source_box"])
    with tempfile.TemporaryDirectory(prefix="vc-pool-chemistry-") as temporary:
        path = Path(temporary) / "source.pdbqt"
        path.write_bytes(files["ligand.pdbqt"])
        try:
            exported = export_vina_carb_ligand(path, mapping["atom_map"], Path(temporary) / "exported.pdbqt")
            residue_count = len({tuple(atom["source_key"][:4]) for atom in mapping["atom_map"]})
            if len(exported["expected_glycosidic_links"]) != residue_count - 1:
                raise ValueError("Connected source glycan does not expose every inter-residue linkage as a supported torsion")
        except ValueError as error:
            return {"status": "unsupported_case", "error": str(error), "error_type": type(error).__name__}
    return {"status": "supported", "expected_glycosidic_links": len(exported["expected_glycosidic_links"]),
            "source_box": row["source_box"], "source_ligand_sha256": row["source_ligand_sha256"]}


def _score_one(root, output, source, evidence, mapping, installation, case_status, method_lock_sha256):
    base = {"schema_version": 1, "source": source, "method_lock_sha256": method_lock_sha256,
            "coordinate_optimization": False, "native_coordinates_used": False}
    if case_status["status"] != "supported":
        return {**base, "status": "unsupported_case", "error": case_status["error"]}
    files = evidence.bundle(source["source_box"])
    raw_pose = (root / source["source_pose_file"]).read_bytes()
    if (hashlib.sha256(raw_pose).hexdigest() != source["source_pose_file_sha256"]
            or hashlib.sha256(model_blocks(raw_pose)[source["model_rank"]]).hexdigest() != source["source_model_sha256"]):
        raise ValueError("Source pose changed after indexing")
    with tempfile.TemporaryDirectory(prefix="vc-common-pose-") as temporary:
        temporary = Path(temporary)
        prepared, receptor, poses, exported = [temporary / name for name in ("source.pdbqt", "receptor.pdbqt", "poses.pdbqt", "exported.pdbqt")]
        prepared.write_bytes(files["ligand.pdbqt"])
        receptor.write_bytes(files["neighborhood.pdbqt"])
        poses.write_bytes(raw_pose)
        try:
            receipt = export_vina_carb_pose(prepared, mapping["atom_map"], poses, source["model_rank"], exported)
            score = _invoke_fixed_pose(receptor, exported, receipt, installation, root=root)
            record = {**base, "status": "ok", "score": score, "exported_ligand_sha256": sha256(exported),
                      "export_map_sha256": sha256(exported.with_suffix(".pdbqt.map.json")),
                      "expected_glycosidic_links": len(receipt["expected_glycosidic_links"])}
            try:
                verify_score_record(record, source, method_lock_sha256)
            except (ValueError, KeyError, TypeError) as error:
                raise NativeScoringFailure(str(error), score) from error
        except NativeScoringFailure as error:
            record = {**base, "status": "pose_scoring_failed", "error": str(error),
                      "error_type": type(error).__name__, "native_attempt": error.attempt,
                      "exported_ligand_sha256": sha256(exported),
                      "export_map_sha256": sha256(exported.with_suffix(".pdbqt.map.json")),
                      "expected_glycosidic_links": len(receipt["expected_glycosidic_links"])}
        except (ValueError, RuntimeError) as error:
            record = {**base, "status": "pose_scoring_failed", "error": str(error), "error_type": type(error).__name__}
    return verify_score_record(record, source, method_lock_sha256)


def rank_complete_pools(index, records, case_status):
    """Rank only complete pools; failed poses are never removed to improve rank."""
    sources = index["poses"]
    result = []
    normalizers = {}
    for record in records.values():
        if record["status"] == "ok":
            normalizers.setdefault(record["source"]["complex_id"], []).append(record["score"]["torsion_normalization_scale"])
    conflicts = {case for case, values in normalizers.items() if max(values) - min(values) > 1e-12}
    for condition in index["conditions"]:
        for budget in (1, 5):
            for seed in (104729, 130363, 155921):
                selected = [row for row in sources if row["condition_id"] == condition["condition_id"] and row["method"] == condition["method"]
                            and row["seed"] == seed and row["candidate_class_rank"] <= budget]
                selected.sort(key=lambda row: (row["original_vina_score"], row["candidate_class_rank"], row["model_rank"]))
                source_pending = any(not any(run["seed"] == seed and run["status"] in {"ok", "failed"} for run in box["runs"])
                                     for box in condition["boxes"] if box["candidate_class_rank"] <= budget and box["status"] != "preparation_failed")
                identity = {key: condition[key] for key in ("condition_id", "complex_id", "method")}
                row = {**identity, "seed": seed, "site_budget": budget, "eligible_workflow_pool": True,
                       "source_retained_pose_count": len(selected), "original_vina_order": [r["pose_id"] for r in selected],
                       "coefficient0_order": [], "coefficient1_order": [], "atomic_accuracy_computed": False}
                if source_pending:
                    status = "pending_source"
                elif not selected:
                    status = "source_workflow_failure_no_poses"
                elif condition["complex_id"] in conflicts:
                    status = "unsupported_case_normalization_conflict"
                elif case_status.get(condition["complex_id"], {}).get("status") == "unsupported_case":
                    status = "unsupported_case"
                elif any(r["pose_id"] not in records for r in selected):
                    status = "pending_scoring"
                elif any(records[r["pose_id"]]["status"] != "ok" for r in selected):
                    status = "pose_scoring_failure_complete_pool_required"
                else:
                    status = "completed"
                    for coefficient, field in ((0, "base_shared0"), (1, "combined_shared0")):
                        order = sorted(selected, key=lambda source: (records[source["pose_id"]]["score"][field], source["candidate_class_rank"], source["model_rank"]))
                        row[f"coefficient{coefficient}_order"] = [r["pose_id"] for r in order]
                row["status"] = status
                result.append(row)
    return result, sorted(conflicts)


def rescore_available(root=ROOT, *, installation, output_root=None, max_poses=None):
    """Process newly terminal poses with one native worker and checked resume."""
    root = Path(root).resolve()
    output = Path(output_root).resolve() if output_root else root / "results/common_pose_rescoring"
    if max_poses is not None and (isinstance(max_poses, bool) or not isinstance(max_poses, int) or max_poses < 1):
        raise ValueError("Execution shard limit must be a positive integer")
    index = index_source_pool(root, output_root=output)
    validate_installation(installation, root=root)
    scorer = _scorer_provenance(root, output, installation, index["method_lock_sha256"])
    evidence, _, _, _, _ = _source_context(root)
    records, manifest = _old_records(output, index, recover_orphans=True)
    write_json(output / "scores_manifest.json", manifest)
    if any(record["status"] == "ok" and record["score"]["executable_sha256"] != scorer["executable_sha256"] for record in records.values()):
        raise ValueError("Previously scored poses used a different native scorer")
    mappings = {case: json.loads((output / row["file"]).read_text()) for case, row in index["source_atom_maps"].items()}
    case_status_path = output / "case_preflight.json"
    case_status = json.loads(case_status_path.read_text()) if case_status_path.exists() else {}
    for case in index["eligible_complexes"]:
        recorded_box = case_status.get(case, {}).get("source_box")
        first = next((row for row in index["poses"] if row["complex_id"] == case
                      and (recorded_box is None or row["source_box"] == recorded_box)), None)
        if first is not None and case not in case_status:
            case_status[case] = _case_preflight(root, evidence, first, mappings[case])
        elif mappings[case]["status"] != "mapped":
            case_status[case] = {"status": "unsupported_case", "error": "No valid full-parent chemical preparation is archived"}
    write_json(case_status_path, case_status)
    processed = 0
    for source in index["poses"]:
        identity = source["pose_id"]
        if identity in records:
            continue
        if max_poses is not None and processed >= max_poses:
            break
        record = _score_one(root, output, source, evidence, mappings[source["complex_id"]], installation,
                            case_status[source["complex_id"]], index["method_lock_sha256"])
        target = _record_path(output, identity)
        if target.exists():
            raise FileExistsError("Refusing to overwrite a prior score artifact")
        write_json(target, record)
        records[identity] = record
        manifest["records"][identity] = sha256(target)
        write_json(output / "scores_manifest.json", manifest)
        processed += 1
        if processed % 50 == 0:
            print(f"Scored {processed} new poses; {len(records)}/{len(index['poses'])} indexed records retained", flush=True)
    pools, conflicts = rank_complete_pools(index, records, case_status)
    result = {"schema_version": 1, "configuration_sha256": index["configuration_sha256"],
              "source_pool_sha256": sha256(output / "source_pool.json"), "score_manifest_sha256": sha256(output / "scores_manifest.json") if (output / "scores_manifest.json").exists() else None,
              "source_global_complete_and_verified": index["source_global_complete_and_verified"],
              "indexed_pose_count": len(index["poses"]), "recorded_pose_count": len(records), "new_records_this_pass": processed,
              "successful_pose_scores": sum(r["status"] == "ok" for r in records.values()),
              "unsupported_or_failed_pose_scores": sum(r["status"] != "ok" for r in records.values()),
              "pending_pose_scores": len(index["poses"]) - len(records), "case_preflight": case_status,
              "normalization_conflict_cases": conflicts, "pools": pools, "atomic_accuracy_computed": False,
              "ready_for_separate_accuracy_evaluation": index["source_global_complete_and_verified"] and len(records) == len(index["poses"])}
    write_json(output / "scoring_progress.json", result)
    return result


def verify_common_scores(root=ROOT, *, output_root=None, require_complete=True):
    """Rebuild source identities and rank all retained raw scores without binaries."""
    root = Path(root).resolve()
    output = Path(output_root).resolve() if output_root else root / "results/common_pose_rescoring"
    index = index_source_pool(root, output_root=output)
    records, manifest = _old_records(output, index)
    if require_complete and (not index["source_global_complete_and_verified"] or len(records) != len(index["poses"])):
        raise ValueError("Frozen source execution or common-pool scoring is incomplete")
    case_status = json.loads((output / "case_preflight.json").read_text()) if (output / "case_preflight.json").exists() else {}
    scorer = json.loads((output / "scorer_provenance.json").read_text()) if records else None
    if scorer:
        lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
        if (scorer["method_lock_sha256"] != index["method_lock_sha256"]
                or scorer["caller_source_sha256"] != lock["fixed_pose_rescorer"]["source_sha256"]
                or scorer["underlying_scoring_source_changed"] is not False):
            raise ValueError("Archived native scorer does not match the frozen source caller")
        for relative, digest in scorer["archived_artifacts"].items():
            if Path(relative).is_absolute() or ".." in Path(relative).parts or sha256(output / relative) != digest:
                raise ValueError("Archived scoring driver/native receipt changed")
        if (sha256(output / "scorer_source/driver_source.py") != scorer["driver_source_sha256"]
                or sha256(output / "scorer_source/native_installation.json") != scorer["installation_receipt_sha256"]):
            raise ValueError("Native scoring source/installation receipt binding changed")
    evidence, _, _, _, _ = _source_context(root)
    mappings = {case: json.loads((output / item["file"]).read_text()) for case, item in index["source_atom_maps"].items()}
    for case in index["eligible_complexes"]:
        recorded_box = case_status.get(case, {}).get("source_box")
        first = next((row for row in index["poses"] if row["complex_id"] == case
                      and (recorded_box is None or row["source_box"] == recorded_box)), None)
        if first is not None and _case_preflight(root, evidence, first, mappings[case]) != case_status.get(case):
            raise ValueError("Static chemical support decision changed")
    export_checks = 0
    for source in index["poses"]:
        record = records.get(source["pose_id"])
        if record is None or "exported_ligand_sha256" not in record:
            continue
        raw_score = record["score"] if record["status"] == "ok" else record["native_attempt"]
        if raw_score["executable_sha256"] != scorer["executable_sha256"]:
            raise ValueError("Native scorer identity changed within the fixed pool")
        files = evidence.bundle(source["source_box"])
        with tempfile.TemporaryDirectory(prefix="vc-export-replay-") as temporary:
            temporary = Path(temporary)
            prepared, poses, exported = [temporary / name for name in ("source.pdbqt", "poses.pdbqt", "exported.pdbqt")]
            prepared.write_bytes(files["ligand.pdbqt"])
            poses.write_bytes((root / source["source_pose_file"]).read_bytes())
            receipt = export_vina_carb_pose(prepared, mappings[source["complex_id"]]["atom_map"], poses, source["model_rank"], exported)
            if (sha256(exported) != record["exported_ligand_sha256"]
                    or sha256(exported.with_suffix(".pdbqt.map.json")) != record["export_map_sha256"]
                    or len(receipt["expected_glycosidic_links"]) != record["expected_glycosidic_links"]):
                raise ValueError("Scored coordinates/identities do not reconstruct from the frozen source pose")
            if record["status"] == "ok":
                recognition = temporary / "VC_log.txt"
                recognition.write_text(record["score"]["recognition_log"])
                verify_chi_recognition(recognition, expected_chi_links(receipt))
        export_checks += 1
    pools, conflicts = rank_complete_pools(index, records, case_status)
    if require_complete and any(row["status"].startswith("pending") for row in pools):
        raise ValueError("A prespecified source/score pool is incomplete")
    return {"schema_version": 1, "status": "verified", "native_executables_used": False,
            "configuration_sha256": index["configuration_sha256"], "source_pool_sha256": sha256(output / "source_pool.json"),
            "source_global_denominator": index["source_global_denominator"], "retained_pose_count": len(index["poses"]),
            "recorded_pose_count": len(records), "successful_pose_scores": sum(r["status"] == "ok" for r in records.values()),
            "unsupported_or_failed_pose_scores": sum(r["status"] != "ok" for r in records.values()),
            "pose_identity_exports_reconstructed": export_checks,
            "normalization_conflict_cases": conflicts, "pools": pools, "atomic_accuracy_computed": False,
            "score_manifest_sha256": sha256(output / "scores_manifest.json") if (output / "scores_manifest.json").exists() else None,
            "verifier_sha256": sha256(Path(__file__))}
