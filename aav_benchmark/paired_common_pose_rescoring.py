"""Source-only union of all retained paired Vina/VC poses, without deduplication.

Fixed scores are prohibited until a separate root freeze. Accuracy is never
computed here; both216 source-method slots must be terminal before evaluation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import resource
import time

from .additional_method import (
    _atom, export_vina_carb_pose, read_carb_poses, expected_chi_links, verify_chi_recognition,
)
from .common_pose_rescoring import (
    canonical_sha256, merge_pose_index, model_blocks, verify_score_record,
    _old_records, _record_path, _invoke_fixed_pose, NativeScoringFailure,
)
from .prospective_vina_carb import CASES, METHODS, PROFILES, SEEDS, normalize_atom_map
from .provenance import sha256, write_json, write_text_atomic

ROOT = Path(__file__).resolve().parents[1]
SOURCE_METHODS = ("modern_vina_1_2_7", "vina_carb_1_0")
EXECUTION_HELPERS = {"aav_benchmark/additional_method.py", "aav_benchmark/common_pose_rescoring.py",
                     "aav_benchmark/provenance.py", "aav_benchmark/prospective_vina_carb.py",
                     "aav_benchmark/paired_method_report.py", "aav_benchmark/heldout_pose_report.py",
                     "aav_benchmark/qualified_pose.py", "aav_benchmark/engine.py",
                     "results/additional_method_reproduction/build/fixed_pose_score.cpp"}
DEFAULT_PLAN = "benchmarks/structural_pose_v3_common_pool.proposal.json"
DEFAULT_OUTPUT = "results/structural_pose_v3/common_pose_union"


def _pin(root, entry):
    path = Path(entry["path"])
    if path.is_absolute() or ".." in path.parts or sha256(Path(root) / path) != entry["sha256"]:
        raise ValueError("A frozen union-pool source changed")
    return Path(root) / path


def validate_union_plan(root=ROOT, *, plan_path=DEFAULT_PLAN, receipt_path=None, require_freeze=False):
    root = Path(root).resolve();plan_path = Path(plan_path)
    if not plan_path.is_absolute():
        plan_path = root / plan_path
    plan = json.loads(plan_path.read_text())
    if require_freeze and plan.get("frozen_before_fixed_pose_scoring") is not True:
        raise ValueError("Separate root freeze required before any union fixed scoring")
    if (plan["cases"] != CASES or plan["profiles"] != PROFILES or plan["arms"] != METHODS or plan["seeds"] != SEEDS
            or plan["eligible_seed_pools"] != 108 or plan["eligible_source_method_slots"] != 216
            or plan["chi_coefficient"] != 1 or plan["paired_control_chi_coefficient"] != 0
            or plan["chi_cutoff_kcal_mol"] != 2 or plan["intramolecular_reference"] != 0
            or plan["deduplicate_or_filter_poses"] is not False):
        raise ValueError("Union cohort, exact-pool policy or scoring convention changed")
    for key in ("source_experiment", "source_freeze", "source_modern_experiment", "source_modern_freeze", "source_inputs", "source_schedule", "method_lock", "adapter_source"):
        _pin(root, plan[key])
    if plan["adapter_source"]["sha256"] != sha256(Path(__file__)):
        raise ValueError("Executing union adapter differs from its planned source")
    if set(plan["helper_sources"]) != EXECUTION_HELPERS:
        raise ValueError("Every unchanged scoring/export/execution helper needs an explicit source pin")
    for path, digest in plan["helper_sources"].items():
        _pin(root, {"path": path, "sha256": digest})
    modern_receipt = json.loads(_pin(root, plan["source_modern_freeze"]).read_text())
    if modern_receipt["configuration_sha256"] != plan["source_modern_experiment"]["sha256"]:
        raise ValueError("Modern source freeze and experiment disagree")
    if require_freeze:
        if receipt_path is None:
            raise ValueError("Separate root freeze required before any union fixed scoring")
        receipt_path = Path(receipt_path)
        if not receipt_path.is_absolute():
            receipt_path = root / receipt_path
        receipt = json.loads(receipt_path.read_text())
        if (plan.get("frozen_before_fixed_pose_scoring") is not True
                or receipt.get("frozen_before_fixed_pose_scoring") is not True
                or plan.get("new_fixed_score_or_accuracy_outcomes_consulted") is not False
                or receipt.get("new_fixed_score_or_accuracy_outcomes_consulted") is not False
                or receipt["configuration_sha256"] != sha256(plan_path)):
            raise ValueError("Separate root freeze required before any union fixed scoring")
    return plan, plan_path


def index_union_pool(root=ROOT, *, plan_path=DEFAULT_PLAN, output_root=None):
    """Index only checked computational source executions; never score accuracy."""
    from .prospective_vina_carb import validate_frozen_plan
    from .paired_method_report import read_execution_pair
    from .qualified_pose import Evidence
    root = Path(root).resolve();plan, actual_plan_path = validate_union_plan(root, plan_path=plan_path)
    experiment, conditions = validate_frozen_plan(root, _pin(root, plan["source_experiment"]), _pin(root, plan["source_freeze"]))
    if (experiment["prepared_inputs"] != plan["source_inputs"] or experiment["paired_modern_vina_policy"] != plan["source_modern_experiment"]
            or experiment["method_lock"] != plan["method_lock"]):
        raise ValueError("Union source experiments do not share the exact frozen36 inputs")
    modern = json.loads(_pin(root, plan["source_modern_experiment"]).read_text())
    output = Path(output_root) if output_root else root / DEFAULT_OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    evidence = Evidence(root);poses=[];pools=[];terminal=0;pending=[]
    for row in conditions:
        paired = read_execution_pair(evidence, row, modern, experiment, plan["source_experiment"], plan["source_freeze"],
                                     executable_sha256=experiment["native_binary_sha256"])
        pending.extend(row["condition_id"] + "/" + p for p in paired["pending"])
        for seed in paired["seeds"]:
            pool_id = row["condition_id"] + "/s" + str(seed["seed"])
            pool = {"pool_id": pool_id, "condition_id": row["condition_id"], "case_id": row["case_id"],
                    "receptor_profile": row["receptor_profile"], "arm": row["method"], "seed": seed["seed"],
                    "eligible_workflow_pool": True, "source_statuses": {}, "source_record_sha256": {},
                    "original_source_orders": {name: [] for name in SOURCE_METHODS}}
            for algorithm_rank, algorithm in enumerate(SOURCE_METHODS):
                result = seed["search_methods"][algorithm]
                pool["source_statuses"][algorithm] = result["status"] if result else "pending"
                if result is None:
                    continue
                terminal += 1
                pool["source_record_sha256"][algorithm] = canonical_sha256(result["record"])
                if result["status"] != "ok":
                    continue
                if algorithm == SOURCE_METHODS[0]:
                    path = Path(row["modern_runs"]).parent / result["record"]["pose_file"]
                    export_map = None
                else:
                    path = Path(row["vina_carb_outputs"][str(seed["seed"])] ) / "poses.pdbqt"
                    map_path = path.parent / "ligand.pdbqt.map.json"
                    export_map = {"path": map_path.as_posix(), "sha256": sha256(root / map_path)}
                raw = evidence.read(path);blocks = model_blocks(raw)
                if len(blocks) != len(result["poses"]):
                    raise ValueError("Raw source-model and validated retained-pose denominators differ")
                for pose in result["poses"]:
                    source = {"pose_id": pool_id + "/" + algorithm + f"/p{pose['rank']:02d}", "pool_id": pool_id,
                              "condition_id": row["condition_id"], "complex_id": row["case_id"], "receptor_profile": row["receptor_profile"],
                              "arm": row["method"], "seed": seed["seed"], "source_algorithm": algorithm,
                              "source_algorithm_rank": algorithm_rank, "candidate_class_rank": row["candidate_class_rank"],
                              "candidate_original_rank": row["candidate_original_rank"], "model_rank": pose["rank"],
                              "original_method_score": pose["method_score"], "source_pose_file": path.as_posix(),
                              "source_pose_file_sha256": hashlib.sha256(raw).hexdigest(),
                              "source_model_sha256": hashlib.sha256(blocks[pose["rank"]]).hexdigest(),
                              "source_record_sha256": canonical_sha256(result["record"]), "source_export_map": export_map,
                              "source_ligand_sha256": row["ligand"]["sha256"], "source_receptor_sha256": row["receptor"]["sha256"],
                              "source_preparation_sha256": row["preparation"]["sha256"],
                              "source_canonical_atom_map_sha256": row["canonical_atom_map"]["sha256"]}
                    source["pose_identity_sha256"] = canonical_sha256(source)
                    poses.append(source);pool["original_source_orders"][algorithm].append(source["pose_id"])
            pools.append(pool)
    if len(pools) != 108 or terminal > 216:
        raise ValueError("Union108-pool/216-source-slot denominator changed")
    schedule = json.loads(_pin(root, plan["source_schedule"]).read_text())
    if (len(schedule["pools"]) != 108 or len({pool["pool_id"] for pool in schedule["pools"]}) != 108
            or {pool["pool_id"] for pool in pools} != {pool["pool_id"] for pool in schedule["pools"]}):
        raise ValueError("Actual union pools differ from the explicit pre-score108-pool schedule")
    target = output / "source_pool.json"
    if target.exists():
        previous = json.loads(target.read_text())
        if previous["source_inputs_sha256"] != plan["source_inputs"]["sha256"]:
            raise ValueError("Union source inputs changed during checked resume")
        poses = merge_pose_index(previous["poses"], poses)
    else:
        poses = merge_pose_index([], poses)
    result = {"schema_version": 1, "configuration_sha256": sha256(actual_plan_path),
              "method_lock_sha256": plan["method_lock"]["sha256"], "source_inputs_sha256": plan["source_inputs"]["sha256"],
              "eligible_seed_pools": 108, "eligible_source_method_slots": 216, "terminal_source_method_slots": terminal,
              "all216_source_slots_terminal_and_verified": terminal == 216 and not pending,
              "pending_source_slots": pending, "retained_pose_count": len(poses), "poses": poses, "pools": pools,
              "deduplicated_poses": 0, "native_coordinates_used": False, "atomic_accuracy_computed": False}
    write_json(target, result)
    return result


def restore_carb_pose_aliases(prepared_ligand, raw_poses, export_map, destination):
    """Restore every heavy/H atom's generic identity while preserving raw xyz.

    Frozen read_carb_poses validates the canonical input/output atom universe
    and typing first. Coordinates are then copied literally, never fitted.
    """
    prepared_ligand, raw_poses, destination = map(Path, (prepared_ligand, raw_poses, destination))
    if destination.exists():
        raise FileExistsError("Generic restored poses require a fresh destination")
    read_carb_poses(raw_poses, export_map)
    source = {}
    for line in prepared_ligand.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = _atom(line)
            if atom["serial"] in source:
                raise ValueError("Duplicate prepared atom serial")
            source[atom["serial"]] = (atom, line)
    mapped = {atom["serial"]: atom for atom in export_map["atoms"]}
    if source.keys() != mapped.keys():
        raise ValueError("Canonical export and common source have different atom universes")
    for serial, (atom, _) in source.items():
        entry = mapped[serial]
        if (atom["name"] != entry["generated_name"] or atom["type"] != entry["atom_type"]
                or atom["charge"] != entry["charge"]):
            raise ValueError("Canonical export map no longer restores exact common atom identity")
    output, copied_atoms = [], 0
    for line in raw_poses.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = _atom(line);original = source[atom["serial"]][1]
            restored = line[:12] + original[12:30] + line[30:]
            generic = _atom(restored)
            if any(generic[key] != atom[key] for key in ("serial", "coordinates", "type", "charge")):
                raise AssertionError("Canonical alias restoration changed pose geometry or typing")
            line = restored;copied_atoms += 1
        output.append(line)
    write_text_atomic(destination, "\n".join(output) + "\n")
    return {"source_raw_pose_sha256": sha256(raw_poses), "source_prepared_ligand_sha256": sha256(prepared_ligand),
            "restored_generic_pose_sha256": sha256(destination), "all_heavy_and_hydrogen_coordinates_preserved": True,
            "copied_atom_records": copied_atoms, "native_coordinates_accessed": False}


def export_source_pose(source, row, root, destination):
    """Export one source model to the unchanged original VC fixed-pose helper."""
    root, destination = Path(root), Path(destination)
    prepared = root / row["ligand"]["path"]
    raw = root / source["source_pose_file"]
    mapping = normalize_atom_map(json.loads((root / row["canonical_atom_map"]["path"]).read_text()))
    if (sha256(prepared) != source["source_ligand_sha256"] or sha256(raw) != source["source_pose_file_sha256"]
            or sha256(root / row["canonical_atom_map"]["path"]) != source["source_canonical_atom_map_sha256"]
            or sha256(root / row["receptor"]["path"]) != source["source_receptor_sha256"]
            or hashlib.sha256(model_blocks(raw.read_bytes())[source["model_rank"]]).hexdigest() != source["source_model_sha256"]):
        raise ValueError("Union source pose/input changed after indexing")
    with tempfile.TemporaryDirectory(prefix="paired-union-alias-") as temporary:
        if source["source_algorithm"] == "vina_carb_1_0":
            export_map_path = root / source["source_export_map"]["path"]
            if sha256(export_map_path) != source["source_export_map"]["sha256"]:
                raise ValueError("Original VC full-search identity map changed")
            generic = Path(temporary) / "generic_poses.pdbqt"
            alias = restore_carb_pose_aliases(prepared, raw, json.loads(export_map_path.read_text()), generic)
            receipt = export_vina_carb_pose(prepared, mapping, generic, source["model_rank"], destination)
            receipt["canonical_alias_roundtrip"] = alias
        elif source["source_algorithm"] == "modern_vina_1_2_7":
            receipt = export_vina_carb_pose(prepared, mapping, raw, source["model_rank"], destination)
        else:
            raise ValueError("Unknown union source algorithm")
    receipt["source_raw_pose_sha256"] = source["source_pose_file_sha256"]
    receipt["source_raw_model_sha256"] = source["source_model_sha256"]
    write_json(destination.with_suffix(".pdbqt.map.json"), receipt)
    return receipt


def rank_union_pools(index, records):
    """Require a complete exact union; every unsupported pose fails its pool."""
    rows = []
    normalizers = {}
    for record in records.values():
        if record["status"] == "ok":
            normalizers.setdefault(record["source"]["complex_id"], []).append(record["score"]["torsion_normalization_scale"])
    conflicts = {case for case, values in normalizers.items() if max(values) - min(values) > 1e-12}
    for pool in index["pools"]:
        selected = [source for source in index["poses"] if source["pool_id"] == pool["pool_id"]]
        result = {**pool, "union_pose_ids": [source["pose_id"] for source in selected],
                  "coefficient0_order": [], "coefficient1_order": [], "atomic_accuracy_computed": False}
        if any(status == "pending" for status in pool["source_statuses"].values()):
            status = "pending_source"
        elif any(status != "ok" for status in pool["source_statuses"].values()):
            status = "source_workflow_failure_whole_pair_required"
        elif pool["case_id"] in conflicts:
            status = "normalization_conflict"
        elif not selected:
            status = "source_workflow_failure_no_poses"
        elif any(source["pose_id"] not in records for source in selected):
            status = "pending_scoring"
        elif any(records[source["pose_id"]]["status"] != "ok" for source in selected):
            status = "pose_scoring_failure_whole_union_required"
        else:
            status = "completed"
            for coefficient, field in ((0, "base_shared0"), (1, "combined_shared0")):
                ordered = sorted(selected, key=lambda source: (records[source["pose_id"]]["score"][field],
                                  source["source_algorithm_rank"], source["candidate_class_rank"], source["model_rank"]))
                result[f"coefficient{coefficient}_order"] = [source["pose_id"] for source in ordered]
        result["status"] = status;rows.append(result)
    return rows, sorted(conflicts)


def score_available_union(root=ROOT, *, plan_path, receipt_path, installation, output_root=None, max_poses=None):
    """One sequential fixed-helper worker, allowed only after the root freeze."""
    root = Path(root).resolve()
    plan, actual_plan = validate_union_plan(root, plan_path=plan_path, receipt_path=receipt_path, require_freeze=True)
    if max_poses is not None and (isinstance(max_poses, bool) or not isinstance(max_poses, int) or max_poses < 1):
        raise ValueError("Execution shard limit must be positive; scientific pool membership remains unchanged")
    output = Path(output_root) if output_root else root / DEFAULT_OUTPUT
    index = index_union_pool(root, plan_path=actual_plan, output_root=output)
    experiment = json.loads(_pin(root, plan["source_experiment"]).read_text())
    rows = {row["condition_id"]: row for row in json.loads(_pin(root, plan["source_inputs"]).read_text())["conditions"]}
    records, manifest = _old_records(output, index, recover_orphans=True)
    write_json(output / "scores_manifest.json", manifest)
    source_archive = output / "union_adapter_source.py"
    if source_archive.exists() and sha256(source_archive) != sha256(Path(__file__)):
        raise ValueError("Union source changed after native scoring began")
    if not source_archive.exists():
        write_text_atomic(source_archive, Path(__file__).read_text())
    processed = 0
    for source in index["poses"]:
        identity = source["pose_id"]
        if identity in records:
            continue
        if max_poses is not None and processed >= max_poses:
            break
        base = {"schema_version": 1, "source": source, "method_lock_sha256": index["method_lock_sha256"],
                "coordinate_optimization": False, "native_coordinates_used": False}
        row = rows[source["condition_id"]]
        with tempfile.TemporaryDirectory(prefix="paired-union-score-") as temporary:
            exported = Path(temporary) / "ligand.pdbqt"
            try:
                export = export_source_pose(source, row, root, exported)
                start = time.monotonic();before = resource.getrusage(resource.RUSAGE_CHILDREN)
                try:
                    score = _invoke_fixed_pose(root / row["receptor"]["path"], exported, export, installation, root=root)
                finally:
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    base["native_wall_seconds"] = time.monotonic() - start
                    base["native_cpu_seconds"] = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                if score["executable_sha256"] != experiment["native_preflight_binary_sha256"]:
                    raise NativeScoringFailure("Union fixed helper differs from the observed frozen native scorer", score)
                record = {**base, "status": "ok", "score": score, "exported_ligand_sha256": sha256(exported),
                          "export_map_sha256": sha256(exported.with_suffix(".pdbqt.map.json")),
                          "expected_glycosidic_links": len(export["expected_glycosidic_links"])}
                try:
                    verify_score_record(record, source, index["method_lock_sha256"])
                except (ValueError, TypeError, KeyError) as error:
                    raise NativeScoringFailure(str(error), score) from error
            except NativeScoringFailure as error:
                record = {**base, "status": "pose_scoring_failed", "error": str(error), "native_attempt": error.attempt,
                          "exported_ligand_sha256": sha256(exported), "export_map_sha256": sha256(exported.with_suffix(".pdbqt.map.json")),
                          "expected_glycosidic_links": len(export["expected_glycosidic_links"])}
            except (ValueError, RuntimeError, OSError) as error:
                record = {**base, "status": "pose_scoring_failed", "error": str(error), "error_type": type(error).__name__}
        verify_score_record(record, source, index["method_lock_sha256"])
        target = _record_path(output, identity)
        if target.exists():
            raise FileExistsError("Refusing to overwrite an earlier union score record")
        write_json(target, record);records[identity] = record;manifest["records"][identity] = sha256(target)
        write_json(output / "scores_manifest.json", manifest);processed += 1
        if processed % 50 == 0:
            print(f"Union helper recorded {processed} new poses; total {len(records)}/{len(index['poses'])}", flush=True)
    pools, conflicts = rank_union_pools(index, records)
    result = {"schema_version": 1, "configuration_sha256": index["configuration_sha256"],
              "source_pool_sha256": sha256(output / "source_pool.json"), "score_manifest_sha256": sha256(output / "scores_manifest.json"),
              "eligible_seed_pools": 108, "eligible_rescore_method_slots": 216,
              "indexed_pose_count": len(index["poses"]), "recorded_pose_count": len(records),
              "successful_pose_scores": sum(record["status"] == "ok" for record in records.values()),
              "failed_pose_scores": sum(record["status"] != "ok" for record in records.values()),
              "new_records_this_pass": processed, "pools": pools, "normalization_conflict_cases": conflicts,
              "all216_source_slots_terminal_and_verified": index["all216_source_slots_terminal_and_verified"],
              "ready_for_separate_accuracy_evaluation": index["all216_source_slots_terminal_and_verified"] and len(records) == len(index["poses"]),
              "atomic_accuracy_computed": False}
    write_json(output / "scoring_progress.json", result)
    return result


def verify_union_scores(root=ROOT, *, plan_path, receipt_path, output_root=None, require_complete=True):
    """Portable raw-output/alias replay; never invoke a native executable."""
    root = Path(root).resolve()
    plan, actual_plan = validate_union_plan(root, plan_path=plan_path, receipt_path=receipt_path, require_freeze=True)
    output = Path(output_root) if output_root else root / DEFAULT_OUTPUT
    index = index_union_pool(root, plan_path=actual_plan, output_root=output)
    records, _ = _old_records(output, index)
    if require_complete and (not index["all216_source_slots_terminal_and_verified"] or len(records) != len(index["poses"])):
        raise ValueError("All216 source slots and every retained union pose must be terminal before metric aggregation")
    rows = {row["condition_id"]: row for row in json.loads(_pin(root, plan["source_inputs"]).read_text())["conditions"]}
    reconstructed = 0
    for source in index["poses"]:
        record = records.get(source["pose_id"])
        if record is None or "exported_ligand_sha256" not in record:
            continue
        with tempfile.TemporaryDirectory(prefix="paired-union-verify-") as temporary:
            destination = Path(temporary) / "ligand.pdbqt"
            export = export_source_pose(source, rows[source["condition_id"]], root, destination)
            if (sha256(destination) != record["exported_ligand_sha256"]
                    or sha256(destination.with_suffix(".pdbqt.map.json")) != record["export_map_sha256"]
                    or len(export["expected_glycosidic_links"]) != record["expected_glycosidic_links"]):
                raise ValueError("Union score cannot reconstruct unchanged source pose geometry/identity")
            if record["status"] == "ok":
                log = Path(temporary) / "VC_log.txt";log.write_text(record["score"]["recognition_log"])
                verify_chi_recognition(log, expected_chi_links(export))
        reconstructed += 1
    pools, conflicts = rank_union_pools(index, records)
    return {"status": "verified", "native_executables_used": False, "eligible_seed_pools": 108,
            "eligible_source_method_slots": 216, "eligible_rescore_method_slots": 216,
            "terminal_source_method_slots": index["terminal_source_method_slots"], "retained_pose_count": len(index["poses"]),
            "recorded_pose_count": len(records), "source_exports_reconstructed": reconstructed,
            "all216_source_slots_terminal_and_verified": index["all216_source_slots_terminal_and_verified"],
            "pools": pools, "normalization_conflict_cases": conflicts, "atomic_accuracy_computed": False}
