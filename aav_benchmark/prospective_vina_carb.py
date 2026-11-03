"""Prospective original Vina-Carb full search on exactly paired prepared inputs.

No native call is permitted without a separately hashed root freeze. This
module never reads experimental ligand coordinates or computes pose accuracy.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time

from .additional_method import (
    export_vina_carb_ligand, expected_chi_links, read_carb_poses,
    validate_installation, verify_chi_recognition,
)
from .common_pose_rescoring import NativeScoringFailure, _invoke_fixed_pose
from .provenance import sha256, write_json, write_text_atomic

ROOT = Path(__file__).resolve().parents[1]
CASES = ["2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U"]
SEEDS = [104729, 130363, 155921]
METHODS = ["fpocket", "p2rank", "oracle"]
PROFILES = ["strict", "repair"]
PARAMETERS = {"cpu": 2, "exhaustiveness": 16, "num_modes": 20,
              "energy_range_kcal_mol": 10, "chi_coeff": 1, "chi_cutoff": 2}
HELPER_SOURCES = {"aav_benchmark/additional_method.py", "aav_benchmark/common_pose_rescoring.py",
                  "aav_benchmark/provenance.py", "results/additional_method_reproduction/build/fixed_pose_score.cpp"}


def _checked(root, entry):
    path = Path(entry["path"])
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("Prepared inputs must use repository-relative paths")
    target = Path(root) / path
    if target.is_symlink() or sha256(target) != entry["sha256"]:
        raise ValueError("A prospective prepared input changed")
    return target


def nominal_cost(secondary_top5=False, receptor_profiles=2):
    if isinstance(receptor_profiles, bool) or not isinstance(receptor_profiles, int) or receptor_profiles < 1:
        raise ValueError("Receptor-profile count must be a positive integer")
    boxes = 11 if secondary_top5 else 3
    jobs = len(CASES) * boxes * len(SEEDS) * receptor_profiles
    return {"maximum_VinaCarb_seed_box_jobs_before_preparation_failures": jobs,
            "maximum_jobs_per_receptor_profile": jobs // receptor_profiles,
            "paired_modern_Vina_additional_jobs_if_all_inputs_prepared": jobs,
            "per_search_cpu_allocation": 2,
            "wall_or_CPU_time_estimate_claimed": False,
            "secondary_top5_enabled": bool(secondary_top5)}


def validate_frozen_plan(root, plan_path, receipt_path):
    root, plan_path, receipt_path = Path(root).resolve(), Path(plan_path), Path(receipt_path)
    plan = json.loads(plan_path.read_text());receipt = json.loads(receipt_path.read_text())
    if (plan.get("status") != "frozen_before_native_full_search"
            or plan.get("frozen_before_native_execution") is not True
            or receipt.get("frozen_before_native_execution") is not True
            or receipt.get("new_method_outcomes_consulted") is not False
            or plan.get("new_method_outcomes_consulted") is not False
            or receipt.get("configuration_sha256") != sha256(plan_path)):
        raise ValueError("A separate root freeze is required before any native Vina-Carb call")
    if (plan["cases"] != CASES or plan["seeds"] != SEEDS or plan["parameters"] != PARAMETERS
            or plan["receptor_profiles"] != PROFILES or plan["methods"] != METHODS
            or plan["method_lock"]["path"] != "benchmarks/additional_method.lock.json"):
        raise ValueError("Prospective method, cohort, paired profiles or effort differs from the reviewed plan")
    if plan["adapter_source"]["path"] != "aav_benchmark/prospective_vina_carb.py":
        raise ValueError("Unexpected prospective native adapter")
    if plan["adapter_source"]["sha256"] != sha256(Path(__file__)):
        raise ValueError("Executing adapter differs from the independently frozen source")
    for key in ("method_lock", "initialization_policy", "initialization_algorithm", "prepared_inputs", "adapter_source", "paired_modern_vina_policy", "native_build_receipt", "native_preflight_receipt"):
        _checked(root, plan[key])
        if receipt.get(key + "_sha256") != plan[key]["sha256"]:
            raise ValueError("Root receipt does not bind all method/initialization/input artifacts")
    if set(plan["helper_sources"]) != HELPER_SOURCES or receipt.get("helper_sources") != plan["helper_sources"]:
        raise ValueError("Every imported native helper and fixed caller needs a root-frozen source pin")
    for path, digest in plan["helper_sources"].items():
        _checked(root, {"path": path, "sha256": digest})
    lock = json.loads(_checked(root, plan["method_lock"]).read_text())
    if lock["fixed_pose_rescorer"]["source_sha256"] != plan["helper_sources"]["results/additional_method_reproduction/build/fixed_pose_score.cpp"]:
        raise ValueError("Frozen method lock and fixed CHI caller disagree")
    native = json.loads(_checked(root, plan["native_build_receipt"]).read_text())
    preflight = json.loads(_checked(root, plan["native_preflight_receipt"]).read_text())
    if (native["executable_sha256"] != plan["native_binary_sha256"]
            or preflight["executable_sha256"] != plan["native_preflight_binary_sha256"]
            or preflight["caller_source_sha256"] != lock["fixed_pose_rescorer"]["source_sha256"]
            or preflight["method_lock_sha256"] != plan["method_lock"]["sha256"]):
        raise ValueError("Frozen native executable/caller receipts disagree")
    if plan.get("reuse_byte_identical_profile_jobs") is not False:
        raise ValueError("This reviewed adapter recomputes profile jobs; reuse requires a separate declared implementation")
    manifest = json.loads(_checked(root, plan["prepared_inputs"]).read_text())
    rows = manifest["conditions"]
    expected = {(case, profile, method) for case in CASES for profile in PROFILES for method in METHODS}
    observed = {(row["case_id"], row["receptor_profile"], row["method"]) for row in rows}
    if observed != expected:
        raise ValueError("Every six-case primary/repair/predictor condition, including failures, must remain eligible")
    identities = [row["condition_id"] for row in rows]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate prospective condition identity")
    for identity in expected:
        selected = [row for row in rows if (row["case_id"], row["receptor_profile"], row["method"]) == identity]
        limit = 5 if plan["secondary_top5_enabled"] and identity[-1] != "oracle" else 1
        prepared = [row for row in selected if row["status"] == "prepared"]
        ranks = [row["candidate_class_rank"] for row in prepared]
        if len(selected) > limit or len(set(ranks)) != len(ranks) or (prepared and len(prepared) != len(selected)):
            raise ValueError("A paired condition duplicated, mixed, or exceeded its candidate budget")
    for row in rows:
        source_condition = json.loads(_checked(root, row["source_condition"]).read_text())
        if source_condition["case_id"] != row["case_id"] or source_condition["method"] != row["method"]:
            raise ValueError("Shared source-condition identity changed")
        if set(row["vina_carb_outputs"]) != {str(seed) for seed in SEEDS}:
            raise ValueError("Every eligible seed needs an explicit future VC output path")
        for path in [row["modern_runs"], *row["vina_carb_outputs"].values()]:
            if Path(path).is_absolute() or ".." in Path(path).parts:
                raise ValueError("Future native outputs must be repository-relative")
        if row["status"] not in {"prepared", "preparation_failed", "predictor_unavailable"}:
            raise ValueError("Unknown prospective preparation status")
        if row["status"] != "prepared":
            if not row.get("failure_reason"):
                raise ValueError("An eligible failed condition needs its explicit mechanism")
            continue
        source_jobs = [job for job in source_condition["jobs"] if job["candidate_class_rank"] == row["candidate_class_rank"]]
        if (len(source_jobs) != 1 or source_jobs[0]["status"] != "prepared"
                or source_jobs[0]["preparation_sha256"] != row["preparation"]["sha256"]
                or source_jobs[0]["candidate_original_rank"] != row["candidate_original_rank"]):
            raise ValueError("Paired candidate rank or shared preparation binding changed")
        for key in ("center", "size"):
            values = row["box"][key]
            if len(values) != 3 or not all(isinstance(v, (float, int)) and not isinstance(v, bool) and math.isfinite(v) for v in values):
                raise ValueError("Invalid paired docking box")
        if min(row["box"]["size"]) <= 0:
            raise ValueError("Invalid paired docking box size")
        for key in ("ligand", "canonical_atom_map", "initializer_receipt", "initializer_sdf", "receptor", "preparation"):
            _checked(root, row[key])
        initializer = json.loads(_checked(root, row["initializer_receipt"]).read_text())
        if (initializer["case_id"] != row["case_id"] or initializer["profile_id"] != plan["initialization_profile_id"]
                or initializer["policy_sha256"] != plan["initialization_policy"]["sha256"]
                or initializer["algorithm_source_sha256"] != plan["initialization_algorithm"]["sha256"]
                or initializer["search"]["status"] != "selected"
                or initializer["serialized_sdf_validation"]["valid"] is not True
                or initializer["selected_sdf_sha256"] != row["initializer_sdf"]["sha256"]
                or initializer["no_docking_or_accuracy_evaluation"] is not True):
            raise ValueError("Paired ligand is not bound to a successful independently validated chemical-ideal initializer")
        paired = json.loads(_checked(root, row["preparation"]).read_text())
        directory = Path(row["modern_vina_prepared_directory"])
        if (row["ligand"]["path"] != (directory / paired["ligand_file"]).as_posix()
                or row["receptor"]["path"] != (directory / "neighborhood.pdbqt").as_posix()
                or row["preparation"]["path"] != (directory / "preparation.json").as_posix()
                or paired["box_center"] != row["box"]["center"] or paired["box_size"] != row["box"]["size"]
                or paired["case_id"] != row["case_id"]
                or paired["configuration_sha256"] != plan["paired_modern_vina_policy"]["sha256"]
                or paired["native_coordinates_in_runner_reference"] is not False):
            raise ValueError("Vina-Carb and modern Vina must use the same prepared ligand, receptor, configuration and box")
        if (row["initializer_profile_id"] != plan["initialization_profile_id"]
                or row["candidate_class_rank"] < 1
                or row["candidate_class_rank"] > (5 if plan["secondary_top5_enabled"] and row["method"] != "oracle" else 1)):
            raise ValueError("Initializer or candidate budget differs from the frozen paired protocol")
    return plan, rows


def normalize_atom_map(value):
    """Convert only independent chemical identities to the frozen exporter API."""
    atoms = value.get("atoms", value.get("canonical_atom_map", value.get("atom_map"))) if isinstance(value, dict) else value
    if not isinstance(atoms, list) or not atoms:
        raise ValueError("An independent canonical atom map is required")
    result = []
    for atom in atoms:
        if "source_key" in atom:
            key = atom["source_key"]
        else:
            number = int(atom["branch_num"])
            if number < 1:
                raise ValueError("Independent chemical branch numbers must be positive")
            key = ["B", number, "", atom["ccd_id"], atom["atom_name"]]
        result.append({"name": atom["name"], "source_key": key})
    return result


def render_configuration(row, seed):
    if seed not in SEEDS:
        raise ValueError("Unscheduled prospective seed")
    values = {"receptor": "receptor.pdbqt", "ligand": "ligand.pdbqt",
              **{f"center_{axis}": row["box"]["center"][i] for i, axis in enumerate("xyz")},
              **{f"size_{axis}": row["box"]["size"][i] for i, axis in enumerate("xyz")},
              "cpu": PARAMETERS["cpu"], "seed": seed, "exhaustiveness": PARAMETERS["exhaustiveness"],
              "num_modes": PARAMETERS["num_modes"], "energy_range": PARAMETERS["energy_range_kcal_mol"],
              "chi_coeff": PARAMETERS["chi_coeff"], "chi_cutoff": PARAMETERS["chi_cutoff"],
              "out": "poses.pdbqt", "log": "search.log"}
    return "".join(f"{key} = {value}\n" for key, value in values.items())


def build_prepared_manifest(*, root=ROOT, prepared_root, initialization_root, modern_config, output, preparation_index=None):
    """Index the ordinary shared prepared boxes; no native or metric calls.

    All36 case/profile/site conditions must have completed preparation records.
    Selecting candidate class rank1 is the declared primary budget, independent
    of any downstream docking result.
    """
    root = Path(root).resolve();base = Path(prepared_root).resolve();initializers = Path(initialization_root).resolve()
    def pin(path):
        path = Path(path).resolve()
        return {"path": path.relative_to(root).as_posix(), "sha256": sha256(path)}
    preparation_index = Path(preparation_index) if preparation_index else base / "prepared_conditions.json"
    source_index = json.loads(preparation_index.read_text())
    if (source_index["status"] not in {"complete", "completed"} or source_index["condition_count"] != 36
            or source_index["expected_condition_count"] != 36
            or source_index["configuration_sha256"] != sha256(modern_config)
            or source_index["no_docking_started"] is not True or source_index["accuracy_evaluated"] is not False):
        raise ValueError("The complete source-only36-condition preparation index is required before native input freeze")
    indexed = {row["condition_path"]: row["condition_sha256"] for row in source_index["conditions"]}
    conditions = []
    for profile in PROFILES:
        for case in CASES:
            for method in METHODS:
                directory = base / profile / case / method / "prepared"
                condition_path = directory / "condition.json"
                if indexed.get(condition_path.relative_to(root).as_posix()) != sha256(condition_path):
                    raise ValueError("Prepared condition differs from the complete source index")
                condition = json.loads(condition_path.read_text())
                if condition["case_id"] != case or condition["method"] != method or condition.get("status") in {"preparing", "pending"}:
                    raise ValueError("Shared source preparation is missing, pending or mismatched")
                row = {"condition_id": f"{profile}/{case}/{method}/box01", "case_id": case, "receptor_profile": profile,
                       "method": method, "source_condition": pin(condition_path), "eligible_workflow_condition": True}
                selected = [job for job in condition["jobs"] if job["candidate_class_rank"] == 1]
                if len(selected) > 1:
                    raise ValueError("A shared source condition has duplicate top-ranked boxes")
                if not selected:
                    row.update(status="predictor_unavailable" if condition["status"] not in {"preparation_failed", "source_preparation_failed"} else "preparation_failed",
                               failure_reason=condition.get("reason", condition["status"]))
                elif selected[0]["status"] != "prepared":
                    row.update(status="preparation_failed", failure_reason=selected[0].get("reason", selected[0]["status"]))
                else:
                    job = selected[0];box = directory / job["id"];preparation = json.loads((box / "preparation.json").read_text())
                    if sha256(box / "preparation.json") != job["preparation_sha256"]:
                        raise ValueError("Shared source job/preparation hashes differ")
                    atom_map = box / "ligand_atom_map.json"
                    if not atom_map.exists():
                        atom_map = box / "reference.json"
                    source_initializer = initializers / case
                    report = json.loads((source_initializer / "report.json").read_text())
                    row.update(status="prepared", candidate_class_rank=job["candidate_class_rank"],
                               candidate_original_rank=job["candidate_original_rank"], modern_vina_prepared_directory=box.relative_to(root).as_posix(),
                               box={"center": preparation["box_center"], "size": preparation["box_size"]},
                               ligand=pin(box / "ligand.pdbqt"), canonical_atom_map=pin(atom_map),
                               receptor=pin(box / "neighborhood.pdbqt"), preparation=pin(box / "preparation.json"),
                               initializer_profile_id=report["profile_id"], initializer_receipt=pin(source_initializer / "report.json"),
                               initializer_sdf=pin(source_initializer / "ligand.sdf"))
                row["modern_runs"] = (directory.parent / "vina/box01/runs.json").relative_to(root).as_posix()
                row["vina_carb_outputs"] = {str(seed): (directory.parent / "vina_carb/box01" / f"s{seed}").relative_to(root).as_posix() for seed in SEEDS}
                conditions.append(row)
    result = {"schema_version": 1, "primary_candidate_budget": 1, "secondary_top5_enabled": False,
              "source_preparation_index": pin(preparation_index),
              "modern_vina_configuration": pin(modern_config), "eligible_conditions": 36,
              "eligible_seed_slots_per_engine": 36 * len(SEEDS), "conditions": conditions,
              "prepared_condition_count": sum(row["status"] == "prepared" for row in conditions),
              "actual_native_search_slots_per_engine": sum(row["status"] == "prepared" for row in conditions) * len(SEEDS),
              "source_only_index": True, "docking_or_accuracy_outcomes_consulted": False}
    write_json(output, result)
    return result


def execute_one(*, root=ROOT, plan_path, receipt_path, condition_id, seed, installation, output):
    """Run one explicitly scheduled job only after an independent root freeze."""
    root = Path(root).resolve()
    plan, rows = validate_frozen_plan(root, plan_path, receipt_path)
    row = next((row for row in rows if row["condition_id"] == condition_id), None)
    if row is None or seed not in SEEDS:
        raise ValueError("Condition/seed is outside the frozen full-search schedule")
    output = Path(output).resolve()
    if output != (root / row["vina_carb_outputs"][str(seed)]).resolve():
        raise ValueError("Native output directory differs from the frozen schedule")
    if output.exists():
        raise FileExistsError("Refusing to overwrite a prospective search attempt")
    output.mkdir(parents=True)
    record = {"schema_version": 1, "condition": row, "seed": seed, "configuration_sha256": sha256(plan_path),
              "freeze_receipt_sha256": sha256(receipt_path), "driver_source_sha256": sha256(Path(__file__)),
              "atomic_accuracy_computed": False, "experimental_ligand_coordinates_read": False,
              "status": "prepared_input_failure" if row["status"] != "prepared" else "preparing",
              "full_search_started": False, "cpu_seconds": 0.0, "wall_seconds": 0.0}
    if row["status"] != "prepared":
        record["error"] = row["failure_reason"];write_json(output / "execution.json", record);return record
    executable = validate_installation(installation, root=root)
    runtime_scorer = json.loads((Path(installation) / "fixed_pose_scorer.json").read_text())
    if (sha256(executable) != plan["native_binary_sha256"]
            or runtime_scorer["executable_sha256"] != plan["native_preflight_binary_sha256"]
            or runtime_scorer["method_lock_sha256"] != plan["method_lock"]["sha256"]
            or runtime_scorer["caller_source_sha256"] != plan["helper_sources"]["results/additional_method_reproduction/build/fixed_pose_score.cpp"]):
        raise ValueError("Installed native binaries differ from the frozen observed build")
    try:
        source = _checked(root, row["ligand"])
        mapping = normalize_atom_map(json.loads(_checked(root, row["canonical_atom_map"]).read_text()))
        export = export_vina_carb_ligand(source, mapping, output / "ligand.pdbqt")
        residue_count = len({tuple(atom["source_key"][:4]) for atom in mapping})
        if residue_count < 2 or len(export["expected_glycosidic_links"]) != residue_count - 1:
            raise ValueError("Every full-parent interresidue linkage must be exposed; zero-link substitution is forbidden")
        (output / "receptor.pdbqt").write_bytes(_checked(root, row["receptor"]).read_bytes())
        write_text_atomic(output / "configuration.txt", render_configuration(row, seed))
        before = resource.getrusage(resource.RUSAGE_CHILDREN);start = time.monotonic()
        try:
            preflight = _invoke_fixed_pose(output / "receptor.pdbqt", output / "ligand.pdbqt", export, installation, root=root)
            if (preflight["glycosidic_links"] != len(export["expected_glycosidic_links"])
                    or preflight["fully_supported_links"] != len(export["expected_glycosidic_links"])):
                raise NativeScoringFailure("Not every source glycosidic link is fully supported", preflight)
            record["preflight"] = preflight
        finally:
            after = resource.getrusage(resource.RUSAGE_CHILDREN)
            record["preflight_wall_seconds"] = time.monotonic() - start
            record["preflight_cpu_seconds"] = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
        record.update(status="running", full_search_started=True, executable_sha256=sha256(executable),
                      exported_ligand_sha256=sha256(output / "ligand.pdbqt"), receptor_sha256=sha256(output / "receptor.pdbqt"),
                      command=[str(executable), "--config", "configuration.txt"], parameters=PARAMETERS)
        write_json(output / "execution.json", record)
        before = resource.getrusage(resource.RUSAGE_CHILDREN);start = time.monotonic()
        with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
            completed = subprocess.run(record["command"], cwd=output,
                                       env={**os.environ, "VINA_CARB": str(Path(installation).resolve() / "VC_1_0")}, stdout=stdout, stderr=stderr)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        record.update(exit_code=completed.returncode, wall_seconds=time.monotonic() - start,
                      cpu_seconds=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime)
        if completed.returncode:
            raise RuntimeError("Original Vina-Carb full search exited unsuccessfully")
        verify_chi_recognition(output / "VC_log.txt", expected_chi_links(export))
        poses = read_carb_poses(output / "poses.pdbqt", export)
        if not 1 <= len(poses) <= PARAMETERS["num_modes"]:
            raise ValueError("Returned pose count differs from the frozen output cap")
        record.update(status="ok", returned_poses=len(poses), generated_names=poses[0]["names"],
                      score_role="Original Vina-Carb full search including CHI; original shared best-mode intra convention")
    except NativeScoringFailure as error:
        record.update(status="unsupported_preflight", error=str(error), native_attempt=error.attempt)
    except (ValueError, RuntimeError, OSError) as error:
        record.update(status="failed", error=str(error), error_type=type(error).__name__)
    record["artifacts"] = {path.name: sha256(path) for path in output.iterdir() if path.is_file() and path.name != "execution.json"}
    record["total_native_wall_seconds"] = record["wall_seconds"] + record.get("preflight_wall_seconds", 0)
    record["total_native_cpu_seconds"] = record["cpu_seconds"] + record.get("preflight_cpu_seconds", 0)
    write_json(output / "execution.json", record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("plan", "receipt", "installation", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--condition-id", required=True);parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args(argv)
    result = execute_one(plan_path=args.plan, receipt_path=args.receipt, condition_id=args.condition_id,
                         seed=args.seed, installation=args.installation, output=args.output)
    print(json.dumps(result, indent=2));return 0 if result["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
