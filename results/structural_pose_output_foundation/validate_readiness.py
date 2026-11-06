#!/usr/bin/env python3
"""Root-authorized metadata/opaque-byte readiness; no scientific readers."""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import socket
import subprocess

PACKAGER_SHA256 = "400cc0cbc4fb3700a93403ae642c1f9dc4b5c0f468892db8f55295278d71b12d"
PLAN_SHA256 = "961e5537b1a48223a35cad10f61720db2c0f08772326a1b6f8892d90874bc600"
DISTRIBUTION = "benchmarks/method_protocol_public_distribution.json"
DISTRIBUTION_SHA256 = "bed3115a6ab2685aa3c9cd790ff5080e931ee813c37c720188b4e44901570beb"
MAIN = "results/structural_pose_v3"
ASSEMBLY = "results/structural_assembly_pose_v3"
HADD = "results/haddock_paired_cohort"
UNION = MAIN + "/common_pose_union"
VC_PLAN = "benchmarks/structural_test_vina_carb.json"
MODERN = "benchmarks/structural_pose_v3_experiment.json"
VC_FREEZE = MAIN + "/vina_carb_freeze_receipt.json"
MAIN_FREEZE = MAIN + "/freeze_receipt.json"
INPUTS = MAIN + "/vina_carb_prepared_inputs.json"
ASSEMBLY_CONFIG = "benchmarks/structural_assembly_pose_v3_experiment.json"
ASSEMBLY_INPUTS = ASSEMBLY + "/prepared_conditions.json"
ASSEMBLY_FREEZE = ASSEMBLY + "/freeze_receipt.json"
ASSEMBLY_DRIVER_FREEZE = ASSEMBLY + "/execution_driver_freeze.json"
HADD_PLAN = HADD + "/prepared_plan.json"
HADD_FREEZE = HADD + "/execution_freeze.json"
RECOVERY_PLAN = HADD + "/orchestration_recovery/recovery_plan.json"
RECOVERY_FREEZE = HADD + "/orchestration_recovery/recovery_freeze.json"
UNION_PLAN = "benchmarks/structural_pose_v3_common_pool.json"
UNION_FREEZE = UNION + "/freeze_receipt.json"
SCHEDULE = UNION + "/pool_schedule.json"
VINA_LOCK = "benchmarks/vina.lock.json"
DRIVERS = {"main": MAIN + "/paired_execution_driver.py", "assembly": ASSEMBLY + "/execution_driver.py",
           "haddock": "aav_benchmark/haddock_cohort_resume.py", "union": "aav_benchmark/paired_common_pose_rescoring.py"}
STATIC_NAMES = [VC_PLAN, MODERN, VC_FREEZE, MAIN_FREEZE, INPUTS, ASSEMBLY_CONFIG, ASSEMBLY_INPUTS,
                ASSEMBLY_FREEZE, ASSEMBLY_DRIVER_FREEZE, HADD_PLAN, HADD_FREEZE, RECOVERY_PLAN,
                RECOVERY_FREEZE, UNION_PLAN, UNION_FREEZE, SCHEDULE, VINA_LOCK,
                *[HADD + "/orchestration_recovery/" + n for n in
                  ["adoption_receipt.json", "orphan_completion_proof.json", "recovery_preparation_receipt.json"]],
                *DRIVERS.values()]
WAIT_SOURCES = {
    "main": [DRIVERS["main"], "aav_benchmark/engine.py", "aav_benchmark/prospective_vina_carb.py", "aav_benchmark/common_pose_rescoring.py"],
    "assembly": [DRIVERS["assembly"], "aav_benchmark/engine.py"],
    "haddock": [DRIVERS["haddock"], "aav_benchmark/haddock_paired_cohort.py"],
    "union": [DRIVERS["union"], "aav_benchmark/common_pose_rescoring.py"],
}
WAIT_CONTRACTS = {
    "main": "sequential_engine_and_VC_subprocess_run_waits_then_synchronous_receipts",
    "assembly": "sequential_engine_subprocess_run_waits_then_synchronous_receipts",
    "haddock": "future_results_and_ProcessPoolExecutor_join_then_synchronous_final_progress",
    "union": "sequential_fixed_helper_subprocess_run_waits_then_synchronous_final_progress",
}
for _role in WAIT_SOURCES:
    WAIT_SOURCES[_role].append("aav_benchmark/provenance.py")
STATIC_NAMES += sorted({n for names in WAIT_SOURCES.values() for n in names} - set(STATIC_NAMES))
SEEDS = [104729, 130363, 155921]
CASES = ["2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U"]
ARMS = ["fpocket", "p2rank", "oracle"]
ALGORITHMS = ["modern_vina_1_2_7", "vina_carb_1_0"]
VC_TERMINAL = {"ok", "failed", "prepared_input_failure", "unsupported_preflight"}


class OpaqueFloat(str):
    """JSON number token retained verbatim; never converted to a float."""


def load_json(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate metadata JSON key")
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError("Nonfinite JSON token: " + value)
    return json.loads(raw, parse_float=OpaqueFloat, parse_constant=nonfinite, object_pairs_hook=unique)


def canonical_bytes(value):
    if isinstance(value, OpaqueFloat):
        require(re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", value) is not None,
                "Invalid opaque JSON number")
        return value.encode()
    if isinstance(value, dict):
        require(all(type(k) is str for k in value), "Canonical keys must be strings")
        return b"{" + b",".join(json.dumps(k).encode() + b":" + canonical_bytes(value[k]) for k in sorted(value)) + b"}"
    if isinstance(value, list):
        return b"[" + b",".join(canonical_bytes(v) for v in value) + b"]"
    require(type(value) in (str, int, bool, type(None)), "Use opaque JSON tokens, not interpreted floats")
    return json.dumps(value, allow_nan=False).encode()


def canonical_sha256(value):
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, *, minimum=0):
    require(type(value) is int and value >= minimum, "Exact integer identity/count required")
    return value


def sha(value):
    require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None, "Invalid SHA256")
    return value


def packager():
    path = Path(__file__).with_name("package_outputs.py")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == PACKAGER_SHA256, "Reviewed packager changed")
    spec = importlib.util.spec_from_file_location("reviewed_structural_packager", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


class Evidence:
    def __init__(self, root, pkg, static):
        self.root, self.pkg, self.static, self.files = root, pkg, static, {}

    def read(self, name, expected=None):
        self.pkg.portable_path(name); self.pkg._source_path_policy(name)
        raw = self.pkg._read(self.root, name)
        pin = {"path": name, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
        require(expected is None or pin["sha256"] == sha(expected), "Checksum changed: " + name)
        require(name not in self.files or self.files[name] == pin, "Evidence changed during validation: " + name)
        self.files[name] = pin
        return raw

    def json(self, name, expected=None):
        return load_json(self.read(name, expected))

    def pinned(self, pin, *, json_value=False):
        require(isinstance(pin, dict) and {"path", "sha256"} <= pin.keys(), "Path/SHA source pin required")
        raw = self.read(pin["path"], pin["sha256"])
        if "bytes" in pin:
            require(integer(pin["bytes"]) == len(raw), "Pinned byte count changed")
        return load_json(raw) if json_value else raw

    def source(self, name):
        require(name in self.static, "Unpinned source metadata: " + name)
        return self.pinned(self.static[name], json_value=True)

    def recheck(self):
        for pin in list(self.files.values()):
            self.pinned(pin)


def unique_rows(rows, field, expected=None):
    require(isinstance(rows, list), "Expected metadata row list")
    result = {}
    for row in rows:
        require(isinstance(row, dict) and field in row, "Missing metadata identity")
        key = row[field]
        require(type(key) in (str, int) and key not in result, "Duplicate/invalid metadata identity")
        result[key] = row
    require(expected is None or set(result) == set(expected), "Missing/extra metadata identities")
    return result


def model_hashes(raw):
    """MODEL identity syntax only; atom/energy bytes remain opaque."""
    rows, current, rank = {}, None, None
    for line in raw.splitlines(keepends=True):
        if line.startswith(b"MODEL"):
            require(current is None and re.fullmatch(rb"MODEL\s+[0-9]+\s*", line) is not None, "Malformed MODEL identity")
            rank = int(line.split()[1]); require(rank == len(rows) + 1, "Missing/duplicate MODEL identity")
            current = [line]
        elif current is not None:
            current.append(line)
            if line.startswith(b"ENDMDL"):
                rows[rank] = hashlib.sha256(b"".join(current)).hexdigest(); current = None
        elif line.startswith((b"ATOM  ", b"HETATM", b"ENDMDL")):
            raise ValueError("Model payload outside MODEL identity")
    require(current is None and rows, "Truncated/missing MODEL identity")
    return rows


def local_file(directory, name):
    require(type(name) is str and Path(name).name == name and name not in {".", ".."}, "Unsafe raw artifact filename")
    return directory + "/" + name


def streams(record):
    for key in ("raw_stdout", "raw_stderr", "recognition_log"):
        require(type(record.get(key)) is str and hashlib.sha256(record[key].encode()).hexdigest() == record.get(key + "_sha256"),
                "Opaque native stream checksum changed")


def modern_runs(e, receipt, path, config_sha, *, row=None):
    require(receipt.get("status") == "complete" and receipt.get("native_search_started") is True
            and receipt.get("runs_path") == path and integer(receipt.get("terminal_seed_runs")) == 3,
            "Modern condition is not terminal")
    runs = e.json(path, sha(receipt.get("runs_sha256"))); by_seed = unique_rows(runs, "seed", SEEDS)
    e.read(str(Path(path).parent / "configuration.json"), config_sha)
    require(receipt.get("seed_statuses") == {str(s): by_seed[s].get("status") for s in SEEDS}, "Modern receipt/run status mismatch")
    output = {}
    for seed, run in by_seed.items():
        require(run.get("status") in {"ok", "failed"} and run.get("fingerprint", {}).get("config_sha256") == config_sha,
                "Modern seed missing, pending or rebound")
        if row is not None:
            expected = {"preparation_sha256": row["preparation"]["sha256"], "ligand_sha256": row["ligand"]["sha256"],
                        "receptor_sha256": row["receptor"]["sha256"], "full_assembly_sha256": None}
            prepared = row["modern_vina_prepared_directory"]
            for filename, key in [("reference.json", "reference_sha256"), ("neighborhood.pdb", "receptor_pdb_sha256")]:
                name = prepared + "/" + filename; e.read(name); expected[key] = e.files[name]["sha256"]
            executables = {v["sha256"] for v in e.source(VINA_LOCK)["assets"].values()}
            require(all(run["fingerprint"].get(k) == v for k, v in expected.items())
                    and run["fingerprint"].get("executable_sha256") in executables
                    and run.get("case_id") == row["case_id"] and run.get("receptor") == "neighborhood"
                    and run.get("vina_version") == "AutoDock Vina v1.2.7"
                    and run.get("information_regime") == ("true-site-local" if row["method"] == "oracle" else "predicted-site-local"),
                    "Modern prepared-input lineage changed")
        directory = Path(path).parent.as_posix(); e.read(local_file(directory, run["log_file"]))
        blocks, pose_path = {}, None
        if run["status"] == "ok":
            require(type(run.get("exit_code")) is int and run["exit_code"] == 0, "Successful modern seed has nonzero/unknown exit")
            pose_path = local_file(directory, run["pose_file"])
            blocks = model_hashes(e.read(pose_path, sha(run.get("pose_sha256"))))
            require(integer(run.get("returned_poses"), minimum=1) == len(blocks) <= 20, "Modern model count changed")
        output[seed] = {"status": run["status"], "record": run, "model_hashes": blocks, "pose_file": pose_path}
    return output


def main_workflow(e):
    experiment, modern = e.source(VC_PLAN), e.source(MODERN)
    freeze, main_freeze, manifest = e.source(VC_FREEZE), e.source(MAIN_FREEZE), e.source(INPUTS)
    require(experiment.get("seeds") == SEEDS and experiment.get("cases") == CASES and experiment.get("methods") == ARMS
            and experiment.get("receptor_profiles") == ["strict", "repair"], "Main schedule changed")
    require(experiment["prepared_inputs"]["sha256"] == e.static[INPUTS]["sha256"]
            and main_freeze["configuration_sha256"] == e.static[MODERN]["sha256"]
            and freeze["configuration_sha256"] == e.static[VC_PLAN]["sha256"], "Main freeze/source pins changed")
    expected = {f"{p}/{c}/{m}/box01" for p in ["strict", "repair"] for c in CASES for m in ARMS}
    rows = unique_rows(manifest["conditions"], "condition_id", expected); outcomes = {}
    for identity, row in rows.items():
        profile, case, arm, _ = identity.split("/")
        require((row["receptor_profile"], row["case_id"], row["method"]) == (profile, case, arm)
                and row.get("eligible_workflow_condition") is True, "Main condition identity changed")
        base = f"{MAIN}/{profile}/{case}/{arm}"; path = base + "/vina/box01/runs.json"
        require(row["modern_runs"] == path and row["source_condition"]["path"] == base + "/prepared/condition.json",
                "Main output/source path changed")
        e.pinned(row["source_condition"]); receipt = e.json(base + "/vina/execution.json")
        require(receipt.get("condition_id") == identity and receipt.get("configuration_sha256") == e.static[MODERN]["sha256"]
                and receipt.get("frozen_input_manifest_sha256") == e.static[INPUTS]["sha256"]
                and receipt.get("source_condition") == row["source_condition"], "Modern condition receipt rebound")
        prepared = row["status"] == "prepared"
        require(row["status"] in {"prepared", "preparation_failed", "predictor_unavailable"}, "Unknown source preparation")
        if prepared:
            for key in ["ligand", "receptor", "preparation", "canonical_atom_map"]:
                e.pinned(row[key])
            values = modern_runs(e, receipt, path, e.static[MODERN]["sha256"], row=row)
        else:
            require(receipt.get("status") == "prepared_input_failure" and receipt.get("native_search_started") is False
                    and receipt.get("eligible_seed_slots") == SEEDS and receipt.get("no_runs_or_poses_invented") is True,
                    "Source-failed modern condition needs actual three-slot failure receipt")
            values = {s: {"status": "prepared_input_failure", "record": receipt, "model_hashes": {}, "pose_file": None} for s in SEEDS}
        for seed in SEEDS:
            outcomes[(identity, seed, ALGORITHMS[0])] = values[seed]
            directory = base + f"/vina_carb/box01/s{seed}"
            require(row["vina_carb_outputs"].get(str(seed)) == directory and set(row["vina_carb_outputs"]) == {str(s) for s in SEEDS},
                    "VC seed-output identity changed")
            record = e.json(directory + "/execution.json")
            require(canonical_bytes(record.get("condition")) == canonical_bytes(row) and record.get("seed") == seed and type(record["seed"]) is int
                    and record.get("configuration_sha256") == e.static[VC_PLAN]["sha256"]
                    and record.get("freeze_receipt_sha256") == e.static[VC_FREEZE]["sha256"]
                    and record.get("driver_source_sha256") == experiment["adapter_source"]["sha256"]
                    and record.get("atomic_accuracy_computed") is False and record.get("experimental_ligand_coordinates_read") is False,
                    "VC receipt source/seed/freeze rebound")
            require(record.get("status") in VC_TERMINAL and type(record.get("full_search_started")) is bool, "VC seed remains pending")
            require(prepared != (record["status"] == "prepared_input_failure"), "VC source-failure status changed")
            if not prepared:
                require(record["full_search_started"] is False, "Unprepared VC job claims execution")
            if record["full_search_started"]:
                require(record.get("executable_sha256") == experiment["native_binary_sha256"]
                        and canonical_bytes(record.get("parameters")) == canonical_bytes(experiment["parameters"]),
                        "VC started-search method metadata changed")
                command = record.get("command")
                require(isinstance(command, list) and len(command) == 3 and command[1:] == ["--config", "configuration.txt"],
                        "VC recorded command changed")
                require(record.get("receptor_sha256") == row["receptor"]["sha256"], "VC copied receptor identity changed")
                e.read(directory + "/receptor.pdbqt", row["receptor"]["sha256"])
                e.read(directory + "/ligand.pdbqt", sha(record.get("exported_ligand_sha256")))
                e.read(directory + "/configuration.txt")
            for field in ["preflight", "native_attempt"]:
                if field in record:
                    streams(record[field])
            for name, digest in record.get("artifacts", {}).items():
                e.read(local_file(directory, name), digest)
            blocks, pose_path = {}, None
            if record["status"] == "ok":
                needed = {"poses.pdbqt", "ligand.pdbqt", "ligand.pdbqt.map.json", "receptor.pdbqt", "configuration.txt",
                          "VC_log.txt", "search.log", "stdout.log", "stderr.log"}
                require(needed <= set(record.get("artifacts", {})), "Successful VC omitted raw artifact identities")
                require(record["full_search_started"] is True and type(record.get("exit_code")) is int and record["exit_code"] == 0,
                        "Successful VC lacks successful actual execution")
                require(record.get("executable_sha256") == experiment["native_binary_sha256"], "VC executable provenance changed")
                pose_path = directory + "/poses.pdbqt"
                blocks = model_hashes(e.read(pose_path, sha(record.get("artifacts", {}).get("poses.pdbqt"))))
                require(integer(record.get("returned_poses"), minimum=1) == len(blocks) <= 20, "VC model count changed")
            outcomes[(identity, seed, ALGORITHMS[1])] = {"status": record["status"], "record": record,
                                                        "model_hashes": blocks, "pose_file": pose_path}
    progress = e.json(MAIN + "/paired_execution_progress.json")
    require(progress.get("status") == "complete" and integer(progress.get("planned_eligible_method_seed_slots")) == 216
            and integer(progress.get("terminal_eligible_seed_slots")) == 216 and progress.get("atomic_accuracy_evaluated") is False,
            "Main progress is not fully terminal")
    expected_records = {(i, "vina_1_2_7", None) for i in rows} | {(i, "vina_carb_1_0", s) for i in rows for s in SEEDS}
    actual = []
    for record in progress["records"]:
        key = (record.get("condition_id"), record.get("method"), record.get("seed"))
        require(key in expected_records, "Unexpected main progress identity")
        identity, method, seed = key
        expected_path = str(Path(rows[identity]["modern_runs"]).parent.parent / "execution.json") if seed is None else rows[identity]["vina_carb_outputs"][str(seed)] + "/execution.json"
        require(record.get("execution_receipt") == expected_path, "Main progress receipt path changed")
        expected_status = ("complete" if rows[identity]["status"] == "prepared" else "prepared_input_failure") if seed is None else outcomes[(identity, seed, ALGORITHMS[1])]["status"]
        require(record.get("status") == expected_status, "Main progress disagrees with receipt")
        actual.append(key)
    require(len(actual) == len(set(actual)) == 144 and set(actual) == expected_records, "Main progress omits/duplicates slots")
    return rows, outcomes, progress


def assembly_workflow(e):
    config, index = e.source(ASSEMBLY_CONFIG), e.source(ASSEMBLY_INPUTS)
    freeze, execution_freeze = e.source(ASSEMBLY_FREEZE), e.source(ASSEMBLY_DRIVER_FREEZE)
    require(config.get("seeds") == SEEDS and freeze["configuration_sha256"] == e.static[ASSEMBLY_CONFIG]["sha256"]
            and execution_freeze["configuration_sha256"] == e.static[ASSEMBLY_CONFIG]["sha256"]
            and execution_freeze["prepared_conditions_sha256"] == e.static[ASSEMBLY_INPUTS]["sha256"],
            "Assembly source/configuration freeze changed")
    expected = {(p, c, a) for p in ["strict_monomer", "repair_monomer"] for c in ["3ZWE", "5HZA"] for a in ARMS}
    rows = index["conditions"]
    require(len(rows) == 12 and {(r["receptor_profile"], r["case_id"], r["method"]) for r in rows} == expected,
            "Assembly condition denominator changed")
    receipts, slots = {}, []
    for row in rows:
        identity = "/".join([row["receptor_profile"], row["case_id"], row["method"], "box01"])
        base = "/".join([ASSEMBLY, row["receptor_profile"], row["case_id"], row["method"]])
        require(identity not in receipts and row["condition_path"] == base + "/prepared/condition.json", "Assembly source identity changed")
        condition = e.json(row["condition_path"], row["condition_sha256"])
        receipt_path = base + "/vina/execution.json"; receipt = e.json(receipt_path)
        require(receipt.get("condition_id") == identity and receipt.get("configuration_sha256") == e.static[ASSEMBLY_CONFIG]["sha256"]
                and receipt.get("prepared_conditions_sha256") == e.static[ASSEMBLY_INPUTS]["sha256"]
                and receipt.get("source_condition") == {"path": row["condition_path"], "sha256": row["condition_sha256"]},
                "Assembly execution source/configuration rebound")
        prepared = integer(row["prepared_jobs"])
        require(prepared in {0, 1}, "Assembly preparation denominator changed")
        if prepared:
            jobs = condition["jobs"]
            require(len(jobs) == 1 and jobs[0].get("status") == "prepared", "Assembly prepared box changed")
            prepared_dir = base + "/prepared/" + jobs[0]["prepared_directory"]
            lineage = {"case_id": row["case_id"], "method": row["method"], "modern_vina_prepared_directory": prepared_dir,
                       "preparation": {"path": prepared_dir + "/preparation.json", "sha256": jobs[0]["preparation_sha256"]}}
            e.pinned(lineage["preparation"])
            for key, filename in [("ligand", "ligand.pdbqt"), ("receptor", "neighborhood.pdbqt")]:
                name = prepared_dir + "/" + filename; e.read(name); lineage[key] = e.files[name]
            modern_runs(e, receipt, base + "/vina/box01/runs.json", e.static[ASSEMBLY_CONFIG]["sha256"], row=lineage)
        else:
            require(receipt.get("status") == "prepared_input_failure" and receipt.get("native_search_started") is False
                    and receipt.get("eligible_seed_slots") == SEEDS and receipt.get("no_runs_or_poses_invented") is True,
                    "Assembly failed source lacks actual three-slot receipt")
        receipts[identity] = (receipt_path, receipt["status"])
        slots.extend(identity + "/s" + str(s) for s in SEEDS)
    progress = e.json(ASSEMBLY + "/execution_progress.json")
    require(progress.get("status") == "complete" and integer(progress.get("planned_eligible_seed_slots")) == 36
            and integer(progress.get("terminal_eligible_seed_slots")) == 36
            and progress.get("atomic_or_site_accuracy_evaluated") is False, "Assembly progress is not terminal")
    observed = unique_rows(progress["records"], "condition_id", receipts)
    for identity, record in observed.items():
        require((record.get("execution_receipt"), record.get("status")) == receipts[identity], "Assembly progress receipt changed")
    return slots, progress


def haddock_artifacts(e, directory, identity, physical):
    """Declared stage/file hashes only; models and energy tables stay opaque."""
    for stage in physical.get("stages", physical.get("retained_partial_stages", [])):
        name = stage["stage"]
        require(type(name) is str and re.fullmatch(r"[0-9]+_[A-Za-z0-9_]+", name), "Unsafe HADDOCK stage identity")
        base = directory + "/run/" + name
        if stage.get("io_sha256") is not None:
            e.read(base + "/io.json", sha(stage["io_sha256"]))
        models = unique_rows(stage.get("models", []), "file_name")
        require(integer(stage.get("actual_physical_model_files")) == len(models), "HADDOCK declared model inventory changed")
        for filename, model in models.items():
            path = local_file(base, filename)
            require(filename.endswith((".pdb", ".pdb.gz")) and model.get("model_id") == identity + "/" + name + "/" + filename.removesuffix(".gz"),
                    "HADDOCK model identity/path changed")
            e.read(path, sha(model.get("sha256")))


def haddock_workflow(e):
    plan, freeze, recovery, recovery_freeze = e.source(HADD_PLAN), e.source(HADD_FREEZE), e.source(RECOVERY_PLAN), e.source(RECOVERY_FREEZE)
    require(plan.get("seeds") == SEEDS and freeze["prepared_plan"]["sha256"] == e.static[HADD_PLAN]["sha256"]
            and recovery_freeze["recovery_plan"]["sha256"] == e.static[RECOVERY_PLAN]["sha256"], "HADDOCK freeze/source rebound")
    for key, name in [("original_plan", HADD_PLAN), ("original_execution_freeze", HADD_FREEZE)]:
        require(recovery.get(key, {}).get("sha256") == e.static[name]["sha256"], "Recovery original source lineage changed")
    for key in ["recovery_wrapper", "recovery_policy", "orphan_completion_proof"]:
        e.pinned(recovery[key])
    for pin in recovery.get("existing_terminal_receipts", []):
        e.pinned(pin)
    proof = e.source(HADD + "/orchestration_recovery/orphan_completion_proof.json")
    for pin in proof.get("raw_artifacts", []):
        require(pin["path"].startswith(HADD + "/"), "Recovered raw proof escaped HADDOCK root")
        e.pinned(pin)
    jobs = unique_rows(plan["jobs"], "job_id")
    expected = {f"{p}/{c}/{a}/box01/s{s}" for p in ["strict", "repair"] for c in CASES for a in ARMS for s in SEEDS}
    require(set(jobs) == expected, "HADDOCK108 identities changed")
    progress = e.json(HADD + "/execution_progress.json")
    for field, name in [("prepared_plan", HADD_PLAN), ("root_freeze", HADD_FREEZE),
                        ("execution_recovery_freeze", RECOVERY_FREEZE), ("recovery_plan", RECOVERY_PLAN)]:
        require(progress.get(field) == e.static[name], "HADDOCK progress transport pin changed")
    require(type(progress.get("schema_version")) is int and progress["schema_version"] == 2
            and integer(progress.get("seed_cell_denominator")) == 108 and integer(progress.get("terminal_cells")) == 108
            and progress.get("all_terminal") is True and progress.get("no_accuracy_evaluation") is True,
            "HADDOCK progress is not fully terminal")
    observed = unique_rows(progress["jobs"], "job_id", jobs)
    ready = 0
    for identity, job in jobs.items():
        row = observed[identity]
        require(job.get("eligible_workflow_pool") is True and job.get("condition_id") == identity.rsplit("/", 1)[0],
                "HADDOCK frozen job eligibility/condition changed")
        require(row.get("terminal") is True and row.get("eligible_workflow_pool") is True, "HADDOCK cell is not terminal/eligible")
        if job["status"] != "ready":
            require(row.get("status") == job["status"] and type(row.get("physical_models")) is int and row["physical_models"] == 0
                    and row.get("no_native_execution") is True and "execution_receipt" not in row,
                    "Nonready HADDOCK cell requires original progress-only failure")
            continue
        ready += 1
        profile, case, arm, _, seed_name = identity.split("/")
        directory = f"{HADD}/{profile}/{case}/{arm}/{seed_name}"
        require(job.get("directory") == directory and integer(job.get("seed")) == int(seed_name[1:]), "HADDOCK native job path/seed changed")
        present = [name for name in ["execution.json", "wrapper_failure.json"] if (e.root / directory / name).exists()]
        require(len(present) == 1, "Missing/conflicting HADDOCK physical terminal receipts")
        pin = row.get("execution_receipt")
        require(isinstance(pin, dict) and pin.get("path") in {directory + "/execution.json", directory + "/wrapper_failure.json"},
                "Ready HADDOCK cell lacks physical terminal receipt")
        physical = e.pinned(pin, json_value=True)
        require(canonical_bytes({k: v for k, v in row.items() if k not in {"execution_receipt", "resource_usage"}}) == canonical_bytes(physical),
                "HADDOCK progress differs from physical receipt")
        require(physical.get("prepared_plan") == e.static[HADD_PLAN]
                and physical.get("config") == job["config"] and physical.get("no_accuracy_evaluation") is True,
                "HADDOCK physical source/configuration rebound")
        e.pinned(job["config"])
        for source in job.get("input_files", []):
            e.pinned(source)
        require(physical.get("status") in {"completed", "native_execution_failed", "complete_outputs_recovered", "orchestration_failed"},
                "HADDOCK native cell has unknown/pending status")
        if physical["status"] == "orchestration_failed":
            require(pin["path"].endswith("/wrapper_failure.json") and physical.get("native_retry") is False
                    and physical.get("process_return_code") is None, "Wrapper failure provenance changed")
        else:
            require(pin["path"].endswith("/execution.json") and physical.get("root_freeze") == e.static[HADD_FREEZE],
                    "Native receipt original freeze changed")
            require(physical.get("helper_source") == plan["driver"] and physical.get("CNS_binary") == plan["CNS_binary"]
                    and integer(physical.get("ncores")) == 1 and integer(physical.get("requested_initial_model_slots")) == 20
                    and physical.get("failed_model_slots_not_replaced") is True and physical.get("clustering_performed") is False,
                    "HADDOCK recorded native source/protocol changed")
        if physical["status"] == "completed":
            require(type(physical.get("process_return_code")) is int and physical["process_return_code"] == 0
                    and physical.get("final_score_metadata_complete") is True, "Completed HADDOCK receipt lacks successful terminal metadata")
        if physical["status"] == "complete_outputs_recovered":
            require(identity == "strict/2RDK/fpocket/box01/s104729" and physical.get("native_rerun") is False and physical.get("process_return_code") is None
                    and all(k in physical and physical[k] is None for k in
                            ["process_return_code", "process_CPU_seconds", "external_measured_wall_seconds"]),
                    "Recovered first job must preserve unknown original telemetry")
            require(physical.get("recovery_receipt") == e.static[HADD + "/orchestration_recovery/adoption_receipt.json"],
                    "Recovered first-job adoption receipt changed")
        haddock_artifacts(e, directory, identity, physical)
        if "resource_usage" in row:
            resource = e.pinned(row["resource_usage"], json_value=True)
            require(row["resource_usage"]["path"] == directory + "/measured_resources.json"
                    and resource.get("job_id") == identity and integer(resource.get("native_cpus_per_job")) == 1
                    and resource.get("original_plan") == e.static[HADD_PLAN]
                    and resource.get("frozen_helper") == plan["driver"], "HADDOCK resource-sidecar identity changed")
    require(ready == 39 and len(jobs) - ready == 69, "HADDOCK ready/progress-only denominator changed")
    return list(jobs), progress


def union_workflow(e, rows, outcomes):
    plan, freeze, schedule = e.source(UNION_PLAN), e.source(UNION_FREEZE), e.source(SCHEDULE)
    plan_sha = e.static[UNION_PLAN]["sha256"]
    require(plan["source_inputs"]["sha256"] == e.static[INPUTS]["sha256"]
            and plan["source_schedule"]["sha256"] == e.static[SCHEDULE]["sha256"], "Union source/schedule rebound")
    # The original freeze names its configuration, not a recomputed report.
    require(freeze.get("configuration_sha256") == plan_sha, "Union execution freeze rebound")
    index = e.json(UNION + "/source_pool.json")
    require(index.get("configuration_sha256") == plan_sha and index.get("source_inputs_sha256") == e.static[INPUTS]["sha256"]
            and index.get("method_lock_sha256") == plan["method_lock"]["sha256"]
            and integer(index.get("eligible_seed_pools")) == 108 and integer(index.get("eligible_source_method_slots")) == 216
            and integer(index.get("terminal_source_method_slots")) == 216
            and index.get("all216_source_slots_terminal_and_verified") is True and index.get("pending_source_slots") == []
            and index.get("native_coordinates_used") is False and index.get("atomic_accuracy_computed") is False,
            "Final union source index is pending or rebound")
    pool_ids = {i + "/s" + str(s) for i in rows for s in SEEDS}
    pools = unique_rows(index["pools"], "pool_id", pool_ids)
    unique_rows(schedule["pools"], "pool_id", pool_ids)
    sources = unique_rows(index["poses"], "pose_id")
    expected_sources = {}
    for identity, row in rows.items():
        for seed in SEEDS:
            pool_id = identity + "/s" + str(seed); pool = pools[pool_id]
            require(pool.get("condition_id") == identity and type(pool.get("seed")) is int and pool["seed"] == seed
                    and pool.get("case_id") == row["case_id"] and pool.get("receptor_profile") == row["receptor_profile"]
                    and pool.get("arm") == row["method"] and pool.get("eligible_workflow_pool") is True,
                    "Union pool identity changed")
            require(set(pool["source_statuses"]) == set(ALGORITHMS) and set(pool["source_record_sha256"]) == set(ALGORITHMS)
                    and set(pool["original_source_orders"]) == set(ALGORITHMS), "Union source-method identity set changed")
            for algorithm_index, algorithm in enumerate(ALGORITHMS):
                outcome = outcomes[(identity, seed, algorithm)]
                source_digest = canonical_sha256(outcome["record"])
                require(pool["source_statuses"][algorithm] == outcome["status"]
                        and pool["source_record_sha256"][algorithm] == source_digest, "Union source status/record lineage changed")
                expected_order = [pool_id + f"/{algorithm}/p{rank:02d}" for rank in outcome["model_hashes"]]
                require(pool["original_source_orders"][algorithm] == expected_order, "Union omitted/reordered retained source MODEL identities")
                for rank, model_digest in outcome["model_hashes"].items():
                    pose_id = pool_id + f"/{algorithm}/p{rank:02d}"
                    expected_sources[pose_id] = True
                    require(pose_id in sources, "Missing retained source pose identity")
                    source = sources[pose_id]
                    expected = {"pool_id": pool_id, "condition_id": identity, "complex_id": row["case_id"],
                                "receptor_profile": row["receptor_profile"], "arm": row["method"], "seed": seed,
                                "source_algorithm": algorithm, "source_algorithm_rank": algorithm_index,
                                "model_rank": rank, "source_pose_file": outcome["pose_file"],
                                "source_pose_file_sha256": e.files[outcome["pose_file"]]["sha256"],
                                "source_model_sha256": model_digest, "source_record_sha256": source_digest,
                                "source_ligand_sha256": row["ligand"]["sha256"], "source_receptor_sha256": row["receptor"]["sha256"],
                                "source_preparation_sha256": row["preparation"]["sha256"],
                                "source_canonical_atom_map_sha256": row["canonical_atom_map"]["sha256"],
                                "candidate_class_rank": row["candidate_class_rank"], "candidate_original_rank": row["candidate_original_rank"]}
                    require(all(type(source.get(k)) is type(v) and source[k] == v for k, v in expected.items()), "Union retained source identity/hash rebound")
                    require(source.get("pose_identity_sha256") == canonical_sha256({k: v for k, v in source.items() if k != "pose_identity_sha256"}),
                            "Union canonical pose identity changed")
                    if algorithm_index == 0:
                        require(source.get("source_export_map") is None, "Modern source export identity changed")
                    else:
                        pin = source.get("source_export_map")
                        require(isinstance(pin, dict) and pin.get("path") == str(Path(outcome["pose_file"]).parent / "ligand.pdbqt.map.json"),
                                "VC export-map source path changed")
                        e.pinned(pin)
    require(set(sources) == set(expected_sources) and integer(index.get("retained_pose_count")) == len(sources),
            "Final union source set omits/duplicates/adds retained models")
    manifest = e.json(UNION + "/scores_manifest.json")
    records = manifest["records"]
    require(isinstance(records, dict) and manifest.get("configuration_sha256") == plan_sha and set(records) == set(sources),
            "Retained pose-ID set must equal fixed-score record keys exactly")
    physical = {p.relative_to(e.root / UNION / "scores").as_posix()[:-5]
                for p in (e.root / UNION / "scores").rglob("*.json")}
    require(physical == set(records), "Orphan/missing fixed-score files")
    successful = 0
    for identity, digest in records.items():
        e.pkg.portable_path(identity)
        record = e.json(UNION + "/scores/" + identity + ".json", sha(digest))
        require(canonical_bytes(record.get("source")) == canonical_bytes(sources[identity]) and record.get("method_lock_sha256") == index["method_lock_sha256"]
                and record.get("coordinate_optimization") is False and record.get("native_coordinates_used") is False,
                "Fixed-score source/method policy rebound")
        require(record.get("status") in {"ok", "pose_scoring_failed"}, "Fixed-score record remains pending or has unknown producer status")
        if record["status"] == "ok":
            require(isinstance(record.get("score"), dict), "Missing opaque score packet")
            streams(record["score"]); successful += 1
        else:
            require(type(record.get("error")) is str and record["error"], "Fixed-score failure lacks recorded reason")
            if "native_attempt" in record:
                streams(record["native_attempt"])
    progress = e.json(UNION + "/scoring_progress.json")
    require(progress.get("configuration_sha256") == plan_sha
            and progress.get("source_pool_sha256") == e.files[UNION + "/source_pool.json"]["sha256"]
            and progress.get("score_manifest_sha256") == e.files[UNION + "/scores_manifest.json"]["sha256"]
            and integer(progress.get("eligible_seed_pools")) == 108 and integer(progress.get("eligible_rescore_method_slots")) == 216
            and integer(progress.get("indexed_pose_count")) == len(sources) and integer(progress.get("recorded_pose_count")) == len(records)
            and integer(progress.get("successful_pose_scores")) == successful and integer(progress.get("failed_pose_scores")) == len(records) - successful
            and progress.get("all216_source_slots_terminal_and_verified") is True
            and progress.get("ready_for_separate_accuracy_evaluation") is True and progress.get("atomic_accuracy_computed") is False,
            "Final fixed-score progress is stale/pending or changed")
    return list(sources), list(pools)


def process_absent(pid):
    """Current-host existence query only; no terminating/stopping signals."""
    integer(pid, minimum=2)
    require(pid != os.getpid(), "Validator cannot certify itself as a stopped producer")
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    except (PermissionError, OverflowError) as error:
        raise ValueError("Cannot establish current producer PID absence") from error
    return False


def benchmark_command_tags(root, command):
    """Conservative writer candidates, not exhaustive historical PID discovery."""
    try:
        words = shlex.split(command)
    except ValueError:
        words = command.split()
    if not words:
        return []
    executable = Path(words[0]).name.casefold()
    runners = ("python", "pypy", "sh", "bash", "zsh", "dash", "env", "nohup", "timeout", "time", "nice")
    native = ("vina", "haddock3", "cns", "fixed_pose_score")
    is_runner = any(executable == n or executable.startswith(n + ".") or
                    (n in {"python", "pypy"} and executable.startswith(n)) for n in runners)
    is_native = any(executable == n or executable.startswith(n + "-") or executable.startswith(n + "_") for n in native)
    is_native = is_native or any(words[0].startswith(str(root / name)) for name in
                                 ["tools/vina-1.2.7", "tools/vina-carb-1.0", "tools/haddock3-2026.10.1"])
    if not (is_runner or is_native):
        return []
    markers = {
        "main_driver": ("paired_execution_driver.py",),
        "assembly_driver": ("structural_assembly_pose_v3/execution_driver.py",),
        "haddock_driver": ("haddock_cohort_resume", "haddock_paired_cohort", "run_resume_pool("),
        "union_driver": ("paired_common_pose_rescoring", "score_available_union("),
        "benchmark_result_path": (MAIN + "/", ASSEMBLY + "/", HADD + "/"),
        "benchmark_native_path": (str(root / "tools/vina-1.2.7"), str(root / "tools/vina-carb-1.0"),
                                  str(root / "tools/haddock3-2026.10.1")),
    }
    tags = [tag for tag, literals in markers.items() if any(literal in command for literal in literals)]
    if is_native:
        tags.append("conservative_native_writer_candidate")
    return sorted(tags)


def current_benchmark_processes(root):
    """Read-only ps query; discard unrelated commands and all raw command text."""
    try:
        result = subprocess.run(["/bin/ps", "-ww", "-axo", "pid=,command="], capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError("Cannot obtain current benchmark-writer command inventory") from error
    require(result.returncode == 0, "Current benchmark-writer command inventory failed")
    candidates = []
    for line in result.stdout.splitlines():
        parts = line.strip().split(None, 1)
        require(len(parts) == 2 and parts[0].isdigit(), "Malformed current process inventory row")
        pid = int(parts[0])
        if pid == os.getpid():
            continue
        tags = benchmark_command_tags(root, parts[1])
        if tags:
            candidates.append({"pid": pid, "tags": tags, "command_sha256": hashlib.sha256(parts[1].encode()).hexdigest()})
    return sorted(candidates, key=lambda row: row["pid"])


def producer_checks(e, context):
    require(type(context.get("schema_version")) is int and context["schema_version"] == 2
            and context.get("scope") == "actual_current_host_producers_only"
            and context.get("repository_root") == str(e.root) and context.get("hostname") == socket.gethostname()
            and context.get("source_review_and_observed_completion_account_for_writer_lifetimes") is True
            and "all_writer_descendants_accounted_for" not in context,
            "Explicit source-bound current-host completion context required; historical PIDs are not portable gates")
    roles = unique_rows(context["producers"], "role", DRIVERS)
    all_pids = []
    for role, record in roles.items():
        require(record.get("driver_source") == e.static[DRIVERS[role]], "Producer launch source identity changed")
        e.pinned(record["driver_source"]); e.pinned(record["launch_provenance"])
        complete = record.get("historical_pid_inventory_complete")
        pids = record.get("actual_known_pids")
        require(type(complete) is bool and "actual_known_pids" in record and "pids" not in record
                and (pids is None or isinstance(pids, list)), "Explicit known/unknown PID inventory required")
        if role in {"main", "assembly"}:
            require(isinstance(pids, list) and pids, "Actual main/assembly driver PID inventory required")
        if complete:
            require(isinstance(pids, list) and pids and record.get("missing_historical_pids") == []
                    and record.get("pid_inventory_reason") == "complete_recorded_inventory", "Complete PID inventory provenance missing")
        else:
            require("missing_historical_pids" in record and record["missing_historical_pids"] is None
                    and record.get("pid_inventory_reason") == "not_recorded_during_execution", "Unknown historical PIDs must stay explicit null")
        completion = e.pinned(record["operator_completion"], json_value=True)
        require(type(completion.get("schema_version")) is int and completion["schema_version"] == 1
                and completion.get("kind") == "operator_observed_synchronous_tool_completion" and completion.get("role") == role
                and type(completion.get("tool_session_id")) is str and completion["tool_session_id"].strip()
                and type(completion.get("exit_code")) is int and completion["exit_code"] == 0
                and completion.get("observed_complete") is True and completion.get("synchronous_wait_observed") is True
                and completion.get("exit_code_scope") == "tool_session_only_not_original_native_process"
                and completion.get("original_native_process_return_code_inferred") is False
                and completion.get("launch_provenance") == record["launch_provenance"]
                and completion.get("driver_source") == record["driver_source"], "Observed synchronous tool completion proof missing/invalid")
        proof = e.pinned(record["source_wait_proof"], json_value=True)
        source_pins = [e.static[name] for name in WAIT_SOURCES[role]]
        require(type(proof.get("schema_version")) is int and proof["schema_version"] == 1 and proof.get("role") == role
                and proof.get("source_pins") == source_pins and proof.get("wait_contract") == WAIT_CONTRACTS[role]
                and proof.get("all_workers_and_child_writers_joined") is True and proof.get("no_unjoined_writer_lifetimes") is True,
                "Source-pinned synchronous writer wait/join proof missing/invalid")
        for pin in source_pins:
            e.pinned(pin)
        for pid in pids or []:
            integer(pid, minimum=2)
            require(process_absent(pid), "Actual known producer PID remains present; no readiness")
        all_pids.extend(pids or [])
    require(len(all_pids) == len(set(all_pids)), "Duplicate actual known producer/descendant PID")
    candidates = current_benchmark_processes(e.root)
    require(not candidates, "Current benchmark-writer command candidate remains; no readiness: " + json.dumps(candidates, sort_keys=True))
    return roles


def validate(root, plan_path, authorization_path, *, allow_synthetic_fixture=False):
    """Read-only certification; production use requires a source-pinned root grant."""
    pkg = packager(); root = pkg.checked_root(root)
    plan_name = pkg._metadata_name(root, plan_path); auth_name = pkg._metadata_name(root, authorization_path)
    plan_raw = pkg._read(root, plan_name); plan = pkg.validate_plan(load_json(plan_raw), allow_synthetic_fixture=allow_synthetic_fixture)
    synthetic = plan["profile"] == "synthetic_test_fixture"
    if synthetic:
        pkg._fixture_root(root)
    else:
        require(hashlib.sha256(plan_raw).hexdigest() == PLAN_SHA256, "Root-reviewed real packaging plan changed")
    auth_raw = pkg._read(root, auth_name); auth = load_json(auth_raw)
    plan_pin = {"path": plan_name, "bytes": len(plan_raw), "sha256": hashlib.sha256(plan_raw).hexdigest()}
    require(type(auth.get("schema_version")) is int and auth["schema_version"] == 1
            and auth.get("status") == "approved_metadata_only_terminal_validation" and auth.get("validation_authorized") is True
            and auth.get("packaging_plan") == plan_pin and auth.get("packager_sha256") == PACKAGER_SHA256
            and auth.get("validator_sha256") == hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "Explicit root validation grant must pin this validator, packager and exact plan")
    for pin in plan["source_contracts"]:
        pkg._checked(root, pin)
    if synthetic:
        pins = auth.get("source_metadata_pins")
        static = unique_rows(pins, "path", STATIC_NAMES)
        for pin in pins:
            pkg._pin(pin)
    else:
        require("source_metadata_pins" not in auth, "Public frozen-source pins cannot be rebound")
        distribution_raw = pkg._read(root, DISTRIBUTION)
        require(hashlib.sha256(distribution_raw).hexdigest() == DISTRIBUTION_SHA256, "Published source map changed")
        published = unique_rows(load_json(distribution_raw)["root_files"], "path")
        require(set(STATIC_NAMES) <= set(published), "Published source map lacks readiness prerequisites")
        static = {name: published[name] for name in STATIC_NAMES}
    e = Evidence(root, pkg, static)
    e.read(plan_name, plan_pin["sha256"]); e.read(auth_name, hashlib.sha256(auth_raw).hexdigest())
    context_pin = auth.get("producer_context"); pkg._pin(context_pin); pkg._metadata_name(root, context_pin["path"])
    context = e.pinned(context_pin, json_value=True)
    roles = producer_checks(e, context)
    pkg._ensure_absent(root, plan["must_be_absent"])
    # No source/raw inventory occurs until the actual-host stopped-writer gate.
    before = pkg.inventory(root, plan); historical = pkg.verify_historical(root, plan, before)
    for name in STATIC_NAMES:
        e.pinned(static[name])
    rows, outcomes, main_progress = main_workflow(e)
    assembly_slots, assembly_progress = assembly_workflow(e)
    haddock_slots, _ = haddock_workflow(e)
    retained, pools = union_workflow(e, rows, outcomes)
    for role, progress in [("main", main_progress), ("assembly", assembly_progress)]:
        require(integer(progress.get("driver_pid"), minimum=2) in (roles[role]["actual_known_pids"] or []),
                "Current producer context does not include actual driver PID")
    require(pkg.inventory(root, plan) == before, "Full included/excluded output snapshot changed")
    e.recheck()
    for pin in plan["source_contracts"]:
        pkg._checked(root, pin)
    pkg._ensure_absent(root, plan["must_be_absent"]); producer_checks(e, context)
    require(pkg._read(root, plan_name) == plan_raw and pkg._read(root, auth_name) == auth_raw, "Root validation grant/plan changed")
    return {"schema_version": 1, "status": "all_structural_sources_terminal_metadata_only",
            "accuracy_computed": False, "numeric_pose_or_energy_payloads_interpreted": False,
            "float_number_tokens_retained_without_numeric_conversion": True,
            "writers_stopped": True, "pending": [], "source_snapshot_sha256": pkg.snapshot_sha256(before),
            "counts": {"main": len(outcomes), "assembly": len(assembly_slots), "haddock": len(haddock_slots),
                       "union_sources": len(outcomes), "pools": len(pools), "retained": len(retained), "recorded": len(retained)},
            "unique_retained_pose_ids_equal_fixed_score_record_keys": True,
            "all_returned_source_MODEL_rank_ids_preserved": True,
            "slot_identity_sha256": {"main": canonical_sha256(sorted([list(k) for k in outcomes])),
                                     "assembly": canonical_sha256(sorted(assembly_slots)), "haddock": canonical_sha256(sorted(haddock_slots)),
                                     "union_pools": canonical_sha256(sorted(pools)), "retained_pose_ids": canonical_sha256(sorted(retained))},
            "historical_first_job": historical, "source_contracts": plan["source_contracts"],
            "metadata_and_opaque_payload_file_pins": sorted(e.files.values(), key=lambda p: p["path"]),
            "producer_context": context_pin, "producer_quiescence_checks": {
                "scope": "actual_current_host_producers_only", "host": socket.gethostname(),
                "checked_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "all_actual_known_pids_absent_before_and_after": True,
                "filtered_current_benchmark_command_candidates_absent_before_and_after": True,
                "synchronous_completion_and_source_wait_proofs_verified": True,
                "original_native_process_return_codes_inferred_from_tool_exit": False,
                "process_scan_is_not_exhaustive_historical_PID_discovery": True,
                "portable_replay_must_not_query_historical_pids": True,
                "actual_known_pids": {role: record["actual_known_pids"] for role, record in roles.items()},
                "historical_pid_inventory_complete": {role: record["historical_pid_inventory_complete"] for role, record in roles.items()},
                "missing_historical_pids": {role: record["missing_historical_pids"] for role, record in roles.items()},
                "operator_completion_receipts": {role: record["operator_completion"] for role, record in roles.items()},
                "source_wait_proofs": {role: record["source_wait_proof"] for role, record in roles.items()}},
            "validation_authorization": {"path": auth_name, "bytes": len(auth_raw), "sha256": hashlib.sha256(auth_raw).hexdigest()},
            "packaging_plan": plan_pin, "validator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "packager_sha256": PACKAGER_SHA256, "synthetic_test_fixture": synthetic,
            "scientific_reader_or_native_calls": False,
            "interpretation": "Terminal identity/byte preservation only; no numerical score validation, geometry validation or accuracy authorization."}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    pkg = packager(); root = pkg.checked_root(args.root)
    output = pkg._metadata_name(root, args.output)
    require(output not in {pkg._metadata_name(root, args.plan), pkg._metadata_name(root, args.authorization)}, "Output cannot overwrite validation inputs")
    require(not (root / output).exists() and not (root / output).is_symlink(), "Readiness output must be fresh")
    result = validate(root, args.plan, args.authorization)
    # Exclusive creation; parent must already exist and contain no symlink.
    target = root / output
    require(target.parent.is_dir() and not any(p.is_symlink() for p in (target.parent, *target.parent.parents)), "Unsafe readiness output parent")
    raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    with target.open("xb") as stream:
        stream.write(raw)
    print(json.dumps({"readiness": output, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
