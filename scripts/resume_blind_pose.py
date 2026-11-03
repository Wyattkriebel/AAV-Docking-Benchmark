#!/usr/bin/env python3
"""Audit and recover the frozen blind-pose execution without changing its plan.

The default is read-only. --execute must run in a retained foreground session.
Each worker owns one fresh seed shard; only this parent merges canonical boxes.
No terminal outcome, including a failure, is rerun. Interrupted costs stay unknown.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


TERMINAL = {"ok", "failed"}
INTERRUPTED = {"running", "pending", "cancelled", "interrupted"}


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def _json(path):
    return json.loads(Path(path).read_text())


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def require_stopped(root):
    """A stale lock alone is insufficient evidence that its workers stopped."""
    result = root / "results/blind_pose"
    for path in (result / "execution_driver.lock", result / "recovery_driver.lock"):
        if path.exists() and _alive(_json(path)["pid"]):
            raise ValueError(f"Execution driver still alive: {path}")
    process_list = subprocess.check_output(["ps", "-Ao", "pid=,command="], text=True)
    executable = str(root / "tools/vina-1.2.7")
    for line in process_list.splitlines():
        fields = line.strip().split(None, 1)
        if len(fields) == 2 and int(fields[0]) != os.getpid():
            command = fields[1]
            if command.startswith(executable + " ") or command == executable:
                raise ValueError(f"Vina process still alive: {fields[0]}")


def validate_rows(root, prepared, directory, config_path, *, canonical):
    """Validate complete atoms, invocation, input hashes and terminal artifacts."""
    from aav_benchmark.engine import read_poses
    from aav_benchmark.provenance import sha256
    from aav_benchmark.qualified_pose import _pose_log, _validate_command
    manifest = directory / "runs.json"
    if not manifest.exists():
        return []
    config, prep = _json(config_path), _json(prepared / "preparation.json")
    reference = _json(prepared / "reference.json")
    if sha256(directory / "configuration.json") != sha256(config_path):
        raise ValueError("Run directory changed the frozen configuration")
    if len(prep["receptors"]) != 1 or prep["receptors"][0]["id"] != "neighborhood":
        raise ValueError("Recovery expects exactly the frozen neighborhood receptor")
    receptor = prep["receptors"][0]
    common = {"config_sha256": sha256(config_path), "preparation_sha256": sha256(prepared / "preparation.json"),
              "reference_sha256": sha256(prepared / "reference.json"),
              "ligand_sha256": sha256(prepared / prep["ligand_file"]),
              "receptor_sha256": sha256(prepared / receptor["pdbqt_file"]),
              "receptor_pdb_sha256": sha256(prepared / receptor["pdb_file"]),
              "full_assembly_sha256": None, "executable_sha256": sha256(root / "tools/vina-1.2.7")}
    rows = _json(manifest)
    seen = set()
    for row in rows:
        identity = (row["receptor"], row["seed"])
        if identity in seen or identity[0] != "neighborhood" or identity[1] not in config["seeds"]:
            raise ValueError("Unknown or duplicate seed identity")
        seen.add(identity)
        if row.get("fingerprint") != common:
            raise ValueError("Run inputs, executable or configuration changed")
        if (row.get("case_id") != prep["complex_id"] or row.get("family_id") != prep["family_id"]
                or row.get("information_regime") != prep["execution_information_regime"]
                or row.get("vina_version") != "AutoDock Vina v1.2.7"):
            raise ValueError("Run identity or information regime changed")
        _validate_command(row, prep, config, prep["condition_id"], prep["method"], prepared.name)
        command = row["command"]
        if Path(command[0]).resolve() != (root / "tools/vina-1.2.7").resolve():
            raise ValueError("Actual command executable differs from the pinned binary")
        actual_output = Path(command[command.index("--out") + 1])
        for name in (row["pose_file"], row["log_file"]):
            if Path(name).name != name or name in {".", ".."}:
                raise ValueError("Unsafe raw artifact filename")
        if actual_output != directory / row["pose_file"]:
            # A canonical merge retains the actual shard command, never rewrites it.
            provenance = row.get("shard_provenance", {})
            archive = directory / provenance.get("archived_run_manifest", "MISSING")
            if not canonical or not archive.is_file() or sha256(archive) != provenance.get("source_run_manifest_sha256"):
                raise ValueError("Relocated command lacks its archived shard manifest")
            source = [r for r in _json(archive) if (r["receptor"], r["seed"]) == identity]
            preserved_fields = ("command", "fingerprint", "status", "pose_sha256", "returned_poses", "exit_code",
                                "wall_seconds", "cpu_seconds", "child_peak_rss_bytes", "error")
            if len(source) != 1 or any(source[0].get(key) != row.get(key) for key in preserved_fields):
                raise ValueError("Archived shard does not support the retained actual command")
        if row["status"] not in TERMINAL | INTERRUPTED:
            raise ValueError("Unknown execution state")
        log = directory / row["log_file"]
        if row["status"] in TERMINAL:
            if not log.is_file() or row.get("log_sha256", sha256(log)) != sha256(log):
                raise ValueError("Terminal execution log missing or changed")
            for field in ("wall_seconds", "cpu_seconds", "child_peak_rss_bytes"):
                if field not in row or row[field] < 0:
                    raise ValueError("Terminal execution lacks measured resource use")
        if row["status"] == "ok":
            pose = directory / row["pose_file"]
            if row.get("exit_code") != 0 or sha256(pose) != row.get("pose_sha256"):
                raise ValueError("Successful pose or exit status changed")
            poses = read_poses(pose, reference["names"])
            if len(poses) != row["returned_poses"] or not 1 <= len(poses) <= config["num_modes"]:
                raise ValueError("Successful pose count changed")
            _pose_log(log.read_bytes(), poses)
    return rows


def audit(root):
    """Read-only preflight; never evaluates atomic docking accuracy."""
    from aav_benchmark.blind_docking import planning_inventory
    from aav_benchmark.provenance import sha256
    from aav_benchmark.qualified_pose import Evidence, _all_execution, _check_atomic_sources
    require_stopped(root)
    result, config_path = root / "results/blind_pose", root / "benchmarks/blind_pose_experiment.json"
    config, snapshot = _json(config_path), result / "software/adapters_at_execution_start"
    original = _json(snapshot / "manifest.json")
    for name, digest in original["files"].items():
        if sha256(snapshot / name) != digest:
            raise ValueError("Original execution code snapshot changed")
        if name != "execution_driver.py" and sha256(root / "aav_benchmark" / name) != digest:
            raise ValueError("Frozen runtime adapter changed")
    for package, version in original["versions"].items():
        if importlib.metadata.version(package) != version:
            raise ValueError("Execution software environment changed")
    before, current = _json(result / "planning_inventory.json"), planning_inventory(root)
    for value in (before, current):
        value.pop("candidate_and_ligand_code_sha256", None)
    if before != current:
        raise ValueError("Frozen source-only planning inventory changed")
    # Stop at mechanical validation even when all runs are complete: audit must
    # never fall through to reference fitting or atomic accuracy evaluation.
    evidence = Evidence(root)
    freeze = evidence.json("results/blind_pose/freeze_receipt.json")
    evidence.read("benchmarks/blind_pose_experiment.json", freeze["configuration_sha256"])
    if (config.get("frozen_before_new_docking") is not True or freeze.get("frozen_before_new_docking") is not True
            or freeze.get("new_pose_outcomes_consulted") is not False
            or config["atomic_reference_ledger_sha256"] != freeze["atomic_reference_ledger_sha256"]):
        raise ValueError("Pre-execution freeze changed")
    ledger = evidence.json(config["atomic_reference_ledger"], freeze["atomic_reference_ledger_sha256"])
    inventory = evidence.json("results/blind_pose/planning_inventory.json", freeze["planning_inventory_sha256"])
    oracle = evidence.json("results/blind_pose/oracle_controls.json", freeze["oracle_controls_sha256"])
    _check_atomic_sources(evidence, ledger)
    _, pending, denominator = _all_execution(evidence, config, freeze, inventory, oracle)
    boxes, conditions = [], []
    counts = {"canonical_ok": 0, "canonical_failed": 0, "stale_running": 0,
              "reusable_shard_terminal": 0, "new_seed_jobs": 0}
    for case in config["cases"]:
        for method in config["methods"]:
            cid = case["condition_id"]
            prepared_condition = root / "data/processed/blind_pose" / cid / method
            portable_condition = result / cid / method / "condition_preparation.json"
            if (prepared_condition / "condition.json").read_bytes() != portable_condition.read_bytes():
                raise ValueError("Processed condition differs from its immutable archive")
            condition = _json(portable_condition)
            conditions.append({"prepared": str(prepared_condition), "output": str(result / cid / method), "condition": condition})
            for job in condition["jobs"]:
                if job["status"] != "prepared":
                    continue
                prepared, output = prepared_condition / job["id"], result / cid / method / job["id"]
                for entry in _json(output / "inputs/manifest.json")["files"]:
                    if sha256(prepared / entry["filename"]) != entry["sha256"]:
                        raise ValueError("Processed box differs from its immutable archive")
                rows = validate_rows(root, prepared, output, config_path, canonical=True)
                terminal = {r["seed"]: r for r in rows if r["status"] in TERMINAL}
                counts["canonical_ok"] += sum(r["status"] == "ok" for r in rows)
                counts["canonical_failed"] += sum(r["status"] == "failed" for r in rows)
                counts["stale_running"] += sum(r["status"] in INTERRUPTED for r in rows)
                sources, stale = {}, []
                pattern = f"*/seed-*/results/blind_pose/{cid}/{method}/{job['id']}/runs.json"
                for path in sorted((result / "_recovery/attempts").glob(pattern)):
                    shard_rows = validate_rows(root, prepared, path.parent, config_path, canonical=False)
                    if len(shard_rows) != 1:
                        raise ValueError("Recovery shard must contain exactly one seed")
                    row = shard_rows[0]
                    if row["seed"] in terminal:
                        continue
                    if row["status"] in TERMINAL:
                        if row["seed"] in sources:
                            raise ValueError("Duplicate unmerged terminal seed shards")
                        sources[row["seed"]] = str(path.parent)
                    else:
                        stale.append({"directory": str(path.parent), "row": row})
                missing = [s for s in config["seeds"] if s not in terminal and s not in sources]
                counts["reusable_shard_terminal"] += len(sources)
                counts["new_seed_jobs"] += len(missing)
                boxes.append({"condition_id": cid, "method": method, "box_id": job["id"],
                              "prepared": str(prepared), "output": str(output), "baseline": rows,
                              "sources": sources, "stale_shards": stale, "missing": missing})
    if len(conditions) != 27 or len(boxes) != 52 or len(boxes) * len(config["seeds"]) != 156:
        raise ValueError("Frozen 27-condition/52-box/156-seed denominator changed")
    return {"root": str(root), "configuration_sha256": sha256(config_path), "config": config,
            "conditions": conditions, "boxes": boxes, "counts": counts,
            "atomic_accuracy_evaluated": False, "portable_validation_status": "pending" if pending else "complete",
            "portable_denominator": denominator}


def preserve(path, target):
    from aav_benchmark.provenance import sha256
    if not path.exists():
        return {"source": str(path), "present": False}
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    return {"source": str(path), "archive": str(target), "present": True, "sha256": sha256(target)}


def worker(task):
    """One worker owns one new seed directory; it never writes canonical files."""
    root, prepared, output, seed = task
    sys.path.insert(0, root)
    from aav_benchmark.engine import dock
    return dock(prepared, output, Path(root) / "tools/vina-1.2.7",
                Path(root) / "benchmarks/blind_pose_experiment.json", seeds=[seed])


def assemble_conditions(plan, *, complete):
    """Publish the same receipts as run_condition without rerunning failures."""
    from aav_benchmark.provenance import sha256, write_json
    for item in plan["conditions"]:
        condition, output = item["condition"], Path(item["output"])
        jobs = []
        for job in condition["jobs"]:
            rows = _json(output / job["id"] / "runs.json") if job["status"] == "prepared" and (output / job["id"] / "runs.json").exists() else []
            terminal = len(rows) == len(plan["config"]["seeds"]) and all(r["status"] in TERMINAL for r in rows)
            status = ("preparation_failed" if job["status"] != "prepared" else
                      "ok" if terminal and all(r["status"] == "ok" for r in rows) else
                      "execution_failed" if terminal else "queued")
            jobs.append({**job, "execution_status": status, "runs": rows})
        pending = any(j["execution_status"] == "queued" for j in jobs)
        if complete and pending:
            raise ValueError("Cannot finalize an incomplete condition")
        write_json(output / "condition_runs.json", {"schema_version": 1, "condition_id": condition["condition_id"],
            "complex_id": condition["complex_id"], "method": condition["method"],
            "preparation_status": condition["status"],
            "condition_manifest_sha256": sha256(Path(item["prepared"]) / "condition.json"),
            "configuration_sha256": plan["configuration_sha256"], "eligible_condition": True, "jobs": jobs,
            "status": "running" if pending else "ok" if jobs and all(j["execution_status"] == "ok" for j in jobs) else "workflow_failed",
            "note": "Missing/failed predictor or chemistry remains in eligible denominator; no inferred biological nonbinding"})


def execute(plan, workers):
    from aav_benchmark.provenance import sha256, write_json
    from aav_benchmark.shards import merge_shards
    root, config_path = Path(plan["root"]), Path(plan["root"]) / "benchmarks/blind_pose_experiment.json"
    result = root / "results/blind_pose"
    lock = result / "recovery_driver.lock"
    if lock.exists():
        if _alive(_json(lock)["pid"]):
            raise ValueError("Recovery driver already running")
        # Old lock is archived with the remaining prior infrastructure state below.
    attempt_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + f"-{os.getpid()}"
    attempt = result / "_recovery/attempts" / attempt_id
    attempt.mkdir(parents=True, exist_ok=False)
    infrastructure = []
    for name in ("execution_driver.lock", "execution_progress.json", "execution_driver.log", "recovery_driver.lock"):
        infrastructure.append(preserve(result / name, attempt / "baseline/infrastructure" / name))
    if lock.exists():
        lock.unlink()
    with lock.open("x") as stream:
        json.dump({"pid": os.getpid(), "attempt": str(attempt), "configuration_sha256": plan["configuration_sha256"]}, stream)
    receipt = {"schema_version": 1, "attempt_id": attempt_id, "audit_completed_at_UTC": utcnow(),
        "configuration_sha256": plan["configuration_sha256"], "frozen_settings_changed": False,
        "atomic_accuracy_evaluated": False, "prior_checkpoint": _json(result / "execution_progress.json"),
        "interruption_reason": "unknown; prior driver and Vina workers absent at recovery audit",
        "interruption_time_UTC": None, "interrupted_attempt_costs": "unknown where not measured; never imputed as zero",
        "preflight_counts": plan["counts"], "infrastructure_archives": infrastructure,
        "interrupted_attempts": [], "workers": workers, "maximum_vina_cpu": workers * plan["config"]["cpu"],
        "software": {}, "started_at_UTC": utcnow()}
    sources = [Path(__file__), root / "aav_benchmark/shards.py", root / "aav_benchmark/provenance.py"]
    for source in sources:
        preserve(source, attempt / "software" / source.name)
        receipt["software"][str(source.relative_to(root))] = sha256(source)
    tasks, task_boxes = [], {}
    for index, box in enumerate(plan["boxes"]):
        output = Path(box["output"])
        if box["baseline"]:
            preserve(output / "runs.json", attempt / "baseline" / box["condition_id"] / box["method"] / box["box_id"] / "runs.json")
        for row in box["baseline"]:
            if row["status"] not in INTERRUPTED:
                continue
            saved = {"condition_id": box["condition_id"], "method": box["method"], "box_id": box["box_id"],
                     "row": row, "unknown_interrupt_reason": True, "artifacts": []}
            for name in (row["pose_file"], row["log_file"]):
                saved["artifacts"].append(preserve(output / name, attempt / "baseline" / box["condition_id"] / box["method"] / box["box_id"] / name))
            receipt["interrupted_attempts"].append(saved)
        for stale in box["stale_shards"]:
            receipt["interrupted_attempts"].append({"condition_id": box["condition_id"], "method": box["method"],
                "box_id": box["box_id"], "directory": stale["directory"], "row": stale["row"],
                "unknown_interrupt_reason": True, "evidence_retained_in_original_unique_directory": True})
        for seed in box["missing"]:
            shard = attempt / f"seed-{seed}" / "results/blind_pose" / box["condition_id"] / box["method"] / box["box_id"]
            task = (str(root), box["prepared"], str(shard), seed)
            tasks.append(task)
            task_boxes[task] = index
    write_json(attempt / "recovery_receipt.json", receipt)
    # Long 2VUZ seed jobs occupy workers early; finer seed shards avoid box stragglers.
    tasks.sort(key=lambda t: (Path(t[1]).parts[-3] != "2VUZ", t[1], t[3]))
    merged = set()

    def merge_ready():
        for index, box in enumerate(plan["boxes"]):
            if index in merged or box["missing"]:
                continue
            baseline = box["baseline"]
            if box["sources"]:
                merge_shards(box["prepared"], box["output"], config_path, list(box["sources"].values()),
                             baseline_runs=baseline if baseline else None)
            elif len(baseline) != len(plan["config"]["seeds"]) or any(r["status"] not in TERMINAL for r in baseline):
                raise ValueError("Ready box lacks its complete terminal denominator")
            merged.add(index)

    def checkpoint(status):
        merged_rows, available = [], {}
        for box in plan["boxes"]:
            output = Path(box["output"])
            rows = _json(output / "runs.json") if (output / "runs.json").exists() else []
            merged_rows.extend(rows)
            for row in rows:
                if row["status"] not in TERMINAL:
                    continue  # stale prior attempts are history, not current workers
                available[(box["condition_id"], box["method"], box["box_id"], row["seed"])] = row
            for directory in box["sources"].values():
                for row in _json(Path(directory) / "runs.json"):
                    available[(box["condition_id"], box["method"], box["box_id"], row["seed"])] = row
        for task in tasks:
            path = Path(task[2]) / "runs.json"
            if path.exists():
                box = plan["boxes"][task_boxes[task]]
                for row in _json(path):
                    available[(box["condition_id"], box["method"], box["box_id"], row["seed"])] = row
        rows = list(available.values())
        progress = {"configuration_sha256": plan["configuration_sha256"], "driver_pid": os.getpid(), "workers": workers,
            "scheduled_box_count": 52, "scheduled_seed_runs": 156, "observed_seed_records": len(rows),
            "successful_seed_runs": sum(r["status"] == "ok" for r in rows),
            "failed_seed_runs": sum(r["status"] == "failed" for r in rows),
            "running_seed_runs": sum(r["status"] in INTERRUPTED for r in rows),
            "complete_boxes": len(merged), "canonical_terminal_seed_runs": sum(r["status"] in TERMINAL for r in merged_rows),
            "status": status, "updated_at_UTC": utcnow(), "recovery_attempt": str(attempt.relative_to(root)),
            "previous_interrupted_attempts": len(receipt["interrupted_attempts"]), "atomic_accuracy_evaluated": False}
        write_json(result / "execution_progress.json", progress)
        write_json(attempt / "progress.json", progress)
        print(json.dumps(progress, sort_keys=True), flush=True)

    try:
        merge_ready()
        assemble_conditions(plan, complete=False)
        checkpoint("running")
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(worker, task): task for task in tasks}
            while futures:
                done, _ = wait(futures, timeout=30, return_when=FIRST_COMPLETED)
                for future in done:
                    task = futures.pop(future)
                    rows = future.result()
                    if len(rows) != 1 or rows[0]["seed"] != task[3] or rows[0]["status"] not in TERMINAL:
                        raise ValueError("Worker did not return one terminal frozen-seed run")
                    box = plan["boxes"][task_boxes[task]]
                    # Previous interrupted shard attempts are retained separately and
                    # explicitly attached, without inventing terminal resource values.
                    stale = [s for s in box["stale_shards"] if s["row"]["seed"] == task[3]]
                    if stale:
                        rows[0]["previous_attempts"] = rows[0].get("previous_attempts", []) + [
                            {**s["row"], "status_before_recovery": s["row"]["status"], "status": "interrupted",
                             "interruption_reason": "unknown", "source_directory": s["directory"],
                             "missing_resource_measurements": [k for k in ("wall_seconds", "cpu_seconds", "child_peak_rss_bytes") if k not in s["row"]]}
                            for s in stale]
                        write_json(Path(task[2]) / "runs.json", rows)
                    validate_rows(root, Path(task[1]), Path(task[2]), config_path, canonical=False)
                    box["sources"][task[3]] = task[2]
                    box["missing"].remove(task[3])
                merge_ready()
                assemble_conditions(plan, complete=False)
                checkpoint("running")
        if len(merged) != 52:
            raise ValueError("Not all 52 boxes reached their full frozen denominator")
        assemble_conditions(plan, complete=True)
        checkpoint("complete")
        receipt["completed_at_UTC"], receipt["status"] = utcnow(), "complete"
        write_json(attempt / "recovery_receipt.json", receipt)
    except BaseException as error:
        receipt["stopped_at_UTC"], receipt["status"], receipt["driver_error"] = utcnow(), "interrupted", repr(error)
        write_json(attempt / "recovery_receipt.json", receipt)
        checkpoint("interrupted")
        raise
    finally:
        lock.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--execute", action="store_true", help="Start workers only after the complete read-only audit passes")
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error("At most four two-CPU Vina workers are authorized")
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    plan = audit(root)
    print(json.dumps({"audit": "passed", "read_only": not args.execute, "configuration_sha256": plan["configuration_sha256"],
                      "counts": plan["counts"], "scheduled_seed_runs": 156, "atomic_accuracy_evaluated": False}, sort_keys=True), flush=True)
    if args.execute:
        execute(plan, args.workers)


if __name__ == "__main__":
    main()
