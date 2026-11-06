#!/usr/bin/env python3
"""Package completed development trees without invoking an empirical reader."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import tarfile


ROOT = Path(__file__).resolve().parents[2]
GROUPS = {
    "old156": "results/blind_pose",
    "corrected9": "results/blind_pose_neighborhood8",
}
EXCLUDED_DERIVED_PATHS = {"results/blind_pose/qualified_pose_summary.json"}
SEEDS = {104729, 130363, 155921}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def regular(root, name):
    path = root / name
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise ValueError("Symlink source: " + str(path))
    if not stat.S_ISREG(path.stat().st_mode):
        raise ValueError("Not an ordinary file: " + str(path))
    return path


def snapshot(root, relative):
    source = root / relative
    if not source.is_dir() or source.is_symlink():
        raise ValueError("Missing ordinary source directory: " + relative)
    records = []
    for directory, dirs, files in os.walk(source, followlinks=False):
        for name in dirs:
            path = Path(directory) / name
            if path.is_symlink() or not stat.S_ISDIR(path.stat().st_mode):
                raise ValueError("Nonordinary source directory: " + str(path))
        for name in files:
            path = Path(directory) / name
            relative_name = path.relative_to(root).as_posix()
            if relative_name in EXCLUDED_DERIVED_PATHS:
                continue
            raw = regular(root, relative_name).read_bytes()
            records.append({"path": relative_name, "bytes": len(raw), "sha256": digest(raw)})
    records.sort(key=lambda row: row["path"])
    folded = [row["path"].casefold() for row in records]
    if not records or len(set(folded)) != len(folded):
        raise ValueError("Empty or case-colliding source inventory")
    return records


def terminal_gates(root):
    progress = json.loads(regular(root, "results/blind_pose/execution_progress.json").read_bytes())
    if (progress.get("status") != "complete"
            or progress.get("canonical_terminal_seed_runs") != 156
            or progress.get("running_seed_runs") != 0
            or progress.get("successful_seed_runs", -1) + progress.get("failed_seed_runs", -1) != 156):
        raise ValueError("The original156 workflow is not terminal")
    validation_path = regular(root, "results/blind_pose/report_release/execution_validation.json")
    validation = json.loads(validation_path.read_bytes())
    if (validation.get("atomic_accuracy_computed") is not False
            or validation.get("pending") != []
            or validation.get("evaluator_sha256") != "c77333c8fbda7fe81f79fbe99211a65f3c91f51382daaa80f2fb89bd12dfc525"
            or validation["denominator"]["terminal_seed_runs"] != 156):
        raise ValueError("Missing workflow-only raw-output validation")
    for name, pin in validation["source_sha256"].items():
        raw = regular(root, name).read_bytes()
        if len(raw) != pin["bytes"] or digest(raw) != pin["sha256"]:
            raise ValueError("Validated raw evidence changed: " + name)
    manifest_path = regular(root, "results/blind_pose_neighborhood8/execution_manifest.json")
    corrected = json.loads(manifest_path.read_bytes())
    boxes = {"aav9_apo/oracle/box01", "aav9_apo/p2rank/box03", "aav9_holo/p2rank/box05"}
    if len(corrected["jobs"]) != 3 or {row["box_key"] for row in corrected["jobs"]} != boxes:
        raise ValueError("Corrected control box denominator changed")
    for job in corrected["jobs"]:
        if len(job["runs"]) != 3 or {row["seed"] for row in job["runs"]} != SEEDS:
            raise ValueError("Corrected control seed denominator changed")
        for row in job["runs"]:
            if row["status"] not in {"ok", "failed"}:
                raise ValueError("Corrected control run is not terminal")
            base = "results/blind_pose_neighborhood8/" + job["box_key"] + "/"
            for key in ("log_file", "pose_file"):
                if Path(row[key]).name != row[key] or row[key] in {".", ".."}:
                    raise ValueError("Unsafe control output name")
            regular(root, base + row["log_file"])
            if row["status"] == "ok":
                if row.get("exit_code") != 0 or digest(regular(root, base + row["pose_file"]).read_bytes()) != row["pose_sha256"]:
                    raise ValueError("Corrected control pose provenance changed")
    return {
        "old156": {key: progress[key] for key in ("status", "canonical_terminal_seed_runs", "successful_seed_runs", "failed_seed_runs", "running_seed_runs", "previous_interrupted_attempts")},
        "corrected9_terminal_runs": 9,
        "workflow_validation_sha256": digest(validation_path.read_bytes()),
        "corrected_execution_manifest_sha256": digest(manifest_path.read_bytes()),
    }


def package(root, destination):
    root = root.resolve()
    destination = destination.resolve()
    if any(destination.is_relative_to(root / value) for value in GROUPS.values()):
        raise ValueError("Archive destination must be outside both source trees")
    gates = terminal_gates(root)
    inventories = {name: snapshot(root, path) for name, path in GROUPS.items()}
    destination.mkdir(parents=True, exist_ok=True)
    created, output_records = [], []
    try:
        for name, records in inventories.items():
            target = destination / (name + "_outputs.tar.gz")
            with target.open("xb") as output:
                created.append(target)
                with gzip.GzipFile(fileobj=output, filename="", mode="wb", mtime=0, compresslevel=9) as compressed:
                    with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                        for row in records:
                            raw = regular(root, row["path"]).read_bytes()
                            if len(raw) != row["bytes"] or digest(raw) != row["sha256"]:
                                raise ValueError("Source changed during packaging: " + row["path"])
                            info = tarfile.TarInfo(row["path"])
                            info.size, info.mode, info.uid, info.gid, info.mtime = len(raw), 0o644, 0, 0, 0
                            info.uname = info.gname = ""
                            archive.addfile(info, io.BytesIO(raw))
            if target.stat().st_size >= 100_000_000:
                raise ValueError("Archive needs a separately reviewed split below100MB")
            archive_pin = {"path": target.name, "bytes": target.stat().st_size, "sha256": digest(target.read_bytes())}
            manifest = {"schema_version": 1, "archive": target.name, "archive_sha256": archive_pin["sha256"],
                        "target_prefix": "", "source_root": GROUPS[name], "files": records,
                        "excluded_derived_paths": sorted(EXCLUDED_DERIVED_PATHS),
                        "complete_ordinary_file_tree": True, "numeric_payloads_interpreted": False,
                        "source_snapshot_sha256": digest(encoded(records))}
            path = destination / (name + "_outputs_manifest.json")
            with path.open("xb") as stream:
                created.append(path)
                stream.write(encoded(manifest))
            output_records.append({"archive": archive_pin, "manifest": {"path": path.name, "bytes": path.stat().st_size, "sha256": digest(path.read_bytes())},
                                   "files": len(records), "uncompressed_bytes": sum(row["bytes"] for row in records)})
        if terminal_gates(root) != gates or any(snapshot(root, GROUPS[name]) != rows for name, rows in inventories.items()):
            raise ValueError("Source trees or terminal gates changed during packaging")
        receipt = {"schema_version": 1, "stage": "Post-execution packaging; not a scientific or contemporaneous producer freeze.",
                   "terminal_gates": gates, "archives": output_records, "all_before_after_source_snapshots_equal": True,
                   "all_failures_partials_recovery_and_interrupted_attempts_retained": True,
                   "resource_measurements_not_invented_or_filled": True, "atomic_accuracy_computed": False,
                   "packager_sha256": digest(Path(__file__).read_bytes())}
        path = destination / "packaging_receipt.json"
        with path.open("xb") as stream:
            created.append(path)
            stream.write(encoded(receipt))
        return receipt
    except BaseException:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(package(args.root, args.destination), indent=2, sort_keys=True))
