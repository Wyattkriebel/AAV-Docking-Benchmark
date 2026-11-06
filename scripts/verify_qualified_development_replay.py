#!/usr/bin/env python3
"""Replay corrected OLD V2 results from a pinned, freshly exported workspace."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aav_benchmark.development_report_replay import compare_development_report_files
from aav_benchmark.report_replay import load_report_json
from aav_benchmark.replay_workspace import relative_path, restore_workspace, validate_plan

PLAN = "docs/examples/development_v2_workspace_plan.json"
FREEZE = "results/development_pose_report_release/evaluation_v2_freeze.json"
# Filled after the source/data commit and its evaluation freeze exist. The
# freeze must not pin this verifier: that would create a source-hash cycle.
FREEZE_SHA256 = "0d529ac316c3f79263f01aedc9d5f8e11343d5395f91d99a601911e671657863"
RELEASE = "results/development_pose_report_release/generation_v2_receipt.json"
SUMMARY = "results/development_pose_report_release/qualified_pose_v2_summary.json"
MARKDOWN = "docs/QUALIFIED_BLIND_POSE_V2_RESULTS.md"
CLI = "scripts/regenerate_qualified_pose_v2_report.py"
EVALUATOR = "aav_benchmark/qualified_pose_v2.py"
MAX_TOLERANCE = 1e-10
REQUIRED_SOURCES = frozenset({EVALUATOR, CLI, "aav_benchmark/metrics.py",
    "aav_benchmark/development_report_replay.py", "aav_benchmark/report_replay.py",
    "aav_benchmark/replay_workspace.py"})
ARCHIVES = (
    ("results/development_pose_output_foundation/old156_outputs.tar.gz",
     "28ee3e3d2d305f285ef4a6c66f8b6f2075c5fc17e1b164d55a921575de8c1321", 19426016,
     "results/development_pose_output_foundation/old156_outputs_manifest.json",
     "01de653fe57f39d745759323a8c527253a6e3571ee6d5b8efe70d611940f86d0", 313460),
    ("results/development_pose_output_foundation/corrected9_outputs.tar.gz",
     "91e2a67a1a876c1c0e831b7b9d375755d0cb610a496a186893423af5d9608798", 3848957,
     "results/development_pose_output_foundation/corrected9_outputs_manifest.json",
     "b9bf002a405aad206381ff601a86c00e7113dedf04c3d65eb1bb3c5dfecacbae", 23132),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _file(root, relative):
    relative_path(relative)
    path = root / relative
    if any(part.is_symlink() for part in (path, *path.parents) if part != root and root in part.parents):
        raise ValueError("Pinned paths must not contain symbolic links: " + relative)
    if not path.is_file() or not path.resolve().is_relative_to(root):
        raise ValueError("Missing ordinary repository file: " + relative)
    return path


def _json(root, relative):
    return load_report_json(_file(root, relative).read_bytes())


def _pin(row):
    if (type(row) is not dict or set(row) != {"path", "sha256", "bytes"}
            or not isinstance(row["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None
            or type(row["bytes"]) is not int or row["bytes"] < 0):
        raise ValueError("Every file pin requires a relative path, SHA256 and integer byte size")
    relative_path(row["path"])
    return row


def _verify_pin(root, row):
    _pin(row)
    path = _file(root, row["path"])
    if path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
        raise ValueError("Pinned file bytes changed: " + row["path"])


def _sources(root, freeze):
    rows = freeze.get("source_files")
    if type(rows) is not list or not rows:
        raise ValueError("V2 freeze requires an explicit source-file list")
    paths = [_pin(row)["path"] for row in rows]
    if len(set(paths)) != len(paths) or not REQUIRED_SOURCES <= set(paths):
        raise ValueError("V2 freeze omits or duplicates required evaluator/replay source files")
    if Path(__file__).relative_to(ROOT).as_posix() in paths:
        raise ValueError("Evaluation freeze must not pin its own freeze-pinning verifier")
    for row in rows:
        _verify_pin(root, row)
    amendment = _pin(freeze["implementation_amendment"])
    _verify_pin(root, amendment)
    return {row["path"]: row for row in rows}, amendment


def _plan(root):
    plan = _json(root, PLAN)
    validate_plan(plan)
    expected = [{"archive": {"path": a, "sha256": ash, "bytes": abytes},
                 "manifest": {"path": m, "sha256": msh, "bytes": mbytes}, "target_prefix": ""}
                for a, ash, abytes, m, msh, mbytes in ARCHIVES]
    if plan["archives"] != expected:
        raise ValueError("Corrected V2 replay requires exactly the two pinned OLD output archives")
    return plan


def verify(repository, numeric_tolerance=None):
    repository = repository.resolve(strict=True)
    if re.fullmatch(r"[0-9a-f]{64}", FREEZE_SHA256) is None:
        raise ValueError("Corrected V2 evaluation freeze is not pinned; replay is disabled")
    if digest(_file(repository, FREEZE)) != FREEZE_SHA256:
        raise ValueError("Corrected V2 evaluation freeze changed")
    freeze = _json(repository, FREEZE)
    if (type(freeze) is not dict or freeze.get("outputs") != {"summary": SUMMARY, "markdown": MARKDOWN}
            or type(freeze.get("planned_cross_platform_absolute_derived_tolerance")) is not float
            or freeze["planned_cross_platform_absolute_derived_tolerance"] != MAX_TOLERANCE):
        raise ValueError("Corrected V2 output roles or frozen tolerance changed")
    if numeric_tolerance is not None and (
            type(numeric_tolerance) not in (int, float) or not 0 <= numeric_tolerance <= MAX_TOLERANCE
            or not math.isfinite(numeric_tolerance)):
        raise ValueError("Explicit tolerance must be finite and within the frozen 1e-10 bound")
    sources, amendment = _sources(repository, freeze)
    plan = _plan(repository)
    release = _json(repository, RELEASE)
    if (type(release) is not dict or type(release.get("process_exit_code")) is not int
            or release["process_exit_code"] != 0 or release.get("atomic_accuracy_generated") is not True
            or release.get("source_commit") != plan["commit"]
            or release.get("workspace_plan_sha256") != digest(_file(repository, PLAN))):
        raise ValueError("Canonical V2 report release does not match the replay plan")
    outputs = release.get("outputs")
    if (type(outputs) is not list or len(outputs) != 2
            or {_pin(row)["path"] for row in outputs} != {SUMMARY, MARKDOWN}):
        raise ValueError("Canonical V2 development report roles changed")
    for row in outputs:
        _verify_pin(repository, row)
    summary_pin = next(row for row in outputs if row["path"] == SUMMARY)
    with tempfile.TemporaryDirectory(prefix="aav-qualified-development-v2-replay-") as temporary:
        workspace = (Path(temporary) / "workspace").resolve()
        provenance = restore_workspace(repository, repository / PLAN, workspace)
        if any(path.exists() or path.is_symlink() for path in (workspace / "tools", workspace / "data/processed")):
            raise ValueError("Replay workspace must contain no installed tools or processed tree")
        if digest(_file(workspace, FREEZE)) != FREEZE_SHA256:
            raise ValueError("Exported V2 evaluation freeze changed")
        _sources(workspace, freeze)
        if any(path.exists() or path.is_symlink() for path in (workspace / SUMMARY, workspace / MARKDOWN)):
            raise ValueError("Replay source commit already contains a derived V2 development report")
        if provenance["commit"] != plan["commit"]:
            raise ValueError("Replay export does not match the frozen source commit")
        env = os.environ.copy()
        env.update({name: "1" for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
        run = subprocess.run(
            [sys.executable, "-I", str(workspace / CLI), "--root", str(workspace),
             "--implementation-amendment", amendment["path"],
             "--implementation-amendment-sha256", amendment["sha256"]],
            cwd=workspace, env=env, capture_output=True, text=True,
        )
        if run.returncode:
            raise ValueError("Fresh corrected V2 development report failed: " + run.stderr.strip())
        for row in outputs:
            _verify_pin(repository, row)
        comparison = compare_development_report_files(
            repository / SUMMARY, _file(workspace, SUMMARY), numeric_tolerance=numeric_tolerance)
        if comparison["expected_file"]["sha256"] != summary_pin["sha256"]:
            raise ValueError("Compared canonical V2 JSON differs from its generation receipt")
        if (repository / MARKDOWN).read_bytes() != _file(workspace, MARKDOWN).read_bytes():
            raise ValueError("Replayed corrected V2 development Markdown differs")
        for row in outputs:
            _verify_pin(repository, row)
        return {
            "schema_version": 1, "status": "verified_qualified_development_v2_replay",
            "source_commit": provenance["commit"], "workspace_plan_sha256": digest(repository / PLAN),
            "evaluation_freeze_sha256": FREEZE_SHA256, "evaluator_source_sha256": sources[EVALUATOR]["sha256"],
            "implementation_amendment": amendment,
            "exported_git_files": provenance["exported_files"],
            "verified_archives": len(plan["archives"]),
            "verified_archive_members": provenance["verified_archive_members"],
            "native_executables_used": False, "tools_absent": True, "processed_tree_absent": True,
            "comparison": comparison, "markdown_bytes_identical": True,
            "canonical_summary_sha256": summary_pin["sha256"],
            "recomputed_summary_sha256": digest(workspace / SUMMARY),
            "new_six_case_structural_metrics_computed": False,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=ROOT)
    parser.add_argument("--numeric-tolerance", type=float,
                        help="Explicit named-derived-field tolerance, at most 1e-10; default requires byte equality")
    parser.add_argument("--output", type=Path, help="Write a new verification receipt without replacing an existing file")
    args = parser.parse_args()
    try:
        result = verify(args.repository, args.numeric_tolerance)
        raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as stream:
                stream.write(raw)
        print(raw, end="")
    except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Qualified development V2 replay failed: " + str(error) + "\n")


if __name__ == "__main__":
    main()
