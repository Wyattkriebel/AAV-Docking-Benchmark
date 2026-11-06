"""Portable union-workflow preflight; no native atomic metrics are permitted."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from aav_benchmark.provenance import sha256, write_json


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "union_paired_public_fixture", ROOT / "tests/test_paired_method_public_portability.py")
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)


def _copy_union_sources(destination):
    """Use only declared source/archive bytes; native outcomes remain absent."""
    omitted, _ = fixture._copy_source_tree(destination)
    return omitted


METRIC_NAMES = {
    "generate_paired_report", "_source_context", "score_pool", "_score_pool",
    "protein_fits", "_protein_alignment", "site_summary", "accuracy_summary",
    "native_metric_controls", "paired_statistics", "site_statistics",
    "comparison_statistics", "conditional_statistics", "macro_summaries",
    "profile_sensitivity_statistics", "_statistics", "paired_family_bootstrap",
    "aggregate_runs", "fit_receptor", "ligand_rmsd", "contact_pairs",
    "contact_metrics", "clash_metrics", "summarize_poses", "site_metrics",
}


def _block_metrics(monkeypatch):
    """Patch aliases by identity, including from-import aliases in the adapter."""
    import aav_benchmark.heldout_pose_report
    import aav_benchmark.metrics
    import aav_benchmark.paired_method_report
    import aav_benchmark.qualified_pose
    import aav_benchmark.statistics
    modules = [module for name, module in tuple(sys.modules.items())
               if name.startswith("aav_benchmark.") and module is not None]
    targets = {id(value) for module in modules for name, value in vars(module).items()
               if name in METRIC_NAMES and callable(value)}
    def forbidden(*args, **kwargs):
        pytest.fail("Portable union pending workflow entered a main report or atomic metric")
    for module in modules:
        for name, value in tuple(vars(module).items()):
            if callable(value) and id(value) in targets:
                monkeypatch.setattr(module, name, forbidden)


def _assert_pending(result):
    assert result["status"] in {"pending", "pending_incomplete_workflow"}
    assert result.get("metrics_eligible", False) is False
    counts = result["counts"]
    assert counts["eligible_seed_pools"] == 108
    assert counts["eligible_source_method_slots"] == 216
    assert counts["terminal_source_method_slots"] == 0
    assert counts["recorded_pose_count"] == 0


def test_union_public_source_only_clone_blocks_every_metric_and_preserves_live_index(tmp_path, monkeypatch):
    from aav_benchmark import paired_common_pose_report as report
    omitted = _copy_union_sources(tmp_path)
    assert len(omitted) == 11 and all(not (tmp_path / path).exists() for path in omitted)
    assert not (tmp_path / "tools").exists()
    assert not (tmp_path / "data/processed").exists()
    _block_metrics(monkeypatch)
    evaluation = report.example_plan(tmp_path)
    evaluation["frozen_before_union_pose_accuracy"] = True
    plan_path = "benchmarks/union_portability_evaluation.json"
    freeze_path = "results/union_portability_evaluation_freeze.json"
    write_json(tmp_path / plan_path, evaluation)
    write_json(tmp_path / freeze_path, {
        "evaluation_plan_sha256": sha256(tmp_path / plan_path),
        "frozen_before_test_pose_accuracy": True,
        "frozen_before_union_pose_accuracy": True,
        "new_test_pose_accuracy_consulted": False,
    })
    union_output = tmp_path / "results/structural_pose_v3/common_pose_union"
    before = {path.relative_to(union_output).as_posix(): path.read_bytes()
              for path in union_output.rglob("*") if path.is_file()}
    original_open = Path.open
    def guarded_open(path, *args, **kwargs):
        resolved = path.resolve()
        if resolved.is_relative_to(tmp_path / "tools") or resolved.is_relative_to(tmp_path / "data/processed"):
            pytest.fail("Portable union reader attempted a native or ignored processed dependency")
        mode = args[0] if args else kwargs.get("mode", "r")
        if resolved.is_relative_to(union_output) and any(flag in mode for flag in "wax+"):
            pytest.fail("Portable union report attempted to mutate its live source/score output")
        return original_open(path, *args, **kwargs)
    monkeypatch.setattr(Path, "open", guarded_open)
    probe = report.probe_workflow(tmp_path, plan_path, freeze_path, require_freeze=True)
    _assert_pending(probe)
    pending = report.generate_report(tmp_path, plan_path, freeze_path, allow_pending=True)
    _assert_pending(pending)
    assert before == {path.relative_to(union_output).as_posix(): path.read_bytes()
                      for path in union_output.rglob("*") if path.is_file()}
    assert not (union_output / "source_pool.json").exists()
    assert not (union_output / "scores_manifest.json").exists()
    assert not (union_output / "scores").exists()
    assert not (tmp_path / "tools").exists()
    assert not (tmp_path / "data/processed").exists()
    with pytest.raises(ValueError, match="incomplete|pending|terminal"):
        report.generate_report(tmp_path, plan_path, freeze_path)
    # Pending source workflows never waive the separate union accuracy freeze.
    original_plan = (tmp_path / plan_path).read_bytes()
    unfrozen = json.loads(original_plan)
    unfrozen["frozen_before_union_pose_accuracy"] = False
    write_json(tmp_path / plan_path, unfrozen)
    try:
        with pytest.raises(ValueError, match="pre-accuracy|freeze"):
            report.probe_workflow(tmp_path, plan_path, freeze_path, require_freeze=True)
        with pytest.raises(ValueError, match="pre-accuracy|freeze"):
            report.generate_report(tmp_path, plan_path, freeze_path, allow_pending=True)
    finally:
        (tmp_path / plan_path).write_bytes(original_plan)
    assert (tmp_path / plan_path).read_bytes() == original_plan

    script = r"""
import json
from pathlib import Path
import sys
import aav_benchmark.paired_common_pose_report as report
import aav_benchmark.heldout_pose_report
import aav_benchmark.metrics
import aav_benchmark.paired_method_report
import aav_benchmark.qualified_pose
import aav_benchmark.statistics
assert Path(report.__file__).resolve().is_relative_to(Path.cwd())
names = set(METRIC_NAMES_LITERAL)
modules = [module for name, module in tuple(sys.modules.items()) if name.startswith('aav_benchmark.') and module is not None]
targets = {id(value) for module in modules for name, value in vars(module).items() if name in names and callable(value)}
def forbidden(*args, **kwargs):
    raise AssertionError('Atomic metric reached from source-only union workflow')
for module in modules:
    for name, value in tuple(vars(module).items()):
        if callable(value) and id(value) in targets:
            setattr(module, name, forbidden)
original_open = Path.open
def guarded_open(path, *args, **kwargs):
    resolved = path.resolve()
    if resolved.is_relative_to(Path.cwd() / 'tools') or resolved.is_relative_to(Path.cwd() / 'data/processed'):
        raise AssertionError('Ignored native/processed dependency reached')
    mode = args[0] if args else kwargs.get('mode', 'r')
    if resolved.is_relative_to(Path.cwd() / 'results/structural_pose_v3/common_pose_union') and any(flag in mode for flag in 'wax+'):
        raise AssertionError('Live union source/score output mutated by report')
    return original_open(path, *args, **kwargs)
Path.open = guarded_open
plan = 'benchmarks/union_portability_evaluation.json'
freeze = 'results/union_portability_evaluation_freeze.json'
probe = report.probe_workflow(Path.cwd(), plan, freeze, require_freeze=True)
pending = report.generate_report(Path.cwd(), plan, freeze, allow_pending=True)
print(json.dumps({'probe_status': probe['status'], 'report_status': pending['status'],
                  'probe_terminal': probe['counts']['terminal_source_method_slots'],
                  'report_terminal': pending['counts']['terminal_source_method_slots'],
                  'recorded': pending['counts']['recorded_pose_count'],
                  'tools_present': Path('tools').exists(),
                  'processed_present': Path('data/processed').exists(),
                  'live_index_present': Path('results/structural_pose_v3/common_pose_union/source_pool.json').exists()}))
""".replace("METRIC_NAMES_LITERAL", repr(sorted(METRIC_NAMES)))
    replay = json.loads(subprocess.check_output([sys.executable, "-c", script], cwd=tmp_path, text=True))
    assert replay == {"probe_status": probe["status"], "report_status": pending["status"],
                      "probe_terminal": 0, "report_terminal": 0, "recorded": 0,
                      "tools_present": False, "processed_present": False, "live_index_present": False}
    # Portability does not waive implementation/source pins.
    target = tmp_path / "aav_benchmark/paired_method_report.py"
    target.write_bytes(target.read_bytes() + b"\n# source tamper\n")
    with pytest.raises(ValueError, match="changed|checksum|pin|source|implementation"):
        report.probe_workflow(tmp_path, plan_path, freeze_path, require_freeze=True)
