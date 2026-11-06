"""Archive-only guided-report preflight; actual native accuracy is forbidden."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

import aav_benchmark.haddock_pose_report as guided
from aav_benchmark.provenance import sha256, write_json


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("paired_public_fixture", ROOT / "tests/test_paired_method_public_portability.py")
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)


def test_guided_report_source_archives_need_no_installation_or_processed_tree(tmp_path, monkeypatch):
    omitted, _ = fixture._copy_source_tree(tmp_path)
    def copy_file(path):
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / path, target)
    for path in ("prepared_plan.json", "preparation_receipt.json", "prepared_inputs_manifest.json", "prepared_inputs.tar.gz",
                 "native_software_manifest.json", "native_score_weights.json", "execution_freeze.json"):
        copy_file("results/haddock_paired_cohort/" + path)
    prepared = json.loads((ROOT / "results/haddock_paired_cohort/prepared_plan.json").read_text())
    for key in ("proposal", "driver", "assembly_protocol_policy", "source_manifest", "preflight_receipt", "release_lock"):
        copy_file(prepared[key]["path"])
    reference_plan = json.loads((ROOT / "benchmarks/structural_paired_method_evaluation.json").read_text())
    for key in ("reference_evaluation_plan", "reference_evaluation_freeze"):
        copy_file(reference_plan[key]["path"])
    assert all(not (tmp_path / path).exists() for path in omitted)
    assert not (tmp_path / "tools").exists() and not (tmp_path / "data/processed").exists()
    evaluation = guided.example_plan(tmp_path)
    evaluation.pop("orchestration_recovery")  # No native outcome files are copied into this source-only fixture.
    plan_path, freeze_path = "benchmarks/guided_test_evaluation.json", "results/guided_test_evaluation_freeze.json"
    write_json(tmp_path / plan_path, evaluation)
    write_json(tmp_path / freeze_path, {"evaluation_plan_sha256": sha256(tmp_path / plan_path),
        "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    def forbidden(*args, **kwargs):
        pytest.fail("Source-only pending validation reached actual guided pose accuracy")
    monkeypatch.setattr(guided, "_source_context", forbidden)
    monkeypatch.setattr(guided, "score_pool", forbidden)
    result = guided.generate_report(tmp_path, plan_path, freeze_path, allow_pending=True)
    assert not result["metrics_eligible"] and result["status"] == "pending_incomplete_workflow"
    assert result["counts"]["eligible_cells"] == 108 and result["counts"]["native_ready_cells"] == 39
    assert result["counts"]["terminal_cells"] == 0 and len(result["pending"]) == 108
    assert result["native_provenance"]["recorded_installed_source_files"] == 370
    assert result["native_provenance"]["local_installation_or_binary_reverified"] is False
    assert result["reference_validation"]["retained_archive_records_verified"] == 163
    assert not any(p.startswith(("tools/", "data/processed/")) for p in result["source_sha256_inventory"])
    # The copied package itself also imports/runs with no HADDOCK installation.
    script = """
import json
from pathlib import Path
import aav_benchmark.haddock_pose_report as report
assert Path(report.__file__).resolve().is_relative_to(Path.cwd())
def forbidden(*args, **kwargs):
    raise AssertionError('Atomic accuracy entered before all108 cells terminal')
report._source_context = forbidden
report.score_pool = forbidden
result = report.generate_report(Path.cwd(), 'benchmarks/guided_test_evaluation.json', 'results/guided_test_evaluation_freeze.json', allow_pending=True)
print(json.dumps({'status': result['status'], 'ready': result['counts']['native_ready_cells'], 'tools': Path('tools').exists()}))
"""
    replay = json.loads(subprocess.check_output([sys.executable, "-c", script], cwd=tmp_path, text=True))
    assert replay == {"status": "pending_incomplete_workflow", "ready": 39, "tools": False}
    with pytest.raises(ValueError, match="remains incomplete"):
        guided.generate_report(tmp_path, plan_path, freeze_path)
