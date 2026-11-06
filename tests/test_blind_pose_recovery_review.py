"""Independent controls for the recovery driver's non-execution contracts."""

import importlib.util
import json
from pathlib import Path
import sys

import pytest

from aav_benchmark.provenance import sha256, write_json


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/resume_blind_pose.py"
SPEC = importlib.util.spec_from_file_location("blind_pose_recovery_review", SCRIPT)
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


def condition_plan(tmp_path, rows):
    prepared, output = tmp_path / "prepared", tmp_path / "output"
    condition = {"condition_id": "case", "complex_id": "case", "method": "fpocket",
                 "status": "prepared", "jobs": [{"id": "box01", "status": "prepared"}]}
    write_json(prepared / "condition.json", condition)
    write_json(output / "box01/runs.json", rows)
    return {"configuration_sha256": "a" * 64, "config": {"seeds": [11, 13, 17]},
            "conditions": [{"prepared": str(prepared), "output": str(output), "condition": condition}]}, output


def test_terminal_failure_is_retained_without_rerunning_or_rewriting_history(tmp_path, monkeypatch):
    import aav_benchmark.engine
    monkeypatch.setattr(aav_benchmark.engine, "dock", lambda *a, **k: pytest.fail("Assembly invoked docking"))
    prior = {"seed": 13, "status": "running", "command": ["historical-vina", "--seed", "13"]}
    rows = [{"seed": 11, "status": "ok", "wall_seconds": 2.0},
            {"seed": 13, "status": "failed", "error": "Measured terminal failure",
             "wall_seconds": 4.0, "previous_attempts": [prior]},
            {"seed": 17, "status": "ok", "wall_seconds": 3.0}]
    plan, output = condition_plan(tmp_path, rows)
    manifest = output / "box01/runs.json"
    original = manifest.read_bytes()
    recovery.assemble_conditions(plan, complete=True)
    receipt = json.loads((output / "condition_runs.json").read_text())
    assert receipt["status"] == "workflow_failed"
    assert receipt["jobs"][0]["execution_status"] == "execution_failed"
    assert receipt["jobs"][0]["runs"] == rows
    assert manifest.read_bytes() == original
    assert "wall_seconds" not in receipt["jobs"][0]["runs"][1]["previous_attempts"][0]
    assert "exit_code" not in receipt["jobs"][0]["runs"][1]["previous_attempts"][0]


def test_partial_seed_denominator_cannot_be_finalized(tmp_path):
    plan, output = condition_plan(tmp_path, [{"seed": 11, "status": "ok"},
                                             {"seed": 13, "status": "running"}])
    recovery.assemble_conditions(plan, complete=False)
    receipt = json.loads((output / "condition_runs.json").read_text())
    assert receipt["status"] == "running"
    assert receipt["jobs"][0]["execution_status"] == "queued"
    original = (output / "condition_runs.json").read_bytes()
    with pytest.raises(ValueError, match="incomplete condition"):
        recovery.assemble_conditions(plan, complete=True)
    assert (output / "condition_runs.json").read_bytes() == original


def test_default_cli_is_check_only_even_when_missing_jobs_exist(tmp_path, monkeypatch, capsys):
    audit_roots = []
    def fake_audit(root):
        audit_roots.append(root)
        return {"configuration_sha256": "a" * 64, "counts": {"canonical_ok": 64,
                "canonical_failed": 0, "stale_running": 4, "new_seed_jobs": 92}}
    monkeypatch.setattr(recovery, "audit", fake_audit)
    monkeypatch.setattr(recovery, "execute", lambda *a, **k: pytest.fail("Default CLI executed workers"))
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--root", str(tmp_path)])
    recovery.main()
    result = json.loads(capsys.readouterr().out)
    assert audit_roots == [tmp_path.resolve()]
    assert result["read_only"] is True
    assert result["scheduled_seed_runs"] == 156
    assert result["atomic_accuracy_evaluated"] is False
    assert result["counts"]["new_seed_jobs"] == 92


def test_preservation_copies_exact_bytes_and_does_not_invent_missing_artifacts(tmp_path):
    original = tmp_path / "interrupted.log"
    original.write_bytes(b"partial Vina output\nwithout measured termination\n")
    destination = tmp_path / "archive/interrupted.log"
    receipt = recovery.preserve(original, destination)
    assert original.read_bytes() == destination.read_bytes()
    assert receipt["sha256"] == sha256(original)
    missing = recovery.preserve(tmp_path / "absent.pdbqt", tmp_path / "archive/absent.pdbqt")
    assert missing["present"] is False
    assert not (tmp_path / "archive/absent.pdbqt").exists()
    assert "sha256" not in missing
