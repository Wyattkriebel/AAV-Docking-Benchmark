"""Replay commands must keep retained evidence immutable and use frozen inputs."""
import importlib.util
import json
from pathlib import Path
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/run_next_experiments.py"
spec = importlib.util.spec_from_file_location("run_next_experiments", MODULE_PATH)
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def _source(tmp_path):
    source = tmp_path / "source"
    (source / "benchmarks").mkdir(parents=True)
    for arm in ("site", "gal9"):
        (source / "benchmarks" / f"{arm}_experiment.json").write_text('{"seeds": [1]}\n')
    return source


def test_replay_rejects_repository_and_resolved_symlink_destination(tmp_path):
    source = _source(tmp_path)
    with pytest.raises(ValueError, match="outside"):
        replay._fresh_workspace(source, source / "results/new", "site")
    link = tmp_path / "linked_source"
    link.symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError, match="outside"):
        replay._fresh_workspace(source, link / "data/replay", "site")
    assert not (source / "results").exists()


def test_replay_preparation_never_overwrites_prior_work(tmp_path):
    source = _source(tmp_path)
    destination = tmp_path / "replay"
    destination.mkdir()
    sentinel = destination / "prior_evidence"
    sentinel.write_bytes(b"original")
    with pytest.raises(FileExistsError, match="new or empty"):
        replay._fresh_workspace(source, destination, "site")
    assert sentinel.read_bytes() == b"original"
    assert list(destination.iterdir()) == [sentinel]


def test_workspace_checks_arm_and_frozen_plan_on_every_stage(tmp_path):
    source = _source(tmp_path)
    work = replay._fresh_workspace(source, tmp_path / "replay", "site")
    assert replay._workspace(source, work, "site") == work
    marker = json.loads((work / "replay.json").read_text())
    assert "not newly prospective" in marker["role"]
    with pytest.raises(ValueError, match="arm or frozen configuration"):
        replay._workspace(source, work, "gal9")
    (work / "benchmarks/site_experiment.json").write_text('{"seeds": [2]}\n')
    with pytest.raises(ValueError, match="arm or frozen configuration"):
        replay._workspace(source, work, "site")
    assert json.loads((source / "benchmarks/site_experiment.json").read_text()) == {"seeds": [1]}


def test_source_restoration_rejects_mutated_archived_bytes(tmp_path):
    import gzip
    archive = tmp_path / "source.gz"
    archive.write_bytes(gzip.compress(b"changed"))
    with pytest.raises(ValueError, match="checksum differs"):
        replay._restore_source(archive, tmp_path / "raw/source", "0" * 64)


def test_duplicate_method_request_rejected_before_execution(tmp_path, monkeypatch):
    def never_run(*args, **kwargs):
        pytest.fail("duplicate method request reached scientific execution")
    monkeypatch.setattr(replay, "run_site", never_run)
    with pytest.raises(SystemExit) as error:
        replay.main(["site-run", "--work-dir", str(tmp_path), "--method", "fpocket", "--method", "fpocket"])
    assert error.value.code == 1


def test_changed_method_lock_rejected_before_scientific_execution(tmp_path):
    source = _source(tmp_path)
    method = source / "benchmarks/method.lock.json"
    method.write_text('{"fixed": true}\n')
    config = source / "benchmarks/site_experiment.json"
    config.write_text(json.dumps({"method_locks": {"benchmarks/method.lock.json": replay.sha256(method)}}))
    work = replay._fresh_workspace(source, tmp_path / "replay", "site")
    assert replay._workspace(source, work, "site") == work
    (work / "benchmarks/method.lock.json").write_text('{"fixed": false}\n')
    with pytest.raises(ValueError, match="method lock differs"):
        replay._workspace(source, work, "site")
