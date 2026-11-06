"""Synthetic work-sharing races and frozen-engine guards; no native launch."""
from pathlib import Path
from types import SimpleNamespace
import json

import pytest

from aav_benchmark import engine
from aav_benchmark.paired_worksharing_model import (
    LogFrontier, PARENT_EVENTS, WORKER_EVENTS, SchedulingStop, SyntheticWorkerBudget,
    interleavings, parse_frontier, postclaim_decision, reconcile_schedule, simulate,
)


def test_postclaim_barrier_prevents_duplicates_in_all70_interleavings():
    schedules = list(interleavings())
    assert len(schedules) == 70
    results = [simulate(schedule) for schedule in schedules]
    assert not any(row["duplicate_search"] for row in results)
    assert any(row["parent_search"] for row in results)
    assert any(row["worker_search"] for row in results)
    assert not any(row["any_process_paused"] for row in results)


def test_exclusive_manifest_alone_has_a_real_check_create_race():
    schedule = (PARENT_EVENTS[0], PARENT_EVENTS[1], *WORKER_EVENTS,
                PARENT_EVENTS[2], PARENT_EVENTS[3])
    assert simulate(schedule, postclaim_barrier=False)["duplicate_search"] is True
    assert simulate(schedule)["duplicate_search"] is False
    with pytest.raises(SchedulingStop, match="eight declared atomic events"):
        simulate((*schedule, "unmodeled_external_writer"))


def test_frontier_ignores_incomplete_print_but_rejects_restarts_skips_and_unknown_ids():
    ids = ["a", "b", "c"]
    result = parse_frontier(b"START_MODERN a\nnoise\nSTART_MODERN b", ids, (1, 2))
    assert result.index == 0
    for raw in [b"START_MODERN b\n", b"START_MODERN a\nSTART_MODERN a\n", b"START_MODERN unknown\n"]:
        with pytest.raises(SchedulingStop):
            parse_frontier(raw, ids, (1, 2))


def decision(**changes):
    args = dict(before=LogFrontier((1, 2), 100, (0,)),
                after=LogFrontier((1, 2), 110, (0,)), index=35,
                prepared=True, exclusive_creation=True, owns_empty=True,
                verified_regular_stdout_binding=True)
    args.update(changes)
    return postclaim_decision(**args)


def test_only_future_prepared_owned_empty_condition_can_pass_the_model():
    assert decision() == "barrier_holds_in_model"
    for change in [dict(index=0), dict(prepared=False), dict(exclusive_creation=False),
                   dict(owns_empty=False), dict(verified_regular_stdout_binding=False)]:
        with pytest.raises(SchedulingStop):
            decision(**change)


def test_caught_up_parent_aborts_without_removing_reservation_or_launching():
    after = LogFrontier((1, 2), 200, tuple(range(36)))
    assert decision(after=after) == "abort_and_retain_reservation"


@pytest.mark.parametrize("after", [LogFrontier((9, 9), 110, (0,)),
                                   LogFrontier((1, 2), 99, (0,)),
                                   LogFrontier((1, 2), 110, ())])
def test_log_rebinding_truncation_and_frontier_regression_stop(after):
    with pytest.raises(SchedulingStop):
        decision(after=after)


def test_complete_denominator_survives_legacy_exit_before_unstarted_work():
    ids = [f"condition{i:02}" for i in range(36)]
    report = reconcile_schedule(ids, ids[:5] + ids[-2:])
    assert report["eligible_method_seed_slots"] == 216
    assert report["verified_terminal_seed_slots"] == 42
    assert report["pending_condition_ids"] == ids[5:-2]
    assert report["requires_preserving_coordinator"] is True
    assert report["complete"] is False
    assert reconcile_schedule(ids, ids)["complete"] is True
    assert reconcile_schedule(ids, ids[:-1], blocked_ids=[ids[-1]])["complete"] is False
    with pytest.raises(SchedulingStop):
        reconcile_schedule(ids, [*ids, "replacement"])
    with pytest.raises(SchedulingStop):
        reconcile_schedule(ids, ids, blocked_ids=[ids[0]])


def test_only_two_extra_workers_after_old156_and_recovery_ready_and_finally_releases_slots():
    for old_done, recovery in [(False, True), (True, False), ("true", True)]:
        with pytest.raises(SchedulingStop):
            SyntheticWorkerBudget(old156_verified_terminal=old_done, recovery_coordinator_ready=recovery)
    budget = SyntheticWorkerBudget(old156_verified_terminal=True, recovery_coordinator_ready=True)
    with pytest.raises(RuntimeError):
        with budget.slot(), budget.slot():
            with pytest.raises(SchedulingStop):
                with budget.slot():
                    pass
            raise RuntimeError("synthetic worker interruption")
    assert budget.active == 0


def test_frozen_modern_engine_stops_fresh_driver_at_reserved_empty_manifest(tmp_path, monkeypatch):
    prepared, output = tmp_path / "prepared", tmp_path / "output"
    prepared.mkdir(); output.mkdir()
    (prepared / "preparation.json").write_text(json.dumps({"receptors": []}))
    (output / "runs.json").write_text("[]\n")
    config, executable = tmp_path / "config.json", tmp_path / "vina"
    config.write_text(json.dumps({"seeds": [104729, 130363, 155921]}))
    executable.write_text("synthetic executable identity\n")
    monkeypatch.setattr(engine.subprocess, "check_output", lambda *args, **kwargs: "AutoDock Vina v1.2.7")
    def forbidden(*args, **kwargs):
        raise AssertionError("A fresh legacy call must stop before a search")
    monkeypatch.setattr(engine.subprocess, "run", forbidden)
    with pytest.raises(FileExistsError, match="Run manifest exists"):
        engine.dock(prepared, output, executable, config)
    assert (output / "runs.json").read_text() == "[]\n"


def test_frozen_engine_can_fill_owned_empty_manifest_with_exactly_three_mocked_seed_jobs(tmp_path, monkeypatch):
    prepared, output = tmp_path / "prepared", tmp_path / "output"
    prepared.mkdir(); output.mkdir()
    metadata = {"box_center": [0, 0, 0], "box_size": [24, 24, 24], "ligand_file": "ligand.pdbqt",
                "receptors": [{"id": "r", "state": "holo", "assembly_condition": "full",
                               "pdb_file": "r.pdb", "pdbqt_file": "r.pdbqt"}]}
    (prepared / "preparation.json").write_text(json.dumps(metadata))
    (prepared / "reference.json").write_text(json.dumps({"names": ["C1"]}))
    for name in ["ligand.pdbqt", "r.pdb", "r.pdbqt"]:
        (prepared / name).write_text("synthetic input identity\n")
    (output / "runs.json").write_text("[]\n")
    config, executable = tmp_path / "config.json", tmp_path / "vina"
    seeds = [104729, 130363, 155921]
    config.write_text(json.dumps(dict(seeds=seeds, cpu=2, exhaustiveness=16, num_modes=20,
                                     min_rmsd_A=.5, energy_range_kcal_mol=10)))
    executable.write_text("synthetic executable identity\n")
    monkeypatch.setattr(engine.platform, "platform", lambda: "synthetic-platform")
    monkeypatch.setattr(engine.subprocess, "check_output", lambda *args, **kwargs: "AutoDock Vina v1.2.7")
    launched = []
    def fake_search(command, **kwargs):
        launched.append(int(command[command.index("--seed") + 1]))
        pose = Path(command[command.index("--out") + 1])
        atom = f"ATOM  {1:5d} {'C1':>4s} LIG A   1    {1.:8.3f}{2.:8.3f}{3.:8.3f}  1.00  0.00     0.000 {'C':>2s}"
        pose.write_text("MODEL 1\nREMARK VINA RESULT: -1 0 0\n" + atom + "\nENDMDL\n")
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(engine.subprocess, "run", fake_search)
    rows = engine.dock(prepared, output, executable, config, resume=True)
    assert launched == seeds
    assert [row["seed"] for row in rows] == seeds
    assert all(row["status"] == "ok" for row in rows)
