"""Frozen-seed recovery controls; no docking executable is ever invoked."""
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.shards import merge_shards


SPEC = importlib.util.spec_from_file_location("blind_pose_recovery", Path(__file__).parents[1] / "scripts/resume_blind_pose.py")
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)
SEEDS = [104729, 130363, 155921]


def setup(tmp_path):
    prepared = tmp_path / "data/processed/blind_pose/C0/fpocket/box01"
    output = tmp_path / "results/blind_pose/C0/fpocket/box01"
    prepared.mkdir(parents=True)
    output.mkdir(parents=True)
    executable = tmp_path / "tools/vina-1.2.7"
    executable.parent.mkdir()
    executable.write_text("fingerprinted only; not an executable")
    config = tmp_path / "benchmarks/blind_pose_experiment.json"
    config.parent.mkdir()
    write_json(config, {"seeds": SEEDS, "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
                        "min_rmsd_A": 0.5, "energy_range_kcal_mol": 10})
    prep = {"complex_id": "C0", "case_id": "C0", "family_id": "family", "condition_id": "C0", "method": "fpocket",
            "execution_information_regime": "predicted-site-local", "ligand_file": "ligand.pdbqt",
            "box_center": [0, 0, 0], "box_size": [24, 24, 24],
            "receptors": [{"id": "neighborhood", "pdbqt_file": "neighborhood.pdbqt", "pdb_file": "neighborhood.pdb"}]}
    write_json(prepared / "preparation.json", prep)
    write_json(prepared / "reference.json", {"names": ["C1"], "pose_reference_eligible": False})
    for name in ("ligand.pdbqt", "neighborhood.pdbqt", "neighborhood.pdb"):
        (prepared / name).write_text("frozen " + name)
    return prepared, output, config


def record(root, prepared, directory, config, seed, status="ok"):
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config, directory / "configuration.json")
    pose, log = f"neighborhood_{seed}.pdbqt", f"neighborhood_{seed}.log"
    row = {"case_id": "C0", "family_id": "family", "receptor": "neighborhood", "seed": seed,
           "information_regime": "predicted-site-local", "vina_version": "AutoDock Vina v1.2.7", "status": status,
           "pose_file": pose, "log_file": log, "fingerprint": {
               "config_sha256": sha256(config), "preparation_sha256": sha256(prepared / "preparation.json"),
               "reference_sha256": sha256(prepared / "reference.json"), "ligand_sha256": sha256(prepared / "ligand.pdbqt"),
               "receptor_sha256": sha256(prepared / "neighborhood.pdbqt"),
               "receptor_pdb_sha256": sha256(prepared / "neighborhood.pdb"), "full_assembly_sha256": None,
               "executable_sha256": sha256(root / "tools/vina-1.2.7")}}
    row["command"] = [str(root / "tools/vina-1.2.7"), "--receptor", str(prepared / "neighborhood.pdbqt"),
        "--ligand", str(prepared / "ligand.pdbqt"), "--out", str(directory / pose), "--scoring", "vina",
        "--cpu", "2", "--seed", str(seed), "--exhaustiveness", "16", "--num_modes", "20", "--min_rmsd", "0.5",
        "--energy_range", "10", "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "24", "--size_y", "24", "--size_z", "24"]
    (directory / log).write_text(" 1 -5.000 0.000 0.000\n" if status == "ok" else "partial interrupted evidence\n")
    if status in {"ok", "failed"}:
        row.update(wall_seconds=3.5, cpu_seconds=5.8, child_peak_rss_bytes=1024)
    if status == "ok":
        atom = f"ATOM  {1:5d} {'C1':>4s} GAL A{1:4d}    {0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00    +0.000 C"
        (directory / pose).write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\n" + atom + "\nENDMDL\n")
        row.update(exit_code=0, returned_poses=1, pose_sha256=sha256(directory / pose))
    return row


def test_original_and_shard_records_merge_without_changing_completed_science(tmp_path):
    prepared, output, config = setup(tmp_path)
    completed = record(tmp_path, prepared, output, config, SEEDS[0])
    interrupted = record(tmp_path, prepared, output, config, SEEDS[1], "running")
    write_json(output / "runs.json", [completed, interrupted])
    baseline = json.loads((output / "runs.json").read_text())
    dirs = []
    for seed in SEEDS[1:]:
        directory = tmp_path / "results/blind_pose/_recovery/attempts/attempt" / f"seed-{seed}" / "results/blind_pose/C0/fpocket/box01"
        write_json(directory / "runs.json", [record(tmp_path, prepared, directory, config, seed)])
        recovery.validate_rows(tmp_path, prepared, directory, config, canonical=False)
        dirs.append(directory)
    rows = merge_shards(prepared, output, config, dirs, baseline_runs=baseline)
    recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)
    for key, value in completed.items():
        assert rows[0][key] == value
    old = rows[1]["previous_attempts"][0]
    assert old["status_before_merge"] == "running"
    assert "wall_seconds" not in old and "cpu_seconds" not in old
    assert (output / old["preserved_artifacts"]["log_file"]).read_text() == "partial interrupted evidence\n"
    assert rows[1]["command"] != interrupted["command"]
    # A relocated command needs immutable source evidence, including measured costs.
    rows[1]["cpu_seconds"] = 0
    write_json(output / "runs.json", rows)
    with pytest.raises(ValueError, match="Archived shard"):
        recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)


def test_pending_attempt_fingerprints_and_commands_are_checked(tmp_path):
    prepared, output, config = setup(tmp_path)
    row = record(tmp_path, prepared, output, config, SEEDS[0], "running")
    row["fingerprint"]["ligand_sha256"] = "0" * 64
    write_json(output / "runs.json", [row])
    with pytest.raises(ValueError, match="inputs"):
        recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)
    row = record(tmp_path, prepared, output, config, SEEDS[0], "running")
    row["command"][row["command"].index("--cpu") + 1] = "4"
    write_json(output / "runs.json", [row])
    with pytest.raises(ValueError, match="frozen option"):
        recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)


def test_completed_pose_atom_loss_and_log_score_disagreement_fail(tmp_path):
    prepared, output, config = setup(tmp_path)
    row = record(tmp_path, prepared, output, config, SEEDS[0])
    write_json(output / "runs.json", [row])
    (output / row["log_file"]).write_text("1 -6.000 0.000 0.000\n")
    with pytest.raises(ValueError, match="score table"):
        recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)
    (output / row["log_file"]).write_text("1 -5.000 0.000 0.000\n")
    (output / row["pose_file"]).write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\nENDMDL\n")
    row["pose_sha256"] = sha256(output / row["pose_file"])
    write_json(output / "runs.json", [row])
    with pytest.raises(ValueError, match="Empty model"):
        recovery.validate_rows(tmp_path, prepared, output, config, canonical=True)
