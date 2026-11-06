"""Execution-shard merging must preserve the full frozen scientific denominator."""

import json
from pathlib import Path
import shutil

import pytest

from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.shards import merge_shards


SEEDS = [104729, 130363, 155921]


def setup_shards(tmp_path):
    prepared, output = tmp_path / "prepared", tmp_path / "merged"
    prepared.mkdir()
    for name in ("ligand.pdbqt", "holo.pdb", "holo.pdbqt"):
        (prepared / name).write_text(name + " frozen bytes\n")
    write_json(prepared / "reference.json", {"names": ["C1"], "coordinates": [[0, 0, 0]]})
    write_json(prepared / "preparation.json", {"case_id": "control", "family_id": "control",
        "box_center": [0, 0, 0], "box_size": [20, 20, 20],
        "ligand_file": "ligand.pdbqt", "receptors": [{"id": "holo", "pdb_file": "holo.pdb", "pdbqt_file": "holo.pdbqt"}]})
    config = tmp_path / "config.json"
    write_json(config, {"seeds": SEEDS, "cpu": 2, "exhaustiveness": 16,
                        "num_modes": 20, "min_rmsd_A": 0.5, "energy_range_kcal_mol": 10})
    fingerprint = {
        "config_sha256": sha256(config), "preparation_sha256": sha256(prepared / "preparation.json"),
        "reference_sha256": sha256(prepared / "reference.json"), "ligand_sha256": sha256(prepared / "ligand.pdbqt"),
        "receptor_sha256": sha256(prepared / "holo.pdbqt"), "receptor_pdb_sha256": sha256(prepared / "holo.pdb"),
        "full_assembly_sha256": None, "executable_sha256": "a" * 64,
    }
    dirs = []
    for seed in SEEDS:
        directory = tmp_path / f"seed-{seed}"
        directory.mkdir()
        shutil.copyfile(config, directory / "configuration.json")
        pose_name, log_name = f"holo_{seed}.pdbqt", f"holo_{seed}.log"
        atom = f"ATOM  {1:5d} {'C1':>4s} GAL A{1:4d}    {0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00    +0.000 C"
        (directory / pose_name).write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\n" + atom + "\nENDMDL\n")
        (directory / log_name).write_text(f"real output seed {seed}\n")
        write_json(directory / "runs.json", [{"case_id": "control", "family_id": "control", "receptor": "holo", "seed": seed,
            "status": "ok", "fingerprint": fingerprint, "pose_file": pose_name, "log_file": log_name,
            "pose_sha256": sha256(directory / pose_name), "returned_poses": 1,
            "command": ["/original/vina", "--out", str(directory / pose_name), "--seed", str(seed),
                        "--receptor", str(prepared / "holo.pdbqt"), "--ligand", str(prepared / "ligand.pdbqt"),
                        "--scoring", "vina", "--cpu", "2", "--exhaustiveness", "16", "--num_modes", "20",
                        "--min_rmsd", "0.5", "--energy_range", "10",
                        "--center_x", "0", "--center_y", "0", "--center_z", "0",
                        "--size_x", "20", "--size_y", "20", "--size_z", "20"],
            "wall_seconds": 2.5, "cpu_seconds": 4.2}])
        dirs.append(directory)
    return prepared, output, config, dirs


def mutate(directory, callback):
    path = directory / "runs.json"
    rows = json.loads(path.read_text())
    callback(rows)
    write_json(path, rows)


def test_full_merge_preserves_commands_resources_config_and_source_manifests(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    source = json.loads((dirs[0] / "runs.json").read_text())[0]
    result = merge_shards(prepared, output, config, dirs)
    assert [row["seed"] for row in result] == SEEDS
    assert result[0]["command"] == source["command"]
    assert result[0]["wall_seconds"] == 2.5 and result[0]["cpu_seconds"] == 4.2
    assert sha256(output / "configuration.json") == sha256(config)
    for row in result:
        p = row["shard_provenance"]
        assert sha256(output / p["archived_run_manifest"]) == p["source_run_manifest_sha256"]
        assert sha256(output / row["pose_file"]) == row["pose_sha256"]
        assert sha256(output / row["log_file"]) == row["log_sha256"]


def test_partial_subset_cannot_be_published_as_full_denominator(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    with pytest.raises(ValueError, match="Incomplete full"):
        merge_shards(prepared, output, config, dirs[:2])
    assert not output.exists()


def test_duplicate_terminal_seed_rejected_even_when_same_bytes(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    duplicate = tmp_path / "duplicate"
    shutil.copytree(dirs[0], duplicate)
    def relocate_invocation(rows):
        command = rows[0]["command"]
        command[command.index("--out") + 1] = str(duplicate / rows[0]["pose_file"])
    mutate(duplicate, relocate_invocation)
    with pytest.raises(ValueError, match="Duplicate terminal"):
        merge_shards(prepared, output, config, dirs + [duplicate])


def test_changed_config_or_prepared_inputs_fail_before_output_mutation(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    (dirs[0] / "configuration.json").write_text('{"seeds":[104729]}\n')
    with pytest.raises(ValueError, match="configuration differs"):
        merge_shards(prepared, output, config, dirs)
    shutil.copyfile(config, dirs[0] / "configuration.json")
    (prepared / "holo.pdb").write_text("changed scoring coordinates")
    with pytest.raises(ValueError, match="prepared inputs"):
        merge_shards(prepared, output, config, dirs)
    assert not output.exists()


def test_pose_corruption_and_missing_log_fail(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    pose = dirs[0] / f"holo_{SEEDS[0]}.pdbqt"
    pose.write_text(pose.read_text() + "changed")
    with pytest.raises(ValueError, match="pose artifact checksum"):
        merge_shards(prepared, output, config, dirs)
    mutate(dirs[0], lambda rows: rows[0].update(pose_sha256=sha256(pose)))
    (dirs[0] / f"holo_{SEEDS[0]}.log").unlink()
    with pytest.raises(ValueError, match="log is missing"):
        merge_shards(prepared, output, config, dirs)


def test_incompatible_executable_and_live_shards_fail(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    mutate(dirs[1], lambda rows: rows[0]["fingerprint"].update(executable_sha256="b" * 64))
    with pytest.raises(ValueError, match="different executable"):
        merge_shards(prepared, output, config, dirs)
    mutate(dirs[1], lambda rows: rows[0]["fingerprint"].update(executable_sha256="a" * 64))
    mutate(dirs[1], lambda rows: rows[0].update(status="running"))
    with pytest.raises(ValueError, match="nonterminal work"):
        merge_shards(prepared, output, config, dirs)


def test_baseline_completed_seed_and_interrupted_attempt_are_preserved(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    shutil.copytree(dirs[0], output)
    original = json.loads((output / "runs.json").read_text())
    original[0]["command"][original[0]["command"].index("--out") + 1] = str(output / original[0]["pose_file"])
    interrupted = json.loads((dirs[1] / "runs.json").read_text())[0]
    interrupted.update(status="running")
    interrupted["command"][interrupted["command"].index("--out") + 1] = str(output / interrupted["pose_file"])
    for field in ("wall_seconds", "cpu_seconds", "returned_poses", "pose_sha256"):
        interrupted.pop(field)
    original.append(interrupted)
    write_json(output / "runs.json", original)
    (output / interrupted["log_file"]).write_text("partial cancelled attempt\n")
    result = merge_shards(prepared, output, config, dirs[1:], baseline_runs=original)
    attempt = result[1]["previous_attempts"][0]
    assert attempt["status"] == "cancelled" and attempt["status_before_merge"] == "running"
    assert "wall_seconds" not in attempt and "cpu_seconds" not in attempt
    assert set(attempt["missing_resource_measurements"]) >= {"wall_seconds", "cpu_seconds"}
    assert (output / attempt["preserved_artifacts"]["log_file"]).read_text() == "partial cancelled attempt\n"
    assert attempt["command"] == interrupted["command"]
    assert result[1]["command"] != interrupted["command"]


def test_existing_different_raw_artifact_is_saved_before_copy(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    output.mkdir()
    name = f"holo_{SEEDS[0]}.pdbqt"
    (output / name).write_text("prior incomplete output")
    result = merge_shards(prepared, output, config, dirs)
    previous = result[0]["preserved_previous_artifacts"]["pose_file"]
    assert (output / previous).read_text() == "prior incomplete output"
    assert sha256(output / name) == result[0]["pose_sha256"]


def test_atom_loss_is_rejected_even_if_pose_checksum_updated(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    pose = dirs[0] / f"holo_{SEEDS[0]}.pdbqt"
    pose.write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\nENDMDL\n")
    mutate(dirs[0], lambda rows: rows[0].update(pose_sha256=sha256(pose)))
    with pytest.raises(ValueError, match="Empty model"):
        merge_shards(prepared, output, config, dirs)


def test_claimed_config_hash_does_not_hide_changed_actual_search_settings(tmp_path):
    prepared, output, config, dirs = setup_shards(tmp_path)
    def change(rows):
        command = rows[0]["command"]
        command[command.index("--exhaustiveness") + 1] = "4"
    mutate(dirs[0], change)
    with pytest.raises(ValueError, match="changed frozen setting"):
        merge_shards(prepared, output, config, dirs)
