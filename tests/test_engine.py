from pathlib import Path

import numpy as np
import pytest

from aav_benchmark.engine import read_poses


def atom(name, xyz=(1, 2, 3), atom_type="C", serial=1):
    return f"ATOM  {serial:5d} {name:>4s} LIG A   1    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00     0.000 {atom_type:>2s}"


def test_pdbqt_retains_score_rank_and_exact_atom_map(tmp_path):
    text = "MODEL 1\nREMARK VINA RESULT: -4.2 0 0\n" + atom("O1", (3, 2, 1), "OA") + "\n" + atom("C1") + "\n" + atom("H1", atom_type="HD") + "\nENDMDL\n"
    path = tmp_path / "pose.pdbqt"
    path.write_text(text)
    result = read_poses(path, ["C1", "O1"])
    assert result[0]["rank"] == 1
    assert result[0]["vina_score"] == -4.2
    np.testing.assert_array_equal(result[0]["coordinates"], [[1, 2, 3], [3, 2, 1]])


@pytest.mark.parametrize("text", [
    "", "MODEL 1\nREMARK VINA RESULT: -1 0 0\n" + atom("C1"),
    "MODEL 1\n" + atom("C1") + "\nENDMDL\n",
    "MODEL 1\nREMARK VINA RESULT: nan 0 0\n" + atom("C1") + "\nENDMDL\n",
    "MODEL 2\nREMARK VINA RESULT: -1 0 0\n" + atom("C1") + "\nENDMDL\n",
    "MODEL 1\nREMARK VINA RESULT: -1 0 0\n" + atom("C1") + "\n" + atom("C1") + "\nENDMDL\n",
])
def test_corrupt_or_incomplete_pose_is_not_a_success(tmp_path, text):
    path = tmp_path / "pose.pdbqt"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_poses(path, ["C1"])


@pytest.mark.parametrize("execution_label,expected", [(None, "true-site-local"),
                                                       ("true-site-local", "true-site-local"),
                                                       ("predicted-site-local", "predicted-site-local")])
def test_execution_information_is_explicit_and_legacy_labels_stay_unchanged(tmp_path, monkeypatch, execution_label, expected):
    import json
    from types import SimpleNamespace
    import aav_benchmark.engine as engine
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    prep = {"information_regime": "legacy descriptive text; not a machine execution label",
            "box_center": [0, 0, 0], "box_size": [24, 24, 24], "ligand_file": "ligand.pdbqt",
            "receptors": [{"id": "r", "state": "holo", "assembly_condition": "full",
                           "pdb_file": "r.pdb", "pdbqt_file": "r.pdbqt"}]}
    if execution_label is not None:
        prep["execution_information_regime"] = execution_label
    (prepared / "preparation.json").write_text(json.dumps(prep))
    (prepared / "reference.json").write_text(json.dumps({"names": ["C1"]}))
    for name in ("ligand.pdbqt", "r.pdb", "r.pdbqt"):
        (prepared / name).write_text("input fixture\n")
    executable = tmp_path / "vina"
    executable.write_text("test executable identity\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"seeds": [1], "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
                                  "min_rmsd_A": .5, "energy_range_kcal_mol": 10}))
    monkeypatch.setattr(engine.platform, "platform", lambda: "test-platform")
    monkeypatch.setattr(engine.subprocess, "check_output", lambda *a, **k: "AutoDock Vina v1.2.7")
    def run(command, **kwargs):
        output = Path(command[command.index("--out") + 1])
        output.write_text("MODEL 1\nREMARK VINA RESULT: -1 0 0\n" + atom("C1") + "\nENDMDL\n")
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(engine.subprocess, "run", run)
    rows = engine.dock(prepared, tmp_path / "results", executable, config)
    assert rows[0]["status"] == "ok"
    assert rows[0]["information_regime"] == expected
