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
