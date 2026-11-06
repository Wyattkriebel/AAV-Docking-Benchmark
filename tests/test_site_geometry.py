import numpy as np
import pytest

from aav_benchmark.site_geometry import distinct_candidates, evaluate_candidates, capsid_symmetry


def test_receptor_symmetry_collapses_ranked_copies_without_native_input():
    rows = [{"rank": 1, "center_A": [3, 0, 0]}, {"rank": 2, "center_A": [-3, 0, 0]},
            {"rank": 3, "center_A": [0, 9, 0]}]
    rotation = np.diag([-1., -1., 1.])
    unique, mapping = distinct_candidates(rows, [(np.eye(3), np.zeros(3)), (rotation, np.zeros(3))], .1)
    assert len(unique) == 2
    assert mapping[1]["candidate_class"] == 1
    assert unique[1]["original_rank"] == 3


def test_distance_and_coverage_use_the_same_native_copy_and_receptor_frame():
    truth = {"evaluation_only": True, "sites": [
        {"id": "one", "coordinates": [[0, 0, 0]], "interface_residue_keys": [["A", 1, "", "ALA"]]},
        {"id": "two", "coordinates": [[20, 0, 0]], "interface_residue_keys": [["B", 1, "", "ALA"]]}]}
    keys = [("A", 1, "", "ALA", "CA"), ("B", 1, "", "ALA", "CA")]
    result = evaluate_candidates([{"rank": 1, "center_A": [19, 0, 0], "residue_keys": [keys[1][:4]]}], truth, keys, np.array([[0, 0, 0], [20, 0, 0]]))
    assert result["top1_success_4A"]
    assert result["candidates"][0]["closest_native_site_id"] == "two"
    assert result["top1_common_interface_recall"] == 1


def test_missing_predictions_retain_zero_recovery():
    result = evaluate_candidates([], {"evaluation_only": True, "sites": [{"coordinates": [[0, 0, 0]]}]},
                                 [("A", 1, "", "ALA", "CA")], [[0, 0, 0]])
    assert result["top1_success_4A"] is False
    assert result["top5_success_4A"] is False


def test_homologous_chain_fit_alone_cannot_certify_capsid():
    with pytest.raises(ValueError, match="60"):
        capsid_symmetry([("A", 1, "", "ALA", "CA")], [[0, 0, 0]])
