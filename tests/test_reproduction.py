import json

import pytest

from aav_benchmark.reproduction import compare_tables


def test_explicit_cross_platform_roundoff_leaves_default_byte_exact(tmp_path):
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    a.write_text("case_id,rank,rmsd_A\nAAV9,1,3.333033140649321\n")
    b.write_text("case_id,rank,rmsd_A\nAAV9,1,3.333033140649322\n")
    with pytest.raises(ValueError, match="bytes"):
        compare_tables(a, b)
    result = compare_tables(a, b, 1e-10)
    assert result["exact_recomputation"] is False
    assert 0 < result["max_numeric_difference"] < 1e-10


@pytest.mark.parametrize("changed", [{"success": False, "count": 3, "rmsd": 1.0},
                                     {"success": True, "count": 4, "rmsd": 1.0},
                                     {"success": True, "count": 3, "rmsd": 1.001}])
def test_numeric_tolerance_cannot_change_endpoints_or_counts(tmp_path, changed):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps({"success": True, "count": 3, "rmsd": 1.0}))
    b.write_text(json.dumps(changed))
    with pytest.raises(ValueError):
        compare_tables(a, b, 1e-10)


def test_numeric_looking_identity_is_never_rounded(tmp_path):
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    a.write_text("case_id,rmsd_A\n1.23,2.0\n")
    b.write_text("case_id,rmsd_A\n1.23000000001,2.0\n")
    with pytest.raises(ValueError, match="identity"):
        compare_tables(a, b, 1e-10)


@pytest.mark.parametrize("field", ["kd_uM", "wall_seconds", "score"])
def test_copied_float_measurements_and_resources_must_match_exactly(tmp_path, field):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps({field: 1.0, "rmsd_A": 2.0}))
    b.write_text(json.dumps({field: 1.0 + 1e-11, "rmsd_A": 2.0}))
    with pytest.raises(ValueError, match="Copied"):
        compare_tables(a, b, 1e-10)
