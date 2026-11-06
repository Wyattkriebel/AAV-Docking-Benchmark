"""Analytic JSON comparison controls only; no scientific report is generated."""
import copy
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

from aav_benchmark.report_replay import (
    ReportReplayError, compare_report_files, compare_report_values, json_pointer, load_report_json,
)

ROOT = Path(__file__).resolve().parents[1]
RANKER = "VinaCarb1.0_Vina1.1.2_core_shared_zero_base_coeff0"


def report():
    pose = {"rank": 1, "pooled_rank": 1, "model_id": "synthetic_model", "global_RMSD_A": 1.25,
            "local_interface_RMSD_A": 1.125, "full_parent_min_protein_distance_A": 2.75,
            "full_observed_assembly_min_protein_distance_A": 2.625,
            "global_contacts": {"native_contact_count": 8, "shared_contact_count": 4,
                                "contact_recall": 0.5, "interface_residue_jaccard": 0.25},
            "geometry_diagnostics": {"fixed_parent_graph_and_heavy_map_preserved": True,
                "max_heavy_bond_length_change_from_input_A": 0.002,
                "max_ring_pair_distance_change_from_input_A": 0.003,
                "glycosidic_links": [{"heavy_indices": [2, 9], "distance_A": 1.4, "input_distance_A": 1.401}]},
            "native_method_score": {"method": "synthetic", "value": -5.125}}
    summary = {"returned_pose_count": 1, "top1_RMSD_A": 1.25, "top5_min_RMSD_A": 1.25,
               "best_returned_RMSD_A": 1.25, "top1_success_2A": True}
    method = {"poses": [pose], "global": summary,
              "resource_use": {"cpu_seconds": 2.25, "wall_seconds": 1.5}}
    control = {"copy_id": "synthetic_copy", "alignments": {"global": {
        "matched_backbone_atoms": 9, "backbone_RMSD_A": 0.002,
        "rotation": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "translation": [0.0, 0.0, 0.0]}},
        "native_controls": {"native_RMSD_A": 0.0, "known_translation_RMSD_A": 5.0, "known_translation_A": [3, 4, 0]}}
    seed = {"seed": 104729, "search_methods": {"modern_vina_1_2_7": copy.deepcopy(method)},
            "ranking_methods": {RANKER: copy.deepcopy(method)}, "poses": [copy.deepcopy(pose)],
            "global": copy.deepcopy(summary), "elapsed_seconds": 4.25,
            "model_reference_controls": [{"model_id": "synthetic_model", "protein_fits_and_native_controls": [copy.deepcopy(control)]}]}
    return {"schema_version": 1, "counts": {"eligible": 2},
        "source_sha256_inventory": {"source.bin": {"sha256": "a" * 64, "bytes": 19}},
        "coordinates": [[1.25, 2.5, 3.75]], "global_RMSD_A": 1.25,
        "conditions": [{"condition_id": "synthetic", "native_controls": [control], "seeds": [seed],
            "site": {"candidates": [{"candidate_class_rank": 1, "center_to_native_atom_A": 2.5,
                                     "qualified_copy_distances_A": {"synthetic_copy": 2.5}}],
                     "1": {"minimum_DCA_A": 2.5, "success_4A": True, "qualified_copy_coverage_4A": 0.5,
                           "qualified_copy_max_box_atom_coverage": [0.5]}}}],
        "paired_family_statistics": [{"confidence": 0.95, "threshold_A": 2,
            "difference_a_minus_b": 0.25, "confidence_interval": [0.125, 0.5],
            "per_family": [{"family_id": "synthetic_family", "recovery_a": 0.5, "recovery_b": 0.25, "difference": 0.25}],
            "eligible_and_conditional_summary": {"overall": [{"method": "synthetic", "family_macro_recovery": 0.5}]}}],
        "per_source_pose_metrics": [copy.deepcopy(pose)]}


def changed_at(value, path, replacement):
    result = copy.deepcopy(value)
    current = result
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = replacement
    return result


def test_exact_structural_comparison_counts_values_and_keeps_unknowns_exact():
    value = report()
    result = compare_report_values(value, copy.deepcopy(value))
    assert result["status"] == "exact_structural_match"
    assert result["compared_float_values"] > result["approved_derived_float_values"] > 0
    assert result["tolerated_difference_count"] == result["max_numeric_difference"] == 0
    assert result["metrics_recomputed"] is result["source_validation_performed"] is False


@pytest.mark.parametrize("path", [
    ("conditions", 0, "seeds", 0, "search_methods", "modern_vina_1_2_7", "poses", 0, "global_RMSD_A"),
    ("conditions", 0, "seeds", 0, "poses", 0, "full_observed_assembly_min_protein_distance_A"),
    ("conditions", 0, "seeds", 0, "poses", 0, "global_contacts", "contact_recall"),
    ("conditions", 0, "seeds", 0, "poses", 0, "geometry_diagnostics", "glycosidic_links", 0, "distance_A"),
    ("conditions", 0, "seeds", 0, "ranking_methods", RANKER, "global", "top1_RMSD_A"),
    ("conditions", 0, "native_controls", 0, "alignments", "global", "rotation", 0, 0),
    ("conditions", 0, "seeds", 0, "model_reference_controls", 0, "protein_fits_and_native_controls", 0, "alignments", "global", "translation", 0),
    ("conditions", 0, "site", "candidates", 0, "qualified_copy_distances_A", "synthetic_copy"),
    ("conditions", 0, "site", "1", "qualified_copy_max_box_atom_coverage", 0),
    ("paired_family_statistics", 0, "confidence_interval", 0),
    ("paired_family_statistics", 0, "eligible_and_conditional_summary", "overall", 0, "family_macro_recovery"),
    ("per_source_pose_metrics", 0, "global_RMSD_A"),
])
def test_explicit_tolerance_accepts_only_positioned_derived_metrics(path):
    expected = report(); current = expected
    for key in path:
        current = current[key]
    actual = changed_at(expected, path, current + 1e-12)
    with pytest.raises(ReportReplayError):
        compare_report_values(expected, actual)
    result = compare_report_values(expected, actual, numeric_tolerance=1e-10)
    assert result["tolerated_difference_count"] == 1
    assert 0 < result["max_numeric_difference"] <= 1e-10
    assert result["max_difference_path"] == result["displayed_differences"][0]["path"] == json_pointer(path)
    assert result["relative_tolerance"] == 0
    with pytest.raises(ReportReplayError, match="exceeds"):
        compare_report_values(expected, actual, numeric_tolerance=1e-13)


@pytest.mark.parametrize("path", [
    ("global_RMSD_A",), ("coordinates", 0, 0),
    ("conditions", 0, "seeds", 0, "elapsed_seconds"),
    ("conditions", 0, "seeds", 0, "search_methods", "modern_vina_1_2_7", "resource_use", "cpu_seconds"),
    ("conditions", 0, "seeds", 0, "search_methods", "modern_vina_1_2_7", "poses", 0, "native_method_score", "value"),
    ("per_source_pose_metrics", 0, "native_method_score", "value"),
    ("paired_family_statistics", 0, "confidence"),
    ("conditions", 0, "native_controls", 0, "native_controls", "native_RMSD_A"),
])
def test_tolerance_never_relaxes_raw_copied_or_unknown_values(path):
    expected = report(); current = expected
    for key in path:
        current = current[key]
    actual = changed_at(expected, path, current + 1e-12)
    with pytest.raises(ReportReplayError, match="Exact floating") as error:
        compare_report_values(expected, actual, numeric_tolerance=1e-8)
    assert error.value.path == json_pointer(path)


@pytest.mark.parametrize("path,new", [
    (("counts", "eligible"), 3),
    (("source_sha256_inventory", "source.bin", "sha256"), "b" * 64),
    (("conditions", 0, "seeds", 0, "global", "top1_success_2A"), False),
    (("conditions", 0, "seeds", 0, "seed"), 104730),
    (("conditions", 0, "seeds", 0, "poses", 0, "rank"), 2),
    (("conditions", 0, "seeds", 0, "poses", 0, "geometry_diagnostics", "fixed_parent_graph_and_heavy_map_preserved"), False),
    (("conditions", 0, "native_controls", 0, "native_controls", "known_translation_A", 0), 4),
])
def test_identities_decisions_counts_and_source_hashes_stay_exact(path, new):
    expected = report()
    with pytest.raises(ReportReplayError):
        compare_report_values(expected, changed_at(expected, path, new), numeric_tolerance=1e-8)


@pytest.mark.parametrize("left,right", [(1, 1.0), (True, 1), (False, 0.0), (None, 0), ([1], (1,))])
def test_type_changes_are_never_numeric_roundoff(left, right):
    with pytest.raises(ReportReplayError):
        compare_report_values({"value": left}, {"value": right}, numeric_tolerance=1e-8)


def test_object_shape_list_lengths_and_order_are_preserved():
    for expected, actual in [({"x": 1}, {"x": 1, "y": 2}), ([1, 2], [1]),
                             ([{"id": "a"}, {"id": "b"}], [{"id": "b"}, {"id": "a"}])]:
        with pytest.raises(ReportReplayError):
            compare_report_values(expected, actual, numeric_tolerance=1e-8)
    assert compare_report_values({"a": 1, "b": 2}, {"b": 2, "a": 1})["status"] == "exact_structural_match"


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_equal_nonfinite_values_are_always_rejected(bad):
    with pytest.raises(ReportReplayError, match="Nonfinite"):
        compare_report_values({"nested": [bad]}, {"nested": [bad]})


@pytest.mark.parametrize("raw", [b'{"value":NaN}', b'{"value":Infinity}', b'{"value":-Infinity}',
                                b'{"value":1e999}', b'{"x":1,"x":1}', b'{"x":}',
                                '{"x":1}'.encode("utf-16"), b'\xff'])
def test_strict_json_rejects_nonfinite_duplicates_and_non_utf8(raw):
    with pytest.raises(ReportReplayError):
        load_report_json(raw)


@pytest.mark.parametrize("bad", [True, -1, float("nan"), float("inf"), 1e-7, 10 ** 1000, "1e-10"])
def test_numeric_tolerance_is_explicit_finite_bounded_and_not_boolean(bad):
    with pytest.raises(ReportReplayError):
        compare_report_values({}, {}, numeric_tolerance=bad)


def test_signed_zero_is_exact_for_raw_floats_and_explicit_for_derived_floats():
    with pytest.raises(ReportReplayError):
        compare_report_values({"raw": -0.0}, {"raw": 0.0}, numeric_tolerance=1e-8)
    expected = report(); path = ("per_source_pose_metrics", 0, "global_RMSD_A")
    expected = changed_at(expected, path, -0.0); actual = changed_at(expected, path, 0.0)
    with pytest.raises(ReportReplayError):
        compare_report_values(expected, actual)
    result = compare_report_values(expected, actual, numeric_tolerance=1e-10)
    assert result["tolerated_difference_count"] == 1 and result["max_numeric_difference"] == 0


def test_diagnostic_limit_does_not_stop_comparison_or_hide_later_failure():
    expected = {"per_source_pose_metrics": [{"global_RMSD_A": 1.0} for _ in range(5)]}
    actual = {"per_source_pose_metrics": [{"global_RMSD_A": 1.0 + 1e-12} for _ in range(5)]}
    result = compare_report_values(expected, actual, numeric_tolerance=1e-10, differences_limit=1)
    assert result["tolerated_difference_count"] == 5 and len(result["displayed_differences"]) == 1
    assert result["differences_display_truncated"] is True
    actual["per_source_pose_metrics"][-1]["global_RMSD_A"] = 1.1
    with pytest.raises(ReportReplayError, match="exceeds"):
        compare_report_values(expected, actual, numeric_tolerance=1e-10, differences_limit=0)
    for bad in (-1, True, 1.5):
        with pytest.raises(ReportReplayError):
            compare_report_values({}, {}, differences_limit=bad)


def test_numeric_difference_overflow_fails_closed():
    expected = {"per_source_pose_metrics": [{"global_RMSD_A": 1e308}]}
    actual = {"per_source_pose_metrics": [{"global_RMSD_A": -1e308}]}
    with pytest.raises(ReportReplayError, match="exceeds"):
        compare_report_values(expected, actual, numeric_tolerance=1e-8)


def test_exact_binary_bound_cannot_be_relaxed_by_rounded_subtraction():
    expected = {"paired_family_statistics": [{"difference_a_minus_b": 1e-10}]}
    actual = {"paired_family_statistics": [{"difference_a_minus_b": -1e-27}]}
    assert abs(1e-10 - (-1e-27)) == 1e-10  # Ordinary subtraction loses the excess.
    with pytest.raises(ReportReplayError, match="exceeds") as error:
        compare_report_values(expected, actual, numeric_tolerance=1e-10)
    assert error.value.exact_difference is not None
    actual["paired_family_statistics"][0]["difference_a_minus_b"] = 0.0
    result = compare_report_values(expected, actual, numeric_tolerance=1e-10)
    assert result["tolerated_difference_count"] == 1
    assert result["max_numeric_difference_exact_ratio"] == result["tolerance_exact_ratio"]


def test_cyclic_python_values_and_excessively_nested_json_fail_as_input_errors():
    cyclic = {}; cyclic["cycle"] = cyclic
    with pytest.raises(ReportReplayError, match="Cyclic or excessively"):
        compare_report_values(cyclic, cyclic)
    with pytest.raises(ReportReplayError, match="Invalid UTF-8 JSON"):
        load_report_json("[" * 2000 + "0" + "]" * 2000)


def test_default_file_mode_requires_bytes_and_optional_structural_mode_is_exact(tmp_path):
    expected, actual = tmp_path / "expected.json", tmp_path / "actual.json"
    expected.write_text('{"value":1,"other":2}\n'); actual.write_text('{"other": 2, "value": 1}\n')
    with pytest.raises(ReportReplayError, match="bytes differ"):
        compare_report_files(expected, actual)
    result = compare_report_files(expected, actual, structural=True)
    assert result["exact_bytes"] is False and result["comparison_mode"] == "exact_structure_and_types"
    assert result["expected_file"]["sha256"] != result["actual_file"]["sha256"]
    assert compare_report_files(expected, expected)["status"] == "exact_byte_match"
    expected.write_text('{"value":NaN}'); actual.write_bytes(expected.read_bytes())
    with pytest.raises(ReportReplayError, match="Nonfinite"):
        compare_report_files(expected, actual)


def test_json_pointer_escapes_unusual_dictionary_keys():
    with pytest.raises(ReportReplayError) as error:
        compare_report_values({"a/b": {"~": 1}}, {"a/b": {"~": 2}})
    assert error.value.path == "/a~1b/~0"


def test_stdlib_only_cli_compares_existing_synthetic_files_without_scientific_imports(tmp_path):
    expected, actual = tmp_path / "expected.json", tmp_path / "actual.json"
    expected.write_text(json.dumps({"per_source_pose_metrics": [{"global_RMSD_A": 1.0}]}))
    actual.write_text(json.dumps({"per_source_pose_metrics": [{"global_RMSD_A": 1.0 + 1e-12}]}))
    command = [sys.executable, "-S", str(ROOT / "scripts/compare_report_replay.py"), str(expected), str(actual)]
    rejected = subprocess.run(command, capture_output=True, text=True)
    assert rejected.returncode == 1 and json.loads(rejected.stdout)["status"] == "comparison_rejected"
    accepted = subprocess.run(command + ["--numeric-tolerance", "1e-10"], capture_output=True, text=True)
    assert accepted.returncode == 0 and accepted.stderr == ""
    result = json.loads(accepted.stdout)
    assert result["tolerated_difference_count"] == 1 and result["metrics_recomputed"] is False
    imported = subprocess.run([sys.executable, "-S", "-c", "import sys;sys.path.insert(0,sys.argv[1]);import aav_benchmark.report_replay;print(sorted(n for n in sys.modules if n.startswith('aav_benchmark.')))", str(ROOT)], capture_output=True, text=True, check=True)
    assert imported.stdout.strip() == "['aav_benchmark.report_replay']"
