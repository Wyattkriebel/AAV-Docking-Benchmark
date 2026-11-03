"""Synthetic old-schema comparison controls; no evaluator or actual data calls."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from aav_benchmark import development_report_replay as replay
from aav_benchmark import report_replay as core

ROOT = Path(__file__).resolve().parents[1]


def synthetic_report():
    pose = {"rank": 1, "pooled_rank": 1, "vina_score": -6.25,
        "candidate_class_rank": 1, "candidate_original_rank": 1, "box_id": "box01",
        "global_RMSD_A": 1.25, "local_interface_RMSD_A": 1.125,
        "global_chemical_atom_map_index": 0, "local_interface_chemical_atom_map_index": 0,
        "full_parent_min_protein_distance_A": 2.75,
        "global_contacts": {"native_contact_count": 6, "contact_recall": 0.5, "contact_precision": 0.75,
                            "interface_residue_recall": 0.5, "interface_residue_precision": 0.75,
                            "interface_residue_jaccard": 0.25},
        "geometry_diagnostics": {"fixed_parent_graph_and_heavy_map_preserved": True,
            "stereochemistry_from_pose_3D_matches_input": True,
            "max_heavy_bond_length_change_from_input_A": 0.002,
            "max_ring_pair_distance_change_from_input_A": 0.003,
            "glycosidic_links": [{"heavy_indices": [2, 9], "distance_A": 1.4, "input_distance_A": 1.401}]}}
    second = copy.deepcopy(pose); second.update(rank=2, pooled_rank=2, vina_score=-6.0)
    summary = {"returned_pose_count": 2, "top1_RMSD_A": 1.25,
        "top5_min_RMSD_A": 1.125, "best_returned_RMSD_A": 1.125, "best_returned_pooled_rank": 2,
        **{field: True for field in replay._SUCCESS_FIELDS}}
    averaged = {field: 0.5 for field in replay._SUCCESS_FIELDS}
    seed = {"seed": 104729, "poses": [pose, second], "global": copy.deepcopy(summary), "local_interface": copy.deepcopy(summary)}
    budget = {"global": {**averaged, "seeds_with_returned_pose": 3},
              "local_interface": {**averaged, "seeds_with_returned_pose": 3}, "seed_evaluations": [seed]}
    align = {"matched_backbone_atom_count": 9, "backbone_RMSD_A": 0.002,
             "rotation": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "translation": [0.0, 0.0, 0.0]}
    aggregate = {"eligible_complexes": 2, "complex_ids": ["synthetic_a", "synthetic_b"], **averaged}
    return {"schema_version": 1, "status": "complete_development_evaluation", "metrics_eligible": True,
        "denominator": {"scheduled_seed_runs": 156, "terminal_seed_runs": 156},
        "evaluator_sha256": replay.EVALUATOR_SHA256,
        "source_sha256": {"source.bin": {"sha256": "a" * 64, "bytes": 19}},
        "resource_use": {"cpu_seconds": 2.25, "wall_seconds": 1.5},
        "conditions": [{"condition_id": "synthetic", "complex_id": "synthetic_a", "method": "fpocket",
            "atom_mapping": {"indices": [2, 4], "coordinates": [[1.25, 2.5, 3.75], [2.0, 3.0, 4.0]]},
            "native_controls": {"native_RMSD_A": 0.0, "native_contact_recall": 1.0,
                "known_translation_A": [3, 4, 0], "known_translation_RMSD_A": 5.0},
            "receptor_alignments": {alignment: copy.deepcopy(align) for alignment in replay._ALIGNMENTS},
            "site_budgets": {budget_id: copy.deepcopy(budget) for budget_id in ("1", "5")}}],
        "aggregates": {scope: {method: {budget_id: {alignment: copy.deepcopy(aggregate) for alignment in replay._ALIGNMENTS}
                           for budget_id in ("1", "5")} for method in ("fpocket", "p2rank", "oracle")}
                       for scope in ("full_ligand", "modeled_core")},
        "paired_contrasts": [{"atomic_scope": "full_ligand", "method": "fpocket", "predicted_site_budget": 5,
            "oracle_site_budget": 1, "alignment": "global", "pose_rank_budget": "top1",
            "complexes": [{"complex_id": "synthetic_a", "oracle_minus_predicted": 0.25}],
            "mean_oracle_minus_predicted": 0.25}]}


def get(value, path):
    for key in path:
        value = value[key]
    return value


def changed_at(value, path, replacement):
    result = copy.deepcopy(value); current = result
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = replacement
    return result


POSE = ("conditions", 0, "site_budgets", "1", "seed_evaluations", 0, "poses", 0)


@pytest.mark.parametrize("path", [
    POSE + ("global_RMSD_A",), POSE + ("local_interface_RMSD_A",),
    POSE + ("full_parent_min_protein_distance_A",),
    POSE + ("global_contacts", "contact_recall"),
    POSE + ("geometry_diagnostics", "max_ring_pair_distance_change_from_input_A"),
    POSE + ("geometry_diagnostics", "glycosidic_links", 0, "input_distance_A"),
    ("conditions", 0, "site_budgets", "5", "seed_evaluations", 0, "global", "top1_RMSD_A"),
    ("conditions", 0, "site_budgets", "1", "global", "top1_success_2A"),
    ("conditions", 0, "site_budgets", "5", "local_interface", "best_returned_success_3A"),
    ("conditions", 0, "receptor_alignments", "global", "rotation", 0, 0),
    ("conditions", 0, "receptor_alignments", "local_interface", "translation", 0),
    ("conditions", 0, "receptor_alignments", "global", "backbone_RMSD_A"),
    ("conditions", 0, "native_controls", "known_translation_RMSD_A"),
    ("aggregates", "full_ligand", "oracle", "1", "global", "top1_success_3A"),
    ("paired_contrasts", 0, "complexes", 0, "oracle_minus_predicted"),
    ("paired_contrasts", 0, "mean_oracle_minus_predicted"),
])
def test_old_schema_derived_values_use_separate_positioned_policy(path):
    expected = synthetic_report(); actual = changed_at(expected, path, get(expected, path) + 1e-12)
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values(expected, actual)
    result = replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-10)
    assert result["policy_id"] == replay.POLICY_ID
    assert result["tolerated_difference_count"] == 1
    assert result["max_difference_path"] == result["displayed_differences"][0]["path"] == core.json_pointer(path)
    assert result["displayed_differences"][0]["approved_path_rule"].startswith("/" + path[0])
    assert result["metrics_recomputed"] is result["global_comparison_policy_modified"] is False
    with pytest.raises(core.ReportReplayError, match="exceeds") as error:
        replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-13)
    assert error.value.path == core.json_pointer(path)


@pytest.mark.parametrize("path", [
    POSE + ("vina_score",),
    ("conditions", 0, "atom_mapping", "coordinates", 0, 0),
    ("conditions", 0, "native_controls", "native_RMSD_A"),
    ("resource_use", "wall_seconds"),
])
def test_native_values_raw_coordinates_and_analytic_zero_remain_exact(path):
    expected = synthetic_report(); actual = changed_at(expected, path, get(expected, path) + 1e-12)
    with pytest.raises(core.ReportReplayError) as error:
        replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-8)
    assert error.value.path == core.json_pointer(path)


@pytest.mark.parametrize("path,new", [
    (POSE + ("rank",), 2), (POSE + ("global_chemical_atom_map_index",), 1),
    (POSE + ("geometry_diagnostics", "stereochemistry_from_pose_3D_matches_input"), False),
    (("denominator", "terminal_seed_runs"), 155),
    (("source_sha256", "source.bin", "sha256"), "b" * 64),
    (("conditions", 0, "native_controls", "known_translation_A", 0), 4),
    (("conditions", 0, "site_budgets", "1", "seed_evaluations", 0, "global", "top1_success_2A"), False),
])
def test_counts_hashes_maps_and_seed_decisions_remain_exact(path, new):
    expected = synthetic_report()
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values(expected, changed_at(expected, path, new), numeric_tolerance=1e-8)


@pytest.mark.parametrize("replacement", [None, True, 1, "1.25"])
def test_mask_projection_cannot_hide_a_derived_float_type_change(replacement):
    expected = synthetic_report(); actual = changed_at(expected, POSE + ("global_RMSD_A",), replacement)
    for left, right in ((expected, actual), (actual, expected)):
        with pytest.raises(core.ReportReplayError, match="path or type"):
            replay.compare_development_report_values(left, right, numeric_tolerance=1e-8)


def test_exact_ratio_limit_and_diagnostic_truncation_are_reused():
    expected = synthetic_report(); path = ("paired_contrasts", 0, "mean_oracle_minus_predicted")
    expected = changed_at(expected, path, 1e-10); actual = changed_at(expected, path, -1e-27)
    with pytest.raises(core.ReportReplayError, match="exceeds") as error:
        replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-10)
    assert error.value.path == core.json_pointer(path) and error.value.exact_difference is not None
    expected = synthetic_report(); actual = changed_at(expected, path, get(expected, path) + 1e-12)
    result = replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-10, differences_limit=0)
    assert result["tolerated_difference_count"] == 1 and result["displayed_differences"] == []
    assert result["differences_display_truncated"] is True


def test_inputs_and_frozen_global_policy_are_never_modified():
    expected = synthetic_report(); original = copy.deepcopy(expected)
    rules_before = {key: tuple(value) for key, value in core._RULES_BY_FINAL.items()}
    core_sha = hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest()
    replay.compare_development_report_values(expected, copy.deepcopy(expected), numeric_tolerance=1e-10)
    assert expected == original
    assert {key: tuple(value) for key, value in core._RULES_BY_FINAL.items()} == rules_before
    assert core.POLICY_ID == "global_report_derived_paths_v1"
    assert hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest() == core_sha == replay.CORE_SHA256
    # The new-report shape is still outside the old development policy.
    left = {"per_source_pose_metrics": [{"global_RMSD_A": 1.0}]}
    right = {"per_source_pose_metrics": [{"global_RMSD_A": 1.0 + 1e-12}]}
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values(left, right, numeric_tolerance=1e-10)
    assert core.compare_report_values(left, right, numeric_tolerance=1e-10)["tolerated_difference_count"] == 1


def test_frozen_core_binding_fails_closed_without_changing_any_file(monkeypatch):
    monkeypatch.setattr(replay, "CORE_SHA256", "0" * 64)
    with pytest.raises(core.ReportReplayError, match="core changed"):
        replay.compare_development_report_values({}, {})


def test_shape_and_model_order_are_exact_and_nonfinite_input_is_rejected():
    expected = synthetic_report(); actual = copy.deepcopy(expected)
    actual["conditions"][0]["site_budgets"]["1"]["seed_evaluations"][0]["poses"].reverse()
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values(expected, actual, numeric_tolerance=1e-8)
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values(expected, changed_at(expected, POSE + ("global_RMSD_A",), float("nan")))
    with pytest.raises(core.ReportReplayError):
        replay.compare_development_report_values({"value": 1}, {"value": 1.0})


def test_projection_depth_limit_is_a_controlled_input_rejection():
    nested = 0
    for _ in range(600):
        nested = {"nested": nested}
    with pytest.raises(core.ReportReplayError, match="excessively nested"):
        replay.compare_development_report_values(nested, nested)


def test_existing_file_cli_is_byte_exact_by_default_and_stdlib_only(tmp_path):
    expected = synthetic_report(); actual = changed_at(expected, POSE + ("global_RMSD_A",), 1.25 + 1e-12)
    left, right = tmp_path / "expected.json", tmp_path / "actual.json"
    left.write_text(json.dumps(expected)); right.write_text(json.dumps(actual))
    command = [sys.executable, "-S", str(ROOT / "scripts/compare_development_report_replay.py"), str(left), str(right)]
    strict = subprocess.run(command, capture_output=True, text=True)
    assert strict.returncode == 1 and strict.stderr == ""
    tolerant = subprocess.run(command + ["--numeric-tolerance", "1e-10"], capture_output=True, text=True)
    assert tolerant.returncode == 0 and tolerant.stderr == ""
    assert json.loads(tolerant.stdout)["tolerated_difference_count"] == 1
    right.write_bytes(left.read_bytes())
    assert replay.compare_development_report_files(left, right)["status"] == "exact_byte_match"
    right.write_text(json.dumps(expected, indent=2))
    with pytest.raises(core.ReportReplayError, match="bytes differ"):
        replay.compare_development_report_files(left, right)
    assert replay.compare_development_report_files(left, right, structural=True)["status"] == "exact_structural_match"
    left.write_text('{"value":NaN}'); right.write_bytes(left.read_bytes())
    with pytest.raises(core.ReportReplayError, match="Nonfinite"):
        replay.compare_development_report_files(left, right)
    imported = subprocess.run([sys.executable, "-S", "-c", "import sys;sys.path.insert(0,sys.argv[1]);import aav_benchmark.development_report_replay;print(sorted(n for n in sys.modules if n.startswith('aav_benchmark.')))", str(ROOT)], check=True, capture_output=True, text=True)
    assert imported.stdout.strip() == "['aav_benchmark.development_report_replay', 'aav_benchmark.report_replay']"
