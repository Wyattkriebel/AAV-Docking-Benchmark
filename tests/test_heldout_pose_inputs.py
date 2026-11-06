import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from aav_benchmark.heldout_pose_inputs import (
    CASE_IDS, independent_box_size, validate_candidate_reuse, select_execution_candidates, export_independent_ligand,
)
from aav_benchmark.provenance import sha256
from scripts.verify_serialization_archives import read_archive

ROOT = Path(__file__).resolve().parents[1]


def selection():
    return {"case_id": "2RDK", "method": "fpocket", "configuration_sha256": "config", "prediction_input_sha256": "input",
            "native_coordinate_access": False, "predictor_status": "completed",
            "candidates": [{"rank": 1, "original_rank": 2, "center_A": [1., 2., 3.]},
                           {"rank": 2, "original_rank": 7, "center_A": [4., 5., 6.]}]}


def test_box_radius_is_translation_invariant_and_ceil_rule():
    side, radius = independent_box_size([[-5.1, 0, 0], [5.1, 0, 0]])
    assert side == 27 and radius == 5.1
    assert independent_box_size([[94.9, -4, 2], [105.1, -4, 2]])[0] == side
    assert independent_box_size([[0, 0, 0]]) == (24, 0.)
    assert independent_box_size([[-13, 0, 0], [13, 0, 0]]) == (42, 13.)


@pytest.mark.parametrize("xyz", [[], [[float("nan"), 0, 0]], [[0, 0]], [[-13.01, 0, 0], [13.01, 0, 0]]])
def test_invalid_or_oversize_independent_box_is_rejected(xyz):
    with pytest.raises(ValueError): independent_box_size(xyz)


def test_candidate_reuse_preserves_original_centers_objects_and_ranks():
    source = selection(); old = copy.deepcopy(source)
    result = validate_candidate_reuse(source, "2RDK", "fpocket", "config", "input")
    assert result is source["candidates"] and source == old
    source["predictor_status"] = "failed"
    with pytest.raises(ValueError, match="Failed predictor"):
        validate_candidate_reuse(source, "2RDK", "fpocket", "config", "input")
    source["candidates"] = []
    assert validate_candidate_reuse(source, "2RDK", "fpocket", "config", "input") == []


def test_top_one_execution_keeps_original_full_selection_unchanged():
    source = selection(); old = copy.deepcopy(source)
    assert select_execution_candidates(source["candidates"], 1) == [source["candidates"][0]]
    assert source == old
    assert select_execution_candidates([], 1) == []
    for budget in (None, True, 0, 2, 5):
        with pytest.raises(ValueError, match="top-one"):
            select_execution_candidates(source["candidates"], budget)


@pytest.mark.parametrize("field,value", [("case_id", "3ZWE"), ("method", "p2rank"), ("configuration_sha256", "other"),
                                        ("prediction_input_sha256", "other"), ("native_coordinate_access", True)])
def test_candidate_identity_cannot_change(field, value):
    source = selection(); source[field] = value
    with pytest.raises(ValueError): validate_candidate_reuse(source, "2RDK", "fpocket", "config", "input")


@pytest.mark.parametrize("changes", [({"rank": 3}), ({"original_rank": 1}), ({"center_A": [math.inf, 0, 0]})])
def test_candidate_class_order_original_order_and_finite_centers(changes):
    source = selection(); source["candidates"][1].update(changes)
    with pytest.raises(ValueError): validate_candidate_reuse(source, "2RDK", "fpocket", "config", "input")


@pytest.mark.parametrize("case_id", CASE_IDS)
def test_actual_v3_export_is_same_full_chemical_parent_and_literal_sdf(case_id, tmp_path):
    # The actual fixed six inputs are mandatory public archive members.
    parent = ROOT / "results/structural_test_v2"
    manifest = json.loads((parent / "source_archive_manifest.json").read_bytes())
    v2 = read_archive((parent / "source_preparation.tar.gz").read_bytes(), manifest, expected_count=len(manifest["files"]))
    parent = ROOT / "results/structural_pose_v3"
    manifest = json.loads((parent / "strict_sources_manifest.json").read_bytes())
    v3 = read_archive((parent / manifest["archive"]).read_bytes(), manifest, expected_count=len(manifest["files"]))
    source, strict = tmp_path / "initializer", tmp_path / "original_map"
    source.mkdir(); strict.mkdir()
    for name, archived_name in (("ligand.sdf", "ligand.sdf"), ("report.json", "initializer_report.json")):
        (source / name).write_bytes(v3[f"{case_id}/{archived_name}"])
    for name in ("ligand_atom_map.json", "chemistry.json"):
        (strict / name).write_bytes(v2[f"{case_id}/{name}"])
    before = sha256(source / "ligand.sdf")
    chemistry = export_independent_ligand(source / "ligand.sdf", source / "report.json", strict / "ligand_atom_map.json",
                                        strict / "chemistry.json", tmp_path / "export")
    assert sha256(tmp_path / "export/ligand.sdf") == before == sha256(source / "ligand.sdf")
    assert (tmp_path / "export/ligand.pdbqt").read_bytes() == v3[f"{case_id}/ligand.pdbqt"]
    assert chemistry["box_size"] == json.loads(v3[f"{case_id}/chemistry.json"])["box_size"]
    assert chemistry["native_coordinates_used_for_embedding"] is False
    assert chemistry["full_parent_canonical_key_and_name_bijection_preserved"] is True
    assert chemistry["no_new_embedding_or_minimization_during_export"] is True
    assert chemistry["serialization"]["rounded_CCD_ideal_validation"]["valid"] is True
    assert chemistry["serialization"]["max_heavy_serialization_displacement_A"] <= .0009
    assert all(row["status"] == "consistent" for row in chemistry["signed_chirality_tests"])
