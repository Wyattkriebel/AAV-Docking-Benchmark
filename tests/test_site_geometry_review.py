"""Independent controls for symmetry, receptor frames and footprint interpretation."""
import gzip
from pathlib import Path

import gemmi
import numpy as np
import pytest

from aav_benchmark.aav import PROTEIN_RESIDUES
from aav_benchmark.site_geometry import capsid_symmetry, distinct_candidates, evaluate_candidates


def test_spatial_merge_preserves_highest_rank_center_and_rotates_about_receptor_center():
    # Symmetry axis is displaced from the coordinate origin: x' = 20-x.
    operations = [(np.eye(3), np.zeros(3)),
                  (np.diag([-1., -1., 1.]), np.array([20., 0., 0.]))]
    rows = [{"rank": 1, "center_A": [12, 2, 0], "score": 9},
            {"rank": 2, "center_A": [8, -2, 0], "score": 8},
            {"rank": 3, "center_A": [12.2, 2, 0], "score": 7},
            {"rank": 4, "center_A": [10, 10, 0], "score": 6}]
    unique, mapping = distinct_candidates(rows, operations, merge_distance_A=.3)
    assert [r["original_rank"] for r in unique] == [1, 4]
    assert unique[0]["center_A"] == [12., 2., 0.]
    assert unique[0]["score"] == 9
    assert [r["candidate_class"] for r in mapping] == [1, 1, 1, 2]


def test_good_pairwise_chain_fits_do_not_certify_an_asymmetric_60_chain_arrangement():
    keys, xyz = [], []
    # Every chain has exactly the same fold, so all individual chain fits are
    # perfect. Their irregular arrangement is not a finite capsid symmetry.
    for c in range(60):
        for n, point in enumerate(((0., 0., 0.), (1., 0., 0.), (0., 1., 1.)), 1):
            keys.append((str(c), n, "", "ALA", "CA"))
            xyz.append(np.array(point) + [c * c * 3., c * 11., 0.])
    with pytest.raises(ValueError):
        capsid_symmetry(keys, np.asarray(xyz))


def test_actual_aav_receptor_operators_reproduce_all_native_gal_copies_posthoc():
    archive = Path(__file__).resolve().parents[1] / "data/archive/7MTZ.cif.gz"
    block = gemmi.cif.read_string(gzip.decompress(archive.read_bytes()).decode()).sole_block()
    model = gemmi.make_structure_from_block(block)[0]
    keys, xyz = [], []
    for chain in model:
        for residue in chain:
            if residue.het_flag != "A" or residue.name not in PROTEIN_RESIDUES:
                continue
            for atom in residue:
                if atom.name == "CA":
                    keys.append((chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, "CA"))
                    xyz.append(list(atom.pos))
    # Only protein C-alpha arrays are supplied to symmetry certification. Native
    # ligand extraction follows its completed return and cannot determine it.
    operations, audit = capsid_symmetry(keys, np.asarray(xyz), tolerance_A=.05)
    assert len(operations) == 60
    assert all(row["whole_capsid_atom_correspondence_bijective"] for row in audit["operators"])
    assert max(row["whole_capsid_max_CA_displacement_A"] for row in audit["operators"]) <= .05
    assert max(row["chain_proposal_whole_capsid_max_A"] for row in audit["operators"]) > .05
    native = []
    for chain in model:
        for residue in chain:
            if residue.name == "GAL":
                native.append(np.asarray([list(atom.pos) for atom in sorted(residue, key=lambda atom: atom.name)
                                          if not atom.element.is_hydrogen]))
    assert len(native) == 60
    assert all(site.shape == (12, 3) for site in native)
    native = np.asarray(native)
    images = set()
    for rotation, translation in operations:
        per_atom = np.linalg.norm(native - (native[0] @ rotation + translation)[None, :, :], axis=2)
        nearest_copy = int(np.argmin(np.mean(per_atom ** 2, axis=1)))
        assert per_atom[nearest_copy].max() <= .05
        images.add(nearest_copy)
    assert len(images) == 60


def test_apo_footprint_uses_input_coordinates_while_center_uses_holo_rotation_and_translation():
    keys = [("A", 1, "", "ALA", "CA"), ("A", 2, "B", "ASN", "CA"),
            ("B", 1, "", "ALA", "CA")]
    xyz = np.array([[0., 0., 0.], [2., 0., 0.], [20., 0., 0.]])
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    translation = np.array([101., 45., 3.])
    center = np.array([1., 0., 0.])
    truth = {"evaluation_only": True, "sites": [
        {"id": "holo", "coordinates": [(center @ rotation + translation).tolist()],
         "interface_residue_keys": [list(k[:4]) for k in keys[:2]]}]}
    result = evaluate_candidates(
        [{"rank": 1, "center_A": center.tolist(), "residue_keys": [list(keys[1][:4])]}],
        truth, keys, xyz, rotation=rotation, translation=translation, footprint_radius_A=8.)
    row = result["candidates"][0]
    assert row["center_to_native_atom_A"] == pytest.approx(0.)
    assert row["evaluated_center_A"] == pytest.approx(center @ rotation + translation)
    assert row["common_radius_footprint"]["recall"] == 1.
    assert row["method_footprint"]["recall"] == .5
    assert row["common_radius_residue_keys"] == [list(k[:4]) for k in keys[:2]]


def test_common_radius_and_reported_method_footprints_are_distinct_diagnostics():
    keys = [("A", 1, "", "ALA", "CA"), ("A", 2, "", "ASN", "CA"),
            ("A", 3, "", "GLU", "CA")]
    xyz = np.array([[0., 0., 0.], [7., 0., 0.], [7.5, 0., 0.]])
    truth = {"evaluation_only": True, "sites": [
        {"id": "native", "coordinates": [[0., 0., 0.]],
         "interface_residue_keys": [list(keys[0][:4]), list(keys[1][:4])]}]}
    result = evaluate_candidates(
        [{"rank": 1, "center_A": [0., 0., 0.], "footprint_residues": [list(keys[0][:4])]}],
        truth, keys, xyz, footprint_radius_A=8.)
    row = result["candidates"][0]
    assert row["common_radius_footprint"] == {"recall": 1., "precision": 2 / 3}
    assert row["method_footprint"] == {"recall": .5, "precision": 1.}


def test_closest_copy_distance_cannot_be_combined_with_other_copy_interface():
    keys = [("A", 1, "", "ALA", "CA"), ("B", 1, "", "ALA", "CA")]
    truth = {"evaluation_only": True, "sites": [
        {"id": "one", "coordinates": [[0., 0., 0.]], "interface_residue_keys": [list(keys[0][:4])]},
        {"id": "two", "coordinates": [[20., 0., 0.]], "interface_residue_keys": [list(keys[1][:4])]}]}
    result = evaluate_candidates(
        [{"rank": 1, "center_A": [19., 0., 0.], "residue_keys": [list(keys[0][:4])]}],
        truth, keys, np.array([[0., 0., 0.], [20., 0., 0.]]))
    row = result["candidates"][0]
    assert row["closest_native_site_id"] == "two"
    assert row["center_to_native_atom_A"] == 1.
    assert row["method_footprint"]["recall"] == 0.
    assert row["common_radius_footprint"]["recall"] == 1.


def test_no_candidates_retains_zero_recovery_for_every_prespecified_threshold():
    truth = {"evaluation_only": True, "sites": [
        {"id": "one", "coordinates": [[0., 0., 0.]], "interface_residue_keys": [["A", 1, "", "ALA"]]}]}
    result = evaluate_candidates([], truth, [("A", 1, "", "ALA", "CA")], [[0., 0., 0.]])
    for radius in (3, 4, 5):
        assert result[f"top1_success_{radius}A"] is False
        assert result[f"top5_success_{radius}A"] is False
    assert result["returned_distinct_candidates"] == 0
    assert result["top1_common_interface_recall"] == 0.
    assert result["top1_method_interface_recall"] == 0.
    assert result["top1_distance_A"] is None
