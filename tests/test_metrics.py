"""Analytic positive/negative controls; these are not docking results."""

import unittest

import numpy as np

from aav_benchmark.metrics import (
    clash_metrics,
    contact_metrics,
    contact_pairs,
    fit_receptor,
    ligand_rmsd,
    mapped_coordinates,
    site_metrics,
    summarize_poses,
)


class PoseMetricControls(unittest.TestCase):
    def setUp(self):
        self.receptor = np.array([[0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=float)
        self.ligand = np.array([[0.2, 0.4, 2.5], [0.8, 0.4, 3.7], [0.2, 1.8, 2.5]], dtype=float)

    def test_native_is_zero_and_contacts_are_complete(self):
        self.assertEqual(ligand_rmsd(self.ligand, self.ligand), 0.0)
        native = contact_pairs(self.ligand, self.receptor)
        result = contact_metrics(self.ligand, self.receptor, native, receptor_residues=["A1", "A2", "A3", "A4"])
        self.assertEqual(result["contact_recall"], 1.0)
        self.assertEqual(result["contact_precision"], 1.0)
        self.assertEqual(result["interface_residue_recall"], 1.0)

    def test_translation_is_not_fitted_away(self):
        translated = self.ligand + [3, 4, 0]
        self.assertAlmostEqual(ligand_rmsd(translated, self.ligand), 5.0)
        self.assertAlmostEqual(ligand_rmsd(translated, self.ligand, receptor_source=self.receptor, receptor_target=self.receptor), 5.0)

    def test_common_receptor_transform_preserves_score(self):
        rotation = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=float)
        shift = np.array([20, -7, 8])
        pose = self.ligand + [0.6, 0.8, 0]
        transformed_receptor = self.receptor @ rotation + shift
        transformed_pose = pose @ rotation + shift
        rmsd = ligand_rmsd(transformed_pose, self.ligand, receptor_source=transformed_receptor, receptor_target=self.receptor)
        self.assertAlmostEqual(rmsd, 1.0)
        r, t = fit_receptor(transformed_receptor, self.receptor)
        np.testing.assert_allclose(transformed_receptor @ r + t, self.receptor, atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(r), 1.0)

    def test_atom_loss_and_duplicates_fail_visibly(self):
        for names, coords in [(["C1", "O1"], self.ligand[:2]), (["C1", "C1", "O2"], self.ligand), (["C1", "O1", "EXTRA"], self.ligand)]:
            with self.assertRaises(ValueError):
                mapped_coordinates(names, coords, ["C1", "O1", "O2"])
        with self.assertRaises(ValueError):
            ligand_rmsd(self.ligand[:2], self.ligand)

    def test_atom_order_can_be_mapped_exactly(self):
        result = mapped_coordinates(["O2", "C1", "O1"], self.ligand[[2, 0, 1]], ["C1", "O1", "O2"])
        np.testing.assert_array_equal(result, self.ligand)

    def test_only_explicit_symmetry_is_used(self):
        translation = [10, 0, 0]
        shifted = self.ligand + translation
        self.assertAlmostEqual(ligand_rmsd(shifted, self.ligand), 10)
        self.assertAlmostEqual(ligand_rmsd(shifted, self.ligand, symmetry_transforms=[(np.eye(3), [-10, 0, 0])]), 0)
        with self.assertRaises(ValueError):
            ligand_rmsd(shifted, self.ligand, symmetry_transforms=[(np.diag([-1, 1, 1]), [0, 0, 0])])

    def test_permutation_requires_complete_bijection(self):
        swapped = self.ligand[[1, 0, 2]]
        self.assertGreater(ligand_rmsd(swapped, self.ligand), 0)
        self.assertAlmostEqual(ligand_rmsd(swapped, self.ligand, permutations=[[1, 0, 2]]), 0)
        for bad in [[0, 0, 2], [0, 1], [0.0, 1.0, 2.0]]:
            with self.assertRaises(ValueError):
                ligand_rmsd(swapped, self.ligand, permutations=[bad])

    def test_alignment_cannot_use_a_degenerate_receptor(self):
        line = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
        with self.assertRaises(ValueError):
            fit_receptor(line, line)
        with self.assertRaises(ValueError):
            ligand_rmsd(self.ligand, self.ligand, receptor_source=self.receptor)

    def test_nonfinite_coordinates_fail(self):
        bad = self.ligand.copy()
        bad[0, 0] = float("nan")
        with self.assertRaises(ValueError):
            ligand_rmsd(bad, self.ligand)

    def test_reflection_is_not_an_alignment(self):
        reflected = self.receptor * [-1, 1, 1]
        rotation, translation = fit_receptor(reflected, self.receptor)
        self.assertAlmostEqual(np.linalg.det(rotation), 1)
        self.assertGreater(np.linalg.norm(reflected @ rotation + translation - self.receptor), 0.1)


class ContactAndRankingControls(unittest.TestCase):
    def test_contact_cutoff_is_strict_and_precision_penalizes_extra_contacts(self):
        native = {(0, 0)}
        ligand = [[0, 0, 0]]
        receptor = [[1, 0, 0], [3, 0, 0], [4, 0, 0]]
        self.assertEqual(contact_pairs(ligand, receptor), {(0, 0), (0, 1)})
        result = contact_metrics(ligand, receptor, native, receptor_residues=["a", "b", "c"])
        self.assertEqual(result["contact_recall"], 1)
        self.assertEqual(result["contact_precision"], 0.5)
        self.assertEqual(result["interface_residue_jaccard"], 0.5)
        missing = contact_metrics([[100, 0, 0]], receptor, native)
        self.assertEqual(missing["contact_recall"], 0)
        self.assertEqual(missing["contact_precision"], 0)

    def test_empty_native_contacts_are_undefined_recall(self):
        result = contact_metrics([[0, 0, 0]], [[0, 0, 0]], [])
        self.assertIsNone(result["contact_recall"])
        with self.assertRaises(ValueError):
            contact_metrics([[0, 0, 0]], [[0, 0, 0]], [(0, 1)])

    def test_apo_missing_atoms_stay_in_holo_recall_denominator(self):
        result = contact_metrics(
            [[0, 0, 0]], [[1, 0, 0]], {(0, 0), (0, 1)},
            receptor_native_indices=[1], native_receptor_count=2,
            receptor_residues=["missing_holo_residue", "present_apo_residue"],
        )
        self.assertEqual(result["native_contact_count"], 2)
        self.assertEqual(result["contact_recall"], 0.5)
        self.assertEqual(result["contact_precision"], 1)
        self.assertEqual(result["interface_residue_recall"], 0.5)
        with self.assertRaises(ValueError):
            contact_metrics([[0, 0, 0]], [[1, 0, 0]], {(0, 0)}, native_receptor_count=2)

    def test_clashes_use_explicit_radius_overlap(self):
        result = clash_metrics([[0, 0, 0]], [[2, 0, 0], [4, 0, 0]], [1.7], [1.7, 1.7])
        self.assertEqual(result["clash_pair_count"], 1)
        self.assertEqual(result["clashing_ligand_atom_count"], 1)
        self.assertAlmostEqual(result["max_overlap_A"], 1.4)
        with self.assertRaises(ValueError):
            clash_metrics([[0, 0, 0]], [[2, 0, 0]], [], [1.7])

    def test_explicit_score_rank_not_input_order_controls_top_five(self):
        rows = [{"rank": rank, "rmsd_A": rmsd} for rank, rmsd in [(6, 0.1), (2, 4), (1, 5), (3, 3), (4, 2), (5, 2.1)]]
        result = summarize_poses(rows)
        self.assertEqual(result["top1_rmsd_A"], 5)
        self.assertEqual(result["top5_min_rmsd_A"], 2)
        self.assertEqual(result["best_returned_rank"], 6)
        self.assertFalse(result["top1_success_2A"])
        self.assertTrue(result["top5_success_2A"])
        self.assertTrue(result["retained_success_1A"])
        self.assertFalse(result["top5_success_1A"])
        self.assertFalse(result["top1_success_given_retained_1A"])

    def test_no_returned_pose_remains_a_failure(self):
        result = summarize_poses([])
        self.assertEqual(result["returned_poses"], 0)
        self.assertFalse(result["top1_success_2A"])
        self.assertFalse(result["retained_success_2A"])
        self.assertIsNone(result["top1_success_given_retained_2A"])
        self.assertIsNone(result["best_returned_rmsd_A"])

    def test_invalid_ranks_and_scores_fail(self):
        for rows in [[{"rank": 2, "rmsd_A": 1}], [{"rank": 1, "rmsd_A": 1}, {"rank": 1, "rmsd_A": 2}], [{"rank": 1, "rmsd_A": float("nan")}]]:
            with self.assertRaises(ValueError):
                summarize_poses(rows)

    def test_site_equivalents_do_not_consume_multiple_top_five_slots(self):
        centers = [[10, 0, 0]] * 5 + [[1, 0, 0]]
        no_merge = site_metrics(centers, [[0, 0, 0]])
        merged = site_metrics(centers, [[0, 0, 0]], equivalent_classes=["a"] * 5 + ["b"])
        self.assertFalse(no_merge["top5_success_4A"])
        self.assertTrue(merged["top5_success_4A"])
        self.assertEqual(merged["distinct_returned_sites"], 2)

    def test_equivalent_site_is_evaluated_against_transformed_native(self):
        result = site_metrics(
            [[100, 0, 0], [0, 0, 0]], [[0, 0, 0]], equivalent_classes=["a", "a"],
            symmetry_transforms=[(np.eye(3), [100, 0, 0])],
        )
        self.assertEqual(result["distinct_returned_sites"], 1)
        self.assertTrue(result["top1_success_4A"])

    def test_malformed_empty_coordinates_fail(self):
        for empty in [[[]], [[], []], np.empty((0, 2))]:
            with self.assertRaises(ValueError):
                site_metrics(empty, [[0, 0, 0]])
        self.assertFalse(site_metrics([], [[0, 0, 0]])["top1_success_4A"])


if __name__ == "__main__":
    unittest.main()
