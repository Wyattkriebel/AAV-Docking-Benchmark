"""Geometry and provenance controls for the receptor-only random baseline."""
import inspect
import json
import unittest

import gemmi
import numpy as np

from aav_benchmark.surface import (
    DEFAULT_POINTS_PER_ATOM, build_surface_geometry, element_radii_A,
    fibonacci_sphere, iter_solvent_accessible_points, random_surface_patches,
    sample_surface_patches,
)


def collect_surface(coordinates, elements, **kwargs):
    batches = list(iter_solvent_accessible_points(coordinates, elements, **kwargs))
    return tuple(np.concatenate([row[i] for row in batches]) for i in range(3))


class ProbeSurfaceGeometry(unittest.TestCase):
    def test_isolated_atom_has_analytic_expanded_sphere_area_and_radius(self):
        center = np.array([7.0, -5.0, 11.0])
        points, owners, areas = collect_surface([center], ["C"], points_per_atom=128)
        radius = float(gemmi.Element("C").vdw_r) + 1.4
        self.assertEqual(len(points), 128)
        np.testing.assert_allclose(np.linalg.norm(points - center, axis=1), radius, atol=1e-12)
        np.testing.assert_array_equal(owners, 0)
        self.assertAlmostEqual(float(areas.sum()), 4 * np.pi * radius ** 2, places=10)

    def test_nested_atom_points_are_buried_and_not_sampled(self):
        # The oxygen expanded sphere lies entirely inside the calcium sphere.
        xyz = np.array([[0., 0., 0.], [0.05, 0., 0.]])
        points, owners, areas = collect_surface(xyz, ["Ca", "O"], points_per_atom=128)
        self.assertEqual(len(points), 128)
        np.testing.assert_array_equal(owners, 0)
        radius = float(gemmi.Element("Ca").vdw_r) + 1.4
        self.assertAlmostEqual(float(areas.sum()), 4 * np.pi * radius ** 2, places=10)
        result = random_surface_patches(xyz, ["Ca", "O"], ["outer", "buried"],
                                        max_sites=1, patch_radius_A=4.)
        self.assertEqual(result["sites"][0]["source_atom_index"], 0)

    def test_overlapping_equal_spheres_match_analytic_union_surface(self):
        separation = 3.
        xyz = np.array([[0., 0., 0.], [separation, 0., 0.]])
        points, _, areas = collect_surface(xyz, ["C", "C"], points_per_atom=2048)
        radius = float(gemmi.Element("C").vdw_r) + 1.4
        expected_area = 4 * np.pi * radius ** 2 + 2 * np.pi * radius * separation
        self.assertAlmostEqual(float(areas.sum()) / expected_area, 1., delta=.005)
        distances = np.linalg.norm(points[:, None, :] - xyz[None, :, :], axis=2)
        self.assertTrue(np.all(distances >= radius - 1e-9))

    def test_coincident_identical_spheres_do_not_double_surface_area(self):
        points, owners, areas = collect_surface([[0., 0., 0.]] * 2, ["C"] * 2)
        radius = float(gemmi.Element("C").vdw_r) + 1.4
        self.assertEqual(len(points), DEFAULT_POINTS_PER_ATOM)
        np.testing.assert_array_equal(owners, 0)
        self.assertAlmostEqual(float(areas.sum()), 4 * np.pi * radius ** 2, places=10)

    def test_different_radius_occluder_is_checked_beyond_nearest_atom(self):
        xyz = np.array([[0., 0., 0.], [3.0, 0., 0.], [2.2, 2.3, 0.]])
        elements = ["O", "Ca", "C"]
        points, _, _ = collect_surface(xyz, elements, points_per_atom=256)
        radii = element_radii_A(elements) + 1.4
        distance = np.linalg.norm(points[:, None, :] - xyz[None, :, :], axis=2)
        self.assertTrue(np.all(distance >= radii[None, :] - 1e-9))

    def test_probe_exposure_includes_unconnected_cavity(self):
        # A closed shell of expanded atom spheres surrounds the origin. Its
        # cavity boundary is probe-exposed even though this method makes no
        # connectivity claim. Both radial sides must remain in the sample.
        xyz = fibonacci_sphere(160) * 12.
        points, _, _ = collect_surface(xyz, ["C"] * len(xyz), points_per_atom=128)
        radii_from_origin = np.linalg.norm(points, axis=1)
        self.assertLess(float(radii_from_origin.min()), 10.)
        self.assertGreater(float(radii_from_origin.max()), 14.)


class RandomSurfacePatches(unittest.TestCase):
    def setUp(self):
        self.xyz = np.column_stack((np.arange(8) * 20., np.zeros(8), np.zeros(8)))
        self.elements = ["C"] * len(self.xyz)
        self.residues = [("A", i + 1, "", "ALA") for i in range(len(self.xyz))]

    def test_frozen_seed_exact_budget_distinct_centers_and_radius_footprints(self):
        kwargs = dict(seed=42, patch_radius_A=4., max_sites=5)
        result = random_surface_patches(self.xyz, self.elements, self.residues, **kwargs)
        repeated = random_surface_patches(self.xyz, self.elements, self.residues, **kwargs)
        self.assertEqual(result, repeated)
        self.assertTrue(result["budget_complete"])
        self.assertEqual(result["returned_site_count"], 5)
        centers = np.array([site["center_A"] for site in result["sites"]])
        distances = np.linalg.norm(centers[:, None] - centers[None, :], axis=2)
        self.assertGreaterEqual(float(distances[np.triu_indices(5, 1)].min()), 8.)
        for site in result["sites"]:
            source = site["source_atom_index"]
            self.assertEqual(site["footprint_residues"], [self.residues[source]])
            self.assertEqual(site["footprint_atom_count"], 1)
            self.assertEqual(site["footprint_radius_A"], 4.)
        json.dumps(result, allow_nan=False)

    def test_batch_size_does_not_change_selection(self):
        small = random_surface_patches(self.xyz, self.elements, self.residues,
                                       patch_radius_A=4., batch_atoms=1)
        large = random_surface_patches(self.xyz, self.elements, self.residues,
                                       patch_radius_A=4., batch_atoms=512)
        self.assertEqual(small["sites"], large["sites"])
        self.assertEqual(small["sampling"]["surface_mask_bytes"], len(self.xyz) * DEFAULT_POINTS_PER_ATOM)

    def test_cached_geometry_is_reusable_frozen_and_equal_to_wrapper(self):
        geometry = build_surface_geometry(self.xyz, self.elements, self.residues)
        original_mask = geometry.exposed_mask.copy()
        cached = sample_surface_patches(geometry, seed=42, patch_radius_A=4.)
        direct = random_surface_patches(self.xyz, self.elements, self.residues, seed=42, patch_radius_A=4.)
        self.assertEqual(cached, direct)
        sample_surface_patches(geometry, seed=43, patch_radius_A=4.)
        np.testing.assert_array_equal(geometry.exposed_mask, original_mask)
        self.assertFalse(geometry.exposed_mask.flags.writeable)
        self.xyz[0, 0] = -100.
        self.assertEqual(geometry.coordinates_A[0, 0], 0.)
        self.assertTrue(self.xyz.flags.writeable)

    def test_seed_changes_random_ranking(self):
        first = random_surface_patches(self.xyz, self.elements, self.residues, seed=1, patch_radius_A=4.)
        second = random_surface_patches(self.xyz, self.elements, self.residues, seed=2, patch_radius_A=4.)
        self.assertNotEqual(first["sites"], second["sites"])
        self.assertEqual(first["parameters"]["seed"], 1)
        self.assertTrue(first["reference_independent"])

    def test_area_weighting_uses_expanded_sphere_area(self):
        # Reproduce the first RNG draw against unequal sphere areas. Uniform
        # atom sampling would choose oxygen for this seed; area weighting picks
        # the larger calcium sphere instead.
        xyz = [[0., 0., 0.], [30., 0., 0.]]
        seed = next(seed for seed in range(1000)
                    if .5 < np.random.default_rng(seed).random() <
                    (float(gemmi.Element("Ca").vdw_r) + 1.4) ** 2 /
                    sum((float(gemmi.Element(symbol).vdw_r) + 1.4) ** 2 for symbol in ("Ca", "O")))
        result = random_surface_patches(xyz, ["Ca", "O"], ["Ca", "O"], seed=seed,
                                        max_sites=1, patch_radius_A=4.)
        self.assertEqual(result["sites"][0]["source_atom_index"], 0)

    def test_source_only_symmetry_excludes_equivalent_centers(self):
        xyz = np.array([[-20., 0., 0.], [20., 0., 0.]])
        rotation = np.diag([-1., -1., 1.])
        result = random_surface_patches(xyz, ["C", "C"], ["A", "B"], max_sites=2,
                                        patch_radius_A=4., symmetry_transforms=[(rotation, np.zeros(3))])
        self.assertEqual(result["returned_site_count"], 1)
        self.assertFalse(result["budget_complete"])
        self.assertEqual(result["parameters"]["symmetry_transform_count"], 2)

    def test_budget_exhaustion_is_explicit(self):
        result = random_surface_patches([[0., 0., 0.]], ["C"], ["A"], max_sites=5, patch_radius_A=4.)
        self.assertFalse(result["budget_complete"])
        self.assertEqual(result["returned_site_count"], 1)
        self.assertEqual(result["requested_site_count"], 5)
        self.assertEqual(result["sampling"]["remaining_eligible_point_count"], 0)

    def test_no_reference_or_ligand_parameter(self):
        parameters = set(inspect.signature(random_surface_patches).parameters)
        self.assertFalse(any("native" in value or "reference" in value or "ligand" in value or "box" in value
                             for value in parameters))

    def test_invalid_receptor_and_configuration_fail(self):
        for elements in (["H"], ["D"], ["X"], ["Carbon"], [""]):
            with self.subTest(elements=elements), self.assertRaises(ValueError):
                random_surface_patches([[0., 0., 0.]], elements, ["A"])
        invalid = [{"seed": -1}, {"seed": True}, {"max_sites": 0}, {"max_sites": 1.5},
                   {"patch_radius_A": 0}, {"probe_radius_A": float("nan")},
                   {"points_per_atom": 4}, {"batch_atoms": 0}, {"min_center_separation_A": 0},
                   {"symmetry_transforms": [(np.diag([-1., 1., 1.]), np.zeros(3))]}]
        for kwargs in invalid:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                random_surface_patches(self.xyz, self.elements, self.residues, **kwargs)
        with self.assertRaises(ValueError):
            random_surface_patches([[np.nan, 0., 0.]], ["C"], ["A"])
        with self.assertRaises(ValueError):
            random_surface_patches([[0., 0., 0.]], ["C"], [])


if __name__ == "__main__":
    unittest.main()
