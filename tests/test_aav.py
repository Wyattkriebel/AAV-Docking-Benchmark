"""Preparation controls: assembly boundaries, atom identity and protein topology."""
import tempfile
import unittest
from pathlib import Path

import gemmi
import numpy as np

from aav_benchmark.aav import (
    _pdb_atoms, _protein_topology, _restore_terminal_oxygen_names, box_distance, fit,
    materialize_assembly, protein_atoms, select_neighborhood, validate_native_stereochemistry,
)


def pdb_atom(serial, name, chain, residue, xyz, element="O", icode=""):
    return (f"ATOM  {serial:5d} {name:>4s} ALA {chain}{residue:4d}{icode:1s}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n")


class PreparationControls(unittest.TestCase):
    def test_neighborhood_covers_box_corners_plus_cutoff(self):
        atoms = {
            ("A", 1, "", "ALA", "CA"): np.array([15., 10., 10.]),
            ("A", 2, "", "ALA", "CA"): np.array([100., 0., 0.]),
            ("B", 1, "", "ALA", "CA"): np.array([17., 10., 10.]),
        }
        selected, omitted = select_neighborhood(atoms, [0, 0, 0], [20, 20, 20], 6)
        self.assertEqual(selected, ["A"])
        self.assertAlmostEqual(omitted, 7.)
        # A radial 20 A shell would miss this atom at the box corner.
        self.assertGreater(np.linalg.norm(atoms[("A", 1, "", "ALA", "CA")]), 20)
        self.assertEqual(box_distance([[10, 10, 10]], [0, 0, 0], [20, 20, 20])[0], 0)

    def test_invalid_box_fails(self):
        with self.assertRaises(ValueError):
            box_distance([[0, 0, 0]], [0, 0, 0], [20, -1, 20])
        with self.assertRaises(ValueError):
            select_neighborhood({("A", 1, "", "ALA", "CA"): [0, 0, 0]}, [0, 0, 0], [20, 20, 20], 0)

    def test_apo_fit_preserves_known_rigid_transform(self):
        source = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
        angle = .7
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        target = source @ rotation + [4, -2, 9]
        recovered_r, recovered_t = fit(source, target)
        np.testing.assert_allclose(source @ recovered_r + recovered_t, target, atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(recovered_r), 1.)

    def test_degenerate_receptor_fit_fails(self):
        with self.assertRaises(ValueError):
            fit([[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]])

    def test_insertion_codes_are_distinct_atom_identities(self):
        text = pdb_atom(1, "O", "A", 42, [0, 0, 0]) + pdb_atom(2, "O", "A", 42, [1, 0, 0], icode="A")
        atoms = _pdb_atoms(text)
        self.assertEqual(len(atoms), 2)
        self.assertIn(("A", 42, "A", "ALA", "O"), atoms)

    def test_duplicate_atom_identity_fails(self):
        text = pdb_atom(1, "O", "A", 42, [0, 0, 0]) + pdb_atom(2, "O", "A", 42, [1, 0, 0])
        with self.assertRaises(ValueError):
            _pdb_atoms(text)

    def test_only_exact_terminal_oxygen_permutation_is_restored(self):
        before_text = pdb_atom(1, "O", "A", 42, [0, 0, 0]) + pdb_atom(2, "OXT", "A", 42, [2, 0, 0])
        swapped = pdb_atom(1, "O", "A", 42, [2, 0, 0]) + pdb_atom(2, "OXT", "A", 42, [0, 0, 0])
        restored, permutations = _restore_terminal_oxygen_names(swapped, _pdb_atoms(before_text))
        self.assertEqual(len(permutations), 1)
        self.assertEqual(_pdb_atoms(restored), _pdb_atoms(before_text))
        moved = pdb_atom(1, "O", "A", 42, [2.1, 0, 0]) + pdb_atom(2, "OXT", "A", 42, [0, 0, 0])
        unchanged, permutations = _restore_terminal_oxygen_names(moved, _pdb_atoms(before_text))
        self.assertEqual(unchanged, moved)
        self.assertEqual(permutations, [])

    def test_close_nonbonded_chains_cannot_create_peptide_link(self):
        records = []
        serial = 1
        for chain, shift in [("A", 0.), ("B", 3.)]:
            for name, xyz, element in [("N", [0, 0, 0], "N"), ("CA", [1.45, 0, 0], "C"),
                                       ("C", [2., 1.3, 0], "C"), ("O", [1.5, 2.3, 0], "O"),
                                       ("CB", [2., -1., 1.], "C")]:
                records.append(pdb_atom(serial, name, chain, 1, np.asarray(xyz) + [shift, 0, 0], element))
                serial += 1
            records.append("TER\n")
        raw, bonds, log = _protein_topology("".join(records))
        self.assertEqual(set(raw), {"A:1", "B:1"})
        self.assertEqual(bonds, {})
        self.assertTrue(log["rejected_distance_inferred_nonbonded_edges"])

    def test_native_stereochemistry_rejects_mirror_image_and_accepts_reordering(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("F[C@](Cl)(Br)I")
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 20261024
        self.assertEqual(AllChem.EmbedMolecule(mol, parameters), 0)
        names = [f"A{i}" for i in range(mol.GetNumAtoms())]
        for atom, name in zip(mol.GetAtoms(), names):
            atom.SetProp("ccd_name", name)
        xyz = np.array(mol.GetConformer().GetPositions())
        report = validate_native_stereochemistry(mol, names[::-1], xyz[::-1])
        self.assertEqual(report["status"], "passed")
        mirrored = xyz.copy()
        mirrored[:, 0] *= -1
        with self.assertRaisesRegex(ValueError, "stereochemistry differs"):
            validate_native_stereochemistry(mol, names, mirrored)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_native_stereochemistry(mol, ["A0"] * len(names), xyz)

    def test_blank_source_cannot_be_treated_as_assembly(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "blank.cif"
            path.write_text("\r\n")
            with self.assertRaises((ValueError, RuntimeError)):
                materialize_assembly(path)

    def test_duplicate_protein_atom_not_overwritten(self):
        text = pdb_atom(1, "CA", "A", 42, [0, 0, 0], "C") + pdb_atom(2, "CA", "A", 42, [1, 0, 0], "C")
        model = gemmi.read_pdb_string(text)[0]
        with self.assertRaises(ValueError):
            protein_atoms(model)


if __name__ == "__main__":
    unittest.main()
