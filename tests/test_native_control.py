"""Native-score controls must preserve deposited coordinates and source inputs."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from aav_benchmark.native_control import _ligand_atoms, build_native_ligand, parse_score_only
from aav_benchmark.provenance import sha256


class NativeLigandControls(unittest.TestCase):
    def setUp(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.source = Path(self.temporary.name) / "source"
        self.output = Path(self.temporary.name) / "control"
        self.source.mkdir()
        molecule = Chem.AddHs(Chem.MolFromSmiles("CO"))
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 20261024
        self.assertEqual(AllChem.EmbedMolecule(molecule, parameters), 0)
        heavy_names, h_number = ["C1", "O1"], 0
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() > 1:
                name = heavy_names[atom.GetIdx()]
            else:
                h_number += 1
                name = f"H{h_number}"
            info = Chem.AtomPDBResidueInfo()
            info.SetName(name.rjust(4))
            info.SetResidueName("GAL")
            info.SetResidueNumber(1)
            info.SetChainId("L")
            atom.SetMonomerInfo(info)
        setup = MoleculePreparation().prepare(molecule)[0]
        text, ok, error = PDBQTWriterLegacy.write_string(setup)
        self.assertTrue(ok, error)
        (self.source / "galactose.pdbqt").write_text(text)
        writer = Chem.SDWriter(str(self.source / "galactose.sdf"))
        writer.write(molecule)
        writer.close()
        self.preparation = {"ligand": {"ccd_atom_names": heavy_names}, "ligand_file": "galactose.pdbqt"}
        (self.source / "preparation.json").write_text(json.dumps(self.preparation))
        # Native coordinates deliberately differ from independently generated input.
        self.reference = {"names": ["O1", "C1"], "coordinates": [[7.4, 7, 8], [6, 7, 8]]}
        (self.source / "reference.json").write_text(json.dumps(self.reference))

    def test_exact_native_coordinates_hydrogen_ledger_and_source_immutability(self):
        originals = {path: sha256(path) for path in self.source.iterdir()}
        result = build_native_ligand(self.source, self.output)
        atoms, torsdof = _ligand_atoms((self.output / "native_ligand.pdbqt").read_text())
        for name, xyz in zip(self.reference["names"], self.reference["coordinates"]):
            np.testing.assert_array_equal(atoms[name]["coordinates"], xyz)
        self.assertEqual(result["native_heavy_rmsd_A"], 0)
        self.assertFalse(result["count_as_docking_success"])
        self.assertFalse(result["native_heavy_atoms_optimized"])
        self.assertFalse(result["hydrogen_optimization"])
        self.assertEqual(result["modeled_hydrogen_count"], 4)
        self.assertEqual(len(result["hydrogen_placement"]), 4)
        self.assertEqual(CounterParents(result["hydrogen_placement"]), {"C1": 3, "O1": 1})
        self.assertEqual(sum(row["retained_in_pdbqt"] for row in result["hydrogen_placement"]), 1)
        self.assertEqual(result["native_torsdof"], result["source_torsdof"])
        self.assertEqual(torsdof, result["native_torsdof"])
        self.assertEqual({path: sha256(path) for path in originals}, originals)

    def test_missing_native_atom_fails(self):
        (self.source / "reference.json").write_text(json.dumps({"names": ["C1"], "coordinates": [[6, 7, 8]]}))
        with self.assertRaisesRegex(ValueError, "native atom map"):
            build_native_ligand(self.source, self.output)

    def test_source_graph_order_mismatch_fails(self):
        self.preparation["ligand"]["ccd_atom_names"] = ["O1", "C1"]
        (self.source / "preparation.json").write_text(json.dumps(self.preparation))
        with self.assertRaisesRegex(ValueError, "source SDF/PDBQT atom map"):
            build_native_ligand(self.source, self.output)

    def test_nonfinite_native_coordinate_fails(self):
        self.reference["coordinates"][0][0] = float("nan")
        (self.source / "reference.json").write_text(json.dumps(self.reference))
        with self.assertRaisesRegex(ValueError, "finite"):
            build_native_ligand(self.source, self.output)


class VinaScoreParsingControls(unittest.TestCase):
    def test_exact_score_only_result(self):
        text = "AutoDock Vina v1.2.7\nEstimated Free Energy of Binding   : -1.858 (kcal/mol) [=(1)+(2)+(3)-(4)]\n"
        self.assertEqual(parse_score_only(text), -1.858)

    def test_missing_duplicate_nonfinite_results_fail(self):
        line = "Estimated Free Energy of Binding : -1.858 (kcal/mol)\n"
        for text in ("", line + line, "Estimated Free Energy of Binding : nan (kcal/mol)"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_score_only(text)


def CounterParents(rows):
    from collections import Counter
    return dict(Counter(row["parent_ccd_atom"] for row in rows))


if __name__ == "__main__":
    unittest.main()
