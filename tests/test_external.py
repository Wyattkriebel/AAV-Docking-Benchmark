"""Offline chemistry and input-integrity controls, not performance evidence."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from aav_benchmark.external import (
    _ccd_ligand,
    build_ligand,
    case_definition,
    cohort,
    fetch_external,
    pdb_heavy_atoms,
    prepare,
    prepare_external_receptor,
    verify_external,
)


FIXTURES = Path(__file__).parent / "fixtures" / "external"


class ExternalChemistryControls(unittest.TestCase):
    def test_released_trisaccharide_has_complete_atom_map_and_independent_conformer(self):
        molecule, reference, atom_map, metadata = build_ligand(FIXTURES, "154L")
        self.assertEqual(len(reference), 43)
        self.assertEqual(len(atom_map), 43)
        self.assertEqual(Chem.RemoveHs(molecule).GetNumAtoms(), 43)
        self.assertEqual(len({atom["name"] for atom in atom_map}), 43)
        self.assertTrue(all(row["handedness_verified"] for row in metadata["stereocenter_checks"]))
        # Independently generated input cannot reproduce the bound conformation
        # merely by translation (all pair distances would then agree).
        generated = molecule.GetConformer().GetPositions()[:43]
        native = np.array([atom["coordinates"] for atom in reference])
        generated_distances = np.linalg.norm(generated[:, None] - generated[None, :], axis=2)
        native_distances = np.linalg.norm(native[:, None] - native[None, :], axis=2)
        self.assertGreater(float(np.max(np.abs(generated_distances-native_distances))), 0.1)

    def test_atom_loss_fails_before_embedding(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory)
            (raw / "NAG.cif").write_bytes((FIXTURES / "NAG.cif").read_bytes())
            lines = (FIXTURES / "154L_l_b.pdb").read_text().splitlines()
            # Drop an actual carbon atom, not a glycosidic leaving group.
            lines = [line for line in lines if not (line.startswith("ATOM")
                     and int(line[22:26]) == 1 and line[12:16].strip() == "C2")]
            (raw / "154L_l_b.pdb").write_text("\n".join(lines)+"\n")
            with self.assertRaisesRegex(ValueError, "Reference atom mismatch"):
                build_ligand(raw, "154L")

    def test_anomer_inconsistency_fails_visibly(self):
        # The released 1QFO GLC reducing end has the opposite handedness to
        # CCD GLC at C1. It must not silently pass using residue-name identity.
        with self.assertRaisesRegex(ValueError, "stereochemistry mismatch 1:C1"):
            build_ligand(FIXTURES, "1QFO")

    def test_ccd_link_replacement_preserves_lacnac_stereochemistry(self):
        links = [{"donor_residue": 2, "donor_atom": "C1", "donor_leaving_atom": "O1",
                  "acceptor_residue": 1, "acceptor_atom": "O4", "acceptor_leaving_atom": "HO4",
                  "bond_order": "SING"}]
        ideal, atom_map, metadata = _ccd_ligand(FIXTURES, ["NAG", "GAL"], links)
        verified, _, _, verified_metadata = build_ligand(FIXTURES, "1KJL")
        self.assertEqual(metadata["isomeric_smiles"], verified_metadata["isomeric_smiles"])
        self.assertEqual(len(atom_map), 26)
        self.assertEqual(Chem.RemoveHs(ideal).GetNumAtoms(), Chem.RemoveHs(verified).GetNumAtoms())

    def test_anomers_keep_composition_and_differ_in_stereochemistry(self):
        alpha, _, alpha_data = _ccd_ligand(FIXTURES, ["GLA"], [])
        beta, _, beta_data = _ccd_ligand(FIXTURES, ["GAL"], [])
        self.assertEqual(alpha_data["heavy_atoms"], beta_data["heavy_atoms"])
        self.assertEqual(alpha_data["formal_charge"], beta_data["formal_charge"])
        self.assertNotEqual(alpha_data["isomeric_smiles"], beta_data["isomeric_smiles"])
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(alpha), isomericSmiles=False),
                         Chem.MolToSmiles(Chem.RemoveHs(beta), isomericSmiles=False))

    def test_duplicate_and_alternate_atom_identities_are_rejected(self):
        line = (FIXTURES / "154L_l_b.pdb").read_text().splitlines()[0]
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            pdb_heavy_atoms(line+"\n"+line)
        with self.assertRaisesRegex(ValueError, "Alternate"):
            pdb_heavy_atoms(line[:16]+"A"+line[17:])

    def test_source_integrity_and_unknown_cases_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "154L_r_b.pdb"
            path.write_text("mutated")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                verify_external(directory, "154L")
        with self.assertRaisesRegex(ValueError, "Unknown"):
            case_definition("INVALID")

    def test_archived_sources_replay_with_network_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch("aav_benchmark.external.urllib.request.urlopen", side_effect=AssertionError("network requested")):
                files = fetch_external(directory, "154L")
            self.assertIn("154L_r_b.pdb", files)
            self.assertIn("NAG.cif", files)
            verify_external(directory, "154L")

    def test_failed_rerun_invalidates_stale_preparation_and_preserves_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "processed"
            output.mkdir()
            (output / "preparation.json").write_text('{"status":"prepared"}')
            with self.assertRaises(FileNotFoundError):
                prepare(root / "missing_raw", output, "154L")
            self.assertFalse((output / "preparation.json").exists())
            attempt = json.loads((output / "preparation_attempt.json").read_text())
            self.assertEqual(attempt["status"], "failed")
            self.assertEqual(attempt["error_type"], "FileNotFoundError")

    def test_native_calcium_retains_identity_coordinates_and_explicit_plus_two_charge(self):
        from aav_benchmark.aav import _pdb_atoms

        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / "raw"
            fetch_external(raw, "1PWB")
            text = (raw / "1PWB_r_b.pdb").read_text()
            original = pdb_heavy_atoms(text)
            stem = Path(directory) / "prepared"
            log = prepare_external_receptor(text, stem, case_definition("1PWB"), raw)
            output = _pdb_atoms(stem.with_suffix(".pdbqt").read_text(), pdbqt=True)
            heavy_keys = {key for key, atom in output.items() if not atom["hydrogen"]}
            self.assertEqual(heavy_keys, {tuple(atom["key"]) for atom in original})
            calcium = next(atom for atom in original if atom["element"] == "CA")
            key = tuple(calcium["key"])
            np.testing.assert_allclose(output[key]["coordinates"], calcium["coordinates"], atol=0.0001)
            self.assertEqual(output[key]["element_or_type"], "Ca")
            self.assertEqual(log["native_cofactors"][0]["formal_charge"], 2)
            self.assertEqual(log["native_cofactors"][0]["pdbqt_charge"], 2.0)
            self.assertEqual(log["input_heavy_atoms"], log["pdbqt_heavy_atoms"])
            self.assertFalse(log["heavy_atom_completion"])

    def test_uncurated_metal_substitution_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / "raw"
            fetch_external(raw, "1PWB")
            lines = []
            for line in (raw / "1PWB_r_b.pdb").read_text().splitlines():
                if line.startswith(("ATOM", "HETATM")) and line[17:20] == "CA2":
                    line = line[:76]+"MG"+line[78:]
                lines.append(line)
            with self.assertRaisesRegex(ValueError, "Unexpected required-calcium"):
                prepare_external_receptor("\n".join(lines)+"\n", Path(directory)/"prepared",
                                          case_definition("1PWB"), raw)

    def test_deposited_numbering_exposes_missing_loop_without_atom_or_bond_invention(self):
        from aav_benchmark.aav import _pdb_atoms

        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / "raw"
            fetch_external(raw, "2VUZ")
            text = (raw / "2VUZ_r_b.pdb").read_text()
            original = pdb_heavy_atoms(text)
            stem = Path(directory) / "prepared"
            log = prepare_external_receptor(text, stem, case_definition("2VUZ"), raw)
            output = _pdb_atoms(stem.with_suffix(".pdbqt").read_text(), pdbqt=True)
            native_keys = {tuple(atom["key"]) for atom in original}
            self.assertEqual(native_keys, {key for key, atom in output.items() if not atom["hydrogen"]})
            self.assertEqual(log["input_heavy_atoms"], 990)
            self.assertEqual(log["pdbqt_heavy_atoms"], 990)
            self.assertEqual(log["max_heavy_displacement_A"], 0.0)
            self.assertEqual(log["receptor_numbering_audit"]["missing_author_residues"], [91])
            self.assertIn(["A:90", "A:92"], log["topology"]["unresolved_sequence_gaps"])
            self.assertFalse(any(set(bond["residues"]) == {"A:90", "A:92"}
                                 for bond in log["topology"]["accepted_inter_residue_bonds"]))
            # These are the released identities after the explicit output map,
            # not silently renumbered receptor atoms used for pose evaluation.
            self.assertIn(("A", 91, "", "ASN", "N"), output)
            self.assertFalse(log["heavy_atom_completion"])

    def test_assay_plan_predeclares_no_pose_or_binary_endpoint(self):
        frozen = cohort()["gal3_panel_frozen_analysis"]
        self.assertIn("pose_RMSD", frozen["excluded_endpoints"])
        self.assertIn("binary_AUROC", frozen["excluded_endpoints"])
        self.assertEqual(set(frozen["reducing_end_forms"]), {"LacNAc", "lactose", "galactose"})
        self.assertTrue(all(set(forms) == {"alpha", "beta"}
                            for forms in frozen["reducing_end_forms"].values()))


if __name__ == "__main__":
    unittest.main()
