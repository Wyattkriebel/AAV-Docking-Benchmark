"""Evaluation truth controls, using pinned archives without network or docking."""
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.spatial import cKDTree

from aav_benchmark.site_truth import DEFAULT_LEDGER, _contact_residues, load_site_truth


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_COUNTS = {"154L": 43, "1GNY": 46, "1KJL": 26, "1QFO": 43,
                   "2J1T": 46, "1PWB": 23, "2VUZ": 90}


class SiteTruthControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.raw = Path(cls.temporary.name)
        (cls.raw / "external").mkdir()
        for case in EXTERNAL_COUNTS:
            for part in ("r", "l"):
                name = f"{case}_{part}_b.pdb"
                archive = ROOT / "data/archive/external" / f"{name}.gz"
                (cls.raw / "external" / name).write_bytes(gzip.decompress(archive.read_bytes()))
        (cls.raw / "7MTZ.cif").write_bytes(
            gzip.decompress((ROOT / "data/archive/7MTZ.cif.gz").read_bytes()))

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_all_seven_external_truths_keep_every_modeled_atom_and_ledger_identity(self):
        for case, expected_count in EXTERNAL_COUNTS.items():
            with self.subTest(case=case):
                truth = load_site_truth(case.lower(), self.raw)
                self.assertEqual(truth["case_id"], case)
                self.assertEqual(truth["site_class_count"], 1)
                self.assertEqual(len(truth["sites"]), 1)
                site = truth["sites"][0]
                self.assertEqual(len(site["coordinates"]), expected_count)
                self.assertEqual(len(site["atom_keys"]), expected_count)
                self.assertEqual(len(set(site["names"])), expected_count)
                self.assertTrue(site["interface_residue_keys"])
                self.assertFalse(truth["atomic_pose_eligible"])
                self.assertFalse(truth["confirmatory_site"])
                self.assertTrue(truth["evaluation_only"])
                self.assertEqual(truth["coordinate_frame"], "pinned_holo_source")
                self.assertTrue(all(s["sha256"] == s["expected_sha256"] for s in truth["source_files"]))
                json.dumps(truth)  # No ndarray/scalar leakage into the archived JSON contract.

    def test_sialyllactose_site_does_not_require_invalid_full_chemistry_to_pass(self):
        with patch("aav_benchmark.external.build_ligand", side_effect=AssertionError("No chemistry preparation")):
            truth = load_site_truth("1QFO", self.raw)
        site = truth["sites"][0]
        self.assertEqual({r["deposited"]["author_chain"]
                          for r in site["deposited_ligand_residue_map"]}, {"E"})
        self.assertIn(["A", 96, "", "ARG"], site["interface_residue_keys"])
        self.assertIn("chirality", truth["reference_limitations"])

    def test_calcium_is_not_counted_as_a_protein_interface_residue(self):
        truth = load_site_truth("1PWB", self.raw)
        self.assertEqual(truth["native_protein_atom_count"], 1154)
        contacts = truth["sites"][0]["interface_residue_keys"]
        self.assertEqual(len(contacts), 6)
        self.assertNotIn("CA2", {key[3] for key in contacts})

    def test_one_tree_full_capsid_truth_keeps_neighbor_contacts_and_one_site_class(self):
        truth = load_site_truth("7mtz", self.raw)
        self.assertEqual(truth["case_id"], "aav9")
        self.assertEqual(truth["native_protein_atom_count"], 247860)
        self.assertEqual(truth["native_protein_chain_count"], 60)
        self.assertEqual(truth["site_class_count"], 1)
        self.assertEqual(truth["symmetry_copy_count"], 60)
        self.assertEqual(len(truth["sites"]), 60)
        self.assertEqual(len({site["id"] for site in truth["sites"]}), 60)
        self.assertTrue(all(len(site["coordinates"]) == 12 for site in truth["sites"]))
        self.assertTrue(all(len(site["interface_residue_keys"]) == 3 for site in truth["sites"]))
        representative = next(site for site in truth["sites"] if site["id"] == "7MTZ:A:801:GAL")
        self.assertEqual(representative["deposited_label_chain"], "IB")
        self.assertEqual({tuple(key) for key in representative["interface_residue_keys"]},
                         {("A", 503, "", "TRP"), ("I", 470, "", "ASN"), ("I", 473, "", "VAL")})
        self.assertFalse(truth["apo_coordinate_transform_applied"])
        json.dumps(truth)

    def test_contacts_use_strict_boundary_and_preserve_insertion_codes(self):
        xyz = np.asarray([[3.999, 0, 0], [4, 0, 0], [0, 0, 0], [3, 0, 0]], float)
        keys = [["A", 9, "B", "ASN"], ["A", 10, "", "GLU"],
                ["C", 9, "", "TYR"], ["A", 9, "B", "ASN"]]
        contacts = _contact_residues(cKDTree(xyz), xyz, keys, [[0, 0, 0]])
        self.assertEqual(contacts, [["A", 9, "B", "ASN"], ["C", 9, "", "TYR"]])

    def test_changed_native_source_bytes_are_refused_before_parsing(self):
        path = self.raw / "external/154L_l_b.pdb"
        payload = path.read_bytes()
        try:
            path.write_bytes(payload + b"REMARK changed source\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                load_site_truth("154L", self.raw)
        finally:
            path.write_bytes(payload)

    def test_frozen_ledger_atom_loss_identity_or_contacts_cannot_be_silent(self):
        for alteration in ("atom_count", "ligand_identity", "contact_identity"):
            with self.subTest(alteration=alteration):
                ledger = json.loads(DEFAULT_LEDGER.read_text())
                case = next(r for r in ledger["cases"] if r["case_id"] == "154L")
                evaluation = case["evaluation_only"]
                if alteration == "atom_count":
                    evaluation["ligand_atoms"] -= 1
                elif alteration == "ligand_identity":
                    evaluation["ligand_residue_map"][0]["prepared_resnum"] = 999
                else:
                    evaluation["native_contact_residues_4A"][0]["prepared_residue_number"] = 999
                path = self.raw / "altered_ledger.json"
                path.write_text(json.dumps(ledger))
                with self.assertRaises(ValueError):
                    load_site_truth("154L", self.raw, ledger_path=path)

    def test_custom_external_directory_and_unknown_case(self):
        truth = load_site_truth("154L", "/unused", external_raw_dir=self.raw / "external")
        self.assertEqual(len(truth["sites"][0]["coordinates"]), 43)
        with self.assertRaisesRegex(ValueError, "Unknown"):
            load_site_truth("not_a_case", self.raw)


if __name__ == "__main__":
    unittest.main()
