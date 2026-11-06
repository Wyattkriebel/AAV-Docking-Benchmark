"""Fixed-graph matched epimers, shared inputs, and conditional assay metrics."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from rdkit import Chem

from aav_benchmark.engine import analyze
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.specificity import (
    _load_panel,
    prepare_specificity_panel,
    summarize_specificity_scores,
)


PANEL = Path(__file__).resolve().parent.parent / "benchmarks/specificity_resolution.json"


class SpecificityPreparationControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.receptor = cls.root / "receptor"
        cls.receptor.mkdir()
        cls.source_panel = cls.root / "panel.json"
        cls.source_panel.write_bytes(PANEL.read_bytes())
        (cls.receptor / "holo_complete.pdb").write_text("receptor PDB fixture\n")
        (cls.receptor / "holo_complete.pdbqt").write_text("receptor PDBQT fixture\n")
        write_json(cls.receptor / "holo_complete_atom_changes.json", {"heavy_atom_completion": False})
        write_json(cls.receptor / "preparation.json", {
            "case_id": "3NV3", "target_id": "O00182", "status": "prepared",
            "box_center": [10, 20, 30], "box_size": [25, 25, 25],
            "receptors": [{"id": "holo_complete", "state": "holo", "assembly_condition": "isolated_CRD",
                            "pdb_file": "holo_complete.pdb", "pdbqt_file": "holo_complete.pdbqt",
                            "atom_changes_file": "holo_complete_atom_changes.json"}]})
        cls.input_hashes = {path: sha256(path) for path in [cls.source_panel, *cls.receptor.iterdir()]}
        cls.output = cls.root / "prepared"
        cls.prepared = prepare_specificity_panel(cls.source_panel, cls.receptor, cls.output)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_all_twelve_fixed_stereo_graphs_and_source_inputs_preserved(self):
        self.assertEqual(len(self.prepared), 12)
        self.assertEqual({row["ligand_id"] for row in self.prepared},
                         {"1", "14", "4", "17", "5", "18", "6", "19", "7", "20", "8", "21"})
        self.assertEqual({path: sha256(path) for path in self.input_hashes}, self.input_hashes)
        for row in self.prepared:
            with self.subTest(ligand=row["ligand_id"]):
                path = self.output / row["case_id"]
                ligand = Chem.SDMolSupplier(str(path / "ligand.sdf"), removeHs=False)[0]
                self.assertIsNotNone(ligand)
                heavy = Chem.RemoveHs(ligand)
                self.assertEqual(Chem.MolToSmiles(heavy, isomericSmiles=True), row["assay_row"]["isomeric_smiles"])
                self.assertEqual(heavy.GetNumAtoms(), row["assay_row"]["expected_heavy_atoms"])
                self.assertEqual(Chem.GetFormalCharge(heavy), row["assay_row"]["formal_charge"])
                self.assertEqual(row["chemistry"]["embedding_seed"], 20261024)
                self.assertTrue(row["chemistry"]["stereochemistry_verified_from_3d"])
                self.assertFalse(row["chemistry"]["native_coordinates_used_for_embedding"])

    def test_same_receptor_box_and_ineligible_pose_reference_for_every_alternative(self):
        for row in self.prepared:
            path = self.output / row["case_id"]
            self.assertEqual(row["box_center"], [10, 20, 30])
            self.assertEqual(row["box_size"], [25, 25, 25])
            for filename in ("holo_complete.pdb", "holo_complete.pdbqt", "holo_complete_atom_changes.json"):
                self.assertEqual(sha256(path / filename), sha256(self.receptor / filename))
            reference = json.loads((path / "reference.json").read_text())
            self.assertFalse(reference["pose_reference_eligible"])
            self.assertEqual(len(reference["names"]), row["chemistry"]["heavy_atoms"])
            np.testing.assert_allclose(np.mean(reference["coordinates"], axis=0), [10, 20, 30], atol=1e-10)
            with self.assertRaisesRegex(ValueError, "no experimental atomic pose reference"):
                analyze(path, self.root / "unused_results")

    def test_source_graphs_and_pair_ids_are_archived(self):
        for row in self.prepared:
            source = json.loads((self.output / row["case_id"] / "source_graph.json").read_text())
            self.assertEqual(source["source_row"]["pair_id"], row["pair_id"])
            self.assertEqual(source["panel_input_sha256"], sha256(self.source_panel))
            self.assertEqual(len(source["atom_map"]), row["chemistry"]["heavy_atoms"])

    def test_preparation_never_reuses_stale_ligand_directories(self):
        with self.assertRaisesRegex(ValueError, "fresh"):
            prepare_specificity_panel(self.source_panel, self.receptor, self.output)


class SpecificityValidationAndMetrics(unittest.TestCase):
    def setUp(self):
        self.document = json.loads(PANEL.read_text())
        self.panel = self.document["primary_specificity_panel"]
        self.rows = [row for row in self.panel["rows"] if row["primary_cohort"]]
        self.scores = {row["ligand_id"]: -7.0 if row["binary_label"] else -6.0 for row in self.rows}
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name) / "panel.json"

    def save(self):
        write_json(self.path, self.document)
        return self.path

    def test_perfect_matched_direction_and_tied_size_charge_baselines(self):
        result = summarize_specificity_scores(self.save(), self.scores)
        metrics = result["classification"]["contexts"][0]
        self.assertEqual(metrics["raw"]["auroc"], 1)
        self.assertEqual(metrics["raw"]["auprc_average_precision"], 1)
        self.assertEqual(metrics["matched"]["pairs"], 6)
        self.assertEqual(metrics["matched"]["pair_concordance"], 1)
        for baseline in result["baselines"].values():
            self.assertEqual(baseline["contexts"][0]["raw"]["auroc"], 0.5)
            self.assertEqual(baseline["contexts"][0]["matched"]["pair_concordance"], 0.5)
        self.assertEqual(result["positive_only_quantitative_association"]["contexts"][0]["n_ligands"], 6)
        self.assertIsNone(result["right_censored_quantitative_contrast"]["estimate"])

    def test_ties_receive_half_credit_and_scores_are_not_fit_to_labels(self):
        self.scores = {key: -6.0 for key in self.scores}
        result = summarize_specificity_scores(self.save(), self.scores)
        metrics = result["classification"]["contexts"][0]
        self.assertEqual(metrics["raw"]["auroc"], 0.5)
        self.assertEqual(metrics["raw"]["auprc_average_precision"], 0.5)
        self.assertEqual(metrics["matched"]["pair_concordance"], 0.5)

    def test_missing_nonfinite_and_extra_ligand_scores_fail(self):
        for scores in [{key: value for key, value in self.scores.items() if key != "1"},
                       {**self.scores, "1": float("nan")}, {**self.scores, "extra": -5}]:
            with self.subTest(scores=scores), self.assertRaises(ValueError):
                summarize_specificity_scores(self.save(), scores)

    def test_exact_kd_invention_for_nb_and_undefined_stereo_are_rejected(self):
        for field, value in [("kd_uM", 2000), ("numeric_kd_lower_bound_uM", 2000),
                             ("isomeric_smiles", self.rows[1]["isomeric_smiles"].replace("@@", ""))]:
            original = copy.deepcopy(self.document)
            self.rows[1][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                _load_panel(self.save())
            self.document = original
            self.panel = self.document["primary_specificity_panel"]
            self.rows = [row for row in self.panel["rows"] if row["primary_cohort"]]

    def test_missing_cohort_member_duplicate_pair_and_changed_charge_fail(self):
        mutations = [lambda panel: panel["rows"].remove(next(row for row in panel["rows"] if row["ligand_id"] == "14")),
                     lambda panel: panel["matched_pairs"].append(copy.deepcopy(panel["matched_pairs"][0])),
                     lambda panel: next(row for row in panel["rows"] if row["ligand_id"] == "14").update(formal_charge=1)]
        for mutation in mutations:
            original = copy.deepcopy(self.document)
            mutation(self.document["primary_specificity_panel"])
            with self.assertRaises(ValueError):
                _load_panel(self.save())
            self.document = original


if __name__ == "__main__":
    unittest.main()
