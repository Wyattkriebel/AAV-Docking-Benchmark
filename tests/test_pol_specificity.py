"""POL source chemistry, complete assemblies and freeze gates; no scientific docking."""
import copy
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem

from aav_benchmark.external import pdb_heavy_atoms
from aav_benchmark.pol_specificity import (
    ARCHIVE, CA_SHA256, COHORT_ID, LEDGER, MODEL_SPECS, ROOT,
    _ligand_models, dock_pol_panel, experiment_template, load_source_qualification,
    planned_cases, preflight, prepare_pol_panel, receptor_sources, summarize_pol_scores, validate_experiment,
)
from aav_benchmark.provenance import sha256, write_json


def source_copy(destination):
    shutil.copytree(ARCHIVE, destination)
    ca = destination / "CA.cif"
    if not ca.exists():
        ca.write_bytes(gzip.decompress((ROOT / "data/archive/external/CA.cif.gz").read_bytes()))
    if sha256(ca) != CA_SHA256:
        raise ValueError("test requires the actual frozen calcium CCD")


class POLSourceControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sources = self.root / "sources"
        source_copy(self.sources)
        self.ledger = self.root / "ledger.json"
        self.ledger.write_bytes(LEDGER.read_bytes())
        self.config = experiment_template(self.ledger, self.sources)

    def test_two_assay_species_four_exact_model_graphs_and_full_run_denominator(self):
        source, graphs = load_source_qualification(self.ledger, self.sources)
        self.assertEqual(source["candidate_id"], COHORT_ID)
        self.assertEqual(set(graphs), {"GLA", "GAL", "GLC", "BGC"})
        cases = planned_cases(self.config)
        self.assertEqual(len(cases), 16)
        self.assertEqual(len({r["case_id"] for r in cases}), 16)
        self.assertEqual(len(cases) * len(self.config["metal_profiles"]) * len(self.config["seeds"]), 96)
        self.assertEqual({r["species_id"] for r in cases}, {"D_galactose", "D_glucose"})

    def test_unfrozen_experiment_blocks_execution_before_engine_or_artifact_access(self):
        path = self.root / "draft.json"
        write_json(path, self.config)
        with patch("aav_benchmark.engine.dock") as engine:
            with self.assertRaisesRegex(ValueError, "separate frozen experiment"):
                dock_pol_panel(path, self.root / "nonexistent", self.root / "unused", self.root / "vina",
                               ledger_path=self.ledger, source_dir=self.sources)
        engine.assert_not_called()

    def test_changed_inventory_seeds_boxes_and_primary_are_rejected(self):
        for key, value in [("seeds", self.config["seeds"][:1]), ("model_ids", self.config["model_ids"][:2]),
                           ("box_size", [20, 20, 20]), ("expected_runs", 12),
                           ("primary_condition", {"state": "apo", "site": "C", "metal_profile": "dry"})]:
            altered = {**self.config, key: value}
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, key):
                validate_experiment(altered, self.ledger, self.sources)

    def test_mutated_raw_source_or_calcium_ccd_is_rejected(self):
        for filename in ["GLA.cif", "6T0Q.cif", "CA.cif"]:
            path = self.sources / filename
            original = path.read_bytes()
            path.write_bytes(original + b"\n# changed\n")
            with self.subTest(source=filename), self.assertRaisesRegex(ValueError, "checksum|pinned"):
                validate_experiment(self.config, self.ledger, self.sources)
            path.write_bytes(original)

    def test_numeric_glucose_kd_and_source_graph_alias_change_fail(self):
        document = json.loads(self.ledger.read_text())
        pol = next(row for row in document["candidates"] if row["candidate_id"] == COHORT_ID)
        glucose = next(row for row in pol["source_tested_inventory"] if row["species_id"] == "D_glucose")
        glucose["Kd_M"] = 0.01
        write_json(self.ledger, document)
        with self.assertRaisesRegex(ValueError, "invented a glucose Kd"):
            load_source_qualification(self.ledger, self.sources)
        self.ledger.write_bytes(LEDGER.read_bytes())
        _, graphs = load_source_qualification(self.ledger, self.sources)
        altered = copy.deepcopy(graphs)
        altered["GLA"]["canonical_smiles"] = altered["GLC"]["canonical_smiles"]
        with self.assertRaisesRegex(ValueError, "exact source graph"):
            _ligand_models(self.sources, altered)

    def test_actual_assemblies_use_whole_residue_altloc_and_backbone_only_mapping(self):
        receptors, centers, alignment = receptor_sources(self.sources)
        for state, expected in [("holo", 5034), ("apo", 5058)]:
            text, audit = receptors[state]
            atoms = pdb_heavy_atoms(text)
            protein = [a for a in atoms if a["key"][3] != "CA2"]
            self.assertEqual(len(protein), expected)
            self.assertEqual(len(atoms) - len(protein), 4)
            self.assertEqual(audit["materializations"], 1)
        apo = receptors["apo"][1]
        self.assertEqual(apo["raw_protein_heavy_atoms_per_chain_before_altloc_selection"], {"A": 2532, "B": 2532})
        self.assertEqual(len(apo["alternate_conformers"]), 2)
        for row in apo["alternate_conformers"]:
            self.assertEqual(row["selected"], "B")
            self.assertEqual(set(row["retained_selected_atom_names"]), {"CA", "CB", "OG"})
            self.assertEqual(set(row["discarded_atom_names"]), {"CA", "CB", "OG"})
            self.assertFalse(row["occupancy_normalized"])
        self.assertEqual(alignment["matched_backbone_atoms"], 2704)
        self.assertFalse(alignment["native_ligand_coordinates_used_for_alignment"])
        self.assertFalse(alignment["receptor_coordinates_transformed"])
        self.assertEqual(len(alignment["candidate_chain_correspondences"]), 2)
        self.assertEqual(alignment["backbone_rmsd_A"], min(r["backbone_rmsd_A"] for r in alignment["candidate_chain_correspondences"]))
        rotation, translation = np.array(alignment["rotation_row_vector"]), np.array(alignment["translation_A"])
        for site in ("N", "C"):
            np.testing.assert_allclose(np.array(centers["holo"][site]) @ rotation + translation, centers["apo"][site])


class POLActualPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.sources = cls.root / "sources"
        source_copy(cls.sources)
        cls.config = experiment_template(source_dir=cls.sources)
        cls.config["analysis_frozen_before_new_docking"] = True
        cls.config_path = cls.root / "experiment.json"
        write_json(cls.config_path, cls.config)
        cls.output = cls.root / "prepared"
        cls.panel = prepare_pol_panel(cls.config_path, cls.output, source_dir=cls.sources)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_every_prepared_ligand_has_its_exact_source_graph_and_no_pose_truth(self):
        self.assertEqual(len(self.panel["cases"]), 16)
        _, graphs = load_source_qualification(source_dir=self.sources)
        for row in self.panel["cases"]:
            directory = self.output / row["prepared_dir"]
            molecule = Chem.SDMolSupplier(str(directory / "ligand.sdf"), removeHs=False)[0]
            self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True), graphs[row["ccd_id"]]["canonical_smiles"])
            self.assertEqual(molecule.GetNumHeavyAtoms(), 12)
            reference = json.loads((directory / "reference.json").read_text())
            metadata = json.loads((directory / "preparation.json").read_text())
            self.assertFalse(reference["pose_reference_eligible"])
            self.assertFalse(metadata["experimental_anomer_label_available"])
            self.assertEqual(metadata["box_size"], [24, 24, 24])
            self.assertEqual(metadata["chemistry"]["embedding_seed"], 20261024)
            self.assertFalse(metadata["chemistry"]["native_coordinates_used_for_embedding"])
            np.testing.assert_allclose(np.mean(reference["coordinates"], axis=0), metadata["box_center"], atol=1e-10)

    def test_ca_profiles_have_exact_ions_and_identical_complete_protein_coordinates(self):
        for state, expected in [("holo", 5034), ("apo", 5058)]:
            wet = pdb_heavy_atoms((self.output / "receptors" / f"{state}_Ca_retained.pdb").read_text())
            dry = pdb_heavy_atoms((self.output / "receptors" / f"{state}_dry.pdb").read_text())
            ions = [a for a in wet if a["key"][3] == "CA2"]
            self.assertEqual(len(ions), 4)
            self.assertTrue(all(a["key"][4] == "CA+2" for a in ions))
            self.assertEqual([a for a in wet if a["key"][3] != "CA2"], dry)
            self.assertEqual(len(dry), expected)
            log = json.loads((self.output / "receptors" / f"{state}_Ca_retained_atom_changes.json").read_text())
            self.assertEqual(len(log["native_cofactors"]), 4)
            self.assertTrue(all(r["pdbqt_charge"] == 2 and r["autodock_atom_type"] == "Ca" for r in log["native_cofactors"]))
            self.assertLess(log["max_heavy_displacement_A"], 0.002)

    def test_wrapper_preserves_96_run_design_and_explicit_execution_subsets(self):
        def fake_engine(*args, **kwargs):
            return [{"status": "ok"} for _ in range(2 * len(kwargs.get("seeds") or self.config["seeds"]))]
        with patch("aav_benchmark.engine.dock", side_effect=fake_engine) as engine:
            rows = dock_pol_panel(self.config_path, self.output, self.root / "mock_execution", self.root / "vina",
                                  source_dir=self.sources)
            self.assertEqual(engine.call_count, 16)
            self.assertEqual(sum(r["runs"] for r in rows), 96)
        with patch("aav_benchmark.engine.dock", side_effect=fake_engine) as engine:
            rows = dock_pol_panel(self.config_path, self.output, self.root / "mock_subset", self.root / "vina",
                                  source_dir=self.sources, case_ids=[self.panel["cases"][0]["case_id"]], seeds=[104729])
            self.assertEqual(engine.call_count, 1)
            self.assertEqual(rows[0]["runs"], 2)

    def test_changed_prepared_ligand_blocks_before_engine(self):
        row = self.panel["cases"][0]
        path = self.output / row["prepared_dir"] / "ligand.pdbqt"
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"\nREMARK changed\n")
            with patch("aav_benchmark.engine.dock") as engine:
                with self.assertRaisesRegex(ValueError, "prepared artifact changed"):
                    dock_pol_panel(self.config_path, self.output, self.root / "unused", self.root / "vina",
                                   source_dir=self.sources, case_ids=[row["case_id"]])
            engine.assert_not_called()
        finally:
            path.write_bytes(original)

    def test_preparation_requires_fresh_output(self):
        with self.assertRaisesRegex(ValueError, "fresh"):
            prepare_pol_panel(self.config_path, self.output, source_dir=self.sources)


class POLAggregationControls(unittest.TestCase):
    def setUp(self):
        self.rows = []
        for case in planned_cases({}):
            for profile in ("Ca_retained", "dry"):
                for seed in (104729, 130363, 155921):
                    # Synthetic values exercise aggregation only.
                    value = -6.0 if case["species_id"] == "D_galactose" else -5.0
                    value += 0.5 if case["anomer"] == "beta" else 0
                    self.rows.append({"case_id": case["case_id"], "metal_profile": profile,
                                      "seed": seed, "top1_score": value})

    def test_seed_then_equal_anomer_mean_and_separate_sensitivities(self):
        result = summarize_pol_scores(self.rows)
        primary = result["primary"]
        self.assertEqual(len(result["conditions"]), 8)
        self.assertEqual((primary["state"], primary["site"], primary["metal_profile"]), ("holo", "N", "Ca_retained"))
        self.assertEqual(primary["species"]["D_galactose"]["equal_anomer_mean_kcal_mol"], -5.75)
        self.assertEqual(primary["species"]["D_galactose"]["minimum_anomer_seed_mean_kcal_mol"], -6)
        self.assertEqual(primary["species"]["D_galactose"]["maximum_anomer_seed_mean_kcal_mol"], -5.5)
        self.assertEqual(primary["descriptive_single_pair_concordance"], 1)
        self.assertEqual(result["independent_families"], 1)
        self.assertEqual(result["experimental_epimer_contrasts"], 1)
        self.assertIsNone(result["auroc"])
        self.assertFalse(result["affinity_calibration"])

    def test_ties_missing_runs_duplicates_and_nonfinite_scores(self):
        for row in self.rows:
            row["top1_score"] = -5
        self.assertEqual(summarize_pol_scores(self.rows)["primary"]["descriptive_single_pair_concordance"], 0.5)
        for invalid in [self.rows[:-1], self.rows + [self.rows[0]],
                        [{**self.rows[0], "top1_score": float("nan")}, *self.rows[1:]]]:
            with self.assertRaises(ValueError):
                summarize_pol_scores(invalid)


if __name__ == "__main__":
    unittest.main()
