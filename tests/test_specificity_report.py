"""Fixed-denominator and evidence-tampering tests for the Gal9C report."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from aav_benchmark.specificity import _load_panel
from aav_benchmark.specificity_report import generate_specificity_report, write_specificity_report


REPOSITORY = Path(__file__).resolve().parents[1]
SEEDS = [104729, 130363, 155921]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def atom_line(name, index, element):
    kind = {"O": "OA", "N": "NA"}.get(element, element)
    return (f"HETATM{index:5d} {name:>4s} LIG L   1    {float(index):8.3f}{0.:8.3f}{0.:8.3f}"
            f"  1.00  0.00    {0.:6.3f} {kind:>2s}\n")


class SpecificityReportEvidence(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        for name in ("specificity_resolution.json", "gal9_experiment.json", "gal9_receptor_audit.json", "vina.lock.json"):
            (self.root / "benchmarks").mkdir(exist_ok=True)
            shutil.copyfile(REPOSITORY / "benchmarks" / name, self.root / "benchmarks" / name)
        source = self.root / "data/archive/specificity/3NV3.cif.gz"
        source.parent.mkdir(parents=True)
        shutil.copyfile(REPOSITORY / "data/archive/specificity/3NV3.cif.gz", source)
        self.prepared = self.root / "data/processed/gal9_panel"
        self.results = self.root / "results/gal9_specificity"
        panel, rows, _, graphs = _load_panel(self.root / "benchmarks/specificity_resolution.json")
        self.rows = rows
        config_path = self.root / "benchmarks/gal9_experiment.json"
        config = json.loads(config_path.read_text())
        panel_hash = digest(self.root / "benchmarks/specificity_resolution.json")
        audit_hash = digest(self.root / "benchmarks/gal9_receptor_audit.json")
        config_hash = digest(config_path)
        executable_hash = json.loads((self.root / "benchmarks/vina.lock.json").read_text())["assets"]["mac_aarch64"]["sha256"]
        save(self.results / "freeze_receipt.json", {"frozen_before_new_docking": True,
             "experiment_sha256": config_hash, "panel_sha256": panel_hash,
             "receptor_audit_sha256": audit_hash, "source_receptor_sha256": config["source_receptor_sha256"]})
        receptor = {"id": "holo_complete", "pdb_file": "holo_complete.pdb", "pdbqt_file": "holo_complete.pdbqt",
                    "atom_changes_file": "holo_complete_atom_changes.json"}
        for index, row in enumerate(rows):
            case = row["case_id"]
            prepared, results = self.prepared / case, self.results / case
            prepared.mkdir(parents=True)
            results.mkdir(parents=True)
            atoms = graphs[row["ligand_id"]].GetAtoms()
            names = [f"L{i + 1:03d}" for i in range(len(atoms))]
            atom_map = [{"name": name, "rdkit_index": i, "source_smiles_atom_index": i,
                         "element": atom.GetSymbol(), "formal_charge": atom.GetFormalCharge()}
                        for i, (name, atom) in enumerate(zip(names, atoms))]
            ligand = "".join(atom_line(name, i + 1, atom.GetSymbol()) for i, (name, atom) in enumerate(zip(names, atoms)))
            (prepared / "ligand.pdbqt").write_text("ROOT\n" + ligand + "ENDROOT\nTORSDOF 0\n")
            # These fixtures test report provenance and atom-map parsing; they
            # are deliberately synthetic, never scientific docking evidence.
            (prepared / "ligand.sdf").write_text("synthetic fixture conformer\n")
            (prepared / "holo_complete.pdb").write_text("synthetic fixed receptor\n")
            (prepared / "holo_complete.pdbqt").write_text("synthetic fixed receptor PDBQT\n")
            save(prepared / "holo_complete_atom_changes.json", {"synthetic": True})
            receptor_files = {name: digest(prepared / name) for name in
                              ("holo_complete.pdb", "holo_complete.pdbqt", "holo_complete_atom_changes.json")}
            receptor_prep = {"source_receptor_sha256": receptor_files["holo_complete.pdb"],
                             "source_provenance": {"source_cif_sha256": config["source_receptor_sha256"],
                                                   "audit_sha256": audit_hash}, "receptors": [receptor]}
            save(prepared / "receptor_preparation.json", receptor_prep)
            save(prepared / "reference.json", {"names": names, "atom_map": atom_map, "pose_reference_eligible": False})
            save(prepared / "source_graph.json", {"source": panel["source"], "source_row": row,
                 "atom_map": atom_map, "panel_input_sha256": panel_hash})
            preparation = {"case_id": case, "ligand_id": row["ligand_id"], "assay_row": row,
                           "status": "prepared", "pose_reference_eligible": False,
                           "frozen_analysis": panel["frozen_analysis"], "box_center": config["box_center"],
                           "box_size": config["box_size"], "input_fingerprint": {"panel_sha256": panel_hash,
                           "receptor_preparation_sha256": digest(prepared / "receptor_preparation.json"),
                           "receptor_files_sha256": receptor_files}, "receptors": [receptor],
                           "ligand_file": "ligand.pdbqt", "ligand_sha256": digest(prepared / "ligand.pdbqt"),
                           "ligand_sdf_sha256": digest(prepared / "ligand.sdf"), "chemistry": {
                           "heavy_atoms": len(names), "formal_charge": row["formal_charge"],
                           "source_isomeric_smiles": row["isomeric_smiles"], "embedding_seed": config["embedding_seed"],
                           "native_coordinates_used_for_embedding": False, "native_pose_reference_available": False,
                           "stereochemistry_verified_from_3d": True}}
            save(prepared / "preparation.json", preparation)
            shutil.copyfile(config_path, results / "configuration.json")
            runs = []
            pair = index // 2
            base_score = -7. + pair * .2 + (.1 if index % 2 and pair else 0.)
            for seed_index, seed in enumerate(SEEDS):
                score = round(base_score + .01 * seed_index, 3)
                pose_name, log_name = f"holo_complete_{seed}.pdbqt", f"holo_complete_{seed}.log"
                (results / pose_name).write_text(f"MODEL 1\nREMARK VINA RESULT: {score:.3f} 0 0\n" + ligand + "ENDMDL\n")
                center = " ".join(f"{axis} {value:.4f}" for axis, value in zip("XYZ", config["box_center"]))
                (results / log_name).write_text("AutoDock Vina v1.2.7\nScoring function : vina\nExhaustiveness: 16\nCPU: 2\n"
                    + f"Grid center: {center}\nGrid size: X 24 Y 24 Z 24\nPerforming docking (random seed: {seed})\n"
                    + f"1 {score:.3f} 0 0\n")
                command = ["/frozen/vina-1.2.7", "--receptor", f"/frozen/{case}/holo_complete.pdbqt",
                           "--ligand", f"/frozen/{case}/ligand.pdbqt", "--out", f"/frozen/{case}/{pose_name}",
                           "--scoring", "vina", "--cpu", "2", "--seed", str(seed), "--exhaustiveness", "16",
                           "--num_modes", "20", "--min_rmsd", str(config["min_rmsd_A"]),
                           "--energy_range", str(config["energy_range_kcal_mol"])]
                for axis, center, size in zip("xyz", config["box_center"], config["box_size"]):
                    command.extend(["--center_" + axis, str(center), "--size_" + axis, str(size)])
                fingerprint = {"config_sha256": config_hash, "executable_sha256": executable_hash,
                    "ligand_sha256": digest(prepared / "ligand.pdbqt"), "preparation_sha256": digest(prepared / "preparation.json"),
                    "reference_sha256": digest(prepared / "reference.json"), "full_assembly_sha256": None,
                    "receptor_sha256": receptor_files["holo_complete.pdbqt"], "receptor_pdb_sha256": receptor_files["holo_complete.pdb"]}
                runs.append({"case_id": case, "receptor": "holo_complete", "seed": seed, "status": "ok",
                             "vina_version": "AutoDock Vina v1.2.7", "fingerprint": fingerprint,
                             "command": command, "pose_file": pose_name, "log_file": log_name,
                             "pose_sha256": digest(results / pose_name), "returned_poses": 1, "exit_code": 0})
            save(results / "runs.json", runs)

    def mutate_runs(self, update):
        path = self.results / "gal9c_1/runs.json"
        runs = json.loads(path.read_text())
        update(runs)
        save(path, runs)

    def test_complete_raw_recomputation_ties_baselines_and_endpoint_scope(self):
        summary, markdown = generate_specificity_report(self.root)
        self.assertEqual(summary["denominator"]["successful_runs"], 36)
        self.assertTrue(summary["denominator"]["metrics_eligible"])
        self.assertEqual(summary["matched_pair_counts"]["exact_aggregate_ties"], 1)
        self.assertEqual(summary["matched_pair_counts"]["measured_binder_lower_energy"], 5)
        self.assertAlmostEqual(summary["matched_pair_counts"]["tie_aware_concordance"], 5.5 / 6)
        for baseline in summary["metrics"]["baselines"].values():
            self.assertEqual(baseline["contexts"][0]["raw"]["auroc"], .5)
        quant = summary["metrics"]["positive_only_quantitative_association"]["contexts"][0]
        self.assertEqual(quant["n_ligands"], 6)
        self.assertEqual(len(quant["ligands"]), 6)
        for row in summary["compounds"]:
            self.assertIsNone(row["numeric_kd_lower_bound_uM"])
            if not row["binary_label"]:
                self.assertIsNone(row["kd_uM"])
        self.assertNotIn("top1_rmsd", json.dumps(summary))
        self.assertIn("NB at 2 mM", markdown)
        self.assertEqual(summary, generate_specificity_report(self.root)[0])

    def test_missing_and_failed_runs_keep_36_denominator_and_disable_all_metrics(self):
        self.mutate_runs(lambda rows: rows.pop())
        summary, _ = generate_specificity_report(self.root)
        self.assertEqual(summary["denominator"]["expected_runs"], 36)
        self.assertEqual(summary["denominator"]["missing_runs"], 1)
        self.assertIsNone(summary["metrics"])
        self.assertIsNone(summary["compounds"][0]["mean_top1_score_kcal_mol"])
        self.mutate_runs(lambda rows: rows[0].update({"status": "failed", "exit_code": 1, "error": "fixture failure"}))
        summary, _ = generate_specificity_report(self.root)
        self.assertEqual(summary["denominator"]["failed_runs"], 1)
        self.assertEqual(summary["denominator"]["complete_ligands"], 11)

    def test_duplicate_or_unexpected_run_cannot_shrink_or_expand_denominator(self):
        self.mutate_runs(lambda rows: rows.append(dict(rows[0])))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            generate_specificity_report(self.root)

    def test_changed_input_fingerprint_and_cpu_command_are_rejected(self):
        self.mutate_runs(lambda rows: rows[0]["fingerprint"].update({"reference_sha256": "0" * 64}))
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            generate_specificity_report(self.root)
        reference = digest(self.prepared / "gal9c_1/reference.json")
        self.mutate_runs(lambda rows: rows[0]["fingerprint"].update({"reference_sha256": reference}))
        self.mutate_runs(lambda rows: rows[0]["command"].__setitem__(rows[0]["command"].index("--cpu") + 1, "3"))
        with self.assertRaisesRegex(ValueError, "--cpu"):
            generate_specificity_report(self.root)

    def test_changed_pose_hash_and_invalid_atom_map_are_rejected(self):
        path = self.results / "gal9c_1/holo_complete_104729.pdbqt"
        path.write_text(path.read_text().replace("L001", "FAIL"))
        with self.assertRaisesRegex(ValueError, "pose hash"):
            generate_specificity_report(self.root)
        self.mutate_runs(lambda rows: rows[0].update({"pose_sha256": digest(path)}))
        with self.assertRaisesRegex(ValueError, "atom identity mismatch"):
            generate_specificity_report(self.root)

    def test_raw_log_score_table_must_match_pose_scores(self):
        path = self.results / "gal9c_1/holo_complete_104729.log"
        path.write_text(path.read_text().replace("1 -7.000 0 0", "1 -7.100 0 0"))
        with self.assertRaisesRegex(ValueError, "score table"):
            generate_specificity_report(self.root)

    def test_prepared_ligand_and_frozen_receipt_changes_are_rejected(self):
        path = self.prepared / "gal9c_1/ligand.sdf"
        original = path.read_bytes()
        path.write_bytes(original + b"tampered")
        with self.assertRaisesRegex(ValueError, "ligand hash"):
            generate_specificity_report(self.root)
        path.write_bytes(original)
        receipt = self.results / "freeze_receipt.json"
        document = json.loads(receipt.read_text())
        document["experiment_sha256"] = "0" * 64
        save(receipt, document)
        with self.assertRaisesRegex(ValueError, "freeze receipt"):
            generate_specificity_report(self.root)

    def test_archived_gzip_inputs_recompute_identical_metrics_and_commands(self):
        original, _ = generate_specificity_report(self.root)
        for directory in (self.prepared, self.results):
            for path in list(directory.rglob("*")):
                if not path.is_file():
                    continue
                path.with_name(path.name + ".gz").write_bytes(gzip.compress(path.read_bytes(), mtime=0))
                path.unlink()
        archived, _ = generate_specificity_report(self.root, prepared_root=self.prepared, results_root=self.results)
        self.assertEqual(original["metrics"], archived["metrics"])
        self.assertEqual(original["runs"], archived["runs"])
        self.assertEqual(original["compounds"], archived["compounds"])
        self.assertIn("gzip_sha256", archived["source_hash_manifest"]["verified_artifacts"]["prepared/gal9c_1/ligand.pdbqt"])

    def test_write_and_check_detects_stale_output(self):
        paths = write_specificity_report(self.root)
        self.assertEqual(paths, write_specificity_report(self.root, check=True))
        paths[0].write_text("stale\n")
        with self.assertRaisesRegex(ValueError, "stale"):
            write_specificity_report(self.root, check=True)

    def test_portable_preparation_bundles_are_default_and_validate_compressed_bytes(self):
        from aav_benchmark.archive import bundle
        original, _ = generate_specificity_report(self.root)
        for directory in self.prepared.iterdir():
            bundle(directory, self.results / directory.name)
        shutil.rmtree(self.prepared)
        archived, _ = generate_specificity_report(self.root)
        self.assertEqual(original["metrics"], archived["metrics"])
        self.assertEqual(original["runs"], archived["runs"])
        corrupted = self.results / "gal9c_1/inputs/ligand.pdbqt.gz"
        corrupted.write_bytes(corrupted.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "archive bytes changed"):
            generate_specificity_report(self.root)

    def test_explicit_float_tolerance_accepts_only_small_derived_correlation_noise(self):
        import contextlib
        import io
        paths = write_specificity_report(self.root)
        original = json.loads(paths[0].read_text())
        original["metrics"]["positive_only_quantitative_association"]["contexts"][0]["spearman_rho"] += 1e-12
        save(paths[0], original)
        with self.assertRaisesRegex(ValueError, "first difference: .*spearman_rho"):
            write_specificity_report(self.root, check=True)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(paths, write_specificity_report(self.root, check=True, numeric_tolerance=1e-10))
        self.assertIn("maximum difference", output.getvalue())
        self.assertIn("explicit tolerance 1e-10", output.getvalue())

    def test_float_tolerance_rejects_copied_data_counts_hashes_large_drift_and_markdown(self):
        paths = write_specificity_report(self.root)
        original = paths[0].read_text()
        changes = [
            lambda x: x["compounds"][0].update(kd_uM=x["compounds"][0]["kd_uM"] + 1),
            lambda x: x["runs"][0].update(top1_score_kcal_mol=x["runs"][0]["top1_score_kcal_mol"] + 1e-12),
            lambda x: x["denominator"].update(successful_runs=35),
            lambda x: x["runs"][0].update(pose_sha256="0" * 64),
            lambda x: x["metrics"]["positive_only_quantitative_association"]["contexts"][0].update(spearman_rho=.12345),
        ]
        for change in changes:
            with self.subTest(change=change):
                document = json.loads(original)
                change(document)
                save(paths[0], document)
                with self.assertRaisesRegex(ValueError, "stale"):
                    write_specificity_report(self.root, check=True, numeric_tolerance=1e-10)
        paths[0].write_text(original)
        paths[1].write_text(paths[1].read_text() + "changed\n")
        with self.assertRaisesRegex(ValueError, "stale"):
            write_specificity_report(self.root, check=True, numeric_tolerance=1e-10)

    def test_float_tolerance_must_be_explicit_finite_and_bounded(self):
        for tolerance in (-1., float("nan"), float("inf"), True, 1.1e-8):
            with self.subTest(tolerance=tolerance), self.assertRaisesRegex(ValueError, "Numeric tolerance"):
                write_specificity_report(self.root, check=True, numeric_tolerance=tolerance)


if __name__ == "__main__":
    unittest.main()
