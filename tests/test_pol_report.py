"""Synthetic 96-run report controls using real source/preparation input graphs.

All result rows, poses and log scores created here are explicitly synthetic and
live only in TemporaryDirectory fixtures. No scientific docking is performed,
and no actual POL scores are copied into these fixtures.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from aav_benchmark.pol_report import generate_pol_report, write_pol_report
from aav_benchmark.pol_specificity import COHORT_ID, MODEL_SPECS, planned_cases
from aav_benchmark.provenance import sha256


REPO = Path(__file__).resolve().parents[1]
LIVE_RESULTS = REPO / "results/pol_specificity"


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def archived_read(case_dir, filename):
    return gzip.decompress((case_dir / "inputs" / (filename + ".gz")).read_bytes())


def synthetic_fixture(root):
    """Clone source evidence and synthesize outputs; never open real run files."""
    result = root / "results/pol_specificity"
    result.mkdir(parents=True)
    for name in ("pol_experiment.json", "independent_specificity_candidates.json", "vina.lock.json"):
        target = root / "benchmarks" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / "benchmarks" / name, target)
    shutil.copytree(REPO / "data/archive/independent_specificity_sources",
                    root / "data/archive/independent_specificity_sources")
    for name in ("freeze_receipt.json", "proposal_before_freeze.json", "pol_experiment.json",
                 "panel_preparation.json", "execution_snapshot.json", "source_preflight.json"):
        shutil.copyfile(LIVE_RESULTS / name, result / name)
    shutil.copytree(LIVE_RESULTS / "execution_snapshot", result / "execution_snapshot")
    config = json.loads((root / "benchmarks/pol_experiment.json").read_text())
    for name in config["implementation_files_sha256"]:
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / name, target)
    panel = json.loads((result / "panel_preparation.json").read_text())
    frozen_sha = sha256(root / "benchmarks/pol_experiment.json")
    executable_sha = json.loads((result / "execution_snapshot.json").read_text())["executable_sha256"]
    for case in panel["cases"]:
        case_id = case["case_id"]
        directory = result / case_id
        shutil.copytree(LIVE_RESULTS / case_id / "inputs", directory / "inputs")
        preparation = json.loads(archived_read(directory, "preparation.json"))
        ligand = archived_read(directory, "ligand.pdbqt").decode()
        reference = archived_read(directory, "reference.json")
        rows = []
        for receptor in preparation["receptors"]:
            for seed_index, seed in enumerate(config["seeds"]):
                profile = receptor["id"]
                score = -6.0 if case["species_id"] == "D_galactose" else -5.0
                score += (0.2 if case["anomer"] == "beta" else 0) + 0.01 * seed_index
                pose = f"{profile}_{seed}.pdbqt"
                log = f"{profile}_{seed}.log"
                pose_text = f"MODEL 1\nREMARK VINA RESULT: {score:.3f} 0.000 0.000\n{ligand.rstrip()}\nENDMDL\n"
                (directory / pose).write_text(pose_text)
                prepared_case = root / "data/processed/pol_specificity" / case["prepared_dir"]
                command = [str(root / "tools/vina-1.2.7"), "--receptor", str(prepared_case / receptor["pdbqt_file"]),
                           "--ligand", str(prepared_case / "ligand.pdbqt"), "--out", str(directory / pose),
                           "--scoring", "vina", "--cpu", str(config["cpu"]), "--seed", str(seed),
                           "--exhaustiveness", str(config["exhaustiveness"]), "--num_modes", str(config["num_modes"]),
                           "--min_rmsd", str(config["min_rmsd_A"]), "--energy_range", str(config["energy_range_kcal_mol"])]
                for axis, center, size in zip("xyz", preparation["box_center"], preparation["box_size"]):
                    command.extend(["--center_" + axis, str(center), "--size_" + axis, str(size)])
                center = preparation["box_center"]
                size = preparation["box_size"]
                log_text = ("AutoDock Vina v1.2.7\nScoring function : vina\n"
                            f"Rigid receptor: {prepared_case / receptor['pdbqt_file']}\n"
                            f"Ligand: {prepared_case / 'ligand.pdbqt'}\n"
                            f"Grid center: X {center[0]:.6g} Y {center[1]:.6g} Z {center[2]:.6g}\n"
                            f"Grid size  : X {size[0]:.6g} Y {size[1]:.6g} Z {size[2]:.6g}\n"
                            "Exhaustiveness: 16\nCPU: 2\n"
                            f"Performing docking (random seed: {seed}) ...\n"
                            "mode |   affinity | dist from best mode\n"
                            "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
                            "-----+------------+----------+----------\n"
                            f"    1      {score:.3f}      0.000      0.000\n")
                (directory / log).write_text(log_text)
                fingerprint = {"executable_sha256": executable_sha, "config_sha256": frozen_sha,
                               "ligand_sha256": digest_bytes(ligand.encode()),
                               "preparation_sha256": digest_bytes(archived_read(directory, "preparation.json")),
                               "reference_sha256": digest_bytes(reference), "full_assembly_sha256": None,
                               "receptor_sha256": digest_bytes(archived_read(directory, receptor["pdbqt_file"])),
                               "receptor_pdb_sha256": digest_bytes(archived_read(directory, receptor["pdb_file"]))}
                rows.append({"case_id": case_id, "family_id": "POL_CBM67", "receptor": profile,
                             "receptor_state": case["state"], "assembly_condition": receptor["assembly_condition"],
                             "seed": seed, "information_regime": "true-site-local", "status": "ok", "exit_code": 0,
                             "fingerprint": fingerprint, "command": command, "pose_file": pose, "log_file": log,
                             "returned_poses": 1, "pose_sha256": sha256(directory / pose),
                             "vina_version": "AutoDock Vina v1.2.7", "wall_seconds": 1.0, "cpu_seconds": 2.0,
                             "software": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "gemmi", "rdkit", "meeko")}})
        dump(directory / "runs.json", rows)
        shutil.copyfile(root / "benchmarks/pol_experiment.json", directory / "configuration.json")
    dump(result / "execution_status.json", {"status": "completed", "expected_runs": 96,
         "completed_cases": 16, "terminal_runs": 96, "successful_runs": 96, "failed_runs": 0,
         "configuration_sha256": frozen_sha})
    (root / "SYNTHETIC_FIXTURE.txt").write_text("Synthetic test outputs only; no scientific docking was performed.\n")


class POLReportControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_temporary = tempfile.TemporaryDirectory()
        cls.base = Path(cls.base_temporary.name) / "base"
        synthetic_fixture(cls.base)

    @classmethod
    def tearDownClass(cls):
        cls.base_temporary.cleanup()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "repo"
        shutil.copytree(self.base, self.root)
        self.results = self.root / "results/pol_specificity"
        self.case = planned_cases({})[0]["case_id"]
        self.case_dir = self.results / self.case

    def report(self):
        return generate_pol_report(self.root)

    def rows(self, case=None):
        return json.loads((self.results / (case or self.case) / "runs.json").read_text())

    def save_rows(self, rows, case=None):
        dump(self.results / (case or self.case) / "runs.json", rows)

    def archive_edit(self, filename, content, case=None):
        """Update local hash chains to exercise semantic, not merely byte, gates."""
        case = case or self.case
        directory = self.results / case
        if isinstance(content, str):
            content = content.encode()
        elif not isinstance(content, bytes):
            content = (json.dumps(content, indent=2, sort_keys=True) + "\n").encode()
        manifest_path = directory / "inputs/manifest.json"
        manifest = json.loads(manifest_path.read_text())
        entry = next(r for r in manifest["files"] if r["filename"] == filename)
        compressed = directory / "inputs" / entry["archive"]
        compressed.write_bytes(gzip.compress(content, mtime=0))
        entry.update(sha256=digest_bytes(content), archive_sha256=sha256(compressed), original_bytes=len(content))
        dump(manifest_path, manifest)
        panel_path = self.results / "panel_preparation.json"
        panel = json.loads(panel_path.read_text())
        item = next(r for r in panel["cases"] if r["case_id"] == case)
        item["prepared_artifact_sha256"][filename] = digest_bytes(content)
        if filename == "preparation.json":
            item["preparation_sha256"] = digest_bytes(content)
        dump(panel_path, panel)
        fingerprint_key = {"preparation.json": "preparation_sha256", "reference.json": "reference_sha256",
                           "ligand.pdbqt": "ligand_sha256"}.get(filename)
        rows = self.rows(case)
        preparation = json.loads(archived_read(directory, "preparation.json"))
        receptors = {r["id"]: r for r in preparation["receptors"]}
        for row in rows:
            if fingerprint_key:
                row["fingerprint"][fingerprint_key] = digest_bytes(content)
            receptor = receptors[row["receptor"]]
            if filename == receptor["pdb_file"]:
                row["fingerprint"]["receptor_pdb_sha256"] = digest_bytes(content)
            if filename == receptor["pdbqt_file"]:
                row["fingerprint"]["receptor_sha256"] = digest_bytes(content)
        self.save_rows(rows, case)

    def test_complete_synthetic_evidence_replays_all_96_and_eight_separate_conditions(self):
        summary, markdown = self.report()
        self.assertTrue(summary["metrics_eligible"])
        self.assertEqual(summary["denominator"]["successful_runs"], 96)
        self.assertEqual(len(summary["runs"]), 96)
        aggregation = summary["aggregation"]
        self.assertEqual(len(aggregation["conditions"]), 8)
        self.assertEqual(aggregation["experimental_species"], 2)
        self.assertEqual(aggregation["experimental_epimer_contrasts"], 1)
        self.assertEqual(aggregation["independent_families"], 1)
        self.assertAlmostEqual(aggregation["primary"]["glucose_minus_galactose_kcal_mol"], 1.0)
        self.assertIsNone(aggregation["auroc"])
        self.assertIsNone(aggregation["auprc"])
        self.assertFalse(aggregation["affinity_calibration"])
        self.assertTrue(summary["source_hash_manifest"])
        self.assertIn("anomer", markdown.lower())
        self.assertIn("one", markdown.lower())

    def test_portable_report_needs_no_original_processed_input_directory(self):
        self.assertFalse((self.root / "data/processed/pol_specificity").exists())
        summary, _ = self.report()
        self.assertTrue(summary["metrics_eligible"])

    def test_missing_run_keeps_96_denominator_and_disables_all_aggregates(self):
        self.save_rows(self.rows()[:-1])
        summary, _ = self.report()
        self.assertFalse(summary["metrics_eligible"])
        self.assertEqual(summary["denominator"]["successful_runs"], 95)
        self.assertEqual(len(summary["runs"]), 96)
        self.assertIsNone(summary["aggregation"])

    def test_explicit_failed_run_is_visible_and_never_averaged_away(self):
        rows = self.rows()
        rows[0].update(status="failed", exit_code=1, error="synthetic intentional failure")
        self.save_rows(rows)
        summary, _ = self.report()
        self.assertFalse(summary["metrics_eligible"])
        self.assertEqual(summary["denominator"]["successful_runs"], 95)
        self.assertIsNone(summary["aggregation"])

    def test_duplicate_unknown_seed_and_wrong_case_run_are_rejected(self):
        original = self.rows()
        variants = [original + [original[0]],
                    [{**original[0], "seed": 9999}, *original[1:]],
                    [{**original[0], "case_id": planned_cases({})[1]["case_id"]}, *original[1:]]]
        for rows in variants:
            self.save_rows(rows)
            with self.subTest(row=rows[0]), self.assertRaises(ValueError):
                self.report()
        self.save_rows(original)

    def test_unexpected_run_bearing_case_directory_is_rejected(self):
        extra = self.results / "unexpected_case"
        extra.mkdir()
        dump(extra / "runs.json", self.rows())
        with self.assertRaises(ValueError):
            self.report()

    def test_source_freeze_and_execution_snapshot_mutations_are_rejected(self):
        paths = [self.root / "data/archive/independent_specificity_sources/GLA.cif",
                 self.root / "data/archive/independent_specificity_sources/CA.cif",
                 self.results / "freeze_receipt.json",
                 self.results / "execution_snapshot/aav_benchmark/pol_specificity.py"]
        for path in paths:
            original = path.read_bytes()
            if path.suffix == ".json":
                changed = json.loads(original)
                changed["configuration_sha256"] = "0" * 64
                dump(path, changed)
            else:
                path.write_bytes(original + b"\n# changed synthetic control\n")
            with self.subTest(path=path.name), self.assertRaises(ValueError):
                self.report()
            path.write_bytes(original)

    def test_prepared_archive_and_raw_pose_byte_hashes_are_enforced(self):
        archive = self.case_dir / "inputs/ligand.pdbqt.gz"
        original = archive.read_bytes()
        archive.write_bytes(original + b"changed")
        with self.assertRaises(ValueError):
            self.report()
        archive.write_bytes(original)
        row = self.rows()[0]
        pose = self.case_dir / row["pose_file"]
        pose.write_text(pose.read_text() + "REMARK changed\n")
        with self.assertRaises(ValueError):
            self.report()

    def test_actual_command_cpu_box_and_duplicate_option_are_checked(self):
        original = self.rows()
        for flag, value in [("--cpu", "1"), ("--center_x", "999")]:
            rows = json.loads(json.dumps(original))
            command = rows[0]["command"]
            command[command.index(flag) + 1] = value
            self.save_rows(rows)
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                self.report()
        rows = json.loads(json.dumps(original))
        rows[0]["command"] += ["--cpu", "2"]
        self.save_rows(rows)
        with self.assertRaises(ValueError):
            self.report()

    def test_raw_log_table_must_match_actual_pose_ranks_and_scores(self):
        row = self.rows()[0]
        log = self.case_dir / row["log_file"]
        content = log.read_text()
        self.assertIn("    1      -6.000", content)
        log.write_text(content.replace("    1      -6.000", "    1      -9.000"))
        with self.assertRaises(ValueError):
            self.report()

    def test_logged_mode_outside_saved_energy_window_is_explicitly_retained_as_log_evidence(self):
        row = self.rows()[0]
        path = self.case_dir / row["log_file"]
        path.write_text(path.read_text() + "    2      20.00       1.000      1.000\n")
        summary, _ = self.report()
        self.assertTrue(summary["metrics_eligible"])
        self.assertEqual(summary["resources"]["energy_window_omitted_modes"], 1)
        verified = next(r for r in summary["runs"] if r["case_id"] == self.case and r["metal_profile"] == row["receptor"] and r["seed"] == row["seed"])
        self.assertEqual(verified["logged_modes"], 2)
        self.assertEqual(verified["returned_poses"], 1)
        self.assertEqual(verified["omitted_log_scores_kcal_mol"], [20.0])

    def test_logged_mode_inside_energy_window_cannot_be_silently_absent_from_poses(self):
        row = self.rows()[0]
        path = self.case_dir / row["log_file"]
        path.write_text(path.read_text() + "    2      -5.000       1.000      1.000\n")
        with self.assertRaisesRegex(ValueError, "inside the frozen output energy window"):
            self.report()

    def test_integer_formatted_log_score_matches_the_same_three_decimal_pose_score(self):
        row = self.rows()[0]
        path = self.case_dir / row["log_file"]
        path.write_text(path.read_text().replace("    1      -6.000", "    1      -6"))
        self.assertTrue(self.report()[0]["metrics_eligible"])

    def test_raw_log_protocol_fields_require_exact_values_not_substring_prefixes(self):
        row = self.rows()[0]
        log = self.case_dir / row["log_file"]
        original = log.read_text()
        for before, after in [("CPU: 2", "CPU: 20"),
                              ("Exhaustiveness: 16", "Exhaustiveness: 160"),
                              ("Scoring function : vina", "Scoring function : vinardo"),
                              ("AutoDock Vina v1.2.7", "AutoDock Vina v1.2.70")]:
            log.write_text(original.replace(before, after))
            with self.subTest(field=before), self.assertRaises(ValueError):
                self.report()
        log.write_text(original)

    def test_actual_executable_command_identity_is_checked(self):
        rows = self.rows()
        rows[0]["command"][0] = str(self.root / "tools/unapproved_program")
        self.save_rows(rows)
        with self.assertRaises(ValueError):
            self.report()

    def test_raw_pose_element_map_and_3d_stereochemistry_are_checked_even_with_rehashed_pose(self):
        rows = self.rows()
        pose = self.case_dir / rows[0]["pose_file"]
        original = pose.read_text()
        lines = original.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM")) and line[77:].strip() == "C")
        lines[index] = lines[index][:77] + "N"
        pose.write_text("\n".join(lines) + "\n")
        rows[0]["pose_sha256"] = sha256(pose)
        self.save_rows(rows)
        with self.assertRaises(ValueError):
            self.report()
        # Reflect all heavy atoms around their centroid: map/elements are intact,
        # but each tetrahedral centre has the opposite handedness.
        lines = original.splitlines()
        heavy = [i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM"))
                 and line[77:].strip() not in {"H", "HD", "HS"}]
        center_x = np.mean([float(lines[i][30:38]) for i in heavy])
        for i in heavy:
            x = 2 * center_x - float(lines[i][30:38])
            lines[i] = lines[i][:30] + f"{x:8.3f}" + lines[i][38:]
        pose.write_text("\n".join(lines) + "\n")
        rows[0]["pose_sha256"] = sha256(pose)
        self.save_rows(rows)
        with self.assertRaises(ValueError):
            self.report()

    def test_duplicate_heavy_atom_map_and_anomer_source_alias_are_rejected_after_rehash(self):
        reference = json.loads(archived_read(self.case_dir, "reference.json"))
        original = archived_read(self.case_dir, "reference.json")
        reference["names"][1] = reference["names"][0]
        self.archive_edit("reference.json", reference)
        with self.assertRaises(ValueError):
            self.report()
        self.archive_edit("reference.json", original)
        source = json.loads(archived_read(self.case_dir, "source_graph.json"))
        beta_case = self.results / "POL_GAL_beta__holo__N"
        beta_source = json.loads(archived_read(beta_case, "source_graph.json"))
        source["ccd_graph"] = beta_source["ccd_graph"]
        self.archive_edit("source_graph.json", source)
        with self.assertRaises(ValueError):
            self.report()

    def test_actual_sdf_graph_cannot_be_replaced_by_chemically_matched_epimer_after_rehash(self):
        other = self.results / "POL_GLC_alpha__holo__N"
        self.archive_edit("ligand.sdf", archived_read(other, "ligand.sdf"))
        with self.assertRaises(ValueError):
            self.report()

    def test_receptor_native_calcium_count_is_checked_after_rehash(self):
        filename = "holo_Ca_retained.pdb"
        lines = archived_read(self.case_dir, filename).decode().splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CA2")
        del lines[index]
        self.archive_edit(filename, "\n".join(lines) + "\n")
        with self.assertRaises(ValueError):
            self.report()

    def test_write_and_check_detect_stale_derived_report(self):
        write_pol_report(self.root)
        write_pol_report(self.root, check=True)
        summary_path = self.results / "summary.json"
        summary = json.loads(summary_path.read_text())
        summary["denominator"]["successful_runs"] = 95
        dump(summary_path, summary)
        with self.assertRaises(ValueError):
            write_pol_report(self.root, check=True)


if __name__ == "__main__":
    unittest.main()
