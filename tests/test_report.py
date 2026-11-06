"""Controls for denominator completeness and reproducible report provenance."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from aav_benchmark.report import generate_report, write_report
from aav_benchmark.statistics import quantitative_specificity_by_context


SEEDS = [104729, 130363, 155921]
CASES = ["154L", "1GNY", "1KJL", "1QFO", "2J1T", "1PWB", "2VUZ"]


def put_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def make_runs(root, case, directory_name, conditions, *, energy=-5, structural=True):
    directory = root / "results/development" / directory_name
    directory.mkdir(parents=True, exist_ok=True)
    config = {"seeds": SEEDS, "information_regime": "true-site-local", "cpu": 2}
    put_json(directory / "configuration.json", config)
    config_hash = hashlib.sha256((directory / "configuration.json").read_bytes()).hexdigest()
    runs, summaries, metrics = [], [], []
    for condition in conditions:
        for seed in SEEDS:
            pose_file = f"{condition}_{seed}.pdbqt"
            content = f"MODEL 1\nREMARK VINA RESULT: {energy} 0 0\nENDMDL\n"
            (directory / pose_file).write_text(content)
            raw = dict(case_id=case, receptor=condition, seed=seed, status="ok",
                       fingerprint={"config_sha256": config_hash}, pose_file=pose_file,
                       pose_sha256=hashlib.sha256(content.encode()).hexdigest(), returned_poses=1,
                       wall_seconds=1.25, cpu_seconds=2.0)
            runs.append(raw)
            if structural:
                summaries.append({**raw, "top1_rmsd_A": 1.0, "top5_min_rmsd_A": 1.0,
                                  "best_returned_rmsd_A": 1.0, "top1_score": energy,
                                  "top1_contact_recall": 0.5, "top1_contact_precision": 0.25,
                                  **{f"{prefix}_success_{threshold}A": True
                                     for threshold in (2, 3) for prefix in ("top1", "top5", "retained")}})
                metrics.append(dict(case_id=case, receptor=condition, seed=seed, rank=1,
                                    rmsd_A=1.0, vina_score=energy, contact_recall=0.5,
                                    contact_precision=0.25, interface_residue_recall=0.5,
                                    interface_residue_precision=0.5, native_contact_count=4,
                                    native_interface_residue_count=2, evaluation_pairs_below_1_5A=0,
                                    evaluation_receptor_scope="provided_receptor_only"))
    put_json(directory / "runs.json", runs)
    if structural:
        put_json(directory / "summary.json", summaries)
        with (directory / "pose_metrics.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
            writer.writeheader()
            writer.writerows(metrics)
    return directory


def fixture(root):
    cases = [dict(case_id=case, status=("ok" if case == "154L" else "failed" if case == "1QFO" else "unsupported"),
                  eligible_preparation_denominator=True) for case in CASES]
    put_json(root / "benchmarks/development_preparation.json",
             dict(eligible_cases=7, successful_preparations=1, failed_preparations=1,
                  unsupported_preparations=5, cases=cases))
    put_json(root / "benchmarks/external_run.json", {"cases": CASES})
    put_json(root / "benchmarks/development.json", {"seeds": SEEDS})
    make_runs(root, "AAV9", "aav9", ["holo_complete", "holo_monomer", "apo_complete", "apo_monomer"])
    make_runs(root, "154L", "154l", ["holo_complete"])
    source = {"panel_id": "panel", "target": {"uniprot_id": "P17931"},
              "rows": [dict(ligand_id=ligand, kd_uM=kd, kd_sem_uM=1.0, expected_heavy_atoms=atoms)
                       for ligand, kd, atoms in (("LacNAc", 24.0, 26), ("lactose", 93.5, 23), ("galactose", 3100.0, 12))]}
    put_json(root / "benchmarks/galectin3_itc_panel.json", source)
    entries, correlations = [], []
    for anomer in ("alpha", "beta"):
        rows = []
        for item in source["rows"]:
            case = f"gal3_{item['ligand_id']}_{anomer}"
            energy = -float(item["expected_heavy_atoms"])
            make_runs(root, case, case, ["holo_complete"], energy=energy, structural=False)
            entries.append(dict(case_id=case, ligand_id=item["ligand_id"], anomer=anomer,
                                top1_scores_kcal_mol=[energy] * 3, seed_mean_top1_score=energy,
                                kd_uM=item["kd_uM"], kd_sem_uM=1.0, heavy_atoms=item["expected_heavy_atoms"],
                                eligible_runs=3, completed_runs=3))
            rows.append(dict(context_id="panel", target_id="P17931", ligand_id=item["ligand_id"],
                             kd=item["kd_uM"], kd_unit="uM", score=-energy))
        correlation = quantitative_specificity_by_context(rows)
        correlations.append(dict(anomer=anomer, vina=correlation, heavy_atom_baseline=correlation))
    put_json(root / "benchmarks/gal3_run.json", {"cases": [item["case_id"] for item in entries]})
    put_json(root / "results/development/gal3_panel_summary.json",
             dict(ligands=entries, correlations=correlations))


def native_fixture(root):
    parent = root / "results/development/aav9"
    directory = parent / "native_control"
    directory.mkdir()
    ligand = b"fixed native ligand\n"
    (directory / "native_ligand.pdbqt").write_bytes(ligand)
    (directory / "native_ligand.sdf").write_text("source SDF\n")
    ligand_hash = hashlib.sha256(ligand).hexdigest()
    put_json(directory / "ligand_preparation.json",
             dict(native_ligand_sha256=ligand_hash, diagnostic_only=True, count_as_docking_success=False,
                  native_heavy_atoms_optimized=False, hydrogen_optimization=False, native_heavy_rmsd_A=0.0))
    runs = json.loads((parent / "runs.json").read_text())
    summaries = json.loads((parent / "summary.json").read_text())
    for run in runs + summaries:
        run["fingerprint"].update({name: "0" * 64 for name in
                                  ("executable_sha256", "receptor_sha256", "preparation_sha256", "reference_sha256")})
    put_json(parent / "runs.json", runs)
    put_json(parent / "summary.json", summaries)
    controls = []
    for condition in sorted({run["receptor"] for run in runs}):
        selected = [run for run in summaries if run["receptor"] == condition]
        log_name = condition + "_score_only.log"
        log = "Estimated Free Energy of Binding : -1.000 (kcal/mol)\n"
        (directory / log_name).write_text(log)
        controls.append(dict(receptor=condition, status="ok", diagnostic_only=True, count_as_docking_success=False,
                             fingerprint={**selected[0]["fingerprint"], "native_ligand_sha256": ligand_hash},
                             native_heavy_rmsd_A=0.0, native_score_only_kcal_mol=-1.0,
                             log_file=log_name, log_sha256=hashlib.sha256(log.encode()).hexdigest(), wall_seconds=0.1,
                             docked_top1_comparisons=[dict(seed=run["seed"], pose_file=run["pose_file"],
                                                          pose_sha256=run["pose_sha256"], docked_top1_rmsd_A=1.0,
                                                          docked_top1_score_kcal_mol=-5.0,
                                                          native_minus_docked_top1_score_kcal_mol=4.0) for run in selected]))
    put_json(directory / "summary.json", controls)


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        fixture(self.root)

    def mutate_json(self, path, callback):
        target = self.root / path
        value = json.loads(target.read_text())
        callback(value)
        put_json(target, value)

    def test_deterministic_report_retains_scope_and_source_hashes(self):
        report = generate_report(self.root)
        self.assertEqual(report, generate_report(self.root))
        self.assertIn("original scientific purpose remains incomplete", report)
        self.assertIn("**7 eligible cases**", report)
        self.assertIn("**33/33 completed runs**", report)
        self.assertIn("no tested nonbinders", report)
        self.assertIn("Size alone ranks this panel perfectly", report)
        source = self.root / "results/development/154l/summary.json"
        self.assertIn(hashlib.sha256(source.read_bytes()).hexdigest(), report)

    def test_missing_configured_seed_cannot_shrink_denominator(self):
        self.mutate_json("results/development/154l/runs.json", lambda value: value.pop())
        with self.assertRaisesRegex(ValueError, "three-seed denominator"):
            generate_report(self.root)

    def test_incomplete_existing_summary_fails(self):
        self.mutate_json("results/development/154l/summary.json", lambda value: value.pop())
        with self.assertRaisesRegex(ValueError, "three-seed structural summary"):
            generate_report(self.root)

    def test_absent_listed_completed_case_fails(self):
        self.mutate_json("benchmarks/development_preparation.json", lambda value: value["cases"][0].update(docking_complete=True))
        (self.root / "results/development/154l/summary.json").unlink()
        with self.assertRaisesRegex(ValueError, "completed docking case has no summary"):
            generate_report(self.root)

    def test_pending_successful_preparation_is_explicit(self):
        def change(value):
            value["cases"][1]["status"] = "ok"
            value["successful_preparations"] += 1
            value["unsupported_preparations"] -= 1
        self.mutate_json("benchmarks/development_preparation.json", change)
        report = generate_report(self.root)
        self.assertIn("| 1GNY | ok | not completed / no summary |", report)
        self.assertIn("**7 eligible cases**", report)

    def test_pose_change_fails_checksum_validation(self):
        path = self.root / "results/development/154l/holo_complete_104729.pdbqt"
        path.write_text(path.read_text() + "tampered\n")
        with self.assertRaisesRegex(ValueError, "Changed pose artifact"):
            generate_report(self.root)

    def test_contact_summary_must_match_per_pose_evidence(self):
        self.mutate_json("results/development/154l/summary.json", lambda value: value[0].update(top1_contact_recall=0.9))
        with self.assertRaisesRegex(ValueError, "top1_contact_recall"):
            generate_report(self.root)

    def test_native_contact_denominator_cannot_change(self):
        path = self.root / "results/development/aav9/pose_metrics.csv"
        rows = list(csv.DictReader(path.read_text().splitlines()))
        rows[0]["native_contact_count"] = "8"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(ValueError, "Native contact denominator changed"):
            generate_report(self.root)

    def test_panel_requires_all_configured_seeds(self):
        self.mutate_json("results/development/gal3_LacNAc_alpha/runs.json", lambda value: value.pop())
        with self.assertRaisesRegex(ValueError, "three-seed denominator"):
            generate_report(self.root)

    def test_native_control_stays_separate_from_docking_success(self):
        native_fixture(self.root)
        report = generate_report(self.root)
        self.assertIn("**never docking successes**", report)
        self.assertIn("**33/33 completed runs**", report)
        self.assertIn("| holo_complete | 0.000000 | -1.000 | -5.000 | 4.000 |", report)
        self.mutate_json("results/development/aav9/native_control/summary.json",
                         lambda value: value[0].update(count_as_docking_success=True))
        with self.assertRaisesRegex(ValueError, "cannot be counted as successful docking"):
            generate_report(self.root)

    def test_native_control_log_change_is_detected(self):
        native_fixture(self.root)
        (self.root / "results/development/aav9/native_control/holo_complete_score_only.log").write_text("changed")
        with self.assertRaisesRegex(ValueError, "log checksum mismatch"):
            generate_report(self.root)

    def test_check_mode_detects_stale_output_without_rewriting(self):
        output = write_report(self.root)
        self.assertEqual(write_report(self.root, check=True), output)
        output.write_text("stale\n")
        with self.assertRaisesRegex(ValueError, "stale or missing"):
            write_report(self.root, check=True)
        self.assertEqual(output.read_text(), "stale\n")


if __name__ == "__main__":
    unittest.main()
