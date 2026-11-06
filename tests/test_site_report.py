"""Fixed-denominator and evidence-tampering checks for site report generation."""
from copy import deepcopy
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.site_report import FIELDS, METHODS, generate_site_report, write_site_report


REPOSITORY = Path(__file__).resolve().parents[1]


class SiteReportEvidence(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "benchmarks").mkdir()
        for name in ("site_experiment.json", "site_reference_ledger.json", "site_methods.lock.json",
                     "p2rank.lock.json", "site_geometry_amendment.json"):
            shutil.copyfile(REPOSITORY / "benchmarks" / name, self.root / "benchmarks" / name)
        self.directory = self.root / "results/site_discovery"
        self.directory.mkdir(parents=True)
        shutil.copyfile(REPOSITORY / "results/site_discovery/freeze_receipt.json", self.directory / "freeze_receipt.json")
        config = json.loads((self.root / "benchmarks/site_experiment.json").read_text())
        sources, rows = {}, []

        def bind(path, value):
            write_json(path, value)
            sources[str(path.relative_to(self.root))] = sha256(path)

        for name in ("site_experiment.json", "site_reference_ledger.json", "site_geometry_amendment.json"):
            path = self.root / "benchmarks" / name
            sources[str(path.relative_to(self.root))] = sha256(path)
        # Synthetic evaluated outcomes exercise fixed denominators and evidence
        # binding. These files are deliberately not scientific predictions.
        for case in config["cases"]:
            condition = case["condition_id"]
            details = {"condition_id": condition, "complex_id": case["complex_id"], "methods": {}}
            write_json(self.directory / condition / "fpocket/execution.json", {"status": "failed", "returncode": 1})
            for method in METHODS:
                baseline = method == "random_surface"
                failed = method == "fpocket" and case["complex_id"] == "aav9"
                status = "failed" if failed else "ok"
                trials = []
                for index in range(100 if baseline else 1):
                    # Both real methods succeed on 154L. The null has one
                    # successful draw in each condition, never 100 observations.
                    hit = int(index == 0) if baseline else int(condition == "154L")
                    trials.append({**{field: float(hit) for field in FIELDS}, "raw_candidate_count": 0 if failed else 5,
                                   "returned_distinct_candidates": 0 if failed else 5,
                                   "top1_distance_A": 1. if hit else None, "top5_min_distance_A": 1. if hit else None,
                                   "first_matching_class_4A": 1 if hit else None})
                row = {"condition_id": condition, "complex_id": case["complex_id"], "state": case["state"],
                       "method": method, "status": status, "eligible_conditions": 1,
                       "baseline_replicates": 100 if baseline else None,
                       **{field: (.01 if baseline else float(condition == "154L")) for field in FIELDS},
                       "raw_candidate_count": 0 if failed else 5, "distinct_candidate_count": 0 if failed else 5,
                       "top1_distance_A": None if baseline else trials[0]["top1_distance_A"],
                       "top5_min_distance_A": None if baseline else trials[0]["top5_min_distance_A"],
                       "first_matching_class_4A": None if baseline else trials[0]["first_matching_class_4A"]}
                relative = {"fpocket": "fpocket_v2/execution.json", "p2rank": "p2rank/result.json",
                            "random_surface": "random_surface_v2/result.json"}[method]
                record = {"status": ("failed" if failed else "completed") if method == "fpocket" else status,
                          "returncode": -5 if failed else 0}
                path = self.directory / condition / relative
                bind(path, record)
                evidence = {"status": status, "source_record_sha256": sha256(path), "evaluations": trials}
                row["merge_sensitivity"] = {distance: deepcopy(row) for distance in ("2.0", "6.0")}
                evidence["merge_sensitivity"] = {distance: {"evaluations": deepcopy(trials)} for distance in ("2.0", "6.0")}
                details["methods"][method] = evidence
                rows.append(row)
            bind(self.directory / condition / "evaluation_details.json", details)
        self.evaluation = {"complete": True, "pending": [], "eligible_complexes": 8,
                           "eligible_conditions": 9, "eligible_method_conditions": 27,
                           "rows": rows, "source_sha256": sources}
        self.save_evaluation()

    def save_evaluation(self):
        write_json(self.directory / "evaluation.json", self.evaluation)

    def mutate_details(self, change):
        path = self.directory / "154L/evaluation_details.json"
        details = json.loads(path.read_text())
        change(details)
        write_json(path, details)
        self.evaluation["source_sha256"][str(path.relative_to(self.root))] = sha256(path)
        self.save_evaluation()

    def test_fixed_denominator_and_monte_carlo_aggregation(self):
        summary, markdown = generate_site_report(self.root)
        self.assertEqual(summary["primary_holo_systems"], 8)
        self.assertEqual(summary["method_conditions"], 27)
        self.assertEqual(summary["primary_macro"]["fpocket"]["top5_success_4A"], 1 / 8)
        self.assertEqual(summary["primary_macro"]["random_surface"]["top5_success_4A"], .01)
        self.assertEqual(summary["fpocket_corrected_failures"], 2)
        self.assertEqual(summary["fpocket_original_reader_failures"], 9)
        self.assertTrue(summary["merge_sensitivity_complete"])
        self.assertIn("conditional Monte Carlo null replicates", markdown)
        self.assertIn("without an additional independent system", markdown)

    def test_duplicate_or_missing_condition_is_rejected(self):
        self.evaluation["rows"][-1] = deepcopy(self.evaluation["rows"][0])
        self.save_evaluation()
        with self.assertRaisesRegex(ValueError, "Duplicate, missing or extra"):
            generate_site_report(self.root)

    def test_failed_execution_cannot_leave_primary_denominator(self):
        self.evaluation["rows"] = [row for row in self.evaluation["rows"] if row["status"] != "failed"]
        self.save_evaluation()
        with self.assertRaisesRegex(ValueError, "Duplicate, missing or extra"):
            generate_site_report(self.root)

    def test_source_tampering_is_rejected(self):
        path = self.directory / "154L/p2rank/result.json"
        path.write_text("{}\n")
        with self.assertRaisesRegex(ValueError, "source changed"):
            generate_site_report(self.root)

    def test_row_aggregate_tampering_is_rejected(self):
        row = next(row for row in self.evaluation["rows"] if row["condition_id"] == "154L" and row["method"] == "p2rank")
        row["top5_success_4A"] = 0.
        self.save_evaluation()
        with self.assertRaisesRegex(ValueError, "disagrees with detailed evidence"):
            generate_site_report(self.root)

    def test_baseline_replication_cannot_inflate_or_drop_draws(self):
        self.mutate_details(lambda details: details["methods"]["random_surface"]["evaluations"].pop())
        with self.assertRaisesRegex(ValueError, "replicate denominator"):
            generate_site_report(self.root)

    def test_unbound_prediction_reference_is_rejected(self):
        self.mutate_details(lambda details: details["methods"]["p2rank"].update(source_record_sha256="0" * 64))
        with self.assertRaisesRegex(ValueError, "Unbound normalized"):
            generate_site_report(self.root)

    def test_missing_merge_sensitivity_is_explicit(self):
        self.evaluation["rows"][0].pop("merge_sensitivity")
        self.save_evaluation()
        summary, markdown = generate_site_report(self.root)
        self.assertFalse(summary["merge_sensitivity_complete"])
        self.assertIn("Those evaluations are pending", markdown)

    def test_check_mode_detects_stale_generated_document(self):
        output = write_site_report(self.root)
        self.assertEqual(write_site_report(self.root, check=True), output)
        output.write_text(output.read_text() + "stale\n")
        with self.assertRaisesRegex(ValueError, "report is stale"):
            write_site_report(self.root, check=True)


if __name__ == "__main__":
    unittest.main()
