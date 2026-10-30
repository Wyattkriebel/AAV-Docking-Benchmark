"""Statistical controls against pseudoreplication, attrition and assay pooling."""
import math
import unittest

from aav_benchmark.statistics import (
    aggregate_runs, binary_ranking_metrics, paired_family_bootstrap,
    specificity_by_context, quantitative_specificity_by_context,
)


def run(family, complex_id, method, seed, success=None, status="ok", near_native=None):
    row = dict(family_id=family, complex_id=complex_id, method=method,
               seed=seed, status=status, success=success)
    if near_native is not None:
        row["near_native_present"] = near_native
    return row


class AggregationTests(unittest.TestCase):
    def test_seeds_then_complexes_then_families(self):
        # Family A has three systems and extra seeds; neither adds biological weight.
        rows = [run("A", "a1", "m", seed, 1) for seed in range(10)]
        rows += [run("A", "a2", "m", 0, 1), run("A", "a3", "m", 0, 1)]
        rows += [run("B", "b1", "m", 0, 0)]
        report = aggregate_runs(rows)["overall"][0]
        self.assertEqual(report["family_macro_recovery"], 0.5)
        self.assertEqual(report["complex_macro_recovery"], 0.75)
        self.assertEqual(report["families"], 2)

    def test_failed_and_unsupported_in_eligible_denominator(self):
        rows = [run("A", "a", "m", 1, True, near_native=True),
                run("A", "a", "m", 2, status="failed"),
                run("A", "a", "m", 3, status="unsupported")]
        report = aggregate_runs(rows, expected_seeds=[1, 2, 3])["overall"][0]
        self.assertAlmostEqual(report["family_macro_recovery"], 1 / 3)
        self.assertEqual(report["conditional_completed_recovery"], 1)
        self.assertEqual(report["conditional_near_native_recovery"], 1)
        self.assertEqual(report["eligible_runs"], 3)
        self.assertEqual(report["failed_runs"], 1)
        self.assertEqual(report["unsupported_runs"], 1)

    def test_uneven_seeds_do_not_change_complex_weight(self):
        rows = [run("A", "a1", "m", 1, 1), run("A", "a1", "m", 2, 0),
                run("A", "a2", "m", 1, 1)]
        report = aggregate_runs(rows)["overall"][0]
        self.assertEqual(report["family_macro_recovery"], 0.75)

    def test_unknown_near_native_status_is_not_zero(self):
        rows = [run("A", "a", "m", 1, 1)]
        report = aggregate_runs(rows)["overall"][0]
        self.assertIsNone(report["near_native_runs"])
        self.assertIsNone(report["conditional_near_native_recovery"])

    def test_missing_duplicate_and_invalid_runs_fail(self):
        valid = run("A", "a", "m", 1, 1)
        with self.assertRaises(ValueError):
            aggregate_runs([valid], expected_seeds=[1, 2])
        with self.assertRaises(ValueError):
            aggregate_runs([valid, valid])
        for bad in (math.nan, 2, "1", None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                aggregate_runs([dict(valid, success=bad)])
        with self.assertRaises(ValueError):
            aggregate_runs([dict(valid, status="failed")])
        with self.assertRaises(ValueError):
            aggregate_runs([dict(valid, near_native_present=False)])
        with self.assertRaises(ValueError):
            aggregate_runs([valid, dict(valid, method="other", family_id="B")])


class PairedBootstrapTests(unittest.TestCase):
    def test_exact_pairing_is_required(self):
        with self.assertRaises(ValueError):
            paired_family_bootstrap([run("A", "a", "x", 1, 1), run("A", "a", "y", 2, 0)], "x", "y")

    def test_whole_family_pairing_and_reproducibility(self):
        rows = []
        for family, wins in (("A", [1, 1, 1]), ("B", [0]), ("C", [1])):
            for index, win in enumerate(wins):
                for seed in (1, 2, 3):
                    rows += [run(family, family + str(index), "x", seed, win),
                             run(family, family + str(index), "y", seed, 0)]
        first = paired_family_bootstrap(rows, "x", "y", expected_seeds=[1, 2, 3])
        second = paired_family_bootstrap(rows, "x", "y", expected_seeds=[1, 2, 3])
        self.assertEqual(first, second)
        self.assertAlmostEqual(first["difference_a_minus_b"], 2 / 3)
        self.assertEqual(first["n_families"], 3)
        self.assertEqual(first["n_complexes"], 5)
        self.assertEqual(first["confidence_interval"], [0, 1])

    def test_one_family_does_not_imply_certain_generalization(self):
        rows = [run("A", "a", "x", 1, 1), run("A", "a", "y", 1, 0)]
        report = paired_family_bootstrap(rows, "x", "y")
        self.assertEqual(report["difference_a_minus_b"], 1)
        self.assertIsNone(report["confidence_interval"])
        self.assertEqual(report["n_resamples"], 0)

    def test_expected_seed_iterators_are_not_exhausted(self):
        rows = [run("A", "a", "x", 1, 1), run("A", "a", "y", 1, 0)]
        report = paired_family_bootstrap(rows, "x", "y", expected_seeds=iter([1]))
        self.assertEqual(report["difference_a_minus_b"], 1)


class SpecificityTests(unittest.TestCase):
    def test_ranking_and_tie_controls(self):
        perfect = binary_ranking_metrics([1, 0, 1, 0], [4, 1, 3, 2])
        self.assertEqual(perfect["auroc"], 1)
        self.assertEqual(perfect["auprc_average_precision"], 1)
        reverse = binary_ranking_metrics([1, 1, 0, 0], [1, 2, 3, 4])
        self.assertEqual(reverse["auroc"], 0)
        ties = binary_ranking_metrics([1, 0, 1, 0], [2, 2, 2, 2])
        self.assertEqual(ties["auroc"], 0.5)
        self.assertEqual(ties["auprc_average_precision"], 0.5)
        self.assertEqual(ties["prevalence"], 0.5)
        with self.assertRaises(ValueError):
            binary_ranking_metrics([1, 1], [1, 2])

    def test_contexts_never_pool_and_untested_not_negative(self):
        rows = [dict(context_id="c1", target_id="t", ligand_id="a", label="binder", score=1),
                dict(context_id="c1", target_id="t", ligand_id="b", label="nonbinder", score=0),
                dict(context_id="c1", target_id="t", ligand_id="c", label="untested"),
                dict(context_id="c2", target_id="t", ligand_id="a", label="binder", score=0),
                dict(context_id="c2", target_id="t", ligand_id="b", label="nonbinder", score=1)]
        report = specificity_by_context(rows)
        self.assertIsNone(report["pooled_assay_metrics"])
        self.assertEqual([c["raw"]["auroc"] for c in report["contexts"]], [1, 0])
        self.assertEqual(report["contexts"][0]["untested_ligands"], 1)
        self.assertEqual(report["contexts"][0]["raw"]["tested_ligands"], 2)

    def test_single_class_returns_no_metrics(self):
        report = specificity_by_context([dict(context_id="c", target_id="t", ligand_id="a", label="binder", score=1)])
        self.assertEqual(report["contexts"][0]["status"], "insufficient_tested_classes")
        self.assertIsNone(report["contexts"][0]["raw"])

    def test_disjoint_prespecified_matching_retains_raw(self):
        rows = [dict(context_id="c", target_id="t", ligand_id="a", label="binder", score=1,
                     heavy_atoms=12, formal_charge=-1),
                dict(context_id="c", target_id="t", ligand_id="b", label="nonbinder", score=2,
                     heavy_atoms=12, formal_charge=-1),
                dict(context_id="c", target_id="t", ligand_id="d", label="nonbinder", score=0,
                     heavy_atoms=24, formal_charge=-2)]
        pair = dict(context_id="c", binder_id="a", nonbinder_id="b")
        context = specificity_by_context(rows, matched_pairs=[pair])["contexts"][0]
        self.assertEqual(context["raw"]["auroc"], 0.5)
        self.assertEqual(context["matched"]["auroc"], 0)
        self.assertEqual(context["matched"]["pair_concordance"], 0)
        with self.assertRaises(ValueError):
            specificity_by_context(rows, matched_pairs=[pair, pair])
        with self.assertRaises(ValueError):
            specificity_by_context(rows, matched_pairs=[dict(pair, nonbinder_id="d")])

    def test_invalid_context_or_missing_tested_score_fails(self):
        row = dict(context_id="c", target_id="t", ligand_id="a", label="binder", score=1)
        for invalid in (dict(row, score=None), dict(row, label="negative")):
            with self.assertRaises(ValueError):
                specificity_by_context([invalid])
        with self.assertRaises(ValueError):
            specificity_by_context([row, dict(row, target_id="other", ligand_id="b")])


class QuantitativeSpecificityTests(unittest.TestCase):
    def rows(self, scores=(1, 2, 3), kd=(100, 10, 1), context="c"):
        return [dict(context_id=context, target_id="t", ligand_id=str(index),
                     kd=affinity, kd_unit="uM", score=score)
                for index, (score, affinity) in enumerate(zip(scores, kd))]

    def test_perfect_and_reversed_affinity_ranking(self):
        report = quantitative_specificity_by_context(self.rows())["contexts"][0]
        self.assertEqual(report["spearman_rho"], 1)
        self.assertEqual(report["kendall_tau_b"], 1)
        self.assertEqual(report["concordant_pairs"], 3)
        reverse = quantitative_specificity_by_context(self.rows(scores=(3, 2, 1)))["contexts"][0]
        self.assertEqual(reverse["spearman_rho"], -1)
        self.assertEqual(reverse["kendall_tau_b"], -1)
        self.assertIsNone(report["p_values"])
        self.assertIsNone(report["confidence_intervals"])

    def test_average_tied_ranks_and_tau_b(self):
        # Only the first score pair ties; no affinity pair ties.
        report = quantitative_specificity_by_context(self.rows(scores=(1, 1, 2)))["contexts"][0]
        self.assertAlmostEqual(report["spearman_rho"], math.sqrt(3) / 2)
        self.assertAlmostEqual(report["kendall_tau_b"], 2 / math.sqrt(6))
        self.assertEqual(report["score_only_ties"], 1)
        self.assertEqual(report["affinity_only_ties"], 0)
        both = quantitative_specificity_by_context(self.rows(scores=(1, 1, 2), kd=(100, 100, 1)))["contexts"][0]
        self.assertEqual(both["spearman_rho"], 1)
        self.assertEqual(both["kendall_tau_b"], 1)
        self.assertEqual(both["both_ties"], 1)

    def test_mixed_explicit_units_are_normalized_before_ranking(self):
        rows = self.rows(kd=(1, 10, 100))
        rows[0]["kd_unit"] = "mM"  # 1 mM is weaker than 10 uM, then 100 nM.
        rows[2]["kd_unit"] = "nM"
        report = quantitative_specificity_by_context(rows)["contexts"][0]
        self.assertEqual(report["spearman_rho"], 1)
        self.assertEqual([r["pkd"] for r in report["ligands"]], [3, 5, 7])

    def test_contexts_are_separate_even_for_identical_ligand_ids(self):
        rows = self.rows(context="a") + self.rows(scores=(3, 2, 1), context="b")
        report = quantitative_specificity_by_context(rows)
        self.assertIsNone(report["pooled_assay_metrics"])
        self.assertEqual([r["spearman_rho"] for r in report["contexts"]], [1, -1])

    def test_single_and_constant_variables_have_no_correlation(self):
        single = quantitative_specificity_by_context(self.rows()[:1])["contexts"][0]
        self.assertEqual(single["status"], "insufficient_ligands")
        self.assertIsNone(single["spearman_rho"])
        self.assertIsNone(single["kendall_tau_b"])
        constant = quantitative_specificity_by_context(self.rows(scores=(1, 1, 1)))["contexts"][0]
        self.assertEqual(constant["status"], "undefined_constant_ranks")
        self.assertIsNone(constant["spearman_rho"])
        self.assertIsNone(constant["kendall_tau_b"])

    def test_invalid_exact_affinities_and_censoring_fail(self):
        row = self.rows()[0]
        for patch in ({"kd": 0}, {"kd": -1}, {"kd": math.inf}, {"kd": math.nan},
                      {"kd": "<1"}, {"kd": None}, {"censored": True},
                      {"censored": "false"}, {"kd_relation": ">"}, {"kd_unit": "RFU"},
                      {"score": math.nan}):
            with self.subTest(patch=patch), self.assertRaises(ValueError):
                quantitative_specificity_by_context([dict(row, **patch)])
        with self.assertRaisesRegex(ValueError, "censored KD values are unsupported"):
            quantitative_specificity_by_context([dict(row, kd_relation="<")])
        with self.assertRaises(ValueError):
            quantitative_specificity_by_context([row, row])
        with self.assertRaises(ValueError):
            quantitative_specificity_by_context([row, dict(row, target_id="other", ligand_id="different")])


if __name__ == "__main__":
    unittest.main()
