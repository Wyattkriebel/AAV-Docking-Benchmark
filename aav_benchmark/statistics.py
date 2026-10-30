"""Family-aware benchmark summaries and assay-context ranking.

Every submitted run is in the eligible denominator. Failed and unsupported runs
count as zero recovery, never as exclusions. Callers must enumerate every frozen
eligible run, including failures; ``expected_seeds`` verifies its seed coverage.
Conditions (information regime, receptor state, etc.) must be filtered separately
or included in the method identifier. No function pools independent assays.
"""
from __future__ import annotations

from collections import defaultdict
import math
import random
from statistics import mean
from typing import Iterable, Mapping, Sequence


_STATUSES = {"ok", "failed", "unsupported"}
_LABELS = {"binder", "nonbinder", "untested"}


def _identifier(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def _probability(value, name):
    if isinstance(value, bool):
        return float(value)
    value = _finite(value, name)
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be between zero and one")
    return value


def _seed(value):
    if isinstance(value, bool) or not isinstance(value, (int, str)) or not str(value).strip():
        raise ValueError("seed must be a nonempty string or integer")
    return str(value)


def _validated_runs(rows, metric, expected_seeds):
    rows = [dict(row) for row in rows]
    if not rows:
        raise ValueError("at least one eligible run is required")
    expected = None
    if expected_seeds is not None:
        original = list(expected_seeds)
        expected = {_seed(s) for s in original}
        if not expected or len(expected) != len(original):
            raise ValueError("expected_seeds must be nonempty and unique")
    seen, families, seeds = set(), {}, defaultdict(set)
    for row in rows:
        for name in ("family_id", "complex_id", "method"):
            row[name] = _identifier(row.get(name), name)
        row["seed"] = _seed(row.get("seed"))
        key = (row["complex_id"], row["method"], row["seed"])
        if key in seen:
            raise ValueError(f"duplicate eligible run: {key}")
        seen.add(key)
        if row["complex_id"] in families and families[row["complex_id"]] != row["family_id"]:
            raise ValueError("a complex cannot belong to multiple families")
        families[row["complex_id"]] = row["family_id"]
        status = row.get("status")
        if status not in _STATUSES:
            raise ValueError(f"status must be one of {sorted(_STATUSES)}")
        if status == "ok":
            row[metric] = _probability(row.get(metric), metric)
        else:
            if row.get(metric) is not None and _probability(row[metric], metric) != 0:
                raise ValueError("failed/unsupported runs cannot have nonzero recovery")
            row[metric] = 0.0
        if row.get("near_native_present") is not None:
            if not isinstance(row["near_native_present"], bool):
                raise ValueError("near_native_present must be boolean")
            if status != "ok" and row["near_native_present"]:
                raise ValueError("a failed run cannot have a returned near-native pose")
            if row[metric] > 0 and not row["near_native_present"]:
                raise ValueError("recovery cannot succeed without a near-native returned pose")
        seeds[(row["complex_id"], row["method"])].add(row["seed"])
    if expected is not None:
        for key, actual in seeds.items():
            if actual != expected:
                raise ValueError(f"seed coverage differs from expected_seeds for {key}: {sorted(actual)}")
    return rows


def _optional_mean(values):
    values = [value for value in values if value is not None]
    return mean(values) if values else None


def aggregate_runs(rows: Iterable[Mapping], metric: str = "success", *, expected_seeds=None):
    """Average seeds within complexes, then complexes within families.

    ``metric`` must be a recovery probability/indicator in [0,1]. Return records
    per complex and family and an overall family-macro rate for each method.
    Conditional rates use only their stated subsets; their denominators are
    included. An unknown near-native flag is never silently treated as false.
    This function verifies supplied rows, not eligibility against an absent
    cohort ledger: callers must join their frozen ledger to execution statuses.
    """
    rows = _validated_runs(rows, metric, expected_seeds)
    groups = defaultdict(list)
    for row in rows:
        groups[(row["method"], row["family_id"], row["complex_id"])].append(row)
    complexes = []
    for (method, family, complex_id), runs in sorted(groups.items()):
        completed = [r for r in runs if r["status"] == "ok"]
        flags_complete = all(r.get("near_native_present") is not None for r in completed)
        available = [r for r in completed if r.get("near_native_present") is True]
        complexes.append({
            "method": method, "family_id": family, "complex_id": complex_id,
            "eligible_runs": len(runs), "completed_runs": len(completed),
            "failed_runs": sum(r["status"] == "failed" for r in runs),
            "unsupported_runs": sum(r["status"] == "unsupported" for r in runs),
            "recovery": mean(r[metric] for r in runs),
            "conditional_completed_recovery": _optional_mean(r[metric] for r in completed),
            "near_native_runs": len(available) if flags_complete else None,
            "conditional_near_native_recovery": (
                _optional_mean(r[metric] for r in available) if flags_complete else None
            ),
        })
    family_groups = defaultdict(list)
    for row in complexes:
        family_groups[(row["method"], row["family_id"])].append(row)
    families = []
    for (method, family), systems in sorted(family_groups.items()):
        families.append({
            "method": method, "family_id": family, "complexes": len(systems),
            "recovery": mean(s["recovery"] for s in systems),
            "conditional_completed_recovery": _optional_mean(s["conditional_completed_recovery"] for s in systems),
            "conditional_near_native_recovery": _optional_mean(s["conditional_near_native_recovery"] for s in systems),
            "complexes_with_completed_runs": sum(s["completed_runs"] > 0 for s in systems),
            "complexes_with_near_native_runs": sum((s["near_native_runs"] or 0) > 0 for s in systems),
            "near_native_flags_complete": all(s["near_native_runs"] is not None for s in systems),
        })
    overall = []
    for method in sorted({r["method"] for r in rows}):
        method_families = [f for f in families if f["method"] == method]
        method_complexes = [s for s in complexes if s["method"] == method]
        flags_complete = all(s["near_native_runs"] is not None for s in method_complexes)
        overall.append({
            "method": method, "families": len(method_families), "complexes": len(method_complexes),
            **{name: sum(s[name] for s in method_complexes) for name in
               ("eligible_runs", "completed_runs", "failed_runs", "unsupported_runs")},
            "family_macro_recovery": mean(f["recovery"] for f in method_families),
            "complex_macro_recovery": mean(s["recovery"] for s in method_complexes),
            "conditional_completed_recovery": _optional_mean(f["conditional_completed_recovery"] for f in method_families),
            "families_with_completed_runs": sum(f["complexes_with_completed_runs"] > 0 for f in method_families),
            "near_native_runs": sum(s["near_native_runs"] for s in method_complexes) if flags_complete else None,
            "conditional_near_native_recovery": (
                _optional_mean(f["conditional_near_native_recovery"] for f in method_families)
                if flags_complete else None
            ),
            "families_with_near_native_runs": (
                sum(f["complexes_with_near_native_runs"] > 0 for f in method_families)
                if flags_complete else None
            ),
        })
    return {"metric": metric, "weighting": "seed -> complex -> family; equal weight at each level",
            "per_complex": complexes, "per_family": families, "overall": overall}


def _quantile(values, probability):
    values = sorted(values)
    position = (len(values) - 1) * probability
    lo = math.floor(position)
    hi = math.ceil(position)
    return values[lo] + (values[hi] - values[lo]) * (position - lo)


def paired_family_bootstrap(rows: Iterable[Mapping], method_a: str, method_b: str,
                            metric: str = "success", *, n_resamples: int = 10000,
                            seed: int = 20261028, confidence: float = 0.95,
                            expected_seeds=None):
    """Paired family-macro difference (A minus B), with percentile bootstrap CI.

    Exact complex and seed pairing is required. A failure remains eligible and
    contributes zero. Resample whole families after averaging seeds and systems;
    never bootstrap individual poses or seeds as biological replicates. One
    family has no estimable between-family uncertainty, so its CI is null.
    """
    _identifier(method_a, "method_a")
    _identifier(method_b, "method_b")
    if method_a == method_b:
        raise ValueError("paired methods must differ")
    if isinstance(n_resamples, bool) or not isinstance(n_resamples, int) or n_resamples < 1:
        raise ValueError("n_resamples must be a positive integer")
    confidence = _finite(confidence, "confidence")
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between zero and one")
    # Validation is repeated by aggregate_runs below; retain generator inputs.
    if expected_seeds is not None:
        expected_seeds = list(expected_seeds)
    rows = _validated_runs(rows, metric, expected_seeds)
    selected = [r for r in rows if r["method"] in (method_a, method_b)]
    keys = {method: {(r["family_id"], r["complex_id"], r["seed"]) for r in selected if r["method"] == method}
            for method in (method_a, method_b)}
    if not keys[method_a] or keys[method_a] != keys[method_b]:
        raise ValueError("paired comparison requires identical complex/seed keys for both methods")
    summary = aggregate_runs(selected, metric, expected_seeds=expected_seeds)
    family_values = defaultdict(dict)
    for row in summary["per_family"]:
        family_values[row["family_id"]][row["method"]] = row["recovery"]
    differences = [values[method_a] - values[method_b] for _, values in sorted(family_values.items())]
    ci = None
    if len(differences) > 1:
        rng = random.Random(seed)
        samples = [mean(rng.choices(differences, k=len(differences))) for _ in range(n_resamples)]
        alpha = (1 - confidence) / 2
        ci = [_quantile(samples, alpha), _quantile(samples, 1 - alpha)]
    return {"method_a": method_a, "method_b": method_b, "metric": metric,
            "difference_a_minus_b": mean(differences), "confidence": confidence,
            "confidence_interval": ci, "n_families": len(differences),
            "n_complexes": len({r["complex_id"] for r in selected}),
            "n_resamples": n_resamples if ci is not None else 0, "bootstrap_seed": seed,
            "interval_method": "family-cluster percentile bootstrap",
            "uncertainty_note": ("Between-family uncertainty is not estimable from one family."
                                 if ci is None else "Descriptive interval; few families can yield unstable coverage."),
            "per_family": [{"family_id": family, "recovery_a": values[method_a],
                            "recovery_b": values[method_b], "difference": values[method_a] - values[method_b]}
                           for family, values in sorted(family_values.items())]}


def binary_ranking_metrics(labels: Sequence, scores: Sequence):
    """AUROC and step-wise average precision (AUPRC), higher score = binder.

    Ties are processed as a single threshold. Both measured classes are required.
    Average precision is the sum of precision times each recall increment; it is
    not trapezoidal interpolation of the precision-recall curve.
    """
    labels, scores = list(labels), list(scores)
    if not labels or len(labels) != len(scores):
        raise ValueError("labels and scores must have equal nonzero length")
    if any(label not in (0, 1, False, True) for label in labels):
        raise ValueError("binary labels must be zero or one")
    labels = [int(label) for label in labels]
    scores = [_finite(s, "score") for s in scores]
    n_positive = sum(labels)
    n_negative = len(labels) - n_positive
    if not n_positive or not n_negative:
        raise ValueError("AUROC/AUPRC require tested binders and tested nonbinders")
    # Mann-Whitney statistic: positive score wins, ties receive half credit.
    auc = sum((a > b) + 0.5 * (a == b)
              for label, a in zip(labels, scores) if label
              for other, b in zip(labels, scores) if not other) / (n_positive * n_negative)
    thresholds = defaultdict(lambda: [0, 0])
    for label, score in zip(labels, scores):
        thresholds[score][0] += label
        thresholds[score][1] += 1
    true_positive = selected = 0
    ap = 0.0
    for _, (new_positive, new_total) in sorted(thresholds.items(), reverse=True):
        true_positive += new_positive
        selected += new_total
        ap += (new_positive / n_positive) * (true_positive / selected)
    return {"auroc": auc, "auprc_average_precision": ap,
            "prevalence": n_positive / len(labels), "tested_ligands": len(labels),
            "binders": n_positive, "nonbinders": n_negative}


def specificity_by_context(rows: Iterable[Mapping], *, matched_pairs=None,
                           heavy_atom_tolerance: int = 0, charge_tolerance: float = 0):
    """Evaluate tested binding labels separately within each assay/target context.

    Rows require context_id, target_id, ligand_id, label and (for tested rows) a
    finite score, oriented so higher predicts stronger binding. Untested ligands
    are counted but excluded from metrics. ``context_id`` must identify the
    complete assay setting (presentation, concentration, valency, salt, etc.).

    Optional prespecified, disjoint matching records have context_id, binder_id,
    nonbinder_id. They require heavy_atoms and formal_charge on matched rows and
    satisfy the explicitly supplied tolerances (default: exact size and charge).
    Both raw and matched metrics are returned; matching never overwrites raw.
    Callers must freeze pairs before inspecting scores; this cannot be inferred
    from the numbers. No threshold or classifier is fitted here.
    """
    if isinstance(heavy_atom_tolerance, bool) or not isinstance(heavy_atom_tolerance, int) or heavy_atom_tolerance < 0:
        raise ValueError("heavy_atom_tolerance must be a nonnegative integer")
    charge_tolerance = _finite(charge_tolerance, "charge_tolerance")
    if charge_tolerance < 0:
        raise ValueError("charge_tolerance cannot be negative")
    rows = [dict(row) for row in rows]
    if not rows:
        raise ValueError("at least one assay record is required")
    groups, lookup, targets = defaultdict(list), {}, {}
    for row in rows:
        for name in ("context_id", "target_id", "ligand_id"):
            row[name] = _identifier(row.get(name), name)
        context, ligand = row["context_id"], row["ligand_id"]
        if context in targets and targets[context] != row["target_id"]:
            raise ValueError("a context cannot mix protein targets")
        targets[context] = row["target_id"]
        if (context, ligand) in lookup:
            raise ValueError("duplicate ligand within assay context")
        if row.get("label") not in _LABELS:
            raise ValueError(f"label must be one of {sorted(_LABELS)}")
        if row["label"] != "untested":
            row["score"] = _finite(row.get("score"), "score")
        lookup[(context, ligand)] = row
        groups[context].append(row)
    matched, used = defaultdict(list), set()
    for pair in ([] if matched_pairs is None else matched_pairs):
        context = _identifier(pair.get("context_id"), "context_id")
        binder_id = _identifier(pair.get("binder_id"), "binder_id")
        nonbinder_id = _identifier(pair.get("nonbinder_id"), "nonbinder_id")
        try:
            binder, nonbinder = lookup[(context, binder_id)], lookup[(context, nonbinder_id)]
        except KeyError as error:
            raise ValueError("matched pair references an absent context/ligand") from error
        if binder["label"] != "binder" or nonbinder["label"] != "nonbinder":
            raise ValueError("matched pairs require one tested binder and one tested nonbinder")
        for row in (binder, nonbinder):
            key = (context, row["ligand_id"])
            if key in used:
                raise ValueError("matched pairs must be disjoint within each context")
            used.add(key)
            atoms = row.get("heavy_atoms")
            if isinstance(atoms, bool) or not isinstance(atoms, int) or atoms < 1:
                raise ValueError("matched ligands require positive integer heavy_atoms")
            row["formal_charge"] = _finite(row.get("formal_charge"), "formal_charge")
        if abs(binder["heavy_atoms"] - nonbinder["heavy_atoms"]) > heavy_atom_tolerance:
            raise ValueError("matched pair exceeds heavy-atom tolerance")
        if abs(binder["formal_charge"] - nonbinder["formal_charge"]) > charge_tolerance:
            raise ValueError("matched pair exceeds formal-charge tolerance")
        matched[context].append((binder, nonbinder))
    result = []
    for context, context_rows in sorted(groups.items()):
        tested = [row for row in context_rows if row["label"] != "untested"]
        classes = {row["label"] for row in tested}
        raw = (binary_ranking_metrics([r["label"] == "binder" for r in tested], [r["score"] for r in tested])
               if len(classes) == 2 else None)
        pairs = matched[context]
        matched_metrics = None
        if pairs:
            subset = [row for pair in pairs for row in pair]
            matched_metrics = binary_ranking_metrics([r["label"] == "binder" for r in subset], [r["score"] for r in subset])
            matched_metrics.update({"pairs": len(pairs), "pair_concordance": mean(
                float(a["score"] > b["score"]) + 0.5 * (a["score"] == b["score"]) for a, b in pairs),
                "heavy_atom_tolerance": heavy_atom_tolerance, "charge_tolerance": charge_tolerance,
                "ligand_ids": [r["ligand_id"] for r in subset]})
        result.append({"context_id": context, "target_id": targets[context],
                       "status": "evaluated" if raw is not None else "insufficient_tested_classes",
                       "tested_ligands": len(tested), "untested_ligands": len(context_rows) - len(tested),
                       "raw": raw, "matched": matched_metrics})
    return {"score_direction": "higher predicts binding", "auprc_definition": "step-wise average precision",
            "pooled_assay_metrics": None, "contexts": result}


_KD_SCALES = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "μM": 1e-6,
              "nM": 1e-9, "pM": 1e-12, "fM": 1e-15}


def _average_ranks(values):
    """Ascending ranks, with each tied group assigned its average rank."""
    ranks = [0.0] * len(values)
    ordered = sorted(range(len(values)), key=values.__getitem__)
    start = 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and values[ordered[stop]] == values[ordered[start]]:
            stop += 1
        average_rank = (start + 1 + stop) / 2
        for position in ordered[start:stop]:
            ranks[position] = average_rank
        start = stop
    return ranks


def _spearman(values_a, values_b):
    ranks_a, ranks_b = _average_ranks(values_a), _average_ranks(values_b)
    center_a, center_b = mean(ranks_a), mean(ranks_b)
    covariance = sum((a - center_a) * (b - center_b) for a, b in zip(ranks_a, ranks_b))
    variance_a = sum((a - center_a) ** 2 for a in ranks_a)
    variance_b = sum((b - center_b) ** 2 for b in ranks_b)
    denominator = math.sqrt(variance_a * variance_b)
    return covariance / denominator if denominator else None


def _kendall_tau_b(values_a, values_b):
    concordant = discordant = score_ties = affinity_ties = both_ties = 0
    for i in range(len(values_a)):
        for j in range(i + 1, len(values_a)):
            delta_a = (values_a[i] > values_a[j]) - (values_a[i] < values_a[j])
            delta_b = (values_b[i] > values_b[j]) - (values_b[i] < values_b[j])
            if not delta_a and not delta_b:
                both_ties += 1
            elif not delta_a:
                score_ties += 1
            elif not delta_b:
                affinity_ties += 1
            elif delta_a == delta_b:
                concordant += 1
            else:
                discordant += 1
    comparable = concordant + discordant
    denominator = math.sqrt((comparable + score_ties) * (comparable + affinity_ties))
    return {"kendall_tau_b": (concordant - discordant) / denominator if denominator else None,
            "concordant_pairs": concordant, "discordant_pairs": discordant,
            "score_only_ties": score_ties, "affinity_only_ties": affinity_ties,
            "both_ties": both_ties}


def quantitative_specificity_by_context(rows: Iterable[Mapping]):
    """Descriptive within-assay score association with exact measured KD.

    Rows require context_id, target_id, ligand_id, score (higher predicts stronger
    binding), kd (strictly positive finite number) and kd_unit. Supported units
    are M, mM, uM/µM/μM, nM, pM and fM; all are normalized before calculating
    pKD = -log10(KD in M). No KD is inferred from a docking energy.

    ``censored=True`` or ``kd_relation`` other than '='/'exact' is explicitly
    unsupported: censored bounds are not treated as exact affinities. The caller
    must use a context_id that preserves presentation, protein construct, assay
    and solution conditions. Duplicate ligands and mixed targets fail visibly.

    Spearman uses average tied ranks and Kendall tau-b corrects for ties. One
    ligand or a constant variable gives null correlations with explicit status.
    No p-values, confidence intervals, classifier labels or cross-assay pooled
    metrics are produced. Small panels are exploratory descriptions.
    """
    rows = [dict(row) for row in rows]
    if not rows:
        raise ValueError("at least one quantitative assay record is required")
    groups, seen, targets = defaultdict(list), set(), {}
    for row in rows:
        for name in ("context_id", "target_id", "ligand_id"):
            row[name] = _identifier(row.get(name), name)
        context = row["context_id"]
        if context in targets and targets[context] != row["target_id"]:
            raise ValueError("a quantitative context cannot mix protein targets")
        targets[context] = row["target_id"]
        key = (context, row["ligand_id"])
        if key in seen:
            raise ValueError("duplicate ligand within quantitative assay context")
        seen.add(key)
        if not isinstance(row.get("censored", False), bool):
            raise ValueError("censored must be a boolean")
        if row.get("censored", False) or row.get("kd_relation", "=") not in ("=", "exact"):
            raise ValueError("censored KD values are unsupported; exact measured affinities are required")
        kd = _finite(row.get("kd"), "kd")
        if kd <= 0:
            raise ValueError("kd must be strictly positive")
        unit = row.get("kd_unit")
        if not isinstance(unit, str) or unit not in _KD_SCALES:
            raise ValueError(f"kd_unit must be one of {sorted(_KD_SCALES)}")
        row["score"] = _finite(row.get("score"), "score")
        # Separate logarithms avoid overflow/underflow during unit conversion.
        row["pkd"] = -math.log10(kd) - math.log10(_KD_SCALES[unit])
        groups[context].append(row)
    results = []
    for context, context_rows in sorted(groups.items()):
        scores = [row["score"] for row in context_rows]
        affinities = [row["pkd"] for row in context_rows]
        spearman = _spearman(scores, affinities)
        kendall = _kendall_tau_b(scores, affinities)
        status = "evaluated_descriptive"
        if len(context_rows) < 2:
            status = "insufficient_ligands"
        elif spearman is None or kendall["kendall_tau_b"] is None:
            status = "undefined_constant_ranks"
        results.append({"context_id": context, "target_id": targets[context],
                        "status": status, "n_ligands": len(context_rows),
                        "spearman_rho": spearman, **kendall,
                        "p_values": None, "confidence_intervals": None,
                        "ligands": [{"ligand_id": row["ligand_id"], "kd": row["kd"],
                                     "kd_unit": row["kd_unit"], "pkd": row["pkd"], "score": row["score"]}
                                    for row in context_rows]})
    return {"score_direction": "higher predicts stronger binding",
            "affinity_transform": "-log10(KD in M), measured affinities only",
            "interpretation": "exploratory within-assay association; no significance or generalization claim",
            "measurement_uncertainty": "correlations use reported point estimates; experimental uncertainty is not propagated",
            "pooled_assay_metrics": None, "contexts": results}
