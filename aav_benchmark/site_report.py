"""Render the site-localization report from the frozen evaluated experiment.

This renderer checks the fixed denominator, source digests and the correspondence
between each reported row and its per-replicate evaluation. Raw prediction
reparsing and receptor/reference geometry recomputation are separate checks.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import mean

from .provenance import sha256, write_text_atomic


METHODS = ("fpocket", "p2rank", "random_surface")
LABELS = {"fpocket": "fpocket", "p2rank": "P2Rank", "random_surface": "Surface baseline"}
FIELDS = [f"top{k}_success_{t}A" for k in (1, 5) for t in (3, 4, 5)] + [
    "top1_common_interface_recall", "top1_method_interface_recall",
    "top5_common_union_interface_recall", "top5_method_union_interface_recall"]


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _safe(root, name):
    path = Path(name)
    _require(not path.is_absolute() and ".." not in path.parts, "Evidence path leaves repository")
    return root / path


def _json(path):
    return json.loads(Path(path).read_text())


def _same(actual, expected, label):
    if isinstance(expected, (float, int)) and not isinstance(expected, bool):
        _require(isinstance(actual, (float, int)) and not isinstance(actual, bool)
                 and math.isfinite(actual) and math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12),
                 f"Evaluated row disagrees with detailed evidence: {label}")
    else:
        _require(actual == expected, f"Evaluated row disagrees with detailed evidence: {label}")


def _check_metrics(row, trials, *, baseline, label):
    _require(len(trials) == (100 if baseline else 1), f"Wrong replicate denominator: {label}")
    for trial in trials:
        for field in ("raw_candidate_count", "returned_distinct_candidates"):
            _require(type(trial[field]) is int and trial[field] >= 0,
                     f"Invalid candidate count: {label}/{field}")
        _require(trial["returned_distinct_candidates"] <= trial["raw_candidate_count"],
                 f"Merged candidate count exceeds input: {label}")
        for field in FIELDS:
            value = trial[field]
            _require(isinstance(value, (float, int)) and math.isfinite(value) and 0 <= value <= 1,
                     f"Invalid evaluated endpoint: {label}/{field}")
            if "success" in field:
                _require(value in (0, 1), f"Individual success must be binary: {label}/{field}")
        for budget in (1, 5):
            values = [trial[f"top{budget}_success_{t}A"] for t in (3, 4, 5)]
            _require(values == sorted(values), f"Threshold sensitivity is nonmonotone: {label}")
        for threshold in (3, 4, 5):
            _require(trial[f"top1_success_{threshold}A"] <= trial[f"top5_success_{threshold}A"],
                     f"Candidate-budget sensitivity is nonmonotone: {label}")
    for field in FIELDS:
        _same(row[field], mean(float(trial[field]) for trial in trials), f"{label}/{field}")
    for field, detailed in (("raw_candidate_count", "raw_candidate_count"),
                            ("distinct_candidate_count", "returned_distinct_candidates")):
        _same(row[field], mean(trial[detailed] for trial in trials), f"{label}/{field}")
    for field in ("top1_distance_A", "top5_min_distance_A", "first_matching_class_4A"):
        if field in row:
            _same(row[field], None if baseline else trials[0][field], f"{label}/{field}")


def _macro(rows):
    return {method: {field: mean(row[field] for row in rows if row["method"] == method)
                     for field in FIELDS} for method in METHODS}


def _pct(value):
    return f"{100 * value:.2f}%"


def _pair(row, threshold=4):
    return f"{_pct(row[f'top1_success_{threshold}A'])} / {_pct(row[f'top5_success_{threshold}A'])}"


def generate_site_report(root):
    """Return a validated summary and deterministic Markdown; never run predictors."""
    root = Path(root).resolve()
    directory = root / "results/site_discovery"
    config_path = root / "benchmarks/site_experiment.json"
    config, evaluation = _json(config_path), _json(directory / "evaluation.json")
    receipt = _json(directory / "freeze_receipt.json")
    _require(sha256(config_path) == receipt["configuration_sha256"], "Frozen site configuration changed")
    _require(receipt["parameters_frozen_before_scientific_predictions"] is True, "Missing pre-execution freeze")
    _require(sha256(_safe(root, config["reference_ledger"])) == config["reference_ledger_sha256"],
             "Frozen site reference eligibility changed")
    for name, digest in config["method_locks"].items():
        _require(sha256(_safe(root, name)) == digest, "Frozen site method lock changed")
    cases = {case["condition_id"]: case for case in config["cases"]}
    _require(len(cases) == len(config["cases"]) == 9 and len({c["complex_id"] for c in cases.values()}) == 8,
             "Expected eight systems and nine paired/unpaired conditions")
    _require(tuple(config["methods"]) == METHODS and config["baseline"]["replicates"] == 100,
             "Method or Monte Carlo denominator changed")
    _require(evaluation["complete"] is True and evaluation["pending"] == []
             and evaluation["eligible_complexes"] == 8 and evaluation["eligible_conditions"] == 9
             and evaluation["eligible_method_conditions"] == 27, "Site evaluation is incomplete")
    sources = evaluation["source_sha256"]
    for name, digest in sources.items():
        _require(sha256(_safe(root, name)) == digest, f"Site evaluation source changed: {name}")
    expected = {(condition, method) for condition in cases for method in METHODS}
    rows = evaluation["rows"]
    actual = {(row["condition_id"], row["method"]) for row in rows}
    _require(len(rows) == len(actual) == 27 and actual == expected, "Duplicate, missing or extra method-condition row")
    by_key = {(row["condition_id"], row["method"]): row for row in rows}
    merge_distances = [str(float(value)) for value in config["merge_sensitivity_A"]]
    merge_complete = True
    for condition, case in cases.items():
        details_name = f"results/site_discovery/{condition}/evaluation_details.json"
        _require(details_name in sources, f"Unbound detailed evidence: {condition}")
        details = _json(root / details_name)
        _require(details["condition_id"] == condition and details["complex_id"] == case["complex_id"],
                 "Detailed evidence identifies a different condition")
        for method in METHODS:
            row, evidence = by_key[condition, method], details["methods"][method]
            label = f"{condition}/{method}"
            _require(row["state"] == case["state"] and row["complex_id"] == case["complex_id"]
                     and row["eligible_conditions"] == 1, f"Changed cohort identity: {label}")
            _require(row["status"] == evidence["status"] and row["status"] in ("ok", "failed"),
                     f"Unfinished or inconsistent method status: {label}")
            relative = f"results/site_discovery/{condition}/" + {
                "fpocket": "fpocket_v2/execution.json", "p2rank": "p2rank/result.json",
                "random_surface": "random_surface_v2/result.json"}[method]
            _require(relative in sources and evidence["source_record_sha256"] == sources[relative],
                     f"Unbound normalized prediction evidence: {label}")
            normalized = _json(root / relative)
            _require(("ok" if normalized["status"] in ("ok", "completed") else normalized["status"]) == row["status"],
                     f"Prediction/evaluation status disagreement: {label}")
            baseline = method == "random_surface"
            _require(row["baseline_replicates"] == (100 if baseline else None), f"Changed replicate count: {label}")
            _check_metrics(row, evidence["evaluations"], baseline=baseline, label=label)
            if row["status"] == "failed":
                _require(all(row[field] == 0 for field in FIELDS) and row["raw_candidate_count"] == 0,
                         f"Failed method has fabricated predictions: {label}")
            for distance in merge_distances:
                if distance not in row.get("merge_sensitivity", {}) or distance not in evidence.get("merge_sensitivity", {}):
                    merge_complete = False
                    continue
                _check_metrics(row["merge_sensitivity"][distance], evidence["merge_sensitivity"][distance]["evaluations"],
                               baseline=baseline, label=f"{label}/merge{distance}")
    primary = [row for row in rows if row["state"] == "holo"]
    _require(len(primary) == 24 and len({row["complex_id"] for row in primary}) == 8,
             "Primary descriptive denominator must be eight holo systems")
    macro = _macro(primary)
    sensitivity = {distance: _macro([{**row, **row["merge_sensitivity"][distance]} for row in primary])
                   for distance in merge_distances} if merge_complete else {}
    original, revised = [], []
    for condition in cases:
        first = _json(directory / condition / "fpocket/execution.json")
        retry = _json(directory / condition / "fpocket_v2/execution.json")
        _require(first["status"] == "failed" and first["returncode"] == 1,
                 "Original serialization failure history changed")
        _require(retry["status"] in ("failed", "completed"), "Incomplete fpocket compatibility attempt")
        _require(("ok" if retry["status"] == "completed" else "failed") == by_key[condition, "fpocket"]["status"],
                 "fpocket execution/evaluation status disagreement")
        if retry["status"] == "failed":
            _require(retry["returncode"] == -5, "Unexpected failure category requires report review")
        original.append(first)
        revised.append(retry)
    failed = sum(record["status"] == "failed" for record in revised)
    amendment = _json(root / "benchmarks/site_geometry_amendment.json")
    _require(amendment["decision_before_any_site_metric_evaluation"] is True
             and amendment["effective_operator_max_CA_displacement_A"] == .1
             and amendment["unchanged_scientific_endpoints"]["primary_site_distance_A"] == 4,
             "Source-only symmetry amendment requires report review")
    summary = {"systems": 8, "conditions": 9, "method_conditions": 27,
               "primary_holo_systems": 8, "primary_macro": macro, "merge_sensitivity": sensitivity,
               "merge_sensitivity_complete": merge_complete, "fpocket_original_reader_failures": len(original),
               "fpocket_corrected_failures": failed, "rows": rows,
               "evaluation_sha256": sha256(directory / "evaluation.json")}
    lines = ["# Published-site localization results", "",
             "The eight-system development comparison found the published modeled site among the Top5 distinct candidate classes in "
             f"{sum(r['top5_success_4A'] for r in primary if r['method'] == 'fpocket'):g}/8 systems with fpocket and "
             f"{sum(r['top5_success_4A'] for r in primary if r['method'] == 'p2rank'):g}/8 with P2Rank. "
             "Neither method recovered the AAV9 site in its Top5. These are descriptive localization endpoints; "
             "the references qualify approximate published sites rather than validated atomic poses.", "",
             "## Fixed denominator and primary endpoint", "",
             "The frozen cohort contains eight systems: AAV9 and seven released glycan–protein complexes. "
             "AAV9 apo is a paired ninth condition, reported separately. Three methods give 27 eligible method–condition rows. "
             "Each primary macro average weights the eight holo systems equally and retains execution failures as zero workflow successes. "
             "A failed program supplies no evidence about biological binding or whether a site exists.", "",
             "Success means a candidate center lies within 4 Å of any heavy atom of a modeled reference ligand. "
             "Candidate ranks are merged at 4 Å using receptor-derived capsid symmetry, with budgets of one and five distinct classes. "
             "All native symmetry copies represent one capsid system. The surface baseline is the mean of 100 area-weighted draws "
             "with five separated patches per draw; these are conditional Monte Carlo null replicates, not biological observations.", "",
             "| Method | Top1 success | Top5 success | Top1 common interface recall | Top5 common union recall |",
             "|---|---:|---:|---:|---:|"]
    for method in METHODS:
        m = macro[method]
        lines.append(f"| {LABELS[method]} | {_pct(m['top1_success_4A'])} | {_pct(m['top5_success_4A'])} | "
                     f"{_pct(m['top1_common_interface_recall'])} | {_pct(m['top5_common_union_interface_recall'])} |")
    lines += ["", "Common interface recall uses protein heavy atoms within 8 Å of the candidate center and native protein–ligand contacts "
              "strictly within 4 Å. Top1 distance and coverage refer to the same closest native copy. For Top5, "
              "the union of symmetry-expanded candidate footprints is compared separately against each single native copy, and the maximum recall is reported. "
              "Method-specific footprints remain in the detailed data; their differing definitions are not pooled into the common comparison.", "",
              "## Per-condition outcomes", "", "Each success cell gives Top1 / Top5 at 4 Å. Baseline cells are replicate means.", "",
              "| Condition | fpocket | P2Rank | Surface baseline |",
              "|---|---:|---:|---:|"]
    for condition in cases:
        cells = []
        for method in METHODS:
            row = by_key[condition, method]
            cells.append("failed (0 / 0)" if row["status"] == "failed" else _pair(row))
        lines.append(f"| {condition} | " + " | ".join(cells) + " |")
    lines += ["", "P2Rank returned zero candidates for 1QFO and 2J1T; those completed executions remain zero-success observations. "
              "Reference quality and assembly limitations apply to every row. No lower-quality case was removed after seeing recovery.", "",
              "## Distance and candidate-merge sensitivity", "",
              "The following macro averages retain all eight holo systems. Each cell gives Top1 / Top5 success.", "",
              "| Method | 3 Å center threshold | 4 Å center threshold | 5 Å center threshold |",
              "|---|---:|---:|---:|"]
    for method in METHODS:
        lines.append(f"| {LABELS[method]} | " + " | ".join(_pair(macro[method], threshold) for threshold in (3, 4, 5)) + " |")
    if merge_complete:
        lines += ["", "The frozen 2/6 Å merge sensitivities use the same 4 Å center endpoint and retain the original candidate order.", "",
                  "| Method | Merge 2 Å | Primary merge 4 Å | Merge 6 Å |", "|---|---:|---:|---:|"]
        for method in METHODS:
            lines.append(f"| {LABELS[method]} | {_pair(sensitivity['2.0'][method])} | {_pair(macro[method])} | {_pair(sensitivity['6.0'][method])} |")
    else:
        lines += ["", "The configuration also declares 2/6 Å candidate-merge sensitivity. Those evaluations are pending; "
                  "the present distance sensitivity does not substitute for that analysis."]
    lines += ["", "## Paired AAV9 conditions", "",
              "P2Rank's first matching distinct class at 4 Å is rank "
              f"{by_key['aav9_holo', 'p2rank']['first_matching_class_4A']} for holo and "
              f"{by_key['aav9_apo', 'p2rank']['first_matching_class_4A']} for apo, beyond the primary Top5 budget. "
              f"The raw counts are {by_key['aav9_holo', 'p2rank']['raw_candidate_count']:g} and "
              f"{by_key['aav9_apo', 'p2rank']['raw_candidate_count']:g}, respectively, reduced to "
              f"{by_key['aav9_holo', 'p2rank']['distinct_candidate_count']:g} and "
              f"{by_key['aav9_apo', 'p2rank']['distinct_candidate_count']:g} symmetry-aware classes. "
              "Apo predictions are aligned to the holo evaluation frame with all matching receptor CA atoms; ligands are absent from that fit. "
              "The two states supply a paired robustness observation, without an additional independent system.", "",
              "## Execution history and development amendments", "",
              f"All {len(original)} initial fpocket attempts failed during mmCIF reading (exit 1). "
              "The compatibility retry appended a terminal `#` delimiter required by its bundled reader. "
              "Atom identities, coordinates, method parameters and the original attempts were preserved. "
              f"The corrected runs completed for {len(revised) - failed} small receptors; both full capsid runs terminated with return code −5 "
              "(SIGTRAP), yielding zero predictions. These are workflow failures. The deposited ordered capsid was never cropped to force completion.", "",
              "A source-only symmetry audit increased the full-capsid maximum CA displacement guard from 0.05 to 0.1 Å for both states "
              "before site endpoint scoring. The apo operators reached 0.073316 Å; the chain-fit guard remained 0.05 Å. "
              "The primary 4 Å site endpoint, candidate merge distance and 8 Å footprint remained fixed. "
              "Original failures and the separately hashed amendment are retained. This development history does not support a pristine prospective confirmatory claim.", "",
              "fpocket uses unmodified upstream wall-clock seeding; one corrected scientific execution per condition is retained. "
              "P2Rank uses fixed seed 42 and one thread. The methods use released configurations rather than equal CPU budgets. "
              "Archived raw outputs support result reconstruction without claiming bitwise equivalence of a new native execution.", "",
              "## Scope of inference", "",
              "Every predictor received the full supplied protein-only heavy-atom receptor, with all glycans, waters and metals removed uniformly. "
              "AAV9 used the full 60-subunit deposited ordered capsid. The seven external receptors are full released files, "
              "which can be isolated domains or incomplete biological assemblies. Holo structures retain binding-induced conformational information. "
              "Metal removal particularly limits interpretation for calcium-dependent lectins.", "",
              "The [reference audit](SITE_REFERENCE_AUDIT.md) records weak local density, incomplete glycans, chemistry discrepancies and assembly caveats. "
              "Coverage concerns modeled residues only. The development collection and P2Rank training overlap have not passed independent sequence/fold and calibration audits. "
              "No family-held-out estimate, experimental nonbinding annotation for other patches, atomic-pose accuracy or ligand specificity follows from these site scores.", "",
              "## Reproduction", "",
              "`results/site_discovery/evaluation.json` contains all 27 rows and hashes their detailed/reference records. "
              "Each condition retains normalized outputs, original and revised fpocket attempts, raw-artifact archives, frozen input archives and baseline draws. "
              "The raw-prediction verifier reparses archived native output; the experiment verifier recomputes reference geometry and scores. "
              "This report separately checks source digests, denominator identities and every displayed aggregate against per-replicate evidence.", "",
              "```sh", "python scripts/regenerate_site_report.py --check", "```", "",
              f"Evaluated artifact SHA256: `{summary['evaluation_sha256']}`.", ""]
    return summary, "\n".join(lines)


def write_site_report(root, *, markdown_output=None, check=False):
    root = Path(root).resolve()
    _, markdown = generate_site_report(root)
    output = Path(markdown_output) if markdown_output is not None else root / "docs/SITE_LOCALIZATION_RESULTS.md"
    if check:
        _require(output.is_file() and output.read_text() == markdown, f"Site report is stale: {output}")
    else:
        write_text_atomic(output, markdown)
    return output
