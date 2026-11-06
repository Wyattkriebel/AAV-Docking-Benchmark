"""Deterministic development report from archived run evidence.

The report describes development diagnostics. It deliberately cannot mark the
repository's original scientific purpose complete. Missing, stale or incomplete
completed-run evidence raises rather than shrinking an eligible denominator.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import mean

from .statistics import quantitative_specificity_by_context
from .provenance import write_text_atomic


_AAV_CONDITIONS = {"holo_complete", "holo_monomer", "apo_complete", "apo_monomer"}


def _number(value, name, *, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} is below {minimum}")
    return float(value)


def _close(actual, expected, name):
    actual = _number(actual, name)
    expected = _number(expected, name)
    if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(f"Stale or inconsistent {name}: {actual} != {expected}")


def _cell(value):
    return str(value).replace("|", "\\|").replace("\n", " ")


def _table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |"] +
                     ["| " + " | ".join(_cell(value) for value in row) + " |" for row in rows])


class _Evidence:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.sources = {}

    def read(self, path):
        path = Path(path)
        if not path.is_absolute():
            path = self.root / path
        path = path.resolve()
        if not path.is_relative_to(self.root):
            raise ValueError("Report evidence must stay within the repository")
        content = path.read_bytes()
        self.sources[path.relative_to(self.root).as_posix()] = hashlib.sha256(content).hexdigest()
        return content

    def json(self, path):
        return json.loads(self.read(path))

    def csv(self, path):
        return list(csv.DictReader(self.read(path).decode().splitlines()))


def _configuration(evidence, directory):
    config_path = directory / "configuration.json"
    config = evidence.json(config_path)
    seeds = config.get("seeds")
    if not isinstance(seeds, list) or len(seeds) != 3 or len(set(seeds)) != 3:
        raise ValueError(f"Exactly three unique configured seeds required: {config_path}")
    if any(isinstance(seed, bool) or not isinstance(seed, int) for seed in seeds):
        raise ValueError("Configured seeds must be integers")
    if config.get("information_regime") != "true-site-local":
        raise ValueError("This development report supports only explicitly oracle-site runs")
    if config.get("cpu") != 2:
        raise ValueError("This frozen development report requires two requested CPUs per run")
    return config, evidence.sources[config_path.relative_to(evidence.root).as_posix()]


def _run_records(evidence, directory, conditions):
    config, config_hash = _configuration(evidence, directory)
    runs = evidence.json(directory / "runs.json")
    if not isinstance(runs, list) or not runs:
        raise ValueError(f"Empty run manifest: {directory}")
    identities = [(run.get("receptor"), run.get("seed")) for run in runs]
    expected = {(condition, seed) for condition in conditions for seed in config["seeds"]}
    if len(identities) != len(set(identities)) or set(identities) != expected:
        raise ValueError(f"Incomplete or duplicate configured three-seed denominator: {directory}")
    for run in runs:
        if run.get("status") not in {"ok", "failed", "unsupported"}:
            raise ValueError(f"Unfinished/unknown run status: {directory}")
        if run.get("fingerprint", {}).get("config_sha256") != config_hash:
            raise ValueError(f"Changed configuration fingerprint: {directory}")
        for resource in ("wall_seconds", "cpu_seconds"):
            _number(run.get(resource), resource, minimum=0)
        if run["status"] == "ok":
            pose_path = directory / run["pose_file"]
            content = evidence.read(pose_path)
            if hashlib.sha256(content).hexdigest() != run.get("pose_sha256"):
                raise ValueError(f"Changed pose artifact: {pose_path}")
            lines = content.decode().splitlines()
            if sum(line.startswith("MODEL") for line in lines) != run.get("returned_poses"):
                raise ValueError(f"Returned-pose count mismatch: {pose_path}")
    return config, runs


def _structural_case(evidence, summary_path):
    directory = summary_path.parent
    summaries = evidence.json(summary_path)
    if not isinstance(summaries, list) or not summaries:
        raise ValueError(f"Empty structural summary: {summary_path}")
    case_ids = {row.get("case_id") for row in summaries}
    if len(case_ids) != 1 or None in case_ids:
        raise ValueError("Each structural directory must contain one case")
    case_id = next(iter(case_ids))
    conditions = _AAV_CONDITIONS if case_id.upper() == "AAV9" else {"holo_complete"}
    config, runs = _run_records(evidence, directory, conditions)
    run_map = {(run["receptor"], run["seed"]): run for run in runs}
    identities = [(row.get("receptor"), row.get("seed")) for row in summaries]
    if len(set(identities)) != len(identities) or set(identities) != set(run_map):
        raise ValueError(f"Incomplete or duplicate three-seed structural summary: {summary_path}")
    metrics = evidence.csv(directory / "pose_metrics.csv")
    poses = defaultdict(list)
    for row in metrics:
        if row.get("case_id") != case_id:
            raise ValueError("Pose metrics mix case identities")
        key = (row["receptor"], int(row["seed"]))
        if key not in run_map:
            raise ValueError("Pose metrics contain an unconfigured run")
        row = {**row, "rank": int(row["rank"])}
        for field in ("rmsd_A", "vina_score", "contact_recall", "contact_precision",
                      "interface_residue_recall", "interface_residue_precision"):
            row[field] = _number(float(row[field]), field)
        poses[key].append(row)
    native_contact_counts, native_residue_counts = set(), set()
    for row in summaries:
        key = (row["receptor"], row["seed"])
        raw = run_map[key]
        for field in ("status", "case_id", "pose_sha256"):
            if row.get(field) != raw.get(field):
                raise ValueError(f"Summary/run mismatch in {field}: {summary_path}")
        for field in ("wall_seconds", "cpu_seconds"):
            _close(row.get(field), raw.get(field), field)
        if row["status"] != "ok":
            if poses[key]:
                raise ValueError("Failed run has successful pose metrics")
            continue
        candidates = sorted(poses[key], key=lambda pose: pose["rank"])
        if [pose["rank"] for pose in candidates] != list(range(1, raw["returned_poses"] + 1)):
            raise ValueError(f"Missing or duplicate retained poses: {summary_path}/{key}")
        top1, top5, best = candidates[0], candidates[:5], min(candidates, key=lambda pose: pose["rmsd_A"])
        for field, actual in (("top1_rmsd_A", top1["rmsd_A"]),
                              ("top5_min_rmsd_A", min(p["rmsd_A"] for p in top5)),
                              ("best_returned_rmsd_A", best["rmsd_A"]),
                              ("top1_score", top1["vina_score"]),
                              ("top1_contact_recall", top1["contact_recall"]),
                              ("top1_contact_precision", top1["contact_precision"])):
            _close(row.get(field), actual, field)
        for threshold in (2, 3):
            for prefix, distance in (("top1", top1["rmsd_A"]),
                                     ("top5", min(p["rmsd_A"] for p in top5)),
                                     ("retained", best["rmsd_A"])):
                field = f"{prefix}_success_{threshold}A"
                if row.get(field) is not (distance <= threshold):
                    raise ValueError(f"Stale or missing threshold flag: {field}")
        row["_top1_pose"] = top1
        for pose in candidates:
            native_contact_counts.add(int(pose["native_contact_count"]))
            native_residue_counts.add(int(pose["native_interface_residue_count"]))
    if len(native_contact_counts) > 1 or len(native_residue_counts) > 1:
        raise ValueError(f"Native contact denominator changed between conditions: {summary_path}")
    groups = []
    for condition in sorted(conditions):
        selected = [row for row in summaries if row["receptor"] == condition]
        completed = [row for row in selected if row["status"] == "ok"]
        groups.append({"case_id": case_id, "condition": condition, "rows": selected,
                       "completed": completed, "eligible": len(selected),
                       "native_contacts": next(iter(native_contact_counts), None),
                       "native_residues": next(iter(native_residue_counts), None)})
    return case_id, groups


def _resource_rows(case, condition, runs):
    return [case, condition, str(len(runs)), str(sum(run["status"] == "ok" for run in runs)),
            str(sum(run.get("returned_poses", 0) for run in runs if run["status"] == "ok")),
            f"{sum(run['wall_seconds'] for run in runs):.3f}",
            f"{sum(run['cpu_seconds'] for run in runs):.3f}",
            str(sum(len(run.get("previous_attempts", [])) for run in runs))]


def _panel(evidence, results_dir):
    summary = evidence.json(results_dir / "gal3_panel_summary.json")
    panel = evidence.json("benchmarks/galectin3_itc_panel.json")
    config = evidence.json("benchmarks/gal3_run.json")
    expected_cases = set(config["cases"])
    entries = summary.get("ligands", [])
    case_ids = [row.get("case_id") for row in entries]
    if len(set(case_ids)) != len(case_ids) or set(case_ids) != expected_cases:
        raise ValueError("Incomplete or duplicate Gal3 panel condition denominator")
    measured = {row["ligand_id"]: row for row in panel["rows"]}
    if len(panel["rows"]) != 3 or len(measured) != 3:
        raise ValueError("This exploratory report requires exactly three measured ligand identities")
    if len(config["cases"]) != 6 or expected_cases != {
            f"gal3_{ligand}_{anomer}" for ligand in measured for anomer in ("alpha", "beta")}:
        raise ValueError("Gal3 frozen panel must contain both anomers of all three ligands")
    output, resource_rows, recomputed = [], [], []
    for entry in sorted(entries, key=lambda row: row["case_id"]):
        case_id = entry["case_id"]
        if entry.get("anomer") not in {"alpha", "beta"} or case_id != f"gal3_{entry['ligand_id']}_{entry['anomer']}":
            raise ValueError("Gal3 case identity disagrees with ligand/anomer identity")
        _, runs = _run_records(evidence, results_dir / case_id, {"holo_complete"})
        if any(run["status"] != "ok" for run in runs):
            raise ValueError("Quantitative panel ranking requires all frozen runs completed")
        if any(run.get("case_id") != case_id for run in runs):
            raise ValueError("Gal3 run manifest has the wrong ligand/anomer case identity")
        scores = []
        for run in runs:
            lines = evidence.read(results_dir / case_id / run["pose_file"]).decode().splitlines()
            first = next((line for line in lines if line.startswith("REMARK VINA RESULT:")), None)
            if first is None:
                raise ValueError("Missing Top1 Vina score in panel pose")
            scores.append(_number(float(first.split()[3]), "Top1 Vina score"))
        if entry.get("eligible_runs") != 3 or entry.get("completed_runs") != 3:
            raise ValueError("Incomplete Gal3 panel summary seed counts")
        if len(entry.get("top1_scores_kcal_mol", [])) != 3:
            raise ValueError("Missing Gal3 panel seed scores")
        for raw, stored in zip(sorted(scores), sorted(entry["top1_scores_kcal_mol"])):
            _close(raw, stored, "panel seed score")
        _close(entry.get("seed_mean_top1_score"), mean(scores), "panel seed-mean score")
        source = measured[entry["ligand_id"]]
        for field in ("kd_uM", "kd_sem_uM"):
            _close(entry.get(field), source[field], field)
        if entry.get("heavy_atoms") != source["expected_heavy_atoms"]:
            raise ValueError("Panel heavy-atom baseline changed")
        output.append(entry)
        resource_rows.append(_resource_rows(case_id, "holo_complete", runs))
    correlations = {row["anomer"]: row for row in summary["correlations"]}
    if len(summary["correlations"]) != 2 or set(correlations) != {"alpha", "beta"}:
        raise ValueError("Missing or duplicate Gal3 anomer correlations")
    for anomer in ("alpha", "beta"):
        selected = [row for row in output if row["anomer"] == anomer]
        if {row["ligand_id"] for row in selected} != set(measured) or len(selected) != len(measured):
            raise ValueError("Incomplete Gal3 anomer ligand panel")
        for method in ("vina", "heavy_atom_baseline"):
            rows = [{"context_id": panel["panel_id"], "target_id": panel["target"]["uniprot_id"],
                     "ligand_id": row["ligand_id"], "kd": row["kd_uM"], "kd_unit": "uM",
                     "score": -row["seed_mean_top1_score"] if method == "vina" else row["heavy_atoms"]}
                    for row in selected]
            correlation = quantitative_specificity_by_context(rows)["contexts"][0]
            archived = correlations[anomer][method]["contexts"]
            if len(archived) != 1 or archived[0]["n_ligands"] != len(selected):
                raise ValueError("Stale Gal3 correlation denominator")
            for field in ("spearman_rho", "kendall_tau_b"):
                if correlation[field] is None or archived[0][field] is None:
                    if correlation[field] != archived[0][field]:
                        raise ValueError("Stale undefined Gal3 correlation")
                else:
                    _close(archived[0][field], correlation[field], field)
            recomputed.append([anomer, method, str(len(selected)),
                               "undefined" if correlation["spearman_rho"] is None else f"{correlation['spearman_rho']:.3f}",
                               "undefined" if correlation["kendall_tau_b"] is None else f"{correlation['kendall_tau_b']:.3f}"])
    return output, recomputed, resource_rows, summary


def _native_control(evidence, structural):
    directory = evidence.root / "results/development/aav9/native_control"
    if not (directory / "summary.json").exists():
        return "No archived native-coordinate scoring control is supplied."
    rows = evidence.json(directory / "summary.json")
    preparation = evidence.json(directory / "ligand_preparation.json")
    ligand = evidence.read(directory / "native_ligand.pdbqt")
    evidence.read(directory / "native_ligand.sdf")
    ligand_hash = hashlib.sha256(ligand).hexdigest()
    if preparation.get("native_ligand_sha256") != ligand_hash:
        raise ValueError("Changed native-control ligand")
    if (preparation.get("diagnostic_only") is not True or preparation.get("count_as_docking_success") is not False
            or preparation.get("native_heavy_atoms_optimized") is not False
            or preparation.get("hydrogen_optimization") is not False):
        raise ValueError("Native-coordinate score control must remain a fixed-coordinate diagnostic")
    conditions = [row.get("receptor") for row in rows]
    if len(conditions) != 4 or set(conditions) != _AAV_CONDITIONS:
        raise ValueError("Incomplete or duplicate native-control model conditions")
    groups = {group["condition"]: group for group in structural if group["case_id"].upper() == "AAV9"}
    output = []
    for row in sorted(rows, key=lambda row: row["receptor"]):
        if (row.get("status") != "ok" or row.get("diagnostic_only") is not True
                or row.get("count_as_docking_success") is not False):
            raise ValueError("Native control cannot be counted as successful docking")
        fingerprint = row.get("fingerprint", {})
        if fingerprint.get("native_ligand_sha256") != ligand_hash:
            raise ValueError("Native-control ligand fingerprint mismatch")
        rmsd = _number(row.get("native_heavy_rmsd_A"), "native heavy-coordinate RMSD", minimum=0)
        if rmsd > 0.001:
            raise ValueError("Native heavy coordinates moved beyond rounding tolerance")
        _close(rmsd, preparation["native_heavy_rmsd_A"], "native heavy-coordinate RMSD")
        log = evidence.read(directory / row["log_file"])
        if hashlib.sha256(log).hexdigest() != row.get("log_sha256"):
            raise ValueError("Native score-only log checksum mismatch")
        values = re.findall(r"Estimated Free Energy of Binding\s*:\s*([-+]?\d+(?:\.\d+)?)\s*\(kcal/mol\)", log.decode())
        if len(values) != 1:
            raise ValueError("Native score-only log has no unique energy")
        score = _number(row.get("native_score_only_kcal_mol"), "native score")
        _close(score, float(values[0]), "native log energy")
        docking = {run["seed"]: run for run in groups[row["receptor"]]["completed"]}
        comparisons = row.get("docked_top1_comparisons", [])
        if len(comparisons) != len(docking) or {item["seed"] for item in comparisons} != set(docking):
            raise ValueError("Native-control comparison has incomplete docking seeds")
        for comparison in comparisons:
            run = docking[comparison["seed"]]
            if comparison.get("pose_sha256") != run.get("pose_sha256") or comparison.get("pose_file") != run.get("pose_file"):
                raise ValueError("Native-control comparator pose changed")
            for name in ("config_sha256", "executable_sha256", "receptor_sha256", "preparation_sha256", "reference_sha256"):
                value = fingerprint.get(name)
                if (not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value)
                        or value != run.get("fingerprint", {}).get(name)):
                    raise ValueError(f"Native-control input mismatch: {name}")
            _close(comparison.get("docked_top1_rmsd_A"), run["top1_rmsd_A"], "native comparator RMSD")
            _close(comparison.get("docked_top1_score_kcal_mol"), run["top1_score"], "native comparator score")
            _close(comparison.get("native_minus_docked_top1_score_kcal_mol"), score - run["top1_score"], "native-minus-docked energy")
        wall = _number(row.get("wall_seconds"), "native score wall time", minimum=0)
        average_docked = mean(run["top1_score"] for run in docking.values())
        output.append([row["receptor"], f"{rmsd:.6f}", f"{score:.3f}", f"{average_docked:.3f}",
                       f"{score - average_docked:.3f}", f"{wall:.3f}"])
    return "\n\n".join([
        "Experimental native heavy coordinates were supplied directly; RDKit added modeled hydrogens without optimization. These four score-only evaluations are diagnostics, **never docking successes**, and are excluded from the docking run/pose totals. Positive native-minus-docked gaps mean the prepared Vina score favors the returned Top1 pose over the supplied reference. This is a score-consistency observation for these inputs, not evidence about internal search sampling or measured affinity.",
        _table(["AAV9 condition", "Native heavy RMSD Å", "Native score kcal/mol", "Mean docked Top1 kcal/mol", "Native minus docked kcal/mol", "Score-only wall s"], output),
    ])


def generate_report(root):
    """Return report Markdown; validate every included completed seed and artifact."""
    evidence = _Evidence(root)
    results_dir = evidence.root / "results/development"
    preparation = evidence.json("benchmarks/development_preparation.json")
    frozen = evidence.json("benchmarks/external_run.json")
    evidence.json("benchmarks/development.json")
    cases = preparation["cases"]
    case_ids = [row["case_id"] for row in cases]
    if (len(case_ids) != 7 or len(set(case_ids)) != 7 or len(frozen["cases"]) != 7
            or set(case_ids) != set(frozen["cases"])):
        raise ValueError("The seven-case frozen external preparation denominator changed")
    if preparation.get("eligible_cases") != len(cases) or any(row.get("eligible_preparation_denominator") is not True for row in cases):
        raise ValueError("Every frozen external case must remain in the eligible denominator")
    statuses = Counter(row["status"] for row in cases)
    if set(statuses) - {"ok", "failed", "unsupported", "pending", "running"}:
        raise ValueError("Unknown preparation status")
    for field, status in (("successful_preparations", "ok"), ("failed_preparations", "failed"),
                          ("unsupported_preparations", "unsupported")):
        if preparation.get(field) != statuses[status]:
            raise ValueError(f"Stale preparation count: {field}")
    structural, seen_cases = [], set()
    for path in sorted(results_dir.glob("*/summary.json")):
        case_id, groups = _structural_case(evidence, path)
        normalized = case_id.upper()
        if normalized in seen_cases or normalized not in set(case_ids) | {"AAV9"}:
            raise ValueError("Duplicate or unknown structural case summary")
        seen_cases.add(normalized)
        structural.extend(groups)
    if "AAV9" not in seen_cases:
        raise ValueError("AAV9 diagnostic summary is required")
    case_by_id = {row["case_id"]: row for row in cases}
    for case_id in seen_cases - {"AAV9"}:
        if case_by_id[case_id]["status"] != "ok":
            raise ValueError("Docking summary exists for a preparation not declared successful")
    for case in cases:
        explicitly_complete = case.get("docking_complete") is True or case.get("docking_status") in {"ok", "complete", "completed"}
        if explicitly_complete and case["case_id"] not in seen_cases:
            raise ValueError(f"Listed completed docking case has no summary: {case['case_id']}")
    panel_rows, correlations, panel_resources, panel_summary = _panel(evidence, results_dir)
    native_control = _native_control(evidence, structural)
    for source in ("aav_benchmark/report.py", "aav_benchmark/statistics.py", "scripts/regenerate_report.py"):
        if (evidence.root / source).exists():
            evidence.read(source)
    pose_rows, threshold_rows, contact_rows, resource_rows = [], [], [], []
    for group in sorted(structural, key=lambda group: (group["case_id"], group["condition"])):
        case, condition = group["case_id"], group["condition"]
        completed, eligible = group["completed"], group["eligible"]
        def avg(field):
            return f"{mean(row[field] for row in completed):.3f}" if completed else "undefined"
        pose_rows.append([case, condition, f"{len(completed)}/{eligible}", avg("top1_rmsd_A"),
                          avg("top5_min_rmsd_A"), avg("best_returned_rmsd_A")])
        threshold_rows.append([case, condition] +
                              [f"{sum(row[f'{prefix}_success_{threshold}A'] for row in completed)}/{eligible}"
                               for threshold in (2, 3) for prefix in ("top1", "top5", "retained")])
        top1 = [row["_top1_pose"] for row in completed]
        def contact_avg(field):
            return f"{mean(row[field] for row in top1):.3f}" if top1 else "undefined"
        contact_rows.append([case, condition, str(group["native_contacts"]), str(group["native_residues"]),
                             contact_avg("contact_recall"), contact_avg("contact_precision"),
                             contact_avg("interface_residue_recall"), contact_avg("interface_residue_precision"),
                             f"{sum(int(row['evaluation_pairs_below_1_5A']) > 0 for row in top1)}/{eligible}",
                             ", ".join(sorted({row["evaluation_receptor_scope"] for row in top1}))])
        resource_rows.append(_resource_rows(case, condition, group["rows"]))
    preparation_rows = []
    preparation_notes = []
    for case in sorted(cases, key=lambda row: row["case_id"]):
        docking = "summary present" if case["case_id"] in seen_cases else "not completed / no summary"
        preparation_rows.append([case["case_id"], case["status"], docking,
                                 case.get("failure_category", "—"), case.get("reason", "—")])
        notes = [case[field] for field in ("reference_scope", "calcium_model", "protein_gap") if case.get(field)]
        if notes:
            preparation_notes.append(f"- **{case['case_id']}**: " + "; ".join(notes) + ".")
    all_resources = resource_rows + panel_resources
    total_runs = sum(int(row[2]) for row in all_resources)
    total_completed = sum(int(row[3]) for row in all_resources)
    total_poses = sum(int(row[4]) for row in all_resources)
    # Preserve full precision in total resource sums, rather than summing rounded table cells.
    total_wall = total_cpu = 0.0
    for path in sorted(results_dir.glob("*/runs.json")):
        runs = evidence.json(path)
        if path.parent.name not in {"gal3_" + ligand + "_" + anomer for ligand in
                                        ("LacNAc", "lactose", "galactose") for anomer in ("alpha", "beta")} and not any(
                                            group["case_id"].lower() == path.parent.name.lower() for group in structural):
            continue
        total_wall += sum(run["wall_seconds"] for run in runs)
        total_cpu += sum(run["cpu_seconds"] for run in runs)
    source_rows = sorted(evidence.sources.items())
    digest = hashlib.sha256(json.dumps(source_rows, separators=(",", ":")).encode()).hexdigest()
    text = [
        "# Development docking results", "",
        "Generated deterministically by `python scripts/regenerate_report.py` from the archived summaries, run manifests, configuration files, per-pose metrics and pose artifacts listed below. `--check` verifies that this document is current. No raw result is modified.", "",
        "**The repository's original scientific purpose remains incomplete.** These are AutoDock Vina 1.2.7 development diagnostics using true-site information. They do not establish blind site recovery, held-out pose reliability, discrimination of tested nonbinders, or reliability of other docking methods. No confirmatory cohort has been evaluated here.", "",
        "## Eligible cases and preparation", "",
        f"The frozen external development selection retains **{len(cases)} eligible cases**: {statuses['ok']} successful preparations, {statuses['failed']} failed and {statuses['unsupported']} unsupported; {statuses['pending'] + statuses['running']} still pending/running. These are preparation outcomes, not molecular binding outcomes. Docking accuracy below is conditional on available completed preparation and execution; absent docking is explicit.", "",
        _table(["Case", "Preparation", "Docking evidence", "Failure category", "Reason"], preparation_rows), "",
        "\n".join(preparation_notes), "",
        "AAV9 is a separate diagnostic. Its apo/holo and complete-neighborhood/monomer conditions are paired model conditions for one biological system, not four independent cases. Monomer conditions are deliberately incomplete controls. External released receptors are not certified complete biological assemblies. Source reproduction alone does not qualify local density, stereochemistry, assemblies or family independence for confirmatory inference.", "",
        "External protein-group labels are provisional descriptive annotations. In particular, 1PWB and 2VUZ share the broad C-type lectin fold; no count of independently qualified families or across-family generalization follows from this selection.", "",
        "## Structural pose diagnostics", "",
        "Each completed condition must contain its three configured seeds. Values are means across completed seeds within a case and condition; Top5 is the lowest RMSD among the five highest-ranked poses, and best retained is the lowest RMSD across all returned poses. Ligands are evaluated in the receptor coordinate frame without independent ligand alignment. No across-family confidence interval is estimated from these unqualified development cases.", "",
        _table(["Case", "Condition", "Completed/eligible seeds", "Mean Top1 RMSD Å", "Mean Top5 minimum RMSD Å", "Mean best retained RMSD Å"], pose_rows), "",
        "Counts below retain all configured seeds in each denominator. Repeated seeds and poses do not supply independent biological replication. A near-native pose absent from retained outputs does not identify internal sampling failure.", "",
        _table(["Case", "Condition", "Top1 ≤2 Å", "Top5 ≤2 Å", "Retained ≤2 Å", "Top1 ≤3 Å", "Top5 ≤3 Å", "Retained ≤3 Å"], threshold_rows), "",
        "## Contacts and physical overlap", "",
        "Native heavy-atom contact identities use the fixed holo reference and the <4 Å cutoff. The generator checks that native contact and interface-group denominators remain constant across seeds and conditions for each case. Contact values are seed means for Top1. Interface groups follow receptor PDB residue records and include retained cofactors such as CA2; these values do not measure protein-residue-only coverage. Severe overlap means at least one evaluated ligand–receptor atom pair below 1.5 Å; this distance screen does not validate ring geometry, linkages or chemical plausibility.", "",
        _table(["Case", "Condition", "Native atom contacts", "Native interface groups", "Contact recall", "Contact precision", "Interface-group recall", "Interface-group precision", "Top1 overlap seeds", "Overlap evaluation scope"], contact_rows), "",
        "## Native-coordinate scoring diagnostic", "",
        native_control, "",
        "## Galectin-3 quantitative assay panel", "",
        "This exploratory panel contains three measured binders from one ITC context. All three compounds bind; there are no tested nonbinders, so binary AUROC/AUPRC, sensitivity and specificity are ineligible. KD values are published fitted estimates (±SEM), not raw ITC data; their uncertainty is not propagated into these descriptive rank correlations. The docking receptor is the 1KJL carbohydrate-recognition domain whereas the assay used full-length galectin-3. Reducing-end alpha and beta models are separate sensitivity conditions; solution ensemble weights are unknown and the six modeled inputs are not six independent assay ligands.", "",
        _table(["Ligand", "Heavy atoms", "Measured KD ±SEM µM", "Alpha mean Top1 energy kcal/mol", "Beta mean Top1 energy kcal/mol"],
               [[ligand, str(next(row for row in panel_rows if row["ligand_id"] == ligand)["heavy_atoms"]),
                 f"{next(row for row in panel_rows if row['ligand_id'] == ligand)['kd_uM']:.1f} ± {next(row for row in panel_rows if row['ligand_id'] == ligand)['kd_sem_uM']:.1f}",
                 f"{next(row for row in panel_rows if row['ligand_id'] == ligand and row['anomer'] == 'alpha')['seed_mean_top1_score']:.3f}",
                 f"{next(row for row in panel_rows if row['ligand_id'] == ligand and row['anomer'] == 'beta')['seed_mean_top1_score']:.3f}"]
                for ligand in sorted({row["ligand_id"] for row in panel_rows})]), "",
        "Correlations compare higher-is-stronger scores (negative Vina energy, or heavy-atom count) with measured −log10(KD in M), separately by anomer. They describe three point estimates; no significance or generalization claim follows.", "",
        _table(["Anomer", "Score", "Ligands", "Spearman ρ", "Kendall τ-b"], correlations), "",
        "Size alone ranks this panel perfectly, so the raw Vina rank does not demonstrate recognition beyond size. All compounds are neutral; a charge-only score is constant. Exact size/charge matching is unavailable in this three-ligand panel. No dissociation constant is inferred from docking energy.", "",
        "## Execution resources", "",
        "Wall and CPU values are sums of recorded final-attempt run times; wall sums are not elapsed study duration when jobs overlap. Prior-attempt counts are shown separately; their time is not included here. Each frozen run used two requested CPUs. Returned poses are retained candidates, not an internal sampling count.", "",
        _table(["Case", "Condition", "Eligible runs", "Completed runs", "Returned poses", "Wall sum s", "CPU sum s", "Prior attempts"], all_resources), "",
        f"Included totals: **{total_completed}/{total_runs} completed runs**, **{total_poses} returned poses**, **{total_wall:.3f} wall-seconds summed** and **{total_cpu:.3f} CPU-seconds summed**. Preparation failures remain in the seven-case table rather than disappearing from scope.", "",
        "## Remaining scientific closure", "",
        "The next evidence gates remain independently qualified structural references and assemblies, frozen family splits and training-overlap audits, blind/predicted-site methods and comparable information budgets, qualified tested alternatives, held-out execution and family-level uncertainty. The present results support implementation diagnostics and bounded failure observations only. They do not complete the original scientific question.", "",
        "## Source provenance", "",
        f"SHA256 of the sorted path/hash manifest below: `{digest}`. Paths are repository-relative. All listed artifacts were read or verified during generation; no timestamp or machine-specific path is inserted into this report.", "",
        _table(["Source", "SHA256"], [[f"`{path}`", f"`{checksum}`"] for path, checksum in source_rows]), "",
    ]
    return "\n".join(text)


def write_report(root, output=None, *, check=False):
    """Write the report, or fail when an existing report differs in check mode."""
    root = Path(root).resolve()
    output = Path(output).resolve() if output is not None else root / "docs/DEVELOPMENT_RESULTS.md"
    result = generate_report(root)
    if check:
        if not output.exists() or output.read_text() != result:
            raise ValueError(f"Generated report is stale or missing: {output}")
    else:
        write_text_atomic(output, result)
    return output
