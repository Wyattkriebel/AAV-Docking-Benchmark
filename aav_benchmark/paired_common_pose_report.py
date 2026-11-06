"""H2 exact-union ranking, reusing frozen main metrics after workflow release."""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile

from . import paired_method_report as main
from . import paired_common_pose_rescoring as union
from .common_pose_rescoring import canonical_sha256, model_blocks, _old_records, _record_path
from .additional_method import expected_chi_links, verify_chi_recognition
from .provenance import sha256, write_text_atomic
from .qualified_pose import Evidence, accuracy_summary
from .statistics import aggregate_runs, paired_family_bootstrap

CASES, PROFILES, ARMS, SEEDS = union.CASES, union.PROFILES, union.METHODS, union.SEEDS
RANK_METHODS = ("VinaCarb1.0_Vina1.1.2_core_shared_zero_base_coeff0",
                "VinaCarb1.0_Vina1.1.2_core_shared_zero_base_plus_CHI_coeff1")
DEFAULT_PLAN = "benchmarks/structural_paired_common_pose_evaluation.proposal.json"
MAIN_PLAN = {"path": "benchmarks/structural_paired_method_evaluation_v2.json", "sha256": "4922b83f6e28b221900bac1ee0b52fd6e8056f97e245c277bbd71ca28f4e807e"}
MAIN_FREEZE = {"path": "results/structural_pose_v3/paired_evaluation_v2_freeze.json", "sha256": "19f45adee79d3dbedcd8b788985e0f0770a55a212bbe51f9deb2a8f38edb15c0"}
UNION_PLAN = {"path": "benchmarks/structural_pose_v3_common_pool.json", "sha256": "11d658e719292e72c8dbb7531cf361411b0134162d8363bfce4605fd9b28fd0e"}
UNION_FREEZE = {"path": "results/structural_pose_v3/common_pose_union/freeze_receipt.json", "sha256": "66ca0a298d01ed1058ae4dd11d7d0f23d5183759217986543a7d1a70a8ada5a8"}
CODE_SOURCES = tuple(sorted(set(main.CODE_SOURCES) | union.EXECUTION_HELPERS | {
    "aav_benchmark/paired_common_pose_rescoring.py", "aav_benchmark/paired_common_pose_report.py"}))
PROTOCOL = {"RMSD_thresholds_A": [1, 2, 3], "primary_RMSD_A": 2, "alignments": ["global", "local_interface"],
            "pose_rank_budgets": [1, 5], "best_retained": True, "independent_ligand_fit": False,
            "full_and_core_separate": True, "geometry_filters_success": False}
STATISTICS = {"n_resamples": 10000, "seed": 20261103, "confidence": .95, "unit": "source_family",
              "difference": "coefficient1_minus_coefficient0", "failed_eligible_outcomes_zero": True,
              "conditional_completed_separate": True}
TIMING = {"fixed_score_records_already_existed": 580, "atomic_accuracy_consulted": False,
          "metric_choices_inherited_from_prior_freezes": True}
METRIC_FIELDS = ("geometry_diagnostics", "full_parent_min_protein_distance_A", "full_parent_protein_pairs_below_1p5A",
    "full_observed_assembly_min_protein_distance_A", "full_observed_assembly_protein_pairs_below_1p5A",
    *(f"{a}_{f}" for a in ("global", "local_interface") for f in ("RMSD_A", "qualified_copy_id", "chemical_atom_map_index", "contacts")))


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _pin(evidence, record):
    return evidence.read(record["path"], record["sha256"])


def validate_report_plan(root, plan_path, freeze_path=None, *, require_freeze=True):
    """Pin inherited choices and code; read no experimental atom coordinates."""
    evidence = Evidence(root); plan = evidence.json(plan_path)
    _require(plan["case_ids"] == CASES and plan["profiles"] == PROFILES and plan["arms"] == ARMS
             and plan["seeds"] == SEEDS and plan["ranking_methods"] == list(RANK_METHODS)
             and plan["eligible_seed_pools"] == 108 and plan["eligible_source_method_slots"] == 216
             and plan["eligible_rescore_method_slots"] == 216 and plan["deduplicate_or_filter_poses"] is False,
             "H2 fixed cohort, exact union or ranking methods changed")
    _require(all(plan[key] == pin for key, pin in (("main_evaluation_plan", MAIN_PLAN), ("main_evaluation_freeze", MAIN_FREEZE),
             ("union_scoring_plan", UNION_PLAN), ("union_scoring_freeze", UNION_FREEZE))), "H2 inherited frozen plan identities changed")
    _require(plan["metric_protocol"] == PROTOCOL and plan["statistics"] == STATISTICS
             and plan["implementation_timing"] == TIMING and plan["global_four_workflow_guard_required"] is True,
             "Inherited metric/statistical choices or implementation timing changed")
    _require(set(plan["metric_code_sha256"]) == set(CODE_SOURCES), "Every H2 runtime helper needs an exact source pin")
    for path, digest in plan["metric_code_sha256"].items():
        evidence.read(path, digest)
    _require(plan["adapter_source"] == {"path": "aav_benchmark/paired_common_pose_report.py", "sha256": sha256(Path(__file__))},
             "Running H2 adapter differs from its frozen source")
    main_plan, main_freeze = json.loads(_pin(evidence, MAIN_PLAN)), json.loads(_pin(evidence, MAIN_FREEZE))
    _require(main_plan.get("common_pool_rescoring_enabled") is False and main_freeze["evaluation_plan_sha256"] == MAIN_PLAN["sha256"]
             and main_freeze.get("frozen_before_test_pose_accuracy") is True and main_freeze.get("new_test_pose_accuracy_consulted") is False,
             "Main comparison's independent pre-accuracy freeze changed")
    _require(all(plan["metric_code_sha256"].get(p) == h for p, h in main_plan["metric_code_sha256"].items()),
             "Main and H2 metric implementations differ")
    reference = json.loads(_pin(evidence, main_plan["reference_evaluation_plan"]))
    _require(set(reference["case_to_family"]) == set(CASES), "Frozen source-family denominator changed")
    for pin in (main_plan["reference_evaluation_freeze"], main_plan["public_source_distribution"], UNION_PLAN, UNION_FREEZE):
        _pin(evidence, pin)
    _require(plan["union_output_root"] == union.DEFAULT_OUTPUT, "Fixed-score output identity changed")
    if require_freeze:
        _require(freeze_path is not None and plan.get("frozen_before_union_pose_accuracy") is True,
                 "H2 atomic report requires a separate pre-accuracy freeze")
        freeze = evidence.json(freeze_path)
        _require(freeze.get("evaluation_plan_sha256") == sha256(evidence.root / plan_path)
                 and freeze.get("frozen_before_test_pose_accuracy") is True and freeze.get("frozen_before_union_pose_accuracy") is True
                 and freeze.get("new_test_pose_accuracy_consulted") is False, "H2 atomic report requires a separate pre-accuracy freeze")
    return evidence, plan, main_plan, reference


def _read_workflow(evidence, plan):
    """Frozen source index in TEMP; checked real score reads; no live writes."""
    root = evidence.root
    union_plan, actual = union.validate_union_plan(root, plan_path=UNION_PLAN["path"], receipt_path=UNION_FREEZE["path"], require_freeze=True)
    with tempfile.TemporaryDirectory(prefix="h2-read-only-index-") as temporary:
        index = union.index_union_pool(root, plan_path=actual, output_root=Path(temporary))
    output = root / plan["union_output_root"]; manifest_path = output / "scores_manifest.json"
    before = manifest_path.read_bytes() if manifest_path.exists() else None
    records, manifest = _old_records(output, index, recover_orphans=False)
    _require(before == (manifest_path.read_bytes() if manifest_path.exists() else None), "Score manifest changed during validation")
    if before is not None:
        evidence.read(manifest_path.relative_to(root), hashlib.sha256(before).hexdigest())
    for identity, digest in manifest["records"].items():
        evidence.read(_record_path(output, identity).relative_to(root), digest)
    rows = {r["condition_id"]: r for r in json.loads(_pin(evidence, union_plan["source_inputs"]))["conditions"]}
    reconstructed = 0
    for source in index["poses"]:
        record = records.get(source["pose_id"])
        if record is None or "exported_ligand_sha256" not in record:
            continue
        with tempfile.TemporaryDirectory(prefix="h2-read-only-export-") as temporary:
            destination = Path(temporary) / "ligand.pdbqt"
            export = union.export_source_pose(source, rows[source["condition_id"]], root, destination)
            _require(sha256(destination) == record["exported_ligand_sha256"]
                     and sha256(destination.with_suffix(".pdbqt.map.json")) == record["export_map_sha256"]
                     and len(export["expected_glycosidic_links"]) == record["expected_glycosidic_links"],
                     "Fixed score cannot reconstruct its unchanged source pose")
            if record["status"] == "ok":
                log = Path(temporary) / "VC_log.txt"; log.write_text(record["score"]["recognition_log"])
                verify_chi_recognition(log, expected_chi_links(export))
        reconstructed += 1
    return index, records, rows, reconstructed


def _workflow_summary(index, records, reconstructed):
    pools, sources = index["pools"], index["poses"]
    _require(len(pools) == len({p["pool_id"] for p in pools}) == 108
             and len(sources) == len({p["pose_id"] for p in sources}) <= 2280, "H2 pool/pose denominator changed")
    _require(set(records) <= {p["pose_id"] for p in sources} and type(index["terminal_source_method_slots"]) is int,
             "Unexpected fixed-score identities or source count type")
    statuses = {"pending", "ok", "failed", "prepared_input_failure", "unsupported_preflight"}
    _require(all(s in statuses for p in pools for s in p["source_statuses"].values()), "Unknown source terminal status")
    terminal = sum(s != "pending" for p in pools for s in p["source_statuses"].values())
    _require(all(set(p["source_statuses"]) == set(union.SOURCE_METHODS) for p in pools)
             and terminal == index["terminal_source_method_slots"] <= 216, "H2 source status inventory changed")
    done = sum(all(s != "pending" for s in p["source_statuses"].values())
        and all(s["pose_id"] in records for s in sources if s["pool_id"] == p["pool_id"]) for p in pools)
    ready = terminal == 216 and index["all216_source_slots_terminal_and_verified"] is True and len(records) == len(sources) and done == 108
    pending = list(index.get("pending_source_slots", [])) + ["fixed_score/" + s["pose_id"] for s in sources if s["pose_id"] not in records]
    return {"schema_version": 1, "status": "ready_for_atomic_evaluation" if ready else "pending_incomplete_workflow",
        "workflow_ready": ready, "metrics_eligible": False,
        "counts": {"eligible_seed_pools": 108, "eligible_source_method_slots": 216, "eligible_rescore_method_slots": 216,
            "source_method_slots": 216, "rescore_method_slots": 216,
            "terminal_source_method_slots": terminal, "terminal_rescore_method_slots": 2 * done,
            "retained_pose_count": len(sources), "recorded_pose_count": len(records), "source_exports_reconstructed": reconstructed,
            "successful_fixed_score_records": sum(r["status"] == "ok" for r in records.values()),
            "failed_fixed_score_records": sum(r["status"] != "ok" for r in records.values())},
        "pending": sorted(pending), "accuracy_computed": False, "atomic_accuracy_computed": False, "experimental_reference_access": False,
        "native_executables_used": False, "live_output_mutated": False}


def probe_workflow(root, plan_path, freeze_path=None, *, require_freeze=True):
    """Only normalized readiness/counts; no metric, fitting or context call."""
    evidence, plan, _, _ = validate_report_plan(root, plan_path, freeze_path, require_freeze=require_freeze)
    index, records, _, reconstructed = _read_workflow(evidence, plan)
    return _workflow_summary(index, records, reconstructed)


def _native_record(evidence, algorithm, seed, row):
    if algorithm == union.SOURCE_METHODS[1]:
        path = Path(row["vina_carb_outputs"][str(seed)]) / "execution.json"
        record = evidence.json(path)
    elif row["status"] != "prepared":
        path = Path(row["modern_runs"]).parent.parent / "execution.json"
        record = evidence.json(path)
    else:
        path = Path(row["modern_runs"])
        matches = [r for r in evidence.json(path) if r["seed"] == seed]
        _require(len(matches) == 1, "Native source seed record is absent or duplicated")
        record = matches[0]
    return path.as_posix(), record


def join_main_pose_metrics(evidence, main_report, index, rows, main_plan):
    """Bind fresh in-memory main metrics to exact independently verified sources.

    A supplied result is produced in memory by the pinned global guard. No
    cached metric file is accepted by the CLI. These checks bind its inputs and
    source tuples; a rank alone is never sufficient provenance.
    """
    _require(isinstance(main_report, dict), "Main metrics must be a freshly generated in-memory report")
    _require(main_report.get("status") == "complete" and main_report.get("metrics_eligible") is True
             and main_report.get("counts", {}).get("eligible_method_outcomes") == 216
             and main_report["counts"].get("terminal_method_outcomes") == 216
             and main_report.get("metric_code_sha256") == main_plan["metric_code_sha256"],
             "H2 requires a fresh complete result from the exact frozen main reader")
    inventory = main_report["source_sha256_inventory"]
    conditions = main_report["conditions"]; by_condition = {r["condition_id"]: r for r in conditions}
    _require(len(conditions) == len(by_condition) == len(rows) == 36 and set(by_condition) == set(rows),
             "Main/union condition identities are missing, duplicate or extra")
    main_poses, source_poses = {}, {s["pose_id"]: s for s in index["poses"]}
    pools = {(p["condition_id"], p["seed"]): p for p in index["pools"]}
    for condition in conditions:
        cid, row = condition["case_id"], rows[condition["condition_id"]]
        _require(cid == row["case_id"] and condition["receptor_profile"] == row["receptor_profile"]
                 and condition["site_condition"] == row["method"]
                 and condition["atomic_scope"] == ("modeled_core" if cid == "4R9F" else "full_ligand")
                 and condition["exact_common_input_sha256"] == {k: row[k]["sha256"] for k in
                     ("ligand", "receptor", "preparation", "canonical_atom_map", "initializer_sdf") if k in row},
                 "Main/union chemical, receptor, atomic-mask or common-input identity changed")
        _require(len(condition["seeds"]) == 3 and {s["seed"] for s in condition["seeds"]} == set(SEEDS),
                 "Main metric seed identities are missing or duplicated")
        for seed in condition["seeds"]:
            pool = pools[condition["condition_id"], seed["seed"]]
            _require(set(seed["search_methods"]) == set(union.SOURCE_METHODS), "Main source algorithm identities changed")
            for algorithm, result in seed["search_methods"].items():
                _require(result["status"] == pool["source_statuses"][algorithm], "Main/union source failure status changed")
                path, record = _native_record(evidence, algorithm, seed["seed"], row)
                _require(canonical_sha256(record) == pool["source_record_sha256"][algorithm]
                         and inventory.get(path) == evidence.sources[path], "Main/union raw source record provenance differs")
                for pose in result["poses"]:
                    identity = pool["pool_id"] + "/" + algorithm + f"/p{pose['rank']:02d}"
                    _require(identity not in main_poses, "Duplicate main per-pose metric identity")
                    main_poses[identity] = pose
    _require(set(main_poses) == set(source_poses), "Main metrics and exact union have different pose identities")
    joined, raw_cache = {}, {}
    for identity, source in source_poses.items():
        pose = main_poses[identity]
        row = rows[source["condition_id"]]; pool = pools[source["condition_id"], source["seed"]]
        _require(source["source_record_sha256"] == pool["source_record_sha256"][source["source_algorithm"]]
                 and source["pose_identity_sha256"] == canonical_sha256({k: v for k, v in source.items() if k != "pose_identity_sha256"})
                 and all(source[field] == row[key]["sha256"] for field, key in (
                     ("source_ligand_sha256", "ligand"), ("source_receptor_sha256", "receptor"),
                     ("source_preparation_sha256", "preparation"), ("source_canonical_atom_map_sha256", "canonical_atom_map"))),
                 "Source pose/record/input provenance changed")
        _require(pose["rank"] == source["model_rank"] and pose["pooled_rank"] == source["model_rank"]
                 and pose["candidate_class_rank"] == source["candidate_class_rank"]
                 and pose["candidate_original_rank"] == source["candidate_original_rank"] and pose["box_id"] == "box01"
                 and pose["native_method_score"] == {"method": source["source_algorithm"],
                     "value": source["original_method_score"], "cross_method_comparable": False},
                 "Main pose rank, source score or candidate identity differs from the union")
        _require(all(field in pose for field in METRIC_FIELDS), "Main pose lacks an inherited metric/diagnostic field")
        for field in ("global_RMSD_A", "local_interface_RMSD_A"):
            value = pose[field]
            _require(not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value) and value >= 0,
                     "Invalid inherited distance metric")
        path = source["source_pose_file"]
        _require(inventory.get(path, {}).get("sha256") == source["source_pose_file_sha256"],
                 "A per-pose source file digest differs from the main inventory")
        if path not in raw_cache:
            raw = evidence.read(path, source["source_pose_file_sha256"])
            _require(inventory.get(path) == evidence.sources[path], "Main per-pose raw file hash differs from the union")
            raw_cache[path] = model_blocks(raw)
        _require(source["model_rank"] in raw_cache[path]
                 and hashlib.sha256(raw_cache[path][source["model_rank"]]).hexdigest() == source["source_model_sha256"],
                 "Source MODEL block hash changed")
        proof = copy.deepcopy(source)
        joined[identity] = {**copy.deepcopy(pose), "pose_id": identity, "source_model_rank": source["model_rank"],
                            "source_algorithm": source["source_algorithm"], "source_provenance": proof}
    return joined, by_condition


def summarize_ranked_pools(index, records, joined, main_conditions):
    """Reuse frozen ordering without filtering geometry, failures or duplicates."""
    pools, conflicts = union.rank_union_pools(index, records)
    conditions = {}
    for pool in pools:
        cid = pool["condition_id"]; original = main_conditions[cid]
        if cid not in conditions:
            conditions[cid] = {k: copy.deepcopy(original[k]) for k in
                ("condition_id", "case_id", "receptor_profile", "site_condition", "atomic_scope", "scored_atoms", "parent_atoms",
                 "common_preparation_status", "executed_candidate_budget", "native_controls", "site", "exact_common_input_sha256")}
            conditions[cid]["seeds"] = []
        methods = {}
        for coefficient, method in enumerate(RANK_METHODS):
            order = pool[f"coefficient{coefficient}_order"]
            if pool["status"] == "completed":
                _require(len(order) == len(pool["union_pose_ids"]) and set(order) == set(pool["union_pose_ids"]),
                         "Coefficient rankings must contain the entire exact union")
            else:
                _require(not order, "A failed whole union unexpectedly has a partial ranking")
            ranked = []
            for rank, identity in enumerate(order, 1):
                item = copy.deepcopy(joined[identity]); item["pooled_rank"] = rank
                item["ranking_score"] = {"backend": method, "value": records[identity]["score"][
                    "base_shared0" if coefficient == 0 else "combined_shared0"], "cross_ligand_comparable": False}
                ranked.append(item)
            statistics_status = "ok" if pool["status"] == "completed" else "unsupported" if (
                any(s in {"prepared_input_failure", "unsupported_preflight"} for s in pool["source_statuses"].values())
                or any(records[i]["status"] == "unsupported_case" for i in pool["union_pose_ids"])) else "failed"
            methods[method] = {"status": pool["status"], "statistics_status": statistics_status, "poses": ranked,
                "global": accuracy_summary(ranked, "global_RMSD_A"), "local_interface": accuracy_summary(ranked, "local_interface_RMSD_A")}
        if pool["status"] == "completed":
            for alignment in ("global", "local_interface"):
                for field in ("best_returned_RMSD_A", *(f"best_returned_success_{t}A" for t in (1, 2, 3))):
                    _require(methods[RANK_METHODS[0]][alignment][field] == methods[RANK_METHODS[1]][alignment][field],
                             "Best-retained recovery changed despite the exact same union")
        conditions[cid]["seeds"].append({"seed": pool["seed"], "pool_id": pool["pool_id"], "pool_status": pool["status"],
            "source_statuses": pool["source_statuses"], "source_record_sha256": pool["source_record_sha256"],
            "original_source_orders": pool["original_source_orders"],
            "retained_union_pose_ids": pool["union_pose_ids"],
            "fixed_score_failure_pose_ids": [i for i in pool["union_pose_ids"] if records[i]["status"] != "ok"], "ranking_methods": methods})
    return list(conditions.values()), pools, conflicts


def ranking_statistics(conditions, case_to_family):
    output = []
    for profile in PROFILES:
        for arm in ARMS:
            for scope in ("full_ligand", "modeled_core"):
                selected = [r for r in conditions if r["receptor_profile"] == profile and r["site_condition"] == arm and r["atomic_scope"] == scope]
                for alignment in ("global", "local_interface"):
                    for rank in ("top1", "top5", "best_returned"):
                        for threshold in (1, 2, 3):
                            rows = [{"family_id": case_to_family[r["case_id"]], "complex_id": r["case_id"], "method": method,
                                "seed": seed["seed"], "status": value["statistics_status"],
                                "success": value[alignment][f"{rank}_success_{threshold}A"],
                                "near_native_present": value[alignment][f"best_returned_success_{threshold}A"]}
                                for r in selected for seed in r["seeds"] for method, value in seed["ranking_methods"].items()]
                            if not rows:
                                continue
                            summary = paired_family_bootstrap(rows, RANK_METHODS[1], RANK_METHODS[0], n_resamples=10000,
                                seed=20261103, confidence=.95, expected_seeds=SEEDS)
                            summary.update(receptor_profile=profile, site_condition=arm, atomic_scope=scope, alignment=alignment,
                                pose_rank=rank, threshold_A=threshold, eligible_and_conditional_summary=aggregate_runs(rows, expected_seeds=SEEDS),
                                interpretation="Conditional source-family-cluster recovery on this purposive sample; coefficient1-minus0 on one exact union, no population/model-unseen claim")
                            output.append(summary)
    return output


def generate_report(root, plan_path, freeze_path, *, main_report=None, allow_pending=False):
    """Enforce local gates; the calling frozen guard enforces global release.

    Provided main_report is the guard's freshly generated in-memory pinned main
    result. A dictionary is not independent proof of freshness: the trusted
    pinned caller supplies that guarantee. There is no metric-file CLI input.
    """
    evidence, plan, main_plan, reference = validate_report_plan(root, plan_path, freeze_path)
    index, records, rows, reconstructed = _read_workflow(evidence, plan)
    workflow = _workflow_summary(index, records, reconstructed)
    if not workflow["workflow_ready"]:
        _require(allow_pending, "All216 source slots and every fixed score must be terminal before H2 accuracy")
        return workflow
    if main_report is None:
        main_report = main.generate_paired_report(evidence.root, MAIN_PLAN["path"], MAIN_FREEZE["path"])
    joined, main_conditions = join_main_pose_metrics(evidence, main_report, index, rows, main_plan)
    conditions, pools, conflicts = summarize_ranked_pools(index, records, joined, main_conditions)
    inventory = copy.deepcopy(main_report["source_sha256_inventory"])
    for path, pin in evidence.sources.items():
        _require(path not in inventory or inventory[path] == pin, "Main and union evidence inventories disagree")
        inventory[path] = pin
    return {"schema_version": 1, "status": "complete", "metrics_eligible": True, "counts": workflow["counts"],
        "conditions": conditions, "per_source_pose_metrics": [joined[s["pose_id"]] for s in index["poses"]],
        "paired_family_statistics": ranking_statistics(conditions, reference["case_to_family"]),
        "workflow_status_counts": {s: sum(p["status"] == s for p in pools) for s in sorted({p["status"] for p in pools})},
        "normalization_conflict_cases": conflicts,
        "main_evaluation_plan": MAIN_PLAN, "main_evaluation_freeze": MAIN_FREEZE,
        "union_scoring_plan": UNION_PLAN, "union_scoring_freeze": UNION_FREEZE,
        "source_sha256_inventory": inventory, "metric_code_sha256": plan["metric_code_sha256"], "implementation_timing": plan["implementation_timing"],
        "limitations": [
            "Only coefficient0 versus1 on this identical retained union isolates CHI ranking under the original Vina1.1.2 core; this is not modern Vina1.2.7 rescoring.",
            "Every source/model identity remains, including duplicate coordinates. Any source/fixed-pose failure invalidates the whole pool for both coefficients and counts zero recovery.",
            "Comparing union and original-source orders changes pool membership. Absence from the retained union cannot prove internal sampling failure rather than output retention/pruning.",
            "Five full ligands and one11-atom modeled core remain separate; profiles/site arms remain separate. No ligand fit, geometry success filter or cross-ligand affinity comparison is introduced.",
            "Inherited per-pose geometry, stereo, contacts and clashes are copied unchanged; they do not certify complete physical validity.",
            "Family-cluster intervals describe this purposive sample; no population or model-training-unseen reliability claim."]}


def render_markdown(report):
    if not report["metrics_eligible"]:
        return "# Common-union ranking pending\n\nThe fixed workflow is incomplete. No atomic accuracy is reported.\n"
    from statistics import mean
    text = ["# CHI ranking on the exact common retained union", "",
        "Case/seed macro recovery at2Å; all eligible failures count zero. Each ranker uses the identical retained union.", "",
        "| Profile | Atomic scope | Site arm | Ranker | Eligible pools | Completed pools | Top1 | Top5 | Best retained |",
        "|---|---|---|---|---:|---:|---:|---:|---:|"]
    for profile in PROFILES:
        for scope in ("full_ligand", "modeled_core"):
            for arm in ARMS:
                rows = [r for r in report["conditions"] if r["receptor_profile"] == profile and r["atomic_scope"] == scope and r["site_condition"] == arm]
                for coefficient, method in enumerate(RANK_METHODS):
                    values = [s["ranking_methods"][method] for r in rows for s in r["seeds"]]
                    if not values:
                        continue
                    rates = [mean(v["global"][f"{rank}_success_2A"] for v in values) for rank in ("top1", "top5", "best_returned")]
                    text.append(f"| {profile} | {scope} | {arm} | coeff{coefficient} | {len(values)} | "
                        + str(sum(v["status"] == "completed" for v in values)) + " | " + " | ".join(f"{v:.3f}" for v in rates) + " |")
    text.extend(["", "JSON retains global/local1/2/3Å endpoints, every source pose's unchanged diagnostics and exact provenance, original source orders, and eligible/conditional family statistics.", "", *report["limitations"]])
    return "\n".join(text) + "\n"


def write_report(root, plan_path, freeze_path, summary_output, markdown_output, *, main_report=None, check=False, allow_pending=False):
    report = generate_report(root, plan_path, freeze_path, main_report=main_report, allow_pending=allow_pending)
    for relative, content in ((summary_output, json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"),
                              (markdown_output, render_markdown(report))):
        path = Path(root) / relative
        if check:
            _require(path.read_text() == content, "H2 report replay changed")
        else:
            write_text_atomic(path, content)
    return report


def example_plan(root):
    root = Path(root)
    return {"schema_version": 1, "case_ids": CASES, "profiles": PROFILES, "arms": ARMS, "seeds": SEEDS,
        "ranking_methods": list(RANK_METHODS), "eligible_seed_pools": 108, "eligible_source_method_slots": 216,
        "eligible_rescore_method_slots": 216, "deduplicate_or_filter_poses": False,
        "main_evaluation_plan": MAIN_PLAN, "main_evaluation_freeze": MAIN_FREEZE, "union_scoring_plan": UNION_PLAN,
        "union_scoring_freeze": UNION_FREEZE, "union_output_root": union.DEFAULT_OUTPUT,
        "adapter_source": {"path": "aav_benchmark/paired_common_pose_report.py", "sha256": sha256(Path(__file__))},
        "metric_code_sha256": {path: sha256(root / path) for path in CODE_SOURCES},
        "metric_protocol": PROTOCOL, "statistics": STATISTICS, "implementation_timing": TIMING,
        "frozen_before_union_pose_accuracy": False, "global_four_workflow_guard_required": True,
        "metric_reuse": "Fresh main-reader output only, exact condition/seed/source/model joined with raw file/MODEL/record/input provenance; no refits or new per-pose metrics.",
        "failure_policy": "Any source/fixed-pose failure invalidates the whole union for both coefficients. All108pools/216outcomes remain; no dropping/deduplication/repair.",
        "actual_execution": "The calling pinned global guard exclusively enforces four-workflow release; this adapter enforces its local union gate. CLI is workflow-only. Fresh in-memory main metrics come from the trusted pinned caller; a dictionary is not independent freshness proof. Cached metric-file input is prohibited."}
