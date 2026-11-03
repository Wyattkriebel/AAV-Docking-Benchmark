"""Bounded contrasts of released structural-report indicators; no pose metrics."""
from __future__ import annotations

from collections import Counter
import copy
import hashlib
import json
from pathlib import Path

from .statistics import aggregate_runs, paired_family_bootstrap

CASES = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
PROFILES = ("strict", "repair")
SITES = ("fpocket", "p2rank", "oracle")
SEEDS = (104729, 130363, 155921)
ENGINES = ("modern_vina_1_2_7", "vina_carb_1_0")
UNION_METHODS = ("VinaCarb1.0_Vina1.1.2_core_shared_zero_base_coeff0", "VinaCarb1.0_Vina1.1.2_core_shared_zero_base_plus_CHI_coeff1")
ROLES = ("main_paired", "haddock", "assembly_control", "fixed_union")
HYPOTHESES = {"H1": "true-site guidance", "H2": "candidate ranking", "H3": "assay specificity", "H4": "assembly"}
DIRECTORY = "results/structural_pose_v3/global_reports"
GENERATION = "results/structural_pose_v3/global_report_release/generation_receipt.json"
BINDING = "results/structural_pose_v3/global_report_release/released_reports_binding.json"
GLOBAL = "benchmarks/global_pose_report_gate.json"
AUTHORIZATION = "results/structural_pose_v3/global_report_release/freeze_receipt.json"
DEFAULT_PLAN = "benchmarks/structural_contrast_evaluation.proposal.json"
CODE = ("aav_benchmark/structural_contrast_report.py", "aav_benchmark/statistics.py")
ENDPOINTS = {"atomic_strata": ["full_ligand", "modeled_core"], "alignments": ["global", "local_interface"],
             "rank_budgets": ["top1", "top5", "best_returned"], "thresholds_A": [1, 2, 3]}
NEW_CONTRASTS = {"main_H1_guidance": 288, "main_preparation_sensitivity": 216, "haddock_preparation_sensitivity": 108}
HADD_COMPLETE = {"completed", "complete_outputs_recovered"}
HADD_FAILED = {"native_execution_failed", "orchestration_failed", "predictor_unavailable", "oracle_box_unavailable",
               "source_parameterization_gate_failed", "assembly_protocol_unsupported"}
EVENTS = tuple((a, r, t) for a in ("global", "local_interface") for r in ("top1", "top5", "best_returned") for t in (1, 2, 3))
COUNTS = {
    "main_paired": {"eligible_conditions": 36, "eligible_seed_pairs": 108, "eligible_method_outcomes": 216, "terminal_method_outcomes": 216},
    "haddock": {"eligible_cells": 108, "native_ready_cells": 39, "terminal_cells": 108},
    "assembly_control": {"control_conditions": 12, "comparator_conditions": 12, "eligible_seed_outcomes": 72, "terminal_seed_outcomes": 72},
    "fixed_union": {"eligible_seed_pools": 108, "source_method_slots": 216, "rescore_method_slots": 216,
                    "terminal_source_method_slots": 216, "terminal_rescore_method_slots": 216},
}


def require(value, message):
    if not value:
        raise ValueError(message)


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON object key")
        result[key] = value
    return result


def _json(raw):
    def invalid(token):
        raise ValueError("Nonfinite JSON token: " + token)
    return json.loads(raw, object_pairs_hook=_pairs, parse_constant=invalid)


def _read(root, path, digest=None):
    path = Path(path)
    require(not path.is_absolute() and ".." not in path.parts, "Input path must be repository-relative")
    full = root / path
    try:
        ordinary = full.resolve() == full and not full.is_symlink()
    except (RuntimeError, OSError) as error:
        raise ValueError("Input path cannot be resolved safely") from error
    require(ordinary, "Input must be an ordinary file inside the declared root")
    raw = full.read_bytes()
    require(digest is None or _sha(raw) == digest, "Input checksum changed: " + path.as_posix())
    return raw


def _pin(root, path):
    raw = _read(root, path)
    return {"path": str(path), "sha256": _sha(raw), "bytes": len(raw)}


def _checked(root, pin):
    raw = _read(root, pin["path"], pin["sha256"])
    require("bytes" not in pin or pin["bytes"] == len(raw), "Input byte count changed")
    return _json(raw)


def _bound(root, pin):
    require(isinstance(pin, dict) and set(pin) == {"path", "sha256", "bytes"}
            and type(pin["bytes"]) is int and pin["bytes"] >= 0
            and isinstance(pin["path"], str) and isinstance(pin["sha256"], str) and len(pin["sha256"]) == 64
            and all(c in "0123456789abcdef" for c in pin["sha256"]), "Bound release pins require exact path/SHA/byte fields")
    return _checked(root, pin)


def validate_plan(root, plan_path, freeze_path, *, require_freeze=True):
    root = Path(root).resolve(); plan_raw = _read(root, plan_path); plan = _json(plan_raw)
    require(plan["case_ids"] == list(CASES) and plan["seeds"] == list(SEEDS)
            and plan["hypothesis_labels"] == HYPOTHESES and plan["repair_is_hypothesis_H3"] is False,
            "Structural cohort or protocol hypothesis labels changed")
    require(plan["report_directory"] == DIRECTORY and plan["generation_receipt"] == GENERATION
            and plan["release_binding"] == BINDING and plan["report_roles"] == list(ROLES), "Fixed report input paths/schema changed")
    require(plan["statistics"] == {"n_resamples": 10000, "seed": 20261103, "confidence": .95,
            "unit": "source_family", "eligible_failures_zero": True, "conditional_CI": False}, "Inherited contrast statistics changed")
    require(plan["endpoints"] == ENDPOINTS and plan["new_contrasts"] == NEW_CONTRASTS
            and plan["new_structural_accuracy_consulted"] is False, "Bounded contrast choices or pre-accuracy declaration changed")
    require(set(plan["source_code_sha256"]) == set(CODE), "Aggregation runtime source pins incomplete")
    for path, digest in plan["source_code_sha256"].items():
        _read(root, path, digest)
    require(plan["source_code_sha256"][CODE[0]] == _sha(Path(__file__).read_bytes()), "Executing aggregation source changed")
    contract, authorization = _checked(root, plan["global_contract"]), _checked(root, plan["global_authorization"])
    require(plan["global_contract"]["path"] == GLOBAL and plan["global_authorization"]["path"] == AUTHORIZATION
            and _same_pin(authorization["contract"], plan["global_contract"])
            and authorization.get("real_reports_when_all_ready_authorized") is True
            and contract.get("all_four_required") is True, "Frozen global release contract changed")
    reference = _checked(root, plan["family_mapping_source"])
    require(plan["case_to_family"] == reference["case_to_family"] and set(plan["case_to_family"]) == set(CASES), "Frozen source-family mapping changed")
    _read(root, plan["protocol_source"]["path"], plan["protocol_source"]["sha256"])
    if require_freeze:
        require(plan.get("frozen_before_any_structural_accuracy") is True and freeze_path is not None,
                "Contrasts require a separate freeze before any structural accuracy")
        freeze = _json(_read(root, freeze_path))
        require(freeze.get("aggregation_plan_sha256") == _sha(plan_raw)
                and freeze.get("frozen_before_any_structural_accuracy") is True
                and freeze.get("new_structural_accuracy_consulted") is False, "Aggregation pre-accuracy freeze changed")
    return root, plan, contract


def _same_pin(a, b):
    return all(a.get(k) == b.get(k) for k in ("path", "sha256"))


def _generation(root, plan, pin):
    require(isinstance(pin, dict) and pin.get("path") == plan["generation_receipt"], "Unexpected global generation receipt pin/path")
    receipt = _bound(root, pin)
    require(receipt.get("schema_version") == 1 and receipt.get("status") == "all_four_reports_generated"
            and receipt.get("all_four_ready") is True and receipt.get("accuracy_computed") is True
            and receipt.get("real_reports_generated") is True and receipt.get("output") == plan["report_directory"]
            and _same_pin(receipt.get("contract", {}), plan["global_contract"])
            and _same_pin(receipt.get("probe_freeze", {}), plan["global_authorization"]), "Global generation was not successfully released")
    readers = receipt.get("readers", [])
    by_role = {r["reader"]: r for r in readers}
    require(len(readers) == len(by_role) == 4 and set(by_role) == set(ROLES), "Global generation receipt lost a reader role")
    for role, row in by_role.items():
        require(row.get("status") == "locally_ready_accuracy_blocked" and row.get("accuracy_computed") is False
                and row.get("source_context_not_entered") is True, "Generation receipt lacks blocked pre-accuracy readiness")
        require(all(type(row["counts"].get(k)) is int and row["counts"][k] == v for k, v in COUNTS[role].items()),
                "Global readiness denominator changed")
        if role == "fixed_union":
            require(row.get("pure_validation_only") is True, "Fixed-union readiness was not a pure validation probe")
    return pin, receipt, by_role


def _report_header(value, role, contract, readiness):
    require(value.get("schema_version") == 1 and value.get("status") == "complete" and value.get("metrics_eligible") is True,
            "A globally released report is incomplete")
    counts = value["counts"]
    require(all(type(counts.get(k)) is int and counts[k] == v for k, v in COUNTS[role].items()), "Released report denominator changed")
    require(all(counts.get(k) == v for k, v in readiness["counts"].items()), "Report counters differ from global readiness")
    spec = contract["pure_union"] if role == "fixed_union" else contract["readers"][role]
    require(value.get("metric_code_sha256") == spec["metric_sources"], "Released report metric-code provenance changed")


def _packet_directory(root, plan):
    directory = root / plan["report_directory"]
    try:
        ordinary = directory.resolve() == directory and not directory.is_symlink()
    except (RuntimeError, OSError) as error:
        raise ValueError("Global report directory cannot be resolved safely") from error
    require(ordinary, "Global report directory must not be symlinked")
    files = list(directory.iterdir())
    require({p.name for p in files} == {role + ".json" for role in ROLES}
            and all(not p.is_symlink() and p.is_file() for p in files), "Global report directory is not the exact four-file packet of ordinary files")


def bind_released_reports(root, plan_path, freeze_path):
    """Create a separate, explicitly post-generation four-file hash binding."""
    root, plan, contract = validate_plan(root, plan_path, freeze_path)
    generation_pin = _pin(root, plan["generation_receipt"])
    _, _, readiness = _generation(root, plan, generation_pin)
    _packet_directory(root, plan)
    pins = {}
    for role in ROLES:
        pin = _pin(root, plan["report_directory"] + "/" + role + ".json")
        _report_header(_checked(root, pin), role, contract, readiness[role]); pins[role] = pin
    binding = {"schema_version": 1, "status": "four_released_reports_bound", "generation_receipt": generation_pin,
        "global_contract": plan["global_contract"], "global_authorization": plan["global_authorization"], "reports": pins,
        "hash_binding_timing": "post_generation; the frozen generation receipt itself does not carry report-file hashes",
        "new_pose_metrics_computed": False, "native_execution": False}
    _write_new(root, plan["release_binding"], binding)
    return _pin(root, plan["release_binding"])


def _released_inputs(root, plan, contract, binding_sha256):
    require(isinstance(binding_sha256, str) and len(binding_sha256) == 64, "The release binding needs its explicit SHA256")
    binding = _json(_read(root, plan["release_binding"], binding_sha256))
    require(binding.get("schema_version") == 1 and binding.get("status") == "four_released_reports_bound"
            and set(binding["reports"]) == set(ROLES) and _same_pin(binding["global_contract"], plan["global_contract"])
            and _same_pin(binding["global_authorization"], plan["global_authorization"])
            and binding.get("new_pose_metrics_computed") is False and binding.get("native_execution") is False
            and binding.get("hash_binding_timing") == "post_generation; the frozen generation receipt itself does not carry report-file hashes",
            "Post-generation report binding changed")
    _, _, readiness = _generation(root, plan, binding["generation_receipt"])
    _packet_directory(root, plan)
    reports = {}
    for role, pin in binding["reports"].items():
        require(pin["path"] == plan["report_directory"] + "/" + role + ".json", "Bound report role/path changed")
        value = _bound(root, pin); _report_header(value, role, contract, readiness[role]); reports[role] = value
    return reports, binding


def _write_new(root, relative, value):
    relative = Path(relative)
    require(not relative.is_absolute() and ".." not in relative.parts, "Output must be repository-relative")
    path = root / relative
    require(path.resolve() == path and path.is_relative_to(root) and not path.is_relative_to(root / DIRECTORY),
            "Output path escaped its root or entered the immutable four-report packet")
    raw = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    if path.exists():
        require(path.read_bytes() == raw, "Refusing to replace an earlier aggregation artifact")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)


def _outcome(value, status, kind):
    if kind == "haddock":
        require(status in HADD_COMPLETE | HADD_FAILED, "Unknown/nonterminal HADDOCK workflow status")
        completed = status in HADD_COMPLETE
        normalized = "ok" if completed else "failed"
    elif kind == "union":
        require(status in {"completed", "source_workflow_failure_whole_pair_required", "normalization_conflict",
                "source_workflow_failure_no_poses", "pose_scoring_failure_whole_union_required"}, "Unknown/nonterminal union status")
        completed = status == "completed"; normalized = value["statistics_status"]
        require(normalized in {"ok", "failed", "unsupported"} and (normalized == "ok") == completed, "Union completion category changed")
    else:
        require(status in {"ok", "failed", "prepared_input_failure", "unsupported_preflight"}, "Unknown/nonterminal main workflow status")
        completed = status == "ok"
        normalized = "ok" if completed else "unsupported" if status in {"prepared_input_failure", "unsupported_preflight"} else "failed"
    poses = value["poses"]
    require(isinstance(poses, list) and bool(poses) == completed, "Completion status and retained scored poses disagree")
    success = {}
    for alignment in ("global", "local_interface"):
        summary = value[alignment]
        require(type(summary["returned_pose_count"]) is int and summary["returned_pose_count"] == len(poses), "Summary pose count differs")
        for rank in ("top1", "top5", "best_returned"):
            flags = [summary[f"{rank}_success_{t}A"] for t in (1, 2, 3)]
            require(all(type(flag) is bool for flag in flags) and flags == sorted(flags), "Invalid/nonmonotone recovery indicators")
            require(completed or not any(flags), "Failed/unsupported workflow has nonzero success")
            for threshold, flag in zip((1, 2, 3), flags):
                success[alignment, rank, threshold] = flag
        for threshold in (1, 2, 3):
            require(summary[f"top1_success_{threshold}A"] <= summary[f"top5_success_{threshold}A"]
                    <= summary[f"best_returned_success_{threshold}A"], "Recovery rank budgets are inconsistent")
    return normalized, success


def _conditions(report, kind):
    conditions = report["conditions"]; output = []; seen = set()
    profiles = ("strict", "repair", "strict_monomer", "repair_monomer") if kind == "assembly" else PROFILES
    cases = ("3ZWE", "5HZA") if kind == "assembly" else CASES
    for condition in conditions:
        profile = condition["profile"] if kind in {"haddock", "assembly"} else condition["receptor_profile"]
        site = condition["method"] if kind in {"haddock", "assembly"} else condition.get("site_condition")
        cid = condition["case_id"]; key = cid, profile, site
        require(key not in seen and cid in cases and profile in profiles and site in SITES, "Missing/duplicate/unknown condition group")
        seen.add(key)
        scope = "modeled_core" if cid == "4R9F" else "full_ligand"
        require(condition["atomic_scope"] == scope and type(condition["scored_atoms"]) is int
                and type(condition["parent_atoms"]) is int and 0 < condition["scored_atoms"] <= condition["parent_atoms"]
                and (condition["scored_atoms"] == 11 if scope == "modeled_core" else condition["scored_atoms"] == condition["parent_atoms"]),
                "Full/core atomic stratum changed")
        seeds = condition["seeds"]
        require(len(seeds) == 3 and all(type(s["seed"]) is int for s in seeds) and {s["seed"] for s in seeds} == set(SEEDS),
                "Missing/duplicate/unknown algorithm seed")
        for seed in seeds:
            if kind == "main":
                methods = seed["search_methods"]; require(set(methods) == set(ENGINES), "Missing/extra main search engine")
            elif kind == "union":
                methods = seed["ranking_methods"]
                require(set(methods) == set(UNION_METHODS) and all(v["status"] == seed["pool_status"] for v in methods.values()),
                        "Union ranker identities/statuses changed")
                ids = [[p["pose_id"] for p in methods[m]["poses"]] for m in UNION_METHODS]
                require(len(ids[0]) == len(set(ids[0])) and len(ids[1]) == len(set(ids[1])) and set(ids[0]) == set(ids[1]),
                        "Union rankers no longer retain identical distinct pose IDs")
                require(all(methods[UNION_METHODS[0]][a][f"best_returned_success_{t}A"] == methods[UNION_METHODS[1]][a][f"best_returned_success_{t}A"]
                        for a in ("global", "local_interface") for t in (1, 2, 3)), "Union best-retained recovery changed")
            else:
                methods = {"haddock3_guided_box_air" if kind == "haddock" else "modern_vina_1_2_7": seed}
            for backend, value in methods.items():
                status = seed["workflow_status"] if kind == "haddock" else value["status"]
                normalized, success = _outcome(value, status, kind)
                output.append({"backend": backend, "case_id": cid, "profile": profile, "site": site,
                    "atomic_scope": scope, "seed": seed["seed"], "status": normalized, "workflow_status": status, "success": success})
    require(seen == {(cid, profile, site) for cid in cases for profile in profiles for site in SITES}, "Incomplete fixed condition universe")
    return output


def _existing_statistics(reports):
    specs = {"main_paired": ("paired_family_statistics", 216), "haddock": ("conditional_family_statistics", 144),
             "fixed_union": ("paired_family_statistics", 216), "assembly_control": ("paired_family_statistics", 108)}
    result = {}
    for role, (field, count) in specs.items():
        values = reports[role][field]
        require(isinstance(values, list) and len(values) == count, "Existing contrast inventory is incomplete")
        keys = []
        for row in values:
            require(row.get("bootstrap_seed") == 20261103 and row.get("confidence") == .95
                    and type(row.get("n_resamples")) is int and row["n_resamples"] in {0, 10000}, "Existing bootstrap convention changed")
            keys.append(tuple(row.get(k) for k in ("profile", "receptor_profile", "site_condition", "atomic_scope", "alignment", "pose_rank", "threshold_A", "method_a", "method_b")))
        if role == "haddock":
            expected = {(p, None, None, scope, a, r, t, "oracle", m) for p in PROFILES for scope in ("full_ligand", "modeled_core") for a, r, t in EVENTS for m in ("fpocket", "p2rank")}
        elif role == "assembly_control":
            expected = {(None, None, site, None, a, r, t, p + "_monomer", p) for p in PROFILES for site in SITES for a, r, t in EVENTS}
        else:
            method_a, method_b = (ENGINES[1], ENGINES[0]) if role == "main_paired" else (UNION_METHODS[1], UNION_METHODS[0])
            expected = {(None, p, site, scope, a, r, t, method_a, method_b) for p in PROFILES for site in SITES for scope in ("full_ligand", "modeled_core") for a, r, t in EVENTS}
        require(len(set(keys)) == count and set(keys) == expected, "Missing/duplicate/unknown existing contrast groups")
        result[role] = {"report_field": field, "entries": count, "recomputed": False}
    return result


def _conditional(rows, method_a, method_b):
    by_pair = {}
    for row in rows:
        by_pair.setdefault((row["complex_id"], row["seed"]), {})[row["method"]] = row["status"] == "ok"
    flags = list(by_pair.values())
    counts = {"eligible_pairs": len(flags), "completed_A": sum(v[method_a] for v in flags),
        "completed_B": sum(v[method_b] for v in flags), "both_completed": sum(v[method_a] and v[method_b] for v in flags),
        "only_A_completed": sum(v[method_a] and not v[method_b] for v in flags),
        "only_B_completed": sum(v[method_b] and not v[method_a] for v in flags),
        "neither_completed": sum(not v[method_a] and not v[method_b] for v in flags)}
    summary = aggregate_runs(rows, expected_seeds=list(SEEDS))
    fields = ("method", "eligible_runs", "completed_runs", "failed_runs", "unsupported_runs", "families_with_completed_runs", "conditional_completed_recovery")
    return {"pair_counts": counts, "arm_specific_rates": [{k: row[k] for k in fields} for row in summary["overall"]],
        "per_case_completion": [{k: row[k] for k in ("method", "family_id", "complex_id", "eligible_runs", "completed_runs", "conditional_completed_recovery")}
                                for row in summary["per_complex"]],
        "denominators_may_differ": True, "conditional_difference_or_CI_computed": False,
        "conditional_weighting": "Completed seeds within case, cases within source family; only groups with completed outcomes contribute",
        "interpretation": "Arm-specific completed-only descriptive rates have different possible denominators; they are not the eligible paired contrast and carry no new conditional interval."}


def _contrasts(outcomes, families, *, kind):
    result = []
    for backend in sorted({row["backend"] for row in outcomes}):
        for scope in ("full_ligand", "modeled_core"):
            selected = [row for row in outcomes if row["backend"] == backend and row["atomic_scope"] == scope]
            if kind == "guidance":
                groups = [(p, predicted, "oracle", predicted) for p in PROFILES for predicted in ("fpocket", "p2rank")]
            else:
                groups = [(site, None, "repair", "strict") for site in SITES]
            for group, predicted, method_a, method_b in groups:
                subset = [r for r in selected if (r["profile"] == group and r["site"] in {method_a, method_b})] if kind == "guidance" else [r for r in selected if r["site"] == group]
                require(len(subset) == (5 if scope == "full_ligand" else 1) * 3 * 2, "Contrast case/seed pair universe changed")
                for alignment, rank, threshold in EVENTS:
                    rows = [{"family_id": families[r["case_id"]], "complex_id": r["case_id"], "seed": r["seed"],
                        "method": r["site"] if kind == "guidance" else r["profile"], "status": r["status"],
                        "success": r["success"][alignment, rank, threshold]} for r in subset]
                    summary = paired_family_bootstrap(rows, method_a, method_b, n_resamples=10000, seed=20261103,
                        confidence=.95, expected_seeds=list(SEEDS))
                    summary.update(backend=backend, atomic_scope=scope, alignment=alignment, pose_rank=rank, threshold_A=threshold,
                        contrast_kind="H1_true_site_guidance" if kind == "guidance" else "preparation_sensitivity",
                        receptor_profile=group if kind == "guidance" else None, site_condition=predicted if kind == "guidance" else group,
                        estimand="Eligible workflow recovery A-minus-B; failed/unsupported outcomes remain zero",
                        conditional_completion=_conditional(rows, method_a, method_b),
                        failure_status_counts=dict(Counter(r["workflow_status"] for r in subset if r["status"] != "ok")))
                    result.append(summary)
    return result


def generate_report(root, plan_path, freeze_path, *, binding_sha256):
    root, plan, contract = validate_plan(root, plan_path, freeze_path)
    reports, binding = _released_inputs(root, plan, contract, binding_sha256)
    normalized = {role: _conditions(reports[role], kind) for role, kind in
                  (("main_paired", "main"), ("haddock", "haddock"), ("assembly_control", "assembly"), ("fixed_union", "union"))}
    existing = _existing_statistics(reports)
    guidance = _contrasts(normalized["main_paired"], plan["case_to_family"], kind="guidance")
    preparation = _contrasts(normalized["main_paired"], plan["case_to_family"], kind="preparation")
    preparation.extend(_contrasts(normalized["haddock"], plan["case_to_family"], kind="preparation"))
    require(len(guidance) == 288 and len(preparation) == 324, "Bounded new contrast count changed")
    return {"schema_version": 1, "status": "complete", "hypothesis_labels": copy.deepcopy(HYPOTHESES),
        "main_H1_guidance_contrasts": guidance, "preparation_sensitivity_contrasts": preparation,
        "existing_statistics_referenced_not_recomputed": existing,
        "H3_assay_specificity": {"scope": "Separate assay reports; preparation sensitivity is not H3", "computed_here": False},
        "input_provenance": {"release_binding": {"path": plan["release_binding"], "sha256": binding_sha256}, **binding},
        "aggregation_code_sha256": plan["source_code_sha256"], "new_pose_metrics_computed": False,
        "native_execution": False, "old156_development_reinterpreted": False,
        "limitations": ["Only already released recovery indicators are aggregated; no pose metrics, contacts, fits or native scores are recalculated.",
            "Engines and HADDOCK remain separate. HADDOCK patch AIR guidance is not equivalent to a hard Cartesian Vina box.",
            "Five full-ligand cases and one11-atom core remain separate. All eligible failures count zero; completed-only rates can use different arm-specific denominators and have no new conditional contrast or interval.",
            "Whole source families are resampled after seed/case averaging. This purposive cohort supplies descriptive intervals, not population or model-training-unseen reliability.",
            "The four file hashes are explicitly post-generation bindings; the frozen generation receipt itself contains no report-file hashes."]}


def example_plan(root):
    root = Path(root).resolve()
    family_source = _pin(root, "benchmarks/structural_pose_reference_evaluation.json")
    plan = {"schema_version": 1, "case_ids": list(CASES), "seeds": list(SEEDS), "hypothesis_labels": HYPOTHESES,
        "repair_is_hypothesis_H3": False, "report_roles": list(ROLES), "report_directory": DIRECTORY,
        "generation_receipt": GENERATION, "release_binding": BINDING,
        "global_contract": _pin(root, GLOBAL), "global_authorization": _pin(root, AUTHORIZATION),
        "family_mapping_source": family_source, "case_to_family": _checked(root, family_source)["case_to_family"],
        "protocol_source": _pin(root, "PROTOCOL.md"), "source_code_sha256": {p: _pin(root, p)["sha256"] for p in CODE},
        "statistics": {"n_resamples": 10000, "seed": 20261103, "confidence": .95, "unit": "source_family", "eligible_failures_zero": True, "conditional_CI": False},
        "endpoints": ENDPOINTS, "new_contrasts": NEW_CONTRASTS,
        "existing_contrasts": {"HADDOCK_H1_guidance": "reference144existing", "H2": "reference216existing", "H4": "reference108existing"},
        "frozen_before_any_structural_accuracy": False, "new_structural_accuracy_consulted": False,
        "release_binding_timing": "Post-generation file hashes, separately bound to the successful frozen global generation result; never represented as contemporaneous producer hashes."}
    return copy.deepcopy(plan)
