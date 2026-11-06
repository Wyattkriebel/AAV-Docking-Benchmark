"""Synthetic indicator tables only; no study reports, fits or native software."""
import copy
import json

import pytest

from aav_benchmark import structural_contrast_report as contrast


def value(completed=True, success=True, status="ok", identity="synthetic"):
    poses = [{"pose_id": identity}] if completed else []
    summary = {"returned_pose_count": len(poses)}
    summary.update({f"{rank}_success_{threshold}A": bool(completed and success)
                    for rank in ("top1", "top5", "best_returned") for threshold in (1, 2, 3)})
    return {"status": status, "poses": poses, "global": copy.deepcopy(summary), "local_interface": copy.deepcopy(summary)}


def conditions(kind):
    result = []
    profiles = (*contrast.PROFILES, "strict_monomer", "repair_monomer") if kind == "assembly" else contrast.PROFILES
    cases = ("3ZWE", "5HZA") if kind == "assembly" else contrast.CASES
    for cid in cases:
        for profile in profiles:
            for site in contrast.SITES:
                row = {"case_id": cid, "profile": profile, "method": site, "receptor_profile": profile,
                       "site_condition": site, "atomic_scope": "modeled_core" if cid == "4R9F" else "full_ligand",
                       "scored_atoms": 11 if cid == "4R9F" else 20, "parent_atoms": 20, "seeds": []}
                for seed in contrast.SEEDS:
                    failed = profile == "strict" and site == "fpocket" and seed == contrast.SEEDS[0]
                    status = "prepared_input_failure" if failed else "ok"
                    outcome = value(not failed, site != "p2rank", status, identity=f"{cid}/{profile}/{site}/{seed}")
                    if kind == "main":
                        item = {"seed": seed, "search_methods": {backend: copy.deepcopy(outcome) for backend in contrast.ENGINES}}
                    elif kind == "haddock":
                        item = {"seed": seed, **outcome, "workflow_status": "source_parameterization_gate_failed" if failed else "completed"}
                    elif kind == "union":
                        outcome.update(status="source_workflow_failure_whole_pair_required" if failed else "completed",
                                       statistics_status="unsupported" if failed else "ok")
                        item = {"seed": seed, "pool_status": outcome["status"],
                                "ranking_methods": {backend: copy.deepcopy(outcome) for backend in contrast.UNION_METHODS}}
                    else:
                        item = {"seed": seed, **outcome}
                    row["seeds"].append(item)
                result.append(row)
    return result


def existing(role):
    rows = []
    for profile in contrast.PROFILES:
        for scope in ("full_ligand", "modeled_core"):
            for alignment, rank, threshold in contrast.EVENTS:
                common = {"confidence": .95, "bootstrap_seed": 20261103, "n_resamples": 0 if scope == "modeled_core" else 10000,
                          "alignment": alignment, "pose_rank": rank, "threshold_A": threshold}
                if role == "haddock":
                    rows.extend({**common, "profile": profile, "atomic_scope": scope, "method_a": "oracle", "method_b": site}
                                for site in ("fpocket", "p2rank"))
                elif role == "assembly_control":
                    if scope == "modeled_core":
                        continue
                    rows.extend({**common, "site_condition": site, "method_a": profile + "_monomer", "method_b": profile}
                                for site in contrast.SITES)
                else:
                    methods = contrast.ENGINES if role == "main_paired" else contrast.UNION_METHODS
                    rows.extend({**common, "receptor_profile": profile, "site_condition": site, "atomic_scope": scope,
                                 "method_a": methods[1], "method_b": methods[0]} for site in contrast.SITES)
    return rows


def packet():
    reports = {role: {"conditions": conditions(kind)} for role, kind in
               (("main_paired", "main"), ("haddock", "haddock"), ("assembly_control", "assembly"), ("fixed_union", "union"))}
    for role, report in reports.items():
        field = "conditional_family_statistics" if role == "haddock" else "paired_family_statistics"
        report[field] = existing(role)
    return reports


def test_normalized_universes_preserve_failures_and_all_fixed_groups():
    reports = packet()
    for role, kind, count in (("main_paired", "main", 216), ("haddock", "haddock", 108),
                              ("assembly_control", "assembly", 72), ("fixed_union", "union", 216)):
        rows = contrast._conditions(reports[role], kind)
        assert len(rows) == count and any(r["status"] != "ok" for r in rows)
    inventory = contrast._existing_statistics(reports)
    assert inventory["haddock"]["entries"] == 144 and inventory["haddock"]["recomputed"] is False


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "seed", "engine", "scope", "failed_success", "count", "nonboolean"])
def test_missing_duplicate_and_failed_success_inconsistencies_are_rejected(mutation):
    report = {"conditions": conditions("main")}; first = report["conditions"][0]
    outcome = first["seeds"][0]["search_methods"][contrast.ENGINES[0]]
    if mutation == "missing": report["conditions"].pop()
    elif mutation == "duplicate": report["conditions"].append(copy.deepcopy(first))
    elif mutation == "seed": first["seeds"][1]["seed"] = first["seeds"][0]["seed"]
    elif mutation == "engine": first["seeds"][0]["search_methods"].pop(contrast.ENGINES[1])
    elif mutation == "scope": first["atomic_scope"] = "modeled_core"
    elif mutation == "failed_success": outcome["global"]["top1_success_2A"] = True
    elif mutation == "count": outcome["global"]["returned_pose_count"] = 1
    else: outcome["global"]["top1_success_2A"] = 0
    with pytest.raises(ValueError): contrast._conditions(report, "main")


def test_haddock_partial_models_do_not_promote_failed_outcomes_and_existing_groups_are_strict():
    reports = packet(); seed = reports["haddock"]["conditions"][0]["seeds"][0]
    seed["physical_initial_models"] = seed["physical_refined_models"] = 20
    rows = contrast._conditions(reports["haddock"], "haddock")
    assert rows[0]["status"] == "failed"
    reports["haddock"]["conditional_family_statistics"][1] = copy.deepcopy(reports["haddock"]["conditional_family_statistics"][0])
    with pytest.raises(ValueError, match="contrast groups"): contrast._existing_statistics(reports)


def test_eligible_family_contrasts_have_correct_direction_and_no_conditional_CI(monkeypatch):
    outcomes = contrast._conditions({"conditions": conditions("main")}, "main")
    outcomes = [r for r in outcomes if r["backend"] == contrast.ENGINES[0]]
    monkeypatch.setattr(contrast, "EVENTS", (("global", "top1", 2),))
    families = {cid: "family_" + cid for cid in contrast.CASES}
    guidance = contrast._contrasts(outcomes, families, kind="guidance")
    strict_fp = next(r for r in guidance if r["receptor_profile"] == "strict" and r["site_condition"] == "fpocket" and r["atomic_scope"] == "full_ligand")
    assert strict_fp["method_a"] == "oracle" and strict_fp["method_b"] == "fpocket"
    assert strict_fp["difference_a_minus_b"] == pytest.approx(1 / 3)
    assert strict_fp["bootstrap_seed"] == 20261103 and strict_fp["n_resamples"] == 10000
    prep = contrast._contrasts(outcomes, families, kind="preparation")
    fp = next(r for r in prep if r["site_condition"] == "fpocket" and r["atomic_scope"] == "full_ligand")
    assert fp["method_a"] == "repair" and fp["method_b"] == "strict"
    assert fp["difference_a_minus_b"] == pytest.approx(1 / 3)
    conditional = fp["conditional_completion"]
    assert conditional["pair_counts"] == {"eligible_pairs": 15, "completed_A": 15, "completed_B": 10,
        "both_completed": 10, "only_A_completed": 5, "only_B_completed": 0, "neither_completed": 0}
    assert all(r["conditional_completed_recovery"] == 1 for r in conditional["arm_specific_rates"])
    assert conditional["conditional_difference_or_CI_computed"] is False
    core = next(r for r in prep if r["atomic_scope"] == "modeled_core")
    assert core["n_families"] == 1 and core["confidence_interval"] is None and core["n_resamples"] == 0


def test_complete_aggregation_is_bounded_keeps_backends_separate_and_preserves_H3(tmp_path, monkeypatch):
    reports = packet(); plan = {"case_to_family": {cid: cid for cid in contrast.CASES}, "release_binding": contrast.BINDING,
                              "source_code_sha256": {"synthetic": "pin"}}
    monkeypatch.setattr(contrast, "validate_plan", lambda *a: (tmp_path, plan, {}))
    monkeypatch.setattr(contrast, "_released_inputs", lambda *a: (reports, {"reports": {"synthetic": "pin"}}))
    calls = []
    def bootstrap(rows, method_a, method_b, **kwargs):
        calls.append((rows, method_a, method_b, kwargs))
        return {"method_a": method_a, "method_b": method_b, "synthetic_test": True}
    monkeypatch.setattr(contrast, "paired_family_bootstrap", bootstrap)
    result = contrast.generate_report(tmp_path, "plan", "freeze", binding_sha256="0" * 64)
    assert len(calls) == 612 and all(k == {"n_resamples": 10000, "seed": 20261103, "confidence": .95, "expected_seeds": list(contrast.SEEDS)} for _, _, _, k in calls)
    assert len(result["main_H1_guidance_contrasts"]) == 288 and len(result["preparation_sensitivity_contrasts"]) == 324
    assert {r["backend"] for r in result["main_H1_guidance_contrasts"]} == set(contrast.ENGINES)
    assert contrast.HYPOTHESES["H3"] == "assay specificity" and result["H3_assay_specificity"]["computed_here"] is False
    assert all(r["contrast_kind"] == "preparation_sensitivity" for r in result["preparation_sensitivity_contrasts"])
    assert result["new_pose_metrics_computed"] is result["native_execution"] is False
    json.dumps(result, allow_nan=False)
