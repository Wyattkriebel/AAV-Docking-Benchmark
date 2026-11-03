"""Synthetic H2 joins/ranks/gates; no experimental inputs or native metrics."""
import copy
import hashlib
from pathlib import Path

import pytest

from aav_benchmark import paired_common_pose_report as report
from aav_benchmark.common_pose_rescoring import canonical_sha256, model_blocks
from aav_benchmark.provenance import write_json
from aav_benchmark.qualified_pose import Evidence


def fixture(tmp_path):
    evidence = Evidence(tmp_path); rows = {}; pools = []; sources = []; records = {}; conditions = []
    for profile in report.PROFILES:
        for cid in report.CASES:
            for arm in report.ARMS:
                condition_id = f"{profile}/{cid}/{arm}/box01"
                directory = Path(condition_id)
                row = {"condition_id": condition_id, "case_id": cid, "receptor_profile": profile, "method": arm,
                       "status": "prepared", "modern_runs": str(directory / "modern/runs.json"),
                       "vina_carb_outputs": {str(s): str(directory / f"vc/s{s}") for s in report.SEEDS}}
                row.update({key: {"sha256": key * 2} for key in ("ligand", "receptor", "preparation", "canonical_atom_map", "initializer_sdf")})
                rows[condition_id] = row
                runs = [{"seed": s, "status": "ok", "source": "modern"} for s in report.SEEDS]
                write_json(tmp_path / row["modern_runs"], runs)
                evidence.read(row["modern_runs"])
                seeds = []
                for seed in report.SEEDS:
                    pool_id = f"{condition_id}/s{seed}"
                    pool = {"pool_id": pool_id, "condition_id": condition_id, "case_id": cid,
                            "receptor_profile": profile, "arm": arm, "seed": seed, "eligible_workflow_pool": True,
                            "source_statuses": {}, "source_record_sha256": {}, "original_source_orders": {}}
                    methods = {}
                    for algorithm_rank, algorithm in enumerate(report.union.SOURCE_METHODS):
                        if algorithm_rank:
                            record_path = Path(row["vina_carb_outputs"][str(seed)]) / "execution.json"
                            record = {"seed": seed, "status": "ok", "source": "vc"}
                            write_json(tmp_path / record_path, record); evidence.read(record_path)
                            pose_path = record_path.parent / "poses.pdbqt"
                        else:
                            record = next(r for r in runs if r["seed"] == seed)
                            pose_path = directory / f"modern/poses_{seed}.pdbqt"
                        raw = b"MODEL 1\nREMARK VINA RESULT: -8.000 0.000 0.000\nENDMDL\n"
                        target = tmp_path / pose_path; target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(raw)
                        evidence.read(pose_path)
                        identity = f"{pool_id}/{algorithm}/p01"
                        source = {"pose_id": identity, "pool_id": pool_id, "condition_id": condition_id,
                            "complex_id": cid, "receptor_profile": profile, "arm": arm, "seed": seed,
                            "source_algorithm": algorithm, "source_algorithm_rank": algorithm_rank,
                            "candidate_class_rank": 1, "candidate_original_rank": 1, "model_rank": 1,
                            "original_method_score": -8., "source_pose_file": str(pose_path),
                            "source_pose_file_sha256": hashlib.sha256(raw).hexdigest(),
                            "source_model_sha256": hashlib.sha256(model_blocks(raw)[1]).hexdigest(),
                            "source_record_sha256": canonical_sha256(record),
                            **{f"source_{name}_sha256": row[key]["sha256"] for name, key in
                               (("ligand", "ligand"), ("receptor", "receptor"), ("preparation", "preparation"), ("canonical_atom_map", "canonical_atom_map"))}}
                        source["pose_identity_sha256"] = canonical_sha256(source); sources.append(source)
                        records[identity] = {"status": "ok", "source": source, "score": {"torsion_normalization_scale": .5,
                            "base_shared0": -8. if algorithm_rank == 0 else -7.,
                            "combined_shared0": -5. if algorithm_rank == 0 else -7.}}
                        pool["source_statuses"][algorithm] = "ok"; pool["source_record_sha256"][algorithm] = canonical_sha256(record)
                        pool["original_source_orders"][algorithm] = [identity]
                        value = .5 if algorithm_rank else 4.
                        pose = {"rank": 1, "pooled_rank": 1, "candidate_class_rank": 1, "candidate_original_rank": 1, "box_id": "box01",
                            "native_method_score": {"method": algorithm, "value": -8., "cross_method_comparable": False},
                            "geometry_diagnostics": {"fixed_parent_graph_and_heavy_map_preserved": True, "stereochemistry_from_pose_3D_matches_input": True},
                            "full_parent_min_protein_distance_A": 2., "full_parent_protein_pairs_below_1p5A": 0,
                            "full_observed_assembly_min_protein_distance_A": 2., "full_observed_assembly_protein_pairs_below_1p5A": 0}
                        for alignment in ("global", "local_interface"):
                            pose.update({alignment + "_RMSD_A": value, alignment + "_qualified_copy_id": "synthetic",
                                         alignment + "_chemical_atom_map_index": 0, alignment + "_contacts": {"contact_recall": .5}})
                        methods[algorithm] = {"status": "ok", "poses": [pose]}
                    pools.append(pool); seeds.append({"seed": seed, "search_methods": methods})
                conditions.append({"condition_id": condition_id, "case_id": cid, "receptor_profile": profile, "site_condition": arm,
                    "atomic_scope": "modeled_core" if cid == "4R9F" else "full_ligand", "scored_atoms": 11 if cid == "4R9F" else 20,
                    "parent_atoms": 20, "common_preparation_status": "prepared", "executed_candidate_budget": 1,
                    "native_controls": ["synthetic"], "site": {"synthetic": True},
                    "exact_common_input_sha256": {k: row[k]["sha256"] for k in ("ligand", "receptor", "preparation", "canonical_atom_map", "initializer_sdf")}, "seeds": seeds})
    index = {"pools": pools, "poses": sources, "terminal_source_method_slots": 216,
             "all216_source_slots_terminal_and_verified": True, "pending_source_slots": []}
    main_plan = {"metric_code_sha256": {"synthetic": "code"}}
    main = {"status": "complete", "metrics_eligible": True, "conditions": conditions,
            "counts": {"eligible_method_outcomes": 216, "terminal_method_outcomes": 216},
            "metric_code_sha256": main_plan["metric_code_sha256"], "source_sha256_inventory": copy.deepcopy(evidence.sources)}
    return evidence, rows, index, records, main_plan, main


def test_exact_source_join_preserves_inherited_metrics_and_distinct_source_identities(tmp_path):
    evidence, rows, index, records, plan, main = fixture(tmp_path)
    before = copy.deepcopy(main)
    joined, conditions = report.join_main_pose_metrics(evidence, main, index, rows, plan)
    assert len(joined) == 216 and main == before
    assert len({s["source_model_sha256"] for s in index["poses"]}) == 1
    evaluated, pools, conflicts = report.summarize_ranked_pools(index, records, joined, conditions)
    assert len(evaluated) == 36 and len(pools) == 108 and not conflicts
    first = evaluated[0]["seeds"][0]["ranking_methods"]
    zero, one = (first[m] for m in report.RANK_METHODS)
    assert not zero["global"]["top1_success_2A"] and one["global"]["top1_success_2A"]
    assert {p["pose_id"] for p in zero["poses"]} == {p["pose_id"] for p in one["poses"]}
    assert zero["global"]["best_returned_RMSD_A"] == one["global"]["best_returned_RMSD_A"] == .5
    for method in first.values():
        for pose in method["poses"]:
            original = joined[pose["pose_id"]]
            assert all(pose[field] == original[field] for field in report.METRIC_FIELDS)
            assert pose["source_model_rank"] == 1


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "score", "candidate", "inventory", "model", "record", "inputs"])
def test_join_rejects_missing_duplicate_or_changed_provenance(tmp_path, mutation):
    evidence, rows, index, _, plan, main = fixture(tmp_path)
    pose = main["conditions"][0]["seeds"][0]["search_methods"][report.union.SOURCE_METHODS[0]]["poses"][0]
    if mutation == "duplicate":
        main["conditions"][0]["seeds"][0]["search_methods"][report.union.SOURCE_METHODS[0]]["poses"].append(copy.deepcopy(pose))
    elif mutation == "missing":
        main["conditions"][0]["seeds"][0]["search_methods"][report.union.SOURCE_METHODS[0]]["poses"] = []
    elif mutation == "score":
        pose["native_method_score"]["value"] = -99.
    elif mutation == "candidate":
        pose["candidate_class_rank"] = 2
    elif mutation == "inventory":
        main["source_sha256_inventory"][index["poses"][0]["source_pose_file"]]["sha256"] = "0" * 64
    else:
        key = {"model": "source_model_sha256", "record": "source_record_sha256", "inputs": "source_ligand_sha256"}[mutation]
        source = index["poses"][0]; source[key] = "0" * 64
        source["pose_identity_sha256"] = canonical_sha256({k: v for k, v in source.items() if k != "pose_identity_sha256"})
    with pytest.raises(ValueError):
        report.join_main_pose_metrics(evidence, main, index, rows, plan)


def test_any_pose_failure_or_case_normalizer_conflict_zeroes_whole_pools(tmp_path):
    evidence, rows, index, records, plan, main = fixture(tmp_path)
    joined, conditions = report.join_main_pose_metrics(evidence, main, index, rows, plan)
    records[index["poses"][0]["pose_id"]]["status"] = "pose_scoring_failed"
    evaluated, pools, _ = report.summarize_ranked_pools(index, records, joined, conditions)
    first = evaluated[0]["seeds"][0]
    assert first["pool_status"] == "pose_scoring_failure_whole_union_required"
    assert len(first["retained_union_pose_ids"]) == 2
    assert all(not r["poses"] and not r["global"]["best_returned_success_3A"] for r in first["ranking_methods"].values())
    assert len(pools) == 108
    records[index["poses"][0]["pose_id"]]["status"] = "ok"
    records[index["poses"][0]["pose_id"]]["score"]["torsion_normalization_scale"] = .6
    _, pools, conflicts = report.summarize_ranked_pools(index, records, joined, conditions)
    assert conflicts == [report.CASES[0]]
    assert sum(p["status"] == "normalization_conflict" for p in pools) == 18


def test_probe_ready_and_pending_paths_never_call_any_metric(tmp_path, monkeypatch):
    evidence, rows, index, records, plan, _ = fixture(tmp_path)
    monkeypatch.setattr(report, "validate_report_plan", lambda *a, **k: (evidence, {}, plan, {}))
    monkeypatch.setattr(report, "_read_workflow", lambda *a: (index, records, rows, 216))
    def forbidden(*args, **kwargs):
        pytest.fail("Workflow validation entered metrics")
    for name in ("join_main_pose_metrics", "summarize_ranked_pools", "ranking_statistics", "accuracy_summary", "paired_family_bootstrap", "aggregate_runs"):
        monkeypatch.setattr(report, name, forbidden)
    monkeypatch.setattr(report.main, "generate_paired_report", forbidden)
    ready = report.probe_workflow(tmp_path, "plan", "freeze")
    assert ready["workflow_ready"] and ready["counts"]["terminal_rescore_method_slots"] == 216
    assert ready["metrics_eligible"] is ready["atomic_accuracy_computed"] is False
    del records[index["poses"][0]["pose_id"]]
    pending = report.generate_report(tmp_path, "plan", "freeze", allow_pending=True)
    assert not pending["workflow_ready"] and pending["counts"]["terminal_rescore_method_slots"] == 214
    with pytest.raises(ValueError, match="terminal"):
        report.generate_report(tmp_path, "plan", "freeze")


def test_actual_union_summary_matches_the_global_guard_contract(tmp_path):
    from aav_benchmark.global_pose_report_gate import pure_union_result
    _, _, index, records, _, _ = fixture(tmp_path)
    ready = pure_union_result(report._workflow_summary(index, records, 216))
    assert ready["status"] == "locally_ready_accuracy_blocked" and ready["accuracy_computed"] is False
    del records[index["poses"][0]["pose_id"]]
    pending = pure_union_result(report._workflow_summary(index, records, 215))
    assert pending["status"] == "pending_incomplete_workflow"


def test_conditional_family_statistics_keep_failures_and_fixed_seed_protocol(tmp_path, monkeypatch):
    evidence, rows, index, records, plan, main = fixture(tmp_path)
    joined, conditions = report.join_main_pose_metrics(evidence, main, index, rows, plan)
    index["pools"][0]["source_statuses"][report.union.SOURCE_METHODS[0]] = "failed"
    evaluated, _, _ = report.summarize_ranked_pools(index, records, joined, conditions)
    selected = [evaluated[0]]; calls = []; original = report.paired_family_bootstrap
    def spy(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)
    monkeypatch.setattr(report, "paired_family_bootstrap", spy)
    output = report.ranking_statistics(selected, {selected[0]["case_id"]: "synthetic_family"})
    assert len(output) == 18 and all(c == {"n_resamples": 10000, "seed": 20261103, "confidence": .95, "expected_seeds": report.SEEDS} for c in calls)
    chosen = next(r for r in output if r["alignment"] == "global" and r["pose_rank"] == "top1" and r["threshold_A"] == 2)
    assert chosen["method_a"] == report.RANK_METHODS[1] and chosen["difference_a_minus_b"] == pytest.approx(2 / 3)
    assert chosen["confidence_interval"] is None and chosen["n_families"] == 1
    one = next(r for r in chosen["eligible_and_conditional_summary"]["overall"] if r["method"] == report.RANK_METHODS[1])
    assert one["eligible_runs"] == 3 and one["completed_runs"] == 2 and one["failed_runs"] == 1
    assert one["conditional_completed_recovery"] == 1.


def test_identical_pose_metrics_stay_distinct_and_unsupported_pools_have_shared_zero_recovery(tmp_path):
    evidence, rows, index, records, plan, main = fixture(tmp_path)
    joined, conditions = report.join_main_pose_metrics(evidence, main, index, rows, plan)
    first, second = (index["poses"][i]["pose_id"] for i in (0, 1))
    for field in report.METRIC_FIELDS:
        joined[second][field] = copy.deepcopy(joined[first][field])
    result, _, _ = report.summarize_ranked_pools(index, records, joined, conditions)
    for value in result[0]["seeds"][0]["ranking_methods"].values():
        assert len(value["poses"]) == 2 and {p["pose_id"] for p in value["poses"]} == {first, second}
    index["pools"][0]["source_statuses"][report.union.SOURCE_METHODS[0]] = "unsupported_preflight"
    result, _, _ = report.summarize_ranked_pools(index, records, joined, conditions)
    for value in result[0]["seeds"][0]["ranking_methods"].values():
        assert value["statistics_status"] == "unsupported" and not value["poses"]
        assert not any(value[a][f"best_returned_success_{t}A"] for a in ("global", "local_interface") for t in (1, 2, 3))


def test_guard_supplied_fresh_main_result_is_reused_after_complete_gate(tmp_path, monkeypatch):
    evidence, rows, index, records, plan, main = fixture(tmp_path)
    monkeypatch.setattr(report, "validate_report_plan", lambda *a, **k: (evidence,
        {"metric_code_sha256": {}, "implementation_timing": report.TIMING}, plan,
        {"case_to_family": {cid: cid for cid in report.CASES}}))
    monkeypatch.setattr(report, "_read_workflow", lambda *a: (index, records, rows, 216))
    monkeypatch.setattr(report, "ranking_statistics", lambda *a: ["synthetic_statistics"])
    def forbidden(*args, **kwargs):
        pytest.fail("Guard-provided main metrics were recomputed")
    monkeypatch.setattr(report.main, "generate_paired_report", forbidden)
    result = report.generate_report(tmp_path, "plan", "freeze", main_report=main)
    assert result["metrics_eligible"] and len(result["conditions"]) == 36 and len(result["per_source_pose_metrics"]) == 216
    assert result["implementation_timing"]["fixed_score_records_already_existed"] == 580


@pytest.mark.parametrize("tamper", [False, True])
def test_read_only_verifier_indexes_and_reconstructs_only_in_temporary_directories(tmp_path, monkeypatch, tamper):
    output = tmp_path / "live"; output.mkdir()
    source = {"pose_id": "pool/modern/p01", "condition_id": "condition"}; index = {"poses": [source]}
    exported, mapping = b"exported", b"mapping"
    record = {"status": "pose_scoring_failed", "exported_ligand_sha256": hashlib.sha256(exported).hexdigest(),
              "export_map_sha256": hashlib.sha256(mapping).hexdigest(), "expected_glycosidic_links": 1}
    score_path = report._record_path(output, source["pose_id"]); write_json(score_path, record)
    manifest = {"records": {source["pose_id"]: hashlib.sha256(score_path.read_bytes()).hexdigest()}}
    write_json(output / "scores_manifest.json", manifest)
    write_json(tmp_path / "inputs.json", {"conditions": [{"condition_id": "condition"}]})
    pin = {"path": "inputs.json", "sha256": hashlib.sha256((tmp_path / "inputs.json").read_bytes()).hexdigest()}
    monkeypatch.setattr(report.union, "validate_union_plan", lambda *a, **k: ({"source_inputs": pin}, "plan"))
    def indexer(root, *, plan_path, output_root):
        assert output_root.resolve() != output.resolve() and not output_root.is_relative_to(output)
        write_json(output_root / "source_pool.json", index)
        return index
    monkeypatch.setattr(report.union, "index_union_pool", indexer)
    def old_records(actual, actual_index, *, recover_orphans):
        assert actual == output and actual_index == index and recover_orphans is False
        return {source["pose_id"]: record}, manifest
    monkeypatch.setattr(report, "_old_records", old_records)
    def exporter(actual_source, row, root, destination):
        assert not destination.is_relative_to(output)
        destination.write_bytes(exported)
        destination.with_suffix(".pdbqt.map.json").write_bytes(b"changed" if tamper else mapping)
        return {"expected_glycosidic_links": ["synthetic"]}
    monkeypatch.setattr(report.union, "export_source_pose", exporter)
    before = {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    if tamper:
        with pytest.raises(ValueError, match="reconstruct"):
            report._read_workflow(Evidence(tmp_path), {"union_output_root": "live"})
    else:
        _, _, _, reconstructed = report._read_workflow(Evidence(tmp_path), {"union_output_root": "live"})
        assert reconstructed == 1
    assert before == {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}
