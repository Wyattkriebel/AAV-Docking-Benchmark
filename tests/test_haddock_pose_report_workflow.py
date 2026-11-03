"""Synthetic guided workflow gates/lineage; no actual model or source outcomes."""
import copy
from pathlib import Path

import numpy as np
import pytest

import aav_benchmark.haddock_pose_report as report
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import Evidence


SEEDS = (104729, 130363, 155921)


def workflow_fixture(tmp_path):
    conditions, jobs, records = [], [], []
    for profile in ("strict", "repair"):
        for number in range(6):
            for method in ("fpocket", "p2rank", "oracle"):
                case = f"synthetic_{number:02d}"
                cid = f"{profile}/{case}/{method}/box01"
                conditions.append({"condition_id": cid, "profile": profile, "case_id": case, "method": method,
                                   "status": "source_parameterization_gate_failed"})
                for seed in SEEDS:
                    jid = cid + f"/s{seed}"
                    jobs.append({"job_id": jid, "condition_id": cid, "profile": profile, "case_id": case,
                                 "method": method, "seed": seed, "status": "source_parameterization_gate_failed"})
                    records.append({"job_id": jid, "status": "source_parameterization_gate_failed", "terminal": True,
                                    "eligible_workflow_pool": True, "physical_models": 0, "no_native_execution": True})
    prepared = {"conditions": conditions, "jobs": jobs}
    evaluation = {"prepared_plan": {"path": "source/prepared.json", "sha256": "a" * 64},
                  "native_execution_freeze": {"path": "source/freeze.json", "sha256": "b" * 64},
                  "execution_progress": "results/progress.json"}
    progress = {"schema_version": 1, "prepared_plan": evaluation["prepared_plan"], "root_freeze": evaluation["native_execution_freeze"],
                "seed_cell_denominator": 108, "no_accuracy_evaluation": True, "terminal_cells": 108, "all_terminal": True, "jobs": records}
    publish(tmp_path, evaluation, progress)
    return prepared, evaluation, progress


def publish(root, evaluation, progress):
    write_json(root / evaluation["execution_progress"], progress)


def test_all108_explicit_source_failures_are_terminal_without_native_or_accuracy_calls(tmp_path, monkeypatch):
    prepared, evaluation, _ = workflow_fixture(tmp_path)
    def forbidden(*args, **kwargs):
        pytest.fail("Unexecuted source failure entered native/artifact or accuracy evaluation")
    monkeypatch.setattr(report, "read_job", forbidden)
    monkeypatch.setattr(report, "_source_context", forbidden)
    monkeypatch.setattr(report, "score_pool", forbidden)
    checked, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    assert len(prepared["conditions"]) == 36
    assert len(checked) == counts["eligible_cells"] == counts["terminal_cells"] == counts["source_failed_cells"] == 108
    assert pending == [] and counts["native_ready_cells"] == counts["completed_native_cells"] == 0
    assert counts["retained_initial_models"] == counts["retained_refined_models"] == 0
    assert all(value["final_stage"] is None for value in checked.values())


def test_missing_one_or_all_progress_records_preserves_fixed108_denominator(tmp_path):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    missing = progress["jobs"].pop(12)["job_id"]
    progress.update(terminal_cells=107, all_terminal=False)
    publish(tmp_path, evaluation, progress)
    checked, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    assert len(checked) == counts["terminal_cells"] == 107 and counts["eligible_cells"] == 108
    assert pending == [missing + ":missing_terminal_cell"]
    (tmp_path / evaluation["execution_progress"]).unlink()
    checked, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    assert checked == {} and len(pending) == 108
    assert counts["eligible_cells"] == 108 and counts["terminal_cells"] == 0


@pytest.mark.parametrize("kind", ["duplicate", "unknown", "terminal_count", "prepared_sha", "freeze_sha", "denominator", "accuracy_flag", "complete_flag"])
def test_ambiguous_progress_or_changed_source_gates_are_rejected(tmp_path, kind):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    if kind == "duplicate":
        progress["jobs"].append(copy.deepcopy(progress["jobs"][0])); progress["terminal_cells"] = 109
    elif kind == "unknown":
        progress["jobs"][0]["job_id"] = "outside/frozen/denominator"
    elif kind == "terminal_count":
        progress["terminal_cells"] = 107
    elif kind == "prepared_sha":
        progress["prepared_plan"] = {"sha256": "0" * 64}
    elif kind == "freeze_sha":
        progress["root_freeze"] = {"sha256": "0" * 64}
    elif kind == "denominator":
        progress["seed_cell_denominator"] = 107
    elif kind == "accuracy_flag":
        progress["no_accuracy_evaluation"] = False
    else:
        progress["all_terminal"] = False
    publish(tmp_path, evaluation, progress)
    with pytest.raises(ValueError):
        report.read_workflows(Evidence(tmp_path), evaluation, prepared)


@pytest.mark.parametrize("field,value", [("status", "completed"), ("physical_models", 1), ("no_native_execution", False)])
def test_unexecuted_failure_proof_cannot_be_replaced_by_a_success_claim(tmp_path, field, value):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    progress["jobs"][0][field] = value; publish(tmp_path, evaluation, progress)
    with pytest.raises(ValueError, match="unexecuted failure"):
        report.read_workflows(Evidence(tmp_path), evaluation, prepared)


@pytest.mark.parametrize("field", ["terminal", "eligible_workflow_pool"])
def test_incomplete_cell_metadata_stays_pending_despite_aggregate_complete_claim(tmp_path, field):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    missing = progress["jobs"][0]["job_id"]
    progress["jobs"][0][field] = False; publish(tmp_path, evaluation, progress)
    _, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    assert pending == [missing + ":not_terminal"] and counts["terminal_cells"] == 107
    assert counts["eligible_cells"] == 108


def test_ready_cells_missing_individual_receipts_stay_pending_even_when_progress_claims_exit0(tmp_path):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    prepared["conditions"][0]["status"] = "ready"
    ready = prepared["jobs"][:3]
    for index, job in enumerate(ready):
        directory = f"results/native/synthetic_job/s{job['seed']}"
        job.update(status="ready", directory=directory,
                   config={"path": directory + "/run.cfg", "sha256": "c" * 64})
        progress["jobs"][index] = {"job_id": job["job_id"], "terminal": True,
            "eligible_workflow_pool": True, "status": "completed", "process_return_code": 0}
    publish(tmp_path, evaluation, progress)
    checked, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    assert len(checked) == counts["terminal_cells"] == counts["source_failed_cells"] == 105
    assert counts["eligible_cells"] == 108 and counts["native_ready_cells"] == 3 and counts["completed_native_cells"] == 0
    assert pending == sorted(job["job_id"] + ":missing_native_execution" for job in ready)


def test_recovery_progress_requires_separate_matching_authorization_pins(tmp_path):
    prepared, evaluation, progress = workflow_fixture(tmp_path)
    recovery = {"plan_pin": {"sha256": "c" * 64}, "freeze_pin": {"sha256": "d" * 64}}
    progress.update(schema_version=2, recovery_plan=recovery["plan_pin"], execution_recovery_freeze=recovery["freeze_pin"])
    publish(tmp_path, evaluation, progress)
    with pytest.raises(ValueError, match="independently frozen"):
        report.read_workflows(Evidence(tmp_path), evaluation, prepared)
    checked, pending, counts = report.read_workflows(Evidence(tmp_path), evaluation, prepared, recovery=recovery)
    assert len(checked) == counts["terminal_cells"] == 108 and pending == []
    wrong = copy.deepcopy(recovery); wrong["freeze_pin"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="independently frozen"):
        report.read_workflows(Evidence(tmp_path), evaluation, prepared, recovery=wrong)


def test_each_model_uses_its_own_protein_frame_and_wrong_ligand_shift_is_not_fitted_away():
    keys = [("A", 1, "", "ALA", name) for name in ("N", "CA", "C", "O")]
    observed = np.asarray([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    protein_map = {"atoms": [{"output_key": list(key), "assembly_xyz_A": point.tolist()} for key, point in zip(keys, observed)]}
    native = np.asarray([[0., 0., 3.]])
    mapping = {"indices": [0], "copies": [{"copy_id": "synthetic", "coordinates": native}]}
    class Graph:
        scored_maps = [[0]]
        def diagnostics(self, full):
            assert np.asarray(full).shape == (1, 3)
            return {"physical_validity_certified": False}
    rotation = np.asarray([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    translation = np.asarray([20., -30., 40.])
    moved_protein = observed @ rotation + translation
    moved_ligand = native @ rotation + translation
    row = {"rank": 1, "vina_score": -1., "candidate_class_rank": 1, "candidate_original_rank": 1,
           "box_id": "synthetic", "pooled_rank": 1, "coordinates": moved_ligand}
    correct, controls = report.score_pool([row], mapping, keys, moved_protein, protein_map, Graph())
    assert correct[0]["global_RMSD_A"] < 1e-10 and correct[0]["local_interface_RMSD_A"] < 1e-10
    assert correct[0]["global_contacts"]["contact_recall"] == correct[0]["global_contacts"]["contact_precision"] == 1
    assert controls[0]["native_controls"]["native_RMSD_A"] == 0
    wrong, _ = report.score_pool([{**row, "coordinates": moved_ligand + [3., 4., 0.]}], mapping, keys, moved_protein, protein_map, Graph())
    assert wrong[0]["global_RMSD_A"] == pytest.approx(5.)
    assert wrong[0]["local_interface_RMSD_A"] == pytest.approx(5.)


def lineage_fixture(tmp_path):
    run_dir = Path("results/synthetic_job/run")
    selected = {"rigidbody": {"stage": "1_rigidbody", "models": [
        {"file_name": "rigidbody_1.pdb", "actual_CNS_seed": 101},
        {"file_name": "rigidbody_3.pdb.gz", "actual_CNS_seed": 103}]},
        "flexref": {"stage": "2_flexref", "models": [{"file_name": "flexref_1.pdb.gz", "actual_CNS_seed": 103}]}}
    metadata = {"input": [
        {"file_name": "../1_rigidbody/rigidbody_3.pdb", "topology": {"file_name": "receptor.psf"}},
        {"file_name": "rigidbody_1.pdb", "topology": {"file_name": "glycan.psf"}}],
        "output": [{"file_name": "flexref_1.pdb", "ori_name": "../1_rigidbody/rigidbody_3.pdb"}]}
    publish_lineage(tmp_path, run_dir, selected, metadata)
    return run_dir, selected, metadata


def publish_lineage(root, run_dir, selected, metadata):
    path = root / run_dir / selected["flexref"]["stage"] / "io.json"
    write_json(path, metadata); selected["flexref"]["io_sha256"] = sha256(path)


def test_lineage_passes_every_input_once_and_allows_failed_refinement_subset_with_renumbering(tmp_path):
    run_dir, selected, _ = lineage_fixture(tmp_path)
    result = report._verify_lineage(Evidence(tmp_path), run_dir, selected)
    assert result["initial_models_passed_to_refinement"] == 2
    assert result["refined_parent_models"] == ["rigidbody_3.pdb"]
    assert result["refinement_failed_slots_not_replaced"] is True
    assert selected["flexref"]["models"][0]["actual_CNS_seed"] == 103


@pytest.mark.parametrize("kind", ["missing_input", "duplicate_input", "unknown_input", "missing_parent", "unknown_parent", "wrong_seed", "duplicate_parent", "duplicate_output", "unknown_output"])
def test_missing_duplicated_or_invented_refinement_lineage_is_rejected(tmp_path, kind):
    run_dir, selected, metadata = lineage_fixture(tmp_path)
    if kind == "missing_input":
        metadata["input"].pop()
    elif kind == "duplicate_input":
        metadata["input"][1] = copy.deepcopy(metadata["input"][0])
    elif kind == "unknown_input":
        metadata["input"][1]["file_name"] = "rigidbody_99.pdb"
    elif kind == "missing_parent":
        metadata["output"][0].pop("ori_name")
    elif kind == "unknown_parent":
        metadata["output"][0]["ori_name"] = "rigidbody_99.pdb"
    elif kind == "wrong_seed":
        selected["flexref"]["models"][0]["actual_CNS_seed"] = 999
    elif kind == "duplicate_parent":
        selected["flexref"]["models"].append({"file_name": "flexref_2.pdb", "actual_CNS_seed": 103})
        metadata["output"].append({"file_name": "flexref_2.pdb", "ori_name": "rigidbody_3.pdb"})
    elif kind == "duplicate_output":
        metadata["output"].append(copy.deepcopy(metadata["output"][0]))
    else:
        metadata["output"].append({"file_name": "flexref_99.pdb", "ori_name": "rigidbody_1.pdb"})
    publish_lineage(tmp_path, run_dir, selected, metadata)
    with pytest.raises(ValueError):
        report._verify_lineage(Evidence(tmp_path), run_dir, selected)
