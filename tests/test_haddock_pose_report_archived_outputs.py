"""Synthetic recovered-model archive proof; no source/native pose accuracy."""
import json
from pathlib import Path
import shutil
import tarfile

import pytest

from aav_benchmark.haddock_paired_cohort import model_inventory
from aav_benchmark.haddock_pose_report import OutputEvidence, read_job
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import Evidence


def fixture(tmp_path):
    root = Path("results/haddock")
    directory = root / "strict/synthetic/fpocket/s104729"
    run = tmp_path / directory / "run"
    model_text = "ATOM      1   CA ALA A   1       0.000   0.000   0.000  1.00  0.00      A    C\nEND\n"
    for stage in ("0_topoaa", "1_rigidbody", "2_flexref"):
        (run / stage).mkdir(parents=True)
        (run / stage / "params.cfg").write_text("# synthetic retained native parameters\n")
    for name in ("receptor", "glycan"):
        (run / "0_topoaa" / (name + "_haddock.pdb")).write_text(model_text)
        (run / "0_topoaa" / (name + "_haddock.psf")).write_text("synthetic retained topology\n")
    write_json(run / "0_topoaa/io.json", {"input": [], "output": []})
    initial, final = [], []
    for number in range(1, 21):
        source, refined = f"rigidbody_{number}.pdb", f"flexref_{number}.pdb"
        for stage, name in (("1_rigidbody", source), ("2_flexref", refined)):
            (run / stage / name).write_text(model_text)
        initial.append({"file_name": source, "score": -1. - number / 100., "seed": 104729 + number, "unw_energies": {"vdw": -1.}})
        final.append({"file_name": refined, "ori_name": source, "score": -2. - number / 100., "seed": 104729 + number, "unw_energies": {"vdw": -2.}})
    write_json(run / "1_rigidbody/io.json", {"input": [], "output": initial})
    write_json(run / "2_flexref/io.json", {"input": initial, "output": final})
    for name in ("stdout.log", "stderr.log"):
        (tmp_path / directory / name).write_text("synthetic complete native log\n")
    job_id = "strict/synthetic/fpocket/box01/s104729"
    driver = {"path": "driver.py", "sha256": "a" * 64}
    cns = {"path": "tools/not-required", "sha256": "b" * 64}
    job = {"job_id": job_id, "directory": directory.as_posix(), "config": {"path": "run.cfg", "sha256": "c" * 64}}
    prepared = {"driver": driver, "CNS_binary": cns}
    prepared_pin, freeze_pin = {"sha256": "d" * 64}, {"sha256": "e" * 64}
    adoption_pin = {"path": "adoption.json", "sha256": "f" * 64}
    recovery = {"adoption_pin": adoption_pin, "adoption": {"job_id": job_id, "native_self_reported_elapsed_seconds": 373}}
    execution = {"job_id": job_id, "terminal": True, "eligible_workflow_pool": True, "prepared_plan": prepared_pin,
                 "root_freeze": freeze_pin, "helper_source": driver, "CNS_binary": cns, "config": job["config"], "ncores": 1,
                 "no_accuracy_evaluation": True, "failed_model_slots_not_replaced": True, "clustering_performed": False,
                 "requested_initial_model_slots": 20, "status": "complete_outputs_recovered", "elapsed_seconds": 373,
                 "process_return_code": None, "process_CPU_seconds": None, "external_measured_wall_seconds": None, "native_rerun": False,
                 "recovery_receipt": adoption_pin, "stages": model_inventory(run, job_id), "actual_initial_models": 20,
                 "actual_refined_models": 20, "unfilled_initial_slots": 0, "final_score_metadata_complete": True}
    write_json(tmp_path / directory / "execution.json", execution)
    return root, directory, job, prepared, prepared_pin, freeze_pin, recovery, execution


def test_recovered_outputs_require_explicit_unknown_telemetry_and_complete_budget(tmp_path):
    root, directory, job, prepared, prepared_pin, freeze_pin, recovery, execution = fixture(tmp_path)
    result, pending = read_job(Evidence(tmp_path), job, {}, prepared, prepared_pin, freeze_pin, recovery=recovery)
    assert pending is None and result["status"] == "complete_outputs_recovered"
    assert result["execution"]["process_return_code"] is result["execution"]["process_CPU_seconds"] is None
    assert len(result["final_stage"]["models"]) == result["lineage"]["initial_models_passed_to_refinement"] == 20
    with pytest.raises(ValueError, match="Recovered output"):
        read_job(Evidence(tmp_path), job, {}, prepared, prepared_pin, freeze_pin)
    execution["process_return_code"] = 0
    write_json(tmp_path / directory / "execution.json", execution)
    with pytest.raises(ValueError, match="unknown return-code"):
        read_job(Evidence(tmp_path), job, {}, prepared, prepared_pin, freeze_pin, recovery=recovery)


def test_verified_archive_replays_models_io_parameters_and_lineage_without_run_tree(tmp_path):
    root, directory, job, prepared, prepared_pin, freeze_pin, recovery, execution = fixture(tmp_path)
    files = [p for p in sorted((tmp_path / directory).rglob("*")) if p.is_file()]
    archive_path, manifest_path = root / "outputs.tar.gz", root / "outputs_manifest.json"
    manifest = {"files": [{"path": p.relative_to(tmp_path).as_posix(), "sha256": sha256(p), "bytes": p.stat().st_size} for p in files]}
    write_json(tmp_path / manifest_path, manifest)
    with tarfile.open(tmp_path / archive_path, "w:gz") as archive:
        for path in files:
            archive.add(path, arcname=path.relative_to(tmp_path / root).as_posix())
    archive_pin = {"path": archive_path.as_posix(), "sha256": sha256(tmp_path / archive_path)}
    manifest_pin = {"path": manifest_path.as_posix(), "sha256": sha256(tmp_path / manifest_path)}
    shutil.rmtree(tmp_path / directory)
    evidence = OutputEvidence(tmp_path)
    evidence.add_output_archive(archive_pin, manifest_pin, root)
    result, pending = read_job(evidence, job, {}, prepared, prepared_pin, freeze_pin, recovery=recovery)
    assert pending is None and result["status"] == "complete_outputs_recovered"
    assert len(result["final_stage"]["individual_score_order"]) == 20
    assert result["lineage"]["refinement_failed_slots_not_replaced"] is True
    assert not (tmp_path / directory).exists() and not (tmp_path / "tools").exists()
    assert any(path.endswith("receptor_haddock.psf") for path in evidence.sources)
    with pytest.raises(ValueError, match="escaped"):
        OutputEvidence(tmp_path).add_output_archive(archive_pin, manifest_pin, Path("other"))
