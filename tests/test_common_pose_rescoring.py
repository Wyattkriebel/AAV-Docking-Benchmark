"""Fixed-pool controls keep failures, source identities and score pairs intact."""
import copy
import hashlib
import json
from pathlib import Path

import pytest

from aav_benchmark.common_pose_rescoring import (
    _old_records, _record_path, canonical_sha256, merge_pose_index, model_blocks,
    rank_complete_pools, validate_plan, verify_score_record,
)
from aav_benchmark.provenance import sha256, write_json

ROOT = Path(__file__).resolve().parents[1]


def source(rank=1, box=1):
    row = {"pose_id": f"1GNY/fpocket/box{box:02d}/s104729/p{rank:02d}", "condition_id": "1GNY",
           "complex_id": "1GNY", "method": "fpocket", "box_id": f"box{box:02d}", "seed": 104729,
           "candidate_class_rank": box, "candidate_original_rank": box, "model_rank": rank,
           "original_vina_score": -8.0 + rank + box / 10, "source_receptor_sha256": "a" * 64}
    row["pose_identity_sha256"] = canonical_sha256(row)
    return row


def successful(row, base=-7.0, chi=0.0, normalization=.5):
    raw = {"base_shared0": base, "chi_coefficient1": chi, "combined_shared0": base + chi,
           "chi_coefficient0": 0, "intramolecular_reference": 0, "nonbonded_cutoff_A": 8,
           "torsion_normalization_scale": normalization, "glycosidic_links": 1, "fully_supported_links": 1,
           "original_score_only_base": base + .5, "pose_intramolecular_energy": -1.}
    text = json.dumps(raw) + "\n"
    recognition = "CHI energy penalties applied to NAG1 and GAL2 linkage.\n"
    score = {**raw, "raw_stdout": text, "raw_stdout_sha256": hashlib.sha256(text.encode()).hexdigest(),
             "raw_stderr": "", "raw_stderr_sha256": hashlib.sha256(b"").hexdigest(),
             "recognition_log": recognition, "recognition_log_sha256": hashlib.sha256(recognition.encode()).hexdigest(),
             "ligand_sha256": "b" * 64, "receptor_sha256": row["source_receptor_sha256"],
             "coordinate_optimization": False, "native_coordinates_accessed": False, "executable_sha256": "c" * 64}
    return {"source": row, "method_lock_sha256": "m" * 64, "coordinate_optimization": False,
            "native_coordinates_used": False, "status": "ok", "score": score,
            "exported_ligand_sha256": "b" * 64, "expected_glycosidic_links": 1}


def index(rows):
    return {"poses": rows, "conditions": [{"condition_id": "1GNY", "complex_id": "1GNY", "method": "fpocket",
            "boxes": [{"box_id": f"box{box:02d}", "candidate_class_rank": box, "status": "prepared",
                       "runs": [{"seed": seed, "status": "ok"} for seed in (104729, 130363, 155921)]}
                      for box in sorted({r["candidate_class_rank"] for r in rows})]}]}


def test_source_model_fingerprints_include_every_atom_and_reject_partial_blocks():
    raw = b"MODEL 1\nREMARK VINA RESULT: -7 0 0\nATOM  source1\nENDMDL\nMODEL 2\nATOM  source2\nENDMDL\n"
    blocks = model_blocks(raw)
    assert list(blocks) == [1, 2]
    assert b"source1" in blocks[1] and b"source2" in blocks[2]
    assert hashlib.sha256(blocks[1]).hexdigest() != hashlib.sha256(blocks[2]).hexdigest()
    for bad in (raw[:-8], raw.replace(b"MODEL 2", b"MODEL 3"), b"ATOM  outside\n" + raw):
        with pytest.raises(ValueError):
            model_blocks(bad)


def test_source_resume_appends_but_never_drops_or_rewrites_prior_pose_identities():
    first, second = source(), source(rank=2)
    assert len(merge_pose_index([first], [second, first])) == 2
    with pytest.raises(ValueError, match="disappeared"):
        merge_pose_index([first], [second])
    with pytest.raises(ValueError, match="changed"):
        merge_pose_index([first], [{**first, "original_vina_score": -999.}])
    with pytest.raises(ValueError, match="Duplicate"):
        merge_pose_index([], [first, first])


def test_raw_score_reparse_cannot_accept_changed_scores_sources_or_coordinate_policy():
    row = source();record = successful(row)
    assert verify_score_record(record, row, "m" * 64) is record
    for change in (lambda x: x["score"].update(base_shared0=-99.),
                   lambda x: x["source"].update(pose_identity_sha256="changed"),
                   lambda x: x.update(native_coordinates_used=True),
                   lambda x: x["score"].update(recognition_log="altered")):
        altered = copy.deepcopy(record);change(altered)
        with pytest.raises(ValueError):
            verify_score_record(altered, row, "m" * 64)


def test_an_unscored_or_failed_pose_never_disappears_to_improve_a_pool():
    a, b = source(), source(box=2)
    records = {a["pose_id"]: successful(a)}
    pools, _ = rank_complete_pools(index([a, b]), records, {"1GNY": {"status": "supported"}})
    one = next(r for r in pools if r["seed"] == 104729 and r["site_budget"] == 1)
    five = next(r for r in pools if r["seed"] == 104729 and r["site_budget"] == 5)
    assert one["status"] == "completed"
    assert five["status"] == "pending_scoring" and five["coefficient1_order"] == []
    records[b["pose_id"]] = {"source": b, "status": "pose_scoring_failed"}
    five = next(r for r in rank_complete_pools(index([a, b]), records, {})[0]
                if r["seed"] == 104729 and r["site_budget"] == 5)
    assert five["status"] == "pose_scoring_failure_complete_pool_required"
    assert len(five["original_vina_order"]) == 2 and five["coefficient0_order"] == []


def test_static_unsupported_and_normalization_conflicts_invalidate_all_case_pools():
    a, b = source(), source(rank=2)
    records = {a["pose_id"]: successful(a), b["pose_id"]: successful(b, normalization=.6)}
    pools, conflicts = rank_complete_pools(index([a, b]), records, {})
    assert conflicts == ["1GNY"]
    assert all(r["status"] == "unsupported_case_normalization_conflict" for r in pools if r["seed"] == 104729)
    records[b["pose_id"]] = successful(b)
    pools, _ = rank_complete_pools(index([a, b]), records, {"1GNY": {"status": "unsupported_case"}})
    assert all(r["status"] == "unsupported_case" for r in pools if r["seed"] == 104729)


def test_paired_ranking_uses_same_pose_set_and_prespecified_class_model_ties():
    a, b, c = source(), source(rank=2), source(box=2)
    records = {a["pose_id"]: successful(a, -10., 5.), b["pose_id"]: successful(b, -10., 0.),
               c["pose_id"]: successful(c, -10., 0.)}
    pools, _ = rank_complete_pools(index([c, b, a]), records, {"1GNY": {"status": "supported"}})
    pool = next(r for r in pools if r["seed"] == 104729 and r["site_budget"] == 5)
    assert pool["coefficient0_order"] == [a["pose_id"], b["pose_id"], c["pose_id"]]
    assert pool["coefficient1_order"] == [b["pose_id"], c["pose_id"], a["pose_id"]]
    assert set(pool["original_vina_order"]) == set(pool["coefficient0_order"]) == set(pool["coefficient1_order"])
    assert pool["atomic_accuracy_computed"] is False


def test_interrupted_record_write_is_recoverable_but_changed_manifested_record_is_not(tmp_path):
    row = source();record = successful(row)
    prepared_index = {"poses": [row], "configuration_sha256": "f" * 64, "method_lock_sha256": "m" * 64}
    path = _record_path(tmp_path, row["pose_id"]);write_json(path, record)
    with pytest.raises(ValueError, match="manifest"):
        _old_records(tmp_path, prepared_index)
    records, manifest = _old_records(tmp_path, prepared_index, recover_orphans=True)
    assert list(records) == [row["pose_id"]]
    write_json(tmp_path / "scores_manifest.json", manifest)
    record["score"]["base_shared0"] = -999.;write_json(path, record)
    with pytest.raises(ValueError, match="changed"):
        _old_records(tmp_path, prepared_index, recover_orphans=True)


def test_failed_native_invocation_keeps_all_raw_streams_and_hashes(tmp_path, monkeypatch):
    import subprocess
    from aav_benchmark import common_pose_rescoring as module
    install = tmp_path / "installation";install.mkdir()
    executable = install / "scorer";executable.write_text("frozen source executable")
    root = tmp_path / "repo";(root / "benchmarks").mkdir(parents=True)
    lock = root / "benchmarks/additional_method.lock.json";lock.write_text("frozen method")
    ligand = tmp_path / "ligand.pdbqt";ligand.write_text("unchanged pose")
    receptor = tmp_path / "receptor.pdbqt";receptor.write_text("unchanged receptor")
    write_json(install / "fixed_pose_scorer.json", {"executable": "scorer", "executable_sha256": sha256(executable),
                                                   "method_lock_sha256": sha256(lock)})
    monkeypatch.setattr(module, "validate_installation", lambda *args, **kwargs: executable)
    def fail(command, **kwargs):
        Path(kwargs["cwd"], "VC_log.txt").write_text("Carbohydrate linkage NOT applied\n")
        return subprocess.CompletedProcess(command, -6, "partial raw score\n", "native validation failed\n")
    monkeypatch.setattr(module.subprocess, "run", fail)
    with pytest.raises(module.NativeScoringFailure) as captured:
        module._invoke_fixed_pose(receptor, ligand,
                                  {"ligand_sha256": sha256(ligand), "native_coordinates_accessed": False}, install, root=root)
    attempt = captured.value.attempt
    assert attempt["raw_stdout"] == "partial raw score\n" and attempt["returncode"] == -6
    assert attempt["recognition_log"] == "Carbohydrate linkage NOT applied\n"
    row = source()
    row["source_receptor_sha256"] = sha256(receptor)
    record = {"source": row, "method_lock_sha256": "m" * 64, "coordinate_optimization": False,
              "native_coordinates_used": False, "status": "pose_scoring_failed", "error": "Native failure",
              "native_attempt": attempt, "exported_ligand_sha256": sha256(ligand)}
    assert verify_score_record(record, row, "m" * 64) is record
    record["native_attempt"]["raw_stderr"] = "edited"
    with pytest.raises(ValueError, match="stream changed"):
        verify_score_record(record, row, "m" * 64)


def test_frozen_common_pool_plan_is_current_and_identity_map_does_not_carry_native_positions():
    assert validate_plan(ROOT)["eligible_complexes"] == ["1GNY", "1KJL", "1QFO", "2J1T", "1PWB", "2VUZ"]
    from aav_benchmark.common_pose_rescoring import source_atom_map
    mapping = source_atom_map(ROOT, "1KJL")
    assert mapping["status"] == "mapped"
    assert all(set(row) == {"name", "source_key"} for row in mapping["atom_map"])
    assert mapping["native_coordinate_fields_forwarded"] is False


def test_independent_pucker_convention_controls_and_planar_limit():
    from aav_benchmark.ring_geometry import analytic_ring, cremer_pople6, self_validation
    assert self_validation()["status"] == "passed"
    flat = cremer_pople6(analytic_ring(0, 0))
    assert flat["Q_A"] < 1e-12 and flat["theta_deg"] is None
    with pytest.raises(ValueError, match="Six finite"):
        cremer_pople6([[0., 0., 0.]])


def test_opposite_chairs_are_not_confused_by_pole_only_diagnostics():
    from aav_benchmark.ring_geometry import analytic_ring, compare_to_ideal, cremer_pople6
    north, south = analytic_ring(0, .6), analytic_ring(0, -.6)
    assert cremer_pople6(north)["chair_pole_distance_deg"] < 1e-10
    assert cremer_pople6(south)["chair_pole_distance_deg"] < 1e-10
    assert compare_to_ideal(north, south)["pucker_vector_angle_to_CCD_ideal_deg"] == pytest.approx(180.)
    assert compare_to_ideal(north, north)["pucker_vector_angle_to_CCD_ideal_deg"] == pytest.approx(0.)
