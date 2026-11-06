"""Qualified maps, receptor-only fits, rank pooling and strict portable evidence."""
import copy
import gzip
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import (
    Evidence, ParentGeometry, _protein_alignment, _score_pool, _validate_command,
    accuracy_summary, aggregate_conditions, generate_qualified_pose_report,
    native_metric_controls, pooled_poses, qualified_atom_map, write_qualified_pose_report,
)


def _map():
    reference = {"names": ["L000", "L001", "L002"], "coordinates": [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                 "atom_map": [{"name": f"L{i:03d}", "source_key": ["B", 1, " ", "GAL", name], "rdkit_index": i} for i, name in enumerate(["C1", "C2", "O1"])]}
    record = {"scored_atom_count": 2, "scored_atom_map": [
        {"prepared_chain": "B", "prepared_resnum": 1, "prepared_resname": "GAL", "atom_name": name, "xyz_A": xyz}
        for name, xyz in [("O1", [0, 1, 0]), ("C1", [0, 0, 0])]]}
    return record, reference


def test_qualified_core_keeps_fixed_named_atoms_and_never_silently_scores_parent_tail():
    record, reference = _map()
    result = qualified_atom_map(record, reference)
    assert result["indices"] == [2, 0] and result["names"] == ["L002", "L000"]
    assert result["parent_atom_count"] == 3
    for mutation in ("lost", "duplicated", "coordinate"):
        bad = copy.deepcopy(reference)
        if mutation == "lost":
            bad["atom_map"].pop()
        elif mutation == "duplicated":
            bad["atom_map"][2]["source_key"] = bad["atom_map"][0]["source_key"]
        else:
            bad["coordinates"][2][0] = .01
        with pytest.raises(ValueError):
            qualified_atom_map(record, bad)


def test_fixed_pool_order_ties_and_nested_site_budgets():
    def pose(rank, score):
        return {"rank": rank, "vina_score": score, "coordinates": [[0, 0, 0]]}
    boxes = [{"id": "box02", "candidate_class_rank": 2, "candidate_original_rank": 5, "poses_by_seed": {1: [pose(1, -5), pose(2, -4)]}},
             {"id": "box01", "candidate_class_rank": 1, "candidate_original_rank": 1, "poses_by_seed": {1: [pose(1, -5), pose(2, -3)]}}]
    pooled = pooled_poses(boxes, 1, 5)
    assert [(x["box_id"], x["rank"]) for x in pooled] == [("box01", 1), ("box02", 1), ("box02", 2), ("box01", 2)]
    assert len(pooled_poses(boxes, 1, 1)) == 2
    assert pooled_poses(boxes, 2, 5) == []


def test_failed_eligible_pool_is_zero_and_best_returned_is_not_ranked_success():
    empty = accuracy_summary([], "RMSD_A")
    assert empty["top1_success_2A"] is False and empty["top5_success_2A"] is False
    assert empty["best_returned_RMSD_A"] is None
    result = accuracy_summary([{"RMSD_A": v} for v in [8, 7, 6, 5, 4, .5]], "RMSD_A")
    assert result["top5_success_2A"] is False and result["best_returned_success_2A"] is True
    assert result["best_returned_pooled_rank"] == 6
    rows = [{"atomic_scope": "modeled_core", "method": "p2rank", "complex_id": case,
             "site_budgets": {"5": {"global": {field: float(value) for field, value in metric.items() if "success_" in field}}}}
            for case, metric in [("one", empty), ("two", result)]]
    aggregate = aggregate_conditions(rows, "modeled_core", "p2rank", 5, "global")
    assert aggregate["eligible_complexes"] == 2 and aggregate["top5_success_2A"] == 0
    assert aggregate["best_returned_success_2A"] == .5
    assert aggregate_conditions(rows, "full_ligand", "p2rank", 5, "global")["top5_success_2A"] is None


def test_native_zero_and_known_translation_controls():
    native = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0]])
    protein = np.array([[0., 0, 2], [1, 2, 0]])
    controls = native_metric_controls(native, protein)
    assert controls["native_RMSD_A"] == 0 and controls["known_translation_RMSD_A"] == 5
    assert controls["native_contact_recall"] == controls["native_contact_precision"] == 1


def test_backbone_fit_does_not_fit_a_wrong_ligand_independently():
    keys = [("A", 1, "", "ALA", name) for name in ["N", "CA", "C", "O"]]
    target = np.array([[0., 0, 2], [1, 0, 2], [0, 1, 2], [1, 1, 2]])
    source = target + [10, 0, 0]
    native = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0]])
    fits, contacts, _ = _protein_alignment(keys, source, keys, target, native)
    class Graph:
        scored_maps = [[0, 1, 2]]
        def diagnostics(self, pose):
            return {"not_a_validity_certificate": True}
    pose = {"rank": 1, "vina_score": -3., "candidate_class_rank": 1, "candidate_original_rank": 1,
            "box_id": "box01", "pooled_rank": 1, "coordinates": native + [13, 0, 0]}
    result = _score_pool([pose], native, source, target, keys, fits, contacts, Graph())[0]
    assert result["global_RMSD_A"] == pytest.approx(3)
    assert result["local_interface_RMSD_A"] == pytest.approx(3)


def test_archives_validate_compressed_and_uncompressed_bytes(tmp_path):
    directory = tmp_path / "case/inputs"
    directory.mkdir(parents=True)
    source = b"source atom map"
    compressed = gzip.compress(source, mtime=0)
    (directory / "map.json.gz").write_bytes(compressed)
    import hashlib
    write_json(directory / "manifest.json", {"files": [{"filename": "map.json", "archive": "map.json.gz",
        "sha256": hashlib.sha256(source).hexdigest(), "archive_sha256": hashlib.sha256(compressed).hexdigest(), "original_bytes": len(source)}]})
    assert Evidence(tmp_path).bundle("case")["map.json"] == source
    (directory / "map.json.gz").write_bytes(compressed + b"changed")
    with pytest.raises(ValueError, match="checksum"):
        Evidence(tmp_path).bundle("case")


def test_graph_stereo_and_rigid_ring_diagnostics_do_not_claim_physical_validity():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(=O)O"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=13) == 0
    sdf = Chem.MolToMolBlock(molecule).encode() + b"\n$$$$\n"
    heavy = Chem.RemoveHs(molecule)
    n = heavy.GetNumAtoms()
    reference = {"names": [f"L{i:03d}" for i in range(n)],
        "atom_map": [{"name": f"L{i:03d}", "rdkit_index": i, "source_key": ["B", 1, "", "LIG", str(i)]} for i in range(n)]}
    graph = ParentGeometry(sdf, reference, {"indices": list(range(n))})
    xyz = np.asarray(heavy.GetConformer().GetPositions())
    good = graph.diagnostics(xyz)
    assert good["stereochemistry_from_pose_3D_matches_input"] is True
    assert good["ring_pucker_validity_certified"] is False
    reflected = xyz.copy()
    reflected[:, 0] *= -1
    assert graph.diagnostics(reflected)["stereochemistry_from_pose_3D_matches_input"] is False
    with pytest.raises(ValueError, match="atom loss"):
        graph.diagnostics(xyz[:-1])


def _pending_fixture(root):
    (root / "benchmarks").mkdir()
    (root / "results/blind_pose").mkdir(parents=True)
    write_json(root / "benchmarks/blind_pose_experiment.json", {"frozen_before_new_docking": True,
        "atomic_reference_ledger": "benchmarks/ledger.json", "atomic_reference_ledger_sha256": "placeholder"})
    write_json(root / "benchmarks/ledger.json", {"source_only": True})
    write_json(root / "results/blind_pose/planning_inventory.json", {})
    write_json(root / "results/blind_pose/oracle_controls.json", {})
    config = json.loads((root / "benchmarks/blind_pose_experiment.json").read_text())
    config["atomic_reference_ledger_sha256"] = sha256(root / "benchmarks/ledger.json")
    write_json(root / "benchmarks/blind_pose_experiment.json", config)
    write_json(root / "results/blind_pose/freeze_receipt.json", {"frozen_before_new_docking": True,
        "new_pose_outcomes_consulted": False, "configuration_sha256": sha256(root / "benchmarks/blind_pose_experiment.json"),
        "atomic_reference_ledger_sha256": sha256(root / "benchmarks/ledger.json"),
        "planning_inventory_sha256": sha256(root / "results/blind_pose/planning_inventory.json"),
        "oracle_controls_sha256": sha256(root / "results/blind_pose/oracle_controls.json")})


def test_pending_workflow_never_computes_or_publishes_partial_atomic_accuracy(tmp_path):
    _pending_fixture(tmp_path)
    with patch("aav_benchmark.qualified_pose._check_atomic_sources"), patch(
            "aav_benchmark.qualified_pose._all_execution", return_value=([], ["case/box:running"], {"scheduled_seed_runs": 156})), patch(
            "aav_benchmark.qualified_pose._native_protein", side_effect=AssertionError("Atomic evaluation before complete run denominator")):
        with pytest.raises(ValueError, match="incomplete"):
            generate_qualified_pose_report(tmp_path)
        summary, text = generate_qualified_pose_report(tmp_path, allow_pending=True)
        assert summary["metrics_eligible"] is False and "conditions" not in summary
        assert "No atomic accuracy aggregate" in text
    (tmp_path / "benchmarks/blind_pose_experiment.json").write_text("{}")
    with pytest.raises(ValueError, match="freeze"):
        generate_qualified_pose_report(tmp_path)


def test_portable_check_tolerates_only_explicit_derived_geometry(tmp_path):
    summary = {"global_RMSD_A": 2.123456789, "vina_score": -5., "success": True,
               "count": 156, "source_sha256": "unchanged", "native": {"coordinates": [[1., 2., 3.]]}}
    with patch("aav_benchmark.qualified_pose.generate_qualified_pose_report", return_value=(summary, "Complete report\n")):
        paths = write_qualified_pose_report(tmp_path)
        rounded = copy.deepcopy(summary)
        rounded["global_RMSD_A"] += 1e-12
        write_json(paths[0], rounded)
        with pytest.raises(ValueError, match="stale"):
            write_qualified_pose_report(tmp_path, check=True)
        write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)
        for field, value in [("vina_score", -5. + 1e-12), ("count", 155), ("success", False),
                             ("source_sha256", "changed"), ("global_RMSD_A", 2.2)]:
            changed = copy.deepcopy(summary)
            changed[field] = value
            write_json(paths[0], changed)
            with pytest.raises(ValueError):
                write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)
        changed = copy.deepcopy(summary)
        changed["native"]["coordinates"][0][0] += 1e-12
        write_json(paths[0], changed)
        with pytest.raises(ValueError):
            write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)
        write_json(paths[0], summary)
        paths[1].write_text("Changed prose\n")
        with pytest.raises(ValueError, match="stale"):
            write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)


def test_native_coordinate_pools_recover_each_accepted_scope_without_promoting_rejected_cases(tmp_path):
    """Analytic native pools, not docking results; raw-run integrity is tested separately."""
    from aav_benchmark.qualified_pose import _old_inputs
    from aav_benchmark.method_protocol_sources import prepare_source_fixture
    prepare_source_fixture(tmp_path, groups=("development",))
    root = tmp_path
    evidence = Evidence(root)
    config = json.loads((root / "benchmarks/blind_pose_experiment.json").read_text())
    ledger = json.loads((root / "benchmarks/atomic_reference_ledger.json").read_text())
    records = {r["case_id"]: r for r in ledger["cases"]}
    conditions = []
    for case in config["cases"]:
        cid = case["complex_id"]
        reference = None
        if cid in records and records[cid]["atomic_pose_eligible"]:
            inputs = _old_inputs(evidence, cid)
            if inputs:
                reference = json.loads(inputs["reference.json"])
        for method in config["methods"]:
            boxes = []
            if reference:
                boxes = [{"id": "box01", "candidate_class_rank": 1, "candidate_original_rank": 1,
                          "poses_by_seed": {seed: [{"rank": 1, "vina_score": -5., "coordinates": np.asarray(reference["coordinates"])}] for seed in config["seeds"]}}]
            conditions.append({"condition_id": case["condition_id"], "complex_id": cid, "state": case["state"],
                               "method": method, "preparation_status": "synthetic_native_control", "boxes": boxes})
    denominator = {"eligible_conditions": 27, "prepared_boxes": 52, "scheduled_seed_runs": 156,
                   "terminal_seed_runs": 156, "successful_seed_runs": 156, "failed_seed_runs": 0}
    with patch("aav_benchmark.qualified_pose._all_execution", return_value=(conditions, [], denominator)):
        summary, markdown = generate_qualified_pose_report(root)
    assert summary["denominator"]["atomic_eligible_method_conditions"] == 18
    for method in config["methods"]:
        for budget in ("1", "5"):
            for fit in ("global", "local_interface"):
                assert summary["aggregates"]["full_ligand"][method][budget][fit]["top1_success_2A"] == 1
                # Three supported cores succeed; unsupported full parent1QFO
                # remains zero in the fixed four-core denominator.
                assert summary["aggregates"]["modeled_core"][method][budget][fit]["top1_success_2A"] == .75
    excluded = [row for row in summary["conditions"] if not row["atomic_pose_eligible"]]
    assert len(excluded) == 9 and all(row["atomic_accuracy"] is None for row in excluded)
    assert "No core result is called whole-ligand recovery" in markdown
