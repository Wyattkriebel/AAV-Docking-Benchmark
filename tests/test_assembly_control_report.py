"""Assembly-control provenance and synthetic metric controls; no actual accuracy."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

import aav_benchmark.assembly_control_report as assembly
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import accuracy_summary
from aav_benchmark.qualified_pose import Evidence
from aav_benchmark.heldout_pose_report import score_pool


ROOT = Path(__file__).resolve().parents[1]


def test_control_receipt_view_aliases_only_designated_receipt_and_preserves_raw_hash(tmp_path):
    path = Path("results/control/vina/execution.json")
    raw = {"prepared_conditions_sha256": "a" * 64, "condition_id": "control", "status": "complete"}
    write_json(tmp_path / path, raw)
    original = (tmp_path / path).read_bytes()
    evidence = Evidence(tmp_path)
    view = assembly.ControlReceiptView(evidence, path, "a" * 64)
    normalized = view.json(path)
    assert normalized["frozen_input_manifest_sha256"] == "a" * 64
    assert normalized["prepared_conditions_sha256"] == "a" * 64
    assert (tmp_path / path).read_bytes() == original and json.loads(original) == raw
    assert evidence.sources[path.as_posix()]["sha256"] == sha256(tmp_path / path)
    other = Path("results/control/unrelated.json"); write_json(tmp_path / other, raw)
    assert "frozen_input_manifest_sha256" not in view.json(other)
    with pytest.raises(ValueError, match="prepared-condition index"):
        assembly.ControlReceiptView(evidence, path, "b" * 64).json(path)
    write_json(tmp_path / path, {**raw, "frozen_input_manifest_sha256": "a" * 64})
    with pytest.raises(ValueError, match="duplicate"):
        view.json(path)


def center_conditions():
    rows = []
    for profile in ("strict", "repair", "strict_monomer", "repair_monomer"):
        for cid in assembly.CASE_IDS:
            for method in assembly.METHODS:
                center = [0., 0., 0.] if "monomer" not in profile or method == "oracle" else [10., 0., 0.]
                rows.append({"profile": profile, "case_id": cid, "method": method,
                             "candidates": [{"rank": 1, "original_rank": 3, "center_A": center}]})
    return rows


def test_center_policy_allows_end_to_end_changes_but_keeps_controlled_pairs_exact():
    rows = center_conditions()
    comparisons = assembly.validate_center_policy(rows, assembly.CENTER_POLICY)
    assert len(comparisons) == 6
    assert all(c["exact_candidate_center_and_rank_match"] for c in comparisons if c["method"] == "oracle")
    assert all(not c["exact_candidate_center_and_rank_match"] for c in comparisons if c["method"] != "oracle")
    assert all("end_to_end" in c["interpretation"] for c in comparisons if c["method"] != "oracle")
    changed = copy.deepcopy(rows)
    next(c for c in changed if c["profile"] == "repair_monomer" and c["method"] == "fpocket")["candidates"][0]["center_A"] = [9., 0., 0.]
    with pytest.raises(ValueError, match="Repair sensitivity"):
        assembly.validate_center_policy(changed, assembly.CENTER_POLICY)
    changed = copy.deepcopy(rows)
    for c in changed:
        if c["profile"].endswith("monomer") and c["method"] == "oracle":
            c["candidates"][0]["center_A"] = [9., 0., 0.]
    with pytest.raises(ValueError, match="Oracle assembly"):
        assembly.validate_center_policy(changed, assembly.CENTER_POLICY)


def source_pairs():
    observed = {}
    sources = {p: {} for p in ("strict", "repair", "strict_monomer", "repair_monomer")}
    for cid in assembly.CASE_IDS:
        atoms = [{"output_key": [chain, 1, "", "ALA", "CA"], "assembly_xyz_A": [float(i), 0, 0]}
                 for i, chain in enumerate(("A", "B"))]
        raw_map = json.dumps({"atoms": atoms}).encode()
        observed[f"{cid}/protein_atom_map.json"] = raw_map
        observed[f"{cid}/receptor.cif"] = b"original full observed assembly"
        for profile, source in sources.items():
            for name in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json"):
                source[f"{cid}/{name}"] = b"same chemical parent and independent initializer"
            if "monomer" in profile:
                source[f"{cid}/protein_atom_map.json"] = json.dumps({"atoms": atoms[:1]}).encode()
                source[f"{cid}/full_protein_atom_map.json"] = raw_map
                source[f"{cid}/full_assembly_receptor.cif"] = observed[f"{cid}/receptor.cif"]
    return sources, observed


def test_full_observed_truth_and_fixed_parent_cannot_be_replaced_by_monomer_or_model():
    sources, observed = source_pairs()
    assembly.validate_source_pairs(sources, observed)
    changed = copy.deepcopy(sources)
    changed["repair_monomer"]["3ZWE/full_protein_atom_map.json"] = changed["repair_monomer"]["3ZWE/protein_atom_map.json"]
    with pytest.raises(ValueError, match="full observed assembly truth"):
        assembly.validate_source_pairs(changed, observed)
    changed = copy.deepcopy(sources)
    changed["repair_monomer"]["3ZWE/ligand.sdf"] = b"different initializer"
    with pytest.raises(ValueError, match="initializer"):
        assembly.validate_source_pairs(changed, observed)
    changed = copy.deepcopy(sources)
    altered = json.loads(changed["strict_monomer"]["3ZWE/protein_atom_map.json"])
    altered["atoms"][0]["assembly_xyz_A"] = [99, 0, 0]
    changed["strict_monomer"]["3ZWE/protein_atom_map.json"] = json.dumps(altered).encode()
    with pytest.raises(ValueError, match="source-A subset"):
        assembly.validate_source_pairs(changed, observed)


def test_synthetic_omitted_chain_contacts_and_full_assembly_clashes_remain_in_truth():
    keys = [("A", 1, "", "ALA", name) for name in ("N", "CA", "C", "O")]
    protein = np.asarray([[0., 0, 0], [1., 0, 0], [0., 1., 0], [0., 0., 1.]])
    full_keys = keys + [("B", 1, "", "ALA", name) for name in ("N", "CA", "C", "O")]
    full_protein = np.vstack([protein, protein + [0, 0, 1]])
    protein_map = {"atoms": [{"output_key": list(k), "assembly_xyz_A": xyz.tolist()} for k, xyz in zip(full_keys, full_protein)]}
    mapping = {"indices": [0], "copies": [{"copy_id": "native", "coordinates": np.asarray([[0., 0., 2.]])}]}
    class Graph:
        scored_maps = [[0]]
        def diagnostics(self, full):
            if np.asarray(full).shape != (1, 3):
                raise ValueError("complete parent atom loss")
            return {"complete_physical_validity_certified": False}
    pose = {"rank": 1, "vina_score": -5., "candidate_class_rank": 1, "candidate_original_rank": 1,
            "box_id": "analytic", "pooled_rank": 1, "coordinates": [[0., 0., 2.]]}
    rows, controls = score_pool([pose], mapping, keys, protein, protein_map, Graph(), allow_missing_observed=True)
    assert rows[0]["global_RMSD_A"] == rows[0]["local_interface_RMSD_A"] == 0
    assert rows[0]["global_contacts"]["contact_recall"] == .5
    assert rows[0]["full_parent_min_protein_distance_A"] == 1
    assert rows[0]["full_observed_assembly_min_protein_distance_A"] == 0
    assert controls[0]["missing_observed_input_protein_atoms"] == 4
    shifted, _ = score_pool([{**pose, "coordinates": [[3., 4., 2.]]}], mapping, keys, protein, protein_map, Graph(), allow_missing_observed=True)
    assert shifted[0]["global_RMSD_A"] == pytest.approx(5)
    with pytest.raises(ValueError, match="atom loss"):
        score_pool([{**pose, "coordinates": []}], mapping, keys, protein, protein_map, Graph(), allow_missing_observed=True)


def test_conditional_statistics_keep_all_failed_cases_and_predeclared_resampling(monkeypatch):
    rows = []
    for profile in ("strict", "repair", "strict_monomer", "repair_monomer"):
        for cid in assembly.CASE_IDS:
            for method in assembly.METHODS:
                success = "monomer" not in profile
                summary = accuracy_summary([{"rmsd": 0.}] if success else [], "rmsd")
                rows.append({"profile": profile, "case_id": cid, "method": method,
                    "seeds": [{"seed": seed, "poses": [{}] if success else [], "global": summary, "local_interface": summary} for seed in assembly.SEEDS]})
    actual = assembly.paired_family_bootstrap
    cache = {}
    def bootstrap(records, a, b, **kwargs):
        assert kwargs == {"n_resamples": 10000, "seed": 20261103, "expected_seeds": list(assembly.SEEDS)}
        assert len(records) == 12
        if (a, b) not in cache:
            cache[a, b] = actual(records, a, b, **kwargs)
        return copy.deepcopy(cache[a, b])
    monkeypatch.setattr(assembly, "paired_family_bootstrap", bootstrap)
    results = assembly.comparison_statistics(rows, {cid: cid for cid in assembly.CASE_IDS})
    assert len(results) == 108 and all(r["difference_a_minus_b"] == -1 for r in results)
    assert all(r["n_complexes"] == 2 and r["n_families"] == 2 for r in results)
    assert all("end_to_end" in r["pairing_center_policy"] for r in results if r["site_condition"] != "oracle")


def _public_full_fixture(root):
    spec = importlib.util.spec_from_file_location("paired_public_fixture", ROOT / "tests/test_paired_method_public_portability.py")
    fixture = importlib.util.module_from_spec(spec); spec.loader.exec_module(fixture)
    fixture._copy_source_tree(root)
    experiment = json.loads((ROOT / "benchmarks/structural_test_vina_carb.json").read_text())
    fixture._evaluation_plan(root, experiment)


def _copy_assembly_sources(root):
    """The shared verified archive fixture already contains all12 controls."""
    index = json.loads((root / "results/structural_assembly_pose_v3/prepared_conditions.json").read_text())
    assert len(index["conditions"]) == 12
    for row in index["conditions"]:
        path = root / row["condition_path"]
        assert sha256(path) == row["condition_sha256"]


def test_public_archive_only_pending_gate_requires_all72_control_and_comparator_outcomes(tmp_path, monkeypatch):
    _public_full_fixture(tmp_path)
    _copy_assembly_sources(tmp_path)
    plan = assembly.example_plan(tmp_path)
    plan_path, freeze_path = "benchmarks/assembly_control_evaluation.json", "results/assembly_control_evaluation_freeze.json"
    write_json(tmp_path / plan_path, plan)
    write_json(tmp_path / freeze_path, {"evaluation_plan_sha256": sha256(tmp_path / plan_path),
        "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    def forbidden(*args, **kwargs): pytest.fail("Actual native accuracy reached before all control/comparator outcomes terminal")
    monkeypatch.setattr(assembly, "_source_context", forbidden)
    monkeypatch.setattr(assembly, "score_pool", forbidden)
    monkeypatch.setattr(assembly, "site_summary", forbidden)
    report = assembly.generate_report(tmp_path, plan_path, freeze_path, allow_pending=True)
    assert report["status"] == "pending_incomplete_workflow" and not report["metrics_eligible"]
    assert report["counts"]["eligible_seed_outcomes"] == 72 and report["counts"]["terminal_seed_outcomes"] == 0
    assert len(report["pending"]) == 24 and "conditions" not in report and "paired_family_statistics" not in report
    assert report["reference_validation"]["retained_archive_records_verified"] == 163
    assert not (tmp_path / "tools").exists() and not (tmp_path / "data/processed").exists()
    assert not any(path.startswith(("tools/", "data/processed/")) for path in report["source_sha256_inventory"])
    with pytest.raises(ValueError, match="remains incomplete"):
        assembly.generate_report(tmp_path, plan_path, freeze_path)
