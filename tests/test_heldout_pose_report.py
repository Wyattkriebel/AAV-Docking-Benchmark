"""Prospective controls use source natives or synthetic poses, never outcomes."""
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np
import pytest

import aav_benchmark.heldout_pose_report as heldout
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import Evidence, accuracy_summary


ROOT = Path(__file__).resolve().parents[1]


def pin(root, path, value):
    write_json(root / path, value)
    return {"path": path, "sha256": sha256(root / path)}


def minimal_case():
    atom = {"canonical_atom_key": "1:GAL:C1", "xyz_A": [0, 0, 0]}
    case = {"inventory_entry": {"expected_heavy_atom_count": 2}, "scored_atom_count": 1,
            "scored_atom_map": [atom], "qualified_copies": [{"copy_id": "copyA", "scored_atom_count": 1, "scored_atom_map": [copy.deepcopy(atom)]}]}
    reference = {"names": ["L000", "L001"], "canonical_atom_map": [
        {"canonical_atom_key": "1:GAL:C1", "name": "L000"}, {"canonical_atom_key": "1:GAL:C2", "name": "L001"}]}
    return case, reference


def test_fixed_mask_requires_full_parent_and_same_mask_in_every_copy():
    case, reference = minimal_case()
    mapping = heldout.canonical_mapping(case, reference)
    assert mapping["indices"] == [0] and mapping["parent_atom_count"] == 2
    lost = copy.deepcopy(reference); lost["names"].pop(); lost["canonical_atom_map"].pop()
    with pytest.raises(ValueError, match="Full-parent atom loss"):
        heldout.canonical_mapping(case, lost)
    altered = copy.deepcopy(case)
    altered["qualified_copies"][0]["scored_atom_map"][0]["canonical_atom_key"] = "1:GAL:C2"
    with pytest.raises(ValueError, match="fixed canonical"):
        heldout.canonical_mapping(altered, reference)
    duplicated = copy.deepcopy(reference)
    duplicated["canonical_atom_map"][1]["canonical_atom_key"] = "1:GAL:C1"
    with pytest.raises(ValueError, match="bijection"):
        heldout.canonical_mapping(case, duplicated)


def test_site_any_qualified_copy_fixed_ranks_and_copy_box_coverage():
    copies = [{"copy_id": "A", "coordinates": np.asarray([[0., 0, 0]])},
              {"copy_id": "B", "coordinates": np.asarray([[20., 0, 0]])}]
    candidates = [{"rank": 1, "original_rank": 7, "center_A": [24., 0, 0]},
                  {"rank": 2, "original_rank": 19, "center_A": [0., 0, 0]}]
    value = heldout.site_summary(candidates, copies, box_sizes=[8, 8, 8])
    assert not value["1"]["success_3A"] and value["1"]["success_4A"] and value["1"]["success_5A"]
    assert value["1"]["qualified_copy_coverage_4A"] == .5
    assert value["5"]["qualified_copy_coverage_4A"] == 1
    assert value["1"]["qualified_copy_max_box_atom_coverage"] == [0, 1]
    assert value["candidates"][0]["candidate_original_rank"] == 7
    failed = heldout.site_summary([], copies)
    assert not failed["1"]["success_4A"] and failed["5"]["qualified_copy_coverage_4A"] == 0


def test_protein_fit_does_not_fit_ligand_and_added_atoms_never_become_truth():
    keys = [("A", 1, "", "ALA", name) for name in ("N", "CA", "C", "O")]
    protein = np.asarray([[0., 0, 0], [1., 0, 0], [0., 1., 0], [0., 0., 1.]])
    pmap = {"atoms": [{"output_key": list(key), "assembly_xyz_A": xyz.tolist()} for key, xyz in zip(keys, protein)]}
    mapping = {"indices": [0], "copies": [{"copy_id": "A", "coordinates": np.asarray([[0., 0., 2.]])}]}

    class Graph:
        scored_maps = [[0]]
        def diagnostics(self, full):
            if np.asarray(full).shape != (1, 3):
                raise ValueError("full parent atom loss")
            return {"physical_validity_certified": False}

    pose = {"rank": 1, "vina_score": -5., "candidate_class_rank": 1, "candidate_original_rank": 9,
            "box_id": "box01", "pooled_rank": 1, "coordinates": [[3., 4., 2.]]}
    scored, controls = heldout.score_pool([pose], mapping, keys, protein, pmap, Graph())
    assert scored[0]["global_RMSD_A"] == pytest.approx(5)
    assert scored[0]["local_interface_RMSD_A"] == pytest.approx(5)
    assert controls[0]["native_controls"]["native_RMSD_A"] == 0
    extended_keys, extended_xyz = keys + [("A", 1, "", "ALA", "CB")], np.vstack([protein, [2, 2, 2]])
    fit = heldout.protein_fits(extended_keys, extended_xyz, pmap, mapping["copies"][0]["coordinates"])
    assert fit["generated_protein_atom_count"] == 1 and len(fit["native_protein"]) == 4
    with pytest.raises(ValueError, match="lost an observed"):
        heldout.protein_fits(keys[:-1], protein[:-1], pmap, mapping["copies"][0]["coordinates"])
    bad = {**pose, "coordinates": []}
    with pytest.raises(ValueError, match="atom loss"):
        heldout.score_pool([bad], mapping, keys, protein, pmap, Graph())
    # First-chain controls may omit other input chains, but all observed native
    # assembly atoms stay in the reference/contact universe.
    other_key = ("B", 1, "", "ALA", "CA")
    full_map = {"atoms": pmap["atoms"] + [{"output_key": list(other_key), "assembly_xyz_A": [0, 0, 2]}]}
    fitted = heldout.protein_fits(keys, protein, full_map, mapping["copies"][0]["coordinates"], allow_missing_observed=True)
    assert len(fitted["native_protein"]) == 5 and fitted["missing_observed_protein_atom_count"] == 1
    native_pose = {**pose, "coordinates": [[0., 0., 2.]]}
    reduced, _ = heldout.score_pool([native_pose], mapping, keys, protein, full_map, Graph(), allow_missing_observed=True)
    assert reduced[0]["global_RMSD_A"] == 0
    assert reduced[0]["global_contacts"]["contact_recall"] == .8
    assert reduced[0]["full_observed_assembly_min_protein_distance_A"] == 0


def test_portable_source_archive_checks_members_and_rejects_tamper(tmp_path):
    content = b"fixed portable source\n"
    archive = tmp_path / "portable/source.tar.gz"; archive.parent.mkdir()
    with tarfile.open(archive, "w:gz") as tar:
        item = tarfile.TarInfo("case/source.json"); item.size = len(content)
        tar.addfile(item, io.BytesIO(content))
    manifest = {"archive": "source.tar.gz", "archive_sha256": sha256(archive),
                "files": [{"path": "case/source.json", "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}]}
    record = pin(tmp_path, "portable/manifest.json", manifest)
    assert heldout.archive_files(Evidence(tmp_path), record) == {"case/source.json": content}
    archive.write_bytes(archive.read_bytes() + b"tamper")
    with pytest.raises(ValueError, match="checksum"):
        heldout.archive_files(Evidence(tmp_path), record)


def test_failure_conditions_keep_all_six_cases_and_all18_methods(tmp_path):
    config_record = pin(tmp_path, "config.json", {"primary_case_ids": list(heldout.CASE_IDS), "methods": list(heldout.METHODS),
                                                   "seeds": [1, 2, 3], "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
                                                   "min_rmsd_A": .5, "energy_range_kcal_mol": 10})
    config = json.loads((tmp_path / "config.json").read_text())
    pin(tmp_path, "benchmarks/vina.lock.json", {"assets": {"test": {"sha256": "a" * 64}}})
    sources, descriptors = {}, []
    for cid in heldout.CASE_IDS:
        sources[f"{cid}/receptor.cif"] = b"frozen source protein"
        for method in heldout.METHODS:
            base = f"results/test/{cid}/{method}"
            condition = {"case_id": cid, "method": method, "configuration_sha256": config_record["sha256"],
                         "fixed_denominator_member": True, "jobs": [], "status": "preparation_failed"}
            descriptor = {"case_id": cid, "method": method, "preparation": base + "/prepared/condition.json", "execution_root": base + "/runs"}
            if method != "oracle":
                prediction = pin(tmp_path, base + "/predictor/execution.json", {"status": "failed", "predictions": []})
                selection = {"case_id": cid, "method": method, "configuration_sha256": config_record["sha256"], "native_coordinate_access": False,
                             "prediction_input_sha256": hashlib.sha256(sources[f"{cid}/receptor.cif"]).hexdigest(),
                             "prediction_record": {"filename": "execution.json", "sha256": prediction["sha256"]}, "predictor_status": "failed", "candidates": []}
                condition.update(status="predictor_failed", candidate_selection=selection)
                descriptor["prediction_selection"] = base + "/predictor/candidate_selection.json"
                pin(tmp_path, descriptor["prediction_selection"], selection)
            pin(tmp_path, descriptor["preparation"], condition)
            descriptors.append(descriptor)
    profile = {"configuration": config_record, "conditions": descriptors, "role": "strict_primary"}
    conditions, pending, counts = heldout.read_workflows(Evidence(tmp_path), profile, config, sources)
    assert len(conditions) == counts["eligible_conditions"] == 18
    assert counts["fixed_case_denominator"] == 6 and counts["conditions_without_prepared_boxes"] == 18
    assert counts["predictor_failed_conditions"] == 12 and counts["scheduled_seed_runs"] == 0 and pending == []
    assert not accuracy_summary([], "global_RMSD_A")["top1_success_2A"]
    # Add one synthetic prepared candidate with two terminal failures and a
    # missing third seed. The frozen denominator grows by exactly three runs;
    # it cannot be reduced to the two terminal records.
    first = descriptors[0]
    base = Path(first["preparation"]).parent / "box01"
    atom_map = [{"name": "L000", "canonical_atom_key": "1:GAL:C1"}]
    sources["2RDK/ligand_atom_map.json"] = json.dumps({"atoms": atom_map}).encode()
    sources["2RDK/chemistry.json"] = json.dumps({"box_size": [24, 24, 24]}).encode()
    prep = {"case_id": "2RDK", "configuration_sha256": config_record["sha256"], "box_center": [0, 0, 0], "box_size": [24, 24, 24],
            "neighborhood": {"interaction_cutoff_A": 8, "closest_omitted_atom_to_box_A": None}}
    prep_pin = pin(tmp_path, (base / "preparation.json").as_posix(), prep)
    pin(tmp_path, (base / "reference.json").as_posix(), {"names": ["L000"], "canonical_atom_map": atom_map, "pose_reference_eligible": False})
    for name in ("ligand.pdbqt", "neighborhood.pdbqt", "neighborhood.pdb"):
        (tmp_path / base / name).write_text("frozen synthetic input")
    prediction = pin(tmp_path, str(Path(first["prediction_selection"]).parent / "execution.json"), {"status": "completed", "predictions": []})
    selection = json.loads((tmp_path / first["prediction_selection"]).read_text())
    selection.update(predictor_status="completed", prediction_record={"filename": "execution.json", "sha256": prediction["sha256"]},
                     candidates=[{"rank": 1, "original_rank": 7, "center_A": [0, 0, 0]}])
    pin(tmp_path, first["prediction_selection"], selection)
    condition = json.loads((tmp_path / first["preparation"]).read_text())
    condition.update(status="prepared", candidate_selection=selection,
                     jobs=[{"id": "box01", "prepared_directory": "box01", "candidate_class_rank": 1, "candidate_original_rank": 7,
                            "center_A": [0, 0, 0], "status": "prepared", "preparation_sha256": prep_pin["sha256"]}])
    pin(tmp_path, first["preparation"], condition)
    run_dir = Path(first["execution_root"]) / "box01"
    (tmp_path / run_dir).mkdir(parents=True)
    (tmp_path / run_dir / "configuration.json").write_bytes((tmp_path / "config.json").read_bytes())
    fingerprint = {"config_sha256": config_record["sha256"], "preparation_sha256": prep_pin["sha256"], "full_assembly_sha256": None,
                   "reference_sha256": sha256(tmp_path / base / "reference.json"), "executable_sha256": "a" * 64,
                   "ligand_sha256": sha256(tmp_path / base / "ligand.pdbqt"), "receptor_sha256": sha256(tmp_path / base / "neighborhood.pdbqt"),
                   "receptor_pdb_sha256": sha256(tmp_path / base / "neighborhood.pdb")}
    runs = []
    for seed in (1, 2):
        pose, log = f"neighborhood_{seed}.pdbqt", f"neighborhood_{seed}.log"
        (tmp_path / run_dir / log).write_text("synthetic terminal execution failure\n")
        command = ["vina", "--receptor", str(tmp_path / base / "neighborhood.pdbqt"), "--ligand", str(tmp_path / base / "ligand.pdbqt"),
                   "--out", str(tmp_path / run_dir / pose), "--scoring", "vina", "--seed", str(seed), "--cpu", "2", "--exhaustiveness", "16",
                   "--num_modes", "20", "--min_rmsd", "0.5", "--energy_range", "10", "--center_x", "0", "--center_y", "0", "--center_z", "0",
                   "--size_x", "24", "--size_y", "24", "--size_z", "24"]
        runs.append({"seed": seed, "status": "failed", "fingerprint": fingerprint, "case_id": "2RDK", "receptor": "neighborhood",
                     "information_regime": "predicted-site-local", "vina_version": "AutoDock Vina v1.2.7", "pose_file": pose, "log_file": log,
                     "command": command, "wall_seconds": 2., "cpu_seconds": 3.})
    pin(tmp_path, (run_dir / "runs.json").as_posix(), runs)
    _, pending, counts = heldout.read_workflows(Evidence(tmp_path), profile, config, sources)
    assert counts["scheduled_seed_runs"] == 3 and counts["terminal_seed_runs"] == counts["failed_seed_runs"] == 2
    assert counts["measured_terminal_cpu_seconds"] == 6 and any("missing_seeds" in p for p in pending)
    # Missing evidence is pending, whereas documented failures are terminal zeros.
    (tmp_path / descriptors[0]["preparation"]).unlink()
    _, pending, _ = heldout.read_workflows(Evidence(tmp_path), profile, config, sources)
    assert pending and "missing_condition" in pending[0]


def test_pending_all_profiles_gate_precedes_any_native_fit_or_site_accuracy(tmp_path, monkeypatch):
    config = pin(tmp_path, "config.json", {})
    scientific_freeze = pin(tmp_path, "execution_freeze.json", {"configuration_sha256": config["sha256"], "frozen_before_new_docking": True})
    source_receipt = pin(tmp_path, "source.json", {})
    family = pin(tmp_path, "family.json", {})
    plan = heldout.example_plan()
    plan.pop("site_prediction")
    plan["family_audit"] = family
    plan["profiles"] = [{"id": "pending", "role": "strict_primary", "configuration": config, "freeze_receipt": scientific_freeze,
                         "source_archive_manifest": {}, "source_preparation_receipt": source_receipt}]
    plan_record = pin(tmp_path, "plan.json", plan)
    pin(tmp_path, "freeze.json", {"evaluation_plan_sha256": plan_record["sha256"], "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    monkeypatch.setattr(heldout, "check_reference_sources", lambda *a: {})
    monkeypatch.setattr(heldout, "archive_files", lambda *a: {})
    monkeypatch.setattr(heldout, "read_workflows", lambda *a, **kw: ([], ["case/seed:running"], {"eligible_conditions": 18}))
    monkeypatch.setattr(heldout, "_source_context", lambda *a: pytest.fail("native context entered before full terminal gate"))
    monkeypatch.setattr(heldout, "site_summary", lambda *a, **k: pytest.fail("site accuracy entered before full terminal gate"))
    with pytest.raises(ValueError, match="remains incomplete"):
        heldout.generate_report(tmp_path, "plan.json", "freeze.json")
    status = heldout.generate_report(tmp_path, "plan.json", "freeze.json", allow_pending=True)
    assert not status["metrics_eligible"] and "conditions" not in status and len(status["pending"]) == 1
    pin(tmp_path, "freeze.json", {"evaluation_plan_sha256": "0" * 64, "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    with pytest.raises(ValueError, match="pre-accuracy"):
        heldout.generate_report(tmp_path, "plan.json", "freeze.json")


def test_all_six_pinned_source_native_controls_without_outcome_reads(monkeypatch):
    evidence = Evidence(ROOT)
    record = {"path": "benchmarks/heldout_reference_ledger.json", "sha256": sha256(ROOT / "benchmarks/heldout_reference_ledger.json")}
    from aav_benchmark.paired_method_report import check_public_reference_sources
    cases, _ = check_public_reference_sources(evidence, record,
        {"path": "data/archive/heldout_reference/source_manifest.json", "sha256": sha256(ROOT / "data/archive/heldout_reference/source_manifest.json")},
        {"path": "benchmarks/source_public_distribution.json", "sha256": sha256(ROOT / "benchmarks/source_public_distribution.json")})
    source = heldout.archive_files(evidence, {"path": "results/structural_test/source_archive_manifest.json", "sha256": sha256(ROOT / "results/structural_test/source_archive_manifest.json")})
    # Only source preparation/native references are loaded; predictions and run artifacts are absent from this path.
    for cid, case in cases.items():
        mapping, graph, keys, protein, pmap = heldout._source_context(source, cid, case)
        assert len(mapping["indices"]) == (11 if cid == "4R9F" else mapping["parent_atom_count"])
        full_native = {f"{a['branch_num']}:{a['ccd_id']}:{a['atom_name']}": a["xyz_A"] for a in case["native_atom_map"]}
        atoms = json.loads(source[f"{cid}/ligand_atom_map.json"])["atoms"]
        native = np.asarray([full_native[a["canonical_atom_key"]] for a in atoms])
        row = {"rank": 1, "vina_score": -1., "candidate_class_rank": 1, "candidate_original_rank": 1,
               "box_id": "analytic_control", "pooled_rank": 1, "coordinates": native}
        scored, controls = heldout.score_pool([row], mapping, keys, protein, pmap, graph)
        assert scored[0]["global_RMSD_A"] < 1e-10 and scored[0]["local_interface_RMSD_A"] < 1e-10
        assert controls[0]["native_controls"]["native_RMSD_A"] == 0
        assert controls[0]["native_controls"]["native_contact_recall"] == 1
        assert scored[0]["global_contacts"]["contact_recall"] == scored[0]["global_contacts"]["contact_precision"] == 1
        assert controls[0]["native_controls"]["known_translation_RMSD_A"] == 5
        shifted, _ = heldout.score_pool([{**row, "coordinates": native + [3, 4, 0]}], mapping, keys, protein, pmap, graph)
        assert shifted[0]["global_RMSD_A"] == pytest.approx(5)


def test_family_weighting_and_conditional_bootstrap_preserve_zero_failures(monkeypatch):
    rows = []
    for index, cid in enumerate(heldout.CASE_IDS):
        for method in heldout.METHODS:
            # Every oracle succeeds; every FP/P2 unsupported pool contributes0.
            success = method == "oracle"
            summary = accuracy_summary([{"distance": 0.}] if success else [], "distance")
            budget = {"global": summary, "local_interface": summary, "status": "ok" if success else "unsupported"}
            rows.append({"case_id": cid, "method": method, "atomic_scope": "modeled_core" if cid == "4R9F" else "full_ligand",
                         "seeds": [{"seed": seed, "site_budgets": {"1": budget, "5": budget}} for seed in (1, 2, 3)],
                         "site": {str(b): {"available_candidates": int(success), **{f"success_{t}A": success for t in (3, 4, 5)}} for b in (1, 5)}})
    families = {**heldout.FAMILIES, "3ZWE": heldout.FAMILIES["2RDK"]}  # explicit conservative merge in a frozen plan
    summary = heldout.macro_summaries(rows, families)
    oracle = next(s for s in summary if s["atomic_scope"] == "full_ligand" and s["method"] == "oracle")
    fpocket = next(s for s in summary if s["atomic_scope"] == "full_ligand" and s["method"] == "fpocket")
    assert oracle["eligible_cases"] == 5 and oracle["source_family_clusters"] == 4
    assert oracle["family_macro"]["top1_success_2A"] == 1 and fpocket["family_macro"]["top1_success_2A"] == 0
    original = heldout.paired_family_bootstrap
    cached, calls = {}, []
    def check_bootstrap(records, a, b, **kwargs):
        assert kwargs["n_resamples"] == 10000 and kwargs["seed"] == 20261103
        assert all(not r["success"] for r in records if r["status"] == "unsupported")
        calls.append(1)
        key = (len(records), a, b)
        if key not in cached:
            cached[key] = original(records, a, b, **kwargs)
        return copy.deepcopy(cached[key])
    monkeypatch.setattr(heldout, "paired_family_bootstrap", check_bootstrap)
    paired = heldout.conditional_statistics(rows, families, [1, 2, 3])
    assert len(paired) == len(calls) == 144
    assert all(p["difference_a_minus_b"] == 1 for p in paired)
    assert all(p["confidence_interval"] is None for p in paired if p["atomic_scope"] == "modeled_core")
    assert all("Conditional" in p["interpretation"] for p in paired)
