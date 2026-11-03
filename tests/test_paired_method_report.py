"""Synthetic paired-reader controls; no native tool or actual accuracy runs."""
import copy
import hashlib
import json
import io
from pathlib import Path
import tarfile

import numpy as np
import pytest

import aav_benchmark.paired_method_report as paired
from aav_benchmark.additional_method import export_vina_carb_ligand, expected_chi_links
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import Evidence, accuracy_summary


def pin(root, path, value, *, raw=False):
    target = root / path; target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(value) if raw else write_json(target, value)
    return {"path": path, "sha256": sha256(target)}


def atom(serial, name, residue, number, xyz, kind="C"):
    return (f"ATOM  {serial:5d} {name:>4s} {residue:>3s} B{number:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00    {0.:+6.3f} {kind}")


def synthetic_vc(tmp_path):
    # Two synthetic complete named pyranose atom sets, linked C1-to-O4. This
    # tests alias/serial/provenance mechanics, not any experimental structure.
    lines, mapping = ["ROOT"], []
    names = ("C1", "C2", "C3", "C4", "C5", "O5", "O4")
    for residue in (1, 2):
        if residue == 2:
            lines.extend(["ENDROOT", "BRANCH 1 14"])
        for local, name in enumerate(names):
            serial = (residue - 1) * 7 + local + 1
            generated = f"L{serial - 1:03d}"
            lines.append(atom(serial, generated, "UNL", 1, [float(serial), 0., 0.], "OA" if name.startswith("O") else "C"))
            mapping.append({"name": generated, "source_key": ["B", residue, "", "MAN", name]})
    lines.extend(["ENDBRANCH 1 14", "TORSDOF 1"])
    source = pin(tmp_path, "inputs/ligand.pdbqt", "\n".join(lines) + "\n", raw=True)
    canonical = pin(tmp_path, "inputs/map.json", {"atoms": mapping})
    receptor = pin(tmp_path, "inputs/receptor.pdbqt", "synthetic paired receptor\n", raw=True)
    row = {"condition_id": "strict/synthetic/fpocket/box01", "case_id": "synthetic", "receptor_profile": "strict", "method": "fpocket",
           "status": "prepared", "ligand": source, "canonical_atom_map": canonical, "receptor": receptor,
           "box": {"center": [0., 0., 0.], "size": [24., 24., 24.]}}
    directory = Path("outputs/s104729"); absolute = tmp_path / directory; absolute.mkdir(parents=True)
    export = export_vina_carb_ligand(tmp_path / source["path"], mapping, absolute / "ligand.pdbqt")
    (absolute / "receptor.pdbqt").write_bytes((tmp_path / receptor["path"]).read_bytes())
    (absolute / "configuration.txt").write_text(paired.render_configuration(row, paired.SEEDS[0]))
    (absolute / "poses.pdbqt").write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\n" + (absolute / "ligand.pdbqt").read_text() + "ENDMDL\n")
    links = expected_chi_links(export)
    (absolute / "VC_log.txt").write_text("".join(f"CHI energy penalties applied to {a[0]} {a[1]} and {b[0]} {b[1]} linkage.\n" for a, b in links))
    (absolute / "search.log").write_text("mode | affinity | chi | affinity | dist from best mode\n1 -5.000 0.000 -5.000 0.000 0.000\n")
    for name in ("stdout.log", "stderr.log"):
        (absolute / name).write_text("")
    experiment = {"adapter_source": {"sha256": "a" * 64}}
    experiment_pin, receipt_pin = {"sha256": "b" * 64}, {"sha256": "c" * 64}
    executable = "d" * 64
    record = {"condition": row, "seed": paired.SEEDS[0], "configuration_sha256": experiment_pin["sha256"],
              "freeze_receipt_sha256": receipt_pin["sha256"], "driver_source_sha256": experiment["adapter_source"]["sha256"],
              "atomic_accuracy_computed": False, "experimental_ligand_coordinates_read": False, "status": "ok", "full_search_started": True,
              "parameters": paired.PARAMETERS, "executable_sha256": executable, "exit_code": 0,
              "command": ["/unexecuted/synthetic/vina-carb", "--config", "configuration.txt"], "exported_ligand_sha256": export["ligand_sha256"],
              "receptor_sha256": receptor["sha256"], "returned_poses": 1, "generated_names": [a["name"] for a in mapping],
              "cpu_seconds": 1.5, "wall_seconds": 1., "preflight": {"glycosidic_links": len(links), "fully_supported_links": len(links)}}
    def publish():
        record["artifacts"] = {p.name: sha256(p) for p in absolute.iterdir() if p.name != "execution.json"}
        write_json(absolute / "execution.json", record)
    publish()
    return row, directory, export, experiment, experiment_pin, receipt_pin, executable, record, publish


def test_vc_aliases_restore_canonical_order_and_atom_loss_fails(tmp_path):
    row, directory, export, *rest = synthetic_vc(tmp_path)
    path = tmp_path / directory / "poses.pdbqt"
    names = [a["generated_name"] for a in export["atoms"]][::-1]
    poses = paired.canonical_carb_poses(path, export, names)
    assert poses[0]["coordinates"][0].tolist() == [14., 0., 0.]
    assert poses[0]["coordinates"][-1].tolist() == [1., 0., 0.]
    assert poses[0]["method_score"] == -5
    path.write_text("\n".join(line for line in path.read_text().splitlines() if not line.startswith("ATOM     14")) + "\n")
    with pytest.raises(ValueError, match="Incomplete"):
        paired.canonical_carb_poses(path, export, names)


def test_successful_vc_reader_replays_export_checks_chi_and_raw_score_ranks(tmp_path):
    row, directory, export, experiment, plan_pin, freeze_pin, executable, record, publish = synthetic_vc(tmp_path)
    names = record["generated_names"]
    result, pending = paired.read_vina_carb(Evidence(tmp_path), directory, row, paired.SEEDS[0], experiment, plan_pin, freeze_pin, names, executable_sha256=executable)
    assert pending is None and result["status"] == "ok" and len(result["poses"]) == 1
    assert result["poses"][0]["method_score"] == -5
    (tmp_path / directory / "VC_log.txt").write_text("CHI NOT applied\n"); publish()
    with pytest.raises(ValueError, match="activate CHI"):
        paired.read_vina_carb(Evidence(tmp_path), directory, row, paired.SEEDS[0], experiment, plan_pin, freeze_pin, names, executable_sha256=executable)
    (tmp_path / directory / "VC_log.txt").write_text("CHI energy penalties applied to MAN 1 and MAN 2 linkage.\n")
    (tmp_path / directory / "search.log").write_text("mode | affinity | chi | affinity | dist from best mode\n1 -6.000 0.000 -6.000 0.000 0.000\n"); publish()
    with pytest.raises(ValueError, match="score table"):
        paired.read_vina_carb(Evidence(tmp_path), directory, row, paired.SEEDS[0], experiment, plan_pin, freeze_pin, names, executable_sha256=executable)


def test_vc_hash_updates_cannot_hide_export_changes_or_wrong_paired_settings(tmp_path):
    row, directory, export, experiment, plan_pin, freeze_pin, executable, record, publish = synthetic_vc(tmp_path)
    names = record["generated_names"]
    ligand = tmp_path / directory / "ligand.pdbqt"
    lines = ligand.read_text().splitlines()
    i = next(i for i, line in enumerate(lines) if line.startswith("ATOM"))
    lines[i] = lines[i][:30] + f"{99.:8.3f}" + lines[i][38:]
    ligand.write_text("\n".join(lines) + "\n")
    record["exported_ligand_sha256"] = sha256(ligand); publish()
    with pytest.raises(ValueError, match="alias export changed"):
        paired.read_vina_carb(Evidence(tmp_path), directory, row, paired.SEEDS[0], experiment, plan_pin, freeze_pin, names, executable_sha256=executable)
    record["status"] = "failed"; record["parameters"] = {**paired.PARAMETERS, "cpu": 4}; publish()
    with pytest.raises(ValueError, match="frozen effort"):
        paired.read_vina_carb(Evidence(tmp_path), directory, row, paired.SEEDS[0], experiment, plan_pin, freeze_pin, names, executable_sha256=executable)


def test_unexecuted_common_failure_remains_all_three_modern_seed_outcomes(tmp_path):
    condition = pin(tmp_path, "condition.json", {"case_id": "2RDK", "method": "fpocket", "jobs": [], "status": "preparation_failed"})
    row = {"case_id": "2RDK", "method": "fpocket", "condition_id": "strict/2RDK/fpocket/box01",
           "status": "preparation_failed", "source_condition": condition, "modern_runs": "outputs/vina/box01/runs.json"}
    config_pin = {"sha256": "a" * 64}
    result, pending = paired.read_modern(Evidence(tmp_path), row, {}, config_pin, ["L000"], input_manifest_sha256="b" * 64)
    assert result == {} and pending == ["missing_condition_execution"]
    receipt = {"condition_id": row["condition_id"], "configuration_sha256": config_pin["sha256"],
               "frozen_input_manifest_sha256": "b" * 64, "source_condition": condition, "status": "prepared_input_failure",
               "native_search_started": False, "no_runs_or_poses_invented": True, "eligible_seed_slots": paired.SEEDS}
    pin(tmp_path, "outputs/vina/execution.json", receipt)
    result, pending = paired.read_modern(Evidence(tmp_path), row, {}, config_pin, ["L000"], input_manifest_sha256="b" * 64)
    assert set(result) == set(paired.SEEDS) and pending == []
    assert all(r["status"] == "prepared_input_failure" and r["poses"] == [] for r in result.values())
    for key, value in (("configuration_sha256", "0" * 64), ("condition_id", "wrong"), ("eligible_seed_slots", paired.SEEDS[:2])):
        wrong = {**receipt, key: value}; pin(tmp_path, "outputs/vina/execution.json", wrong)
        with pytest.raises(ValueError, match="identity|three-seed"):
            paired.read_modern(Evidence(tmp_path), row, {}, config_pin, ["L000"], input_manifest_sha256="b" * 64)
    row["source_condition"] = pin(tmp_path, "condition.json", {"case_id": "2RDK", "method": "fpocket", "jobs": [{"status": "prepared"}]})
    pin(tmp_path, "outputs/vina/execution.json", {**receipt, "source_condition": row["source_condition"]})
    with pytest.raises(ValueError, match="hides a prepared"):
        paired.read_modern(Evidence(tmp_path), row, {}, config_pin, ["L000"], input_manifest_sha256="b" * 64)


def test_complete_modern_receipt_raw_poses_and_exact_input_pairing(tmp_path):
    row, _, export, *_ = synthetic_vc(tmp_path)
    names = [a["generated_name"] for a in export["atoms"]]
    config = {"seeds": paired.SEEDS, "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
              "energy_range_kcal_mol": 10, "min_rmsd_A": .5}
    config_pin = pin(tmp_path, "modern_config.json", config)
    row["modern_vina_prepared_directory"] = "inputs"
    row["receptor"] = pin(tmp_path, "inputs/neighborhood.pdbqt", "paired frozen receptor", raw=True)
    pin(tmp_path, "inputs/neighborhood.pdb", "paired observed receptor", raw=True)
    reference = pin(tmp_path, "inputs/reference.json", {"names": names, "pose_reference_eligible": False})
    prep = {"case_id": row["case_id"], "ligand_file": "ligand.pdbqt", "box_center": row["box"]["center"], "box_size": row["box"]["size"], "configuration_sha256": config_pin["sha256"]}
    row["preparation"] = pin(tmp_path, "inputs/preparation.json", prep)
    row["candidate_class_rank"] = 1
    row["source_condition"] = pin(tmp_path, "inputs/condition.json", {"case_id": row["case_id"], "method": "fpocket",
        "jobs": [{"candidate_class_rank": 1, "status": "prepared", "preparation_sha256": row["preparation"]["sha256"], "center_A": row["box"]["center"]}]})
    row["modern_runs"] = "outputs/modern/vina/box01/runs.json"
    directory = tmp_path / Path(row["modern_runs"]).parent; directory.mkdir(parents=True)
    (directory / "configuration.json").write_bytes((tmp_path / config_pin["path"]).read_bytes())
    pin(tmp_path, "benchmarks/vina.lock.json", {"assets": {"synthetic": {"sha256": "e" * 64}}})
    fp = {"config_sha256": config_pin["sha256"], "preparation_sha256": row["preparation"]["sha256"], "ligand_sha256": row["ligand"]["sha256"],
          "receptor_sha256": row["receptor"]["sha256"], "receptor_pdb_sha256": sha256(tmp_path / "inputs/neighborhood.pdb"),
          "reference_sha256": reference["sha256"], "full_assembly_sha256": None, "executable_sha256": "e" * 64}
    runs = []
    for seed in paired.SEEDS:
        pose, log = f"neighborhood_{seed}.pdbqt", f"neighborhood_{seed}.log"
        (directory / pose).write_text("MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\n" + (tmp_path / row["ligand"]["path"]).read_text() + "ENDMDL\n")
        (directory / log).write_text("1 -5.000 0.000 0.000\n")
        command = ["synthetic-vina", "--ligand", str(tmp_path / row["ligand"]["path"]), "--receptor", str(tmp_path / row["receptor"]["path"]),
                   "--out", str(directory / pose), "--scoring", "vina", "--seed", str(seed), "--cpu", "2", "--exhaustiveness", "16", "--num_modes", "20",
                   "--energy_range", "10", "--min_rmsd", "0.5", "--center_x", "0", "--center_y", "0", "--center_z", "0", "--size_x", "24", "--size_y", "24", "--size_z", "24"]
        runs.append({"seed": seed, "fingerprint": fp, "case_id": row["case_id"], "receptor": "neighborhood", "vina_version": "AutoDock Vina v1.2.7",
                     "information_regime": "predicted-site-local", "command": command, "status": "ok", "exit_code": 0, "pose_file": pose, "log_file": log,
                     "pose_sha256": sha256(directory / pose), "returned_poses": 1, "cpu_seconds": 2., "wall_seconds": 1.})
    receipt = {"condition_id": row["condition_id"], "configuration_sha256": config_pin["sha256"], "frozen_input_manifest_sha256": "f" * 64,
               "source_condition": row["source_condition"], "status": "complete", "native_search_started": True, "terminal_seed_runs": 3, "runs_path": row["modern_runs"]}
    def publish():
        write_json(directory / "runs.json", runs)
        receipt.update(runs_sha256=sha256(directory / "runs.json"), seed_statuses={str(r["seed"]): r["status"] for r in runs})
        write_json(directory.parent / "execution.json", receipt)
    publish()
    result, pending = paired.read_modern(Evidence(tmp_path), row, config, config_pin, names, input_manifest_sha256="f" * 64)
    assert pending == [] and len(result) == 3 and all(r["poses"][0]["method_score"] == -5 for r in result.values())
    wrong = copy.deepcopy(row); wrong["ligand"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="exact paired input"):
        paired.read_modern(Evidence(tmp_path), wrong, config, config_pin, names, input_manifest_sha256="f" * 64)
    runs[0]["command"][runs[0]["command"].index("--cpu") + 1] = "4"; publish()
    with pytest.raises(ValueError, match="frozen paired"):
        paired.read_modern(Evidence(tmp_path), row, config, config_pin, names, input_manifest_sha256="f" * 64)
    runs.pop(); publish()
    with pytest.raises(ValueError, match="three actual terminal"):
        paired.read_modern(Evidence(tmp_path), row, config, config_pin, names, input_manifest_sha256="f" * 64)


def test_portable_v2_archive_requires_explicit_pinned_filename(tmp_path):
    content = b"portable observed source"
    archive = tmp_path / "archive.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        member = tarfile.TarInfo("case/receptor.cif"); member.size = len(content)
        tar.addfile(member, io.BytesIO(content))
    record = pin(tmp_path, "manifest.json", {"archive_sha256": sha256(archive), "files": [
        {"path": "case/receptor.cif", "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}]})
    with pytest.raises(ValueError, match="explicit frozen"):
        paired.portable_archive(Evidence(tmp_path), record)
    assert paired.portable_archive(Evidence(tmp_path), record, archive_path="archive.tar.gz") == {"case/receptor.cif": content}
    archive.write_bytes(archive.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="checksum"):
        paired.portable_archive(Evidence(tmp_path), record, archive_path="archive.tar.gz")


def test_global_paired_pending_gate_keeps216_outcomes_and_never_enters_native_metrics(tmp_path, monkeypatch):
    modern = {"seeds": paired.SEEDS, **{k: paired.PARAMETERS[k] for k in ("cpu", "exhaustiveness", "num_modes", "energy_range_kcal_mol")}}
    entry = lambda value: {"path": "unused", "sha256": "a" * 64, "value": value}
    ref = {"source_ledger": {}, "reference_source_manifest": {}, "family_audit": {}, "case_to_family": {},
           "site_prediction": {"configuration": entry({"source_preparation_archive": {"manifest": "unused", "manifest_sha256": "a" * 64, "path": "unused.tar.gz"}}),
                               "freeze_receipt": entry({"configuration_sha256": "a" * 64, "frozen_before_new_predictions": True}), "source_archive_manifest": entry({})},
           "profiles": [{"id": profile, "source_archive_manifest": {}} for profile in paired.PROFILES]}
    experiment = {"adapter_source": entry({}), "secondary_top5_enabled": False, "native_binary_sha256": "d" * 64,
                  "helper_sources": {}, "paired_modern_vina_policy": entry(modern),
                  **{k: entry({}) for k in ("method_lock", "initialization_policy", "initialization_algorithm", "prepared_inputs", "native_build_receipt", "native_preflight_receipt")}}
    plan = {"case_ids": list(paired.CASE_IDS), "seeds": paired.SEEDS, "executed_candidate_budget": 1, "search_methods": list(paired.SEARCH_METHODS),
            "common_pool_rescoring_enabled": False, "metric_code_sha256": {}, "reference_evaluation_plan": entry(ref),
            "public_source_distribution": entry({}),
            "reference_evaluation_freeze": entry({"evaluation_plan_sha256": "a" * 64, "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False}),
            "vina_carb_experiment": entry(experiment), "vina_carb_freeze": entry({}), "vina_carb_executable_sha256": "d" * 64,
            "source_profiles": {profile: profile for profile in paired.PROFILES}}
    plan_record = pin(tmp_path, "plan.json", plan)
    pin(tmp_path, "freeze.json", {"evaluation_plan_sha256": plan_record["sha256"], "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    rows = [{"case_id": cid, "receptor_profile": profile, "method": method, "condition_id": f"{profile}/{cid}/{method}",
             "status": "preparation_failed", "source_condition": entry({}), "vina_carb_outputs": {str(s): f"output/{profile}/{cid}/{method}/{s}" for s in paired.SEEDS}}
            for profile in paired.PROFILES for cid in paired.CASE_IDS for method in paired.METHODS]
    sources = {f"{cid}/{name}": json.dumps({"atoms": [{"name": "L000"}]}).encode() if name == "ligand_atom_map.json" else b"fixed source"
               for cid in paired.CASE_IDS for name in ("ligand_atom_map.json", "ligand.sdf", "ligand.pdbqt", "protein_atom_map.json", "receptor.cif")}
    monkeypatch.setattr(paired, "CODE_SOURCES", ())
    monkeypatch.setattr(paired, "_pin", lambda evidence, record: json.dumps(record.get("value", {})).encode())
    monkeypatch.setattr(paired, "check_public_reference_sources", lambda *a: ({cid: {} for cid in paired.CASE_IDS}, {"all_consumed_numerical_sources_verified": True}))
    monkeypatch.setattr(paired, "archive_files", lambda *a: sources)
    monkeypatch.setattr(paired, "portable_archive", lambda *a, **kw: sources)
    monkeypatch.setattr(paired, "validate_frozen_plan", lambda *a: (experiment, rows))
    monkeypatch.setattr(paired, "canonical_mapping", lambda *a: {})
    monkeypatch.setattr(paired, "common_selection", lambda *a: ({}, []))
    monkeypatch.setattr(paired, "read_modern", lambda *a, **kw: ({s: {"status": "prepared_input_failure", "poses": [], "record": None} for s in paired.SEEDS}, []))
    monkeypatch.setattr(paired, "read_vina_carb", lambda *a, **kw: (None, "missing_execution"))
    monkeypatch.setattr(paired, "_source_context", lambda *a: pytest.fail("Native fitting entered before all paired outcomes terminal"))
    with pytest.raises(ValueError, match="remains incomplete"):
        paired.generate_paired_report(tmp_path, "plan.json", "freeze.json")
    result = paired.generate_paired_report(tmp_path, "plan.json", "freeze.json", allow_pending=True)
    assert result["counts"]["eligible_method_outcomes"] == 216 and result["counts"]["terminal_method_outcomes"] == 108
    assert result["counts"]["eligible_seed_pairs"] == 108 and len(result["pending"]) == 108
    assert not result["metrics_eligible"] and "conditions" not in result and "paired_family_statistics" not in result


def test_paired_family_statistics_compare_recovery_not_native_score_scales(monkeypatch):
    conditions = []
    for profile in paired.PROFILES:
        for site in paired.METHODS:
            for cid in paired.CASE_IDS:
                methods = {}
                for method in paired.SEARCH_METHODS:
                    succeeds = method == paired.SEARCH_METHODS[1] and cid != "4R9F"
                    metric = accuracy_summary([{"rmsd": 0.}] if succeeds else [], "rmsd")
                    methods[method] = {"global": metric, "local_interface": metric, "poses": [{"native_score": 1e9}] if succeeds else []}
                conditions.append({"case_id": cid, "receptor_profile": profile, "site_condition": site,
                    "atomic_scope": "modeled_core" if cid == "4R9F" else "full_ligand",
                    "seeds": [{"seed": s, "search_methods": methods} for s in paired.SEEDS]})
    from aav_benchmark.statistics import paired_family_bootstrap
    cache = {}
    def bootstrap(records, a, b, **kwargs):
        assert kwargs == {"n_resamples": 10000, "seed": 20261103, "expected_seeds": paired.SEEDS}
        assert all("native_score" not in r for r in records)
        key = len(records)
        if key not in cache:
            cache[key] = paired_family_bootstrap(records, a, b, **kwargs)
        return copy.deepcopy(cache[key])
    monkeypatch.setattr(paired, "paired_family_bootstrap", bootstrap)
    result = paired.paired_statistics(conditions, {cid: cid for cid in paired.CASE_IDS})
    assert len(result) == 216
    assert all(r["difference_a_minus_b"] == 1 for r in result if r["atomic_scope"] == "full_ligand")
    assert all(r["difference_a_minus_b"] == 0 and r["confidence_interval"] is None for r in result if r["atomic_scope"] == "modeled_core")
    assert all(r["executed_candidate_budget"] == 1 for r in result)
