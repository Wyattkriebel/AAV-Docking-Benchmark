"""Prospective full search must be separately frozen and retain failures."""
import json
from pathlib import Path
import subprocess
import hashlib

import pytest

from aav_benchmark import prospective_vina_carb as vc
from aav_benchmark.provenance import sha256, write_json


def freeze_all_failed(tmp_path):
    files = {"method_lock": "benchmarks/additional_method.lock.json", "initialization_policy": "inputs/policy.json",
             "initialization_algorithm": "inputs/algorithm.py", "prepared_inputs": "inputs/prepared.json",
             "adapter_source": "aav_benchmark/prospective_vina_carb.py", "paired_modern_vina_policy": "inputs/vina.json",
             "native_build_receipt": "inputs/native_build.json", "native_preflight_receipt": "inputs/native_preflight.json"}
    rows = [{"case_id": case, "receptor_profile": profile, "method": method,
             "condition_id": f"{profile}/{case}/{method}", "status": "preparation_failed",
             "failure_reason": "Preserved original source/MSE preparation failure"}
            for case in vc.CASES for profile in vc.PROFILES for method in vc.METHODS]
    for row in rows:
        path = Path("inputs/source_conditions") / (row["condition_id"].replace("/", "_") + ".json")
        write_json(tmp_path / path, {"case_id": row["case_id"], "method": row["method"], "jobs": [], "status": "preparation_failed"})
        row["source_condition"] = {"path": path.as_posix(), "sha256": sha256(tmp_path / path)}
        row["modern_runs"] = "outputs/" + row["condition_id"] + "/vina/box01/runs.json"
        row["vina_carb_outputs"] = {str(seed): "outputs/" + row["condition_id"] + f"/vina_carb/box01/s{seed}" for seed in vc.SEEDS}
    for key, relative in files.items():
        target = tmp_path / relative;target.parent.mkdir(parents=True, exist_ok=True)
        if key == "adapter_source":
            target.write_bytes(Path(vc.__file__).read_bytes())
        elif key == "method_lock":
            target.write_bytes((vc.ROOT / relative).read_bytes())
        elif key == "native_build_receipt":
            write_json(target, {"executable_sha256": hashlib.sha256(b"test executable").hexdigest()})
        elif key == "native_preflight_receipt":
            target.write_bytes((vc.ROOT / "results/common_pose_rescoring/scorer_source/fixed_pose_scorer.json").read_bytes())
        else:
            write_json(target, {"conditions": rows} if key == "prepared_inputs" else {})
    plan = {"status": "frozen_before_native_full_search", "frozen_before_native_execution": True,
            "new_method_outcomes_consulted": False, "cases": vc.CASES, "seeds": vc.SEEDS,
            "methods": vc.METHODS, "receptor_profiles": vc.PROFILES, "parameters": vc.PARAMETERS,
            "secondary_top5_enabled": False, "initialization_profile_id": "ccd_ideal_pyranose_v3",
            "reuse_byte_identical_profile_jobs": False,
            "native_binary_sha256": hashlib.sha256(b"test executable").hexdigest(),
            "native_preflight_binary_sha256": json.loads((vc.ROOT / "results/common_pose_rescoring/scorer_source/fixed_pose_scorer.json").read_text())["executable_sha256"]}
    plan.update({key: {"path": relative, "sha256": sha256(tmp_path / relative)} for key, relative in files.items()})
    for relative in vc.HELPER_SOURCES:
        target = tmp_path / relative;target.parent.mkdir(parents=True, exist_ok=True);target.write_bytes((vc.ROOT / relative).read_bytes())
    plan["helper_sources"] = {relative: sha256(tmp_path / relative) for relative in vc.HELPER_SOURCES}
    path = tmp_path / "plan.json";write_json(path, plan)
    receipt = {"configuration_sha256": sha256(path), "frozen_before_native_execution": True,
               "new_method_outcomes_consulted": False,
               "helper_sources": plan["helper_sources"],
               **{key + "_sha256": plan[key]["sha256"] for key in files}}
    receipt_path = tmp_path / "receipt.json";write_json(receipt_path, receipt)
    return path, receipt_path, rows


def rebind_prepared(plan, receipt, rows):
    config = json.loads(plan.read_text());root = plan.parent
    write_json(root / config["prepared_inputs"]["path"], {"conditions": rows})
    config["prepared_inputs"]["sha256"] = sha256(root / config["prepared_inputs"]["path"])
    write_json(plan, config)
    frozen = json.loads(receipt.read_text());frozen["prepared_inputs_sha256"] = config["prepared_inputs"]["sha256"]
    frozen["configuration_sha256"] = sha256(plan);write_json(receipt, frozen)


def one_prepared(tmp_path):
    plan, receipt, rows = freeze_all_failed(tmp_path);config = json.loads(plan.read_text())
    directory = Path("results/structural_pose_v3/strict/2RDK/fpocket/prepared/box01")
    absolute = tmp_path / directory;absolute.mkdir(parents=True)
    def file(name, value, json_value=False):
        path = absolute / name
        if json_value:
            write_json(path, value)
        else:
            path.write_text(value)
        return {"path": (directory / name).as_posix(), "sha256": sha256(path)}
    ligand = file("ligand.pdbqt", "paired unchanged ligand")
    receptor = file("neighborhood.pdbqt", "paired unchanged receptor")
    mapping = file("ligand_atom_map.json", {"atoms": [{"name": "L000", "branch_num": "1", "ccd_id": "MAN", "atom_name": "C1"},
                                                     {"name": "L001", "branch_num": "2", "ccd_id": "MAN", "atom_name": "O4"}]}, True)
    sdf = file("initializer.sdf", "independent serialized conformer")
    initializer = file("initializer.json", {"case_id": "2RDK", "profile_id": config["initialization_profile_id"],
                       "policy_sha256": config["initialization_policy"]["sha256"],
                       "algorithm_source_sha256": config["initialization_algorithm"]["sha256"],
                       "search": {"status": "selected"}, "serialized_sdf_validation": {"valid": True},
                       "selected_sdf_sha256": sdf["sha256"], "no_docking_or_accuracy_evaluation": True}, True)
    box = {"center": [1., 2., 3.], "size": [24., 24., 24.]}
    preparation = file("preparation.json", {"case_id": "2RDK", "ligand_file": "ligand.pdbqt", "box_center": box["center"],
                       "box_size": box["size"], "configuration_sha256": config["paired_modern_vina_policy"]["sha256"],
                       "native_coordinates_in_runner_reference": False}, True)
    rows[0] = {**rows[0], "status": "prepared", "candidate_class_rank": 1, "candidate_original_rank": 2, "initializer_profile_id": config["initialization_profile_id"],
               "ligand": ligand, "receptor": receptor, "canonical_atom_map": mapping, "initializer_sdf": sdf,
               "initializer_receipt": initializer, "preparation": preparation, "box": box,
               "modern_vina_prepared_directory": directory.as_posix()}
    source_path = tmp_path / rows[0]["source_condition"]["path"]
    write_json(source_path, {"case_id": "2RDK", "method": "fpocket", "jobs": [{"candidate_class_rank": 1,
               "candidate_original_rank": 2, "preparation_sha256": preparation["sha256"], "status": "prepared"}], "status": "prepared"})
    rows[0]["source_condition"]["sha256"] = sha256(source_path)
    rebind_prepared(plan, receipt, rows)
    return plan, receipt, rows


def test_unfrozen_proposal_cannot_start_any_native_call(tmp_path, monkeypatch):
    proposal = tmp_path / "proposal.json";write_json(proposal, {"status": "proposal", "frozen_before_native_execution": False})
    receipt = tmp_path / "receipt.json";write_json(receipt, {})
    monkeypatch.setattr(vc, "validate_installation", lambda *a, **k: pytest.fail("Unfrozen native validation reached"))
    with pytest.raises(ValueError, match="separate root freeze"):
        vc.execute_one(root=tmp_path, plan_path=proposal, receipt_path=receipt, condition_id="unused", seed=104729,
                       installation=tmp_path / "unused", output=tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("profile", ["strict", "repair"])
def test_all_six_preparation_failures_remain_scheduled_without_native_calls(tmp_path, monkeypatch, profile):
    plan, receipt, rows = freeze_all_failed(tmp_path)
    _, checked = vc.validate_frozen_plan(tmp_path, plan, receipt)
    assert len(checked) == 36 and {row["case_id"] for row in checked} == set(vc.CASES)
    monkeypatch.setattr(vc, "validate_installation", lambda *a, **k: pytest.fail("Failed preparation invoked native tool"))
    selected = next(row for row in rows if row["receptor_profile"] == profile)
    output = tmp_path / selected["vina_carb_outputs"]["104729"]
    result = vc.execute_one(root=tmp_path, plan_path=plan, receipt_path=receipt, condition_id=selected["condition_id"],
                            seed=104729, installation=tmp_path / "unused", output=output)
    assert result["status"] == "prepared_input_failure" and result["full_search_started"] is False
    assert result["atomic_accuracy_computed"] is False
    assert json.loads((output / "execution.json").read_text())["condition"]["case_id"] == vc.CASES[0]


def test_freeze_rejects_missing_cases_and_duplicate_budget_attempts(tmp_path):
    plan, receipt, rows = freeze_all_failed(tmp_path)
    for altered, message in ((rows[:-1], "Every six-case"), (rows + [{**rows[0], "condition_id": "duplicate_rank"}], "candidate budget")):
        rebind_prepared(plan, receipt, altered)
        with pytest.raises(ValueError, match=message):
            vc.validate_frozen_plan(tmp_path, plan, receipt)


def test_imported_helper_source_change_is_rejected_before_native_execution(tmp_path):
    plan, receipt, _ = freeze_all_failed(tmp_path)
    helper = tmp_path / "aav_benchmark/common_pose_rescoring.py"
    helper.write_text(helper.read_text() + "\n# changed after freeze\n")
    with pytest.raises(ValueError, match="prepared input changed"):
        vc.validate_frozen_plan(tmp_path, plan, receipt)


def test_changed_ligand_box_or_nonfinite_values_cannot_be_presented_as_paired_inputs(tmp_path):
    plan, receipt, rows = one_prepared(tmp_path)
    assert vc.validate_frozen_plan(tmp_path, plan, receipt)
    original = json.loads(json.dumps(rows))
    alternate = tmp_path / Path(rows[0]["ligand"]["path"]).with_name("alternate.pdbqt")
    alternate.write_text("different ligand")
    rows[0]["ligand"] = {"path": str(alternate.relative_to(tmp_path)), "sha256": sha256(alternate)}
    rebind_prepared(plan, receipt, rows)
    with pytest.raises(ValueError, match="same prepared ligand"):
        vc.validate_frozen_plan(tmp_path, plan, receipt)
    rows = json.loads(json.dumps(original));rows[0]["box"]["center"][0] += 1
    rebind_prepared(plan, receipt, rows)
    with pytest.raises(ValueError, match="same prepared ligand"):
        vc.validate_frozen_plan(tmp_path, plan, receipt)
    rows[0]["box"]["size"] = [float("nan"), 24., 24.]
    config = json.loads(plan.read_text());path = tmp_path / config["prepared_inputs"]["path"]
    path.write_text(json.dumps({"conditions": rows}))
    config["prepared_inputs"]["sha256"] = sha256(path);write_json(plan, config)
    frozen = json.loads(receipt.read_text());frozen.update(configuration_sha256=sha256(plan), prepared_inputs_sha256=sha256(path));write_json(receipt, frozen)
    with pytest.raises(ValueError, match="Invalid paired docking box"):
        vc.validate_frozen_plan(tmp_path, plan, receipt)


def test_native_exit_zero_with_zero_returned_models_is_an_explicit_failure(tmp_path, monkeypatch):
    plan, receipt, rows = one_prepared(tmp_path)
    executable = tmp_path / "binary";executable.write_text("test executable")
    installation = tmp_path / "fake-installation";installation.mkdir()
    (installation / "fixed_pose_scorer.json").write_bytes((vc.ROOT / "results/common_pose_rescoring/scorer_source/fixed_pose_scorer.json").read_bytes())
    monkeypatch.setattr(vc, "validate_installation", lambda *a, **k: executable)
    def export(source, mapping, target):
        target.write_text("test canonical pose")
        value = {"expected_glycosidic_links": [{"branch_serials": [1, 2]}],
                 "atoms": [{"serial": 1, "carb_resname": "MAN", "carb_resnum": 1, "hydrogen": False},
                           {"serial": 2, "carb_resname": "MAN", "carb_resnum": 2, "hydrogen": False}]}
        write_json(target.with_suffix(".pdbqt.map.json"), value);return value
    monkeypatch.setattr(vc, "export_vina_carb_ligand", export)
    monkeypatch.setattr(vc, "_invoke_fixed_pose", lambda *a, **k: {"glycosidic_links": 1, "fully_supported_links": 1})
    def run(command, **kwargs):
        out = Path(kwargs["cwd"]);(out / "poses.pdbqt").write_text("")
        (out / "VC_log.txt").write_text("CHI energy penalties applied to MAN 1 and MAN 2 linkage.\n")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(vc.subprocess, "run", run)
    result = vc.execute_one(root=tmp_path, plan_path=plan, receipt_path=receipt, condition_id=rows[0]["condition_id"],
                            seed=104729, installation=installation, output=tmp_path / rows[0]["vina_carb_outputs"]["104729"])
    assert result["full_search_started"] is True and result["exit_code"] == 0
    assert result["status"] == "failed" and "Missing or truncated" in result["error"]
    assert "poses.pdbqt" in result["artifacts"] and result["atomic_accuracy_computed"] is False


def test_config_uses_original_carbohydrate_full_search_parameters():
    row = {"box": {"center": [1., 2., 3.], "size": [24., 28., 32.]}}
    text = vc.render_configuration(row, 130363)
    assert "chi_coeff = 1\n" in text and "chi_cutoff = 2\n" in text
    assert "cpu = 2\n" in text and "exhaustiveness = 16\n" in text
    assert "num_modes = 20\n" in text and "energy_range = 10\n" in text
    assert "min_rmsd" not in text  # Original implementation hardcodes its separation.
    assert "score_only" not in text and "local_only" not in text
    with pytest.raises(ValueError, match="Unscheduled"):
        vc.render_configuration(row, 0)


def test_atom_adapter_forwards_chemical_identities_only():
    mapped = vc.normalize_atom_map({"canonical_atom_map": [{"name": "L000", "branch_num": "2", "ccd_id": "MAN",
                                                            "atom_name": "C1", "coordinates": [99, 99, 99]}]})
    assert mapped == [{"name": "L000", "source_key": ["B", 2, "", "MAN", "C1"]}]
    assert "coordinates" not in mapped[0]
    with pytest.raises(ValueError, match="canonical atom map"):
        vc.normalize_atom_map({"atoms": []})


def test_cost_counts_boxes_and_profiles_without_equal_time_claims():
    primary, secondary = vc.nominal_cost(False), vc.nominal_cost(True)
    assert primary["maximum_VinaCarb_seed_box_jobs_before_preparation_failures"] == 108
    assert primary["maximum_jobs_per_receptor_profile"] == 54
    assert secondary["maximum_VinaCarb_seed_box_jobs_before_preparation_failures"] == 396
    assert primary["wall_or_CPU_time_estimate_claimed"] is False
    with pytest.raises(ValueError):
        vc.nominal_cost(receptor_profiles=0)
