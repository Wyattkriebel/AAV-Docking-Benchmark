"""Writer-window I/O validation uses score tokens, never pose geometry."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from aav_benchmark.vina_output_audit import (
    audit_main_vina_outputs, score_token_bounds, validate_vina_1_2_7_output,
)

ROOT = Path(__file__).resolve().parents[1]


def raw_output(log_scores, saved_scores):
    log = ("AutoDock Vina v1.2.7\n\n"
           "mode |   affinity | dist from best mode\n"
           "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
           "-----+------------+----------+----------\n")
    log += "".join(f" {rank} {score} 0 0\n" for rank, score in enumerate(log_scores, 1))
    poses = "".join(f"MODEL {rank}\nREMARK VINA RESULT: {score} 0.000 0.000\nENDMDL\n"
                    for rank, score in enumerate(saved_scores, 1))
    return log.encode(), poses.encode()


@pytest.mark.parametrize("directory,stem,pose_sha,log_sha,lower", [
    ("pol_specificity/POL_GLC_beta__holo__C", "dry_155921",
     "2aeaa84c973063a54a3d3bbd9c5e5959a634d5390e6ca9026d8d56fbd8be2836",
     "0684ed2f2f12e4c6ed8be57b2e4356193a36012f84bcd20368f8bbb094302974", "17.6705"),
    ("development/aav9", "holo_monomer_155921",
     "ff878d026ffb6bfcf4b8196c111f9ee2540fb67be64b45c47aef216f713f1f7a",
     "5904eac9a710d9a06aabc98fd043fabf9d60b3d9d4f8ebf761beb57a69862299", "49.0805"),
    ("development/aav9", "apo_monomer_130363",
     "f4836af143cfe121566256e4abaeef375ad6b759a83b55984f37a6593ca131fe",
     "de7d845a10db57e0ed2470b4237a506c8d26d35c375feaa5f5abb0c9df8f6226", "25.8725"),
])
def test_three_actual_twenty_logged_nineteen_saved_controls(directory, stem, pose_sha, log_sha, lower):
    directory = ROOT / "results" / directory
    log, poses = (directory / f"{stem}.log").read_bytes(), (directory / f"{stem}.pdbqt").read_bytes()
    report = validate_vina_1_2_7_output(log, poses, expected_returned_poses=19)
    assert report["log"]["sha256"] == log_sha == hashlib.sha256(log).hexdigest()
    assert report["poses"]["sha256"] == pose_sha == hashlib.sha256(poses).hexdigest()
    assert report["logged_modes"] == 20 and report["saved_models"] == 19
    assert report["proven_outside_window_omissions"] == 1
    assert report["boundary_ambiguous_omissions"] == 0
    assert report["omitted_suffix"][0]["rank"] == 20
    assert report["omitted_suffix"][0]["delta_from_best_bounds_kcal_mol"]["lower"] == lower
    assert report["experimental_reference_access"] is report["accuracy_computed"] is False


def test_power_of_ten_general_format_bounds_are_asymmetric_and_scientific_zero_is_exact():
    def bounds(token):
        return score_token_bounds(token, format_kind="log_four_significant_digits")
    assert bounds("10")["lower"] == "9.9995"
    assert bounds("10")["upper"] == "10.005"
    assert bounds("-10")["lower"] == "-10.005"
    assert bounds("-10")["upper"] == "-9.9995"
    assert bounds("1e-05")["lower"] == "0.0000099995"
    assert bounds("1e-05")["upper"] == "0.000010005"
    assert bounds("-0")["lower"] in {"0", "-0"}
    assert Decimal(bounds("0")["upper"]) == 0
    fixed = score_token_bounds("0.000", format_kind="pdbqt_three_decimal_places")
    assert fixed["lower"] == "-0.0005" and fixed["upper"] == "0.0005"


def test_suffix_boundary_ambiguity_is_writer_compatible_but_not_proven_outside():
    log, poses = raw_output(["-5", "5"], ["-5.000"])
    report = validate_vina_1_2_7_output(log, poses)
    assert report["omitted_modes"] == report["boundary_ambiguous_omissions"] == 1
    assert report["proven_outside_window_omissions"] == 0
    assert report["omitted_suffix"][0]["window_classification"] == "boundary_ambiguous_writer_compatible"
    assert report["shared_writer_compatible_best_bounds_kcal_mol"]["upper_inclusive"] is False


@pytest.mark.parametrize("omitted", ["4.99", "4.999"])
def test_inside_window_omission_including_upper_bound_exactly_ten_is_rejected(omitted):
    log, poses = raw_output(["-5", omitted], ["-5.000"])
    with pytest.raises(ValueError, match="inside the output energy window"):
        validate_vina_1_2_7_output(log, poses)


def test_saved_boundary_is_reported_and_saved_pose_proven_outside_is_rejected():
    log, poses = raw_output(["0", "10"], ["0.000", "10.000"])
    report = validate_vina_1_2_7_output(log, poses)
    assert report["saved_boundary_ambiguous_ranks"] == [2]
    log, poses = raw_output(["0", "10.02"], ["0.000", "10.020"])
    with pytest.raises(ValueError, match="Saved model is proven outside"):
        validate_vina_1_2_7_output(log, poses)


def test_log_and_saved_intervals_must_intersect_before_window_checks():
    log, poses = raw_output(["-5.123"], ["-5.124"])
    with pytest.raises(ValueError, match="incompatible rounding intervals"):
        validate_vina_1_2_7_output(log, poses)
    log, poses = raw_output(["-5.123"], ["-5.125"])
    with pytest.raises(ValueError, match="incompatible rounding intervals"):
        validate_vina_1_2_7_output(log, poses)
    log, poses = raw_output(["0.0625"], ["0.062"])
    assert validate_vina_1_2_7_output(log, poses)["saved_models"] == 1


@pytest.mark.parametrize("mutation", ["log_gap", "saved_gap", "duplicate_rank", "malformed_last_row", "more_saved", "duplicate_score", "unterminated"])
def test_rank_prefix_structure_and_malformed_final_rows_fail_closed(mutation):
    log, poses = raw_output(["-5", "-4"], ["-5.000", "-4.000"])
    if mutation == "log_gap":
        log = log.replace(b" 2 -4", b" 3 -4")
    elif mutation == "saved_gap":
        poses = poses.replace(b"MODEL 2", b"MODEL 3")
    elif mutation == "duplicate_rank":
        log = log.replace(b" 2 -4", b" 1 -4")
    elif mutation == "malformed_last_row":
        log = log.replace(b" 2 -4 0 0", b" 2 -4 missing")
    elif mutation == "more_saved":
        log = raw_output(["-5"], ["-5.000"])[0]
    elif mutation == "duplicate_score":
        poses = poses.replace(b"ENDMDL", b"REMARK VINA RESULT: -5.000 0.000 0.000\nENDMDL", 1)
    else:
        poses = poses.removesuffix(b"ENDMDL\n")
    with pytest.raises(ValueError):
        validate_vina_1_2_7_output(log, poses)


def fixture_manifest(tmp_path, *, status="ok", verbosity=None):
    directory = tmp_path / "results/development/fixture"
    directory.mkdir(parents=True)
    log, poses = raw_output(["-5", "5.02"], ["-5.000"])
    (directory / "raw.log").write_bytes(log)
    (directory / "raw.pdbqt").write_bytes(poses)
    command = ["vina", "--scoring", "vina", "--energy_range", "10", "--num_modes", "20"]
    if verbosity is not None:
        command += ["--verbosity", verbosity]
    row = {"status": status, "vina_version": "AutoDock Vina v1.2.7", "seed": 104729,
           "receptor": "fixture", "command": command, "exit_code": 0,
           "pose_file": "raw.pdbqt", "log_file": "raw.log", "returned_poses": 1,
           "pose_sha256": hashlib.sha256(poses).hexdigest()}
    (directory / "runs.json").write_text(json.dumps([row]))
    return directory, row


def test_audit_keeps_pending_rows_hashes_every_raw_file_and_never_accesses_reference_files(tmp_path):
    directory, row = fixture_manifest(tmp_path)
    pending = {**row, "status": "running", "pose_file": "absent.pdbqt"}
    (directory / "runs.json").write_text(json.dumps([row, pending]))
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert report["status"] == "all_audited_outputs_writer_compatible"
    assert report["totals"]["audited_records"] == report["totals"]["source_status_running"] == 1
    assert report["records"][1]["status"] == "pending_source"
    assert len(report["records"][0]["log"]["sha256"]) == 64
    assert len(report["records"][0]["manifest"]["sha256"]) == 64
    assert report["experimental_reference_access"] is report["native_execution"] is False
    with pytest.raises(ValueError, match="existing main/control"):
        audit_main_vina_outputs(tmp_path, scopes=("structural_assembly_pose_v3",))


def test_cli_is_read_only_and_rejects_unsupported_verbosity_or_changed_raw_hash(tmp_path):
    directory, row = fixture_manifest(tmp_path, verbosity="2")
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert "verbosity 1" in report["records"][0]["error"]
    row["command"] = row["command"][:-2]
    row["pose_sha256"] = "0" * 64
    (directory / "runs.json").write_text(json.dumps([row]))
    before = {str(path.relative_to(tmp_path)): path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}
    process = subprocess.run([sys.executable, str(ROOT / "scripts/audit_vina_outputs.py"),
                              "--root", str(tmp_path), "--scope", "development"], capture_output=True, check=False)
    assert process.returncode == 1
    report = json.loads(process.stdout)
    assert "hash mismatch" in report["records"][0]["error"]
    assert before == {str(path.relative_to(tmp_path)): path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}


def test_injected_row_before_separator_and_nonobject_manifest_record_fail_closed(tmp_path):
    log, poses = raw_output(["-5"], ["-5.000"])
    log = log.replace(b"-----+", b" 1 -5 0 0\n-----+", 1)
    with pytest.raises(ValueError, match="separator"):
        validate_vina_1_2_7_output(log, poses)
    directory, _ = fixture_manifest(tmp_path)
    (directory / "runs.json").write_text("[null]")
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert report["status"] == "invalid_output_present"
    assert report["records"][0]["status"] == "invalid_record"


def test_exponent_overflow_fails_cleanly_and_exponent_zero_is_normalized():
    with pytest.raises(ValueError, match="finite score token"):
        score_token_bounds("1e999999999999999999999999", format_kind="log_four_significant_digits")
    bounds = score_token_bounds("0e-999999999", format_kind="log_four_significant_digits")
    assert bounds["lower"] == bounds["upper"] == "0"


def test_indented_extra_control_block_is_not_silently_ignored():
    log, poses = raw_output(["-5"], ["-5.000"])
    poses += b" MODEL 2\n REMARK VINA RESULT: -4.000 0.000 0.000\n ENDMDL\n"
    with pytest.raises(ValueError, match="Indented reserved"):
        validate_vina_1_2_7_output(log, poses)


@pytest.mark.parametrize("value", [-100., -10., -9.9995, -.0625, -.00001, 0., .00001, .0625, 9.9995, 10., 9999.5, 1e100])
def test_independent_binary_float_formatting_agrees_with_decimal_token_intervals(value):
    log, poses = raw_output([format(value, ".4g")], [format(value, ".3f")])
    assert validate_vina_1_2_7_output(log, poses)["saved_models"] == 1


def test_source_independent_software_artifact_replays_without_atom_records():
    fixture_path = ROOT / "tests/fixtures/vina_output_io_controls.json"
    fixtures = json.loads(fixture_path.read_bytes())
    artifact = json.loads((ROOT / "results/vina_output_io_audit/software_validation.json").read_bytes())
    for pin in artifact["software_files"]:
        raw = (ROOT / pin["path"]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == pin["sha256"] and len(raw) == pin["bytes"]
    actual = []
    for case in fixtures["cases"]:
        log, poses = case["log_text"].encode(), case["saved_text"].encode()
        assert "ATOM" not in case["saved_text"] and "HETATM" not in case["saved_text"]
        try:
            report = validate_vina_1_2_7_output(log, poses)
            result = {"id": case["id"], "status": "valid_writer_compatible",
                      "logged_modes": report["logged_modes"], "saved_models": report["saved_models"],
                      "proven_outside_window_omissions": report["proven_outside_window_omissions"],
                      "boundary_ambiguous_omissions": report["boundary_ambiguous_omissions"]}
        except ValueError as error:
            result = {"id": case["id"], "status": "invalid_output", "error": str(error)}
        assert result["status"] == case["expected_status"]
        result.update(log_sha256=hashlib.sha256(log).hexdigest(), saved_sha256=hashlib.sha256(poses).hexdigest())
        actual.append(result)
    assert actual == artifact["synthetic_control_results"]
    assert artifact["native_execution"] is artifact["experimental_reference_access"] is artifact["accuracy_computed"] is False


def test_empty_audit_cannot_claim_success_and_symlinked_manifest_is_not_read(tmp_path, monkeypatch):
    empty = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert empty["status"] == "no_auditable_outputs_found"
    directory = tmp_path / "results/development/fixture"
    directory.mkdir(parents=True)
    outside = tmp_path / "must_not_be_opened.json"
    outside.write_text("[]")
    (directory / "runs.json").symlink_to(outside)
    original = Path.read_bytes
    def guarded(path):
        if path.resolve() == outside:
            pytest.fail("Followed a manifest symlink outside the audit scope")
        return original(path)
    monkeypatch.setattr(Path, "read_bytes", guarded)
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert report["status"] == "invalid_output_present"
    assert "Symlinked" in report["records"][0]["error"]


def test_raw_symlink_and_nonraw_filename_are_rejected_before_open(tmp_path):
    directory, row = fixture_manifest(tmp_path)
    (directory / "raw.log").unlink()
    (directory / "raw.log").symlink_to(directory / "raw.pdbqt")
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert "symlinked" in report["records"][0]["error"]
    row["log_file"] = "reference.json"
    (directory / "runs.json").write_text(json.dumps([row]))
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert "filename or suffix" in report["records"][0]["error"]


@pytest.mark.parametrize("counts", [{"num_modes": True}, {"expected_returned_poses": True}])
def test_boolean_counts_are_not_accepted_as_integer_metadata(counts):
    log, poses = raw_output(["-5"], ["-5.000"])
    with pytest.raises(ValueError):
        validate_vina_1_2_7_output(log, poses, **counts)


@pytest.mark.parametrize("kind", ["scope", "manifest", "raw"])
def test_symlink_loops_return_structured_invalid_results(tmp_path, kind):
    if kind == "scope":
        (tmp_path / "results").mkdir()
        path = tmp_path / "results/development"
    else:
        directory, _ = fixture_manifest(tmp_path)
        path = directory / ("runs.json" if kind == "manifest" else "raw.log")
        path.unlink()
    path.symlink_to(path.name)
    report = audit_main_vina_outputs(tmp_path, scopes=("development",))
    assert report["status"] == "invalid_output_present"
    assert report["records"][0]["status"].startswith("invalid")
