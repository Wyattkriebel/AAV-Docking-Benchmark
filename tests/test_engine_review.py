"""Independent regression checks for run provenance and complete denominators.

Subprocesses are simulated: these checks exercise orchestration and analysis,
and do not present mocked coordinates as experimental docking evidence.
"""

import json
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

import gemmi
import numpy as np
import pytest

from aav_benchmark import cli, engine
from aav_benchmark.provenance import fetch_verified, sha256, write_json


def _receptor(path, coordinates):
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    chain = gemmi.Chain("A")
    for number, coordinate in enumerate(coordinates, 1):
        residue = gemmi.Residue()
        residue.name = "ALA"
        residue.seqid = gemmi.SeqId(number, " ")
        residue.het_flag = "A"
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(*coordinate)
        residue.add_atom(atom)
        chain.add_residue(residue)
    structure[0].add_chain(chain)
    structure.write_pdb(str(path))


def _pose_text(x=2.0):
    atom = f"ATOM  {1:5d} {'C1':>4s} GAL A{1:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00    +0.000 C"
    return "MODEL 1\nREMARK VINA RESULT: -5.000 0.000 0.000\n" + atom + "\nENDMDL\n"


@pytest.fixture
def fixture_files(tmp_path):
    prepared = tmp_path / "prepared"
    output = tmp_path / "output"
    prepared.mkdir()
    _receptor(prepared / "holo.pdb", [[0, 0, 0]])
    (prepared / "holo.pdbqt").write_text("receptor preparation bytes\n")
    (prepared / "ligand.pdbqt").write_text("ligand preparation bytes\n")
    reference = {
        "names": ["C1"], "coordinates": [[2, 0, 0]],
        "full_assembly_file": "full_assembly.npz",
        "receptor_atom_keys": [["A", 1, "", "ALA", "CA"]],
        "receptor_coordinates": [[0, 0, 0]],
    }
    write_json(prepared / "reference.json", reference)
    preparation = {
        "case_id": "analytic_control", "family_id": "analytic_control",
        "ligand_file": "ligand.pdbqt", "box_center": [0, 0, 0], "box_size": [20, 20, 20],
        "receptors": [{"id": "holo", "state": "holo", "assembly_condition": "complete-neighborhood",
                       "pdb_file": "holo.pdb", "pdbqt_file": "holo.pdbqt"}],
    }
    write_json(prepared / "preparation.json", preparation)
    np.savez(prepared / "full_assembly.npz", holo_coordinates=[[0., 0., 0.], [2.1, 0., 0.]])
    config = tmp_path / "config.json"
    write_json(config, {"seeds": [104729, 130363], "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
                        "min_rmsd_A": 0.5, "energy_range_kcal_mol": 10})
    executable = tmp_path / "vina"
    executable.write_text("frozen executable bytes\n")
    return prepared, output, executable, config


def _run(fixture_files, *, resume=False, exit_code=0, write_output=True, pose_x=2.0):
    prepared, output, executable, config = fixture_files

    def subprocess_run(command, **kwargs):
        if write_output:
            Path(command[command.index("--out") + 1]).write_text(_pose_text(pose_x))
        return CompletedProcess(command, exit_code)

    with patch.object(engine.subprocess, "check_output", return_value="AutoDock Vina v1.2.7\n"), \
         patch.object(engine.subprocess, "run", side_effect=subprocess_run), \
         patch.object(engine.platform, "platform", return_value="analytic-test-platform"), \
         patch.object(engine.importlib.metadata, "version", return_value="test-version"):
        return engine.dock(prepared, output, executable, config, resume=resume)


def test_analysis_rejects_truncated_completed_manifest(fixture_files):
    prepared, output, _, _ = fixture_files
    runs = _run(fixture_files)
    write_json(output / "runs.json", runs[:1])
    with pytest.raises(ValueError):
        engine.analyze(prepared, output)


def test_analysis_rejects_duplicate_run_identity(fixture_files):
    prepared, output, _, _ = fixture_files
    runs = _run(fixture_files)
    write_json(output / "runs.json", runs + [runs[0]])
    with pytest.raises(ValueError):
        engine.analyze(prepared, output)


def test_blank_insertion_codes_match_native_contact_identities(fixture_files):
    prepared, output, _, _ = fixture_files
    reference_path = prepared / "reference.json"
    reference = json.loads(reference_path.read_text())
    reference["receptor_atom_keys"][0][2] = " "
    write_json(reference_path, reference)
    _run(fixture_files)
    summary = engine.analyze(prepared, output)
    assert all(row["top1_contact_recall"] == 1.0 for row in summary)
    assert all(row["top1_contact_precision"] == 1.0 for row in summary)
    assert (output / "analysis.json").is_file()


def test_panel_cli_retains_later_cases_and_reports_any_failure(tmp_path):
    (tmp_path / "benchmarks").mkdir()
    write_json(tmp_path / "benchmarks/gal3_run.json", {"cases": ["first", "second"]})
    with patch.object(engine, "dock", side_effect=[[{"status": "failed"}], [{"status": "ok"}]]) as mocked:
        status = cli.main(["--root", str(tmp_path), "panel", "dock", "--results-root", str(tmp_path / "replay")])
    assert status == 1
    assert mocked.call_count == 2
    assert mocked.call_args_list[1].args[1] == tmp_path / "replay/second"


def test_analysis_rejects_changed_receptor_scoring_coordinates(fixture_files):
    prepared, output, _, _ = fixture_files
    _run(fixture_files)
    _receptor(prepared / "holo.pdb", [[20, 0, 0]])
    with pytest.raises(ValueError):
        engine.analyze(prepared, output)


def test_analysis_rejects_changed_full_assembly(fixture_files):
    prepared, output, _, _ = fixture_files
    _run(fixture_files)
    np.savez(prepared / "full_assembly.npz", holo_coordinates=[[200., 0., 0.]])
    with pytest.raises(ValueError):
        engine.analyze(prepared, output)


def test_missing_full_assembly_cannot_be_labeled_full_assembly(fixture_files):
    prepared, output, _, _ = fixture_files
    _run(fixture_files)
    (prepared / "full_assembly.npz").unlink()
    with pytest.raises((ValueError, FileNotFoundError)):
        engine.analyze(prepared, output)


def test_undeclared_full_assembly_is_explicitly_labeled_provided_receptor_only(fixture_files):
    prepared, output, _, _ = fixture_files
    reference_path = prepared / "reference.json"
    reference = json.loads(reference_path.read_text())
    reference.pop("full_assembly_file")
    write_json(reference_path, reference)
    _run(fixture_files)
    engine.analyze(prepared, output)
    import csv
    with (output / "pose_metrics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert all(row["evaluation_receptor_scope"] == "provided_receptor_only" for row in rows)
    assert all(float(row["evaluation_min_distance_A"]) == pytest.approx(2) for row in rows)


def test_retry_cannot_accept_output_left_by_a_failed_attempt(fixture_files):
    first = _run(fixture_files, exit_code=1)
    assert all(run["status"] == "failed" for run in first)
    retried = _run(fixture_files, resume=True, write_output=False)
    assert all(run["status"] == "failed" for run in retried)


def test_pose_translation_remains_error_and_clashes_use_full_assembly(fixture_files):
    prepared, output, _, _ = fixture_files
    np.savez(prepared / "full_assembly.npz", holo_coordinates=[[0., 0., 0.], [7.1, 0., 0.]])
    _run(fixture_files, pose_x=7)
    summaries = engine.analyze(prepared, output)
    assert all(run["top1_rmsd_A"] == pytest.approx(5) for run in summaries)
    import csv
    with (output / "pose_metrics.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert all(row["evaluation_receptor_scope"] == "full_biological_assembly" for row in rows)
    assert all(float(row["evaluation_min_distance_A"]) == pytest.approx(0.1) for row in rows)
    assert all(int(row["evaluation_pairs_below_1_5A"]) == 1 for row in rows)


def test_failed_runs_remain_in_summary_denominator(fixture_files):
    prepared, output, _, _ = fixture_files
    runs = _run(fixture_files, exit_code=1)
    summaries = engine.analyze(prepared, output)
    assert len(summaries) == len(runs) == 2
    assert all(run["returned_poses"] == 0 for run in summaries)
    assert all(run["top1_success_2A"] is False for run in summaries)


def test_resume_rejects_changed_frozen_configuration(fixture_files):
    _, _, _, config_path = fixture_files
    _run(fixture_files)
    config = json.loads(config_path.read_text())
    config["exhaustiveness"] = 64
    write_json(config_path, config)
    with pytest.raises(ValueError):
        _run(fixture_files, resume=True)


def test_cli_dock_returns_failure_when_any_run_fails(tmp_path):
    executable = tmp_path / "vina"
    executable.touch()
    with patch.object(engine, "dock", return_value=[{"status": "failed"}]):
        result = cli.main(["--root", str(tmp_path), "dock", "aav9", "--vina", str(executable)])
    assert result == 1


def test_nonfinite_manifest_cannot_replace_previous_valid_record(tmp_path):
    path = tmp_path / "record.json"
    write_json(path, {"value": 1})
    previous = path.read_bytes()
    with pytest.raises(ValueError):
        write_json(path, {"value": float("nan")})
    assert path.read_bytes() == previous
    assert sorted(p.name for p in tmp_path.iterdir()) == ["record.json"]


def test_verified_fetch_rejects_mutated_existing_source(tmp_path):
    path = tmp_path / "source"
    path.write_text("original")
    original_digest = sha256(path)
    path.write_text("changed")
    with pytest.raises(ValueError):
        fetch_verified("https://example.org/source", path, original_digest)
    assert path.read_text() == "changed"
