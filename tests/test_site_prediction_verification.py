"""Input integrity and raw prediction replay without installed predictors."""
import copy
import gzip
import json
from pathlib import Path
import shutil

import pytest

from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.site_archive import pack_run, restore_run
from aav_benchmark.site_prediction_verification import (
    _input_evidence, _same_predictions, _verify_attempt, restore_site_inputs,
)


ROOT = Path(__file__).resolve().parents[1]


def _input_archive(tmp_path):
    directory = tmp_path / "repo/results/site_discovery/example/inputs"
    directory.mkdir(parents=True)
    data = b"data_test\n#\n"
    archive = directory / "receptor.gz"
    archive.write_bytes(gzip.compress(data, mtime=0))
    plain = tmp_path / "source.cif"
    plain.write_bytes(data)
    manifest = {"condition_id": "example", "files": [{"revision": "original", "filename": "receptor.cif",
                "archive": archive.name, "sha256": sha256(plain), "archive_sha256": sha256(archive)}]}
    write_json(directory / "manifest.json", manifest)
    return tmp_path / "repo", directory, manifest, data


def test_restore_inputs_is_portable_and_requires_fresh_destination(tmp_path):
    root, directory, manifest, data = _input_archive(tmp_path)
    destination = tmp_path / "different/workstation/input"
    receipt = restore_site_inputs(root, "example", "original", destination)
    assert (destination / "receptor.cif").read_bytes() == data
    assert receipt["files"] == manifest["files"]
    with pytest.raises(FileExistsError, match="fresh"):
        restore_site_inputs(root, "example", "original", destination)
    assert (destination / "receptor.cif").read_bytes() == data


@pytest.mark.parametrize("field,value", [("filename", "../escape"), ("filename", "/absolute"),
                                        ("archive", "../outside.gz"), ("archive", "a\\b"),
                                        ("sha256", "F" * 64)])
def test_restore_rejects_unsafe_names_and_bad_digests(tmp_path, field, value):
    root, directory, manifest, _ = _input_archive(tmp_path)
    manifest["files"][0][field] = value
    write_json(directory / "manifest.json", manifest)
    destination = tmp_path / "restored"
    with pytest.raises(ValueError):
        restore_site_inputs(root, "example", "original", destination)
    assert not destination.exists()
    assert not (tmp_path / "escape").exists()


@pytest.mark.parametrize("attack", ["condition", "duplicate", "compressed", "uncompressed", "symlink"])
def test_restore_rejects_corrupted_or_ambiguous_inputs_atomically(tmp_path, attack):
    root, directory, manifest, _ = _input_archive(tmp_path)
    if attack == "condition":
        manifest["condition_id"] = "different"
    elif attack == "duplicate":
        manifest["files"].append(copy.deepcopy(manifest["files"][0]))
    elif attack == "compressed":
        (directory / "receptor.gz").write_bytes(b"changed")
    elif attack == "uncompressed":
        manifest["files"][0]["sha256"] = "0" * 64
    elif attack == "symlink":
        archive = directory / "receptor.gz"
        outside = tmp_path / "outside.gz"
        archive.rename(outside)
        archive.symlink_to(outside)
    write_json(directory / "manifest.json", manifest)
    destination = tmp_path / "restored"
    with pytest.raises(ValueError):
        restore_site_inputs(root, "example", "original", destination)
    assert not destination.exists()
    assert not list(tmp_path.glob(".site-input-restore.*"))


def test_only_centers_receive_a_floating_tolerance():
    observed = [{"center": [1., 2., 3.], "score": 0.5, "rank": 1, "residues": ["A:2"]}]
    within = copy.deepcopy(observed)
    within[0]["center"][0] += 1e-11
    _same_predictions(within, observed)
    for key, value in (("center", [1. + 1e-8, 2., 3.]), ("center", [float("nan"), 2., 3.]),
                       ("score", 0.5 + 1e-15), ("rank", 2), ("residues", ["B:2"])):
        changed = copy.deepcopy(observed)
        changed[0][key] = value
        with pytest.raises(ValueError):
            _same_predictions(changed, observed)


@pytest.fixture
def archived_case(tmp_path):
    original, corrected = tmp_path / "original", tmp_path / "corrected"
    restore_site_inputs(ROOT, "154L", "original", original)
    restore_site_inputs(ROOT, "154L", "fpocket_compatibility", corrected)
    config = json.loads((ROOT / "benchmarks/site_experiment.json").read_text())
    case = next(case for case in config["cases"] if case["condition_id"] == "154L")
    evidence = _input_evidence(case, original, corrected)
    return original, corrected, evidence


@pytest.mark.parametrize("directory,method,initial,count", [("fpocket_v2", "fpocket", False, 4),
                                                         ("p2rank", "p2rank", False, 2),
                                                         ("fpocket", "fpocket", True, 0)])
def test_replay_real_archives_after_relocation(tmp_path, archived_case, directory, method, initial, count):
    original, corrected, evidence = archived_case
    source = ROOT / "results/site_discovery/154L" / directory
    portable = tmp_path / "portable" / directory
    portable.mkdir(parents=True)
    record = "execution.json" if method == "fpocket" else "result.json"
    for name in (record, "raw_archive.json", "raw_artifacts.tar.gz"):
        shutil.copy2(source / name, portable / name)
    lock_path = ROOT / ("benchmarks/site_methods.lock.json" if method == "fpocket" else "benchmarks/p2rank.lock.json")
    row = _verify_attempt(portable, tmp_path / "replayed", corrected if directory == "fpocket_v2" else original,
                          method=method, lock=json.loads(lock_path.read_text()), lock_sha=sha256(lock_path),
                          initial=initial, original_sha=evidence["original_sha256"])
    assert row["prediction_count"] == count
    assert row["raw_predictions_reparsed"] is (not initial)
    if initial:
        assert row["returncode"] == 1
        assert "no negative-site biological inference" in row["interpretation"]


@pytest.mark.parametrize("attack", ["score", "center", "rank", "command", "input_count"])
def test_replay_detects_rewritten_normalization_even_after_repacking(tmp_path, archived_case, attack):
    original, _, evidence = archived_case
    source = ROOT / "results/site_discovery/154L/p2rank"
    modified = tmp_path / "modified"
    restore_run(source, modified)
    record_path = modified / "result.json"
    record = json.loads(record_path.read_text())
    if attack == "score":
        record["predictions"][0]["score"] += 1e-12
    elif attack == "center":
        record["predictions"][0]["center"][0] += 1e-7
    elif attack == "rank":
        record["predictions"][0]["rank"] = 2
    elif attack == "command":
        record["command"][-1] = "different_model"
    elif attack == "input_count":
        record["input_atom_count"] -= 1
    write_json(record_path, record)
    pack_run(modified, "result.json")
    lock_path = ROOT / "benchmarks/p2rank.lock.json"
    with pytest.raises(ValueError):
        _verify_attempt(modified, tmp_path / "replayed", original, method="p2rank",
                        lock=json.loads(lock_path.read_text()), lock_sha=sha256(lock_path),
                        original_sha=evidence["original_sha256"])
