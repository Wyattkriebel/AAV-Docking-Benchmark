"""Portable fixed inputs and semantic tamper controls; no native outcomes."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import socket
import subprocess
import sys

import pytest

from scripts.verify_serialization_archives import read_archive


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepared_input_verifier", ROOT / "scripts/verify_prepared_inputs.py")
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


@pytest.fixture(scope="module")
def inputs():
    manifest = json.loads((ROOT / verifier.MANIFEST).read_bytes())
    directory = (ROOT / verifier.MANIFEST).parent
    files = read_archive((directory / manifest["archive"]).read_bytes(), manifest, expected_count=698)
    profiles = {}
    for profile, records in zip(verifier.PROFILE_CASES, manifest["source_archives"], strict=True):
        record = json.loads((ROOT / records["manifest"]["path"]).read_bytes())
        profiles[profile] = read_archive((ROOT / records["archive"]["path"]).read_bytes(), record, expected_count=len(record["files"]))
    return files, profiles, manifest


def first_condition(files, manifest):
    index_pin = manifest["indexes"][0]
    index = json.loads(files[index_pin["path"]])
    row = index["conditions"][0]
    return index_pin, index, row, json.loads(files[row["condition_path"]])


def reseal_condition(files, index_pin, index, row, condition):
    files[row["condition_path"]] = encoded(condition)
    row["condition_sha256"] = verifier.sha(files[row["condition_path"]])
    files[index_pin["path"]] = encoded(index)
    index_pin["sha256"] = verifier.sha(files[index_pin["path"]])


@pytest.mark.parametrize("mutation", ["center", "original_candidates", "ligand", "receptor_atom_loss", "omit_failure", "config_pin", "profile_pin", "profile_pin_removed", "ledger_pin", "reference_coordinates"])
def test_internally_resealed_inputs_cannot_move_candidates_or_lose_atoms_or_failures(inputs, mutation):
    base_files, profiles, base_manifest = inputs
    files, manifest = dict(base_files), copy.deepcopy(base_manifest)
    index_pin, index, row, condition = first_condition(files, manifest)
    box = Path(row["condition_path"]).parent / "box01"
    if mutation == "center":
        condition["jobs"][0]["center_A"][0] += 1
        reseal_condition(files, index_pin, index, row, condition)
    elif mutation == "original_candidates":
        condition["candidate_selection"]["candidates"][0]["center_A"][0] += 1
        reseal_condition(files, index_pin, index, row, condition)
    elif mutation == "ligand":
        files[(box / "ligand.pdbqt").as_posix()] += b"\n"
    elif mutation == "receptor_atom_loss":
        path = (box / "neighborhood.pdbqt").as_posix()
        lines = files[path].decode().splitlines(keepends=True)
        i = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM")) and not line[77:].strip().startswith("H"))
        del lines[i]
        files[path] = "".join(lines).encode()
    elif mutation in ("config_pin", "profile_pin", "ledger_pin"):
        key = {"config_pin": "configuration_sha256", "profile_pin": "source_profile_manifest_sha256", "ledger_pin": "source_ledger_sha256"}[mutation]
        condition[key] = "0" * 64
        reseal_condition(files, index_pin, index, row, condition)
    elif mutation == "profile_pin_removed":
        del condition["source_profile_manifest_sha256"]
        reseal_condition(files, index_pin, index, row, condition)
    elif mutation == "reference_coordinates":
        path = (box / "reference.json").as_posix()
        reference = json.loads(files[path]); reference["coordinates"] = [[0, 0, 0]]
        files[path] = encoded(reference)
    else:
        i = next(i for i, r in enumerate(index["conditions"]) if r["status"] == "preparation_failed")
        del index["conditions"][i]
        index_pin["conditions"] -= 1
        files[index_pin["path"]] = encoded(index)
        index_pin["sha256"] = verifier.sha(files[index_pin["path"]])
    with pytest.raises(ValueError):
        verifier.verify_contents(files, profiles, manifest)


@pytest.mark.parametrize("mutation", ["chemical_parent", "observed_map", "omitted_assembly"])
def test_source_profiles_cannot_replace_shared_parent_or_observed_assembly(inputs, mutation):
    files, original, manifest = inputs
    profiles = {k: dict(v) for k, v in original.items()}
    if mutation == "chemical_parent":
        profiles["repair"]["4R9F/ligand_atom_map.json"] += b"\n"
    elif mutation == "observed_map":
        profiles["repair"]["2RDK/protein_atom_map.json"] += b"\n"
    else:
        profiles["strict_monomer"]["3ZWE/full_assembly_receptor.cif"] = profiles["strict_monomer"]["3ZWE/receptor.cif"]
    with pytest.raises(ValueError):
        verifier.verify_contents(files, profiles, manifest)


def package_closure(root, names):
    """Copy only statically imported local library sources, including lazy imports."""
    pending, seen = list(names), set()
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        source = ROOT / name
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        for node in ast.walk(ast.parse(source.read_text())):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                dependency = "aav_benchmark/" + node.module.split(".")[0] + ".py"
                if (ROOT / dependency).is_file():
                    pending.append(dependency)
    return seen


def test_public_input_bundle_requires_no_ignored_tree_network_or_native_calls(inputs, tmp_path, monkeypatch):
    _, _, manifest = inputs
    names = {verifier.MANIFEST, str(Path(verifier.MANIFEST).parent / manifest["archive"]),
             "scripts/verify_prepared_inputs.py", "scripts/verify_serialization_archives.py", "aav_benchmark/__init__.py", "benchmarks/heldout_reference_ledger.json"}
    names.update(r["path"] for r in manifest["frozen_documents"])
    for records in manifest["source_archives"]:
        names.update((records["manifest"]["path"], records["archive"]["path"]))
    libraries = package_closure(tmp_path, manifest["frozen_producers"])
    assert not any(p in libraries for p in ("aav_benchmark/qualified_pose.py", "aav_benchmark/heldout_pose_report.py",
                                            "aav_benchmark/paired_method_report.py", "aav_benchmark/global_pose_report_gate.py"))
    for name in names:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, path)
    before = {p.relative_to(tmp_path).as_posix(): verifier.sha(p.read_bytes()) for p in tmp_path.rglob("*") if p.is_file()}
    def forbidden(*args, **kwargs):
        pytest.fail("Input-only replay reached network or a native subprocess")
    with monkeypatch.context() as guard:
        guard.setattr(socket, "create_connection", forbidden)
        guard.setattr(subprocess, "run", forbidden)
        guard.setattr(subprocess, "Popen", forbidden)
        result = verifier.verify(tmp_path)
    assert result == {"conditions": 48, "full_assembly_conditions": 36, "monomer_conditions": 12,
        "prepared_boxes": 28, "preparation_failed_conditions": 18, "no_candidate_conditions": 2,
        "original_predictor_streams": 16, "source_profiles": 4, "archived_input_files": 698}
    assert not (tmp_path / "tools").exists() and not (tmp_path / "data/processed").exists()
    assert before == {p.relative_to(tmp_path).as_posix(): verifier.sha(p.read_bytes()) for p in tmp_path.rglob("*") if p.is_file()}
    completed = subprocess.run([sys.executable, "scripts/verify_prepared_inputs.py"], cwd=tmp_path,
                               capture_output=True, text=True, check=True)
    assert json.loads(completed.stdout) == result
    path = tmp_path / "results/structural_pose_v3/strict_sources.tar.gz"
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="SHA256"):
        verifier.verify(tmp_path)
