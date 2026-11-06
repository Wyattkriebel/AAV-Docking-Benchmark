"""Source-pinned Vina-Carb preparation and execution checks.

The original program uses residue numbers and literal CCD ring-atom names to
recognize glycosidic torsions. Generic atom names or one ligand residue silently
disable that mechanism. These adapters preserve a separate unique identity map.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tarfile
import tempfile
import time

import numpy as np

from .provenance import fetch_verified, sha256, write_json


ROOT = Path(__file__).resolve().parents[1]
HYDROGEN_TYPES = frozenset({"H", "HD", "HS"})


def _atom(line):
    if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 78:
        raise ValueError("Expected a complete PDBQT atom record")
    try:
        serial = int(line[6:11])
        xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        charge = float(line[68:76])
    except ValueError as error:
        raise ValueError("Invalid PDBQT atom serial, coordinates or charge") from error
    if not all(math.isfinite(x) for x in [*xyz, charge]):
        raise ValueError("Nonfinite PDBQT atom")
    return {"serial": serial, "name": line[12:16].strip(), "resname": line[17:20].strip(),
            "resnum": int(line[22:26]), "chain": line[21:22], "coordinates": xyz,
            "charge": charge, "type": line[77:].strip()}


def export_vina_carb_ligand(source_pdbqt, atom_map, destination):
    """Restore CCD identities without changing coordinates, charges or torsions.

    ``atom_map`` must come from independently generated chemical preparation;
    only ``name`` and ``source_key`` identities are read, never native positions.
    Source residues receive unique sequential numbers, even across author chains.
    Hydrogens retain their serials and get an isolated non-sugar residue label.
    """
    source_pdbqt, destination = Path(source_pdbqt), Path(destination)
    if destination.exists():
        raise FileExistsError("Vina-Carb ligand destination already exists")
    entries = list(atom_map)
    if not entries or any(set(("name", "source_key")) - e.keys() for e in entries):
        raise ValueError("Complete source atom identities are required")
    mapping = {e["name"]: e for e in entries}
    if len(mapping) != len(entries):
        raise ValueError("Duplicate generated atom identity")
    keys = [tuple(e["source_key"]) for e in entries]
    if any(len(k) != 5 for k in keys) or len(set(keys)) != len(keys):
        raise ValueError("Duplicate or incomplete source atom identity")
    residues = {key: i + 1 for i, key in enumerate(sorted({key[:4] for key in keys}))}
    if len(residues) > 998:
        raise ValueError("Too many source residues for the legacy PDBQT format")
    output, atoms, heavy_seen, serials = [], [], set(), []
    branches = []
    for line in source_pdbqt.read_text().splitlines():
        if line.startswith(("MODEL", "ENDMDL")):
            raise ValueError("Vina-Carb ligand must be a single input, not pose output")
        if line.startswith("BRANCH"):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError("Invalid torsion-tree branch")
            branches.append(tuple(map(int, fields[1:])))
        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue
        atom = _atom(line)
        serials.append(atom["serial"])
        if atom["type"] in HYDROGEN_TYPES:
            source_key = None
            name, residue_name, residue_number = atom["name"], "HDP", 9999
            if not name or len(name) > 3:
                # Hydrogen identity is carried by the unchanged atom serial.
                name = "H"
        else:
            if atom["name"] not in mapping or atom["name"] in heavy_seen:
                raise ValueError("Missing or duplicated heavy atom in source map")
            heavy_seen.add(atom["name"])
            source_key = tuple(mapping[atom["name"]]["source_key"])
            name, residue_name = source_key[-1], source_key[-2]
            residue_number = residues[source_key[:4]]
            if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9]{0,2}", name):
                raise ValueError("CCD atom name is incompatible with Vina-Carb's three-character parser")
            if not isinstance(residue_name, str) or not 1 <= len(residue_name) <= 3:
                raise ValueError("Invalid source residue name")
        rewritten = (line[:12] + f" {name:<3}" + line[16:17] + f"{residue_name:>3}" + line[20:21]
                     + "A" + f"{residue_number:4d}" + " " + line[27:])
        parsed = _atom(rewritten)
        if any(parsed[field] != atom[field] for field in ("serial", "coordinates", "charge", "type")):
            raise AssertionError("Identity export changed molecular input")
        atoms.append({"serial": atom["serial"], "generated_name": atom["name"],
                      "source_key": list(source_key) if source_key is not None else None,
                      "carb_name": name, "carb_resname": residue_name, "carb_resnum": residue_number,
                      "hydrogen": atom["type"] in HYDROGEN_TYPES,
                      "atom_type": atom["type"], "charge": atom["charge"]})
        output.append(rewritten)
    if heavy_seen != mapping.keys():
        raise ValueError("Source map and ligand have different heavy-atom universes")
    # The upstream parser dereferences ligand_info[serial - 1].
    if serials != list(range(1, len(serials) + 1)):
        raise ValueError("Vina-Carb requires contiguous atom serials in input-record order")
    indexed = {a["serial"]: a for a in atoms}
    glycosidic = []
    for a, b in branches:
        if a not in indexed or b not in indexed:
            raise ValueError("Torsion branch references an absent atom")
        aa, bb = indexed[a], indexed[b]
        if aa["hydrogen"] or bb["hydrogen"]:
            raise ValueError("Hydrogen torsion branches are unsupported")
        if aa["carb_resnum"] != bb["carb_resnum"]:
            if "C1" not in (aa["carb_name"], bb["carb_name"]):
                raise ValueError("Vina-Carb 1.0 does not support this non-C1 glycosidic donor")
            acceptor = bb if aa["carb_name"] == "C1" else aa
            if acceptor["carb_name"] not in {"O2", "O3", "O4", "O6"}:
                raise ValueError("Vina-Carb 1.0 does not support this glycosidic acceptor")
            for endpoint in (aa, bb):
                names = {a["carb_name"] for a in atoms if a["carb_resnum"] == endpoint["carb_resnum"]}
                if not {"C1", "C2", "C3", "C4", "C5", "O5"}.issubset(names):
                    raise ValueError("Vina-Carb requires a complete supported pyranose ring")
            glycosidic.append({"branch_serials": [a, b], "source_keys": [aa["source_key"], bb["source_key"]]})
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(output) + "\n")
    receipt = {"schema_version": 1, "source_sha256": sha256(source_pdbqt), "ligand_sha256": sha256(destination),
               "native_coordinates_accessed": False, "coordinates_charges_types_and_torsion_tree_preserved": True,
               "heavy_atoms": len(entries), "atoms": atoms, "expected_glycosidic_links": glycosidic}
    write_json(destination.with_suffix(destination.suffix + ".map.json"), receipt)
    return receipt


def verify_chi_recognition(log_path, expected_links):
    """Require every expected inter-residue linkage to activate both CHI terms."""
    text = Path(log_path).read_text()
    applied = re.findall(r"CHI energy penalties applied to\s+(\S+)\s+(\d+)\s+and\s+(\S+)\s+(\d+)\s+linkage\.", text)
    observed = [frozenset(((a, int(b)), (c, int(d)))) for a, b, c, d in applied]
    expected = []
    for link in expected_links:
        expected.append(frozenset((tuple(x) for x in link)))
    if "NOT applied" in text or len(observed) != len(expected) or sorted(map(sorted, observed)) != sorted(map(sorted, expected)):
        raise ValueError("Vina-Carb did not activate CHI for every declared glycosidic linkage")
    return {"recognized_links": len(observed), "all_expected_links_fully_supported": True}


def read_carb_poses(path, export_receipt):
    """Parse poses by invariant serials, restoring unique generated atom names."""
    atoms = {a["serial"]: a for a in export_receipt["atoms"]}
    required = sorted(serial for serial, a in atoms.items() if not a["hydrogen"])
    output, current, score = [], None, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("MODEL"):
            if current is not None or int(line.split()[1]) != len(output) + 1:
                raise ValueError("Vina-Carb model ranks must be contiguous")
            current, score = {}, None
        elif line.startswith("REMARK VINA RESULT:"):
            if current is None or score is not None:
                raise ValueError("Duplicate or misplaced pose score")
            score = float(line.split()[3])
            if not math.isfinite(score):
                raise ValueError("Nonfinite Vina-Carb score")
        elif line.startswith(("ATOM  ", "HETATM")):
            if current is None:
                raise ValueError("Atom outside a pose model")
            atom = _atom(line)
            if atom["serial"] not in atoms or atom["serial"] in current:
                raise ValueError("Unexpected or duplicate pose atom serial")
            expected = atoms[atom["serial"]]
            if (atom["name"], atom["resname"], atom["resnum"]) != (expected["carb_name"], expected["carb_resname"], expected["carb_resnum"]):
                raise ValueError("Pose changed the exported atom/residue identity")
            if "atom_type" in expected and (atom["type"], atom["charge"]) != (expected["atom_type"], expected["charge"]):
                raise ValueError("Pose changed atom typing or charges")
            current[atom["serial"]] = atom
        elif line.startswith("ENDMDL"):
            if current is None or score is None or set(current) != atoms.keys():
                raise ValueError("Incomplete Vina-Carb model")
            output.append({"rank": len(output) + 1, "vina_carb_score": score,
                           "names": [atoms[i]["generated_name"] for i in required],
                           "coordinates": [current[i]["coordinates"] for i in required]})
            current, score = None, None
    if current is not None or not output:
        raise ValueError("Missing or truncated Vina-Carb models")
    if any(a["vina_carb_score"] > b["vina_carb_score"] + .002 for a, b in zip(output, output[1:])):
        raise ValueError("Vina-Carb models are not score-ranked")
    return output


def validate_chi_probe(record, *, require_links=True):
    """Check a fixed-coordinate coefficient control, not docking accuracy."""
    fields = ("initial_chi_kcal_mol", "coefficient_0", "coefficient_1", "coefficient_2")
    if any(isinstance(record.get(k), bool) or not isinstance(record.get(k), (int, float)) or not math.isfinite(record[k]) for k in fields):
        raise ValueError("CHI control has missing or nonfinite energies")
    if record["coefficient_0"] != 0 or abs(record["coefficient_2"] - 2 * record["coefficient_1"]) > 1e-9:
        raise ValueError("CHI coefficient control failed")
    if require_links and (record.get("glycosidic_links", 0) <= 0
                          or record.get("fully_supported_links") != record["glycosidic_links"]
                          or record["coefficient_1"] <= record["initial_chi_kcal_mol"]
                          or record.get("torsion_perturbed", -1) < 0):
        raise ValueError("CHI terms are not demonstrably active")
    return {"status": "verified", "scope": "original fixed-coordinate CHI routine, not native-pose accuracy",
            "coefficient_scaling_exact_to_1e-9": True}


def expected_chi_links(export_receipt):
    atoms = {a["serial"]: a for a in export_receipt["atoms"]}
    return [[(atoms[i]["carb_resname"], atoms[i]["carb_resnum"]) for i in link["branch_serials"]]
            for link in export_receipt["expected_glycosidic_links"]]


def export_vina_carb_pose(source_input, atom_map, pose_file, model_rank, destination):
    """Put retained pose coordinates into the fixed source chemical/torsion tree."""
    if isinstance(model_rank, bool) or not isinstance(model_rank, int) or model_rank < 1:
        raise ValueError("Pose rank must be a positive integer")
    source_input, pose_file = Path(source_input), Path(pose_file)
    selected, current, seen = None, None, []
    for line in pose_file.read_text().splitlines():
        if line.startswith("MODEL"):
            if current is not None:
                raise ValueError("Nested pose model")
            current = int(line.split()[1])
            seen.append(current)
            if current == model_rank:
                selected = {}
        elif line.startswith(("ATOM  ", "HETATM")) and current == model_rank:
            atom = _atom(line)
            if atom["serial"] in selected:
                raise ValueError("Duplicate retained pose atom")
            selected[atom["serial"]] = (atom, line)
        elif line.startswith("ENDMDL"):
            if current is None:
                raise ValueError("Unexpected pose terminator")
            current = None
    if current is not None or seen != list(range(1, len(seen) + 1)) or selected is None:
        raise ValueError("Missing, noncontiguous or truncated retained pose model")
    reconstructed, expected_serials = [], []
    for line in source_input.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = _atom(line)
            expected_serials.append(atom["serial"])
            if atom["serial"] not in selected:
                raise ValueError("Retained pose lost a prepared atom")
            moved, pose_line = selected[atom["serial"]]
            if any(atom[key] != moved[key] for key in ("name", "resname", "resnum", "type", "charge")):
                raise ValueError("Retained pose changed prepared atom identity, typing or charges")
            line = line[:30] + pose_line[30:54] + line[54:]
        reconstructed.append(line)
    if set(expected_serials) != selected.keys():
        raise ValueError("Retained pose has extra atoms")
    with tempfile.TemporaryDirectory() as temporary:
        assembled = Path(temporary) / "pose.pdbqt"
        assembled.write_text("\n".join(reconstructed) + "\n")
        receipt = export_vina_carb_ligand(assembled, atom_map, destination)
    receipt.update({"prepared_input_sha256": sha256(source_input), "pose_file_sha256": sha256(pose_file),
                    "model_rank": model_rank, "source": "retained computational pose; not experimental coordinates"})
    write_json(Path(destination).with_suffix(Path(destination).suffix + ".map.json"), receipt)
    return receipt


def validate_installation(installation_dir, *, root=ROOT):
    root, installation_dir = Path(root).resolve(), Path(installation_dir).resolve()
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    installation = json.loads((installation_dir / "installation.json").read_text())
    if (installation["source_sha256"] != lock["source"]["sha256"]
            or installation["boost_sha256"] != lock["boost"]["sha256"]
            or installation["algorithm_source_changes"] != []):
        raise ValueError("Native installation does not match the pinned algorithm source")
    relative = Path(installation["executable"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Native executable must stay inside its isolated installation")
    executable = installation_dir / relative
    if sha256(executable) != installation["executable_sha256"]:
        raise ValueError("Native executable checksum mismatch")
    source = installation_dir / "VC_1_0/autodock_vina_1_1_2"
    for relative, digest in lock["unchanged_algorithm_source_sha256"].items():
        if sha256(source / relative) != digest:
            raise ValueError("Original Vina-Carb scoring/search source changed")
    return executable


def build_fixed_pose_scorer(installation_dir, *, root=ROOT):
    """Link a separate caller of the unchanged source scoring functions."""
    root, installation_dir = Path(root).resolve(), Path(installation_dir).resolve()
    validate_installation(installation_dir, root=root)
    source = installation_dir / "VC_1_0/autodock_vina_1_1_2"
    build, prefix, glylib = source / "build/native-arm64", installation_dir / "boost-native", source / "src/lib/glylib"
    cpp = root / "results/additional_method_reproduction/build/fixed_pose_score.cpp"
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    if sha256(cpp) != lock["fixed_pose_rescorer"]["source_sha256"]:
        raise ValueError("Fixed-pose scoring caller differs from the frozen method lock")
    executable = build / "vc-fixed-pose-score"
    flags = ["clang++", "-std=c++98", "-O3", "-DNDEBUG", "-Wno-deprecated-declarations", "-Wno-writable-strings",
             "-I" + str(prefix / "include"), "-I" + str(glylib / "inc"), "-I" + str(source / "src/lib")]
    _run(flags + ["-o", str(executable), str(cpp)] + [str(p) for p in sorted(build.glob("*.o")) if p.name != "main.o"]
         + ["-L" + str(prefix / "lib"), "-L" + str(glylib / "lib"), "-lboost_system", "-lboost_thread", "-lboost_serialization", "-lboost_filesystem", "-lboost_program_options", "-lglylib", "-lpthread"],
         build, installation_dir / "fixed_pose_build_logs", "link")
    receipt = {"executable": str(executable.relative_to(installation_dir)), "executable_sha256": sha256(executable),
               "caller_source_sha256": sha256(cpp), "method_lock_sha256": sha256(root / "benchmarks/additional_method.lock.json"),
               "underlying_scoring_source_changed": False, "intramolecular_reference": "shared zero; within-fixed-ligand pool ranking only"}
    write_json(installation_dir / "fixed_pose_scorer.json", receipt)
    return executable


def score_fixed_pose(receptor, ligand, export_receipt, installation_dir, *, root=ROOT):
    """Score one computational pose; no search, minimization or native fit occurs."""
    root, installation_dir = Path(root).resolve(), Path(installation_dir).resolve()
    validate_installation(installation_dir, root=root)
    scorer = json.loads((installation_dir / "fixed_pose_scorer.json").read_text())
    executable = installation_dir / scorer["executable"]
    if (sha256(executable) != scorer["executable_sha256"]
            or sha256(ligand) != export_receipt["ligand_sha256"]
            or export_receipt["native_coordinates_accessed"] is not False
            or scorer["method_lock_sha256"] != sha256(root / "benchmarks/additional_method.lock.json")):
        raise ValueError("Fixed-pose scoring caller, exported ligand or method plan changed")
    with tempfile.TemporaryDirectory(prefix="vc-fixed-pose-") as temporary:
        command = [str(executable), str(Path(receptor).resolve()), str(Path(ligand).resolve())]
        result = subprocess.run(command, cwd=temporary, env={**os.environ, "VINA_CARB": str(installation_dir / "VC_1_0")},
                                capture_output=True, text=True, timeout=120)
        if result.returncode:
            raise RuntimeError(f"Vina-Carb fixed-pose scoring failed: {result.stderr[:300]}")
        output = json.loads(result.stdout.strip().splitlines()[-1])
        recognition_path = Path(temporary) / "VC_log.txt"
        recognition = verify_chi_recognition(recognition_path, expected_chi_links(export_receipt))
        recognition_text = recognition_path.read_text()
    expected_count = len(export_receipt["expected_glycosidic_links"])
    if (output["glycosidic_links"] != expected_count or output["fully_supported_links"] != expected_count
            or output["chi_coefficient0"] != 0 or output["nonbonded_cutoff_A"] != 8
            or output["intramolecular_reference"] != 0 or not 0 < output["torsion_normalization_scale"] <= 1
            or any(isinstance(value, float) and not math.isfinite(value) for value in output.values())
            or abs(output["combined_shared0"] - output["base_shared0"] - output["chi_coefficient1"]) > 1e-9):
        raise ValueError("Fixed-pose score/CHI identity or linkage denominator failed")
    return {**output, "recognition": recognition, "command": command, "executable_sha256": sha256(executable),
            "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand),
            "coordinate_optimization": False, "native_coordinates_accessed": False,
            "raw_stdout": result.stdout, "raw_stderr": result.stderr, "recognition_log": recognition_text,
            "raw_stdout_sha256": hashlib.sha256(result.stdout.encode()).hexdigest(),
            "recognition_log_sha256": hashlib.sha256(recognition_text.encode()).hexdigest(),
            "score_role": "VC1.0 original scoring functions with common zero intramolecular reference; rank only within one fixed ligand/torsion-tree pool"}


def run_released_example(installation_dir, destination, *, root=ROOT):
    """Execute the released tutorial configuration without overwriting evidence."""
    root, installation_dir, destination = Path(root).resolve(), Path(installation_dir).resolve(), Path(destination).resolve()
    executable = validate_installation(installation_dir, root=root)
    if destination.exists():
        raise FileExistsError("Tutorial replay requires a fresh output directory")
    inventory = json.loads((root / "results/additional_method_reproduction/example_inputs/manifest.json").read_text())
    sources = {entry["filename"]: entry for entry in inventory["files"]}
    destination.mkdir(parents=True)
    for name in ("config_vc.txt", "receptor.pdbqt", "LeY-xray.pdbqt"):
        original = root / "results/additional_method_reproduction/example_inputs" / name
        if sha256(original) != sources[name]["sha256"]:
            raise ValueError("Tutorial replay source bytes changed")
        shutil.copyfile(original, destination / name)
    command = [str(executable), "--config", "config_vc.txt", "--log", "vina-carb.log"]
    record = {"command": command, "status": "running", "configuration_modified": False,
              "information_regime": "released known-site tutorial; native-like initial conformer; software validation only",
              "input_sha256": {name: sources[name]["sha256"] for name in ("config_vc.txt", "receptor.pdbqt", "LeY-xray.pdbqt")},
              "executable_sha256": sha256(executable)}
    write_json(destination / "execution.json", record)
    start = time.perf_counter()
    with (destination / "stdout.log").open("w") as out, (destination / "stderr.log").open("w") as err:
        completed = subprocess.run(command, cwd=destination, env={**os.environ, "VINA_CARB": str(installation_dir / "VC_1_0")},
                                   stdout=out, stderr=err, timeout=1800)
    record.update(status="completed" if completed.returncode == 0 else "failed", exit_code=completed.returncode,
                  wall_seconds=time.perf_counter() - start,
                  artifact_sha256={p.name: sha256(p) for p in destination.iterdir() if p.is_file() and p.name != "execution.json"})
    write_json(destination / "execution.json", record)
    return record


def verify_released_example(root=ROOT):
    """Verify retained tutorial output and an independent fixed-coordinate CHI control."""
    root = Path(root).resolve()
    base = root / "results/additional_method_reproduction"
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    inventory = json.loads((base / "example_inputs/manifest.json").read_text())
    if inventory["source_commit"] != lock["tutorial"]["commit"]:
        raise ValueError("Tutorial source commit differs")
    for entry in inventory["files"]:
        if Path(entry["filename"]).name != entry["filename"] or sha256(base / "example_inputs" / entry["filename"]) != entry["sha256"]:
            raise ValueError("Tutorial source checksum mismatch")
    record = json.loads((base / "example_exact/execution.json").read_text())
    if record["status"] != "completed" or record["exit_code"] != 0 or record["configuration_modified"] is not False:
        raise ValueError("Published tutorial workflow did not finish unchanged")
    for filename, digest in record["artifact_sha256"].items():
        if Path(filename).name != filename or sha256(base / "example_exact" / filename) != digest:
            raise ValueError("Tutorial execution artifact checksum mismatch")
    input_atoms = [_atom(line) for line in (base / "example_inputs/LeY-xray.pdbqt").read_text().splitlines()
                   if line.startswith(("ATOM  ", "HETATM"))]
    receipt = {"atoms": [{"serial": a["serial"], "generated_name": f"{a['resnum']}:{a['resname']}:{a['name']}",
                          "carb_name": a["name"], "carb_resname": a["resname"], "carb_resnum": a["resnum"],
                          "hydrogen": a["type"] in HYDROGEN_TYPES} for a in input_atoms]}
    poses = read_carb_poses(base / "example_exact/LeY-xray_vc_out.pdbqt", receipt)
    published = read_carb_poses(base / "example_inputs/LeY-xray_vc_out.pdbqt", receipt)
    expected = [[("NDG", 1), ("GAL", 2)], [("GAL", 2), ("FUC", 3)], [("NDG", 1), ("FUC", 4)]]
    recognition = verify_chi_recognition(base / "example_exact/VC_log.txt", expected)
    probe = json.loads((base / "chi_control/stdout.log").read_text().strip().splitlines()[-1])
    probe_check = validate_chi_probe(probe)
    if len(poses) != 10 or len(published) != 10 or len(poses[0]["coordinates"]) != 46:
        raise ValueError("Tutorial pose or atom denominator differs")
    actual_seed = re.findall(r"Using random seed:\s*(-?\d+)", (base / "example_exact/stdout.log").read_text())
    if len(actual_seed) != 1:
        raise ValueError("Actual random seed was not recorded")
    return {"schema_version": 1, "status": "verified", "method": "Vina-Carb 1.0", "poses": len(poses),
            "heavy_atoms_per_pose": len(poses[0]["coordinates"]), "chi_recognition": recognition,
            "fixed_coordinate_chi_control": {**probe, **probe_check}, "actual_random_seed": int(actual_seed[0]),
            "replayed_top1_score_kcal_mol": poses[0]["vina_carb_score"],
            "published_top1_score_kcal_mol": published[0]["vina_carb_score"],
            "released_artifact_coordinate_rmsd_by_rank_A": [float(np.sqrt(np.mean(np.sum((np.asarray(a['coordinates']) - np.asarray(b['coordinates'])) ** 2, axis=1)))) for a, b in zip(poses, published)],
            "published_coordinates_reproduced": all(a["coordinates"] == b["coordinates"] for a, b in zip(poses, published)),
            "interpretation": "Execution of the published tutorial, atom completeness and active CHI mechanism verified. Explicit seed0 is fixed in this original version; cross-build coordinate equality and paper-level accuracy are not claimed. Reported RMSDs compare computational artifacts, not experimental pose truth.",
            "source_execution_sha256": sha256(base / "example_exact/execution.json"),
            "method_lock_sha256": sha256(root / "benchmarks/additional_method.lock.json"),
            "verifier_sha256": sha256(Path(__file__))}


def verify_fixed_pose_validation(root=ROOT):
    """Reparse exact native scorer output and recompute the tutorial controls."""
    root = Path(root).resolve()
    base = root / "results/additional_method_reproduction"
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    if sha256(base / "build/fixed_pose_score.cpp") != lock["fixed_pose_rescorer"]["source_sha256"]:
        raise ValueError("Native scorer caller source changed")
    summary = json.loads((base / "fixed_pose_validation_v2/summary.json").read_text())
    rows = summary["rows"]
    if [row["model_rank"] for row in rows] != list(range(1, 11)):
        raise ValueError("Tutorial scorer validation denominator changed")
    offsets, scales = [], []
    for row in rows:
        if hashlib.sha256(row["raw_stdout"].encode()).hexdigest() != row["raw_stdout_sha256"]:
            raise ValueError("Raw scorer output checksum changed")
        if hashlib.sha256(row["recognition_log"].encode()).hexdigest() != row["recognition_log_sha256"]:
            raise ValueError("Raw CHI recognition output changed")
        raw = json.loads(row["raw_stdout"].strip().splitlines()[-1])
        if any(row[key] != value for key, value in raw.items()):
            raise ValueError("Parsed native score differs from raw output")
        if (raw["glycosidic_links"] != 3 or raw["fully_supported_links"] != 3
                or raw["chi_coefficient0"] != 0 or raw["intramolecular_reference"] != 0
                or raw["nonbonded_cutoff_A"] != 8
                or abs(raw["combined_shared0"] - raw["base_shared0"] - raw["chi_coefficient1"]) > 1e-9):
            raise ValueError("Original fixed-pose scorer control changed")
        offsets.append(raw["combined_shared0"] - row["reported_fullsearch_score"])
        scales.append(raw["torsion_normalization_scale"])
    if max(scales) - min(scales) > 1e-12 or not 0 < min(scales) <= 1:
        raise ValueError("Shared-reference validation mixes torsion normalization contexts")
    order = [row["model_rank"] for row in sorted(rows, key=lambda row: row["combined_shared0"])]
    if order != list(range(1, 11)) or max(offsets) - min(offsets) > .1:
        raise ValueError("Fixed-pose score does not reproduce source ranks/common offset to printed precision")
    return {"status": "verified", "poses": 10, "source_fullsearch_rank_order_preserved": True,
            "offset_range_kcal_mol": [min(offsets), max(offsets)], "fullsearch_printed_score_unit_kcal_mol": .1,
            "common_reference_offset_admits_all_printed_score_intervals": True,
            "torsion_normalization_scale": scales[0], "source_sha256": sha256(base / "fixed_pose_validation_v2/summary.json"),
            "scope": "Retained tutorial poses and unchanged scoring objects; within-fixed-ligand ranking only, no experimental pose adjudication"}


def _extract(archive, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as source:
        for member in source.getmembers():
            p = Path(member.name)
            if p.is_absolute() or ".." in p.parts or member.isdev():
                raise ValueError("Unsafe source archive member")
        source.extractall(destination, filter="data")


def _run(command, cwd, log_directory, name):
    log_directory.mkdir(parents=True, exist_ok=True)
    with (log_directory / (name + ".stdout.log")).open("w") as stdout, (log_directory / (name + ".stderr.log")).open("w") as stderr:
        result = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, timeout=3600)
    write_json(log_directory / (name + ".command.json"), {"command": list(map(str, command)), "exit_code": result.returncode})
    if result.returncode:
        raise RuntimeError(f"Native build step failed: {name}; see {log_directory}")


def install_vina_carb(destination, *, root=ROOT, archives_dir=None):
    """Compile the checksum-pinned original release in an isolated fresh directory."""
    root, destination = Path(root).resolve(), Path(destination).resolve()
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise ValueError("This audited build recipe supports macOS arm64 only")
    if destination.exists():
        raise FileExistsError("Native installation requires a fresh destination")
    lock = json.loads((root / "benchmarks/additional_method.lock.json").read_text())
    destination.mkdir(parents=True)
    logs = destination / "build_logs"
    source_tar = destination / "vina-carb_1.0.tgz"
    boost_tar = destination / "boost_1_75_0.tar.bz2"
    for archive, entry in ((source_tar, lock["source"]), (boost_tar, lock["boost"])):
        if archives_dir is None:
            fetch_verified(entry["url"], archive, entry["sha256"])
        else:
            cached = Path(archives_dir) / archive.name
            if sha256(cached) != entry["sha256"]:
                raise ValueError("Cached native source archive checksum mismatch")
            shutil.copyfile(cached, archive)
            if sha256(archive) != entry["sha256"]:
                raise ValueError("Copied native source archive checksum mismatch")
        _extract(archive, destination)
    source = destination / "VC_1_0/autodock_vina_1_1_2"
    for relative, changes in lock["portability_replacements"].items():
        path = source / relative
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != lock["original_patched_file_sha256"][relative]:
            raise ValueError("Portability patch source differs from pinned original")
        for before, after in changes:
            if before.encode() not in raw:
                raise ValueError("Portability patch precondition failed")
            raw = raw.replace(before.encode(), after.encode())
        path.write_bytes(raw)
    boost = destination / "boost_1_75_0"
    prefix = destination / "boost-native"
    _run(["./bootstrap.sh", "--with-toolset=clang", "--with-libraries=system,thread,serialization,filesystem,program_options"], boost, logs, "boost-bootstrap")
    _run(["./b2", "-j2", "toolset=clang", "link=static", "runtime-link=shared", "threading=multi", "variant=release", "cxxflags=-std=c++98", "--prefix=" + str(prefix), "install"], boost, logs, "boost-build")
    glylib = source / "src/lib/glylib"
    (glylib / "lib/libglylib.a").unlink(missing_ok=True)
    _run(["make", "justC"], glylib / "lib", logs, "glylib")
    build = source / "build/native-arm64"
    build.mkdir()
    flags = ["clang++", "-std=c++98", "-O3", "-DNDEBUG", "-Wno-deprecated-declarations", "-Wno-writable-strings",
             "-I" + str(prefix / "include"), "-I" + str(glylib / "inc"), "-I" + str(source / "src/lib")]
    files = sorted((source / "src/lib").glob("*.cpp")) + [source / "src/main/main.cpp"]
    def compile_one(path):
        _run(flags + ["-c", str(path), "-o", str(build / (path.stem + ".o"))], build, logs, path.stem)
    with ThreadPoolExecutor(max_workers=2) as workers:
        list(workers.map(compile_one, files))
    executable = build / "vina-carb"
    _run(flags + ["-o", str(executable)] + [str(build / (p.stem + ".o")) for p in files]
         + ["-L" + str(prefix / "lib"), "-L" + str(glylib / "lib"), "-lboost_system", "-lboost_thread", "-lboost_serialization", "-lboost_filesystem", "-lboost_program_options", "-lglylib", "-lpthread"], build, logs, "link")
    version = subprocess.check_output([str(executable), "--version"], text=True).strip()
    if version != lock["version_banner"]:
        raise ValueError("Unexpected Vina-Carb version banner")
    installation = {"source_sha256": lock["source"]["sha256"], "boost_sha256": lock["boost"]["sha256"],
                    "method_lock_sha256": sha256(root / "benchmarks/additional_method.lock.json"),
                    "executable": str(executable.relative_to(destination)), "executable_sha256": sha256(executable),
                    "version": version, "platform": platform.platform(),
                    "algorithm_source_changes": [], "portability_patch_files": sorted(lock["portability_replacements"]),
                    "build_log_sha256": {p.name: sha256(p) for p in sorted(logs.iterdir()) if p.is_file()}}
    write_json(destination / "installation.json", installation)
    return executable


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["install", "build-scorer", "replay-example", "verify-example"])
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--archives-dir", type=Path, help="optional directory containing the two exact source archives")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--installation", type=Path)
    args = parser.parse_args()
    if args.stage == "install":
        if args.destination is None:
            parser.error("install requires --destination")
        print(install_vina_carb(args.destination, archives_dir=args.archives_dir))
    elif args.stage == "build-scorer":
        if args.installation is None:
            parser.error("build-scorer requires --installation")
        print(build_fixed_pose_scorer(args.installation))
    elif args.stage == "replay-example":
        if args.installation is None or args.destination is None:
            parser.error("replay-example requires --installation and --destination")
        print(json.dumps(run_released_example(args.installation, args.destination), indent=2, sort_keys=True))
    else:
        verification = verify_released_example()
        verification["fixed_pose_validation"] = verify_fixed_pose_validation()
        if args.output:
            write_json(args.output, verification)
        print(json.dumps(verification, indent=2, sort_keys=True))
