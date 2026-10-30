"""Replay one frozen AAV9 diagnostic with independent pose/contact arithmetic.

Run from the repository root with its installed dependencies:
    .venv/bin/python scripts/reproduce_aav9_one_seed.py

This isolated check never changes the primary three-seed denominator. It reuses
checksummed prepared inputs, rather than independently validating preparation.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import resource
import subprocess
import time

import gemmi
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def pose_models(path, required):
    """Independent fixed-column parser; atom identity is an exact bijection."""
    result, atoms, score, rank = [], None, None, None
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if atoms is not None:
                raise ValueError("Unterminated independent pose model")
            rank, atoms, score = int(line.split()[1]), {}, None
            if rank != len(result) + 1:
                raise ValueError("Noncontiguous independent pose ranks")
        elif line.startswith("REMARK VINA RESULT:"):
            if atoms is None or score is not None:
                raise ValueError("Missing model or duplicated score")
            score = float(line.split()[3])
        elif line.startswith(("ATOM  ", "HETATM")):
            if atoms is None:
                raise ValueError("Independent pose atom outside MODEL")
            if line.split()[-1] in ("H", "HD", "HS"):
                continue
            name = line[12:16].strip()
            if name in atoms:
                raise ValueError("Duplicate independent pose atom")
            atoms[name] = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        elif line.startswith("ENDMDL"):
            if atoms is None or score is None or set(atoms) != set(required):
                raise ValueError("Incomplete independent pose atom map/score")
            if not all(math.isfinite(x) for xyz in atoms.values() for x in xyz) or not math.isfinite(score):
                raise ValueError("Nonfinite independent pose data")
            result.append({"rank": rank, "score": score, "coordinates": [atoms[name] for name in required]})
            atoms = None
    if atoms is not None or not result:
        raise ValueError("Empty/truncated independent pose artifact")
    return result


def rmsd(a, b):
    if len(a) != len(b) or not a:
        raise ValueError("Independent RMSD requires equal, nonempty atom maps")
    return math.sqrt(sum(sum((x - y)**2 for x, y in zip(pa, pb)) for pa, pb in zip(a, b)) / len(a))


def contacts(ligand, receptor, cutoff=4.0):
    return {(i, j) for i, pa in enumerate(ligand) for j, pb in enumerate(receptor)
            if sum((x - y)**2 for x, y in zip(pa, pb)) < cutoff**2}


def deposited_preparation_controls(preparation, reference):
    """Check the native atom map and recorded apo transform against source CIFs."""
    raw = ROOT / "data/raw"
    source_hashes = {name: digest(raw / name) for name in preparation["input_sha256"]}
    if source_hashes != preparation["input_sha256"]:
        raise ValueError("Independent deposited source checksums changed")
    holo = gemmi.read_structure(str(raw / "7MTZ.cif"))
    apo = gemmi.read_structure(str(raw / "7MT0.cif"))
    ligand_ref = reference["reference_ligand"]
    deposited = next(residue for chain in holo[0] if chain.name == ligand_ref["author_chain"]
                     for residue in chain if residue.name == "GAL" and residue.seqid.num == ligand_ref["author_residue"])
    heavy = [(atom.name, list(atom.pos)) for atom in deposited if not atom.element.is_hydrogen]
    deposited_atoms = dict(heavy)
    if len(heavy) != len(deposited_atoms) or set(deposited_atoms) != set(reference["names"]):
        raise ValueError("Independent deposited native atom map differs")
    if rmsd([deposited_atoms[name] for name in reference["names"]], reference["coordinates"]) != 0:
        raise ValueError("Independent deposited native coordinates differ")

    def ca_atoms(structure):
        atoms = {}
        for chain in structure[0]:
            for residue in chain:
                if residue.het_flag != "A":
                    continue
                for atom in residue:
                    if atom.name == "CA":
                        key = (chain.name, residue.seqid.num, residue.seqid.icode, residue.name)
                        if key in atoms:
                            raise ValueError("Independent CA atom map contains a duplicate")
                        atoms[key] = list(atom.pos)
        return atoms

    holo_ca, apo_ca = ca_atoms(holo), ca_atoms(apo)
    matched = sorted(holo_ca.keys() & apo_ca.keys())
    rotation = np.asarray(preparation["rotation"])
    translation = np.asarray(preparation["translation"])
    moved = np.asarray([apo_ca[key] for key in matched]) @ rotation + translation
    target = np.asarray([holo_ca[key] for key in matched])
    ca_rmsd = float(np.sqrt(np.mean(np.sum((moved - target)**2, axis=1))))
    if len(matched) != preparation["alignment_CA_count"] or not math.isclose(ca_rmsd, preparation["apo_to_holo_CA_rmsd_A"], abs_tol=1e-10):
        raise ValueError("Independent source-CA transform check differs")
    return {"deposited_source_sha256": source_hashes, "deposited_native_ligand_atom_names_and_coordinates_exact": True,
            "source_apo_to_holo_CA_count": len(matched), "source_apo_to_holo_CA_rmsd_A": ca_rmsd,
            "source_CA_rmsd_absolute_difference_from_preparation_A": abs(ca_rmsd - preparation["apo_to_holo_CA_rmsd_A"])}


def main():
    prepared = ROOT / "data/processed/aav9"
    primary = ROOT / "results/development/aav9"
    output = ROOT / "results/independent_reproduction/aav9_one_seed"
    output.mkdir(parents=True, exist_ok=True)
    pose_path = output / "holo_complete_104729.pdbqt"
    if pose_path.exists():
        raise FileExistsError("Independent artifact already exists; preserve it before rerunning")
    preparation = json.loads((prepared / "preparation.json").read_text())
    reference = json.loads((prepared / "reference.json").read_text())
    config = json.loads((ROOT / "benchmarks/development.json").read_text())
    config["seeds"] = [104729]
    write_json(output / "config.json", config)
    selected = next(row for row in json.loads((primary / "summary.json").read_text())
                    if row["receptor"] == "holo_complete" and row["seed"] == 104729)
    inputs = {
        "executable_sha256": ROOT / "tools/vina-1.2.7",
        "receptor_sha256": prepared / "holo_complete_receptor.pdbqt",
        "receptor_pdb_sha256": prepared / "holo_complete_receptor.pdb",
        "ligand_sha256": prepared / "galactose.pdbqt",
        "preparation_sha256": prepared / "preparation.json",
        "reference_sha256": prepared / "reference.json",
        "full_assembly_sha256": prepared / reference["full_assembly_file"],
    }
    fingerprints = {key: digest(path) for key, path in inputs.items()}
    if any(selected["fingerprint"][key] != value for key, value in fingerprints.items()):
        raise ValueError("Prepared input differs from the primary run")
    executable = ROOT / "tools/vina-1.2.7"
    version = subprocess.check_output([str(executable), "--version"], text=True).strip()
    if version != "AutoDock Vina v1.2.7":
        raise ValueError("Unexpected independent executable version")
    command = [str(executable), "--receptor", str(prepared / "holo_complete_receptor.pdbqt"),
               "--ligand", str(prepared / "galactose.pdbqt"), "--out", str(pose_path),
               "--scoring", "vina", "--cpu", str(config["cpu"]), "--seed", "104729",
               "--exhaustiveness", str(config["exhaustiveness"]), "--num_modes", str(config["num_modes"]),
               "--min_rmsd", str(config["min_rmsd_A"]), "--energy_range", str(config["energy_range_kcal_mol"])]
    for axis, center, size in zip("xyz", preparation["box_center"], preparation["box_size"]):
        command.extend(["--center_" + axis, str(center), "--size_" + axis, str(size)])
    primary_comparison = list(selected["command"])
    primary_comparison[primary_comparison.index("--out") + 1] = str(pose_path)
    if command != primary_comparison:
        raise ValueError("Independent command differs beyond output destination")
    started = datetime.now(timezone.utc).isoformat()
    start, usage = time.perf_counter(), resource.getrusage(resource.RUSAGE_CHILDREN)
    with (output / "holo_complete_104729.log").open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
    elapsed = time.perf_counter() - start
    end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    if completed.returncode:
        raise RuntimeError("Independent Vina invocation failed; inspect archived log")
    names, native = reference["names"], reference["coordinates"]
    if len(names) != len(set(names)) or len(names) != len(native):
        raise ValueError("Independent reference atom map failed")
    models = pose_models(pose_path, names)
    primary_models = pose_models(primary / selected["pose_file"], names)
    if len(models) != len(primary_models):
        raise ValueError("Independent returned-pose count differs")
    rows = [{"rank": model["rank"], "vina_score": model["score"], "rmsd_A": rmsd(model["coordinates"], native)}
            for model in models]
    primary_rows = [{"rank": model["rank"], "vina_score": model["score"], "rmsd_A": rmsd(model["coordinates"], native)}
                    for model in primary_models]
    max_coordinate_difference = max(abs(x - y) for a, b in zip(models, primary_models)
                                    for pa, pb in zip(a["coordinates"], b["coordinates"]) for x, y in zip(pa, pb))
    max_score_difference = max(abs(a["vina_score"] - b["vina_score"]) for a, b in zip(rows, primary_rows))
    max_rmsd_difference = max(abs(a["rmsd_A"] - b["rmsd_A"]) for a, b in zip(rows, primary_rows))
    receptor = reference["receptor_coordinates"]
    native_pairs = contacts(native, receptor)
    archived_pairs = set(map(tuple, reference["native_contact_pairs_4A"]))
    if native_pairs != archived_pairs:
        raise ValueError("Independent native contacts differ from prepared reference")
    native_recall = len(native_pairs & native_pairs) / len(native_pairs)
    predicted = contacts(models[0]["coordinates"], receptor)
    top1_recall = len(predicted & native_pairs) / len(native_pairs)
    top1_precision = len(predicted & native_pairs) / len(predicted) if predicted else 0.
    translated_native = [[x + 3, y + 4, z] for x, y, z in native]
    translated_error = rmsd(translated_native, native)
    # An explicit common receptor-frame rigid transform must preserve pose error
    # and native contacts; no ligand-only fit is used anywhere in this script.
    def transform(points):
        return [[-y + 17, x - 9, z + 5] for x, y, z in points]
    transformed_error = rmsd(transform(models[0]["coordinates"]), transform(native))
    transformed_native_pairs = contacts(transform(native), transform(receptor))
    rotation = np.asarray(preparation["rotation"])
    rotation_orthogonality_error = float(np.max(np.abs(rotation.T @ rotation - np.eye(3))))
    rotation_determinant = float(np.linalg.det(rotation))
    if not math.isclose(rotation_determinant, 1, abs_tol=1e-10) or rotation_orthogonality_error > 1e-10:
        raise ValueError("Archived apo-to-holo transform is not a proper rigid rotation")
    if transformed_native_pairs != native_pairs or not math.isclose(transformed_error, rows[0]["rmsd_A"], abs_tol=1e-10):
        raise ValueError("Independent common receptor-frame transform control failed")
    if rmsd(native, native) != 0 or not math.isclose(translated_error, 5, abs_tol=1e-12):
        raise ValueError("Independent native/translation RMSD controls failed")
    source_controls = deposited_preparation_controls(preparation, reference)
    summary_fields = {
        "returned_poses": len(rows), "top1_score": rows[0]["vina_score"], "top1_rmsd_A": rows[0]["rmsd_A"],
        "top5_min_rmsd_A": min(row["rmsd_A"] for row in rows[:5]),
        "best_returned_rmsd_A": min(row["rmsd_A"] for row in rows),
        "best_returned_rank": min(rows, key=lambda row: row["rmsd_A"])["rank"],
        "top1_contact_recall": top1_recall, "top1_contact_precision": top1_precision,
    }
    summary_differences = {key: abs(value - selected[key]) for key, value in summary_fields.items()}
    with (output / "independent_pose_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "scope": "Independent agent replay of one prepared AAV9 oracle-site Vina diagnostic on the same host; not an external laboratory or fresh preparation",
        "status": "passed" if max_score_difference == max_rmsd_difference == max_coordinate_difference == 0 and max(summary_differences.values()) < 1e-10 else "mismatch",
        "started_utc": started, "case_id": "AAV9", "receptor": "holo_complete", "seed": 104729,
        "command": command, "vina_version": version, "platform": platform.platform(),
        "wall_seconds": elapsed,
        "cpu_seconds": end_usage.ru_utime + end_usage.ru_stime - usage.ru_utime - usage.ru_stime,
        "input_sha256": fingerprints, "independent_config_sha256": digest(output / "config.json"),
        "analysis_script_sha256": digest(Path(__file__)),
        "primary_config_sha256": selected["fingerprint"]["config_sha256"],
        "pose_sha256": digest(pose_path), "primary_pose_sha256": selected["pose_sha256"],
        "primary_summary_sha256": digest(primary / "summary.json"),
        "summary": summary_fields, "absolute_differences_from_primary_summary": summary_differences,
        "all_rank_comparison": {"max_score_difference_kcal_mol": max_score_difference,
                                "max_rmsd_difference_A": max_rmsd_difference, "max_coordinate_difference_A": max_coordinate_difference},
        "controls": {"native_rmsd_A": 0., "translation_3_4_0_rmsd_A": translated_error,
                     "native_contact_pairs_4A": len(native_pairs), "native_contact_recall": native_recall,
                     "native_contact_precision": 1., "native_contact_identity_matches_preparation": native_pairs == archived_pairs,
                     "common_receptor_frame_transform_rmsd_A": transformed_error,
                     "common_frame_transform_preserves_native_contacts": transformed_native_pairs == native_pairs,
                     **source_controls,
                     "apo_to_holo_rotation_determinant": rotation_determinant,
                     "apo_to_holo_rotation_max_orthogonality_error": rotation_orthogonality_error,
                     "atom_map_heavy_atoms": len(names)},
        "limitations": ["Shares prepared inputs, executable, host and reference with primary run",
                        "Does not independently adjudicate local density or reference chemistry",
                        "Does not validate six development systems, held-out accuracy, blind site recovery or specificity",
                        "Single deterministic seed replay is not an additional biological replicate"],
    }
    write_json(output / "reproduction.json", result)
    print(json.dumps({"status": result["status"], **summary_fields, "all_rank_comparison": result["all_rank_comparison"]}, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
