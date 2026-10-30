"""Execute frozen Vina diagnostics and evaluate every retained candidate."""

from __future__ import annotations

import csv
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import time

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .metrics import contact_metrics, contact_pairs, ligand_rmsd, mapped_coordinates, summarize_poses
from .provenance import sha256, write_json


def read_poses(path, required_names):
    """Strict PDBQT parser: require complete models, scores and heavy atom maps."""
    models, current, score, model_number = [], None, None, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("MODEL"):
            if current is not None:
                raise ValueError("Nested/unterminated PDBQT MODEL")
            model_number = int(line.split()[1])
            if model_number != len(models) + 1:
                raise ValueError("PDBQT models must be contiguous score ranks")
            current, score = [], None
        elif line.startswith("REMARK VINA RESULT:"):
            if current is None or score is not None:
                raise ValueError("Missing model or duplicate Vina score")
            score = float(line.split()[3])
            if not np.isfinite(score):
                raise ValueError("Nonfinite Vina score")
        elif line.startswith(("ATOM  ", "HETATM")):
            if current is None:
                raise ValueError("Pose atoms outside a MODEL")
            atom_type = line[77:].strip()
            if atom_type in {"H", "HD", "HS"}:
                continue
            current.append((line[12:16].strip(), [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        elif line.startswith("ENDMDL"):
            if current is None or score is None or not current:
                raise ValueError("Empty model or missing Vina score")
            names, coordinates = zip(*current)
            pose = mapped_coordinates(names, coordinates, required_names)
            models.append({"rank": model_number, "vina_score": score, "coordinates": pose})
            current = None
    if current is not None:
        raise ValueError("Truncated PDBQT model")
    if not models:
        raise ValueError("No returned docking poses")
    if any(a["vina_score"] > b["vina_score"] + 0.002 for a, b in zip(models, models[1:])):
        raise ValueError("Vina models are not ordered by score")
    return models


def receptor_atoms(path):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("Expected exactly one receptor model")
    keys, coordinates = [], []
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.element.is_hydrogen:
                    continue
                if atom.altloc not in ("\x00", " "):
                    raise ValueError("Resolve alternate conformations before scoring")
                keys.append((chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, atom.name))
                coordinates.append(list(atom.pos))
    if not keys or len(keys) != len(set(keys)):
        raise ValueError("Empty receptor or duplicate atom identities")
    return keys, np.asarray(coordinates)


def dock(prepared_dir, output_dir, executable, config_path, *, resume=False, seeds=None):
    prepared_dir, output_dir = Path(prepared_dir).resolve(), Path(output_dir).resolve()
    executable, config_path = Path(executable).resolve(), Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    scheduled_seeds = list(config["seeds"] if seeds is None else seeds)
    if not scheduled_seeds or len(set(scheduled_seeds)) != len(scheduled_seeds) or any(seed not in config["seeds"] for seed in scheduled_seeds):
        raise ValueError("Scheduled seeds must be a nonempty unique subset of the frozen configuration")
    preparation = json.loads((prepared_dir / "preparation.json").read_text())
    version = subprocess.check_output([str(executable), "--version"], text=True).strip()
    if version != "AutoDock Vina v1.2.7":
        raise ValueError(f"Unexpected docking version: {version}")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "runs.json"
    if manifest_path.exists() and not resume:
        raise FileExistsError("Run manifest exists; use --resume to verify and reuse matching completed runs")
    old = json.loads(manifest_path.read_text()) if manifest_path.exists() else []
    frozen_config = output_dir / "configuration.json"
    if frozen_config.exists() and sha256(frozen_config) != sha256(config_path):
        raise ValueError("Output configuration changed; use a new output directory")
    if not frozen_config.exists():
        shutil.copyfile(config_path, frozen_config)
    previous = {(r["receptor"], r["seed"]): r for r in old}
    if len(previous) != len(old):
        raise ValueError("Duplicate run identities in existing manifest")
    expected = {(r["id"], seed) for r in preparation["receptors"] for seed in scheduled_seeds}
    if set(previous) - expected:
        raise ValueError("Existing manifest contains runs outside this frozen configuration")
    ligand = prepared_dir / preparation.get("ligand_file", "galactose.pdbqt")
    reference = json.loads((prepared_dir / "reference.json").read_text())
    assembly_file = reference.get("full_assembly_file")
    if assembly_file and not (prepared_dir / assembly_file).is_file():
        raise ValueError("Declared full assembly evaluation coordinates are missing")
    common = {"executable_sha256": sha256(executable), "config_sha256": sha256(config_path),
              "ligand_sha256": sha256(ligand), "preparation_sha256": sha256(prepared_dir / "preparation.json"),
              "reference_sha256": sha256(prepared_dir / "reference.json"),
              "full_assembly_sha256": sha256(prepared_dir / assembly_file) if assembly_file else None}
    runs = []
    for receptor in preparation["receptors"]:
        receptor_path = prepared_dir / receptor["pdbqt_file"]
        fingerprint = {**common, "receptor_sha256": sha256(receptor_path),
                       "receptor_pdb_sha256": sha256(prepared_dir / receptor["pdb_file"])}
        for seed in scheduled_seeds:
            identity = (receptor["id"], seed)
            name = f"{receptor['id']}_{seed}"
            pose_path, log_path = output_dir / (name + ".pdbqt"), output_dir / (name + ".log")
            prior = previous.get(identity)
            if prior and prior["fingerprint"] != fingerprint:
                raise ValueError(f"Inputs/configuration changed for existing run {identity}; use a new output directory")
            if prior and prior["status"] == "ok":
                if not pose_path.exists() or sha256(pose_path) != prior["pose_sha256"]:
                    raise ValueError(f"Completed pose artifact changed: {pose_path}")
                runs.append(prior)
                continue
            command = [str(executable), "--receptor", str(receptor_path), "--ligand", str(ligand), "--out", str(pose_path),
                       "--scoring", "vina", "--cpu", str(config["cpu"]), "--seed", str(seed),
                       "--exhaustiveness", str(config["exhaustiveness"]), "--num_modes", str(config["num_modes"]),
                       "--min_rmsd", str(config["min_rmsd_A"]), "--energy_range", str(config["energy_range_kcal_mol"])]
            for axis, center, size in zip("xyz", preparation["box_center"], preparation["box_size"]):
                command.extend(["--center_" + axis, str(center), "--size_" + axis, str(size)])
            row = {"receptor": receptor["id"], "receptor_state": receptor["state"],
                   "assembly_condition": receptor["assembly_condition"], "seed": seed,
                   "case_id": preparation.get("case_id", "AAV9"), "family_id": preparation.get("family_id", "AAV"),
                   "information_regime": "true-site-local", "status": "running", "fingerprint": fingerprint,
                   "command": command, "pose_file": pose_path.name, "log_file": log_path.name,
                   "platform": platform.platform(), "machine": platform.machine(), "available_logical_cpus": os.cpu_count(),
                   "software": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "gemmi", "rdkit", "meeko")},
                   "python": platform.python_version(), "vina_version": version}
            if prior:
                row["previous_attempts"] = prior.get("previous_attempts", []) + [{k: v for k, v in prior.items() if k != "previous_attempts"}]
            if seeds is not None:
                row["execution_seed_subset"] = scheduled_seeds
            runs.append(row)
            write_json(manifest_path, runs + [r for k, r in previous.items() if k not in {(x["receptor"], x["seed"]) for x in runs}])
            start, usage = time.perf_counter(), resource.getrusage(resource.RUSAGE_CHILDREN)
            try:
                # Keep failed output evidence, but never let it masquerade as a
                # new successful invocation that produced no pose file.
                for artifact in (pose_path, log_path):
                    if artifact.exists():
                        attempt = 1
                        while artifact.with_name(artifact.name + f".previous-{attempt}").exists():
                            attempt += 1
                        artifact.rename(artifact.with_name(artifact.name + f".previous-{attempt}"))
                with log_path.open("w") as log:
                    completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
                row["exit_code"] = completed.returncode
                if completed.returncode:
                    raise RuntimeError(f"Vina returned {completed.returncode}; see {log_path}")
                names = json.loads((prepared_dir / "reference.json").read_text())["names"]
                row["returned_poses"] = len(read_poses(pose_path, names))
                row["pose_sha256"] = sha256(pose_path)
                row["status"] = "ok"
            except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
                row["status"], row["error"] = "failed", str(error)
            finally:
                end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
                row["wall_seconds"] = time.perf_counter() - start
                row["cpu_seconds"] = (end_usage.ru_utime + end_usage.ru_stime) - (usage.ru_utime + usage.ru_stime)
                row["child_peak_rss_bytes"] = int(end_usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024))
                row["memory_measurement"] = "process children high-water RSS; includes earlier runs in this process"
                write_json(manifest_path, runs + [r for k, r in previous.items() if k not in {(x["receptor"], x["seed"]) for x in runs}])
            print(f"{row['case_id']} {row['receptor']} seed {seed}: {row['status']} ({row['wall_seconds']:.1f}s)", flush=True)
    write_json(manifest_path, runs)
    return runs


def analyze(prepared_dir, output_dir):
    prepared_dir, output_dir = Path(prepared_dir), Path(output_dir)
    preparation = json.loads((prepared_dir / "preparation.json").read_text())
    reference = json.loads((prepared_dir / "reference.json").read_text())
    if reference.get("pose_reference_eligible") is False:
        raise ValueError("This assay ligand has no experimental atomic pose reference; use score-panel analysis")
    runs = json.loads((output_dir / "runs.json").read_text())
    if not runs or any(r["status"] == "running" for r in runs):
        raise ValueError("Cannot analyze an empty or unfinished run manifest")
    config_path = output_dir / "configuration.json"
    config = json.loads(config_path.read_text())
    expected = {(r["id"], s) for r in preparation["receptors"] for s in config["seeds"]}
    actual = [(r["receptor"], r["seed"]) for r in runs]
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise ValueError("Incomplete or duplicate receptor-by-seed run denominator")
    native = np.asarray(reference["coordinates"])
    native_keys = [(k[0], k[1], k[2].strip(), k[3], k[4]) for k in reference["receptor_atom_keys"]]
    if len(native_keys) != len(set(native_keys)):
        raise ValueError("Duplicate native receptor atom identity")
    native_receptor = np.asarray(reference["receptor_coordinates"])
    native_contacts = contact_pairs(native, native_receptor)
    receptors = {r["id"]: r for r in preparation["receptors"]}
    full_path = prepared_dir / reference["full_assembly_file"] if reference.get("full_assembly_file") else None
    full = np.load(full_path, allow_pickle=False) if full_path else None
    rows, summaries = [], []
    for run in runs:
        if run["fingerprint"]["reference_sha256"] != sha256(prepared_dir / "reference.json") or run["fingerprint"]["preparation_sha256"] != sha256(prepared_dir / "preparation.json"):
            raise ValueError("Reference/preparation changed since docking")
        if run["fingerprint"]["config_sha256"] != sha256(config_path):
            raise ValueError("Frozen configuration changed since docking")
        receptor = receptors[run["receptor"]]
        if run["fingerprint"].get("receptor_pdb_sha256") != sha256(prepared_dir / receptor["pdb_file"]):
            raise ValueError("Scoring receptor PDB changed since docking")
        if run["fingerprint"].get("full_assembly_sha256") != (sha256(full_path) if full_path else None):
            raise ValueError("Full assembly evaluation coordinates changed since docking")
        if run["status"] != "ok":
            summaries.append({**run, **summarize_poses([])})
            continue
        path = output_dir / run["pose_file"]
        if sha256(path) != run["pose_sha256"]:
            raise ValueError("Pose artifact checksum mismatch")
        keys, coords = receptor_atoms(prepared_dir / receptor["pdb_file"])
        # Extend the fixed holo universe for apo-only atoms: they can add false
        # positive contacts, but cannot remove any native contact denominator.
        universe = dict((k, i) for i, k in enumerate(native_keys))
        for key in keys:
            universe.setdefault(key, len(universe))
        indices = [universe[k] for k in keys]
        residue_ids = [None] * len(universe)
        for k, index in universe.items():
            residue_ids[index] = k[:-1]
        full_coords = full[f"{run['receptor_state']}_coordinates"] if full is not None else coords
        tree = cKDTree(full_coords)
        run_rows = []
        for model in read_poses(path, reference["names"]):
            pose = model.pop("coordinates")
            distance, _ = tree.query(pose)
            severe_pairs = sum(len(v) for v in tree.query_ball_point(pose, 1.5))
            row = {"case_id": run["case_id"], "receptor": run["receptor"], "seed": run["seed"], **model,
                   "rmsd_A": ligand_rmsd(pose, native),
                   **contact_metrics(pose, coords, native_contacts, receptor_residues=residue_ids,
                                     receptor_native_indices=indices, native_receptor_count=len(universe)),
                   "evaluation_receptor_scope": "full_biological_assembly" if full is not None else "provided_receptor_only",
                   "evaluation_min_distance_A": float(distance.min()), "evaluation_pairs_below_1_5A": severe_pairs}
            rows.append(row)
            run_rows.append(row)
        summary = {**run, **summarize_poses(run_rows)}
        summary["top1_score"] = run_rows[0]["vina_score"]
        summary["top1_contact_recall"] = run_rows[0]["contact_recall"]
        summary["top1_contact_precision"] = run_rows[0]["contact_precision"]
        summaries.append(summary)
    write_json(output_dir / "summary.json", summaries)
    if rows:
        with (output_dir / "pose_metrics.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    if full is not None:
        full.close()
    write_json(output_dir / "analysis.json", {
        "schema_version": 1,
        "engine_sha256": sha256(Path(__file__)),
        "metrics_sha256": sha256(Path(__file__).with_name("metrics.py")),
        "run_manifest_sha256": sha256(output_dir / "runs.json"),
        "configuration_sha256": sha256(config_path),
        "reference_sha256": sha256(prepared_dir / "reference.json"),
        "summary_sha256": sha256(output_dir / "summary.json"),
        "pose_metrics_sha256": sha256(output_dir / "pose_metrics.csv") if rows else None,
        "coordinate_frame": "prepared receptor frame; no ligand-only superposition",
        "contact_identity_policy": "full holo atom universe; insertion-code whitespace normalized",
        "contact_cutoff_A": 4.0,
        "overlap_screen_distance_A": 1.5,
    })
    return summaries
