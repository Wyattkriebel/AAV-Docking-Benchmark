"""Native-coordinate scoring controls, kept separate from docking success.

Experimental heavy coordinates are held fixed. RDKit adds modeled hydrogens
without minimization, and the frozen Vina binary evaluates the resulting pose.
This cannot diagnose internal search sampling or predict a dissociation constant.
"""
from __future__ import annotations

from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np

from .aav import validate_native_stereochemistry
from .provenance import sha256, write_json


def _ligand_atoms(text):
    atoms = {}
    torsdof = []
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            if name in atoms:
                raise ValueError(f"duplicate ligand atom name {name}")
            kind = line[77:].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if not kind or not np.isfinite(xyz).all():
                raise ValueError("invalid ligand coordinates/type")
            atoms[name] = {"coordinates": xyz, "type": kind,
                           "hydrogen": kind in {"H", "HD", "HS"}}
        elif line.startswith("TORSDOF"):
            torsdof.append(int(line.split()[1]))
    if not atoms or len(torsdof) != 1 or torsdof[0] < 0:
        raise ValueError("expected one complete ligand PDBQT with TORSDOF")
    return atoms, torsdof[0]


def build_native_ligand(prepared_dir, output_dir):
    """Create native-heavy-coordinate ligand and a complete H-placement ledger."""
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    prepared_dir, output_dir = Path(prepared_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    preparation_path, reference_path = prepared_dir / "preparation.json", prepared_dir / "reference.json"
    preparation = json.loads(preparation_path.read_text())
    reference = json.loads(reference_path.read_text())
    source_pdbqt = prepared_dir / preparation.get("ligand_file", "galactose.pdbqt")
    source_sdf = prepared_dir / "galactose.sdf"
    sources = {"preparation": preparation_path, "reference": reference_path,
               "source_ligand_sdf": source_sdf, "source_ligand_pdbqt": source_pdbqt}
    fingerprints = {name: sha256(path) for name, path in sources.items()}
    molecules = list(Chem.SDMolSupplier(str(source_sdf), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("expected one valid frozen source ligand SDF")
    heavy = Chem.RemoveHs(molecules[0])
    graph_names = preparation["ligand"]["ccd_atom_names"]
    if len(graph_names) != heavy.GetNumAtoms() or len(set(graph_names)) != len(graph_names):
        raise ValueError("source graph/name ledger is incomplete or ambiguous")
    source_atoms, source_torsdof = _ligand_atoms(source_pdbqt.read_text())
    source_heavy = {name: row for name, row in source_atoms.items() if not row["hydrogen"]}
    if set(source_heavy) != set(graph_names):
        raise ValueError("source PDBQT and ligand graph atom identities differ")
    source_xyz = heavy.GetConformer().GetPositions()
    source_mapping_error = max(float(np.linalg.norm(source_xyz[i] - source_heavy[name]["coordinates"]))
                               for i, name in enumerate(graph_names))
    # SDF and PDBQT round independent coordinates to 4 and 3 decimals.
    if source_mapping_error > 0.0011:
        raise ValueError("source SDF/PDBQT atom map is inconsistent")
    names = reference["names"]
    xyz = np.asarray(reference["coordinates"], float)
    if len(set(names)) != len(names) or set(names) != set(graph_names) or xyz.shape != (len(names), 3):
        raise ValueError("native atom map is incomplete or ambiguous")
    if not np.isfinite(xyz).all():
        raise ValueError("native ligand coordinates must be finite")
    native_by_name = dict(zip(names, xyz))
    for atom, name in zip(heavy.GetAtoms(), graph_names):
        atom.SetProp("ccd_name", name)
        heavy.GetConformer().SetAtomPosition(atom.GetIdx(), native_by_name[name])
    chemistry = validate_native_stereochemistry(heavy, names, xyz)
    molecule = Chem.AddHs(heavy, addCoords=True)
    placements = []
    hydrogen_number = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            name = atom.GetProp("ccd_name")
        else:
            hydrogen_number += 1
            name = f"H{hydrogen_number}"
            neighbors = list(atom.GetNeighbors())
            if len(neighbors) != 1 or neighbors[0].GetAtomicNum() <= 1:
                raise ValueError("hydrogen does not have exactly one heavy parent")
            parent = neighbors[0].GetProp("ccd_name")
            placements.append({"name": name, "parent_ccd_atom": parent,
                               "coordinates_A": list(molecule.GetConformer().GetAtomPosition(atom.GetIdx())),
                               "method": "RDKit Chem.AddHs(addCoords=True); no optimization"})
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4))
        info.SetResidueName("GAL")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        atom.SetMonomerInfo(info)
    setup = MoleculePreparation().prepare(molecule)[0]
    pdbqt, success, error = PDBQTWriterLegacy.write_string(setup)
    if not success:
        raise ValueError(error)
    control_atoms, control_torsdof = _ligand_atoms(pdbqt)
    control_heavy = {name: row for name, row in control_atoms.items() if not row["hydrogen"]}
    if set(control_heavy) != set(native_by_name) or control_torsdof != source_torsdof:
        raise ValueError("native ligand atom map or torsional freedom changed during preparation")
    displacement = np.array([control_heavy[name]["coordinates"] - native_by_name[name] for name in names])
    max_error = float(np.max(np.linalg.norm(displacement, axis=1)))
    if max_error > 0.0009:
        raise ValueError("native heavy coordinates moved beyond PDBQT rounding tolerance")
    for placement in placements:
        retained = control_atoms.get(placement["name"])
        placement["retained_in_pdbqt"] = retained is not None
        placement["pdbqt_coordinates_A"] = retained["coordinates"].tolist() if retained else None
        placement["pdbqt_type"] = retained["type"] if retained else "merged nonpolar H"
    pdbqt_path, sdf_path = output_dir / "native_ligand.pdbqt", output_dir / "native_ligand.sdf"
    pdbqt_path.write_text(pdbqt)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(molecule)
    writer.close()
    result = {"status": "prepared", "diagnostic_only": True, "count_as_docking_success": False,
              "native_reference": reference.get("reference_ligand"),
              "heavy_atom_names_in_graph_order": graph_names,
              "native_heavy_atom_count": len(names), "modeled_hydrogen_count": len(placements),
              "max_native_heavy_coordinate_error_A": max_error,
              "native_heavy_rmsd_A": float(np.sqrt(np.mean(np.sum(displacement**2, axis=1)))),
              "source_graph_map_max_rounding_error_A": source_mapping_error,
              "source_torsdof": source_torsdof, "native_torsdof": control_torsdof,
              "native_heavy_atoms_optimized": False, "hydrogen_optimization": False,
              "hydrogen_placement": placements, "chemistry": chemistry,
              "fingerprint": fingerprints, "native_ligand_sha256": sha256(pdbqt_path),
              "versions": {name: importlib.metadata.version(name) for name in ("rdkit", "meeko", "numpy")},
              "interpretation": "Native-pose score diagnostic with modeled H; reference-informed and never a docking success"}
    for name, path in sources.items():
        if sha256(path) != fingerprints[name]:
            raise ValueError("frozen source inputs changed while constructing native control")
    write_json(output_dir / "ligand_preparation.json", result)
    return result


def parse_score_only(text):
    values = re.findall(r"Estimated Free Energy of Binding\s*:\s*([-+]?\d+(?:\.\d+)?)\s*\(kcal/mol\)", text)
    if len(values) != 1:
        raise ValueError("Vina score-only output has no unique binding-energy result")
    value = float(values[0])
    if not np.isfinite(value):
        raise ValueError("Vina native-pose score must be finite")
    return value


def run_native_control(prepared_dir, output_dir, executable, docking_results_dir, receptor_ids=None):
    """Evaluate frozen native coordinates with the same Vina score and box."""
    prepared_dir, output_dir = Path(prepared_dir).resolve(), Path(output_dir).resolve()
    executable, docking_results_dir = Path(executable).resolve(), Path(docking_results_dir).resolve()
    preparation = json.loads((prepared_dir / "preparation.json").read_text())
    docking = json.loads((docking_results_dir / "summary.json").read_text())
    config_path = docking_results_dir / "configuration.json"
    config = json.loads(config_path.read_text())
    version = subprocess.check_output([str(executable), "--version"], text=True).strip()
    if version != "AutoDock Vina v1.2.7":
        raise ValueError("native control requires the frozen official Vina 1.2.7 executable")
    receptor_ids = set(receptor_ids) if receptor_ids is not None else {r["id"] for r in preparation["receptors"]}
    selected = [r for r in preparation["receptors"] if r["id"] in receptor_ids]
    if {r["id"] for r in selected} != receptor_ids:
        raise ValueError("unknown native-control receptor ID")
    ligand = build_native_ligand(prepared_dir, output_dir)
    executable_hash = sha256(executable)
    records = []
    for receptor in selected:
        receptor_path = prepared_dir / receptor["pdbqt_file"]
        previous = [row for row in docking if row["receptor"] == receptor["id"] and row["status"] == "ok"]
        if not previous:
            raise ValueError("native-score comparator needs completed docking results")
        expected_fingerprint = {"executable_sha256": executable_hash, "receptor_sha256": sha256(receptor_path),
                                "config_sha256": sha256(config_path),
                                "preparation_sha256": sha256(prepared_dir / "preparation.json"),
                                "reference_sha256": sha256(prepared_dir / "reference.json")}
        if any(any(row["fingerprint"].get(name) != value for name, value in expected_fingerprint.items()) for row in previous):
            raise ValueError("native control inputs do not match frozen docking comparator")
        command = [str(executable), "--receptor", str(receptor_path), "--ligand", str(output_dir / "native_ligand.pdbqt"),
                   "--scoring", "vina", "--cpu", str(config["cpu"]), "--score_only"]
        for axis, center, size in zip("xyz", preparation["box_center"], preparation["box_size"]):
            command.extend(["--center_" + axis, str(center), "--size_" + axis, str(size)])
        log_path = output_dir / f"{receptor['id']}_score_only.log"
        start = time.perf_counter()
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
        log_path.write_text(completed.stdout)
        if completed.returncode:
            raise RuntimeError(f"native score-only failed; inspect {log_path}")
        energy = parse_score_only(completed.stdout)
        comparison = [{"seed": row["seed"], "docked_top1_score_kcal_mol": row["top1_score"],
                       "docked_top1_rmsd_A": row["top1_rmsd_A"],
                       "native_minus_docked_top1_score_kcal_mol": energy - row["top1_score"],
                       "pose_file": row["pose_file"], "pose_sha256": row["pose_sha256"]} for row in previous]
        record = {"receptor": receptor["id"], "receptor_state": receptor["state"],
                  "assembly_condition": receptor["assembly_condition"], "status": "ok",
                  "diagnostic_only": True, "count_as_docking_success": False,
                  "native_score_only_kcal_mol": energy, "docked_top1_comparisons": comparison,
                  "hydrogen_model": "RDKit AddHs coordinates; no hydrogen or heavy-atom optimization",
                  "native_heavy_rmsd_A": ligand["native_heavy_rmsd_A"],
                  "command": command, "log_file": log_path.name, "log_sha256": sha256(log_path),
                  "fingerprint": {**expected_fingerprint, "native_ligand_sha256": ligand["native_ligand_sha256"]},
                  "wall_seconds": time.perf_counter() - start, "vina_version": version,
                  "interpretation": "Scoring diagnostic; native reference was supplied, no docking recovery or internal-sampling inference"}
        records.append(record)
        write_json(output_dir / "summary.json", records)
    return records
