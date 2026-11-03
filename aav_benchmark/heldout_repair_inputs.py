"""Audit and archive separately modeled repair inputs without changing truth.

This post-generation packaging step never changes the frozen strict inputs or
the frozen PDBFixer worker.  Observed protein maps and ligand graphs are copied
verbatim; added protein atoms have a distinct, explicitly modeled atom ledger.
"""
from __future__ import annotations

import copy
import gzip
import io
import json
from pathlib import Path
import shutil
import tarfile

import numpy as np

from .aav import _pdb_atoms
from .heldout_repair import CASE_IDS, audit_repair, parse_pdb_heavy
from .heldout_serialization import atom_semantics, compare_atom_semantics, finalize_cif_metadata
from .provenance import sha256, write_json


def _load(path):
    return json.loads(Path(path).read_text())


def _check_parameterization(source_text, pdbqt_path):
    expected = parse_pdb_heavy(source_text)
    actual = {key: row for key, row in _pdb_atoms(Path(pdbqt_path).read_text(), pdbqt=True).items() if not row["hydrogen"]}
    if set(expected) != set(actual):
        raise ValueError("Modeled receptor parameterization changed heavy atom identities")
    maximum = max(float(np.linalg.norm(np.asarray(expected[key]["coordinates"]) - actual[key]["coordinates"])) for key in expected)
    if maximum > .0009:
        raise ValueError("Modeled receptor PDBQT moved heavy coordinates beyond serialization rounding")
    return {"status": "passed", "heavy_atom_count": len(expected), "max_heavy_serialization_displacement_A": maximum,
            "pdbqt_sha256": sha256(pdbqt_path)}


def audit_and_archive_repaired(root, primary, generated, staging, output, *, config_path="benchmarks/structural_repair_experiment.json",
                               freeze_path="results/structural_repair/freeze_receipt.json"):
    """Create a portable all-six source archive, preserving unsupported MSE."""
    root, primary, generated, staging, output = map(Path, (root, primary, generated, staging, output))
    config_path, freeze_path = root / config_path, root / freeze_path
    configuration, freeze = _load(config_path), _load(freeze_path)
    source_receipt = _load(primary / "source_preparation.json")
    generation_path = generated / "repair_preparation.json"
    generation = _load(generation_path)
    if (configuration["case_ids"] != list(CASE_IDS) or generation["fixed_case_ids"] != list(CASE_IDS)
            or generation["configuration_sha256"] != sha256(config_path)
            or generation["freeze_sha256"] != sha256(freeze_path)
            or freeze["configuration_sha256"] != sha256(config_path)
            or freeze.get("secondary_profile_frozen_before_new_outcomes") is not True
            or generation.get("primary_strict_inputs_changed") is not False
            or generation.get("modeled_atoms_are_atomic_truth") is not False):
        raise ValueError("Modeled inputs are not bound to the separate frozen sensitivity")
    if staging.exists() or staging.is_symlink(): raise FileExistsError("Preserve earlier repaired source collections")
    archive_path, manifest_path, audit_path = (output / name for name in (
        "prepared_sources_v2.tar.gz", "prepared_sources_v2_manifest.json", "prepared_sources_v2_audit.json"))
    if any(path.exists() for path in (archive_path, manifest_path, audit_path)):
        raise FileExistsError("Preserve earlier repaired source archives")
    original_files = {str(p.relative_to(primary)): sha256(p) for p in primary.rglob("*") if p.is_file()}
    generated_files = {str(p.relative_to(generated)): sha256(p) for p in generated.rglob("*") if p.is_file()}
    by_case = {row["case_id"]: row for row in generation["cases"]}
    if set(by_case) != set(CASE_IDS): raise ValueError("Modeled source cohort is incomplete")
    staging.mkdir(parents=True); output.mkdir(parents=True, exist_ok=True)
    rows = []
    for case_id in CASE_IDS:
        strict = primary / case_id; modeled = generated / case_id; destination = staging / case_id
        destination.mkdir()
        source = (strict / "receptor.pdb").read_text()
        expected = configuration["planned_sidechain_atom_keys"][case_id]
        record = by_case[case_id]
        supported = record["status"] == "prepared"
        if supported:
            if record.get("receptor_parameterization", {}).get("status") != "prepared":
                raise ValueError("Supported repair source lacks validated chemical parameterization")
            audit = audit_repair(source, (modeled / "receptor.pdb").read_text(), expected)
            if record["new_atom_count"] != len(expected) or audit["added_heavy_atoms"] != len(expected):
                raise ValueError("Actual modeled atom count differs from the predeclared allowlist")
            if record["determinism"]["max_new_atom_displacement_A"] > configuration["max_new_atom_repeat_displacement_A"]:
                raise ValueError("Fixed-seed repeat exceeded the frozen modeling tolerance")
            parameterization = _check_parameterization((modeled / "receptor.pdb").read_text(), modeled / "parameterized.pdbqt")
            for path in modeled.iterdir():
                if path.is_file(): shutil.copy2(path, destination / path.name)
            copied_cif = destination / "receptor.cif"; copied_cif.unlink()
            serialization = finalize_cif_metadata(modeled / "receptor.cif", copied_cif)
            write_json(destination / "serialization_v2.json", serialization)
        else:
            if record["status"] != "unsupported" or case_id != "4R9F" or expected:
                raise ValueError("Unexpected repair feasibility/unsupported identity")
            # Preserve the unconverted MSE receptor as an explicit failed workflow
            # member; no modeled coordinate or unsupported charge is invented.
            for name in ("receptor.pdb", "receptor.cif"):
                shutil.copy2(strict / name, destination / name)
            for path in modeled.iterdir():
                if path.is_file(): shutil.copy2(path, destination / path.name)
            audit = audit_repair(source, (destination / "receptor.pdb").read_text(), [])
            parameterization = {"status": "unsupported", "reason": record["reason"]}
            serialization = compare_atom_semantics(atom_semantics(strict / "receptor.cif"), atom_semantics(destination / "receptor.cif"))
        for name in ("ligand_atom_map.json", "ligand.sdf", "ligand.pdbqt", "chemistry.json", "protein_atom_map.json", "source_symmetry.json"):
            if supported and name in ("ligand_atom_map.json", "ligand.sdf", "ligand.pdbqt", "chemistry.json"):
                if sha256(modeled / name) != sha256(strict / name):
                    raise ValueError("Repair generation changed the independent ligand graph/conformer")
            shutil.copy2(strict / name, destination / name)
        # This observed map remains the only protein target for receptor fits.
        observed_map = _load(destination / "protein_atom_map.json")
        observed_keys = {tuple(row["output_key"]) for row in observed_map["atoms"]}
        if observed_keys != set(parse_pdb_heavy(source)):
            raise ValueError("Observed source protein map is incomplete")
        modeled_map = {"schema_version": 1, "atoms": [{"output_key": row["key"], "element": row["element"],
                           "modeled_xyz_A": row["coordinates_A"], "source_atom_site_id": None,
                           "observed_atom": False, "atomic_truth_eligible": False,
                           "method": "frozen PDBFixer1.12/OpenMM8.4.0 Reference seed20261024"} for row in audit["added_atoms"]],
                       "observed_protein_atom_map_sha256": sha256(destination / "protein_atom_map.json"),
                       "modeled_atoms_are_atomic_truth": False}
        write_json(destination / "modeled_atom_map.json", modeled_map)
        manifest = copy.deepcopy(_load(strict / "protein_preparation.json"))
        manifest.update(receptor_sha256=sha256(destination / "receptor.cif"), profile="standard_sidechain_completion_sensitivity",
                        protein_heavy_atom_count=audit["observed_heavy_atoms"]+audit["added_heavy_atoms"],
                        observed_protein_heavy_atom_count=audit["observed_heavy_atoms"], modeled_protein_heavy_atom_count=audit["added_heavy_atoms"],
                        observed_atom_map_sha256=sha256(destination / "protein_atom_map.json"), modeled_atom_map_sha256=sha256(destination / "modeled_atom_map.json"),
                        model="uniform dry protein; separately modeled missing standard sidechains; all observed atoms fixed; no loops/backbone/terminal additions or MSE conversion",
                        native_ligand_access=False, distance_based_selection=False, atomic_truth_changed=False)
        write_json(destination / "receptor.cif.manifest.json", manifest)
        write_json(destination / "protein_preparation.json", manifest)
        row = {"case_id": case_id, "status": "prepared" if supported else "unsupported", "fixed_denominator_member": True,
               "repair_audit": audit, "parameterization_audit": parameterization, "serialization_audit": serialization,
               "observed_protein_map_byte_identical": True, "independent_ligand_files_byte_identical": True,
               "worker_receipt_sha256": sha256(destination / "worker_receipt.json"),
               "source_receptor_sha256": sha256(destination / "receptor.cif"), "native_glycan_coordinate_access": False,
               "modeled_atoms_are_atomic_truth": False}
        write_json(destination / "source_preparation.json", row); rows.append(row)
    if any(sha256(primary / name) != digest for name, digest in original_files.items()):
        raise AssertionError("Strict primary input tree changed")
    if any(sha256(generated / name) != digest for name, digest in generated_files.items()):
        raise AssertionError("Original frozen repair generation changed")
    collection = {"schema_version": 1, "profile": "standard_sidechain_completion_sensitivity", "fixed_case_ids": list(CASE_IDS),
                  "configuration_sha256": sha256(config_path), "freeze_sha256": sha256(freeze_path),
                  "original_generation_receipt_sha256": sha256(generation_path),
                  "primary_source_preparation_sha256": sha256(primary / "source_preparation.json"),
                  "observed_maps_and_independent_ligands_unchanged": True, "primary_and_generation_trees_unchanged": True,
                  "modeled_atoms_are_atomic_truth": False, "native_glycan_coordinate_access": False, "cases": rows,
                  "packaging_code_sha256": sha256(Path(__file__))}
    write_json(staging / "source_preparation.json", collection)
    files = []
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in sorted(staging.rglob("*")):
                    if not path.is_file(): continue
                    name = str(path.relative_to(staging)); content = path.read_bytes()
                    info = tarfile.TarInfo(name); info.size = len(content); info.mode = 0o644; info.mtime = 0
                    archive.addfile(info, io.BytesIO(content))
                    files.append({"path": name, "sha256": sha256(path), "bytes": len(content)})
    manifest = {"schema_version": 1, "archive": archive_path.name, "archive_sha256": sha256(archive_path),
                "files": files, "fixed_case_ids": list(CASE_IDS), "profile": collection["profile"],
                "modeled_atoms_are_atomic_truth": False, "native_glycan_coordinate_access": False}
    write_json(manifest_path, manifest); write_json(audit_path, collection)
    return {"archive_manifest": str(manifest_path), "archive_manifest_sha256": sha256(manifest_path),
            "archive": str(archive_path), "archive_sha256": sha256(archive_path), "file_count": len(files), "cases": rows}
