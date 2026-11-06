"""Separate top-one pose preparation for paired first-chain assembly controls."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil

import numpy as np

from .aav import prepare_receptor
from .blind_docking import translate_ligand
from .heldout_assembly_control import CASE_IDS
from .heldout_docking import LEDGER_SHA256, neighborhood_protein
from .heldout_pose_inputs import _portable_archive, select_execution_candidates, validate_candidate_reuse
from .provenance import sha256, write_json, write_text_atomic

METHODS = ("fpocket", "p2rank", "oracle")


def _load(path): return json.loads(Path(path).read_text())


def materialize_strict_profile(prior, output, archive_output):
    """Copy the already fixed strict monomers, binding their actual prediction input."""
    prior, output = Path(prior), Path(output)
    if output.exists(): raise FileExistsError("Preserve previous strict monomer profiles")
    output.mkdir(parents=True); rows = []
    for case_id in CASE_IDS:
        source, destination = prior / case_id, output / case_id
        shutil.copytree(source, destination)
        old = {str(path.relative_to(source)): sha256(path) for path in source.rglob("*") if path.is_file()}
        if any(sha256(destination / name) != digest for name, digest in old.items()): raise AssertionError("Strict monomer copy changed input bytes")
        row = {"case_id": case_id, "profile": "strict_monomer", "selected_author_chain": "A",
               "receptor_sha256": sha256(destination / "receptor.cif"), "receptor_PDB_sha256": sha256(destination / "receptor.pdb"),
               "original_prediction_input_sha256": sha256(destination / "receptor.cif"),
               "ligand_sha256": sha256(destination / "ligand.pdbqt"), "ligand_sdf_sha256": sha256(destination / "ligand.sdf"),
               "canonical_atom_map_sha256": sha256(destination / "ligand_atom_map.json"), "chemistry_sha256": sha256(destination / "chemistry.json"),
               "initializer_report_sha256": sha256(destination / "initializer_report.json"),
               "observed_protein_atom_map_sha256": sha256(destination / "protein_atom_map.json"),
               "full_experimental_assembly_sha256": sha256(destination / "full_assembly_receptor.cif"),
               "full_observed_protein_map_sha256": sha256(destination / "full_protein_atom_map.json"),
               "source_inputs_byte_identical": True, "receptor_support": "prepared" if case_id == "3ZWE" else "preparation_failed",
               "receptor_support_basis": "same observed source residues as fixed strict full profile; strict5HZA incomplete sidechains are retained",
               "native_ligand_coordinate_access": False, "modeled_atoms_are_atomic_truth": False}
        rows.append(row)
    manifest = {"schema_version": 1, "profile": "strict_monomer", "case_ids": list(CASE_IDS), "cases": rows,
                "chain_selection": "first lexical prepared author chain A fixed before monomer predictions",
                "source_inputs_byte_identical": True, "native_ligand_coordinate_access": False,
                "code_sha256": sha256(Path(__file__))}
    write_json(output / "profile_manifest.json", manifest)
    return _portable_archive(output, archive_output, "strict_monomer_sources")


def validate_control_sources(source, profile):
    """Check fixed first-A protein inputs and the shared independent full ligand."""
    source = Path(source); manifest = _load(source / "profile_manifest.json")
    if profile not in ("strict_monomer", "repair_monomer") or manifest["profile"] != profile or manifest["case_ids"] != list(CASE_IDS):
        raise ValueError("Assembly profile must retain both fixed first-chain controls")
    for row in manifest["cases"]:
        if row["case_id"] not in CASE_IDS or row["selected_author_chain"] != "A":
            raise ValueError("Assembly control changed case or first lexical author chain")
        directory = source / row["case_id"]
        for filename, field in (("receptor.cif", "receptor_sha256"), ("receptor.pdb", "receptor_PDB_sha256"),
                                ("ligand.pdbqt", "ligand_sha256"), ("ligand.sdf", "ligand_sdf_sha256"),
                                ("ligand_atom_map.json", "canonical_atom_map_sha256"), ("chemistry.json", "chemistry_sha256"),
                                ("initializer_report.json", "initializer_report_sha256"),
                                ("full_assembly_receptor.cif", "full_experimental_assembly_sha256"),
                                ("full_protein_atom_map.json", "full_observed_protein_map_sha256")):
            if sha256(directory / filename) != row[field]: raise ValueError(f"Assembly source changed: {row['case_id']}/{filename}")
    return manifest


def prepare_control_condition(root, source, profile, case_id, method, predictions, output, config_path, freeze_path, *, oracle_center=None):
    """Reuse strict source-A site predictions for strict and modeled monomer poses."""
    root, source, output = Path(root), Path(source), Path(output)
    config, freeze = _load(config_path), _load(freeze_path)
    if (freeze.get("configuration_sha256") != sha256(config_path) or freeze.get("frozen_before_new_docking") is not True
            or freeze.get("input_adapter_sha256") != sha256(Path(__file__))):
        raise ValueError("Assembly pose preparation requires a separate pre-pose freeze")
    if config["case_ids"] != list(CASE_IDS) or method not in METHODS or case_id not in CASE_IDS:
        raise ValueError("Unknown fixed assembly control/method")
    select_execution_candidates([], config.get("primary_candidate_budget"))
    descriptor = config["source_profiles"][profile]
    if descriptor["manifest"]["sha256"] != sha256(source / "profile_manifest.json"): raise ValueError("Assembly profile manifest changed after freeze")
    manifest = validate_control_sources(source, profile); row = next(r for r in manifest["cases"] if r["case_id"] == case_id)
    if output.exists(): raise FileExistsError("Preserve previous assembly candidate preparations")
    output.mkdir(parents=True)
    condition = {"schema_version": 1, "case_id": case_id, "complex_id": case_id, "condition_id": case_id,
                 "profile": profile, "method": method, "configuration_sha256": sha256(config_path),
                 "fixed_denominator_member": True, "fixed_control_denominator_member": True, "status": "preparing", "jobs": [],
                 "executed_primary_candidate_budget": 1, "all_prediction_candidates_retained_in_selection": True,
                 "native_coordinate_access": method == "oracle", "source_ledger_sha256": LEDGER_SHA256,
                 "evaluation_reference": "full experimental biological assembly, including omitted chains and all qualified source copies"}
    if method == "oracle":
        center = np.asarray(oracle_center, float)
        if center.shape != (3,) or not np.isfinite(center).all() or predictions is not None: raise ValueError("Oracle requires only the common full-assembly native-informed center")
        candidates = [{"rank": 1, "original_rank": 1, "center_A": center.tolist()}]
    else:
        if oracle_center is not None: raise ValueError("Native center forbidden for predicted monomer boxes")
        predictions = Path(predictions); selection_path = predictions / "candidate_selection.json"; selection = _load(selection_path)
        prediction_config = config["prediction_configuration"]
        if sha256(root / prediction_config["path"]) != prediction_config["sha256"]: raise ValueError("Fixed source-monomer prediction config changed")
        candidates = validate_candidate_reuse(selection, case_id, method, prediction_config["sha256"], row["original_prediction_input_sha256"])
        record = selection["prediction_record"]
        if Path(record["filename"]).name != record["filename"] or sha256(predictions / record["filename"]) != record["sha256"]: raise ValueError("Original source-monomer predictor record changed")
        condition.update(candidate_selection=selection, primary_candidate_selection={"path": str(selection_path), "sha256": sha256(selection_path)}, candidate_centers_and_ranks_exactly_reused=True)
    condition["available_prediction_candidates"] = len(candidates); candidates = select_execution_candidates(candidates, 1)
    condition["scheduled_candidate_jobs"] = len(candidates)
    directory = source / case_id; chemistry, atom_map = _load(directory / "chemistry.json"), _load(directory / "ligand_atom_map.json")["atoms"]
    for candidate in candidates:
        box_id = f"box{candidate['rank']:02d}"; box = output / box_id; box.mkdir()
        job = {"id": box_id, "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"], "center_A": candidate["center_A"], "prepared_directory": box_id}
        try:
            protein, neighborhood = neighborhood_protein(directory / "receptor.cif", candidate["center_A"], chemistry["box_size"])
            receptor_audit = prepare_receptor(protein, box / "neighborhood")
            ligand, translation = translate_ligand((directory / "ligand.pdbqt").read_text(), candidate["center_A"])
            write_text_atomic(box / "ligand.pdbqt", ligand); shutil.copy2(directory / "ligand_atom_map.json", box / "ligand_atom_map.json")
            write_json(box / "reference.json", {"names": [a["name"] for a in atom_map], "canonical_atom_map": atom_map,
                                                 "pose_reference_eligible": False, "purpose": "canonical parent identities only; all atomic truth remains in full experimental assembly"})
            prep = {"case_id": case_id, "family_id": case_id, "method": method, "profile": profile, "status": "prepared",
                    "configuration_sha256": sha256(config_path), "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"],
                    "information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                    "execution_information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                    "box_center": candidate["center_A"], "box_size": chemistry["box_size"], "ligand_file": "ligand.pdbqt",
                    "receptors": [{"id": "neighborhood", "state": "holo", "assembly_condition": "first_chain_A_8A_box_neighborhood", "pdb_file": "neighborhood.pdb", "pdbqt_file": "neighborhood.pdbqt"}],
                    "neighborhood": neighborhood, "ligand_translation": translation, "receptor_preparation": receptor_audit,
                    "source_receptor_sha256": row["receptor_sha256"], "independent_ligand_sha256": row["ligand_sha256"],
                    "native_coordinates_in_runner_reference": False, "modeled_atoms_are_atomic_truth": False,
                    "full_assembly_evaluation_archive_manifest": config["evaluation_source_archive_manifest"]}
            write_json(box / "preparation.json", prep); job.update(status="prepared", preparation_sha256=sha256(box / "preparation.json"))
        except (ValueError, RuntimeError, OSError) as error:
            job.update(status="preparation_failed", error_type=type(error).__name__, reason=str(error)); write_json(box / "preparation_attempt.json", job)
        condition["jobs"].append(job); write_json(output / "condition.json", condition)
    condition["status"] = ("prepared" if candidates and all(j["status"] == "prepared" for j in condition["jobs"])
                           else selection["status"] if method != "oracle" and not candidates else "preparation_failed")
    write_json(output / "condition.json", condition); return condition
