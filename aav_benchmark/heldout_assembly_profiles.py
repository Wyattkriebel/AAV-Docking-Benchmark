"""Paired strict/modeled monomer inputs from the same fixed full assemblies.

Only the predeclared first lexical author chain A is retained.  Modeled monomer
atoms are exact subsets of already generated full-receptor repair coordinates;
no atom is modeled here and full-assembly experimental truth remains separate.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil

from .aav import prepare_receptor
from .heldout_assembly_control import CASE_IDS, subset_monomer_input, _subset_pdb
from .heldout_pose_inputs import _portable_archive
from .heldout_repair import parse_pdb_heavy
from .heldout_serialization import atom_semantics
from .provenance import sha256, write_json


def _load(path): return json.loads(Path(path).read_text())


def build_repair_monomers(strict_monomer, full_repair, output, archive_output):
    """Subset the exact repaired full source and retain original observed maps."""
    strict_monomer, full_repair, output, archive_output = map(Path, (strict_monomer, full_repair, output, archive_output))
    if output.exists(): raise FileExistsError("Preserve previous paired assembly controls")
    full_manifest = _load(full_repair / "profile_manifest.json")
    if full_manifest["profile"] != "repair": raise ValueError("Modeled monomers must be subsets of the same full repair profile")
    full_cases = {row["case_id"]: row for row in full_manifest["cases"]}
    output.mkdir(parents=True); rows = []
    for case_id in CASE_IDS:
        strict, full, destination = strict_monomer / case_id, full_repair / case_id, output / case_id
        source_record = full_cases[case_id]
        if sha256(full / "receptor.cif") != source_record["receptor_sha256"] or source_record["receptor_support"] != "prepared":
            raise ValueError("Full modeled source changed or is unsupported")
        before_hashes = {str(path.relative_to(full)): sha256(path) for path in full.rglob("*") if path.is_file()}
        shutil.copytree(strict, destination)
        for name in ("receptor.cif", "receptor.pdb"): (destination / name).unlink()
        proof = subset_monomer_input(full / "receptor.cif", destination / "receptor.cif")
        if proof["selected_author_chain"] != "A": raise ValueError("The predeclared monomer must remain first lexical author chain A")
        expected = {key: row for key, row in atom_semantics(full / "receptor.cif").items() if key[1] == "A"}
        _subset_pdb(full / "receptor.pdb", destination / "receptor.pdb", "A", expected)
        observed = _load(strict / "protein_atom_map.json")
        modeled = copy.deepcopy(_load(full / "modeled_atom_map.json"))
        modeled["atoms"] = [row for row in modeled["atoms"] if row["output_key"][0] == "A"]
        observed_keys = {tuple(row["output_key"]) for row in observed["atoms"]}
        modeled_keys = {tuple(row["output_key"]) for row in modeled["atoms"]}
        actual_keys = set(parse_pdb_heavy((destination / "receptor.pdb").read_text()))
        if observed_keys & modeled_keys or observed_keys | modeled_keys != actual_keys:
            raise ValueError("Monomer observed and modeled maps do not partition the complete source subset")
        if any(row.get("atomic_truth_eligible") is not False or row.get("source_atom_site_id") is not None for row in modeled["atoms"]):
            raise ValueError("Modeled monomer atoms cannot become experimental truth")
        write_json(destination / "modeled_atom_map.json", modeled)
        # Original full experimental assembly and observed map are inherited
        # unchanged from strict monomer controls. The modeled full source is
        # retained separately solely as preparation provenance.
        shutil.copy2(full / "receptor.cif", destination / "full_modeled_assembly_receptor.cif")
        for name in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json", "chemistry.json", "initializer_report.json"):
            if sha256(full / name) != sha256(strict / name): raise ValueError("Paired assembly profiles changed the independent ligand")
        parameterization = prepare_receptor((destination / "receptor.pdb").read_text(), destination / "parameterized")
        manifest = _load(strict / "protein_preparation.json")
        manifest.update(receptor_sha256=sha256(destination / "receptor.cif"), no_heavy_atom_repair=False,
                        scope="physically incomplete first-chain monomer; exact subset of separately modeled full receptor",
                        protein_heavy_atom_count=len(actual_keys), observed_protein_heavy_atom_count=len(observed_keys),
                        modeled_protein_heavy_atom_count=len(modeled_keys),
                        modeled_atom_map_sha256=sha256(destination / "modeled_atom_map.json"),
                        original_observed_map_byte_identical=True, modeled_atoms_are_atomic_truth=False,
                        no_new_modeling_for_monomer=True)
        write_json(destination / "protein_preparation.json", manifest)
        write_json(destination / "receptor.cif.manifest.json", manifest)
        write_json(destination / "receptor_parameterization.json", {"status": "prepared", "audit": parameterization})
        row = {"case_id": case_id, "profile": "repair_monomer", "selected_author_chain": "A", "receptor_support": "prepared",
               "receptor_sha256": sha256(destination / "receptor.cif"), "receptor_PDB_sha256": sha256(destination / "receptor.pdb"),
               "original_prediction_input_sha256": sha256(strict / "receptor.cif"),
               "observed_heavy_atoms": len(observed_keys), "modeled_heavy_atoms": len(modeled_keys),
               "exact_full_repair_source_subset": proof, "no_new_modeling": True,
               "ligand_sha256": sha256(destination / "ligand.pdbqt"), "ligand_sdf_sha256": sha256(destination / "ligand.sdf"),
               "canonical_atom_map_sha256": sha256(destination / "ligand_atom_map.json"),
               "chemistry_sha256": sha256(destination / "chemistry.json"),
               "initializer_report_sha256": sha256(destination / "initializer_report.json"),
               "full_experimental_assembly_sha256": sha256(destination / "full_assembly_receptor.cif"),
               "full_observed_protein_map_sha256": sha256(destination / "full_protein_atom_map.json"),
               "full_modeled_assembly_source_sha256": sha256(destination / "full_modeled_assembly_receptor.cif"),
               "observed_atoms_and_full_experimental_truth_unchanged": True, "modeled_atoms_are_atomic_truth": False,
               "native_ligand_coordinate_access": False}
        write_json(destination / "assembly_profile_preparation.json", row); rows.append(row)
        if any(sha256(full / name) != digest for name, digest in before_hashes.items()): raise AssertionError("Full modeled source changed")
    manifest = {"schema_version": 1, "profile": "repair_monomer", "case_ids": list(CASE_IDS), "cases": rows,
                "selection_rule": "first lexical author chain A fixed for both cases before any monomer predictions; no ligand/site/interface/outcome selection",
                "source_profile_manifest_sha256": sha256(full_repair / "profile_manifest.json"),
                "same_strict_monomer_predictions_reused_for_modeled_monomer": True,
                "no_new_modeling": True, "native_ligand_coordinate_access": False, "modeled_atoms_are_atomic_truth": False,
                "code_sha256": sha256(Path(__file__))}
    write_json(output / "profile_manifest.json", manifest)
    archive = _portable_archive(output, archive_output, "repair_monomer_sources")
    proposal = {"schema_version": 1, "case_ids": list(CASE_IDS), "receptor_profiles": ["strict_monomer", "repair_monomer"],
                "prediction_source": "original unmodeled strict first-chain monomer inputs only; reuse exactly for modeled monomer docking",
                "full_comparators": "same-case strict full assembly vs strict monomer and modeled full assembly vs modeled monomer",
                "chain_selection": "A, first lexical prepared author chain for both fixed cases; independent of native ligands, interfaces and outcomes",
                "site_methods": ["fpocket", "p2rank"], "site_candidate_budgets": [1, 5], "primary_candidate_budget": 1,
                "pose_methods": ["fpocket", "p2rank", "oracle"], "pose_rank_budgets": [1, 5],
                "site_prediction_parameters": {"p2rank_seed": 42, "threads": 1, "candidate_merge_distance_A": 4., "top_k": 5,
                                               "orbit_deduplication": "no_orbit_dedupe"},
                "ligand_initializer": "exact same independently validated CCD-ideal v3 bytes and full canonical parent mapping as full-assembly comparators",
                "receptor_geometry": "strict monomer exact strict-full atom subset; modeled monomer exact preexisting modeled-full atom subset; no new modeling",
                "repair_monomer_source_archive_manifest": archive,
                "strict_monomer_profile_manifest": {"path": str(strict_monomer / "profile_manifest.json"), "sha256": sha256(strict_monomer / "profile_manifest.json")},
                "evaluation": "unchanged full experimental assemblies, qualified native copies and canonical masks, native contacts and full-assembly clashes; omitted chains/sites remain in truth",
                "interpretation": "two-family paired assembly-ablation sensitivity; physically incomplete monomer controls are not independent biological samples",
                "failure_policy": "retain strict5HZA incomplete-residue preparation failures; no new atoms beyond the separately frozen full repair profile; all methods and uncovered known sites retained",
                "cpu": 2, "exhaustiveness": 16, "num_modes": 20, "energy_range_kcal_mol": 10, "min_rmsd_A": .5,
                "seeds": [104729, 130363, 155921], "source_only_before_monomer_prediction": True, "new_modeling_performed": False,
                "code_sha256": sha256(Path(__file__))}
    write_json(archive_output / "assembly_profiles.proposal.json", proposal)
    return {"manifest": manifest, "archive": archive, "proposal_sha256": sha256(archive_output / "assembly_profiles.proposal.json")}
