"""Source-only monomer controls for the fixed 3ZWE and 5HZA assemblies.

The first lexical prepared author chain is selected before any new prediction.
No ligand coordinates, sites, protein interfaces, or outcomes choose the chain.
These physically incomplete receptors are controls, not new biological samples.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import gemmi

from .heldout_serialization import (
    atom_semantics, compare_atom_semantics, finalize_cif_metadata,
)

CASE_IDS = ("3ZWE", "5HZA")
LIGAND_FILES = ("ligand.pdbqt", "ligand.sdf", "chemistry.json", "ligand_atom_map.json")
REQUIRED_FILES = ("receptor.cif", "receptor.pdb", "protein_atom_map.json", *LIGAND_FILES)
LOCK_FILES = ("benchmarks/site_methods.lock.json", "benchmarks/p2rank.lock.json", "benchmarks/vina.lock.json")
ROOT = Path(__file__).resolve().parents[1]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text())


def _write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _blank(value):
    return "" if value in (None, False, ".", "?", " ") else str(value)


def subset_monomer_input(source_cif, output_cif):
    """Select the first lexical author chain and preserve its exact atom rows."""
    source_cif, output_cif = Path(source_cif), Path(output_cif)
    if output_cif.exists() or output_cif.is_symlink():
        raise FileExistsError("Preserve previous assembly controls; use a fresh output")
    source_hash = _sha(source_cif)
    before = atom_semantics(source_cif)
    chains = sorted({key[1] for key in before})
    if len(chains) < 2 or any(not chain for chain in chains):
        raise ValueError("Assembly control requires multiple unambiguous protein author chains")
    selected_chain = chains[0]
    block = gemmi.cif.read_file(str(source_cif)).sole_block()
    atoms = block.get_mmcif_category("_atom_site.")
    selected = [i for i, chain in enumerate(atoms["auth_asym_id"]) if _blank(chain) == selected_chain]
    subset = {field: [values[i] for i in selected] for field, values in atoms.items()}
    document = gemmi.cif.Document()
    document.add_new_block("source_only_assembly_control").set_mmcif_category("_atom_site.", subset)
    output_cif.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="assembly_subset_", dir=output_cif.parent) as temporary:
        intermediate = Path(temporary) / "source_subset.cif"
        intermediate.write_text(document.as_string())
        metadata = finalize_cif_metadata(intermediate, output_cif)
    expected = {key: row for key, row in before.items() if key[1] == selected_chain}
    proof = compare_atom_semantics(expected, atom_semantics(output_cif))
    if _sha(source_cif) != source_hash:
        raise AssertionError("Source assembly was changed during subset preparation")
    metadata = {key: value for key, value in metadata.items() if key not in ("source_receptor", "source_receptor_sha256")}
    return {"schema_version": 1, "selected_author_chain": selected_chain,
            "source_author_chains": chains, "omitted_author_chains": chains[1:],
            "source_assembly_atom_count": len(before), "selected_atom_count": len(expected),
            "source_receptor_sha256": source_hash, "receptor_sha256": _sha(output_cif),
            "selection_rule": "lexically first prepared author chain; no ligand, site, distance, interface, or outcome selection",
            "native_ligand_coordinate_access": False, "distance_based_selection": False,
            "source_assembly_unchanged": True, "exact_source_atom_subset": True,
            "physically_incomplete_assembly_control": True, "serialization": metadata, **proof}


def _subset_pdb(source, output, selected_chain, expected_atoms):
    """Copy original selected-chain PDB atom lines without coordinate rewriting."""
    from .heldout_repair import parse_pdb_heavy

    if len(selected_chain) != 1:
        raise ValueError("An explicit PDB chain map is required for a multicharacter author chain")
    text = Path(source).read_text()
    lines = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")) and line[21:22] == selected_chain]
    result = "\n".join([*lines, "TER", "END"]) + "\n"
    parsed = parse_pdb_heavy(result)
    expected = {}
    for key, row in expected_atoms.items():
        _, chain, number, insertion, resname, name, _ = key
        expected[(chain, int(number), insertion, resname, name)] = row
    if set(parsed) != set(expected):
        raise ValueError("Selected PDB and CIF source atom identities differ")
    for key, row in parsed.items():
        reference = expected[key]
        if row["element"].upper() != str(reference["type_symbol"]).upper():
            raise ValueError("Selected PDB and CIF elements differ")
        if row["coordinates"] != [float(reference[f"Cartn_{axis}"]) for axis in "xyz"]:
            raise ValueError("Selected PDB and CIF coordinates differ")
    Path(output).write_text(result)
    return {"source_atom_lines_copied_exactly": True, "atom_count": len(parsed), "sha256": _sha(output)}


def materialize_control_case(case_id, primary_case, output):
    """Copy independent chemistry and retain full-assembly and subset maps."""
    if case_id not in CASE_IDS:
        raise ValueError("Only the two predeclared assembly-control cases are eligible")
    primary_case, output = Path(primary_case), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError("Preserve previous control inputs; choose a fresh directory")
    original = {name: _sha(primary_case / name) for name in REQUIRED_FILES}
    output.mkdir(parents=True)
    subset_proof = subset_monomer_input(primary_case / "receptor.cif", output / "receptor.cif")
    selected = subset_proof["selected_author_chain"]
    if selected != "A":
        raise ValueError("Current fixed source selection must remain lexical author chain A for both cases")
    expected_atoms = {key: row for key, row in atom_semantics(primary_case / "receptor.cif").items() if key[1] == selected}
    pdb_proof = _subset_pdb(primary_case / "receptor.pdb", output / "receptor.pdb", selected, expected_atoms)
    full_map = _load(primary_case / "protein_atom_map.json")
    subset_map = copy.deepcopy(full_map)
    subset_map["atoms"] = [row for row in full_map["atoms"] if row["output_key"][0] == selected]
    subset_map["chains"] = [row for row in full_map["chains"] if row["output_chain"] == selected]
    mapped = {tuple(row["output_key"]) for row in subset_map["atoms"]}
    actual = {(key[1], int(key[2]), key[3], key[4], key[5]) for key in expected_atoms}
    if mapped != actual or len(subset_map["chains"]) != 1:
        raise ValueError("Subset source map is not a complete bijection to monomer atom identities")
    shutil.copy2(primary_case / "protein_atom_map.json", output / "full_protein_atom_map.json")
    shutil.copy2(primary_case / "receptor.cif", output / "full_assembly_receptor.cif")
    _write(output / "protein_atom_map.json", subset_map)
    _write(output / "subset_atom_map.json", {"case_id": case_id, "selected_author_chain": selected,
                                           "full_protein_map_sha256": original["protein_atom_map.json"],
                                           "full_assembly_atom_count": subset_proof["source_assembly_atom_count"],
                                           "subset_atom_count": len(mapped), **subset_map})
    for name in LIGAND_FILES:
        shutil.copy2(primary_case / name, output / name)
        if _sha(output / name) != original[name]:
            raise AssertionError("Independent ligand input changed during copying")
    chemistry = _load(output / "chemistry.json")
    manifest = {"schema_version": 1, "status": "prepared", "case_id": case_id,
                "receptor_sha256": _sha(output / "receptor.cif"), "source_receptor_sha256": original["receptor.cif"],
                "atom_map_sha256": _sha(output / "protein_atom_map.json"),
                "protein_heavy_atom_count": len(mapped), "protein_chain_count": 1,
                "scope": "physically incomplete monomer assembly ablation control",
                "selected_author_chain": selected, "native_ligand_access": False, "distance_based_selection": False,
                "no_interface_selection": True, "no_heavy_atom_repair": True,
                "full_assembly_evaluation_receptor": "full_assembly_receptor.cif",
                "full_assembly_evaluation_receptor_sha256": _sha(output / "full_assembly_receptor.cif"),
                "full_assembly_evaluation_map": "full_protein_atom_map.json",
                "prediction_scope_does_not_change_atomic_truth": True}
    _write(output / "receptor.cif.manifest.json", manifest)
    _write(output / "protein_preparation.json", manifest)
    _write(output / "source_symmetry.json", {"orbit_deduplication": "no_orbit_dedupe", "native_ligand_access": False,
                                             "reason": "single retained source chain; no additional symmetry operators proposed"})
    proof = {"schema_version": 1, "status": "prepared", "case_id": case_id,
             "source_only": True, "prediction_outcomes_consulted": False,
             "native_ligand_coordinate_access": False, "fixed_control_denominator_member": True,
             "source_inputs_sha256": original, "subset_proof": subset_proof,
             "pdb_subset_proof": pdb_proof, "protein": manifest,
             "independent_ligand_bytes_and_map_unchanged": True,
             "box_size": chemistry["box_size"], "full_parent_heavy_atoms": chemistry["heavy_atoms"],
             "artifact_sha256": {str(path.relative_to(output)): _sha(path) for path in sorted(output.rglob("*")) if path.is_file()}}
    _write(output / "source_preparation.json", proof)
    if any(_sha(primary_case / name) != expected for name, expected in original.items()):
        raise AssertionError("Source case inputs changed during control preparation")
    return proof


def build_assembly_control(primary, output, *, root=ROOT):
    """Materialize both fixed controls and a pre-outcome proposal; run no models."""
    primary, output, root = Path(primary), Path(output), Path(root)
    if output.exists() or output.is_symlink():
        raise FileExistsError("Preserve previous assembly-control materialization")
    source_receipt = _load(primary / "source_preparation.json")
    if source_receipt.get("serialization_revision") != 2:
        raise ValueError("Assembly controls require the completed source metadata V2 inputs")
    source_cases = {row["case_id"]: row for row in source_receipt["cases"]}
    for case_id in CASE_IDS:
        for name in REQUIRED_FILES:
            if _sha(primary / case_id / name) != source_cases[case_id]["artifact_sha256"][name]:
                raise ValueError(f"Pinned V2 source input changed: {case_id}/{name}")
    source_receipt_sha = _sha(primary / "source_preparation.json")
    output.mkdir(parents=True)
    cases = [materialize_control_case(case_id, primary / case_id, output / case_id) for case_id in CASE_IDS]
    receipt = {"schema_version": 1, "purpose": "source-only assembly ablation mechanical preparation; no predictions, docking, modeling, or metrics",
               "case_ids": list(CASE_IDS), "fixed_control_denominator": 2,
               "primary_source_preparation_sha256": source_receipt_sha,
               "source_ledger_sha256": source_receipt["source_ledger_sha256"],
               "selection_rule": "first lexical prepared author chain; A for both fixed cases",
               "source_only": True, "prediction_outcomes_consulted": False,
               "native_ligand_coordinate_access": False, "model_family_heldout_claim": False,
               "code_sha256": _sha(Path(__file__)),
               "serialization_code_sha256": _sha(Path(__file__).with_name("heldout_serialization.py")), "cases": cases}
    _write(output / "source_preparation.json", receipt)
    proposal = {"schema_version": 1, "status": "proposal_requires_separate_coordinator_freeze_before_new_predictions",
                "purpose": "paired source-only assembly completeness ablation: intact assembly versus first lexical-chain monomer",
                "case_ids": list(CASE_IDS), "source_profile": "current completed metadata V2 full biological assemblies",
                "source_preparation_sha256": _sha(output / "source_preparation.json"),
                "primary_source_preparation_sha256": source_receipt_sha,
                "source_ledger_sha256": source_receipt["source_ledger_sha256"],
                "chain_selection": "A for both cases, chosen as first lexical prepared author chain without ligand/site/interface/outcome access",
                "methods": ["fpocket", "p2rank", "oracle"], "candidate_budgets": [1, 5],
                "prediction": {"scope": "retained complete source chain only; dry protein, no repair",
                               "p2rank_seed": 42, "p2rank_threads": 1, "p2rank_timeout_seconds": 600,
                               "fpocket_timeout_seconds": 3600, "candidate_merge_distance_A": 4., "top_k": 5,
                               "orbit_deduplication": "no_orbit_dedupe"},
                "receptor_policy": {"interaction_cutoff_A": 8., "whole_chain_selection": True,
                                    "no_heavy_atom_repair": True, "no_reference_based_selection_in_predicted_arm": True},
                "ligand_policy": "copy original independently embedded full-parent CCD ligand PDBQT/SDF/chemistry/map bytes exactly; no re-embedding or native coordinates",
                "box_rule": {"side_A": "ceil(max(24,2R+16)) from unchanged independent full-parent conformer", "maximum_side_A": 42,
                             "predicted_center": "predictor-only score-ranked spatially distinct center", "oracle_center": "separate explicitly native-informed matched full-parent reference centroid"},
                "cpu": 2, "exhaustiveness": 16, "num_modes": 20, "min_rmsd_A": .5,
                "energy_range_kcal_mol": 10, "seeds": [104729, 130363, 155921],
                "maximum_method_conditions": 6, "maximum_box_conditions": 22, "maximum_seed_runs": 66,
                "evaluation": {"atomic_reference": "unchanged full-assembly source-qualified copies and canonical scored masks",
                               "clash_denominator": "unchanged full biological assembly, including omitted protein chains",
                               "known_sites": "all uncovered qualified and unqualified known source sites remain annotations; uncovered is not biological nonbinding",
                               "full_vs_core": "preserve original qualification; no score-mask or truth replacement"},
                "inference_scope": "two paired assembly controls in two families; physically incomplete monomers are not independent biological samples",
                "failure_policy": "retain both controls, all methods, unsupported preparation and predictor failures without replacement",
                "locks": {name: _sha(root / name) for name in LOCK_FILES},
                "code_sha256": _sha(Path(__file__)),
                "cases": [{"case_id": case["case_id"], "selected_author_chain": case["subset_proof"]["selected_author_chain"],
                           "source_protein_heavy_atoms": case["subset_proof"]["source_assembly_atom_count"],
                           "control_protein_heavy_atoms": case["subset_proof"]["selected_atom_count"],
                           "full_parent_heavy_atoms": case["full_parent_heavy_atoms"], "box_size": case["box_size"],
                           "receptor_sha256": case["protein"]["receptor_sha256"],
                           "source_inputs_sha256": case["source_inputs_sha256"]} for case in cases]}
    _write(output / "proposal.json", proposal)
    if _sha(primary / "source_preparation.json") != source_receipt_sha:
        raise AssertionError("Source preparation receipt changed")
    return receipt


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, default=ROOT / "data/processed/heldout_structural_v2")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = build_assembly_control(arguments.primary, arguments.output)
    for case in result["cases"]:
        print(f"{case['case_id']}: chain{case['subset_proof']['selected_author_chain']}, {case['subset_proof']['selected_atom_count']} of {case['subset_proof']['source_assembly_atom_count']} protein atoms")
