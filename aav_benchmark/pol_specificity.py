"""Source-pinned POL known-site specificity preparation and execution.

The two free-sugar assay species have four explicit pyranose model states.
Neither independently embedded coordinates nor the melibiose site proxy are
atomic pose truth for those species. No docking is authorized by preparation.
"""
from __future__ import annotations

from collections import Counter
import itertools
import json
import math
from pathlib import Path
import shutil
import tempfile

import gemmi
import numpy as np

from .aav import PROTEIN_RESIDUES, fit
from .external import _ccd_ligand, pdb_heavy_atoms, prepare_external_receptor
from .provenance import sha256, write_json


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "benchmarks/independent_specificity_candidates.json"
ARCHIVE = ROOT / "data/archive/independent_specificity_sources"
COHORT_ID = "perduca_2020_POL_free_hexose_ITC"
EMBEDDING_SEED = 20261024
CA_SHA256 = "5713e078d7604837c4169d281dd40da4dfa21e96340f239776cf32bd3946d65d"
MODEL_SPECS = (
    ("POL_GAL_alpha", "D_galactose", "alpha", "GLA"),
    ("POL_GAL_beta", "D_galactose", "beta", "GAL"),
    ("POL_GLC_alpha", "D_glucose", "alpha", "GLC"),
    ("POL_GLC_beta", "D_glucose", "beta", "BGC"),
)
STATES = ("holo", "apo")
SITES = ("N", "C")
METAL_PROFILES = ("Ca_retained", "dry")
SEEDS = (104729, 130363, 155921)
HEAVY_COUNTS = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    (5, 11, 8, 8, 6, 9, 9, 4, 10, 8, 8, 9, 8, 11, 7, 6, 7, 14, 12, 7)))


def _canonical(molecule):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)


def load_source_qualification(ledger_path=LEDGER, source_dir=ARCHIVE):
    """Use the actual candidate_id; reject source mutation and chemical aliases."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    ledger_path, source_dir = Path(ledger_path), Path(source_dir)
    ledger = json.loads(ledger_path.read_text())
    matches = [r for r in ledger["candidates"] if r["candidate_id"] == COHORT_ID]
    if len(matches) != 1 or matches[0]["qualification"].get("source_qualification_frozen") is not True:
        raise ValueError("POL requires its unique frozen source qualification")
    source = matches[0]
    species = source["qualification"]["recommended_primary_inventory"]
    if species != ["D_galactose", "D_glucose"]:
        raise ValueError("POL source species inventory changed")
    manifest = json.loads((source_dir / "manifest.json").read_text())
    entries = {r["path"]: r for r in manifest["source_files"]}
    required = {"6T1D.cif", "6T0Q.cif", "GLA.cif", "GAL.cif", "GLC.cif", "BGC.cif"}
    if not required <= entries.keys():
        raise ValueError("POL source archive lacks a required coordinate/CCD source")
    for name in required:
        if sha256(source_dir / name) != entries[name]["sha256"]:
            raise ValueError(f"POL source checksum mismatch: {name}")
    audit_entries = [r for r in manifest["derived_files"] if r["path"] == "POL_site_audit.json"]
    if len(audit_entries) != 1 or sha256(source_dir / "POL_site_audit.json") != audit_entries[0]["sha256"]:
        raise ValueError("POL source site-audit checksum mismatch")
    graphs = {r["ccd"]: r for r in source["chemistry"]["model_graphs"]}
    if set(graphs) != {r[3] for r in MODEL_SPECS}:
        raise ValueError("POL requires exactly the four frozen anomer model graphs")
    for ccd, row in graphs.items():
        if sha256(source_dir / f"{ccd}.cif") != row["cif_sha256"]:
            raise ValueError("POL graph CCD checksum differs from source ledger")
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if (mol is None or len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumHeavyAtoms() != 12
                or Chem.GetFormalCharge(mol) != 0 or rdMolDescriptors.CalcMolFormula(mol) != "C6H12O6"):
            raise ValueError("POL model graph changed formula, size, charge or connectivity")
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(centers) != 5 or any(s == "?" for _, s in centers):
            raise ValueError("POL model needs five assigned stereocentres")
    measured = {r["species_id"]: r for r in source["source_tested_inventory"]}
    if (measured["D_galactose"]["ITC_outcome"] != "measured_binding"
            or measured["D_glucose"]["ITC_outcome"] != "no_binding_detected"
            or measured["D_glucose"].get("Kd_M") is not None):
        raise ValueError("POL operational assay labels changed or invented a glucose Kd")
    return source, graphs


def experiment_template(ledger_path=LEDGER, source_dir=ARCHIVE):
    """Concrete reviewable design; deliberately not an execution permission."""
    load_source_qualification(ledger_path, source_dir)
    source_dir = Path(source_dir)
    return {
        "schema_version": 1, "experiment_id": "pol_specificity_v1", "cohort_id": COHORT_ID,
        "analysis_frozen_before_new_docking": False,
        "source_ledger_sha256": sha256(ledger_path),
        "source_manifest_sha256": sha256(source_dir / "manifest.json"),
        "site_audit_sha256": sha256(source_dir / "POL_site_audit.json"),
        "calcium_ccd_sha256": CA_SHA256,
        "implementation_files_sha256": {str(path.relative_to(ROOT)): sha256(path)
                                         for path in (Path(__file__), Path(__file__).with_name("external.py"),
                                                      Path(__file__).with_name("aav.py"), Path(__file__).with_name("specificity.py"),
                                                      Path(__file__).with_name("engine.py"), Path(__file__).with_name("metrics.py"),
                                                      Path(__file__).with_name("provenance.py"))},
        "model_ids": [s[0] for s in MODEL_SPECS], "states": list(STATES), "sites": list(SITES),
        "metal_profiles": list(METAL_PROFILES), "seeds": list(SEEDS),
        "embedding_seed": EMBEDDING_SEED, "box_size": [24.0, 24.0, 24.0],
        "alternate_conformer_policy": "whole-residue highest mean occupancy, ties lexical; apo Ser218 selects B (0.58), retains blank atoms, clears selected altloc and keeps source occupancy",
        "apo_alignment_policy": "whole dimer matched N/CA/C/O; choose A/B chain correspondence by lowest backbone RMSD, ties lexical; transform centers only",
        "water_policy": "all explicit waters omitted in both profiles; dry additionally removes native calcium",
        "cpu": 2, "exhaustiveness": 16, "num_modes": 20,
        "min_rmsd_A": 0.5, "energy_range_kcal_mol": 10,
        "primary_condition": {"state": "holo", "site": "N", "metal_profile": "Ca_retained"},
        "primary_aggregation": "seed mean within anomer, then equal mean of the two anomers per species",
        "anomer_sensitivity": "minimum and maximum anomer seed means per species",
        "tie_rule": "exact equal aggregated scores receive half concordance credit; lower Vina score is favorable",
        "anomer_population_interpretation": "equal model weighting is a heuristic; solution fractions, furanose and open-chain populations are unknown",
        "information_regime": "true-site-local; ligand-informed prior sites; not blind",
        "expected_runs": 96,
    }


def validate_experiment(config, ledger_path=LEDGER, source_dir=ARCHIVE, *, execution=False):
    """Do not silently narrow denominators or alter the reviewed design."""
    expected = experiment_template(ledger_path, source_dir)
    for key, value in expected.items():
        if key == "analysis_frozen_before_new_docking":
            continue
        if config.get(key) != value:
            raise ValueError(f"POL experiment differs from reviewed design: {key}")
    if type(config.get("analysis_frozen_before_new_docking")) is not bool:
        raise ValueError("POL freeze flag must be an explicit boolean")
    if execution and config["analysis_frozen_before_new_docking"] is not True:
        raise ValueError("POL docking requires a separate frozen experiment before execution")
    ca = Path(source_dir) / "CA.cif"
    if not ca.is_file() or sha256(ca) != CA_SHA256:
        raise ValueError("POL preparation requires the pinned CA.cif source")
    return config


def planned_cases(config):
    """Sixteen engine-compatible cases, each with six receptor/seed runs."""
    return [{"case_id": f"{model}__{state}__{site}", "model_id": model, "species_id": species,
             "anomer": anomer, "ccd_id": ccd, "state": state, "site": site}
            for (model, species, anomer, ccd), state, site in itertools.product(MODEL_SPECS, STATES, SITES)]


def _build_receptor(source_dir, state):
    """Expand declared assembly once; rename copied protein chains explicitly."""
    pdb_id = {"holo": "6T1D", "apo": "6T0Q"}[state]
    path = Path(source_dir) / f"{pdb_id}.cif"
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("POL structure needs exactly one deposited model")
    assemblies = [a for a in structure.assemblies if a.name == "1"]
    if len(assemblies) != 1:
        raise ValueError("POL requires declared biological assembly 1")
    assembly = assemblies[0]
    operator_count = sum(len(g.operators) for g in assembly.generators)
    if operator_count != (1 if state == "holo" else 2):
        raise ValueError("POL biological assembly operators changed")
    model = gemmi.make_assembly(assembly, structure[0], gemmi.HowToNameCopiedChain.AddNumber)
    protein_chains = [c for c in model if any(r.entity_type == gemmi.EntityType.Polymer for r in c)]
    if len(protein_chains) != 2 or len({c.name for c in protein_chains}) != 2:
        raise ValueError("POL assembly must contain exactly two protein chains")
    mapping = dict(zip(sorted(c.name for c in protein_chains), ("A", "B")))
    new_model = gemmi.Model("1")
    for target in ("A", "B"):
        new_model.add_chain(gemmi.Chain(target))
    atom_map, omitted, backbone, alternatives = [], Counter(), {}, []
    residues_per_chain, heavy_per_chain, ions_per_chain, raw_heavy_per_chain = Counter(), Counter(), Counter(), Counter()
    for chain in model:
        for residue in chain:
            protein = residue.entity_type == gemmi.EntityType.Polymer
            calcium = residue.name == "CA"
            if not (protein or calcium):
                omitted[residue.name] += sum(not a.element.is_hydrogen for a in residue)
                continue
            if chain.name not in mapping:
                raise ValueError("POL calcium is not assigned to a selected protein chain")
            target_chain = mapping[chain.name]
            new_residue = residue.clone()
            for i in reversed(range(len(new_residue))):
                if new_residue[i].element.is_hydrogen:
                    del new_residue[i]
            if protein:
                raw_heavy_per_chain[target_chain] += len(new_residue)
            altlocs = sorted({a.altloc for a in new_residue if a.altloc not in ("\x00", " ")})
            if altlocs:
                if state != "apo" or residue.name != "SER" or residue.seqid.num != 218 or altlocs != ["A", "B"]:
                    raise ValueError("POL uncurated alternate conformer")
                occupancies = {alt: [float(a.occ) for a in new_residue if a.altloc == alt] for alt in altlocs}
                selected = min(altlocs, key=lambda alt: (-float(np.mean(occupancies[alt])), alt))
                if selected != "B" or not np.allclose(occupancies["B"], [0.58] * 3):
                    raise ValueError("POL frozen apo alternative occupancy changed")
                alternatives.append({"assembly_copy_chain": chain.name, "author_residue": 218,
                                     "selected": selected, "source_occupancies": occupancies,
                                     "rule": "one whole residue conformer by highest mean occupancy; lexical tie break",
                                     "retained_shared_atom_names": [a.name for a in new_residue if a.altloc in ("\x00", " ")],
                                     "retained_selected_atom_names": [a.name for a in new_residue if a.altloc == selected],
                                     "discarded_atom_names": [a.name for a in new_residue if a.altloc not in ("\x00", " ", selected)],
                                     "selected_altloc_cleared": True, "occupancy_normalized": False})
                for i in reversed(range(len(new_residue))):
                    if new_residue[i].altloc not in ("\x00", " ", selected):
                        del new_residue[i]
                for atom in new_residue:
                    if atom.altloc == selected:
                        atom.altloc = "\x00"
            if any(a.altloc not in ("\x00", " ") for a in new_residue):
                raise ValueError("POL alternate conformer was not resolved")
            if not altlocs and any(a.occ != 1 for a in new_residue):
                raise ValueError("POL uncurated partial occupancy atom")
            if protein:
                if residue.name not in PROTEIN_RESIDUES:
                    raise ValueError("POL unsupported protein chemistry")
                names = {a.name for a in new_residue}
                expected_heavy = HEAVY_COUNTS[residue.name] + (1 if "OXT" in names and residue.seqid.num == 373 else 0)
                if len(names) != len(new_residue) or len(names) != expected_heavy:
                    raise ValueError(f"POL incomplete/duplicate heavy atoms: {chain.name}:{residue.seqid}")
                residues_per_chain[target_chain] += 1
                heavy_per_chain[target_chain] += len(new_residue)
                for a in new_residue:
                    if a.name in {"N", "CA", "C", "O"}:
                        backbone[(target_chain, residue.seqid.num, residue.name, a.name)] = list(a.pos)
            else:
                if len(new_residue) != 1 or new_residue[0].element.name != "Ca" or new_residue[0].name != "CA":
                    raise ValueError("POL native calcium identity changed")
                ions_per_chain[target_chain] += 1
                source_key = [chain.name, residue.seqid.num, residue.seqid.icode.strip(), "CA", "CA"]
                new_residue.name, new_residue[0].name = "CA2", "CA+2"
                atom_map.append({"assembly_copy_source_key": source_key,
                                 "prepared_key": [target_chain, residue.seqid.num, residue.seqid.icode.strip(), "CA2", "CA+2"],
                                 "coordinates_A": list(new_residue[0].pos), "formal_charge": 2,
                                 "mapping_reason": "native CA/CA mapped to explicit standalone CA2/CA+2 for the audited Meeko ion helper"})
            new_model[target_chain].add_residue(new_residue)
    nres, natoms = ((338, 2517) if state == "holo" else (340, 2529))
    if (dict(residues_per_chain) != {"A": nres, "B": nres}
            or dict(heavy_per_chain) != {"A": natoms, "B": natoms}
            or dict(ions_per_chain) != {"A": 2, "B": 2}):
        raise ValueError("POL assembly residue/heavy-atom/calcium denominator changed")
    result = gemmi.Structure()
    result.add_model(new_model)
    text = result.make_pdb_string()
    # Parse our actual serialized artifact; no lost/duplicated atoms at conversion.
    parsed = pdb_heavy_atoms(text)
    if len(parsed) != 2 * natoms + 4:
        raise ValueError("POL mmCIF-to-PDB conversion lost heavy atoms")
    audit = {"state": state, "pdb_id": pdb_id, "source_sha256": sha256(path),
             "assembly_id": "1", "materializations": 1, "operator_count": operator_count,
             "chain_copy_map": mapping, "protein_residues_per_chain": dict(residues_per_chain),
             "assembly_operators": [{"name": op.name, "rotation": np.asarray(op.transform.mat).tolist(),
                                     "translation_A": list(op.transform.vec)}
                                    for generator in assembly.generators for op in generator.operators],
             "chain_operator_correspondence": ([{"prepared_chain": "A", "source_author_chain": "A", "operator": "1"},
                                                 {"prepared_chain": "B", "source_author_chain": "B", "operator": "1"}]
                                                if state == "holo" else
                                                [{"prepared_chain": "A", "source_author_chain": "A", "operator": "1"},
                                                 {"prepared_chain": "B", "source_author_chain": "A", "operator": "2"}]),
             "raw_protein_heavy_atoms_per_chain_before_altloc_selection": dict(raw_heavy_per_chain),
             "protein_heavy_atoms_per_chain": dict(heavy_per_chain), "native_calcium_count": 4,
             "source_ion_to_prepared_label_map": atom_map, "removed_nonprotein_heavy_atoms_by_resname": dict(omitted),
             "alternate_conformers": alternatives,
             "observed_residue_heavy_atom_completion": False, "protein_heavy_atom_minimization": False,
             "unmodeled_N_terminal_residues_per_chain": 18 if state == "holo" else 13,
             "apo_context_caveat": "native glycosylated apo protein; covalent glycans removed explicitly" if state == "apo" else None}
    return text, backbone, audit


def receptor_sources(source_dir=ARCHIVE):
    """Keep native coordinate frames; transfer centers using protein alone."""
    source_dir = Path(source_dir)
    holo, hb, ha = _build_receptor(source_dir, "holo")
    apo, ab, aa = _build_receptor(source_dir, "apo")
    candidates = []
    for chain_map in ({"A": "A", "B": "B"}, {"A": "B", "B": "A"}):
        common = sorted(k for k in hb if (chain_map[k[0]], *k[1:]) in ab)
        if len(common) != 2 * 338 * 4:
            raise ValueError("POL holo/apo backbone identity correspondence changed")
        source_xyz = np.array([hb[k] for k in common])
        target_xyz = np.array([ab[(chain_map[k[0]], *k[1:])] for k in common])
        rotation, translation = fit(source_xyz, target_xyz)
        rmsd = float(np.sqrt(np.mean(np.sum((source_xyz @ rotation + translation - target_xyz) ** 2, axis=1))))
        candidates.append((rmsd, tuple(chain_map.items()), rotation, translation))
    rmsd, selected_map, rotation, translation = min(candidates, key=lambda row: (row[0], row[1]))
    audit = json.loads((source_dir / "POL_site_audit.json").read_text())
    selected = {row["glycan_label_asym_id"]: row for row in audit["assembly_1_sites"]}
    centers = {"holo": {"N": selected["G"]["centroid_GLAnonH_A"],
                         "C": selected["H"]["centroid_GLAnonH_A"]}}
    if selected["G"]["protein_author_chain"] != "A" or selected["H"]["protein_author_chain"] != "A":
        raise ValueError("POL known-site chain assignment changed")
    centers["apo"] = {site: (np.asarray(point) @ rotation + translation).tolist()
                      for site, point in centers["holo"].items()}
    alignment = {"direction": "holo center to independently solved native apo coordinate frame",
                 "atom_selection": "N/CA/C/O of exact chain, author residue number, residue identity matches in A/B",
                 "matched_backbone_atoms": len(common), "rotation_row_vector": rotation.tolist(),
                 "holo_to_apo_chain_correspondence": dict(selected_map),
                 "candidate_chain_correspondences": [{"chain_map": dict(row[1]), "backbone_rmsd_A": row[0]}
                                                     for row in candidates],
                 "translation_A": translation.tolist(), "backbone_rmsd_A": rmsd,
                 "native_ligand_coordinates_used_for_alignment": False,
                 "receptor_coordinates_transformed": False, "centers": centers}
    return {"holo": (holo, ha), "apo": (apo, aa)}, centers, alignment


def _ligand_models(source_dir, graphs):
    from .specificity import _check_stereo
    result = {}
    for model, species, anomer, ccd in MODEL_SPECS:
        molecule, atom_map, chemistry = _ccd_ligand(source_dir, [ccd], [], EMBEDDING_SEED)
        expected = graphs[ccd]["canonical_smiles"]
        if _canonical(molecule) != expected:
            raise ValueError("POL CCD chemistry differs from the frozen exact source graph")
        _check_stereo(molecule, expected)
        chemistry.update({"model_id": model, "species_id": species, "anomer": anomer,
                          "ccd_id": ccd, "native_coordinates_used_for_embedding": False,
                          "stereochemistry_verified_from_3d": True, "pose_reference_eligible": False,
                          "source_ccd_sha256": sha256(Path(source_dir) / f"{ccd}.cif"),
                          "experimental_anomer_fraction": None, "other_solution_forms": "unknown"})
        result[model] = (molecule, atom_map, chemistry)
    return result


def _write_ligand(molecule, atom_map, center, output):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from .specificity import _check_stereo

    molecule = Chem.Mol(molecule)
    expected = _canonical(molecule)
    positions = molecule.GetConformer().GetPositions()
    indices = [row["rdkit_index"] for row in atom_map]
    positions += np.asarray(center) - positions[indices].mean(0)
    for i, xyz in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(i, xyz)
    _check_stereo(molecule, expected)
    setup = MoleculePreparation().prepare(molecule)
    if len(setup) != 1:
        raise ValueError("POL ligand unexpectedly has multiple preparation setups")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup[0])
    if not ok:
        raise ValueError(error)
    names = [row["name"] for row in atom_map]
    parsed = {line[12:16].strip(): [float(line[30:38]), float(line[38:46]), float(line[46:54])]
              for line in pdbqt.splitlines() if line.startswith(("ATOM  ", "HETATM"))
              and line[77:].strip() not in {"H", "HD", "HS"}}
    if len(parsed) != 12 or set(parsed) != set(names):
        raise ValueError("POL ligand PDBQT changed its heavy atom map")
    rounded = Chem.Mol(molecule)
    for row in atom_map:
        rounded.GetConformer().SetAtomPosition(row["rdkit_index"], parsed[row["name"]])
    _check_stereo(rounded, expected)
    (output / "ligand.pdbqt").write_text(pdbqt)
    writer = Chem.SDWriter(str(output / "ligand.sdf"))
    writer.write(molecule)
    writer.close()
    loaded = list(Chem.SDMolSupplier(str(output / "ligand.sdf"), removeHs=False))
    if len(loaded) != 1 or loaded[0] is None:
        raise ValueError("POL SDF failed round-trip")
    _check_stereo(loaded[0], expected)
    return names, positions[indices].tolist()


def _prepare_receptors(source_dir, output):
    sources, centers, alignment = receptor_sources(source_dir)
    result = {}
    for state, (text, audit) in sources.items():
        for profile in METAL_PROFILES:
            stem = output / f"{state}_{profile}"
            if profile == "Ca_retained":
                definition = {"receptor_cofactors": [{"ccd_id": "CA", "source_resname": "CA2", "expected_count": 4}]}
                actual = text
            else:
                definition = {}
                actual = "\n".join(line for line in text.splitlines()
                                   if not (line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CA2")) + "\n"
            log = prepare_external_receptor(actual, stem, definition, source_dir)
            expected = sum(audit["protein_heavy_atoms_per_chain"].values()) + (4 if profile == "Ca_retained" else 0)
            if log["pdbqt_heavy_atoms"] != expected:
                raise ValueError("POL receptor PDBQT heavy atom denominator changed")
            result[(state, profile)] = {"stem": stem, "log": log, "source_audit": audit}
    return result, centers, alignment


def preflight(config_path, *, ledger_path=LEDGER, source_dir=ARCHIVE):
    """Exercise actual four ligand and four receptor preparations; never dock."""
    config = validate_experiment(json.loads(Path(config_path).read_text()), ledger_path, source_dir)
    _, graphs = load_source_qualification(ledger_path, source_dir)
    models = _ligand_models(source_dir, graphs)
    with tempfile.TemporaryDirectory(prefix="pol-preflight-") as temporary:
        temp = Path(temporary)
        receptors, centers, alignment = _prepare_receptors(source_dir, temp)
        ligands = []
        for model, (molecule, atom_map, chemistry) in models.items():
            directory = temp / model
            directory.mkdir()
            _write_ligand(molecule, atom_map, centers["holo"]["N"], directory)
            ligands.append(chemistry)
        receptor_audits = [{"state": state, "metal_profile": profile,
                            "source_audit": row["source_audit"],
                            "prepared_heavy_atoms": row["log"]["pdbqt_heavy_atoms"],
                            "prepared_calcium_count": len(row["log"].get("native_cofactors", [])),
                            "max_heavy_displacement_A": row["log"]["max_heavy_displacement_A"]}
                           for (state, profile), row in receptors.items()]
    return {"status": "preflight_passed_no_docking", "cohort_id": COHORT_ID,
            "experiment_sha256": sha256(config_path), "cases": planned_cases(config),
            "model_states": ligands, "receptors": receptor_audits, "alignment": alignment,
            "cost_denominators": {"experimental_species": 2, "independent_families": 1,
                                  "independent_epimer_contrasts": 1, "independently_embedded_model_states": 4,
                                  "prepared_receptor_profiles": 4, "engine_case_directories": 16,
                                  "seed_runs": 96, "primary_runs": 12, "sensitivity_runs": 84,
                                  "maximum_returned_poses": 1920, "cpus_per_run": 2,
                                  "per_run_timeout_seconds": 3600,
                                  "wall_time_estimate": "unmeasured before docking; 96 times observed mean run duration for serial execution"}}


def prepare_pol_panel(config_path, output_root, *, ledger_path=LEDGER, source_dir=ARCHIVE):
    """Prepare all 16 cases; enforce a fresh directory and immutable source pins."""
    config_path, output = Path(config_path).resolve(), Path(output_root).resolve()
    config = validate_experiment(json.loads(config_path.read_text()), ledger_path, source_dir)
    source, graphs = load_source_qualification(ledger_path, source_dir)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("POL preparation output must be fresh")
    write_json(output / "preparation_attempt.json", {"status": "running", "experiment_sha256": sha256(config_path)})
    try:
        cache = output / "receptors"
        cache.mkdir()
        receptors, centers, alignment = _prepare_receptors(source_dir, cache)
        models = _ligand_models(source_dir, graphs)
        rows = []
        for case in planned_cases(config):
            directory = output / "cases" / case["case_id"]
            directory.mkdir(parents=True)
            center = centers[case["state"]][case["site"]]
            molecule, atom_map, chemistry = models[case["model_id"]]
            names, coordinates = _write_ligand(molecule, atom_map, center, directory)
            entries = []
            for profile in METAL_PROFILES:
                row = receptors[(case["state"], profile)]
                stem = row["stem"]
                for suffix in (".pdb", ".pdbqt", "_atom_changes.json"):
                    filename = stem.name + suffix
                    shutil.copyfile(stem.with_name(filename), directory / filename)
                entries.append({"id": profile, "state": case["state"],
                                "assembly_condition": "complete_biological_dimer_" + profile,
                                "pdb_file": stem.name + ".pdb", "pdbqt_file": stem.name + ".pdbqt",
                                "atom_changes_file": stem.name + "_atom_changes.json"})
            reference = {"names": names, "coordinates": coordinates, "atom_map": atom_map,
                         "pose_reference_eligible": False,
                         "coordinate_role": "independently generated input only; not experimental atomic pose truth"}
            metadata = {"schema_version": 1, **case, "family_id": "POL_CBM67", "cohort_id": COHORT_ID,
                        "status": "prepared", "execution_information_regime": "true-site-local",
                        "information_regime": config["information_regime"], "box_center": center,
                        "box_size": config["box_size"], "ligand_file": "ligand.pdbqt",
                        "reference_file": "reference.json", "receptors": entries, "chemistry": chemistry,
                        "experiment_sha256": sha256(config_path), "source_ledger_sha256": sha256(ledger_path),
                        "assay_species_id": case["species_id"], "experimental_anomer_label_available": False,
                        "operational_glucose_NB_is_numeric_Kd_bound": False}
            write_json(directory / "reference.json", reference)
            write_json(directory / "chemistry.json", chemistry)
            write_json(directory / "source_graph.json", {"source_candidate_id": source["candidate_id"],
                       "ccd_graph": graphs[case["ccd_id"]], "source_ledger_sha256": sha256(ledger_path), "atom_map": atom_map})
            write_json(directory / "preparation.json", metadata)
            rows.append({**case, "prepared_dir": str(directory.relative_to(output)),
                         "preparation_sha256": sha256(directory / "preparation.json"),
                         "prepared_artifact_sha256": {p.name: sha256(p) for p in sorted(directory.iterdir()) if p.is_file()}})
        result = {"schema_version": 1, "status": "prepared", "cohort_id": COHORT_ID,
                  "experiment_sha256": sha256(config_path), "source_ledger_sha256": sha256(ledger_path),
                  "cases": rows, "expected_runs": 96, "alignment": alignment,
                  "source_receptor_audits": {state: receptors[(state, "Ca_retained")]["source_audit"] for state in STATES}}
        write_json(output / "panel_preparation.json", result)
    except Exception as error:
        write_json(output / "preparation_attempt.json", {"status": "failed", "reason": str(error), "error_type": type(error).__name__})
        raise
    write_json(output / "preparation_attempt.json", {"status": "ok", "experiment_sha256": sha256(config_path)})
    return result


def dock_pol_panel(config_path, prepared_root, output_root, executable, *, ledger_path=LEDGER,
                   source_dir=ARCHIVE, resume=False, case_ids=None, seeds=None):
    """Delegate real Vina invocations only after an explicit complete freeze."""
    from .engine import dock

    config_path, prepared, output = Path(config_path).resolve(), Path(prepared_root).resolve(), Path(output_root).resolve()
    config = validate_experiment(json.loads(config_path.read_text()), ledger_path, source_dir, execution=True)
    panel = json.loads((prepared / "panel_preparation.json").read_text())
    if panel["experiment_sha256"] != sha256(config_path) or panel["cohort_id"] != COHORT_ID:
        raise ValueError("POL prepared panel differs from the frozen experiment")
    expected = {r["case_id"] for r in planned_cases(config)}
    if {r["case_id"] for r in panel["cases"]} != expected or len(panel["cases"]) != 16:
        raise ValueError("POL prepared case denominator changed")
    selected = expected if case_ids is None else set(case_ids)
    if not selected or not selected <= expected or (case_ids is not None and len(selected) != len(case_ids)):
        raise ValueError("POL execution subset contains missing, duplicate or unknown cases")
    result = []
    for row in panel["cases"]:
        if row["case_id"] not in selected:
            continue
        directory = prepared / row["prepared_dir"]
        if sha256(directory / "preparation.json") != row["preparation_sha256"]:
            raise ValueError("POL case preparation changed")
        for name, expected_hash in row["prepared_artifact_sha256"].items():
            if Path(name).name != name or sha256(directory / name) != expected_hash:
                raise ValueError("POL prepared artifact changed before execution")
        runs = dock(directory, output / row["case_id"], executable, config_path, resume=resume, seeds=seeds)
        result.append({"case_id": row["case_id"], "runs": len(runs),
                       "successful_runs": sum(r["status"] == "ok" for r in runs)})
    return result


def summarize_pol_scores(run_scores):
    """Describe exactly 96 completed top-1 scores without inventing assay units.

    Input rows contain case_id, metal_profile, seed and top1_score. This pure
    aggregation does not read results or calibrate affinity. The caller must
    separately verify raw runs, pose maps and frozen input fingerprints.
    """
    cases = {r["case_id"]: r for r in planned_cases({})}
    expected = set(itertools.product(cases, METAL_PROFILES, SEEDS))
    observed = {}
    for row in run_scores:
        identity = (row["case_id"], row["metal_profile"], row["seed"])
        if identity not in expected or identity in observed or row.get("status", "ok") != "ok":
            raise ValueError("POL aggregation needs one successful score for each frozen run")
        score = row["top1_score"]
        if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
            raise ValueError("POL aggregation requires finite numeric top-1 scores")
        observed[identity] = float(score)
    if set(observed) != expected:
        raise ValueError("POL aggregation requires the complete 96-run denominator")
    conditions = []
    for state, site, profile in itertools.product(STATES, SITES, METAL_PROFILES):
        species = {}
        for species_id in ("D_galactose", "D_glucose"):
            anomers = []
            for model, measured_species, anomer, ccd in MODEL_SPECS:
                if measured_species != species_id:
                    continue
                case_id = f"{model}__{state}__{site}"
                scores = [observed[(case_id, profile, seed)] for seed in SEEDS]
                anomers.append({"model_id": model, "anomer": anomer, "ccd_id": ccd,
                                "seeds": list(SEEDS), "seed_top1_scores_kcal_mol": scores,
                                "seed_mean_kcal_mol": math.fsum(scores) / len(scores)})
            means = [a["seed_mean_kcal_mol"] for a in anomers]
            species[species_id] = {"anomers": anomers,
                                  "equal_anomer_mean_kcal_mol": math.fsum(means) / len(means),
                                  "minimum_anomer_seed_mean_kcal_mol": min(means),
                                  "maximum_anomer_seed_mean_kcal_mol": max(means)}
        delta = (species["D_glucose"]["equal_anomer_mean_kcal_mol"]
                 - species["D_galactose"]["equal_anomer_mean_kcal_mol"])
        conditions.append({"state": state, "site": site, "metal_profile": profile,
                           "species": species, "glucose_minus_galactose_kcal_mol": delta,
                           "measured_binder_lower_score": delta > 0,
                           "exact_aggregate_tie": delta == 0,
                           "descriptive_single_pair_concordance": 1.0 if delta > 0 else 0.0 if delta < 0 else 0.5})
    primary = next(r for r in conditions if (r["state"], r["site"], r["metal_profile"]) == ("holo", "N", "Ca_retained"))
    return {"schema_version": 1, "cohort_id": COHORT_ID, "seed_runs": 96,
            "experimental_species": 2, "experimental_epimer_contrasts": 1, "independent_families": 1,
            "primary": primary, "conditions": conditions,
            "primary_aggregation": "seed mean within anomer, then equal mean of two anomers per species",
            "interpretation": "equal-anomer weighting is a computational heuristic; unknown solution populations, furanose and open-chain fractions are not estimated",
            "sensitivity_scope": "all site/state/calcium conditions remain separate; no pooled assay replicates",
            "affinity_calibration": False, "auroc": None, "auprc": None,
            "generalization_scope": "one descriptive natural epimer contrast in one additional protein family"}
