"""New structural test inputs, prepared directly from pinned source chemistry.

Source preparation is deliberately separate from scientific execution.  No old
site predictions, prepared ligands, docking poses, or accuracy tables are read.
The fixed cohort is a structural test; it is not a model-family-heldout cohort.
"""
from __future__ import annotations

from collections import defaultdict
from decimal import Decimal
import gzip
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

import gemmi
import numpy as np

from .aav import PROTEIN_RESIDUES, _pdb_atoms, box_distance, prepare_receptor
from .blind_docking import ligand_atoms, translate_ligand
from .provenance import sha256, write_json, write_text_atomic
from .site_geometry import distinct_candidates, read_receptor

CASE_IDS = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
LEDGER = "benchmarks/heldout_reference_ledger.json"
LEDGER_SHA256 = "5e1142c49191ee9881efc1e918630abc0f57cffb5f2fb0041350cd6df0ac7c4e"
INVENTORY = "benchmarks/heldout_candidate_inventory.json"
INVENTORY_SHA256 = "11c185ea93bf637ddc4fe4818aa4c2114a87f83e34d7d40c2cd01759f9280458"
EMBEDDING_SEED = 20261024
INTERACTION_CUTOFF_A = 8.0
METHODS = ("fpocket", "p2rank", "oracle")
LOCKS = ("benchmarks/site_methods.lock.json", "benchmarks/p2rank.lock.json", "benchmarks/vina.lock.json")


def _load(path):
    return json.loads(Path(path).read_text())


def _blank(value):
    return "" if value in (None, False, ".", "?", " ") else str(value)


def _fresh(path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Preserve prior experiment artifacts; choose a fresh path: {path}")
    path.mkdir(parents=True)
    return path


def _pin(root, record):
    relative = record.get("archive_path", record.get("archive"))
    path = (Path(root) / relative).resolve()
    if not path.is_relative_to(Path(root).resolve()) or not path.is_file():
        raise ValueError(f"Missing or escaped pinned source: {relative}")
    if sha256(path) != record.get("archive_sha256", record["sha256"]):
        raise ValueError(f"Pinned source changed: {relative}")
    if record.get("archive_encoding", "").startswith("gzip mtime"):
        if hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() != record["sha256"]:
            raise ValueError(f"Uncompressed pinned source changed: {relative}")
    return path


def source_cases(root):
    """Verify the complete fixed denominator and return source-only records."""
    root = Path(root)
    for relative, expected in ((LEDGER, LEDGER_SHA256), (INVENTORY, INVENTORY_SHA256)):
        if sha256(root / relative) != expected:
            raise ValueError(f"Frozen structural source changed: {relative}")
    ledger = _load(root / LEDGER)
    if not ledger.get("source_only") or ledger.get("prediction_outcomes_consulted") is not False:
        raise ValueError("Structural ledger is not source-only")
    by_id = {row["case_id"]: row for row in ledger["cases"]}
    if len(by_id) != 12 or ledger["cohort_denominator"] != 12:
        raise ValueError("Keep the original twelve screened candidates in the denominator")
    cases = []
    for case_id in CASE_IDS:
        case = by_id[case_id]
        if not case["source_conditioned_preparation_eligible"] or not case["atomic_pose_eligible"]:
            raise ValueError(f"Fixed runnable case lost source eligibility: {case_id}")
        entry = case["inventory_entry"]
        _pin(root, entry["sources"]["coordinates"])
        for record in entry["ccd_sources"].values():
            _pin(root, record)
        cases.append(case)
    return cases


def select_residue_altloc(rows):
    """Select a single whole-residue alternate, sharing blank-alt atoms.

    Occupancy means include the common blank-alt atoms for each candidate.
    Ties use the lexical alternate ID; atoms from different alternatives are
    never combined.  A partial chosen alternative is not repaired.
    """
    common = [row for row in rows if not _blank(row.get("label_alt_id"))]
    alternatives = sorted({_blank(row.get("label_alt_id")) for row in rows} - {""})
    candidates = alternatives or [""]
    materialized = {alt: common + [r for r in rows if _blank(r.get("label_alt_id")) == alt]
                    if alt else common for alt in candidates}
    exact_means = {alt: sum((Decimal(str(r["occupancy"])) for r in selected), Decimal(0)) / len(selected)
             for alt, selected in materialized.items()}
    chosen = min(candidates, key=lambda alt: (-exact_means[alt], alt))
    means = {alt: float(value) for alt, value in exact_means.items()}
    selected = materialized[chosen]
    names = [r["label_atom_id"] for r in selected]
    if len(names) != len(set(names)):
        raise ValueError("Whole-residue alternative contains duplicate atom names")
    if any(not math.isfinite(float(r["occupancy"])) or float(r["occupancy"]) <= 0 for r in selected):
        raise ValueError("Protein has nonpositive or nonfinite occupancy")
    return selected, {"selected_alt_id": chosen, "candidate_mean_occupancies": means,
                      "shared_blank_atom_count": len(common), "source_rows": len(rows),
                      "selected_rows": len(selected), "tie_rule": "lexical alternate ID"}


def _category_rows(category):
    return [dict(zip(category, values)) for values in zip(*category.values())]


def prepare_source_protein(root, case, output):
    """Materialize every protein copy in the selected deposited assembly.

    This retains protein chemical modifications, strips all nonprotein species,
    and records a bijection to source author/label identities and rigid copies.
    It neither consults ligand coordinates nor selects a native neighborhood.
    """
    output = Path(output)
    entry = case["inventory_entry"]
    source = _pin(root, entry["sources"]["coordinates"])
    block = gemmi.cif.read_file(str(source))[0]
    poly = block.get_mmcif_category("_entity_poly.")
    protein_entities = {e for e, kind in zip(poly["entity_id"], poly["type"])
                        if str(kind).startswith("polypeptide")}
    atoms = _category_rows(block.get_mmcif_category("_atom_site."))
    model_ids = {_blank(row["pdbx_PDB_model_num"]) for row in atoms}
    if model_ids != {"1"}:
        raise ValueError("An explicit single deposited model is required")
    protein = [r for r in atoms if r["label_entity_id"] in protein_entities and r["type_symbol"] not in ("H", "D")]
    grouped = defaultdict(list)
    for row in protein:
        key = tuple(_blank(row.get(k)) for k in ("label_asym_id", "label_seq_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_comp_id"))
        grouped[key].append(row)
    selected, alt_audit = {}, []
    for key, rows in grouped.items():
        picked, audit = select_residue_altloc(rows)
        selected[key] = picked
        if audit["source_rows"] != audit["selected_rows"]:
            alt_audit.append({"source_residue_key": list(key), **audit})
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    atom_map, chain_map = [], []
    chain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    operations = case["receptor_scope"]["expanded_assembly_operators"]
    seen_copies = set()
    for operation in operations:
        matrix = np.asarray(operation["matrix"], float)
        vector = np.asarray(operation["vector"], float)
        if (matrix.shape != (3, 3) or vector.shape != (3,) or not np.isfinite(matrix).all()
                or not np.isfinite(vector).all() or not np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-6)
                or not np.isclose(np.linalg.det(matrix), 1, atol=1e-6)):
            raise ValueError("Assembly protein operator is not a proper rigid transform")
        for label_chain in operation["subchains"]:
            residues = [(key, rows) for key, rows in selected.items() if key[0] == label_chain]
            if not residues:
                continue
            copy_key = (label_chain, tuple(matrix.ravel()), tuple(vector))
            if copy_key in seen_copies:
                continue
            seen_copies.add(copy_key)
            if len(chain_map) >= len(chain_alphabet):
                raise ValueError("Assembly exceeds the explicit PDB chain-map capacity; no truncation allowed")
            output_chain = chain_alphabet[len(chain_map)]
            chain = gemmi.Chain(output_chain)
            chain_map.append({"output_chain": output_chain, "source_label_chain": label_chain,
                              "source_author_chains": sorted({key[2] for key, _ in residues}),
                              "operator_name": operation["operator_name"], "matrix": matrix.tolist(), "vector": vector.tolist()})
            for key, rows in residues:
                _, label_number, author_chain, author_number, insertion, resname = key
                number = int(label_number or author_number)
                residue = gemmi.Residue()
                residue.name, residue.seqid, residue.het_flag = resname, gemmi.SeqId(number, " "), "A"
                residue.subchain = output_chain
                for row in rows:
                    source_xyz = np.asarray([float(row[f"Cartn_{axis}"]) for axis in "xyz"])
                    xyz = matrix @ source_xyz + vector
                    if not np.isfinite(xyz).all():
                        raise ValueError("Nonfinite protein coordinate")
                    atom = gemmi.Atom()
                    atom.name, atom.element = row["label_atom_id"], gemmi.Element(row["type_symbol"])
                    atom.pos, atom.occ, atom.b_iso = gemmi.Position(*xyz), float(row["occupancy"]), float(row["B_iso_or_equiv"])
                    residue.add_atom(atom)
                    atom_map.append({"output_key": [output_chain, number, "", resname, atom.name],
                                     "source_atom_site_id": str(row["id"]), "source_model": "1",
                                     "source_label_chain": label_chain, "source_label_resnum": label_number,
                                     "source_author_chain": author_chain, "source_author_resnum": author_number,
                                     "source_insertion_code": insertion, "source_resname": resname,
                                     "source_atom_name": atom.name, "source_alt_id": _blank(row["label_alt_id"]),
                                     "element": atom.element.name, "operator_name": operation["operator_name"],
                                     "source_xyz_A": source_xyz.tolist(), "assembly_xyz_A": xyz.tolist()})
                chain.add_residue(residue)
            structure[0].add_chain(chain)
    if not atom_map or len({tuple(r["output_key"]) for r in atom_map}) != len(atom_map):
        raise ValueError("Empty or ambiguous assembly protein identity map")
    if len(atom_map) > 99999:
        raise ValueError("Full assembly exceeds PDB atom capacity; no truncation allowed")
    pdb = structure.make_pdb_string()
    parsed = _pdb_atoms(pdb)
    for row in atom_map:
        if tuple(row["output_key"]) not in parsed:
            raise ValueError("Assembly PDB serialization lost a protein atom")
    delta = max(float(np.linalg.norm(np.asarray(parsed[tuple(row["output_key"])]["coordinates"]) - row["assembly_xyz_A"])) for row in atom_map)
    if delta > 0.0009 or len(parsed) != len(atom_map):
        raise ValueError("Assembly PDB serialization changed protein coordinates or identities")
    output.mkdir(parents=True, exist_ok=True)
    write_text_atomic(output / "receptor.pdb", pdb)
    structure.make_mmcif_document().write_file(str(output / "receptor.cif"))
    write_json(output / "protein_atom_map.json", {"atoms": atom_map, "chains": chain_map})
    modifications = sorted({row["source_resname"] for row in atom_map} - PROTEIN_RESIDUES)
    manifest = {"schema_version": 1, "case_id": case["case_id"], "source_sha256": sha256(source),
                "selected_assembly_id": case["receptor_scope"]["selected_assembly_id"],
                "scope": "all protein subchains and rigid copies of selected deposited biological assembly",
                "protein_heavy_atom_count": len(atom_map), "protein_chain_count": len(chain_map),
                "whole_residue_alternate_selection": alt_audit, "protein_modifications_preserved": modifications,
                "model": "uniform dry protein only; no heavy atom repair, mutation, or native crop",
                "all_nonprotein_components_removed": True, "native_ligand_access": False,
                "distance_based_selection": False, "receptor_sha256": sha256(output / "receptor.cif"),
                "atom_map_sha256": sha256(output / "protein_atom_map.json"),
                "max_PDB_serialization_displacement_A": delta}
    write_json(output / "receptor.cif.manifest.json", manifest)
    write_json(output / "protein_preparation.json", manifest)
    return manifest


def _ccd(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = {}
    for row in _category_rows(block.get_mmcif_category("_chem_comp_atom.")):
        if row["type_symbol"] in ("H", "D"):
            continue
        atoms[row["atom_id"]] = {"element": row["type_symbol"], "charge": int(row["charge"]),
                                  "stereo": row["pdbx_stereo_config"],
                                  "aromatic": row.get("pdbx_aromatic_flag") == "Y",
                                  "ideal": [float(row[f"pdbx_model_Cartn_{a}_ideal"]) for a in "xyz"]}
    bonds = [(r["atom_id_1"], r["atom_id_2"], r["value_order"])
             for r in _category_rows(block.get_mmcif_category("_chem_comp_bond."))
             if r["atom_id_1"] in atoms and r["atom_id_2"] in atoms]
    if not atoms or not bonds:
        raise ValueError("Empty pinned CCD chemical graph")
    return {"atoms": atoms, "bonds": bonds}


def _signed_volume(center, neighbors):
    vectors = np.asarray(neighbors) - np.asarray(center)
    matrix = vectors if len(vectors) == 3 else vectors[:3] - vectors[3]
    return float(np.linalg.det(matrix))


def check_signed_chirality(atom_map, coordinates, stereo_tests):
    """Check named-neighbor handedness, including declared leaving replacements."""
    xyz = dict(zip([row["canonical_atom_key"] for row in atom_map], np.asarray(coordinates, float)))
    results = []
    for test in stereo_tests:
        volume = _signed_volume(xyz[test["center"]], [xyz[n] for n in test["neighbors"]])
        passed = math.isfinite(volume) and volume * test["ideal_signed_volume"] > 0
        results.append({**test, "observed_signed_volume": volume, "status": "consistent" if passed else "inverted_or_degenerate"})
    if not results or any(row["status"] != "consistent" for row in results):
        raise ValueError("CCD signed named-neighbor chirality mismatch")
    return results


def build_source_ligand(root, case, *, embedding_seed=EMBEDDING_SEED):
    """Build the full parent graph from CCD and declared chemical branch links.

    Only pinned CCD ideal coordinates establish tetrahedral identity.  No
    experimental ligand coordinates are read by this chemical constructor.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if embedding_seed != EMBEDDING_SEED:
        raise ValueError("Use the fixed independent embedding seed")
    entry = case["inventory_entry"]
    if entry.get("extra_covalent_edges"):
        raise ValueError("External covalent parent requires an explicit source chemical schema")
    components = {name: _ccd(_pin(root, source)) for name, source in entry["ccd_sources"].items()}
    residues = {str(row["branch_num"]): row["ccd_id"] for row in entry["residues"]}
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    complete = Chem.RWMol()
    full_indices, ideal = {}, {}
    for number, component_id in residues.items():
        component = components[component_id]
        for name, definition in component["atoms"].items():
            key = f"{number}:{component_id}:{name}"
            atom = Chem.Atom(definition["element"])
            atom.SetFormalCharge(definition["charge"])
            atom.SetIsAromatic(definition["aromatic"])
            full_indices[key] = complete.AddAtom(atom)
            ideal[key] = definition["ideal"]
        for left, right, order in component["bonds"]:
            complete.AddBond(full_indices[f"{number}:{component_id}:{left}"], full_indices[f"{number}:{component_id}:{right}"], orders[order.upper()])
    complete = complete.GetMol()
    Chem.SanitizeMol(complete)
    conformer = Chem.Conformer(complete.GetNumAtoms())
    for key, index in full_indices.items():
        conformer.SetAtomPosition(index, ideal[key])
    complete.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(complete, replaceExistingTags=True)
    full_keys = {index: key for key, index in full_indices.items()}
    replacements, links = {}, []
    for row in entry["branch_links"]:
        if row["value_order"].upper() != "SING":
            raise ValueError("Unsupported declared branch bond order")
        endpoints = [f"{row[f'entity_branch_list_num_{i}']}:{row[f'comp_id_{i}']}:{row[f'atom_id_{i}']}" for i in (1, 2)]
        links.append(endpoints)
        for i, endpoint in enumerate(endpoints, 1):
            leaving = f"{row[f'entity_branch_list_num_{i}']}:{row[f'comp_id_{i}']}:{row[f'leaving_atom_id_{i}']}"
            if leaving in full_indices:
                replacements[leaving] = endpoints[2-i]
    editable, indices, atom_map = Chem.RWMol(), {}, []
    for key, original_index in full_indices.items():
        if key in replacements:
            continue
        number, component_id, atom_name = key.split(":")
        atom = Chem.Atom(complete.GetAtomWithIdx(original_index))
        index = editable.GetNumAtoms()
        name = f"L{index:03d}"
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4)); info.SetResidueName(component_id)
        info.SetResidueNumber(int(number)); info.SetChainId("L")
        atom.SetMonomerInfo(info)
        atom.SetProp("canonical_atom_key", key)
        indices[key] = editable.AddAtom(atom)
        definition = components[component_id]["atoms"][atom_name]
        atom_map.append({"name": name, "rdkit_index": index, "canonical_atom_key": key,
                         "branch_num": number, "ccd_id": component_id, "atom_name": atom_name,
                         "element": definition["element"], "formal_charge": definition["charge"],
                         "ccd_stereo_config": definition["stereo"]})
    for bond in complete.GetBonds():
        left, right = full_keys[bond.GetBeginAtomIdx()], full_keys[bond.GetEndAtomIdx()]
        if left in indices and right in indices:
            editable.AddBond(indices[left], indices[right], bond.GetBondType())
    for left, right in links:
        editable.AddBond(indices[left], indices[right], Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    final_keys = {index: key for key, index in indices.items()}
    tests = []
    for key, index in indices.items():
        original = complete.GetAtomWithIdx(full_indices[key])
        number, component_id, atom_name = key.split(":")
        stereo = components[component_id]["atoms"][atom_name]["stereo"]
        if stereo not in ("R", "S"):
            continue
        desired = [replacements.get(full_keys[a.GetIdx()], full_keys[a.GetIdx()]) for a in original.GetNeighbors()]
        actual = [final_keys[a.GetIdx()] for a in molecule.GetAtomWithIdx(index).GetNeighbors()]
        if set(desired) != set(actual) or len(actual) not in (3, 4):
            raise ValueError("Declared linkage changed a CCD stereocenter unexpectedly")
        permutation = [desired.index(neighbor) for neighbor in actual]
        inversions = sum(permutation[i] > permutation[j] for i in range(len(permutation)) for j in range(i+1, len(permutation)))
        if inversions % 2:
            molecule.GetAtomWithIdx(index).InvertChirality()
        neighbors = sorted(full_keys[a.GetIdx()] for a in original.GetNeighbors())
        volume = _signed_volume(ideal[key], [ideal[n] for n in neighbors])
        tests.append({"center": key, "ccd_stereo_config": stereo,
                      "neighbors": [replacements.get(n, n) for n in neighbors],
                      "ideal_neighbor_keys": neighbors, "ideal_signed_volume": volume,
                      "declared_leaving_replacements": {n: replacements[n] for n in neighbors if n in replacements}})
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Full parent ligand graph is disconnected")
    if any(s == "?" for _, s in Chem.FindMolChiralCenters(molecule, includeUnassigned=True, useLegacyImplementation=False)):
        raise ValueError("Full parent ligand has unassigned stereochemistry")
    smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    molecule.RemoveAllConformers()
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3(); parameters.randomSeed = embedding_seed; parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise ValueError("Independent full-parent ETKDG embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(molecule) or AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000) != 0:
        raise ValueError("Independent full-parent MMFF optimization unsupported or unconverged")
    xyz = molecule.GetConformer().GetPositions()[:len(atom_map)]
    stereo_result = check_signed_chirality(atom_map, xyz, tests)
    radius = float(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1).max())
    side = math.ceil(max(24, 2 * radius + 16))
    if side > 42:
        raise ValueError("Independent full-parent conformer exceeds the frozen 42 A box limit")
    chemistry = {"heavy_atoms": len(atom_map), "isomeric_smiles": smiles,
                 "formal_charge": Chem.GetFormalCharge(molecule), "embedding_seed": embedding_seed,
                 "embedding": "ETKDGv3 one thread then MMFF maxIters1000",
                 "native_coordinates_used_for_embedding": False, "full_parent_graph": True,
                 "protonation": "pinned CCD formal charges plus RDKit valence hydrogens",
                 "stereochemistry": "pinned CCD ideal named-neighbor parity with declared branch replacements",
                 "signed_chirality_tests": stereo_result, "box_size": [side] * 3,
                 "heavy_atom_centroid_radius_A": radius,
                 "box_rule": "cube side ceil(max(24,2R+16)) A from independent conformer, maximum42 A",
                 "ccd_sources": entry["ccd_sources"], "branch_links": entry["branch_links"],
                 "versions": {p: importlib.metadata.version(p) for p in ("rdkit", "meeko", "gemmi")}}
    return molecule, atom_map, chemistry


def _serialize_ligand(molecule, atom_map, output, stem="ligand"):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1:
        raise ValueError("Expected one full-parent ligand setup")
    text, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(error)
    observed = {a["name"]: a for a in ligand_atoms(text) if not a["hydrogen"]}
    expected = {r["name"] for r in atom_map}
    if set(observed) != expected:
        raise ValueError("Ligand PDBQT serialization changed canonical heavy atom identities")
    xyz = molecule.GetConformer().GetPositions()
    maximum = max(float(np.linalg.norm(xyz[row["rdkit_index"]] - observed[row["name"]]["coordinates"])) for row in atom_map)
    if maximum > .0009:
        raise ValueError("Ligand PDBQT serialization moved heavy atoms")
    output = Path(output)
    write_text_atomic(output / f"{stem}.pdbqt", text)
    writer = Chem.SDWriter(str(output / f"{stem}.sdf")); writer.write(molecule); writer.close()
    return {"pdbqt_sha256": sha256(output / f"{stem}.pdbqt"), "sdf_sha256": sha256(output / f"{stem}.sdf"),
            "max_heavy_serialization_displacement_A": maximum}


def preflight_receptor(directory):
    """Attempt full-chain chemical parameterization without launching a predictor."""
    directory = Path(directory)
    output = directory / "receptor_parameterization"
    output.mkdir(exist_ok=True)
    protein = _load(directory / "protein_preparation.json")
    if protein["protein_modifications_preserved"]:
        return {"status": "unsupported", "modifications": protein["protein_modifications_preserved"],
                "reason": "No pinned modified-residue Meeko template and finite validated partial-charge/atom-type model; source elements and coordinates preserved without conversion",
                "no_heavy_atom_repair": True}
    try:
        result = prepare_receptor((directory / "receptor.pdb").read_text(), output / "receptor")
        return {"status": "prepared", "parameterization": result,
                "pdb_sha256": sha256(output / "receptor.pdb"), "pdbqt_sha256": sha256(output / "receptor.pdbqt")}
    except (ValueError, RuntimeError, OSError) as error:
        return {"status": "preparation_failed", "reason": str(error), "error_type": type(error).__name__,
                "no_heavy_atom_repair": True}


def finalize_preflight(root, prepared):
    """Add full-protein support, source symmetry and native serialization audits."""
    prepared = Path(prepared)
    receipt = _load(prepared / "source_preparation.json")
    if receipt["source_ledger_sha256"] != LEDGER_SHA256 or receipt["case_ids"] != list(CASE_IDS):
        raise ValueError("Preflight inputs differ from fixed structural sources")
    for row in receipt["cases"]:
        directory = prepared / row["case_id"]
        row["receptor_parameterization"] = preflight_receptor(directory)
        write_json(directory / "receptor_parameterization.json", row["receptor_parameterization"])
        mapping = _load(directory / "protein_atom_map.json")
        _, certificate = certified_protein_operators(directory / "receptor.cif", mapping["chains"])
        write_json(directory / "source_symmetry.json", certificate)
        native = directory / "native_control"
        if native.exists():
            raise FileExistsError("Preserve native controls; finalize preflight only once")
        row["native_control"] = native_ligand_control(root, row["case_id"], directory, native)
        row["artifact_sha256"] = {str(p.relative_to(directory)): sha256(p) for p in sorted(directory.rglob("*"))
                                   if p.is_file() and p.name != "source_preparation.json"}
        write_json(directory / "source_preparation.json", row)
        write_json(prepared / "source_preparation.json", receipt)
    return receipt


def prepare_sources(root, output):
    """Perform only mechanical source/chemical preparation for all fixed six."""
    output = _fresh(output)
    rows = []
    for case in source_cases(root):
        directory = output / case["case_id"]; directory.mkdir()
        row = {"case_id": case["case_id"], "status": "preparing", "qualification": case["qualification"],
               "scored_atom_count": case["scored_atom_count"], "fixed_denominator_member": True}
        try:
            row["protein"] = prepare_source_protein(root, case, directory)
            molecule, atom_map, chemistry = build_source_ligand(root, case)
            chemistry["serialization"] = _serialize_ligand(molecule, atom_map, directory)
            write_json(directory / "ligand_atom_map.json", {"atoms": atom_map})
            write_json(directory / "chemistry.json", chemistry)
            row.update(status="prepared", full_parent_heavy_atoms=len(atom_map), box_size=chemistry["box_size"],
                       docking_receptor_support="standard_protein" if not row["protein"]["protein_modifications_preserved"] else "requires_audited_modified_residue_parameters")
        except (ValueError, RuntimeError, OSError, KeyError) as error:
            row.update(status="preparation_failed", reason=str(error), error_type=type(error).__name__)
        write_json(directory / "source_preparation.json", row)
        rows.append(row)
        receipt = {"schema_version": 1, "purpose": "mechanical source and independent chemical preparation; no scientific predictions or docking",
                   "case_ids": list(CASE_IDS), "screened_denominator": 12, "runnable_structural_denominator": 6,
                   "model_family_heldout_claim": False, "source_ledger_sha256": LEDGER_SHA256,
                   "inventory_sha256": INVENTORY_SHA256, "prediction_outcomes_consulted": False,
                   "code_sha256": sha256(Path(__file__)), "cases": rows}
        write_json(output / "source_preparation.json", receipt)
    return receipt


def native_ligand_control(root, case_id, prepared, output):
    """Preserve deposited heavy coordinates and audit stereo/serialization.

    This is explicitly native-informed and separate from prediction inputs.
    No scoring or optimization is performed; hydrogens alone are modeled.
    """
    from rdkit import Chem

    case = next(c for c in source_cases(root) if c["case_id"] == case_id)
    prepared, output = Path(prepared), _fresh(output)
    molecules = list(Chem.SDMolSupplier(str(prepared / "ligand.sdf"), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Expected one independently prepared full-parent SDF")
    heavy = Chem.RemoveHs(molecules[0])
    atom_map = _load(prepared / "ligand_atom_map.json")["atoms"]
    chemistry = _load(prepared / "chemistry.json")
    native = {f"{r['branch_num']}:{r['ccd_id']}:{r['atom_name']}": r["xyz_A"] for r in case["native_atom_map"]}
    if set(native) != {r["canonical_atom_key"] for r in atom_map}:
        raise ValueError("Native full-parent atom map is incomplete or ambiguous")
    for row in atom_map:
        heavy.GetConformer().SetAtomPosition(row["rdkit_index"], native[row["canonical_atom_key"]])
    checks = check_signed_chirality(atom_map, heavy.GetConformer().GetPositions(), chemistry["signed_chirality_tests"])
    molecule = Chem.AddHs(heavy, addCoords=True)
    # SDF does not retain PDB monomer names: restore names explicitly by map.
    hydrogen_number = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            row = atom_map[atom.GetIdx()]; name, comp, num = row["name"], row["ccd_id"], int(row["branch_num"])
        else:
            hydrogen_number += 1; name, comp, num = f"H{hydrogen_number}", "LIG", 1
        info = Chem.AtomPDBResidueInfo(); info.SetName(name.rjust(4)); info.SetResidueName(comp)
        info.SetResidueNumber(num); info.SetChainId("L"); atom.SetMonomerInfo(info)
    serialization = _serialize_ligand(molecule, atom_map, output, "native_ligand")
    receipt = {"schema_version": 1, "case_id": case_id, "native_coordinate_access": True,
               "purpose": "native-heavy-coordinate stereochemistry and serialization control; not docking success",
               "source_ledger_sha256": LEDGER_SHA256, "heavy_atoms": len(atom_map),
               "heavy_coordinates_optimized": False, "hydrogens": "RDKit AddHs(addCoords=True); no optimization",
               "signed_chirality_tests": checks, "serialization": serialization}
    write_json(output / "native_control.json", receipt)
    return receipt


def planning_inventory(root, prepared=None):
    cases = source_cases(root)
    rows = []
    for case in cases:
        entry = case["inventory_entry"]
        row = {"case_id": case["case_id"], "qualification": case["qualification"],
               "full_parent_heavy_atoms": entry["expected_heavy_atom_count"],
               "scored_heavy_atoms": case["scored_atom_count"],
               "qualified_assembly_copies": len(case["qualified_copies"]),
               "assembly_id": case["receptor_scope"]["selected_assembly_id"],
               "coordinate_source": entry["sources"]["coordinates"], "ccd_sources": entry["ccd_sources"]}
        if prepared is not None:
            directory = Path(prepared) / case["case_id"]
            row["mechanical_preparation"] = _load(directory / "source_preparation.json")
            row["prepared_source_sha256"] = {p.name: sha256(p) for p in directory.iterdir() if p.is_file()}
        rows.append(row)
    return {"schema_version": 1, "screened_denominator": 12, "fixed_structural_test_denominator": 6,
            "fixed_case_ids": list(CASE_IDS), "model_family_heldout_claim": False,
            "source_ledger_sha256": LEDGER_SHA256, "inventory_sha256": INVENTORY_SHA256,
            "prediction_outcomes_consulted": False, "methods": list(METHODS),
            "maximum_predicted_boxes_per_method": 5, "oracle_boxes_per_case": 1,
            "maximum_method_conditions": 18, "maximum_box_conditions": 66,
            "interaction_cutoff_A": INTERACTION_CUTOFF_A, "cases": rows,
            "locks": {relative: sha256(Path(root) / relative) for relative in LOCKS}}


def verify_freeze(root, config_path, freeze_path, prepared=None):
    """Require the coordinator's explicit pre-outcome freeze before execution."""
    config_path, freeze_path = Path(config_path), Path(freeze_path)
    config, freeze = _load(config_path), _load(freeze_path)
    source_cases(root)
    if (freeze.get("configuration_sha256") != sha256(config_path)
            or freeze.get("frozen_before_new_docking") is not True
            or freeze.get("frozen_before_new_predictions") is not True
            or freeze.get("source_ledger_sha256") != LEDGER_SHA256):
        raise ValueError("Execution requires a pre-prediction/pre-docking freeze bound to the source ledger")
    ids = config.get("primary_case_ids", config.get("cohort", {}).get("case_ids"))
    if ids != list(CASE_IDS):
        raise ValueError("Configuration must retain exactly the fixed six structural cases")
    if (config.get("prediction", {}).get("p2rank_seed", 42) != 42
            or config.get("prediction", {}).get("candidate_merge_distance_A", 4.) != 4.
            or config.get("prediction", {}).get("top_k", 5) != 5
            or config.get("ligand", {}).get("embedding_seed", EMBEDDING_SEED) != EMBEDDING_SEED
            or config.get("receptor_policy", {}).get("interaction_cutoff_A") != INTERACTION_CUTOFF_A):
        raise ValueError("Configuration differs from fixed predictor, embedding, merge, or 8A interaction rules")
    for relative in LOCKS:
        if freeze.get("lock_sha256", {}).get(relative) != sha256(Path(root) / relative):
            raise ValueError(f"Method lock changed after freeze: {relative}")
    if prepared is not None:
        path = Path(prepared) / "source_preparation.json"
        if freeze.get("source_preparation_sha256") != sha256(path):
            raise ValueError("Mechanical source preparation changed after freeze")
        for row in _load(path)["cases"]:
            directory = Path(prepared) / row["case_id"]
            for relative, expected in row.get("artifact_sha256", {}).items():
                artifact = (directory / relative).resolve()
                if not artifact.is_relative_to(directory.resolve()) or sha256(artifact) != expected:
                    raise ValueError(f"Mechanical prepared artifact changed after freeze: {row['case_id']}/{relative}")
    return config


def certified_protein_operators(receptor, chain_map, tolerance_A=.05):
    """Certify source-declared operations against the entire protein atom set.

    Only operations actually declared by assembly construction are candidates;
    no ligand positions or fitted homologous chains propose extra symmetry.
    Atom chemical identity and a bijective whole-assembly map are required.
    """
    from scipy.spatial import cKDTree

    keys, xyz, elements = read_receptor(receptor)
    classes = defaultdict(list)
    for i, (key, element) in enumerate(zip(keys, elements)):
        classes[(key[1:], element)].append(i)
    operations, audits = [], []
    proposals = {(tuple(np.asarray(c["matrix"]).ravel()), tuple(c["vector"])) for c in chain_map}
    for flat, vector in sorted(proposals):
        matrix, vector = np.asarray(flat).reshape(3, 3), np.asarray(vector)
        transformed = xyz @ matrix.T + vector
        maximum, valid = 0., True
        for indices in classes.values():
            distances, nearest = cKDTree(xyz[indices]).query(transformed[indices])
            maximum = max(maximum, float(np.max(distances)))
            if maximum > tolerance_A or len(set(nearest.tolist())) != len(indices):
                valid = False; break
        audits.append({"matrix": matrix.tolist(), "vector": vector.tolist(),
                       "whole_protein_atom_bijection": valid, "max_displacement_A": maximum})
        if valid:
            operations.append((matrix.T, vector))
    nonidentity = any(not np.allclose(r, np.eye(3)) or not np.allclose(t, 0) for r, t in operations)
    return operations, {"source": "declared assembly operators certified using the complete protein atom map",
                        "native_ligand_access": False, "tolerance_A": tolerance_A, "operations": audits,
                        "orbit_deduplication": "protein_certified" if nonidentity else "no_orbit_dedupe"}


def select_candidates(predictions, receptor, protein_map):
    operations, certificate = certified_protein_operators(receptor, protein_map["chains"])
    candidates, mapping = distinct_candidates(predictions, operations, merge_distance_A=4.)
    return {"raw_candidate_count": len(predictions), "distinct_candidate_count": len(candidates),
            "candidates": candidates[:5], "candidate_mapping": mapping,
            "candidate_merge_distance_A": 4., "maximum_candidates": 5,
            "symmetry_certificate": certificate, "native_coordinate_access": False}


def run_predictor(root, prepared, case_id, method, output, config_path, freeze_path, *, executable):
    """Run one frozen predictor on the full prepared biological assembly."""
    verify_freeze(root, config_path, freeze_path, prepared)
    if case_id not in CASE_IDS or method not in ("fpocket", "p2rank"):
        raise ValueError("Unknown fixed case or predictor")
    prepared = Path(prepared) / case_id
    receptor = prepared / "receptor.cif"
    if method == "p2rank":
        from .p2rank import run_p2rank
        result = run_p2rank(receptor, output, executable)
    else:
        result = _run_fpocket_full_protein(root, receptor, output, executable)
    completed = result["status"] in ("ok", "completed")
    selection = select_candidates(result.get("predictions", []) if completed else [], receptor,
                                  _load(prepared / "protein_atom_map.json"))
    record_name = "result.json" if method == "p2rank" else "execution.json"
    selection.update(case_id=case_id, method=method, predictor_status=result["status"],
                     prediction_input_sha256=sha256(receptor), configuration_sha256=sha256(config_path),
                     prediction_record={"filename": record_name, "sha256": sha256(Path(output) / record_name)},
                     selection_code_sha256=sha256(Path(__file__)),
                     status=("ok" if selection["candidates"] else "completed_no_candidates") if completed else "predictor_failed")
    write_json(Path(output) / "candidate_selection.json", selection)
    return selection


def _run_fpocket_full_protein(root, receptor, output, executable):
    """Pinned fpocket with source polymer modifications retained and audited."""
    from .sites import audit_fpocket_receptor, parse_fpocket

    output = _fresh(output).resolve(); executable = Path(executable).resolve()
    lock_path = Path(root) / LOCKS[0]
    config = _load(lock_path)["methods"]["fpocket"]
    installation = _load(executable.parent.parent / "installation.json")
    if (any(installation.get(k) != config[k] for k in ("release", "commit", "source_sha256"))
            or installation["executable_sha256"] != sha256(executable)):
        raise ValueError("fpocket installation differs from frozen source/build")
    local = output / "receptor.cif"; shutil.copy2(receptor, local)
    command = [str(executable), "-f", local.name, *config["parameters"]]
    result = {"schema_version": 1, "method": "fpocket", "status": "failed", "predictions": [],
              "command": command, "input_sha256": sha256(local), "method_lock_sha256": sha256(lock_path),
              "native_coordinate_access": False, "scope": "entire source protein biological assembly; modifications retained"}
    start = time.perf_counter()
    try:
        with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
            execution = subprocess.run(command, cwd=output, stdout=stdout, stderr=stderr, timeout=3600, check=False)
        result["returncode"] = execution.returncode
        if execution.returncode:
            raise RuntimeError(f"fpocket exited {execution.returncode}")
        raw = output / "receptor_out"
        info = raw / "receptor_info.txt"
        no_pockets = "no pockets found" in (output / "stdout.log").read_text().splitlines()
        if not info.is_file() and not no_pockets:
            raise RuntimeError("fpocket returned zero without descriptor or explicit no-pockets evidence")
        if info.is_file():
            result["full_receptor_audit"] = audit_fpocket_receptor(raw / "receptor_out.cif", local)
            result["predictions"] = parse_fpocket(raw, local)
            described = len(re.findall(r"(?m)^Pocket\s+\d+\s*:", info.read_text()))
            if described != len(result["predictions"]):
                raise ValueError("fpocket normalized predictions and descriptor count differ")
        else:
            result["full_receptor_audit"] = {"status": "unavailable", "reason": "upstream writes no full receptor when no pockets are found"}
        if sha256(receptor) != result["input_sha256"] or sha256(local) != result["input_sha256"]:
            raise ValueError("Full assembly input changed during fpocket execution")
        result["no_pockets_reported"] = no_pockets
        result["status"] = "completed"
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        result.update(error_type=type(error).__name__, reason=str(error))
    result["runtime_seconds"] = time.perf_counter() - start
    result["artifacts"] = {str(p.relative_to(output)): {"sha256": sha256(p), "bytes": p.stat().st_size}
                           for p in sorted(output.rglob("*")) if p.is_file() and p.name != "execution.json"}
    write_json(output / "execution.json", result)
    return result


def neighborhood_protein(receptor, center, size):
    """Retain complete chains within the Vina 8 A box interaction cutoff."""
    keys, xyz, _ = read_receptor(receptor)
    distances = box_distance(xyz, center, size)
    chains = sorted({key[0] for key, d in zip(keys, distances) if d <= INTERACTION_CUTOFF_A})
    if not chains:
        raise ValueError("No complete protein chain intersects box interaction neighborhood")
    selected = np.asarray([key[0] in chains for key in keys])
    omitted = distances[~selected]
    closest = float(omitted.min()) if len(omitted) else None
    if closest is not None and closest <= INTERACTION_CUTOFF_A:
        raise ValueError("Omitted protein atom lies within the Vina interaction cutoff")
    structure = gemmi.read_structure(str(receptor))
    out = gemmi.Structure(); out.add_model(gemmi.Model("1"))
    for chain in structure[0]:
        if chain.name in chains:
            out[0].add_chain(chain.clone())
    return out.make_pdb_string(), {"selected_chains": chains, "selected_heavy_atoms": int(selected.sum()),
                                  "source_heavy_atoms": len(keys), "interaction_cutoff_A": INTERACTION_CUTOFF_A,
                                  "closest_omitted_atom_to_box_A": closest,
                                  "rule": "complete protein chains with any heavy atom within8 A of closed search box"}


def prepare_condition(root, prepared, case_id, method, output, config_path, freeze_path, *, predictions=None, oracle_center=None):
    """Prepare top-five predicted boxes or the explicitly supplied oracle box."""
    verify_freeze(root, config_path, freeze_path, prepared)
    if case_id not in CASE_IDS or method not in METHODS:
        raise ValueError("Unknown fixed structural condition")
    directory, output = Path(prepared) / case_id, _fresh(output)
    manifest = {"schema_version": 1, "case_id": case_id, "condition_id": case_id, "complex_id": case_id,
                "method": method, "configuration_sha256": sha256(config_path), "status": "preparing", "jobs": [],
                "native_coordinate_access": method == "oracle", "fixed_denominator_member": True}
    if method == "oracle":
        center = np.asarray(oracle_center, float)
        if center.shape != (3,) or not np.isfinite(center).all() or predictions is not None:
            raise ValueError("Oracle requires only an explicit native-informed center")
        candidates = [{"rank": 1, "original_rank": 1, "center_A": center.tolist()}]
    else:
        if oracle_center is not None:
            raise ValueError("Native center forbidden in predicted condition")
        selection = _load(Path(predictions) / "candidate_selection.json")
        if (selection["case_id"] != case_id or selection["method"] != method
                or selection["configuration_sha256"] != sha256(config_path)
                or selection["prediction_input_sha256"] != sha256(directory / "receptor.cif")):
            raise ValueError("Frozen predictor selection does not match case/config/full protein")
        record = selection["prediction_record"]
        if Path(record["filename"]).name != record["filename"] or sha256(Path(predictions) / record["filename"]) != record["sha256"]:
            raise ValueError("Predictor record changed after candidate selection")
        predictor = _load(Path(predictions) / record["filename"])
        replay = select_candidates(predictor.get("predictions", []) if predictor["status"] in ("ok", "completed") else [],
                                   directory / "receptor.cif", _load(directory / "protein_atom_map.json"))
        if replay["candidates"] != selection["candidates"] or replay["candidate_mapping"] != selection["candidate_mapping"]:
            raise ValueError("Candidate centers/ranks differ from source-only selection replay")
        candidates = selection["candidates"]
        manifest["candidate_selection"] = selection
    chemistry = _load(directory / "chemistry.json")
    atom_map = _load(directory / "ligand_atom_map.json")["atoms"]
    for candidate in candidates:
        box_id = f"box{candidate['rank']:02d}"; box = output / box_id; box.mkdir()
        job = {"id": box_id, "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"], "center_A": candidate["center_A"], "prepared_directory": box_id}
        try:
            protein, neighborhood = neighborhood_protein(directory / "receptor.cif", candidate["center_A"], chemistry["box_size"])
            audit = prepare_receptor(protein, box / "neighborhood")
            moved, translation = translate_ligand((directory / "ligand.pdbqt").read_text(), candidate["center_A"])
            write_text_atomic(box / "ligand.pdbqt", moved)
            write_json(box / "reference.json", {"names": [r["name"] for r in atom_map], "canonical_atom_map": atom_map,
                                                 "pose_reference_eligible": False, "purpose": "canonical atom identities only; atomic truth evaluated separately"})
            prep = {"case_id": case_id, "family_id": case_id, "status": "prepared", "cohort_role": "new_structural_test",
                    "information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                    "execution_information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                    "box_center": candidate["center_A"], "box_size": chemistry["box_size"], "ligand_file": "ligand.pdbqt",
                    "receptors": [{"id": "neighborhood", "state": "holo", "assembly_condition": "complete_8A_box_neighborhood", "pdb_file": "neighborhood.pdb", "pdbqt_file": "neighborhood.pdbqt"}],
                    "configuration_sha256": sha256(config_path), "source_ledger_sha256": LEDGER_SHA256,
                    "neighborhood": neighborhood, "ligand_translation": translation, "receptor_preparation": audit,
                    "native_coordinates_in_runner_reference": False, "model_family_heldout_claim": False}
            write_json(box / "preparation.json", prep)
            job.update(status="prepared", preparation_sha256=sha256(box / "preparation.json"))
        except (ValueError, OSError, RuntimeError) as error:
            job.update(status="preparation_failed", reason=str(error), error_type=type(error).__name__)
            write_json(box / "preparation_attempt.json", job)
        manifest["jobs"].append(job)
        write_json(output / "condition.json", manifest)
    manifest["status"] = ("prepared" if candidates and all(j["status"] == "prepared" for j in manifest["jobs"])
                          else selection["status"] if method != "oracle" and not candidates else "preparation_failed")
    write_json(output / "condition.json", manifest)
    return manifest


def run_condition(root, prepared, output, executable, config_path, freeze_path, source_prepared, *, resume=False):
    from .blind_docking import run_condition as execute

    verify_freeze(root, config_path, freeze_path, source_prepared)
    lock = _load(Path(root) / LOCKS[2])
    if sha256(executable) not in {asset["sha256"] for asset in lock["assets"].values()}:
        raise ValueError("Vina binary differs from the frozen official release")
    return execute(prepared, output, executable, config_path, resume=resume)
