"""Conservative Gal9C matched-epimer preparation and descriptive score metrics.

Source graphs define chemistry; independently generated input coordinates never
become pose ground truth. NB at 2 mM is an operational class, not a numeric KD.
"""
from __future__ import annotations

from collections import Counter
import importlib.metadata
import json
import math
from pathlib import Path
import re
import shutil

import numpy as np

from .provenance import sha256, write_json
from .statistics import quantitative_specificity_by_context, specificity_by_context


EMBEDDING_SEED = 20261024


def _canonical(molecule, stereo=True):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=stereo)


def _source_graph(row):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    smiles = row.get("isomeric_smiles", row.get("curated_isomeric_smiles"))
    molecule = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("specificity ligand must have one valid connected source SMILES graph")
    molecule = Chem.RemoveHs(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True, useLegacyImplementation=False)
    if len(centers) != 5 or any(value == "?" for _, value in centers):
        raise ValueError("specificity ligand must have the five fixed, assigned source stereocenters")
    if molecule.GetNumHeavyAtoms() != row["expected_heavy_atoms"]:
        raise ValueError("specificity source heavy-atom count differs from the panel")
    if Chem.GetFormalCharge(molecule) != row["formal_charge"]:
        raise ValueError("specificity source formal charge differs from the panel")
    if row.get("molecular_formula") != rdMolDescriptors.CalcMolFormula(molecule):
        raise ValueError("specificity source formula differs from the panel")
    if row.get("curated_isomeric_smiles") and _canonical(molecule) != _canonical(
            Chem.MolFromSmiles(row["curated_isomeric_smiles"])):
        raise ValueError("specificity source stereochemical aliases disagree")
    return molecule


def _load_panel(panel_path):
    from rdkit import Chem

    document = json.loads(Path(panel_path).read_text())
    panel = document.get("primary_specificity_panel", document)
    rows = [dict(row) for row in panel["rows"] if row.get("primary_cohort") is True]
    frozen = panel["frozen_analysis"]
    if frozen.get("embedding_seed") != EMBEDDING_SEED or frozen.get("analysis_frozen_before_new_docking") is not True:
        raise ValueError("specificity chemistry/analysis must be frozen before docking")
    ids = [row["ligand_id"] for row in rows]
    if len(rows) != 12 or len(set(ids)) != 12 or set(ids) != set(frozen["ligand_ids"]):
        raise ValueError("specificity preparation requires all 12 conservative frozen ligands")
    if panel["target"]["uniprot_id"] != "O00182":
        raise ValueError("specificity panel requires the qualified Gal9C target")
    lookup, graphs = {}, {}
    for row in rows:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", row.get("case_id", "")):
            raise ValueError("specificity case_id must be a safe explicit directory name")
        if row.get("measured") is not True or type(row.get("binary_label")) is not int or row["binary_label"] not in (0, 1):
            raise ValueError("specificity classes require measured finite affinity or explicit conditional NB")
        if type(row.get("expected_heavy_atoms")) is not int or row["expected_heavy_atoms"] <= 0 or type(row.get("formal_charge")) is not int:
            raise ValueError("specificity source size and formal charge require explicit integers")
        if row["binary_label"] == 0:
            if row.get("kd_uM") is not None or row.get("numeric_kd_lower_bound_uM") is not None:
                raise ValueError("NB-at-2mM cannot be converted to an exact or bounded numeric KD")
            if row.get("source_measurement") != "NB_at_highest_tested_2_mM":
                raise ValueError("specificity negative must be the source conditional NB class")
        else:
            kd = row.get("kd_uM")
            if isinstance(kd, bool) or not isinstance(kd, (int, float)) or not math.isfinite(kd) or kd <= 0:
                raise ValueError("specificity positive requires a finite positive measured KD")
        graphs[row["ligand_id"]] = _source_graph(row)
        lookup[row["ligand_id"]] = row
    if len({row["case_id"] for row in rows}) != len(rows):
        raise ValueError("specificity case IDs must be unique")
    pairs = [dict(pair) for pair in panel["matched_pairs"] if pair.get("primary_cohort") is True]
    used = []
    for pair in pairs:
        left, right = pair["binder_id"], pair["nonbinder_id"]
        if left not in lookup or right not in lookup:
            raise ValueError("specificity pair references a ligand outside the frozen cohort")
        a, b = lookup[left], lookup[right]
        if a["binary_label"] != 1 or b["binary_label"] != 0:
            raise ValueError("specificity pair must contain one measured binder and one conditional NB")
        if a["pair_id"] != pair["pair_id"] or b["pair_id"] != pair["pair_id"]:
            raise ValueError("specificity ligand pair IDs disagree")
        if (a["expected_heavy_atoms"], a["formal_charge"]) != (b["expected_heavy_atoms"], b["formal_charge"]):
            raise ValueError("specificity pair is not exactly size and charge matched")
        mol_a, mol_b = graphs[left], graphs[right]
        if _canonical(mol_a, False) != _canonical(mol_b, False):
            raise ValueError("specificity pair stereo-stripped source graphs differ")
        centers_a = dict(Chem.FindMolChiralCenters(mol_a, useLegacyImplementation=False))
        centers_b = dict(Chem.FindMolChiralCenters(mol_b, useLegacyImplementation=False))
        mappings = mol_b.GetSubstructMatches(mol_a, useChirality=False, uniquify=False, maxMatches=1000)
        differences = [sum(value != centers_b.get(mapping[index]) for index, value in centers_a.items())
                       for mapping in mappings]
        if not differences or min(differences) != 1:
            raise ValueError("specificity pair must differ at exactly one assigned stereocenter")
        used.extend((left, right))
    if len(pairs) != 6 or len({pair["pair_id"] for pair in pairs}) != 6 or Counter(used) != Counter(ids):
        raise ValueError("specificity panel requires six disjoint epimer pairs covering all 12 ligands")
    return panel, rows, pairs, graphs


def _check_stereo(molecule, expected):
    from rdkit import Chem

    if _canonical(molecule) != expected:
        raise ValueError("specificity ligand graph/stereochemistry changed")
    derived = Chem.RemoveHs(Chem.Mol(molecule))
    xyz = derived.GetConformer().GetPositions()
    if not np.isfinite(xyz).all():
        raise ValueError("specificity input coordinates must be finite")
    Chem.AssignAtomChiralTagsFromStructure(derived, replaceExistingTags=True)
    Chem.AssignStereochemistry(derived, cleanIt=True, force=True)
    if _canonical(derived) != expected:
        raise ValueError("specificity generated coordinates do not preserve source stereochemistry")


def _prepare_ligand(graph, row, center, output):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    expected = _canonical(graph)
    molecule = Chem.AddHs(graph)
    names, atom_map = [], []
    heavy_number = hydrogen_number = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            heavy_number += 1
            name = f"L{heavy_number:03d}"
            names.append(name)
            atom_map.append({"name": name, "rdkit_index": atom.GetIdx(), "element": atom.GetSymbol(),
                             "formal_charge": atom.GetFormalCharge(),
                             "source_smiles_atom_index": atom.GetIdx()})
        else:
            hydrogen_number += 1
            name = f"H{hydrogen_number:03d}"
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name)
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        atom.SetMonomerInfo(info)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed, parameters.numThreads, parameters.enforceChirality = EMBEDDING_SEED, 1, True
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise ValueError("independent specificity ETKDG embedding failed; no native-coordinate fallback")
    _check_stereo(molecule, expected)
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError("unsupported specificity MMFF chemistry")
    if AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000) != 0:
        raise ValueError("independent specificity MMFF optimization did not converge")
    _check_stereo(molecule, expected)
    positions = molecule.GetConformer().GetPositions()
    heavy_indices = [row["rdkit_index"] for row in atom_map]
    positions += np.asarray(center, float) - positions[heavy_indices].mean(axis=0)
    for index, point in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(index, point)
    _check_stereo(molecule, expected)
    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1:
        raise ValueError("specificity ligand produced multiple Meeko setups")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(error)
    records = {}
    for line in pdbqt.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[77:].strip() in {"H", "HD", "HS"}:
            continue
        name = line[12:16].strip()
        if name in records:
            raise ValueError("specificity PDBQT contains duplicate heavy atom names")
        records[name] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    if set(records) != set(names) or len(records) != row["expected_heavy_atoms"]:
        raise ValueError("specificity Meeko preparation lost heavy atoms")
    smiles_rows = [line[len("REMARK SMILES "):].strip() for line in pdbqt.splitlines()
                   if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX")]
    if len(smiles_rows) != 1 or _canonical(Chem.MolFromSmiles(smiles_rows[0])) != expected:
        raise ValueError("specificity PDBQT source graph/stereochemistry differs")
    rounded = Chem.Mol(molecule)
    for mapping in atom_map:
        rounded.GetConformer().SetAtomPosition(mapping["rdkit_index"], records[mapping["name"]])
    _check_stereo(rounded, expected)
    max_rounding_error = max(float(np.linalg.norm(np.array(records[mapping["name"]]) -
                                                 positions[mapping["rdkit_index"]])) for mapping in atom_map)
    if max_rounding_error > 0.0009 or Chem.GetFormalCharge(molecule) != row["formal_charge"]:
        raise ValueError("specificity coordinate/charge integrity check failed")
    (output / "ligand.pdbqt").write_text(pdbqt)
    writer = Chem.SDWriter(str(output / "ligand.sdf"))
    writer.write(molecule)
    writer.close()
    loaded = list(Chem.SDMolSupplier(str(output / "ligand.sdf"), removeHs=False))
    if len(loaded) != 1 or loaded[0] is None:
        raise ValueError("specificity generated SDF did not round-trip")
    _check_stereo(loaded[0], expected)
    chemistry = {"source_isomeric_smiles": row["isomeric_smiles"], "isomeric_smiles": expected,
                 "heavy_atoms": len(names), "formal_charge": Chem.GetFormalCharge(molecule),
                 "defined_stereocenters": 5, "source_graph": "hand-curated published fixed stereochemical SMILES",
                 "initial_conformer": "independent ETKDGv3 then MMFF94, maximum 1000 iterations; convergence required",
                 "embedding_seed": EMBEDDING_SEED, "embedding_threads": 1,
                 "native_coordinates_used_for_embedding": False, "native_pose_reference_available": False,
                 "stereochemistry_verified_from_3d": True, "pdbqt_max_rounding_error_A": max_rounding_error,
                 "hydrogens": "RDKit valence hydrogens; Meeko merges nonpolar hydrogens",
                 "protonation": "source formal charges retained; no outcome-fitted protonation",
                 "versions": {name: importlib.metadata.version(name) for name in ("rdkit", "meeko", "numpy")}}
    return names, positions[heavy_indices].tolist(), atom_map, chemistry


def _box(metadata):
    center, size = np.asarray(metadata["box_center"], float), np.asarray(metadata["box_size"], float)
    if center.shape != (3,) or size.shape != (3,) or not np.isfinite(center).all() or not np.isfinite(size).all() or np.any(size <= 0):
        raise ValueError("specificity receptor must supply one finite fixed positive-size box")
    return center.tolist(), size.tolist()


def prepare_specificity_receptor(receptor_pdb, output_dir, *, box_center, box_size, provenance):
    """Prepare a caller-selected Gal9C receptor with the audited rigid helper.

    The caller freezes source/chain/site choices. Unsupported atoms, missing
    heavy atoms and cofactor preparation failures propagate; no repair or silent
    removal is attempted here. Galectins do not require calcium substitution.
    """
    from .aav import prepare_receptor

    source, output = Path(receptor_pdb).resolve(), Path(output_dir).resolve()
    if provenance.get("target_id") != "O00182" or provenance.get("source_pdb_id") != "3NV3":
        raise ValueError("specificity receptor requires the reviewed 3NV3 Gal9C provenance")
    center, size = _box({"box_center": box_center, "box_size": box_size})
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("specificity receptor output must be fresh")
    digest = sha256(source)
    log = prepare_receptor(source.read_text(), output / "holo_complete")
    if sha256(source) != digest:
        raise ValueError("specificity receptor input changed during preparation")
    result = {"schema_version": 1, "status": "prepared", "case_id": "3NV3", "family_id": "galectin_9",
              "target_id": "O00182", "source_provenance": provenance, "source_receptor_sha256": digest,
              "box_center": center, "box_size": size,
              "receptors": [{"id": "holo_complete", "state": "holo", "assembly_condition": "isolated_C_terminal_CRD",
                              "pdb_file": "holo_complete.pdb", "pdbqt_file": "holo_complete.pdbqt",
                              "atom_changes_file": "holo_complete_atom_changes.json"}],
              "information_regime": "known Gal9C site from independent prior glycan structure; shared by all alternatives",
              "receptor_preparation": log, "heavy_atom_completion": False, "heavy_atom_minimization": False}
    write_json(output / "preparation.json", result)
    return result


def prepare_specificity_panel(panel_path, receptor_prepared_dir, output_root):
    """Write all 12 conservative ligand cases with one frozen receptor and box."""
    panel_path, receptor_dir, output_root = map(Path, (panel_path, receptor_prepared_dir, output_root))
    panel, rows, _, graphs = _load_panel(panel_path)
    receptor_manifest = receptor_dir / "preparation.json"
    metadata = json.loads(receptor_manifest.read_text())
    if metadata.get("case_id") != "3NV3" or metadata.get("target_id") != panel["target"]["uniprot_id"] or len(metadata["receptors"]) != 1:
        raise ValueError("specificity panel requires one reviewed prepared 3NV3 Gal9C receptor")
    center, size = _box(metadata)
    files = {}
    for receptor in metadata["receptors"]:
        for key in ("pdb_file", "pdbqt_file", "atom_changes_file"):
            if key not in receptor:
                continue
            relative = Path(receptor[key])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("specificity receptor files must stay within the prepared directory")
            files[str(relative)] = sha256(receptor_dir / relative)
    panel_hash, receptor_hash = sha256(panel_path), sha256(receptor_manifest)
    for row in rows:
        output = output_root / row["case_id"]
        if output.exists() and any(output.iterdir()):
            raise ValueError("specificity case directories must be fresh; stale preparations cannot be reused")
    prepared = []
    for row in rows:
        output = output_root / row["case_id"]
        output.mkdir(parents=True, exist_ok=True)
        names, coordinates, atom_map, chemistry = _prepare_ligand(graphs[row["ligand_id"]], row, center, output)
        for relative, expected in files.items():
            if sha256(receptor_dir / relative) != expected:
                raise ValueError("specificity source receptor changed between alternatives")
            destination = output / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(receptor_dir / relative, destination)
            if sha256(destination) != expected:
                raise ValueError("specificity receptor copy differs from its frozen input")
        write_json(output / "source_graph.json", {"panel_id": panel["panel_id"], "source": panel["source"],
                                                 "source_row": row, "atom_map": atom_map,
                                                 "panel_input_sha256": panel_hash})
        write_json(output / "receptor_preparation.json", metadata)
        write_json(output / "reference.json", {"names": names, "coordinates": coordinates, "atom_map": atom_map,
                                                "pose_reference_eligible": False,
                                                "coordinate_source": "independently generated docking input only",
                                                "scope": "atom-identity validation only; input RMSD is not a pose endpoint"})
        result = {"schema_version": 1, "status": "prepared", "case_id": row["case_id"], "family_id": "galectin_9",
                  "panel_id": panel["panel_id"], "ligand_id": row["ligand_id"], "pair_id": row["pair_id"],
                  "cohort_role": "exploratory_assay_development", "pose_reference_eligible": False,
                  "information_regime": "known_site_independent_prior_glycan_structure_shared_by_all_alternatives",
                  "box_center": center, "box_size": size, "receptors": metadata["receptors"],
                  "ligand_file": "ligand.pdbqt", "reference_file": "reference.json", "chemistry": chemistry,
                  "assay_endpoint": "finite measured KD versus conditional NB at source-tested 2mM",
                  "assay_row": row, "frozen_analysis": panel["frozen_analysis"],
                  "input_fingerprint": {"panel_sha256": panel_hash, "receptor_preparation_sha256": receptor_hash,
                                        "receptor_files_sha256": files},
                  "ligand_sha256": sha256(output / "ligand.pdbqt"), "ligand_sdf_sha256": sha256(output / "ligand.sdf"),
                  "limitations": ["one target domain and synthetic glycoside series; exploratory development",
                                  "conditional NB does not establish a numeric KD bound or universal nonbinding",
                                  "single independently generated ring/conformer input; no source atomic pose truth",
                                  panel["target"].get("construct_caveat", "assay construct approximation requires review")]}
        write_json(output / "preparation.json", result)
        prepared.append(result)
    if sha256(panel_path) != panel_hash or sha256(receptor_manifest) != receptor_hash:
        raise ValueError("specificity source manifest changed during preparation")
    if any(sha256(receptor_dir / relative) != expected for relative, expected in files.items()):
        raise ValueError("specificity source receptor changed during preparation")
    return prepared


def summarize_specificity_scores(panel_path, scores_by_ligand):
    """Describe all 12 preaggregated Vina scores; lower score predicts binding.

    The caller supplies one frozen seed-mean top-1 score per ligand and audits
    its complete receptor-by-seed execution denominator before calling this.
    No classifier, threshold, docking parameter, or protonation choice is fit.
    """
    panel, rows, pairs, _ = _load_panel(panel_path)
    if set(scores_by_ligand) != {row["ligand_id"] for row in rows}:
        raise ValueError("specificity score denominator must contain all 12 frozen ligands exactly")
    for score in scores_by_ligand.values():
        if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
            raise ValueError("specificity requires a finite score for every frozen ligand")
    context, target = panel["panel_id"], panel["target"]["uniprot_id"]
    measured = [{"context_id": context, "target_id": target, "ligand_id": row["ligand_id"],
                 "label": "binder" if row["binary_label"] else "nonbinder",
                 "score": -scores_by_ligand[row["ligand_id"]], "heavy_atoms": row["expected_heavy_atoms"],
                 "formal_charge": row["formal_charge"]} for row in rows]
    matching = [{"context_id": context, "binder_id": pair["binder_id"], "nonbinder_id": pair["nonbinder_id"]}
                for pair in pairs]
    classification = specificity_by_context(measured, matched_pairs=matching)
    baselines = {name: specificity_by_context([{**row, "score": row[field]} for row in measured], matched_pairs=matching)
                 for name, field in (("heavy_atom_count", "heavy_atoms"), ("formal_charge", "formal_charge"))}
    quantified = [{"context_id": context, "target_id": target, "ligand_id": row["ligand_id"],
                   "score": -scores_by_ligand[row["ligand_id"]], "kd": row["kd_uM"], "kd_unit": "uM"}
                  for row in rows if row["binary_label"] == 1]
    pair_results = []
    for pair in pairs:
        binder, negative = scores_by_ligand[pair["binder_id"]], scores_by_ligand[pair["nonbinder_id"]]
        pair_results.append({"pair_id": pair["pair_id"], "binder_id": pair["binder_id"],
                             "conditional_nonbinder_id": pair["nonbinder_id"],
                             "binder_score_kcal_mol": binder, "conditional_nonbinder_score_kcal_mol": negative,
                             "binder_minus_conditional_nonbinder_kcal_mol": binder - negative,
                             "concordance": float(binder < negative) + 0.5 * (binder == negative)})
    return {"panel_id": context, "classification": classification, "baselines": baselines,
            "per_pair": pair_results, "score_direction": "lower Vina seed-mean top1 score predicts binding",
            "tie_policy": "exact equal aggregate scores receive half concordance",
            "positive_only_quantitative_association": quantitative_specificity_by_context(quantified),
            "right_censored_quantitative_contrast": {"status": "ineligible_no_justified_numeric_KD_bounds",
                                                      "estimate": None,
                                                      "reason": "NB at tested 2mM is an operational class, not a measured KD lower bound"},
            "interpretation": "descriptive conditional same-assay specificity; one target/scaffold; no held-out generalization"}
