"""AAV9 oracle-site diagnostic preparation with audited atom identities.

This module uses the deposited biological assembly once. Its neighborhood and
20 A box depend on the reference ligand, so its outputs are never blind inputs.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import string

import gemmi
import numpy as np

PROTEIN_RESIDUES = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def _dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def atom_key(chain, residue, atom):
    return (chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, atom.name)


def protein_atoms(model):
    """Return unambiguous protein heavy atoms; insertion codes are identities."""
    result = {}
    for chain in model:
        for residue in chain:
            if residue.het_flag != "A":
                continue
            if residue.name not in PROTEIN_RESIDUES:
                raise ValueError(f"unsupported protein residue {chain.name}:{residue.seqid}:{residue.name}")
            for atom in residue:
                if atom.element.is_hydrogen:
                    continue
                if atom.altloc not in ("\x00", " "):
                    raise ValueError("alternate conformers require explicit curation before preparation")
                key = atom_key(chain, residue, atom)
                if key in result:
                    raise ValueError(f"duplicate protein atom identity {key}")
                result[key] = np.array(list(atom.pos), dtype=float)
    if not result:
        raise ValueError("empty protein model")
    return result


def materialize_assembly(path):
    """Materialize assembly 1 exactly once and verify the expected AAV 60-mer.

    Current 7MTZ/7MT0 deposit all 60 subunits with an identity operator. We use
    the declaration instead of applying an additional icosahedral expansion.
    Unexpected operators or chain multiplicity require a new curation review.
    """
    try:
        structure = gemmi.read_structure(str(path))
    except (RuntimeError, IndexError) as error:
        raise ValueError(f"invalid or empty structure source: {path}") from error
    if len(structure) != 1:
        raise ValueError("AAV diagnostic requires one deposited model")
    assemblies = [a for a in structure.assemblies if a.name == "1"]
    if len(assemblies) != 1:
        raise ValueError("missing or ambiguous biological assembly 1")
    assembly = assemblies[0]
    if len(assembly.generators) != 1 or len(assembly.generators[0].operators) != 1:
        raise ValueError("assembly definition changed: explicit independent review required")
    transform = assembly.generators[0].operators[0].transform
    if not np.allclose(np.array(transform.mat), np.eye(3)) or not np.allclose(list(transform.vec), 0):
        raise ValueError("nonidentity AAV assembly requires curated author-chain mapping")
    model = gemmi.make_assembly(assembly, structure[0], gemmi.HowToNameCopiedChain.Dup)
    chain_names = [c.name for c in model if any(r.het_flag == "A" for r in c)]
    if len(chain_names) != 60 or len(set(chain_names)) != 60:
        raise ValueError("AAV assembly must contain exactly 60 unique protein subunits")
    if any(len(name) != 1 or name not in string.ascii_letters + string.digits for name in chain_names):
        raise ValueError("legacy PDB needs a reviewed unique one-character chain map")
    return model, {"assembly_id": "1", "deposited_model": 1, "materializations": 1,
                   "operator": "identity", "protein_subunits": 60,
                   "author_chain_ids_preserved": True}


def fit(source, target):
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3 or len(source) < 3:
        raise ValueError("alignment needs >=3 matched receptor coordinates")
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError("alignment coordinates must be finite")
    if np.linalg.matrix_rank(source - source.mean(0)) < 2:
        raise ValueError("receptor alignment is geometrically underdetermined")
    u, _, vt = np.linalg.svd((source - source.mean(0)).T @ (target - target.mean(0)))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return rotation, target.mean(0) - source.mean(0) @ rotation


def box_distance(coordinates, center, box_size):
    """Euclidean distance from each atom to the closed docking box."""
    xyz = np.asarray(coordinates, dtype=float)
    center, box_size = np.asarray(center, float), np.asarray(box_size, float)
    if center.shape != (3,) or box_size.shape != (3,) or np.any(box_size <= 0):
        raise ValueError("box needs a finite center and positive dimensions")
    if not np.isfinite(xyz).all() or not np.isfinite(center).all() or not np.isfinite(box_size).all():
        raise ValueError("box coordinates must be finite")
    return np.linalg.norm(np.maximum(np.abs(xyz - center) - box_size / 2, 0), axis=-1)


def select_neighborhood(atoms, center, box_size, cutoff=6.0):
    """Retain whole subunits intersecting the box plus interaction cutoff."""
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("interaction cutoff must be positive")
    selected = sorted({key[0] for key, coordinate in atoms.items()
                       if box_distance([coordinate], center, box_size)[0] <= cutoff})
    if not selected:
        raise ValueError("no subunits cover the docking box")
    omitted = [xyz for key, xyz in atoms.items() if key[0] not in selected]
    nearest = float(np.min(box_distance(omitted, center, box_size))) if omitted else None
    if nearest is not None and nearest <= cutoff:
        raise AssertionError("neighborhood omitted atoms inside the interaction region")
    return selected, nearest


def _protein_pdb(model, selected, rotation=None, translation=None):
    output = gemmi.Structure()
    output.add_model(gemmi.Model("1"))
    for original in model:
        if original.name not in selected:
            continue
        chain = original.clone()
        for i in reversed(range(len(chain))):
            if chain[i].het_flag != "A":
                del chain[i]
        for residue in chain:
            for i in reversed(range(len(residue))):
                if residue[i].element.is_hydrogen:
                    del residue[i]
            if rotation is not None:
                for atom in residue:
                    atom.pos = gemmi.Position(*(np.array(list(atom.pos)) @ rotation + translation))
        output[0].add_chain(chain)
    if not len(output[0]):
        raise ValueError("empty selected receptor")
    return output.make_pdb_string()


def _pdb_atoms(text, pdbqt=False):
    atoms = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(),
               line[17:20].strip(), line[12:16].strip())
        if key in atoms:
            raise ValueError(f"duplicate output atom identity {key}")
        element = line[77:].strip().split()[0] if pdbqt else line[76:78].strip()
        if not element:
            raise ValueError(f"missing atom element/type for {key}")
        is_h = element in {"H", "HD", "HS", "D"}
        atoms[key] = {"coordinates": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                      "element_or_type": element, "hydrogen": is_h}
    if not atoms:
        raise ValueError("empty receptor atom output")
    return atoms


def _protein_topology(pdb_text):
    """Build residue topology from chemical atom names, never surface proximity.

    Meeko 0.8's distance heuristic can join close nonbonded capsid atoms. All
    inferred edges are audited; peptide and disulfide connectivity is explicit.
    """
    from meeko import Polymer
    from meeko.polymer import find_inter_mols_bonds

    parsed = Polymer._pdb_to_residue_mols(pdb_text)
    raw = {}
    for key, (mol, resname, missed, needed) in parsed.items():
        if mol is None or missed or needed or resname not in PROTEIN_RESIDUES:
            raise ValueError(f"uncurated receptor residue {key}:{resname}")
        raw[key] = (mol, resname)
    inferred = find_inter_mols_bonds(raw)
    bonds, rejected, accepted = {}, [], []
    sg_partners = Counter()
    def identity(resid, index):
        atom = raw[resid][0].GetAtomWithIdx(index)
        return atom.GetPDBResidueInfo().GetName().strip()
    def position(resid, index):
        return np.array(list(raw[resid][0].GetConformer().GetAtomPosition(index)))
    for pair, edges in inferred.items():
        left, right = pair
        lc, ln = left.split(":")
        rc, rn = right.split(":")
        # Insertion codes need a separately curated sequence ordering.
        consecutive = ln.lstrip("-").isdigit() and rn.lstrip("-").isdigit() and abs(int(ln)-int(rn)) == 1 and lc == rc
        for li, ri in edges:
            lname, rname = identity(left, li), identity(right, ri)
            distance = float(np.linalg.norm(position(left, li)-position(right, ri)))
            peptide = consecutive and ((int(ln) < int(rn) and lname == "C" and rname == "N") or
                                       (int(rn) < int(ln) and lname == "N" and rname == "C"))
            disulfide = raw[left][1] == raw[right][1] == "CYS" and lname == rname == "SG" and 1.8 <= distance <= 2.25
            detail = {"residues": [left, right], "atoms": [lname, rname], "distance_A": distance}
            if peptide:
                if not 1.1 <= distance <= 1.7:
                    raise ValueError(f"implausible peptide bond {detail}")
                bonds.setdefault(pair, []).append((li, ri))
                accepted.append({**detail, "kind": "peptide"})
            elif disulfide:
                sg_partners.update([left, right])
                bonds.setdefault(pair, []).append((li, ri))
                accepted.append({**detail, "kind": "disulfide"})
            else:
                rejected.append(detail)
    if any(count > 1 for count in sg_partners.values()):
        raise ValueError("ambiguous disulfide pairing")
    # A missing peptide link within a continuously numbered chain is an error,
    # not an instruction to create a covalent bond through an unresolved gap.
    by_chain = {}
    for resid, (mol, _) in raw.items():
        ch, number = resid.split(":")
        if not number.lstrip("-").isdigit():
            raise ValueError("insertion-coded peptide topology needs a curated sequence map")
        by_chain.setdefault(ch, []).append((int(number), resid))
    accepted_peptides = {frozenset(v["residues"]) for v in accepted if v["kind"] == "peptide"}
    gaps = []
    for chain, entries in by_chain.items():
        entries.sort()
        for (ln, left), (rn, right) in zip(entries, entries[1:]):
            if rn == ln + 1 and frozenset((left, right)) not in accepted_peptides:
                raise ValueError(f"missing peptide connection {left}--{right}")
            if rn != ln + 1:
                gaps.append([left, right])
    return raw, bonds, {"policy": "same-chain consecutive C-N peptide bonds 1.1-1.7 A; CYS SG-SG 1.8-2.25 A",
                        "accepted_inter_residue_bonds": accepted,
                        "rejected_distance_inferred_nonbonded_edges": rejected,
                        "unresolved_sequence_gaps": gaps}


def _restore_terminal_oxygen_names(text, before, pdbqt=False):
    """Audit Meeko's graph-equivalent terminal O/OXT naming permutation."""
    after = _pdb_atoms(text, pdbqt=pdbqt)
    replacements, permutations = {}, []
    for key in before:
        if key[-1] != "OXT":
            continue
        other = key[:-1] + ("O",)
        if other not in before or key not in after or other not in after:
            continue
        displacement = [np.linalg.norm(np.array(before[k]["coordinates"])-after[k]["coordinates"]) for k in (key, other)]
        swapped = [np.linalg.norm(np.array(before[k]["coordinates"])-after[j]["coordinates"])
                   for k, j in ((key, other), (other, key))]
        if max(displacement) > 0.002 and max(swapped) <= 0.002:
            replacements[key], replacements[other] = "O", "OXT"
            permutations.append({"residue": list(key[:-1]), "permutation": {"O": "OXT", "OXT": "O"},
                                 "reason": "terminal carboxylate graph-equivalent oxygens; exact coordinate match"})
    lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip(), line[12:16].strip())
            if key in replacements:
                line = line[:12] + replacements[key].rjust(4) + line[16:]
        lines.append(line)
    return "".join(lines), permutations


def prepare_receptor(pdb_text, output_stem):
    """Prepare a rigid receptor with Meeko and audit every atom change.

    No heavy-atom completion, residue deletion or coordinate minimization is
    permitted. Meeko's frozen standard residue templates assign protonation;
    all histidines use neutral HIE. Nonpolar H are merged for PDBQT as usual.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates

    output_stem = Path(output_stem)
    before = _pdb_atoms(pdb_text)
    histidines = {f"{key[0]}:{key[1]}{key[2]}": "HIE" for key in before if key[3] == "HIS"}
    raw, bonds, topology_log = _protein_topology(pdb_text)
    polymer = Polymer(raw, bonds, ResidueChemTemplates.create_from_defaults(),
                      mk_prep=MoleculePreparation(), set_template=histidines)
    ignored = polymer.get_ignored_monomers()
    if ignored:
        raise ValueError(f"Meeko ignored receptor residues: {list(ignored)}")
    prepared_pdb, terminal_permutations = _restore_terminal_oxygen_names(polymer.to_pdb(), before)
    after = _pdb_atoms(prepared_pdb)
    rigid_pdbqt, flexible = PDBQTWriterLegacy.write_from_polymer(polymer)
    if flexible:
        raise ValueError("rigid receptor unexpectedly has flexible components")
    rigid_pdbqt, pdbqt_permutations = _restore_terminal_oxygen_names(rigid_pdbqt, before, pdbqt=True)
    dock_atoms = _pdb_atoms(rigid_pdbqt, pdbqt=True)
    before_heavy = {key for key, value in before.items() if not value["hydrogen"]}
    after_heavy = {key for key, value in after.items() if not value["hydrogen"]}
    dock_heavy = {key for key, value in dock_atoms.items() if not value["hydrogen"]}
    if before_heavy != after_heavy or before_heavy != dock_heavy:
        raise ValueError(f"receptor heavy atom change: added={sorted(after_heavy-before_heavy)}; "
                         f"deleted={sorted(before_heavy-after_heavy)}; pdbqt_missing={sorted(before_heavy-dock_heavy)}")
    displacement = max(float(np.linalg.norm(np.array(before[key]["coordinates"]) - after[key]["coordinates"]))
                       for key in before_heavy)
    if displacement > 0.002:
        raise ValueError("Meeko moved heavy atoms")
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    log = {"status": "ok", "preparation": "Meeko Polymer; Gasteiger charges; rigid receptor",
           "meeko_version": importlib.metadata.version("meeko"),
           "protonation_policy": "Meeko standard amino-acid templates; all histidines neutral HIE; no pKa prediction",
           "heavy_atom_completion": False, "coordinate_minimization": False,
           "input_atoms": len(before), "prepared_atoms": len(after), "pdbqt_atoms": len(dock_atoms),
           "input_heavy_atoms": len(before_heavy), "prepared_heavy_atoms": len(after_heavy),
           "pdbqt_heavy_atoms": len(dock_heavy), "max_heavy_displacement_A": displacement,
           "added_atoms": [list(key) for key in added], "removed_atoms": [list(key) for key in removed],
           "pdbqt_merged_atoms": [list(key) for key in sorted(set(after)-set(dock_atoms))],
           "templates": {key: value.residue_template_key for key, value in polymer.get_valid_monomers().items()},
           "meeko_log": polymer.log, "topology": topology_log,
           "terminal_oxygen_name_permutations": {"pdb": terminal_permutations, "pdbqt": pdbqt_permutations}}
    output_stem.with_suffix(".pdb").write_text(pdb_text)
    output_stem.with_suffix(".pdbqt").write_text(rigid_pdbqt)
    _dump(output_stem.with_name(output_stem.name + "_atom_changes.json"), log)
    return log


def validate_native_stereochemistry(mol, reference_names, reference_coordinates):
    """Validate native coordinates against the CCD-defined stereochemical graph.

    Atom labels supply the graph mapping. Native coordinates only validate the
    reference; independently embedded docking conformers remain unchanged.
    """
    from rdkit import Chem

    reference_names = list(reference_names)
    native_xyz = np.asarray(reference_coordinates, dtype=float)
    if len(set(reference_names)) != len(reference_names):
        raise ValueError("duplicate native ligand atom names")
    if native_xyz.shape != (len(reference_names), 3) or not np.isfinite(native_xyz).all():
        raise ValueError("invalid native ligand coordinate array")
    heavy = Chem.RemoveHs(mol)
    graph_names = [atom.GetProp("ccd_name") for atom in heavy.GetAtoms()]
    if len(set(graph_names)) != len(graph_names) or set(graph_names) != set(reference_names):
        raise ValueError("native ligand and CCD graph atom identities differ")
    expected_smiles = Chem.MolToSmiles(heavy, isomericSmiles=True)
    native_by_name = dict(zip(reference_names, native_xyz))
    conformer = Chem.Conformer(heavy.GetNumAtoms())
    for atom in heavy.GetAtoms():
        conformer.SetAtomPosition(atom.GetIdx(), native_by_name[atom.GetProp("ccd_name")])
    heavy.RemoveAllConformers()
    heavy.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(heavy, replaceExistingTags=True)
    Chem.AssignStereochemistry(heavy, cleanIt=True, force=True)
    native_smiles = Chem.MolToSmiles(heavy, isomericSmiles=True)
    if native_smiles != expected_smiles:
        raise ValueError(f"native stereochemistry differs from CCD graph: {native_smiles} != {expected_smiles}")
    return {"status": "passed", "heavy_atoms": len(graph_names),
            "ccd_isomeric_smiles": expected_smiles, "native_isomeric_smiles": native_smiles,
            "method": "RDKit stereochemistry assigned from mapped native heavy coordinates and compared with CCD ideal graph",
            "density_validation": "not performed; matching stereochemistry alone does not establish local map support"}


def prepare_ligand(ccd_path, reference_names, center, output_dir, reference_coordinates=None):
    """Generate GAL independently from the CCD graph and ideal stereochemistry."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    block = gemmi.cif.read_file(str(ccd_path)).sole_block()
    editable = Chem.RWMol()
    indices = {}
    for row in block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]):
        if row[0] in indices:
            raise ValueError("duplicate CCD atom name")
        atom = Chem.Atom(row[1])
        atom.SetFormalCharge(int(row[2]))
        atom.SetProp("ccd_name", row[0])
        indices[row[0]] = editable.AddAtom(atom)
    for row in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        if row[2] not in {"SING", "DOUB", "TRIP"}:
            raise ValueError(f"unsupported CCD bond order {row[2]}")
        editable.AddBond(indices[row[0]], indices[row[1]],
                         {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}[row[2]])
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for row in block.find("_chem_comp_atom.", ["atom_id", "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal"]):
        conformer.SetAtomPosition(indices[row[0]], tuple(float(x) for x in list(row)[1:]))
    mol.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    heavy_names = [atom.GetProp("ccd_name") for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if set(heavy_names) != set(reference_names) or len(heavy_names) != len(reference_names):
        raise ValueError("native GAL heavy atoms differ from explicit CCD graph")
    if reference_coordinates is not None:
        validation = validate_native_stereochemistry(mol, reference_names, reference_coordinates)
        _dump(Path(output_dir) / "native_chemistry_validation.json", validation)
    mol.RemoveAllConformers()
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20261024
    if AllChem.EmbedMolecule(mol, parameters) != 0:
        raise ValueError("independent ligand conformer embedding failed")
    if AllChem.MMFFOptimizeMolecule(mol) != 0:
        raise ValueError("independent ligand MMFF optimization did not converge")
    heavy_positions = np.array([list(mol.GetConformer().GetAtomPosition(atom.GetIdx()))
                                for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1])
    shift = np.asarray(center) - heavy_positions.mean(0)
    for atom in mol.GetAtoms():
        position = np.array(list(mol.GetConformer().GetAtomPosition(atom.GetIdx()))) + shift
        mol.GetConformer().SetAtomPosition(atom.GetIdx(), position)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(atom.GetProp("ccd_name").rjust(4))
        info.SetResidueName("GAL")
        info.SetResidueNumber(1)
        info.SetChainId("L")
        atom.SetMonomerInfo(info)
    setup = MoleculePreparation().prepare(mol)[0]
    pdbqt, success, message = PDBQTWriterLegacy.write_string(setup)
    if not success:
        raise ValueError(message)
    output_names = [line[12:16].strip() for line in pdbqt.splitlines()
                    if line.startswith(("ATOM", "HETATM")) and line[77:].strip().split()[0] not in {"H", "HD", "HS"}]
    if Counter(output_names) != Counter(reference_names):
        raise ValueError("ligand heavy atom identity lost in PDBQT")
    (Path(output_dir) / "galactose.pdbqt").write_text(pdbqt)
    writer = Chem.SDWriter(str(Path(output_dir) / "galactose.sdf"))
    writer.write(mol)
    writer.close()
    return {"isomeric_smiles": smiles, "heavy_atoms": len(heavy_names), "formal_charge": Chem.GetFormalCharge(mol),
            "graph_source": "CCD GAL explicit atom/bond tables", "stereochemistry_source": "CCD ideal coordinates",
            "initial_conformer": "independent ETKDGv3 then MMFF94; centered in oracle box",
            "embedding_seed": 20261024, "ccd_atom_names": heavy_names,
            "glycosidic_linkages": 0, "native_coordinates_used_for_embedding": False}


def prepare(raw_dir, output_dir):
    """Prepare two receptor states and an explicitly incomplete monomer control."""
    raw_dir, output_dir = Path(raw_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    holo, holo_assembly = materialize_assembly(raw_dir / "7MTZ.cif")
    apo, apo_assembly = materialize_assembly(raw_dir / "7MT0.cif")
    ligands = [(chain, residue) for chain in holo for residue in chain
               if chain.name == "A" and residue.name == "GAL" and residue.seqid.num == 801]
    if len(ligands) != 1:
        raise ValueError("expected exactly reference author A GAL 801")
    chain, ligand = ligands[0]
    native_atoms = [atom for atom in ligand if not atom.element.is_hydrogen]
    names = [atom.name for atom in native_atoms]
    if len(names) != len(set(names)) or any(atom.altloc not in ("\x00", " ") for atom in native_atoms):
        raise ValueError("ambiguous reference ligand atom identity")
    native_xyz = np.array([list(atom.pos) for atom in native_atoms])
    center, box_size = native_xyz.mean(0), np.array([20.0, 20.0, 20.0])
    holo_atoms, apo_atoms = protein_atoms(holo), protein_atoms(apo)
    # Residue identities are checked globally rather than inferring a correspondence from ligand geometry.
    holo_residues = {(k[0], k[1], k[2]): k[3] for k in holo_atoms}
    apo_residues = {(k[0], k[1], k[2]): k[3] for k in apo_atoms}
    substitutions = [key for key in holo_residues.keys() & apo_residues.keys()
                     if holo_residues[key] != apo_residues[key]]
    if substitutions:
        raise ValueError(f"apo/holo sequence identity mismatch: {substitutions}")
    ca_keys = sorted(key for key in holo_atoms.keys() & apo_atoms.keys() if key[-1] == "CA")
    rotation, translation = fit([apo_atoms[key] for key in ca_keys], [holo_atoms[key] for key in ca_keys])
    aligned_apo_atoms = {key: xyz @ rotation + translation for key, xyz in apo_atoms.items()}
    alignment_rmsd = float(np.sqrt(np.mean(np.sum((np.array([aligned_apo_atoms[k] for k in ca_keys]) -
                                                     np.array([holo_atoms[k] for k in ca_keys]))**2, axis=1))))
    selected_holo, omitted_holo = select_neighborhood(holo_atoms, center, box_size)
    selected_apo, omitted_apo = select_neighborhood(aligned_apo_atoms, center, box_size)
    selected = sorted(set(selected_holo) | set(selected_apo))
    if "A" not in selected:
        raise ValueError("reference protein subunit missing from complete neighborhood")
    # The complete receptor has the union of the independently checked apo/holo neighborhoods.
    receptor_keys = sorted(key for key in holo_atoms if key[0] in selected)
    receptor_xyz = np.array([holo_atoms[key] for key in receptor_keys])
    contacts = np.argwhere(np.linalg.norm(native_xyz[:, None] - receptor_xyz[None, :], axis=2) < 4.0)
    ligand_metadata = prepare_ligand(raw_dir / "GAL.cif", names, center, output_dir, reference_coordinates=native_xyz)
    receptors = []
    for state, model in [("holo", holo), ("apo", apo)]:
        for condition, chains in [("complete", selected), ("monomer", ["A"])]:
            receptor_id = f"{state}_{condition}"
            text = _protein_pdb(model, chains, rotation if state == "apo" else None,
                                translation if state == "apo" else None)
            log = prepare_receptor(text, output_dir / f"{receptor_id}_receptor")
            receptors.append({"id": receptor_id, "state": state, "assembly_condition": condition,
                              "pdb_file": f"{receptor_id}_receptor.pdb",
                              "pdbqt_file": f"{receptor_id}_receptor.pdbqt",
                              "atom_log_file": f"{receptor_id}_receptor_atom_changes.json",
                              "author_chains": chains, "heavy_atoms": log["prepared_heavy_atoms"],
                              "is_incomplete_control": condition == "monomer"})
    holo_keys, apo_keys = sorted(holo_atoms), sorted(aligned_apo_atoms)
    np.savez_compressed(output_dir / "full_assembly.npz",
                        holo_keys=np.array(holo_keys, dtype=str),
                        holo_coordinates=np.array([holo_atoms[k] for k in holo_keys]),
                        apo_keys=np.array(apo_keys, dtype=str),
                        apo_coordinates=np.array([aligned_apo_atoms[k] for k in apo_keys]))
    reference = {"names": names, "coordinates": native_xyz.tolist(),
                 "reference_ligand": {"pdb_id": "7MTZ", "model": 1, "author_chain": "A", "author_residue": 801,
                                      "label_asym_id": ligand.subchain, "component_id": "GAL"},
                 "receptor_atom_keys": [list(k) for k in receptor_keys], "receptor_coordinates": receptor_xyz.tolist(),
                 "native_contact_pairs_4A": contacts.tolist(), "contact_cutoff_A": 4.0,
                 "full_assembly_file": "full_assembly.npz", "ligand_occupancies": [float(a.occ) for a in native_atoms],
                 "ligand_B_iso": [float(a.b_iso) for a in native_atoms],
                 "reference_quality": "development diagnostic; local map/density adjudication not completed"}
    _dump(output_dir / "reference.json", reference)
    metadata = {"status": "prepared", "information_regime": "true-site-local; reference-derived neighborhood; not blind",
                "ligand_file": "galactose.pdbqt", "reference_file": "reference.json", "receptors": receptors,
                "box_center": center.tolist(), "box_size": box_size.tolist(), "interaction_cutoff_A": 6.0,
                "neighborhood_rule": "union of apo/holo whole subunits within 6 A of any point in the 20 A box",
                "selected_author_chains": selected, "pdb_chain_mapping": {name: name for name in selected},
                "closest_omitted_holo_atom_to_box_A": omitted_holo, "closest_omitted_apo_atom_to_box_A": omitted_apo,
                "holo_assembly": holo_assembly, "apo_assembly": apo_assembly,
                "apo_to_holo_CA_rmsd_A": alignment_rmsd, "alignment_CA_count": len(ca_keys),
                "alignment_scope": "all corresponding whole-assembly CA atoms; author chain/residue identity",
                "rotation": rotation.tolist(), "translation": translation.tolist(),
                "missing_holo_residues_in_apo": [list(k) for k in sorted(holo_residues.keys()-apo_residues.keys())],
                "missing_apo_residues_in_holo": [list(k) for k in sorted(apo_residues.keys()-holo_residues.keys())],
                "ligand": ligand_metadata,
                "dry_baseline": True, "repaired_heavy_atoms": False,
                "input_sha256": {name: sha256((raw_dir / name).read_bytes()).hexdigest() for name in ["7MTZ.cif", "7MT0.cif", "GAL.cif"]},
                "versions": {name: importlib.metadata.version(name) for name in ["gemmi", "rdkit", "meeko", "numpy", "scipy"]}}
    _dump(output_dir / "preparation.json", metadata)
    return metadata
