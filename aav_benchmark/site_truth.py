"""Evaluation-only published-site truth, separate from site prediction inputs.

This module reads modeled native glycans and must be called only after the
ligand-free prediction inputs and candidate ranks have been frozen. It neither
generates ligands nor tests their full chemical/atomic-pose eligibility.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .aav import PROTEIN_RESIDUES
from .external import case_definition, pdb_heavy_atoms


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEDGER = PROJECT_ROOT / "benchmarks/site_reference_ledger.json"
CONTACT_CUTOFF_A = 4.0


def _read_pinned(path, expected, *, url, role):
    path = Path(path)
    payload = path.read_bytes()
    actual = sha256(payload).hexdigest()
    if actual != expected:
        raise ValueError(f"Site truth source checksum mismatch: {path}")
    return payload, {"path": str(path.resolve()), "sha256": actual,
                     "expected_sha256": expected, "url": url, "role": role}


def _contact_residues(tree, protein_coordinates, residue_keys, ligand_coordinates):
    """Query one cached protein tree; enforce strict <4 A after its inclusive query."""
    ligand_coordinates = np.asarray(ligand_coordinates, dtype=float)
    indices = {i for group in tree.query_ball_point(ligand_coordinates, CONTACT_CUTOFF_A)
               for i in group}
    if not indices:
        return []
    selected = np.asarray(sorted(indices), dtype=int)
    distances = np.linalg.norm(
        protein_coordinates[selected, None, :] - ligand_coordinates[None, :, :], axis=2)
    contacts = {tuple(residue_keys[i]) for i, distance in zip(selected, distances.min(axis=1))
                if distance < CONTACT_CUTOFF_A}
    return [list(key) for key in sorted(contacts)]


def _site(identifier, ligand_atoms, tree, protein_coordinates, residue_keys):
    keys = [list(atom["key"]) for atom in ligand_atoms]
    coordinates = [list(atom["coordinates"]) for atom in ligand_atoms]
    names = [f"{key[0]}:{key[1]}{key[2]}:{key[3]}:{key[4]}" for key in keys]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate site truth atom names")
    if not coordinates or not np.isfinite(np.asarray(coordinates, float)).all():
        raise ValueError("Empty or nonfinite site truth coordinates")
    return {"id": identifier, "coordinates": coordinates, "names": names,
            "atom_keys": keys,
            "interface_residue_keys": _contact_residues(
                tree, protein_coordinates, residue_keys, coordinates)}


def _check_contacts(site, expected):
    observed = {tuple(key) for key in site["interface_residue_keys"]}
    if observed != expected:
        raise ValueError(f"Site truth contact identities differ from frozen ledger: {site['id']}")


def _external_truth(case_id, raw_dir, record):
    definition = case_definition(case_id)
    payloads, source_files = {}, []
    for role, suffix in (("native_ligand", "l"), ("native_receptor", "r")):
        name = f"{case_id}_{suffix}_b.pdb"
        source = next(s for s in definition["sources"] if s["filename"] == name)
        payload, provenance = _read_pinned(
            Path(raw_dir) / name, source["sha256"], url=source["url"], role=role)
        payloads[suffix] = payload
        source_files.append(provenance)
    ligand = pdb_heavy_atoms(payloads["l"].decode())
    receptor = pdb_heavy_atoms(payloads["r"].decode())
    evaluation = record["evaluation_only"]
    if len(ligand) != evaluation["ligand_atoms"]:
        raise ValueError("Native ligand atom count differs from frozen site ledger")
    if evaluation["ligand_exact_deposited_atom_matches"] != len(ligand):
        raise ValueError("Frozen ledger lacks complete native ligand correspondence")
    expected_residues = {
        (r["prepared_chain"], r["prepared_resnum"], "", r["prepared_resname"])
        for r in evaluation["ligand_residue_map"]}
    if {tuple(atom["key"][:4]) for atom in ligand} != expected_residues:
        raise ValueError("Native ligand residue identity differs from frozen site ledger")
    if len(receptor) != evaluation["receptor_atom_count"]:
        raise ValueError("Native receptor atom count differs from frozen site ledger")
    # Canonical protein residues only: retained Ca2+ is not an interface residue.
    protein = [a for a in receptor if a["key"][3] in PROTEIN_RESIDUES]
    if not protein:
        raise ValueError("No canonical protein atoms in native receptor")
    xyz = np.asarray([a["coordinates"] for a in protein], dtype=float)
    keys = [a["key"][:4] for a in protein]
    site = _site(f"{case_id}:released_modeled_glycan", ligand, cKDTree(xyz), xyz, keys)
    expected_contacts = {
        (r["prepared_chain"], r["prepared_residue_number"], "", r["prepared_resname"])
        for r in evaluation["native_contact_residues_4A"]}
    _check_contacts(site, expected_contacts)
    site["deposited_ligand_residue_map"] = evaluation["ligand_residue_map"]
    return [site], source_files, len(protein), len({a["key"][0] for a in protein})


def _aav_truth(raw_dir, record):
    source = record["coordinate_source"]
    payload, provenance = _read_pinned(
        Path(raw_dir) / "7MTZ.cif", source["uncompressed_sha256"],
        url=source["url"], role="native_holo_complete_deposited_capsid")
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(payload.decode()).sole_block())
    if len(structure) != 1:
        raise ValueError("AAV site truth requires exactly one deposited model")
    assemblies = [assembly for assembly in structure.assemblies if assembly.name == "1"]
    if len(assemblies) != 1 or len(assemblies[0].generators) != 1:
        raise ValueError("AAV site truth assembly declaration changed")
    operators = assemblies[0].generators[0].operators
    if (len(operators) != 1 or not np.allclose(np.array(operators[0].transform.mat), np.eye(3))
            or not np.allclose(list(operators[0].transform.vec), 0)):
        raise ValueError("AAV site truth requires an already expanded identity assembly")
    protein_coordinates, residue_keys, protein_chains = [], [], set()
    ligands = []
    seen = set()
    for chain in structure[0]:
        for residue in chain:
            is_protein = residue.het_flag == "A" and residue.name in PROTEIN_RESIDUES
            if residue.het_flag == "A" and not is_protein:
                raise ValueError(f"Unsupported native AAV protein residue: {residue.name}")
            if not is_protein and residue.name != "GAL":
                continue
            atoms = []
            for atom in residue:
                if atom.element.is_hydrogen:
                    continue
                if atom.altloc not in ("\x00", " "):
                    raise ValueError("Native AAV alternate conformers require explicit curation")
                key = (chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, atom.name)
                if key in seen:
                    raise ValueError(f"Duplicate native AAV atom identity: {key}")
                seen.add(key)
                xyz = list(atom.pos)
                if not np.isfinite(xyz).all():
                    raise ValueError("Nonfinite native AAV coordinates")
                if is_protein:
                    protein_coordinates.append(xyz)
                    residue_keys.append(key[:4])
                    protein_chains.add(chain.name)
                else:
                    atoms.append({"key": list(key), "coordinates": xyz})
            if residue.name == "GAL":
                ligands.append((chain.name, residue.seqid.num, residue.subchain, atoms))
    if len(protein_chains) != record["receptor_scope"]["deposited_protein_chain_count"]:
        raise ValueError("Native AAV protein chain count differs from frozen site ledger")
    evaluation = record["evaluation_only"]
    expected = {(r["deposited_ligand_author_chain"], int(r["deposited_ligand_author_resnum"]),
                 r["deposited_ligand_label_chain"]): r for r in evaluation["site_copies"]}
    observed = {(chain, number, label) for chain, number, label, atoms in ligands}
    if len(ligands) != evaluation["symmetry_copy_count"] or observed != set(expected):
        raise ValueError("Native AAV GAL site identities/count differ from frozen site ledger")
    xyz = np.asarray(protein_coordinates, dtype=float)
    tree = cKDTree(xyz)  # Build once for the entire ~248k-atom deposited capsid.
    sites = []
    for chain, number, label, atoms in ligands:
        reference = expected[(chain, number, label)]
        if len(atoms) != reference["native_heavy_atom_count"]:
            raise ValueError("Native AAV GAL atom count differs from frozen site ledger")
        if {a["key"][4] for a in atoms} != {f"{e}{n}" for e in ("C", "O") for n in range(1, 7)}:
            raise ValueError("Native AAV GAL atom identities are incomplete")
        site = _site(f"7MTZ:{chain}:{number}:GAL", atoms, tree, xyz, residue_keys)
        site["deposited_label_chain"] = label
        expected_contacts = {(r["author_chain"], int(r["author_resnum"]), "", r["resname"])
                             for r in reference["native_contact_residues_4A"]}
        _check_contacts(site, expected_contacts)
        sites.append(site)
    return sites, [provenance], len(xyz), len(protein_chains)


def load_site_truth(case_id, raw_dir="data/raw", external_raw_dir=None, *, ledger_path=DEFAULT_LEDGER):
    """Return JSON-compatible native-site evaluation data from checksum-pinned sources.

    External files default to ``raw_dir/external``. AAV9 aliases use all 60 GAL
    copies in 7MTZ without further assembly expansion. Truth always remains in
    the holo source coordinate frame; apo predictions need an independently
    computed receptor alignment before evaluation. No prediction preparation,
    docking chemistry, download, repair, or coordinate transformation occurs.
    """
    ledger_path = Path(ledger_path)
    ledger_bytes = ledger_path.read_bytes()
    ledger = json.loads(ledger_bytes)
    canonical = str(case_id).upper()
    if canonical in {"AAV9", "7MTZ", "AAV9_7MTZ_7MT0"}:
        canonical = "aav9"
    matches = [r for r in ledger["cases"] if r["case_id"] == canonical]
    if len(matches) != 1:
        raise ValueError(f"Unknown or ambiguous site truth case: {case_id!r}")
    record = matches[0]
    if not record["eligible_site_only"] or record["atomic_pose_eligible"]:
        raise ValueError("Site truth requires the frozen site-only eligibility decision")
    if canonical == "aav9":
        sites, sources, atom_count, chain_count = _aav_truth(raw_dir, record)
    else:
        external_dir = Path(external_raw_dir) if external_raw_dir is not None else Path(raw_dir) / "external"
        sites, sources, atom_count, chain_count = _external_truth(canonical, external_dir, record)
    return {"schema_version": 1, "case_id": canonical, "pdb_id": record["pdb_id"],
            "sites": sites, "site_class_count": 1, "symmetry_copy_count": len(sites),
            "atomic_pose_eligible": False, "eligible_site_only": True,
            "confirmatory_site": False, "evaluation_only": True,
            "coordinate_frame": "pinned_holo_source",
            "apo_coordinate_transform_applied": False,
            "contact_cutoff_A": CONTACT_CUTOFF_A,
            "contact_rule": "canonical protein heavy atoms strictly less than 4.0 A from modeled native ligand heavy atoms",
            "native_protein_atom_count": atom_count, "native_protein_chain_count": chain_count,
            "site_reference_ledger": str(ledger_path.resolve()),
            "site_reference_ledger_sha256": sha256(ledger_bytes).hexdigest(),
            "source_files": sources, "site_quality_stratum": record["site_quality_stratum"],
            "reference_limitations": record["atomic_and_context_cautions"]}
