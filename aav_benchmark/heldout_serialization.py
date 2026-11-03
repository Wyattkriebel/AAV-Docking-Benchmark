"""Metadata-only reader compatibility amendment for structural test inputs.

The first twelve native reader attempts are retained.  This module changes no
protein or ligand geometry, author atom identity, element, occupancy, or B
factor.  It creates a fresh version with explicit polymer metadata and the CIF
comment delimiter required by the pinned fpocket/VMD reader.
"""
from __future__ import annotations

from collections import OrderedDict
import copy
import hashlib
import json
from pathlib import Path
import shutil

import gemmi

CASE_IDS = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
ALLOWED_ATOM_METADATA = frozenset(("label_asym_id", "label_entity_id", "label_seq_id"))
SEQUENCE_SCOPE = "observed residue names only; not full unresolved construct"
REPLACED_CATEGORIES = frozenset(("_entity.", "_entity_poly.", "_entity_poly_seq.", "_struct_asym.", "_pdbx_poly_seq_scheme."))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text())


def _write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _blank(value):
    return "" if value in (False, None, ".", "?", " ") else str(value)


def _rows(category):
    if not category: return []
    return [dict(zip(category, row)) for row in zip(*category.values())]


def _author_key(row):
    return tuple(_blank(row.get(field)) for field in ("pdbx_PDB_model_num", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code")) + (
        _blank(row.get("auth_comp_id", row["label_comp_id"])),
        _blank(row.get("auth_atom_id", row["label_atom_id"])), _blank(row.get("label_alt_id")))


def atom_semantics(path):
    """Index every atom-column value by its complete author identity."""
    atoms = gemmi.cif.read_file(str(path)).sole_block().get_mmcif_category("_atom_site.")
    result = {}
    for row in _rows(atoms):
        key = _author_key(row)
        if key in result: raise ValueError(f"Duplicate author atom identity: {key}")
        result[key] = row
    if not result: raise ValueError("Empty mmCIF atom category")
    return result


def compare_atom_semantics(before, after):
    """Permit only label metadata changes; all physical atom values are exact."""
    if set(before) != set(after): raise ValueError("Serialization changed author atom identities")
    changes = {field: 0 for field in ALLOWED_ATOM_METADATA}
    for key, old in before.items():
        new = after[key]
        if set(old) != set(new): raise ValueError("Serialization changed the atom-column schema")
        for field in old:
            if old[field] != new[field]:
                if field not in ALLOWED_ATOM_METADATA:
                    raise ValueError(f"Serialization changed protected atom field {field} at {key}")
                changes[field] += 1
    canonical = [[list(key), {field: value for field, value in before[key].items() if field not in ALLOWED_ATOM_METADATA}]
                 for key in sorted(before)]
    semantic_hash = hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"atom_count": len(before), "author_atom_identities_unchanged": True,
            "elements_unchanged": True, "coordinates_exactly_unchanged": True,
            "occupancies_exactly_unchanged": True, "B_factors_exactly_unchanged": True,
            "all_nonlabel_atom_fields_exactly_unchanged": True,
            "max_coordinate_displacement_A": 0.0, "atom_field_change_counts": changes,
            "protected_atom_semantics_sha256": semantic_hash}


def finalize_cif_metadata(source_path, output_path):
    """Write explicit observed-polymer metadata and a final atom-loop delimiter.

    Each prepared protein chain receives a numeric entity.  Label sequence
    numbers enumerate only the observed residues, preserving author numbers,
    insertion codes, names and source maps.  MSE remains the literal MSE
    monomer with its selenium element and original coordinates.
    """
    source_path, output_path = Path(source_path), Path(output_path)
    if output_path.exists() or output_path.is_symlink(): raise FileExistsError("Preserve prior serialization; choose a fresh CIF")
    original_hash = _sha(source_path)
    old = gemmi.cif.read_file(str(source_path)).sole_block()
    category = old.get_mmcif_category("_atom_site.")
    rows = _rows(category)
    if not rows: raise ValueError("Empty source protein")
    if {_blank(row.get("pdbx_PDB_model_num")) for row in rows} != {"1"}:
        raise ValueError("Serialization amendment requires the existing single model")
    if any(row["type_symbol"].upper() in ("H", "D") or _blank(row.get("label_alt_id")) for row in rows):
        raise ValueError("Use the frozen heavy-atom, whole-residue-alternate protein input")
    before = atom_semantics(source_path)
    chains, residues, label_owners = OrderedDict(), {}, {}
    for row in rows:
        chain = _blank(row["auth_asym_id"])
        if not chain: raise ValueError("Author chain identity is required")
        label = _blank(row.get("label_asym_id")) or chain
        if label in label_owners and label_owners[label] != chain:
            raise ValueError("Source label subchain ambiguously maps to author chains")
        label_owners[label] = chain
        if chain not in chains:
            chains[chain] = {"entity": str(len(chains)+1), "label": label, "residues": OrderedDict()}
        if chains[chain]["label"] != label:
            raise ValueError("Prepared author chain has multiple subchains; explicit curation required")
        residue_key = (_blank(row["auth_seq_id"]), _blank(row.get("pdbx_PDB_ins_code")), row["label_comp_id"])
        if residue_key not in chains[chain]["residues"]:
            chains[chain]["residues"][residue_key] = str(len(chains[chain]["residues"])+1)
        residues[(chain, residue_key)] = chains[chain]["residues"][residue_key]
        row["label_asym_id"] = label
        row["label_entity_id"] = chains[chain]["entity"]
        row["label_seq_id"] = residues[(chain, residue_key)]
    atom_category = {field: [row[field] for row in rows] for field in category}
    entities = {"id": [], "type": []}
    polymers = {"entity_id": [], "type": [], "pdbx_strand_id": [], "pdbx_seq_one_letter_code": []}
    sequences = {"entity_id": [], "num": [], "mon_id": [], "hetero": []}
    subchains = {"id": [], "entity_id": []}
    scheme = {field: [] for field in ("asym_id", "entity_id", "seq_id", "mon_id", "ndb_seq_num", "pdb_seq_num", "auth_seq_num", "pdb_mon_id", "auth_mon_id", "pdb_strand_id", "pdb_ins_code", "hetero")}
    chain_metadata = []
    for chain, definition in chains.items():
        entity, label = definition["entity"], definition["label"]
        names = [key[2] for key in definition["residues"]]
        # Literal nonstandard component names avoid silently relabeling MSE as MET.
        one_letter = "".join(gemmi.find_tabulated_residue(name).one_letter_code
                             if name in ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL") else f"({name})" for name in names)
        entities["id"].append(entity); entities["type"].append("polymer")
        for field, value in zip(polymers, (entity, "polypeptide(L)", chain, one_letter)): polymers[field].append(value)
        subchains["id"].append(label); subchains["entity_id"].append(entity)
        for (number, insertion, name), seq in definition["residues"].items():
            for field, value in zip(sequences, (entity, seq, name, "n")): sequences[field].append(value)
            values = (label, entity, seq, name, seq, number, number, name, name, chain, insertion or False, "n")
            for field, value in zip(scheme, values): scheme[field].append(value)
        chain_metadata.append({"author_chain": chain, "label_asym_id": label, "entity_id": entity,
                               "observed_residue_count": len(names), "observed_residue_names": names})
    document = gemmi.cif.Document(); block = document.add_new_block(old.name)
    for name in old.get_mmcif_category_names():
        if name != "_atom_site." and name not in REPLACED_CATEGORIES:
            block.set_mmcif_category(name, old.get_mmcif_category(name))
    for name, values in (("_entity.", entities), ("_entity_poly.", polymers), ("_entity_poly_seq.", sequences),
                         ("_struct_asym.", subchains), ("_pdbx_poly_seq_scheme.", scheme)):
        block.set_mmcif_category(name, values)
    # Keep atom_site last: fpocket's pinned VMD reader counts rows until '#'.
    block.set_mmcif_category("_atom_site.", atom_category)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document.as_string().rstrip() + "\n#\n")
    proof = compare_atom_semantics(before, atom_semantics(output_path))
    if _sha(source_path) != original_hash: raise AssertionError("Prior serialization was changed")
    return {"schema_version": 1, "serialization_revision": 2, "metadata_scope_only": True,
            "source_receptor_sha256": original_hash, "receptor_sha256": _sha(output_path),
            "source_receptor": str(source_path), "receptor": str(output_path),
            "serialization_sequence": SEQUENCE_SCOPE, "entity_count": len(chains), "chains": chain_metadata,
            "terminal_CIF_comment_delimiter": True,
            "changed_categories": sorted(REPLACED_CATEGORIES | {"_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id"}),
            "reason": "numeric polymer entities and complete observed subchain/label sequence metadata for BioJava; terminal atom-loop comment for fpocket/VMD",
            "native_ligand_coordinate_access": False, "chemical_identity_or_geometry_changed": False,
            "gemmi_version": gemmi.__version__, **proof}


def build_source_v2(primary, output):
    """Copy all strict artifacts, amending only six CIFs and their receipts."""
    primary, output = Path(primary), Path(output)
    if output.exists() or output.is_symlink(): raise FileExistsError("Preserve previous V2 artifacts")
    receipt = _load(primary / "source_preparation.json")
    if receipt["case_ids"] != list(CASE_IDS): raise ValueError("Keep the frozen six-case source cohort")
    original_files = {str(p.relative_to(primary)): _sha(p) for p in sorted(primary.rglob("*")) if p.is_file()}
    for case in receipt["cases"]:
        for relative, expected in case["artifact_sha256"].items():
            if original_files.get(f"{case['case_id']}/{relative}") != expected:
                raise ValueError("Primary source artifact differs from frozen preparation")
    shutil.copytree(primary, output)
    amended = copy.deepcopy(receipt); cases, permitted_changes = [], {"source_preparation.json"}
    for case in amended["cases"]:
        case_id = case["case_id"]; directory = output / case_id
        # Remove only the copied CIF, after verifying its identity to the prior file.
        cif = directory / "receptor.cif"
        if _sha(cif) != original_files[f"{case_id}/receptor.cif"]: raise AssertionError("Copy changed source bytes")
        cif.unlink()
        proof = finalize_cif_metadata(primary / case_id / "receptor.cif", cif)
        _write(directory / "serialization_v2.json", proof)
        for filename in ("receptor.cif.manifest.json", "protein_preparation.json"):
            manifest = _load(directory / filename)
            manifest.update(receptor_sha256=proof["receptor_sha256"], serialization_revision=2,
                            superseded_input_sha256=proof["source_receptor_sha256"],
                            serialization_revision_reason=proof["reason"],
                            serialization_revision_receipt=f"{case_id}/serialization_v2.json",
                            serialization_revision_receipt_sha256=_sha(directory / "serialization_v2.json"))
            _write(directory / filename, manifest)
        case["protein"] = _load(directory / "protein_preparation.json")
        case["serialization_revision"] = proof
        case["artifact_sha256"] = {str(p.relative_to(directory)): _sha(p) for p in sorted(directory.rglob("*"))
                                     if p.is_file() and p.name != "source_preparation.json"}
        _write(directory / "source_preparation.json", case)
        permitted_changes.update(f"{case_id}/{name}" for name in ("receptor.cif", "receptor.cif.manifest.json", "protein_preparation.json", "source_preparation.json"))
        cases.append({"case_id": case_id, **proof})
    amended.update(serialization_revision=2, prior_source_preparation_sha256=_sha(primary / "source_preparation.json"),
                   serialization_amendment_code_sha256=_sha(Path(__file__)),
                   prior_scientific_inputs_and_failures_retained=True, geometry_and_chemical_graphs_changed=False)
    _write(output / "source_preparation.json", amended)
    changed_files, unchanged_files = [], []
    for relative, old_hash in original_files.items():
        new_hash = _sha(output / relative)
        if new_hash != old_hash:
            if relative not in permitted_changes: raise ValueError(f"Unexpected nonmetadata change: {relative}")
            changed_files.append({"path": relative, "prior_sha256": old_hash, "sha256": new_hash})
        else: unchanged_files.append({"path": relative, "sha256": old_hash})
    if any(_sha(primary / relative) != old_hash for relative, old_hash in original_files.items()):
        raise AssertionError("Original source tree changed")
    proof = {"schema_version": 1, "amendment": "serialization_metadata_only", "fixed_case_ids": list(CASE_IDS),
             "prior_source_preparation_sha256": _sha(primary / "source_preparation.json"),
             "source_preparation_v2_sha256": _sha(output / "source_preparation.json"),
             "original_file_count": len(original_files), "changed_original_files": changed_files,
             "unchanged_original_files": unchanged_files, "prior_tree_exactly_unchanged": True,
             "all_PDB_PDBQT_SDF_files_byte_identical": all(_sha(output / relative) == old_hash for relative, old_hash in original_files.items()
                                                           if Path(relative).suffix.lower() in (".pdb", ".pdbqt", ".sdf")),
             "all_author_atoms_elements_coordinates_exactly_unchanged": True,
             "native_ligand_graphs_and_controls_byte_identical": True,
             "native_ligand_coordinate_access": False, "cases": cases,
             "code_sha256": _sha(Path(__file__))}
    _write(output / "serialization_amendment.json", proof)
    return proof
