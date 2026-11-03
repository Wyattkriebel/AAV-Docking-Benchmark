"""Source-only audit and replay for the fixed twelve-candidate structural inventory.

The selected deposited copy is the native source. No prepared ligand or predicted
pose is read or invented. Original branch chemistry, external covalent aglycones,
alternate rows, validation records, and average occupancy semantics are retained.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import numpy as np

from .atomic_reference import (
    AtomicReferenceError, DEFAULT_ROOT, DEFAULT_DENSITY_TOLERANCE_SIGMA,
    _ccd, _chirality, _compare_replay, _path, _read_grid, _residue_failures,
    _same, _sample, _sha, _xyz,
)

DEFAULT_INPUTS = Path("data/archive/heldout_reference/audit_inputs.json")
DEFAULT_LEDGER = Path("benchmarks/heldout_reference_ledger.json")
CASE_IDS = ("1I3H", "2G7C", "2GUP", "2RDK", "2WAB", "2WNB", "3A2S", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
IDENTITY_FIELDS = ("pdbx_PDB_model_num", "label_asym_id", "auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id", "label_alt_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy")


def _load(path):
    return json.loads(Path(path).read_text())


def _rows(category):
    if not category:
        return []
    return [{key: values[i] for key, values in category.items()} for i in range(len(next(iter(category.values()))))]


def _pin(root, source):
    relative = source.get("archive", source.get("archive_path"))
    path = _path(root, relative)
    digest = source.get("archive_sha256", source.get("sha256"))
    if not path.is_file() or _sha(path) != digest:
        raise AtomicReferenceError(f"checksum mismatch: {relative}")
    size = source.get("archive_bytes", source.get("bytes"))
    if size is not None and path.stat().st_size != size:
        raise AtomicReferenceError(f"byte count mismatch: {relative}")
    # Candidate archives distinguish original source bytes from gzip storage bytes.
    if source.get("archive_encoding", "").startswith("gzip mtime"):
        if hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() != source["sha256"]:
            raise AtomicReferenceError(f"uncompressed checksum mismatch: {relative}")
    return path


def _blank(value):
    return "" if value in [None, False, ".", "?", " "] else str(value)


def _identity(category, indices):
    rows = [[category[key][i] for key in IDENTITY_FIELDS] for i in indices]
    rows.sort(key=lambda row: tuple(map(str, row)))
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def _key(category, i):
    return f"{category['label_asym_id'][i]}:{category['auth_seq_id'][i]}:{category['label_comp_id'][i]}:{category['label_atom_id'][i]}:{_blank(category['label_alt_id'][i])}"


def audit_candidate(entry: dict, root: Path | str, policy: dict) -> dict:
    """Compute density, complete native graph, validation and stereo for one copy.

    Source pins are checked on every call. Unknown metrics remain unknown, never
    passing. Every alternate atom row is retained. This function makes numerical
    source gates only; the curator separately freezes paper/assembly/completeness
    judgments and scored whole-residue components.
    """
    root = Path(root)
    coordinate_path = _pin(root, entry["sources"]["coordinates"])
    validation_path = _pin(root, entry["sources"]["validation"])
    for source in entry["ccd_sources"].values():
        _pin(root, source)
    for source in entry["audit_sources"].values():
        if source.get("status", "retrieved") == "retrieved":
            _pin(root, source)
    block = gemmi.cif.read_file(str(coordinate_path))[0]
    category = block.get_mmcif_category("_atom_site.")
    selected = [i for i in range(len(category["id"]))
                if category["label_asym_id"][i] == entry["selected_label_asym_id"]
                and category["type_symbol"][i] not in ["H", "D"]]
    source_identity = _identity(category, selected)
    branch_links = [row for row in _rows(block.get_mmcif_category("_pdbx_entity_branch_link.")) if row["entity_id"] == entry["branch_entity_id"]]
    _same(branch_links, entry["branch_links"], entry["pdb_id"] + " branch graph inventory")
    connections = [row for row in _rows(block.get_mmcif_category("_struct_conn.")) if row["conn_type_id"] == "covale"]
    # Include every explicitly inventoried covalent aglycone residue as part of
    # the chemical parent, rather than silently treating its attachment as OH.
    specs = []
    for residue in entry["residues"]:
        specs.append(dict(branch_num=residue["branch_num"], ccd_id=residue["ccd_id"],
                          author_chain=residue["current_atom_site_auth_chain"],
                          author_resnum=residue["current_atom_site_auth_seq_id"],
                          label_chain=entry["selected_label_asym_id"], external_aglycone=False,
                          expected_heavy_atoms=residue["expected_heavy_atoms_after_leaving_atom_removal"]))
    for connection in entry["extra_covalent_edges"]:
        for partner in ["ptnr1", "ptnr2"]:
            chain = connection[partner + "_label_asym_id"]
            if chain == entry["selected_label_asym_id"]:
                continue
            spec = dict(branch_num="external:" + chain + ":" + connection[partner + "_auth_seq_id"],
                        ccd_id=connection[partner + "_label_comp_id"],
                        author_chain=connection[partner + "_auth_asym_id"],
                        author_resnum=connection[partner + "_auth_seq_id"], label_chain=chain,
                        external_aglycone=True, expected_heavy_atoms=None)
            if spec not in specs:
                specs.append(spec)
    parent_indices = [i for i in range(len(category["id"])) if category["type_symbol"][i] not in ["H", "D"]
                      and any(category["label_asym_id"][i] == spec["label_chain"]
                              and category["auth_seq_id"][i] == spec["author_resnum"]
                              and category["label_comp_id"][i] == spec["ccd_id"] for spec in specs)]
    glycan_identity = source_identity
    source_identity = _identity(category, parent_indices)
    _same(source_identity, entry["atom_coordinate_identity_sha256"], entry["pdb_id"] + " native parent coordinate identity")
    if len({_key(category, i) for i in parent_indices}) != len(parent_indices):
        raise AtomicReferenceError(f"Duplicate native atom key: {entry['pdb_id']}")
    required = policy["required_per_scored_residue"]
    maps, map_limits = {}, []
    for kind, source_key in [("2fo-fc", "2fofc_coefficients"), ("fo-fc", "fofc_coefficients")]:
        source = entry["audit_sources"].get(source_key)
        if source and source.get("status", "retrieved") == "retrieved":
            maps[kind] = _read_grid(_path(root, source.get("archive", source.get("archive_path"))), kind, required["map_fft_sample_rate"])
        else:
            map_limits.append(kind + "_map_unavailable")
    step, radius = required["local_density_sampling_grid_A"], required["local_density_search_radius_A"]
    nstep = round(radius / step)
    offsets = [np.array(offset) * step for offset in itertools.product(range(-nstep, nstep + 1), repeat=3)
               if sum(value * value for value in offset) * step * step <= radius * radius]
    xml = ET.fromstring(gzip.open(validation_path, "rt").read())
    xml_entry = xml.find("Entry").attrib
    stereo_links_raw = [dict(left=[c["ptnr1_auth_asym_id"], c["ptnr1_auth_seq_id"], c["ptnr1_label_atom_id"]],
                         right=[c["ptnr2_auth_asym_id"], c["ptnr2_auth_seq_id"], c["ptnr2_label_atom_id"]],
                         distance_A=float(c["pdbx_dist_value"])) for c in connections]
    # Alternate-specific LINK records can repeat the same named residue edge.
    # They represent one stereochemical substitution, not two oxygen candidates.
    stereo_links = []
    for link in stereo_links_raw:
        if not any(link["left"] == old["left"] and link["right"] == old["right"] for old in stereo_links):
            stereo_links.append(link)
    residues, native_atoms, native_bonds = [], [], []
    for spec in specs:
        indices = [i for i in parent_indices if category["label_asym_id"][i] == spec["label_chain"]
                   and category["auth_seq_id"][i] == spec["author_resnum"] and category["label_comp_id"][i] == spec["ccd_id"]]
        ccd_path = _path(root, entry["ccd_sources"][spec["ccd_id"]].get("archive_path", entry["ccd_sources"][spec["ccd_id"]].get("archive")))
        ccd_atoms, ideal, edges = _ccd(ccd_path)
        atom_info = {name: i for i, name in enumerate(ccd_atoms["atom_id"])}
        expected = spec["expected_heavy_atoms"]
        if expected is None:
            expected = sorted(name for i, name in enumerate(ccd_atoms["atom_id"]) if ccd_atoms["type_symbol"][i] not in ["H", "D"])
        observed = sorted(set(category["label_atom_id"][i] for i in indices))
        atom_set_match = observed == sorted(expected)
        alternates = sorted(set(_blank(category["label_alt_id"][i]) for i in indices) - {""})
        stereo = []
        for alternate in alternates or [""]:
            coordinates = {category["label_atom_id"][i]: _xyz(category, i) for i in indices
                           if _blank(category["label_alt_id"][i]) in ["", alternate]}
            values = _chirality(ccd_atoms, ideal, edges, coordinates,
                                dict(author_chain=spec["author_chain"], author_resnum=spec["author_resnum"]),
                                stereo_links, category)
            stereo.extend(dict(alt_id=alternate, **value) for value in values)
        records = [x for x in xml.findall("ModelledSubgroup") if x.attrib.get("chain") == spec["author_chain"]
                   and x.attrib.get("resnum") == spec["author_resnum"] and x.attrib.get("resname") == spec["ccd_id"]]
        validation_records = [dict(attributes=dict(record.attrib), flags=[dict(kind=x.tag, **x.attrib) for x in record]) for record in records]
        validation = {}
        for key in ["rscc", "rsr", "avgoccu", "owab", "EDIAm", "OPIA"]:
            values = [float(record.attrib[key]) for record in records if key in record.attrib]
            if len(values) == len(records) and values:
                validation[key] = (max(values) if key in ["rsr", "owab"] else min(values))
        flags = [flag for record in validation_records for flag in record["flags"]]
        atoms = []
        for i in indices:
            name, position = category["label_atom_id"][i], _xyz(category, i)
            ci = atom_info[name]
            density = dict(two_mFo_DFc_atom_sigma=None, two_mFo_DFc_max_within_0p5A_sigma=None, mFo_DFc_atom_sigma=None)
            if len(maps) == 2:
                density = dict(two_mFo_DFc_atom_sigma=_sample(maps["2fo-fc"][0], position),
                               two_mFo_DFc_max_within_0p5A_sigma=max(_sample(maps["2fo-fc"][0], position + offset) for offset in offsets),
                               mFo_DFc_atom_sigma=_sample(maps["fo-fc"][0], position))
            atoms.append(dict(native_atom_key=_key(category, i), branch_num=spec["branch_num"], ccd_id=spec["ccd_id"],
                              atom_name=name, element=category["type_symbol"][i], formal_charge=int(ccd_atoms["charge"][ci]), aromatic=ccd_atoms["pdbx_aromatic_flag"][ci] == "Y",
                              ccd_stereo_config=ccd_atoms["pdbx_stereo_config"][ci], model_num=category["pdbx_PDB_model_num"][i],
                              deposited_label_chain=category["label_asym_id"][i], deposited_label_resnum=category["label_seq_id"][i],
                              deposited_author_chain=category["auth_asym_id"][i], deposited_author_resnum=category["auth_seq_id"][i],
                              deposited_atom_site_id=category["id"][i], alt_id=_blank(category["label_alt_id"][i]),
                              xyz_A=position.tolist(), occupancy=float(category["occupancy"][i]), B_iso_A2=float(category["B_iso_or_equiv"][i]),
                              **density))
        if len(maps) == 2 and atoms:
            # The old shared gate also has a released-source correspondence test.
            # Here the source is the deposition itself; no prepared ligand exists.
            gate_atoms = [dict(atom, released_coordinate_delta_A=0.0) for atom in atoms]
            failures, local = _residue_failures(validation, gate_atoms, flags, stereo, policy)
        else:
            failures, local = list(map_limits), None
            # Missing maps do not hide independently inspectable geometry/fit gaps.
            for key, bound, low in [("rscc", required["RSCC_min"], True), ("rsr", required["RSR_max"], False), ("avgoccu", required["occupancy_min"], True)]:
                if key not in validation or (validation[key] < bound if low else validation[key] > bound):
                    failures.append(key + ("_min_gate" if low else "_max_gate"))
        if not atom_set_match:
            failures.append("CCD_link_atom_set_mismatch")
        if len(alternates) > 1:
            failures.append("multiple_native_alternates_no_unique_atomic_target")
        if not records:
            failures.append("validation_record_unavailable")
        full_ccd = {name for i, name in enumerate(ccd_atoms["atom_id"]) if ccd_atoms["type_symbol"][i] not in ["H", "D"]}
        residue = dict(branch_num=spec["branch_num"], ccd_id=spec["ccd_id"], external_aglycone=spec["external_aglycone"],
                       deposited=dict(author_chain=spec["author_chain"], author_resnum=spec["author_resnum"], label_chain=spec["label_chain"], resname=spec["ccd_id"]),
                       heavy_atom_row_count=len(atoms), unique_heavy_atom_count=len(observed), expected_heavy_atoms=sorted(expected),
                       observed_heavy_atoms=observed, atom_set_match=atom_set_match, missing_CCD_heavy_atoms=sorted(full_ccd - set(observed)),
                       alternate_ids=alternates, validation=validation, validation_records=validation_records, validation_flags=flags,
                       validation_occupancy_rule="Minimum reported avgoccu across the residue's alternate records; no all-atom occupancy assertion.",
                       atom_occupancy_range=[min(a["occupancy"] for a in atoms), max(a["occupancy"] for a in atoms)] if atoms else [],
                       stereochemistry=stereo, local_density=local, source_gate_pass=not failures, source_gate_failures=failures,
                       atom_map=atoms, CCD_internal_heavy_bonds=[list(edge) for edge in edges if edge[0] in observed and edge[1] in observed])
        residues.append(residue)
        native_atoms.extend(atoms)
        ccd_bond_rows = _rows(gemmi.cif.read_file(str(ccd_path))[0].get_mmcif_category("_chem_comp_bond."))
        for edge in ccd_bond_rows:
            for left in atoms:
                if left["atom_name"] != edge["atom_id_1"]:
                    continue
                for right in atoms:
                    if right["atom_name"] != edge["atom_id_2"] or (left["alt_id"] and right["alt_id"] and left["alt_id"] != right["alt_id"]):
                        continue
                    native_bonds.append(dict(atom1=left["native_atom_key"], atom2=right["native_atom_key"], order=edge["value_order"].upper(),
                                             aromatic=edge.get("pdbx_aromatic_flag") == "Y", kind="CCD_internal", ccd_id=spec["ccd_id"],
                                             distance_A=float(np.linalg.norm(np.array(left["xyz_A"]) - right["xyz_A"]))))
    selected_connections = []
    for connection in connections:
        hits = []
        for partner in ["ptnr1", "ptnr2"]:
            required_alt = _blank(connection.get("pdbx_" + partner + "_label_alt_id"))
            hits.append([atom for atom in native_atoms if (not required_alt or atom["alt_id"] == required_alt)
                         and atom["deposited_label_chain"] == connection[partner + "_label_asym_id"]
                         and atom["deposited_author_resnum"] == connection[partner + "_auth_seq_id"]
                         and atom["atom_name"] == connection[partner + "_label_atom_id"]])
        if all(hits):
            selected_connections.append(connection)
            for left in hits[0]:
                for right in hits[1]:
                    if left["alt_id"] and right["alt_id"] and left["alt_id"] != right["alt_id"]:
                        continue
                    native_bonds.append(dict(atom1=left["native_atom_key"], atom2=right["native_atom_key"], order=(connection.get("pdbx_value_order") or "SING").upper(),
                                             aromatic=False, kind="struct_conn", struct_conn_id=connection["id"],
                                             struct_conn_distance_A=float(connection["pdbx_dist_value"]),
                                             distance_A=float(np.linalg.norm(np.array(left["xyz_A"]) - right["xyz_A"]))))
    return dict(case_id=entry["pdb_id"], pdb_id=entry["pdb_id"], source_only=True, prediction_outcomes_consulted=False,
                selected_label_asym_id=entry["selected_label_asym_id"], selected_copy_identity_sha256=source_identity, selected_glycan_only_identity_sha256=glycan_identity,
                maps={kind: value[1] for kind, value in maps.items()}, map_limits=map_limits, local_sphere_sample_count=len(offsets),
                entry_validation={key: xml_entry.get(key) for key in ["PDB-resolution", "PDB-R", "EDS_R", "DCC_R", "PDB-Rfree", "DCC_Rfree", "XMLcreationDate", "PDB-revision-date"]},
                branch_links=branch_links, covalent_links=selected_connections, residues=residues,
                native_atom_map=native_atoms, native_bond_graph=native_bonds,
                native_heavy_atom_row_count=len(native_atoms), native_unique_heavy_atom_count=len({tuple(a[key] for key in ["deposited_label_chain", "deposited_author_resnum", "ccd_id", "atom_name"]) for a in native_atoms}),
                graph_rule="Complete deposited selected glycan plus inventoried covalent aglycones; standard CCD names and bonds, declared branch links; every alternate row retained. No prepared graph is inferred from a truncated core.")


def derive_copy_entries(entry: dict, root: Path | str) -> dict:
    """Derive full-chemistry ASU copies and retain other carbohydrate site rows.

    Membership expressions and operator matrices are preserved verbatim. Rigid
    assembly instances must be frozen separately; crystallographic copies from
    different biological assembly IDs are never merged by this helper.
    """
    from copy import deepcopy

    root = Path(root)
    block = gemmi.cif.read_file(str(_pin(root, entry["sources"]["coordinates"])))[0]
    atoms = block.get_mmcif_category("_atom_site.")
    scheme = _rows(block.get_mmcif_category("_pdbx_branch_scheme."))
    branch_lists = _rows(block.get_mmcif_category("_pdbx_entity_branch_list."))
    all_links = _rows(block.get_mmcif_category("_pdbx_entity_branch_link."))
    connections = [row for row in _rows(block.get_mmcif_category("_struct_conn.")) if row["conn_type_id"] == "covale"]
    generators = _rows(block.get_mmcif_category("_pdbx_struct_assembly_gen."))

    def extra_connections(label):
        return [row for row in connections if ((row["ptnr1_label_asym_id"] == label) != (row["ptnr2_label_asym_id"] == label))]

    def signature(entity, label):
        monomers = sorted((row["num"], row["comp_id"]) for row in branch_lists if row["entity_id"] == entity)
        links = []
        for row in all_links:
            if row["entity_id"] != entity:
                continue
            ends = sorted([(row["entity_branch_list_num_1"], row["comp_id_1"], row["atom_id_1"], row["leaving_atom_id_1"]),
                           (row["entity_branch_list_num_2"], row["comp_id_2"], row["atom_id_2"], row["leaving_atom_id_2"])])
            links.append((tuple(ends), row["value_order"].upper()))
        external = []
        for row in extra_connections(label):
            this, other = ("ptnr1", "ptnr2") if row["ptnr1_label_asym_id"] == label else ("ptnr2", "ptnr1")
            residue = next(r for r in scheme if r["asym_id"] == label and r["pdb_seq_num"] == row[this + "_auth_seq_id"])
            external.append((residue["num"], row[this + "_label_comp_id"], row[this + "_label_atom_id"],
                             row[other + "_label_comp_id"], row[other + "_label_atom_id"]))
        return monomers, sorted(links), sorted(external)

    primary_signature = signature(entry["branch_entity_id"], entry["selected_label_asym_id"])
    copies, annotations = [], []
    branch_labels = sorted(set(row["asym_id"] for row in scheme))
    full_labels = set()
    for label in branch_labels:
        entity = next(row["entity_id"] for row in scheme if row["asym_id"] == label)
        if signature(entity, label) != primary_signature:
            continue
        copied = deepcopy(entry)
        copied["selected_label_asym_id"] = label
        copied["branch_entity_id"] = entity
        copied["branch_links"] = [row for row in all_links if row["entity_id"] == entity]
        copied["extra_covalent_edges"] = extra_connections(label)
        copied["copy_id"] = entry["pdb_id"] + ":" + label
        copied["same_full_chemical_graph_as_primary"] = True
        copied["assembly_memberships"] = [row for row in generators if label in row["asym_id_list"].split(",")]
        for residue in copied["residues"]:
            row = next(r for r in scheme if r["asym_id"] == label and r["num"] == residue["branch_num"])
            residue["current_atom_site_auth_chain"] = row["pdb_asym_id"]
            residue["current_atom_site_auth_seq_id"] = row["pdb_seq_num"]
            residue["original_depositor_auth_chain_from_branch_scheme"] = row["auth_asym_id"]
            residue["original_depositor_auth_seq_id_from_branch_scheme"] = row["auth_seq_num"]
        labels = {label}
        for row in copied["extra_covalent_edges"]:
            labels.update([row["ptnr1_label_asym_id"], row["ptnr2_label_asym_id"]])
        indices = [i for i in range(len(atoms["id"])) if atoms["label_asym_id"][i] in labels and atoms["type_symbol"][i] not in ["H", "D"]]
        copied["atom_coordinate_identity_sha256"] = _identity(atoms, indices)
        copied["copy_identity_rule"] = "Derived exact source rows before outcome review; canonical branch CCD/link graph plus external covalent identities equal primary."
        copies.append(copied)
        full_labels.update(labels)
    component_types = {row["id"]: row["type"] for row in _rows(block.get_mmcif_category("_chem_comp."))}
    other_labels = sorted({atoms["label_asym_id"][i] for i in range(len(atoms["id"]))
                           if atoms["label_asym_id"][i] not in full_labels
                           and (atoms["label_asym_id"][i] in branch_labels or "saccharide" in str(component_types.get(atoms["label_comp_id"][i], "")).lower())})
    for label in other_labels:
        rows = [i for i in range(len(atoms["id"])) if atoms["label_asym_id"][i] == label and atoms["type_symbol"][i] not in ["H", "D"]]
        annotations.append(dict(label_asym_id=label, same_full_chemical_graph_as_primary=False, scored=False,
                                reason="Partial or different chemical carbohydrate copy; known site annotation only.",
                                assembly_memberships=[row for row in generators if label in row["asym_id_list"].split(",")],
                                atom_map=[dict(native_atom_key=_key(atoms, i), atom_name=atoms["label_atom_id"][i],
                                               ccd_id=atoms["label_comp_id"][i], author_chain=atoms["auth_asym_id"][i],
                                               author_resnum=atoms["auth_seq_id"][i], alt_id=_blank(atoms["label_alt_id"][i]),
                                               xyz_A=_xyz(atoms, i).tolist()) for i in rows]))
    return dict(source_only=True, prediction_outcomes_consulted=False, same_full_chemical_entries=copies,
                other_carbohydrate_sites=annotations, assembly_generators=generators,
                assembly_operators=_rows(block.get_mmcif_category("_pdbx_struct_oper_list.")))


def audit_all(inputs_path=DEFAULT_INPUTS, root=DEFAULT_ROOT):
    """Return all twelve diagnostics without writing them."""
    root = Path(root)
    inputs = _load(_path(root, inputs_path))
    _same(tuple(e["pdb_id"] for e in inputs["candidates"]), CASE_IDS, "Fixed twelve-candidate denominator")
    _same(_sha(_path(root, inputs["inventory_path"])), inputs["inventory_sha256"], "Inventory checksum")
    _same(_sha(_path(root, inputs["policy_path"])), inputs["policy_sha256"], "Frozen policy checksum")
    _same(inputs["source_only"], True, "Source-only audit input")
    _same(inputs["prediction_outcomes_consulted"], False, "Audit input outcome flag")
    policy = _load(_path(root, inputs["policy_path"]))
    return dict(schema_version=1, source_only=True, prediction_outcomes_consulted=False,
                audit_inputs_sha256=_sha(_path(root, inputs_path)), policy_sha256=inputs["policy_sha256"],
                software=dict(gemmi=gemmi.__version__, numpy=np.__version__),
                cases=[audit_candidate(entry, root, policy) for entry in inputs["candidates"]])


def qualified_component(diagnostic: dict, original_reagent_complete: bool):
    """Apply unchanged whole-sugar connected-core rule; no curation is invented."""
    residues = diagnostic["residues"]
    passed = {r["branch_num"] for r in residues if r["source_gate_pass"] and not r["external_aglycone"]}
    counts = {r["branch_num"]: r["heavy_atom_row_count"] for r in residues}
    by_atom = {a["native_atom_key"]: a["branch_num"] for a in diagnostic["native_atom_map"]}
    edges = [(by_atom[b["atom1"]], by_atom[b["atom2"]]) for b in diagnostic["native_bond_graph"]]
    remaining, components = set(passed), []
    while remaining:
        component, queue = {min(remaining)}, [min(remaining)]
        while queue:
            residue = queue.pop()
            for left, right in edges:
                neighbor = right if left == residue else left if right == residue else None
                if neighbor in passed and neighbor not in component:
                    component.add(neighbor)
                    queue.append(neighbor)
        remaining -= component
        components.append(sorted(component))
    components.sort(key=lambda values: (-sum(counts[r] for r in values), tuple(values)))
    full = bool(residues) and all(r["source_gate_pass"] for r in residues) and original_reagent_complete
    selected = [r["branch_num"] for r in residues] if full else (components[0] if components else [])
    return selected, "full_ligand" if full else "modeled_core" if selected else "rejected"


def canonical_atom_key(atom: dict) -> str:
    return f"{atom['branch_num']}:{atom['ccd_id']}:{atom['atom_name']}"


def _scored_atoms(diagnostic, selected):
    atoms = [dict(atom, canonical_atom_key=canonical_atom_key(atom)) for atom in diagnostic["native_atom_map"] if atom["branch_num"] in selected]
    if len({a["canonical_atom_key"] for a in atoms}) != len(atoms):
        raise AtomicReferenceError("Scored native mask contains unresolved alternate atoms")
    return atoms


def _assembly_operators(entry, root, assembly_id):
    structure = gemmi.read_structure(str(_pin(root, entry["sources"]["coordinates"])))
    assembly = next((item for item in structure.assemblies if item.name == assembly_id), None)
    if assembly is None:
        raise AtomicReferenceError(f"Missing source assembly {entry['pdb_id']}:{assembly_id}")
    return [dict(operator_name=operator.name, subchains=list(generator.subchains),
                 matrix=operator.transform.mat.tolist(), vector=operator.transform.vec.tolist())
            for generator in assembly.generators for operator in generator.operators]


def verify_heldout_references(root=DEFAULT_ROOT, *, ledger_path=DEFAULT_LEDGER,
                              density_tolerance_sigma=DEFAULT_DENSITY_TOLERANCE_SIGMA):
    """Replay all twelve native graphs, fixed masks and qualified assembly copies."""
    root = Path(root)
    ledger_file = _path(root, ledger_path)
    ledger = _load(ledger_file)
    if not np.isfinite(density_tolerance_sigma) or density_tolerance_sigma < 0:
        raise ValueError("Density tolerance must be finite and nonnegative")
    for prefix in ["policy", "manual_adjudications", "inventory", "audit_inputs", "source_diagnostics",
                   "copy_inputs", "all_copy_diagnostics", "assembly_context"]:
        _same(_sha(_path(root, ledger[prefix + "_archive"])), ledger[prefix + "_sha256"], prefix + " checksum")
    _same(tuple(c["case_id"] for c in ledger["cases"]), CASE_IDS, "Frozen ledger denominator/order")
    _same(ledger["cohort_denominator"], 12, "Frozen twelve-case denominator")
    _same(ledger["source_only"], True, "Ledger source-only flag")
    _same(ledger["prediction_outcomes_consulted"], False, "Ledger outcome flag")
    manifest = _load(_path(root, ledger["source_record"]))
    _same(manifest["ledger_sha256"], _sha(ledger_file), "Frozen ledger checksum")
    _same(manifest["source_only"], True, "Manifest source-only flag")
    _same(manifest["prediction_outcomes_consulted"], False, "Manifest outcome flag")
    for source in manifest["archived_files"]:
        _pin(root, source)
    inputs = _load(_path(root, ledger["audit_inputs_archive"]))
    policy = _load(_path(root, ledger["policy_archive"]))
    actual = audit_all(ledger["audit_inputs_archive"], root)
    frozen = _load(_path(root, ledger["source_diagnostics_archive"]))
    deltas = []
    _compare_replay(actual, frozen, "heldout primary source diagnostics", density_tolerance_sigma, deltas)
    manual = _load(_path(root, ledger["manual_adjudications_archive"]))
    _same(manual["source_only"], True, "Manual source-only flag")
    _same(manual["prediction_outcomes_consulted"], False, "Manual outcome flag")
    assembly_context = _load(_path(root, ledger["assembly_context_archive"]))
    contexts = {case["case_id"]: case for case in assembly_context["cases"]}
    _same(assembly_context["source_only"], True, "Assembly context source-only flag")
    _same(assembly_context["prediction_outcomes_consulted"], False, "Assembly context outcome flag")
    copy_inputs = _load(_path(root, ledger["copy_inputs_archive"]))
    actual_copy_entries, actual_copy_cases = [], []
    for entry in inputs["candidates"]:
        derived = derive_copy_entries(entry, root)
        actual_copy_entries.extend(derived["same_full_chemical_entries"])
        copies = [dict(copy_id=copy["copy_id"], assembly_memberships=copy["assembly_memberships"],
                       diagnostics=audit_candidate(copy, root, policy)) for copy in derived["same_full_chemical_entries"]]
        actual_copy_cases.append(dict(case_id=entry["pdb_id"], copies=copies,
                                      other_carbohydrate_sites=derived["other_carbohydrate_sites"],
                                      assembly_generators=derived["assembly_generators"], assembly_operators=derived["assembly_operators"]))
    _same(copy_inputs, dict(schema_version=1, source_only=True, prediction_outcomes_consulted=False, copies=actual_copy_entries), "Derived full-chemical copy inputs")
    actual_all_copies = dict(schema_version=1, source_only=True, prediction_outcomes_consulted=False, cases=actual_copy_cases)
    frozen_all_copies = _load(_path(root, ledger["all_copy_diagnostics_archive"]))
    _compare_replay(actual_all_copies, frozen_all_copies, "all full-chemical source copies", density_tolerance_sigma, deltas)
    reports, total_scored, total_qualified = [], 0, 0
    for case, diagnostic, entry, copy_case in zip(ledger["cases"], actual["cases"], inputs["candidates"], actual_copy_cases):
        label = case["case_id"]
        adjudication = manual["decisions"][label]
        _same(case["inventory_entry"], entry, label + " pinned source inventory entry")
        _same(case["source_only"], True, label + " source-only flag")
        _same(case["prediction_outcomes_consulted"], False, label + " outcome flag")
        _same(case["decision_reason"], adjudication["reason"], label + " pinned manual reason")
        _same(case["original_reagent_complete"], adjudication["full_original_reagent_complete"], label + " original reagent completeness judgment")
        _same(case["source_preparation_blocks"], adjudication["source_preparation_blocks"], label + " source preparation judgments")
        _compare_replay(diagnostic["native_atom_map"], case["native_atom_map"], label + " native atom map", density_tolerance_sigma, deltas)
        _compare_replay(diagnostic["native_bond_graph"], case["native_bond_graph"], label + " native bond graph", density_tolerance_sigma, deltas)
        _same(case["native_heavy_atom_row_count"], diagnostic["native_heavy_atom_row_count"], label + " native atom row count")
        _same(case["native_unique_heavy_atom_count"], diagnostic["native_unique_heavy_atom_count"], label + " unique native atom count")
        selected, scope = qualified_component(diagnostic, adjudication["full_original_reagent_complete"])
        _same(selected, case["scored_residue_ids"], label + " whole-residue scored component")
        _same(scope, case["qualification"], label + " qualification scope")
        _same(bool(selected), case["atomic_pose_eligible"], label + " atomic eligibility")
        _same(scope == "full_ligand", case["full_ligand_eligible"], label + " full-ligand eligibility")
        _same(scope == "modeled_core", case["modeled_core_eligible"], label + " modeled-core eligibility")
        _same(bool(selected) and not adjudication["source_preparation_blocks"], case["source_conditioned_preparation_eligible"], label + " source-conditioned preparation eligibility")
        _same(case["parent_ligand_graph_preserved"], True, label + " complete parent graph preserved")
        _same(case["chemical_repair_authorized"], False, label + " no implicit native chemistry repair")
        scored = _scored_atoms(diagnostic, selected)
        _compare_replay(scored, case["scored_atom_map"], label + " primary scored map", density_tolerance_sigma, deltas)
        _same(len(scored), case["scored_atom_count"], label + " primary scored atom count")
        context = contexts[label]
        _same(context, case["receptor_scope"], label + " pinned receptor scope")
        operators = _assembly_operators(entry, root, context["selected_assembly_id"])
        _same(operators, context["expanded_assembly_operators"], label + " source biological-assembly operators")
        _same(copy_case["assembly_generators"], context["assembly_generators"], label + " source assembly generators")
        _same(copy_case["other_carbohydrate_sites"], case["other_known_carbohydrate_sites"], label + " unscored partial/different known sites")
        source_decisions, qualified = [], []
        wanted = {atom["canonical_atom_key"] for atom in scored}
        for copy in copy_case["copies"]:
            source = copy["diagnostics"]
            in_assembly = any(membership["assembly_id"] == context["selected_assembly_id"] for membership in copy["assembly_memberships"])
            residues = {residue["branch_num"]: residue for residue in source["residues"]}
            mapping = {canonical_atom_key(atom): atom for atom in source["native_atom_map"]}
            unique_mask = all(sum(canonical_atom_key(atom) == key for atom in source["native_atom_map"]) == 1 for key in wanted)
            mask_pass = bool(selected) and unique_mask and all(residues[number]["source_gate_pass"] for number in selected) and wanted.issubset(mapping)
            accepted = bool(mask_pass and in_assembly)
            source_decisions.append(dict(copy_id=copy["copy_id"], source_label_chain=source["selected_label_asym_id"],
                                         assembly_memberships=copy["assembly_memberships"], in_selected_assembly=in_assembly,
                                         same_full_chemical_graph_as_primary=True, qualified_for_primary_mask=accepted,
                                         source_gate_residues=[dict(branch_num=residue["branch_num"], source_gate_pass=residue["source_gate_pass"], source_gate_failures=residue["source_gate_failures"]) for residue in source["residues"]],
                                         known_site_status="qualified_atomic_reference" if accepted else "outside_selected_assembly" if not in_assembly else "unscored_known_same_chemical_site_fails_fixed_mask",
                                         diagnostic_copy_id=copy["copy_id"]))
            if not accepted:
                continue
            for operator in operators:
                if source["selected_label_asym_id"] not in operator["subchains"]:
                    continue
                transformed = []
                for primary_atom in scored:
                    key = primary_atom["canonical_atom_key"]
                    atom = dict(mapping[key], canonical_atom_key=key)
                    atom["source_xyz_A"] = atom["xyz_A"]
                    atom["xyz_A"] = (np.array(operator["matrix"]) @ np.array(atom["xyz_A"]) + np.array(operator["vector"])).tolist()
                    transformed.append(atom)
                qualified.append(dict(copy_id=copy["copy_id"] + "@assembly" + context["selected_assembly_id"] + ":" + operator["operator_name"],
                                      source_copy_id=copy["copy_id"], assembly_id=context["selected_assembly_id"], operator_name=operator["operator_name"],
                                      matrix=operator["matrix"], vector=operator["vector"], scored_atom_count=len(transformed), scored_atom_map=transformed,
                                      coordinate_frame="Deposited biological-assembly receptor world frame; rigid source operator applied to complete native reference.",
                                      density_quality_inheritance="Same experimental ASU copy quality under rigid assembly operator; do not sample transformed coordinates against an untransformed crystal map."))
        _same(source_decisions, case["source_copies"], label + " fixed-mask source-copy decisions")
        _compare_replay(qualified, case["qualified_copies"], label + " transformed qualified copies", density_tolerance_sigma, deltas)
        total_scored += len(scored)
        total_qualified += len(qualified)
        reports.append(dict(case_id=label, qualification=scope, native_atom_rows=diagnostic["native_heavy_atom_row_count"], native_unique_heavy_atoms=diagnostic["native_unique_heavy_atom_count"],
                            scored_atoms=len(scored), source_copies=len(source_decisions), qualified_assembly_copies=len(qualified), status="verified"))
    _same(ledger["qualification_counts"], {scope: sum(c["qualification"] == scope for c in reports) for scope in ["full_ligand", "modeled_core", "rejected"]}, "Qualification counts")
    _same(ledger["source_preparation_eligible_count"], sum(c["source_conditioned_preparation_eligible"] for c in ledger["cases"]), "Source preparation count")
    _same(ledger["qualified_assembly_copy_count"], total_qualified, "Qualified biological-assembly copy count")
    return dict(schema_version=1, status="verified", source_only=True, prediction_outcomes_consulted=False,
                ledger_sha256=_sha(ledger_file), case_count=len(actual["cases"]),
                source_atom_row_count=sum(c["native_heavy_atom_row_count"] for c in actual["cases"]),
                source_unique_heavy_atom_count=sum(c["native_unique_heavy_atom_count"] for c in actual["cases"]),
                primary_scored_atom_count=total_scored, full_chemical_source_copy_count=len(actual_copy_entries), qualified_assembly_copy_count=total_qualified,
                checked_archive_count=len(manifest["archived_files"]), density_tolerance_sigma=density_tolerance_sigma,
                maximum_absolute_density_difference_sigma=max(deltas, default=0),
                tolerance_scope="Only normalized map samples. Native identities, coordinates, graphs, cutoff decisions, canonical masks and rigid assembly transforms are exact. Chirality compares determinant signs, not informative magnitudes.",
                manual_adjudications="Pinned judgments retained; replay does not independently prove paper/reagent/functional-site/assembly/family/training judgments.", cases=reports)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--density-tolerance-sigma", type=float, default=DEFAULT_DENSITY_TOLERANCE_SIGMA)
    parser.add_argument("--audit-only", action="store_true", help="Emit source diagnostics before a qualification ledger is frozen")
    args = parser.parse_args(argv)
    try:
        result = audit_all(root=args.root) if args.audit_only else verify_heldout_references(args.root, ledger_path=args.ledger, density_tolerance_sigma=args.density_tolerance_sigma)
    except (ValueError, KeyError, StopIteration, RuntimeError, OSError) as exc:
        print(json.dumps(dict(status="failed", error=str(exc))))
        return 1
    print(json.dumps(result, indent=2))
    return 0
