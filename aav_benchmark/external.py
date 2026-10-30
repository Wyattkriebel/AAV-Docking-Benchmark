"""Explicit-chemistry preparation of released protein–glycan development cases.

These source-reproduction cases exercise the implementation.  They are not an
independently curated, held-out estimate of biological docking reliability.
"""

from __future__ import annotations

import hashlib
import gzip
import json
import math
import pathlib
import urllib.request


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
COHORT_PATH = PROJECT_ROOT / "benchmarks" / "external_cohort.json"


def cohort():
    return json.loads(COHORT_PATH.read_text())


def case_definition(case_id):
    case_id = case_id.upper()
    matches = [row for row in cohort()["cases"] if row["case_id"] == case_id]
    if len(matches) != 1:
        raise ValueError(f"Unknown external development case {case_id!r}")
    return matches[0]


def _case_sources(case_id):
    manifest = cohort()
    definition = case_definition(case_id)
    components = sorted({row["comp_id"] for row in definition["residues"]})
    components += [row["ccd_id"] for row in definition.get("receptor_cofactors", [])]
    return (definition["sources"]
            + [manifest["ccd_sources"][comp] for comp in components]
            + [manifest["license_source"]])


def verify_external(raw_dir, case_id="154L"):
    """Refuse changed or missing source bytes before any preparation."""
    raw_dir = pathlib.Path(raw_dir)
    for source in _case_sources(case_id):
        path = raw_dir / source["filename"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {path}")


def fetch_external(raw_dir, case_id="154L"):
    """Download only this case's small, checksum-pinned development sources."""
    raw_dir = pathlib.Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    for source in _case_sources(case_id):
        destination = raw_dir / source["filename"]
        if destination.exists():
            actual = hashlib.sha256(destination.read_bytes()).hexdigest()
            if actual != source["sha256"]:
                raise ValueError(f"Existing source checksum mismatch: {destination}")
            continue
        archived = PROJECT_ROOT / source.get("archive", "no-archive")
        if archived.is_file():
            payload = gzip.decompress(archived.read_bytes())
        else:
            with urllib.request.urlopen(source["url"], timeout=60) as response:
                payload = response.read()
        if hashlib.sha256(payload).hexdigest() != source["sha256"]:
            raise ValueError(f"Downloaded source checksum mismatch: {source['url']}")
        destination.write_bytes(payload)
    verify_external(raw_dir, case_id)
    return [source["filename"] for source in _case_sources(case_id)]


def pdb_heavy_atoms(text):
    """Read a single-model, unambiguous PDB; retain full residue identities."""
    result = []
    seen = set()
    models = 0
    for line in text.splitlines():
        if line.startswith("MODEL"):
            models += 1
            if models > 1:
                raise ValueError("Multiple PDB models require an explicit selection")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 78:
            raise ValueError("PDB atom line lacks an explicit element field")
        element = line[76:78].strip()
        if element in {"H", "D"}:
            continue
        if not element:
            raise ValueError("Missing PDB element")
        if line[16:17] not in {" ", ""}:
            raise ValueError("Alternate conformers require an explicit selection")
        identity = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(),
                    line[17:20].strip(), line[12:16].strip())
        if identity in seen:
            raise ValueError(f"Duplicate PDB atom identity: {identity}")
        seen.add(identity)
        coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError("Nonfinite PDB coordinates")
        result.append({"key": list(identity), "element": element,
                       "coordinates": coordinates})
    if not result:
        raise ValueError("PDB contains no heavy atoms")
    return result


def _ccd_component(path):
    import gemmi

    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = {}
    columns = ["atom_id", "type_symbol", "charge", "pdbx_stereo_config",
               "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal",
               "pdbx_model_Cartn_z_ideal"]
    for row in block.find("_chem_comp_atom.", columns):
        if row[1] in {"H", "D"}:
            continue
        atoms[row[0]] = {"element": row[1], "charge": int(row[2]),
                         "stereo": row[3], "ideal": list(map(float, list(row)[4:]))}
    bonds = []
    for row in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        if row[0] in atoms and row[1] in atoms:
            bonds.append((row[0], row[1], row[2]))
    if not atoms or (not bonds and len(atoms) > 1):
        raise ValueError(f"Empty CCD heavy-atom graph: {path}")
    return {"atoms": atoms, "bonds": bonds}


def build_ligand(raw_dir, case_id="154L", embedding_seed=20261024):
    """Build from CCD bonds and declared branch links, then embed independently.

    Reference coordinates are used only to verify chemical handedness against
    CCD ideal coordinates and assign tetrahedral identity.  Every reference
    conformer is then discarded before ETKDG embedding and MMFF optimization.
    Glycosidic bonds are never inferred from spatial proximity.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    raw_dir = pathlib.Path(raw_dir)
    definition = case_definition(case_id)
    reference = pdb_heavy_atoms((raw_dir / f"{case_id}_l_b.pdb").read_text())
    # Released inputs use one chain B and consecutive branch numbers.
    by_identity = {}
    for atom in reference:
        chain, residue, insertion, comp, name = atom["key"]
        if chain != "B" or insertion.strip():
            raise ValueError("External ligand residue mapping is not the declared B/branch map")
        by_identity[(residue, comp, name)] = atom
    components = {row["comp_id"]: _ccd_component(raw_dir / f"{row['comp_id']}.cif")
                  for row in definition["residues"]}
    residue_components = {row["resnum"]: row["comp_id"] for row in definition["residues"]}
    leaving = set()
    replacements = {}
    links = []
    for row in definition["links"]:
        if row["value_order"].lower() != "sing":
            raise ValueError("Unsupported glycosidic bond order")
        left = (int(row["entity_branch_list_num_1"]), row["comp_id_1"], row["atom_id_1"])
        right = (int(row["entity_branch_list_num_2"]), row["comp_id_2"], row["atom_id_2"])
        links.append((left, right))
        # Heavy leaving atoms must be absent; hydrogen removal follows valence.
        for key, other, leave_name in [(left, right, row["leaving_atom_id_1"]),
                                      (right, left, row["leaving_atom_id_2"])]:
            comp = components[key[1]]
            if leave_name in comp["atoms"]:
                leave_key = (key[0], key[1], leave_name)
                leaving.add(leave_key)
                replacements[leave_key] = other
    expected = {(resnum, comp, name) for resnum, comp in residue_components.items()
                for name in components[comp]["atoms"]} - leaving
    if expected != set(by_identity):
        raise ValueError(f"Reference atom mismatch: missing={sorted(expected-set(by_identity))}; "
                         f"unexpected={sorted(set(by_identity)-expected)}")
    # Verify each tetrahedral CCD stereocenter using named-neighbor handedness.
    stereo_checks = []
    for resnum, comp_name in residue_components.items():
        component = components[comp_name]
        for name, atom in component["atoms"].items():
            if atom["stereo"] not in {"R", "S"}:
                continue
            neighbors = sorted({b if a == name else a for a, b, _ in component["bonds"]
                                if a == name or b == name})
            if len(neighbors) not in {3, 4}:
                raise ValueError(f"Cannot validate tetrahedral CCD center {resnum}:{name}")
            ideal_center = np.asarray(atom["ideal"])
            observed_center = np.asarray(by_identity[(resnum, comp_name, name)]["coordinates"])
            ideal_vectors, observed_vectors = [], []
            for neighbor in neighbors:
                neighbor_key = (resnum, comp_name, neighbor)
                observed_key = replacements.get(neighbor_key, neighbor_key)
                ideal_vectors.append(np.asarray(component["atoms"][neighbor]["ideal"]) - ideal_center)
                observed_vectors.append(np.asarray(by_identity[observed_key]["coordinates"]) - observed_center)
            if len(neighbors) == 4:
                ideal_vectors = [v - ideal_vectors[3] for v in ideal_vectors[:3]]
                observed_vectors = [v - observed_vectors[3] for v in observed_vectors[:3]]
            expected_volume = float(np.linalg.det(ideal_vectors))
            observed_volume = float(np.linalg.det(observed_vectors))
            if abs(expected_volume) < 0.05 or abs(observed_volume) < 0.05:
                raise ValueError(f"Ambiguous/flat chemical stereocenter {resnum}:{name}")
            if expected_volume * observed_volume <= 0:
                raise ValueError(f"CCD/reference stereochemistry mismatch {resnum}:{name}")
            stereo_checks.append({"resnum": resnum, "comp_id": comp_name, "atom": name,
                                  "ccd_stereo": atom["stereo"], "handedness_verified": True})
    editable = Chem.RWMol()
    indices = {}
    atom_map = []
    for number, reference_atom in enumerate(reference):
        chain, resnum, insertion, comp, name = reference_atom["key"]
        key = (resnum, comp, name)
        template = components[comp]["atoms"][name]
        if reference_atom["element"].upper() != template["element"].upper():
            raise ValueError(f"CCD/reference element mismatch at {key}")
        atom = Chem.Atom(template["element"])
        atom.SetFormalCharge(template["charge"])
        unique_name = f"L{number:03d}"
        atom.SetProp("source_atom_id", name)
        atom.SetProp("benchmark_name", unique_name)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(unique_name.rjust(4))
        info.SetResidueName(comp)
        info.SetResidueNumber(resnum)
        info.SetChainId(chain)
        atom.SetMonomerInfo(info)
        indices[key] = editable.AddAtom(atom)
        atom_map.append({"name": unique_name, "source_key": reference_atom["key"],
                         "rdkit_index": number, "element": atom.GetSymbol()})
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for resnum, comp in residue_components.items():
        for left, right, order in components[comp]["bonds"]:
            key1, key2 = (resnum, comp, left), (resnum, comp, right)
            if key1 in indices and key2 in indices:
                editable.AddBond(indices[key1], indices[key2], orders[order])
    for left, right in links:
        editable.AddBond(indices[left], indices[right], Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    Chem.SanitizeMol(molecule)
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Declared ligand graph is disconnected")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, atom in enumerate(reference):
        conformer.SetAtomPosition(index, atom["coordinates"])
    molecule.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(molecule, replaceExistingTags=True)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(molecule, includeUnassigned=True, useLegacyImplementation=False)
    if any(stereo == "?" for _, stereo in centers):
        raise ValueError("Unassigned ligand stereochemistry")
    smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    molecule.RemoveAllConformers()
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = embedding_seed
    parameters.numThreads = 1
    embedding_attempts = [{"algorithm": "ETKDGv3", "seed": embedding_seed,
                           "use_random_coordinates": False}]
    embedded = AllChem.EmbedMolecule(molecule, parameters)
    embedding_attempts[-1]["return_code"] = embedded
    if embedded != 0:
        # A fixed second development attempt, recorded even when successful.
        molecule.RemoveAllConformers()
        parameters.useRandomCoords = True
        parameters.maxIterations = 1000
        embedding_attempts.append({"algorithm": "ETKDGv3", "seed": embedding_seed,
                                   "use_random_coordinates": True, "max_iterations": 1000})
        embedded = AllChem.EmbedMolecule(molecule, parameters)
        embedding_attempts[-1]["return_code"] = embedded
        if embedded != 0:
            raise ValueError(f"Independent ligand embedding failed: {embedding_attempts}")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError("Unsupported MMFF ligand chemistry")
    optimization = AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000)
    if optimization != 0:
        raise ValueError(f"Independent ligand MMFF optimization failed ({optimization})")
    return molecule, reference, atom_map, {
        "isomeric_smiles": smiles, "heavy_atoms": len(reference),
        "formal_charge": Chem.GetFormalCharge(molecule),
        "embedding_seed": embedding_seed, "embedding": "ETKDGv3; one thread; no reference coordinates",
        "embedding_attempts": embedding_attempts,
        "optimization": "MMFF94; maximum 1000 iterations; convergence required",
        "stereochemistry": "CCD ideal named-neighbor handedness checked against reference; reference coordinates discarded before embedding",
        "stereocenter_checks": stereo_checks, "declared_glycosidic_links": definition["links"],
        "protonation": "CCD formal charges, RDKit valence hydrogens; not an assay-specific pH model",
    }


def _receptor_numbering_for_preparation(pdb_text, definition):
    """Expose independently mapped missing residues without inventing atoms."""
    mapping_path = definition.get("receptor_numbering_map")
    if not mapping_path:
        return pdb_text, {}, None
    path = PROJECT_ROOT / mapping_path
    if hashlib.sha256(path.read_bytes()).hexdigest() != definition["receptor_numbering_map_sha256"]:
        raise ValueError("Receptor author-numbering map checksum mismatch")
    mapping = json.loads(path.read_text())
    forward, reverse = {}, {}
    for row in mapping["residues"]:
        source = (row["prepared_chain"], row["prepared_residue_number"], "", row["prepared_resname"])
        author = (row["deposited_author_chain"], row["deposited_author_number"],
                  row["deposited_insertion_code"], row["deposited_resname"])
        if source[3] != author[3] or source in forward or author in reverse:
            raise ValueError("Nonbijective or chemically changing receptor numbering map")
        forward[source] = author
        reverse[author] = source
    observed = {tuple(atom["key"][:4]) for atom in pdb_heavy_atoms(pdb_text)}
    if observed != set(forward):
        raise ValueError("Receptor residues do not match the frozen author-numbering map")
    transformed = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM", "TER   ")):
            source = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())
            if source in forward:
                chain, number, insertion, residue = forward[source]
                line = line[:17]+f"{residue:3s}"+line[20:21]+f"{chain:1s}{number:4d}{insertion or ' ':1s}"+line[27:]
        transformed.append(line)
    return "\n".join(transformed)+"\n", reverse, mapping


def _restore_prepared_source_numbering(text, reverse):
    if not reverse:
        return text
    result = []
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            author = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())
            if author not in reverse:
                raise ValueError("Prepared protein atom lacks a source numbering identity")
            chain, number, insertion, residue = reverse[author]
            line = line[:17]+f"{residue:3s}"+line[20:21]+f"{chain:1s}{number:4d}{insertion or ' ':1s}"+line[27:]
        result.append(line)
    return "\n".join(result)+"\n"


def prepare_external_receptor(pdb_text, output_stem, definition, raw_dir):
    """Keep curated calcium ions as standard fixed [Ca+2] Vina atoms.

    Protein preparation remains the same audited Meeko Polymer procedure.
    Calcium is noncovalent, parameterized independently, and retained at its
    exact native coordinate and source identity. No protein atom is repaired.
    """
    import numpy as np
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from .aav import prepare_receptor, _pdb_atoms
    from .provenance import write_json

    output_stem = pathlib.Path(output_stem)
    definitions = definition.get("receptor_cofactors", [])
    if not definitions:
        return prepare_receptor(pdb_text, output_stem)
    if len(definitions) != 1 or definitions[0]["ccd_id"] != "CA":
        raise ValueError("Only explicitly curated native calcium is supported")
    ion_definition = definitions[0]
    ccd = _ccd_component(pathlib.Path(raw_dir) / "CA.cif")
    if len(ccd["atoms"]) != 1:
        raise ValueError("Calcium CCD is not a single ion")
    ccd_atom = next(iter(ccd["atoms"].values()))
    if ccd_atom["element"].upper() != "CA" or ccd_atom["charge"] != 2 or ccd["bonds"]:
        raise ValueError("Curated calcium does not match CCD [Ca+2]")
    original = pdb_heavy_atoms(pdb_text)
    ions = [atom for atom in original if atom["key"][3] == ion_definition["source_resname"]]
    if len(ions) != ion_definition["expected_count"]:
        raise ValueError("Required calcium count differs from the frozen source")
    ion_keys = {tuple(atom["key"]) for atom in ions}
    if any(atom["element"].upper() != "CA" or atom["key"][4] != "CA+2" for atom in ions):
        raise ValueError("Unexpected required-calcium chemical identity")
    protein_lines = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(),
                   line[17:20].strip(), line[12:16].strip())
            if key in ion_keys:
                continue
        protein_lines.append(line)
    protein_text, reverse_numbering, numbering_audit = _receptor_numbering_for_preparation(
        "\n".join(protein_lines)+"\n", definition)
    protein_log = prepare_receptor(protein_text, output_stem)
    receptor_pdbqt = _restore_prepared_source_numbering(
        output_stem.with_suffix(".pdbqt").read_text(), reverse_numbering)
    serials = [int(line[6:11]) for line in receptor_pdbqt.splitlines()
               if line.startswith(("ATOM  ", "HETATM"))]
    next_serial = max(serials) + 1
    ion_logs = []
    appended = []
    for ion in ions:
        molecule = Chem.MolFromSmiles("[Ca+2]")
        conformer = Chem.Conformer(1)
        conformer.SetAtomPosition(0, ion["coordinates"])
        molecule.AddConformer(conformer)
        chain, number, insertion, residue, name = ion["key"]
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4))
        info.SetResidueName(residue)
        info.SetResidueNumber(number)
        info.SetChainId(chain)
        info.SetInsertionCode(insertion)
        molecule.GetAtomWithIdx(0).SetMonomerInfo(info)
        setup = MoleculePreparation().prepare(molecule)
        if len(setup) != 1:
            raise ValueError("Expected one standalone calcium setup")
        text, ok, error = PDBQTWriterLegacy.write_string(setup[0])
        if not ok:
            raise ValueError(error)
        lines = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
        if len(lines) != 1 or float(lines[0][70:76]) != 2 or lines[0][77:].strip() != "Ca":
            raise ValueError("Meeko did not produce standard fixed +2 calcium parameters")
        line = lines[0][:6] + f"{next_serial:5d}" + lines[0][11:]
        next_serial += 1
        appended.append(line)
        ion_logs.append({"source_atom_key": ion["key"], "coordinates_A": ion["coordinates"],
                         "formal_charge": 2, "pdbqt_charge": float(line[70:76]),
                         "autodock_atom_type": "Ca", "ccd_id": "CA", "source_identity_preserved": True,
                         "parameterization": "Meeko MoleculePreparation on standalone [Ca+2]"})
    combined = receptor_pdbqt.rstrip()+"\n"+"\n".join(appended)+"\n"
    dock_atoms = _pdb_atoms(combined, pdbqt=True)
    before = _pdb_atoms(pdb_text)
    before_heavy = {key for key, atom in before.items() if not atom["hydrogen"]}
    dock_heavy = {key for key, atom in dock_atoms.items() if not atom["hydrogen"]}
    if before_heavy != dock_heavy:
        raise ValueError("Protein/calcium preparation changed the full native heavy-atom identity set")
    displacement = max(float(np.linalg.norm(np.asarray(before[key]["coordinates"])-dock_atoms[key]["coordinates"]))
                       for key in before_heavy)
    if displacement > 0.002:
        raise ValueError("Protein/calcium preparation moved native heavy atoms")
    output_stem.with_suffix(".pdb").write_text(pdb_text)
    output_stem.with_suffix(".pdbqt").write_text(combined)
    log = {**protein_log, "preparation": "Audited Meeko protein Polymer plus native standalone [Ca+2] cofactor",
           "input_atoms": protein_log["input_atoms"]+len(ions),
           "prepared_atoms": protein_log["prepared_atoms"]+len(ions),
           "pdbqt_atoms": protein_log["pdbqt_atoms"]+len(ions),
           "input_heavy_atoms": len(before_heavy), "prepared_heavy_atoms": len(before_heavy),
           "pdbqt_heavy_atoms": len(dock_heavy), "max_heavy_displacement_A": displacement,
           "native_cofactors": ion_logs,
           "cofactor_limitation": "Fixed +2 Vina calcium model; no specialized coordination or polarization; not assay-validated"}
    if numbering_audit:
        log["receptor_numbering_audit"] = numbering_audit
        log["protein_preparation_identity_frame"] = "deposited author numbering; output PDBQT restored to upstream prepared identifiers"
        log["unresolved_region_policy"] = "No loop atoms or cross-gap peptide bond added; existing Meeko residue templates; native heavy atoms unchanged"
    write_json(output_stem.with_name(output_stem.name+"_atom_changes.json"), log)
    return log


def _prepare(raw_dir, output_dir, case_id="154L"):
    """Prepare a released bound example in the generic runner's schema."""
    import numpy as np
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    case_id = case_id.upper()
    raw_dir, output_dir = pathlib.Path(raw_dir), pathlib.Path(output_dir)
    definition = case_definition(case_id)
    verify_external(raw_dir, case_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    molecule, reference, atom_map, chemistry = build_ligand(raw_dir, case_id)
    (output_dir / "chemistry.json").write_text(json.dumps(chemistry, indent=2) + "\n")
    reference_coordinates = np.asarray([atom["coordinates"] for atom in reference])
    center = reference_coordinates.mean(axis=0)
    # Translation contains site information in this explicitly oracle-only arm.
    positions = molecule.GetConformer().GetPositions()
    positions += center - positions[:len(reference)].mean(axis=0)
    for index, point in enumerate(positions):
        molecule.GetConformer().SetAtomPosition(index, point)
    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1:
        raise ValueError("Expected one Meeko ligand setup")
    pdbqt_text, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(error)
    prepared_names = [line[12:16].strip() for line in pdbqt_text.splitlines()
                      if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip().startswith("L")]
    expected_names = [atom["name"] for atom in atom_map]
    if len(prepared_names) != len(expected_names) or set(prepared_names) != set(expected_names):
        raise ValueError("Meeko changed or lost a ligand heavy-atom identity")
    (output_dir / "ligand.pdbqt").write_text(pdbqt_text)
    writer = Chem.SDWriter(str(output_dir / "ligand.sdf"))
    writer.write(molecule)
    writer.close()
    receptor_text = (raw_dir / f"{case_id}_r_b.pdb").read_text()
    receptor_atoms = pdb_heavy_atoms(receptor_text)
    preparation_changes = prepare_external_receptor(receptor_text, output_dir / "holo_complete", definition, raw_dir)
    box_size = np.maximum(20.0, np.ptp(reference_coordinates, axis=0) + 16.0)
    if np.any(box_size > definition.get("reviewed_max_box_side_A", 40) + 1e-8):
        raise ValueError("Oracle box exceeds the protocol's 40 Å review limit")
    reference_data = {
        "names": expected_names, "coordinates": reference_coordinates.tolist(),
        "atom_map": atom_map,
        "receptor_atom_keys": [atom["key"] for atom in receptor_atoms],
        "receptor_coordinates": [atom["coordinates"] for atom in receptor_atoms],
        "receptor_residue_keys": [atom["key"][:4] for atom in receptor_atoms],
        "coordinate_frame": "released_bound_protein_ligand_pair",
    }
    (output_dir / "reference.json").write_text(json.dumps(reference_data, indent=2) + "\n")
    metadata = {
        "case_id": case_id, "family_id": definition["family_id"],
        "status": "prepared", "cohort_role": "development_only",
        "information_regime": "true_site_local; not blind",
        "reference_quality": definition["reference_quality"],
        "source_commit": cohort()["source_commit"],
        "box_center": center.tolist(), "box_size": box_size.tolist(),
        "box_review_reason": definition.get("box_review_reason"),
        "ligand_file": "ligand.pdbqt", "reference_file": "reference.json",
        "receptors": [{"id": "holo_complete", "state": "holo",
                       "assembly_condition": "released_bound_receptor",
                       "pdb_file": "holo_complete.pdb", "pdbqt_file": "holo_complete.pdbqt"}],
        "chemistry": chemistry, "receptor_preparation": preparation_changes,
        "scope": "Reuses a released HADDOCK example with Vina; not a reproduction of published HADDOCK performance or a confirmatory cohort.",
    }
    (output_dir / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def prepare(raw_dir, output_dir, case_id="154L"):
    """Prepare one case and preserve a visible record if any stage fails."""
    from .provenance import write_json

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # A failed rerun cannot leave a previous success manifest authoritative.
    (output_dir / "preparation.json").unlink(missing_ok=True)
    attempt = {"case_id": case_id.upper(), "status": "in_progress"}
    write_json(output_dir / "preparation_attempt.json", attempt)
    try:
        result = _prepare(raw_dir, output_dir, case_id)
    except Exception as error:
        write_json(output_dir / "preparation_attempt.json",
                   {**attempt, "status": "failed", "error_type": type(error).__name__, "reason": str(error)})
        raise
    write_json(output_dir / "preparation_attempt.json", {**attempt, "status": "ok"})
    return result


def prepare_development(raw_dir, output_root):
    """Retain every frozen case, including chemistry/preparation failures."""
    from .provenance import write_json

    output_root = pathlib.Path(output_root)
    statuses = []
    for definition in cohort()["cases"]:
        case_id = definition["case_id"]
        row = {"case_id": case_id, "family_id": definition["family_id"],
               "reference_quality": definition["reference_quality"]}
        try:
            fetch_external(raw_dir, case_id)
            metadata = prepare(raw_dir, output_root / case_id, case_id)
            statuses.append({**row, "status": "ok", "prepared_dir": str(output_root / case_id),
                             "heavy_atoms": metadata["chemistry"]["heavy_atoms"]})
        except Exception as error:
            unsupported = "uncurated receptor residue" in str(error) and ":CA2" in str(error)
            category = ("unsupported_required_calcium_cofactor" if unsupported else
                        "chemistry_reference_mismatch" if "stereochemistry mismatch" in str(error) else
                        "preparation_failure")
            statuses.append({**row, "status": "unsupported" if unsupported else "failed",
                             "failure_category": category,
                             "error_type": type(error).__name__, "reason": str(error)})
        write_json(output_root / "development_preparation.json", {
            "cohort_role": "implementation_development_only", "eligible_cases": len(cohort()["cases"]),
            "completed_preparation_attempts": len(statuses), "cases": statuses})
    return statuses


def _ccd_ligand(raw_dir, component_ids, branch_links, embedding_seed=20261024):
    """Assemble CCD stereochemical graphs without experimental coordinates.

    Copy the tetrahedral orientation from each ideal CCD monomer. Replacing a
    leaving oxygen by a glycosidic oxygen preserves its neighbor-order parity;
    this avoids treating absolute CIP labels as invariant after substitution.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    raw_dir = pathlib.Path(raw_dir)
    components = [_ccd_component(raw_dir / f"{name}.cif") for name in component_ids]
    complete = Chem.RWMol()
    full_indices, full_coordinates = {}, []
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for resnum, component in enumerate(components, 1):
        for name, definition in component["atoms"].items():
            atom = Chem.Atom(definition["element"])
            atom.SetFormalCharge(definition["charge"])
            full_indices[(resnum, name)] = complete.AddAtom(atom)
            full_coordinates.append(definition["ideal"])
        for left, right, order in component["bonds"]:
            complete.AddBond(full_indices[(resnum, left)], full_indices[(resnum, right)], orders[order])
    complete = complete.GetMol()
    Chem.SanitizeMol(complete)
    conf = Chem.Conformer(complete.GetNumAtoms())
    for index, position in enumerate(full_coordinates):
        conf.SetAtomPosition(index, position)
    complete.AddConformer(conf)
    Chem.AssignAtomChiralTagsFromStructure(complete, replaceExistingTags=True)
    index_keys = {index: key for key, index in full_indices.items()}
    replacements = {}
    links = []
    for link in branch_links:
        if link["bond_order"] != "SING":
            raise ValueError("Unsupported assay glycosidic bond order")
        donor = (link["donor_residue"], link["donor_atom"])
        acceptor = (link["acceptor_residue"], link["acceptor_atom"])
        links.append((donor, acceptor))
        leaving = (donor[0], link["donor_leaving_atom"])
        if leaving not in full_indices:
            raise ValueError(f"CCD donor leaving atom missing: {leaving}")
        replacements[leaving] = acceptor
    editable = Chem.RWMol()
    indices, atom_map = {}, []
    for key, full_index in full_indices.items():
        if key in replacements:
            continue
        atom = Chem.Atom(complete.GetAtomWithIdx(full_index))
        number = editable.GetNumAtoms()
        name = f"L{number:03d}"
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4))
        info.SetResidueName(component_ids[key[0]-1])
        info.SetResidueNumber(key[0])
        info.SetChainId("B")
        atom.SetMonomerInfo(info)
        indices[key] = editable.AddAtom(atom)
        atom_map.append({"name": name, "source_key": ["B", key[0], "", component_ids[key[0]-1], key[1]],
                         "rdkit_index": number, "element": atom.GetSymbol()})
    for bond in complete.GetBonds():
        left, right = index_keys[bond.GetBeginAtomIdx()], index_keys[bond.GetEndAtomIdx()]
        if left in indices and right in indices:
            editable.AddBond(indices[left], indices[right], bond.GetBondType())
    for left, right in links:
        editable.AddBond(indices[left], indices[right], Chem.BondType.SINGLE)
    molecule = editable.GetMol()
    final_keys = {index: key for key, index in indices.items()}
    for key, index in indices.items():
        original = complete.GetAtomWithIdx(full_indices[key])
        if original.GetChiralTag() not in {Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                                          Chem.ChiralType.CHI_TETRAHEDRAL_CW}:
            continue
        desired = [replacements.get(index_keys[a.GetIdx()], index_keys[a.GetIdx()])
                   for a in original.GetNeighbors()]
        actual = [final_keys[a.GetIdx()] for a in molecule.GetAtomWithIdx(index).GetNeighbors()]
        if set(desired) != set(actual):
            raise ValueError("CCD linkage changed a stereocenter's neighbor identity unexpectedly")
        permutation = [desired.index(neighbor) for neighbor in actual]
        inversions = sum(permutation[i] > permutation[j] for i in range(len(permutation))
                         for j in range(i+1, len(permutation)))
        if inversions % 2:
            molecule.GetAtomWithIdx(index).InvertChirality()
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Assay ligand graph is disconnected")
    if any(stereo == "?" for _, stereo in Chem.FindMolChiralCenters(molecule, includeUnassigned=True,
                                                                  useLegacyImplementation=False)):
        raise ValueError("Unassigned CCD assay ligand stereochemistry")
    smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    molecule.RemoveAllConformers()
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = embedding_seed
    parameters.numThreads = 1
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise ValueError("Independent assay ligand embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError("Unsupported assay ligand MMFF chemistry")
    if AllChem.MMFFOptimizeMolecule(molecule, maxIters=1000) != 0:
        raise ValueError("Independent assay ligand MMFF optimization did not converge")
    return molecule, atom_map, {"isomeric_smiles": smiles, "heavy_atoms": len(atom_map),
                               "formal_charge": Chem.GetFormalCharge(molecule),
                               "embedding_seed": embedding_seed,
                               "stereochemistry": "CCD ideal stereochemistry; explicit neighbor-order parity preserved across branch links",
                               "protonation": "CCD formal charges; RDKit valence hydrogens"}


def fetch_gal3_panel(raw_dir):
    """Fetch the frozen exploratory panel's CCD monomer sources."""
    raw_dir = pathlib.Path(raw_dir)
    fetch_external(raw_dir, "1KJL")
    manifest = cohort()
    for name in ["NDG", "BGC", "GLA", "GLC"]:
        source = manifest["ccd_sources"][name]
        path = raw_dir / source["filename"]
        if path.exists():
            payload = path.read_bytes()
        else:
            archived = PROJECT_ROOT / source.get("archive", "no-archive")
            if archived.is_file():
                payload = gzip.decompress(archived.read_bytes())
            else:
                with urllib.request.urlopen(source["url"], timeout=60) as response:
                    payload = response.read()
        if hashlib.sha256(payload).hexdigest() != source["sha256"]:
            raise ValueError(f"Assay CCD checksum mismatch: {name}")
        path.write_bytes(payload)


def prepare_gal3_panel(raw_dir, receptor_dir, output_root):
    """Prepare six ligand identities/forms for descriptive affinity ranking.

    The same LacNAc-informed box and isolated Gal3 CRD receptor are supplied to
    every alternative. Generated coordinates are explicitly ineligible as
    atomic reference truth.  Free-sugar anomer populations are not assumed.
    """
    import numpy as np
    import shutil
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    raw_dir, receptor_dir, output_root = map(pathlib.Path, (raw_dir, receptor_dir, output_root))
    fetch_gal3_panel(raw_dir)
    metadata = json.loads((receptor_dir / "preparation.json").read_text())
    if metadata["case_id"] != "1KJL" or len(metadata["receptors"]) != 1:
        raise ValueError("Gal3 panel requires the frozen 1KJL bound receptor")
    panel = json.loads((PROJECT_ROOT / "benchmarks" / "galectin3_itc_panel.json").read_text())
    frozen = cohort()["gal3_panel_frozen_analysis"]
    prepared = []
    for row in panel["rows"]:
        for anomer, component_ids in frozen["reducing_end_forms"][row["ligand_id"]].items():
            molecule, atom_map, chemistry = _ccd_ligand(raw_dir, component_ids, row["branch_links"])
            if chemistry["heavy_atoms"] != row["expected_heavy_atoms"] or chemistry["formal_charge"] != row["formal_charge"]:
                raise ValueError("Assay ligand graph does not match frozen measured chemical identity")
            positions = molecule.GetConformer().GetPositions()
            positions += np.asarray(metadata["box_center"]) - positions[:len(atom_map)].mean(axis=0)
            for index, point in enumerate(positions):
                molecule.GetConformer().SetAtomPosition(index, point)
            setup = MoleculePreparation().prepare(molecule)
            if len(setup) != 1:
                raise ValueError("Expected one assay ligand Meeko setup")
            pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup[0])
            if not ok:
                raise ValueError(error)
            names = [atom["name"] for atom in atom_map]
            written = [line[12:16].strip() for line in pdbqt.splitlines()
                       if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip().startswith("L")]
            if len(written) != len(names) or set(written) != set(names):
                raise ValueError("Assay ligand atom loss in PDBQT preparation")
            case_id = f"gal3_{row['ligand_id']}_{anomer}"
            output_dir = output_root / case_id
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "ligand.pdbqt").write_text(pdbqt)
            writer = Chem.SDWriter(str(output_dir / "ligand.sdf"))
            writer.write(molecule)
            writer.close()
            for filename in ["holo_complete.pdb", "holo_complete.pdbqt", "holo_complete_atom_changes.json"]:
                shutil.copyfile(receptor_dir / filename, output_dir / filename)
            reference = {"names": names, "coordinates": positions[:len(names)].tolist(),
                         "atom_map": atom_map, "pose_reference_eligible": False,
                         "coordinate_source": "independently generated input; no experimental pose truth",
                         "scope": "Used only for atom identity verification; RMSD against this input is not a scientific endpoint"}
            (output_dir / "reference.json").write_text(json.dumps(reference, indent=2) + "\n")
            result = {**{key: metadata[key] for key in ["box_center", "box_size", "receptors", "ligand_file", "reference_file"]},
                      "case_id": case_id, "family_id": "galectin_3", "status": "prepared",
                      "cohort_role": "exploratory_assay_development", "information_regime": "1KJL_LacNAc_oracle_box_shared_by_all_alternatives",
                      "assay_endpoint": "within_assay_rank_only", "pose_reference_eligible": False,
                      "panel_id": panel["panel_id"], "ligand_id": row["ligand_id"], "anomer": anomer,
                      "assay_kd_uM": row["kd_uM"], "assay_kd_sem_uM": row["kd_sem_uM"],
                      "chemistry": chemistry, "frozen_analysis": frozen,
                      "construct_mismatch": panel["target"]["construct_mismatch"],
                      "scope": "All compounds are measured binders; no tested-negative classification; three identities on one target; reducing-end equilibrium populations unknown"}
            (output_dir / "preparation.json").write_text(json.dumps(result, indent=2) + "\n")
            prepared.append(result)
    return prepared
