"""Independently audit raw POL evidence before descriptive aggregation.

Archive relocation is supported. Commands retain their historical paths; their
case/artifact identities and all actual input bytes are verified. The complete
96-run denominator is required for every aggregate, including the primary.
"""
from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import tempfile

import numpy as np

from .aav import _pdb_atoms
from .engine import read_poses
from .external import _ccd_component, pdb_heavy_atoms
from .pol_specificity import (COHORT_ID, METAL_PROFILES, MODEL_SPECS, SEEDS,
                              load_source_qualification, planned_cases,
                              receptor_sources, summarize_pol_scores,
                              validate_experiment)
from .provenance import sha256, write_text_atomic
from .specificity import _check_stereo


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _name(name):
    _require(isinstance(name, str) and Path(name).name == name and name not in {"", ".", ".."},
             "Artifact filename must be simple and local")
    return name


def _same_derived(left, right, tolerance=1e-9):
    """Allow source-only linear-algebra roundoff; raw artifacts stay exact."""
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(_same_derived(left[k], right[k], tolerance) for k in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(_same_derived(a, b, tolerance) for a, b in zip(left, right))
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isfinite(left) and math.isfinite(right) and math.isclose(left, right, rel_tol=0, abs_tol=tolerance)
    return type(left) is type(right) and left == right


class _Evidence:
    def __init__(self):
        self.sources = {}

    def read(self, path, key):
        path = Path(path)
        if not path.is_file() and path.with_name(path.name + ".gz").is_file():
            path = path.with_name(path.name + ".gz")
        raw = path.read_bytes()
        content = gzip.decompress(raw) if path.suffix == ".gz" else raw
        row = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
        if path.suffix == ".gz":
            row.update({"gzip_sha256": hashlib.sha256(raw).hexdigest(), "gzip_bytes": len(raw)})
        self.sources[key] = row
        return content

    def json(self, path, key):
        return json.loads(self.read(path, key))

    def digest(self, path, key):
        self.read(path, key)
        return self.sources[key]["sha256"]


def _prepared_files(evidence, root, row):
    case = row["case_id"]
    candidates = [root / row["prepared_dir"], root / case, root / case / "inputs"]
    direct = next((p for p in candidates if (p / "preparation.json").is_file()), None)
    expected = row["prepared_artifact_sha256"]
    _require(isinstance(expected, dict) and expected, f"Missing prepared artifact pins: {case}")
    if direct is not None:
        content = {name: evidence.read(direct / _name(name), f"prepared/{case}/{name}") for name in expected}
    else:
        directory = root / case / "inputs"
        manifest = evidence.json(directory / "manifest.json", f"prepared/{case}/archive_manifest.json")
        entries = manifest["files"]
        names = [r["filename"] for r in entries]
        archives = [r["archive"] for r in entries]
        _require(len(names) == len(set(names)) and len(archives) == len(set(archives))
                 and set(names) == set(expected), f"Prepared archive denominator/duplicates differ: {case}")
        content = {}
        for entry in entries:
            name, archive = _name(entry["filename"]), _name(entry["archive"])
            compressed = directory / archive
            _require(sha256(compressed) == entry["archive_sha256"], f"Prepared archive bytes changed: {case}/{name}")
            raw = evidence.read(compressed, f"prepared/{case}/{name}")
            _require(hashlib.sha256(raw).hexdigest() == entry["sha256"]
                     and len(raw) == entry["original_bytes"], f"Prepared archive contents changed: {case}/{name}")
            content[name] = raw
    for name, digest in expected.items():
        _require(hashlib.sha256(content[name]).hexdigest() == digest,
                 f"Prepared artifact hash changed: {case}/{name}")
    _require(hashlib.sha256(content["preparation.json"]).hexdigest() == row["preparation_sha256"],
             f"Prepared metadata hash changed: {case}")
    return content


def _canonical(molecule):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)


def _atom_type_element(atom_type):
    types = {"C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N",
             "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S", "Ca": "CA",
             "H": "H", "HD": "H", "HS": "H"}
    _require(atom_type in types, f"Unexpected PDBQT atom type: {atom_type}")
    return types[atom_type]


def _ligand(content, case, graph, source_dir):
    from rdkit import Chem

    reference = json.loads(content["reference.json"])
    metadata = json.loads(content["preparation.json"])
    chemistry = metadata["chemistry"]
    source = json.loads(content["source_graph.json"])
    names, atom_map = reference["names"], reference["atom_map"]
    ccd = _ccd_component(source_dir / f"{case['ccd_id']}.cif")
    expected_keys = [["B", 1, "", case["ccd_id"], name] for name in ccd["atoms"]]
    _require(len(names) == 12 and len(set(names)) == 12
             and names == [r["name"] for r in atom_map]
             and [r["rdkit_index"] for r in atom_map] == list(range(12))
             and [r["source_key"] for r in atom_map] == expected_keys,
             f"Incomplete or incorrect CCD heavy atom map: {case['case_id']}")
    _require(reference.get("pose_reference_eligible") is False
             and chemistry.get("pose_reference_eligible") is False,
             "Generated POL input cannot become atomic pose truth")
    _require(source.get("source_candidate_id") == COHORT_ID and source.get("ccd_graph") == graph
             and source.get("atom_map") == atom_map
             and source.get("source_ledger_sha256") == metadata.get("source_ledger_sha256"),
             f"Prepared source graph identity changed: {case['case_id']}")
    _require(chemistry == json.loads(content["chemistry.json"])
             and chemistry.get("heavy_atoms") == 12 and chemistry.get("formal_charge") == 0
             and chemistry.get("ccd_id") == case["ccd_id"]
             and chemistry.get("model_id") == case["model_id"]
             and chemistry.get("species_id") == case["species_id"]
             and chemistry.get("anomer") == case["anomer"]
             and chemistry.get("embedding_seed") == 20261024
             and chemistry.get("native_coordinates_used_for_embedding") is False
             and chemistry.get("stereochemistry_verified_from_3d") is True,
             f"Prepared chemistry provenance changed: {case['case_id']}")
    molecule = Chem.MolFromMolBlock(content["ligand.sdf"].decode(), removeHs=False, sanitize=True)
    _require(molecule is not None and _canonical(molecule) == graph["canonical_smiles"]
             and chemistry["isomeric_smiles"] == graph["canonical_smiles"],
             f"Actual SDF chemical graph/stereochemistry differs: {case['case_id']}")
    _check_stereo(molecule, graph["canonical_smiles"])
    elements = {}
    for row, ccd_atom in zip(atom_map, ccd["atoms"].values()):
        element = molecule.GetAtomWithIdx(row["rdkit_index"]).GetSymbol().upper()
        _require(element == row["element"].upper() == ccd_atom["element"].upper(),
                 "POL atom map elements differ from CCD/SDF")
        elements[row["name"]] = element
    _require(Counter(elements.values()) == Counter({"C": 6, "O": 6}), "POL ligand element denominator changed")
    xyz = np.asarray(reference["coordinates"], float)
    _require(xyz.shape == (12, 3) and np.isfinite(xyz).all(), "Invalid POL generated input coordinates")
    _require(np.allclose(xyz.mean(0), metadata["box_center"], rtol=0, atol=1e-9),
             "POL generated ligand centroid differs from its frozen condition box")
    sdf_xyz = molecule.GetConformer().GetPositions()[[r["rdkit_index"] for r in atom_map]]
    _require(np.max(np.abs(sdf_xyz - xyz)) <= 0.000051, "POL SDF/reference input coordinates differ")
    pdbqt = content["ligand.pdbqt"].decode()
    smiles = [line[len("REMARK SMILES "):].strip() for line in pdbqt.splitlines()
              if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX")]
    _require(len(smiles) == 1 and Chem.MolFromSmiles(smiles[0]) is not None
             and _canonical(Chem.MolFromSmiles(smiles[0])) == graph["canonical_smiles"],
             "POL actual PDBQT chemical graph differs from the source anomer")
    coords = {}
    for line in pdbqt.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = _atom_type_element(line[77:].strip())
        if element == "H":
            continue
        name = line[12:16].strip()
        _require(name not in coords and elements.get(name) == element, "POL PDBQT heavy map/type changed")
        coords[name] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    _require(set(coords) == set(names), "POL PDBQT lost heavy atoms")
    _require(np.max(np.abs(np.array([coords[n] for n in names]) - xyz)) <= 0.000501,
             "POL actual PDBQT coordinates differ from generated input")
    return names, elements, molecule, atom_map, metadata


def _receptor(content, descriptor, expected_text, state, profile):
    original = pdb_heavy_atoms(expected_text)
    actual = pdb_heavy_atoms(content[descriptor["pdb_file"]].decode())
    _require(actual == original, f"POL actual receptor differs from native complete assembly: {state}/{profile}")
    text = content[descriptor["pdbqt_file"]].decode()
    atoms = _pdb_atoms(text, pdbqt=True)
    heavy = {key: row for key, row in atoms.items() if not row["hydrogen"]}
    before = {tuple(row["key"]): row for row in original}
    _require(set(heavy) == set(before), f"POL PDBQT receptor heavy atom identity changed: {state}/{profile}")
    for key, row in heavy.items():
        _require(_atom_type_element(row["element_or_type"]) == before[key]["element"].upper(),
                 "POL receptor PDBQT element changed")
        _require(np.linalg.norm(np.array(row["coordinates"]) - before[key]["coordinates"]) <= 0.002,
                 "POL receptor PDBQT moved native heavy atoms")
    ions = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CA2"]
    expected_count = 4 if profile == "Ca_retained" else 0
    _require(len(ions) == expected_count and all(line[12:16].strip() == "CA+2"
             and line[77:].strip() == "Ca" and float(line[70:76]) == 2 for line in ions),
             "POL native calcium count/charge/type changed")
    log = json.loads(content[descriptor["atom_changes_file"]])
    _require(log.get("status") == "ok" and log.get("heavy_atom_completion") is False
             and log.get("coordinate_minimization") is False
             and log.get("pdbqt_heavy_atoms") == len(before)
             and log.get("max_heavy_displacement_A", math.inf) <= 0.002,
             "POL receptor preparation audit changed")
    _require(len(log.get("native_cofactors", [])) == expected_count, "POL calcium preparation audit denominator changed")


def _command(row, case, metadata, descriptor):
    command = row.get("command")
    _require(isinstance(command, list) and len(command) > 1 and all(isinstance(v, str) for v in command)
             and len(command[1:]) % 2 == 0, "Missing/malformed actual POL command")
    _require(Path(command[0]).name == "vina-1.2.7", "POL command executable identity changed")
    options = dict(zip(command[1::2], command[2::2]))
    _require(len(options) * 2 == len(command) - 1, "Duplicate actual POL command flag")
    numeric = {"--seed": row["seed"], "--cpu": 2, "--exhaustiveness": 16,
               "--num_modes": 20, "--min_rmsd": 0.5, "--energy_range": 10}
    for axis, center, size in zip("xyz", metadata["box_center"], metadata["box_size"]):
        numeric["--center_" + axis], numeric["--size_" + axis] = center, size
    _require(set(options) == set(numeric) | {"--scoring", "--receptor", "--ligand", "--out"},
             "Unexpected/missing actual POL command options")
    for flag, expected in numeric.items():
        value = float(options[flag])
        _require(math.isfinite(value) and value == expected, f"Actual POL command changed {flag}")
    _require(options["--scoring"] == "vina", "POL scoring function changed")
    for flag, filename in (("--receptor", descriptor["pdbqt_file"]), ("--ligand", "ligand.pdbqt"),
                           ("--out", row["pose_file"])):
        _require(Path(options[flag]).parts[-2:] == (case, filename), f"POL command {flag} changed case/artifact")


def _raw_poses(raw, names, elements, molecule, atom_map, expected_smiles):
    from rdkit import Chem

    model_count, element_counts = 0, []
    current = None
    for line in raw.decode().splitlines():
        if line.startswith("MODEL"):
            current = Counter()
        elif line.startswith(("ATOM  ", "HETATM")):
            element = _atom_type_element(line[77:].strip())
            if element != "H":
                _require(current is not None and elements.get(line[12:16].strip()) == element,
                         "Raw POL pose heavy element/map differs")
                current[element] += 1
        elif line.startswith("ENDMDL"):
            _require(current == Counter({"C": 6, "O": 6}), "Raw POL pose lost a ligand heavy element")
            element_counts.append(dict(current))
            model_count += 1
            current = None
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "poses.pdbqt"
        path.write_bytes(raw)
        poses = read_poses(path, names)
    _require(len(poses) == model_count, "POL raw MODEL denominator changed")
    for pose in poses:
        derived = Chem.Mol(molecule)
        for row, point in zip(atom_map, pose["coordinates"]):
            derived.GetConformer().SetAtomPosition(row["rdkit_index"], point)
        _check_stereo(derived, expected_smiles)
    return poses


def _log(raw, row, metadata, poses):
    text = raw.decode()
    fields = ((r"^\s*AutoDock Vina\s+(\S+)\s*$", "v1.2.7", "version"),
              (r"^\s*Scoring function\s*:\s*(\S+)\s*$", "vina", "scoring function"),
              (r"^\s*Exhaustiveness\s*:\s*(\S+)\s*$", "16", "exhaustiveness"),
              (r"^\s*CPU\s*:\s*(\S+)\s*$", "2", "CPU"),
              (r"^\s*Performing docking \(random seed:\s*(\d+)\)\s*\.\.\.\s*$", str(row["seed"]), "random seed"))
    for pattern, expected, label in fields:
        _require(re.findall(pattern, text, re.MULTILINE) == [expected],
                 f"Raw POL log does not certify exact {label}")
    for label, filename in (("Rigid receptor", f"{metadata['state']}_{row['receptor']}.pdbqt"),
                            ("Ligand", "ligand.pdbqt")):
        paths = re.findall(r"^\s*" + label + r"\s*:\s*(.+?)\s*$", text, re.MULTILINE)
        _require(len(paths) == 1 and Path(paths[0]).parts[-2:] == (metadata["case_id"], filename),
                 f"Raw POL log changed {label} input identity")
    for label, values in (("Grid center", metadata["box_center"]), ("Grid size", metadata["box_size"])):
        found = re.findall(r"^\s*" + label + r"\s*:\s*X\s+(\S+)\s+Y\s+(\S+)\s+Z\s+(\S+)\s*$", text, re.MULTILINE)
        _require(len(found) == 1 and all(math.isclose(float(a), b, rel_tol=0, abs_tol=0.0001)
                                       for a, b in zip(found[0], values)), f"Raw POL log changed {label}")
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    energies = [(int(rank), float(score)) for rank, score in re.findall(
        r"^\s*(\d+)\s+(" + number + r")\s+\S+\s+\S+\s*$", text, re.MULTILINE)]
    _require(1 <= len(energies) <= 20 and [rank for rank, _ in energies] == list(range(1, len(energies) + 1))
             and all(math.isfinite(score) for _, score in energies), "Raw POL log mode denominator/ranks differ")
    _require(all(a[1] <= b[1] + 0.002 for a, b in zip(energies, energies[1:])), "Raw POL log modes are not score ordered")
    # Vina 1.2.7 global_search logs all modes at four significant digits;
    # get_poses writes only the energy-window prefix at three decimal places.
    # https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp
    def log_roundoff(score):
        return 0.5 * 10 ** (math.floor(math.log10(abs(score))) - 3) if score else 0.0005
    _require(len(poses) <= len(energies), "Raw POL pose count exceeds its log table")
    for pose, (rank, score) in zip(poses, energies):
        _require(pose["rank"] == rank
                 and abs(pose["vina_score"] - score) <= 0.0005 + log_roundoff(score) + 1e-12,
                 "Raw POL log/pose score tables differ")
    best = energies[0][1]
    for index, (_, score) in enumerate(energies):
        uncertainty = log_roundoff(best) + log_roundoff(score) + 1e-12
        if index < len(poses):
            _require(score - best <= 10 + uncertainty, "Raw POL pose exceeds the frozen output energy window")
        else:
            _require(score - best >= 10 - uncertainty, "Raw POL omitted logged mode is inside the frozen output energy window")
    return {"logged_modes": len(energies), "energy_window_omitted_modes": len(energies) - len(poses),
            "omitted_log_scores_kcal_mol": [score for _, score in energies[len(poses):]],
            "log_pose_score_comparison": "four-significant-digit log / three-decimal-place PDBQT rounding compatibility; complete retained rank prefix"}


def generate_pol_report(root, *, prepared_root=None, results_root=None):
    root = Path(root).resolve()
    results = Path(results_root).resolve() if results_root else root / "results/pol_specificity"
    prepared = Path(prepared_root).resolve() if prepared_root else results
    archive = root / "data/archive/independent_specificity_sources"
    ledger = root / "benchmarks/independent_specificity_candidates.json"
    config_path = root / "benchmarks/pol_experiment.json"
    evidence = _Evidence()
    config = evidence.json(config_path, "benchmarks/pol_experiment.json")
    config_hash = sha256(config_path)
    _require(evidence.digest(results / "pol_experiment.json", "results/pol_experiment.json") == config_hash,
             "Portable POL experiment copy changed")
    validate_experiment(config, ledger, archive, execution=True)
    source, graphs = load_source_qualification(ledger, archive)
    evidence.read(ledger, "benchmarks/independent_specificity_candidates.json")
    manifest = evidence.json(archive / "manifest.json", "source/manifest.json")
    for row in manifest["source_files"] + manifest["derived_files"]:
        _require(evidence.digest(archive / _name(row["path"]), "source/" + row["path"]) == row["sha256"],
                 "POL source archive content changed")
    _require(evidence.digest(archive / "CA.cif", "source/CA.cif") == config["calcium_ccd_sha256"], "Calcium CCD changed")
    calcium = evidence.json(archive / "calcium_source.json", "source/calcium_source.json")
    _require(calcium.get("filename") == "CA.cif" and calcium.get("sha256") == config["calcium_ccd_sha256"]
             and calcium.get("bytes") == (archive / "CA.cif").stat().st_size,
             "Separate calcium source pin differs")
    freeze = evidence.json(results / "freeze_receipt.json", "results/freeze_receipt.json")
    proposal = evidence.json(results / "proposal_before_freeze.json", "results/proposal_before_freeze.json")
    expected_proposal = dict(config, analysis_frozen_before_new_docking=False)
    _require(proposal == expected_proposal and sha256(results / "proposal_before_freeze.json") == freeze["proposal_sha256"],
             "Original POL proposal/freeze transition changed")
    _require(freeze.get("configuration_sha256") == config_hash
             and freeze.get("analysis_frozen_before_new_docking") is True
             and freeze.get("implementation_files_sha256") == config["implementation_files_sha256"]
             and freeze.get("source_ledger_sha256") == config["source_ledger_sha256"]
             and freeze.get("source_manifest_sha256") == config["source_manifest_sha256"]
             and freeze.get("new_scores_consulted") is False
             and freeze.get("experimental_species") == 2
             and freeze.get("independent_assay_contexts") == 1
             and freeze.get("expected_runs") == 96, "POL pre-run freeze receipt changed")
    snapshot = evidence.json(results / "execution_snapshot.json", "results/execution_snapshot.json")
    snapshot_rows = snapshot["files"]
    expected_files = set(config["implementation_files_sha256"]) | {"aav_benchmark/archive.py", "scripts/run_pol_specificity.py"}
    _require(snapshot.get("configuration_sha256") == config_hash and snapshot.get("workers") == 1
             and snapshot.get("cpu_per_run") == 2 and len(snapshot_rows) == len(expected_files)
             and {r["path"] for r in snapshot_rows} == expected_files, "POL execution source snapshot denominator changed")
    for row in snapshot_rows:
        relative = Path(row["snapshot_path"])
        _require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe execution snapshot path")
        digest = evidence.digest(results / relative, "snapshot/" + row["path"])
        _require(digest == row["sha256"] and (row["path"] not in config["implementation_files_sha256"]
                 or digest == config["implementation_files_sha256"][row["path"]]), "POL frozen execution code snapshot changed")
    _require(evidence.digest(results / "execution_snapshot/execution_driver.py", "snapshot/execution_driver.py") == snapshot["driver_sha256"],
             "POL execution driver snapshot changed")
    lock = evidence.json(root / "benchmarks/vina.lock.json", "benchmarks/vina.lock.json")
    approved_binaries = {r["sha256"] for r in lock["assets"].values()}
    _require(snapshot["executable_sha256"] in approved_binaries, "POL executable is not an official pinned Vina asset")
    panel = evidence.json(results / "panel_preparation.json", "results/panel_preparation.json")
    cases = planned_cases(config)
    canonical = {r["case_id"]: r for r in cases}
    rows = panel["cases"]
    _require(panel.get("status") == "prepared" and panel.get("cohort_id") == COHORT_ID
             and panel.get("experiment_sha256") == config_hash and panel.get("expected_runs") == 96
             and panel.get("source_ledger_sha256") == config["source_ledger_sha256"]
             and len(rows) == 16 and {r["case_id"] for r in rows} == set(canonical), "POL prepared case denominator/freeze differs")
    for row in rows:
        _require(all(row.get(k) == v for k, v in canonical[row["case_id"]].items())
                 and row.get("prepared_dir") == "cases/" + row["case_id"], "POL portable prepared case identity differs")
    unexpected = [p.name for p in results.iterdir() if p.is_dir() and (p / "runs.json").is_file() and p.name not in canonical]
    _require(not unexpected, "Unexpected POL run-bearing case directory")
    native, centers, alignment = receptor_sources(archive)
    _require(_same_derived(panel.get("alignment"), alignment), "POL prepared protein-only alignment changed")
    for state in ("holo", "apo"):
        _require(panel["source_receptor_audits"][state] == native[state][1], "POL assembly/alternate-conformer source audit changed")
    shared_receptors, outputs, scores = {}, [], []
    for row in rows:
        case_id, case = row["case_id"], canonical[row["case_id"]]
        files = _prepared_files(evidence, prepared, row)
        names, elements, molecule, atom_map, metadata = _ligand(files, case, graphs[case["ccd_id"]], archive)
        _require(all(metadata.get(k) == v for k, v in case.items())
                 and metadata.get("cohort_id") == COHORT_ID and metadata.get("family_id") == "POL_CBM67"
                 and metadata.get("status") == "prepared"
                 and metadata.get("experiment_sha256") == config_hash
                 and metadata.get("source_ledger_sha256") == config["source_ledger_sha256"]
                 and _same_derived(metadata.get("box_center"), centers[case["state"]][case["site"]])
                 and metadata.get("box_size") == [24.0] * 3
                 and metadata.get("execution_information_regime") == "true-site-local"
                 and metadata.get("experimental_anomer_label_available") is False
                 and metadata.get("operational_glucose_NB_is_numeric_Kd_bound") is False,
                 "POL prepared condition/box/source identity differs")
        receptors = metadata["receptors"]
        _require(len(receptors) == 2 and [r["id"] for r in receptors] == list(METAL_PROFILES), "POL prepared metal profile denominator changed")
        descriptors = {}
        for descriptor in receptors:
            profile, state = descriptor["id"], case["state"]
            stem = f"{state}_{profile}"
            _require(descriptor == {"id": profile, "state": state, "assembly_condition": "complete_biological_dimer_" + profile,
                     "pdb_file": stem + ".pdb", "pdbqt_file": stem + ".pdbqt", "atom_changes_file": stem + "_atom_changes.json"},
                     "POL prepared receptor descriptor differs")
            expected_text = native[state][0]
            if profile == "dry":
                expected_text = "\n".join(line for line in expected_text.splitlines()
                    if not (line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "CA2")) + "\n"
            _receptor(files, descriptor, expected_text, state, profile)
            hashes = tuple(hashlib.sha256(files[descriptor[k]]).hexdigest() for k in ("pdb_file", "pdbqt_file", "atom_changes_file"))
            _require((state, profile) not in shared_receptors or shared_receptors[(state, profile)] == hashes,
                     "POL ligands/sites did not share identical receptor preparation")
            shared_receptors[(state, profile)] = hashes
            descriptors[profile] = descriptor
        directory = results / case_id
        records = evidence.json(directory / "runs.json", f"results/{case_id}/runs.json") if (directory / "runs.json").is_file() else []
        _require(isinstance(records, list), "POL run manifest must be a list")
        identities = [(r.get("receptor"), r.get("seed")) for r in records]
        _require(len(identities) == len(set(identities)), "Duplicate POL run identity")
        _require(set(identities) <= set((p, s) for p in METAL_PROFILES for s in SEEDS), "Unexpected POL run identity")
        if records:
            _require(evidence.digest(directory / "configuration.json", f"results/{case_id}/configuration.json") == config_hash,
                     "POL case frozen configuration changed")
        lookup = {(r["receptor"], r["seed"]): r for r in records}
        for profile in METAL_PROFILES:
            descriptor = descriptors[profile]
            expected_fp = {"config_sha256": config_hash, "executable_sha256": snapshot["executable_sha256"],
                           "ligand_sha256": hashlib.sha256(files["ligand.pdbqt"]).hexdigest(),
                           "preparation_sha256": hashlib.sha256(files["preparation.json"]).hexdigest(),
                           "reference_sha256": hashlib.sha256(files["reference.json"]).hexdigest(),
                           "receptor_sha256": hashlib.sha256(files[descriptor["pdbqt_file"]]).hexdigest(),
                           "receptor_pdb_sha256": hashlib.sha256(files[descriptor["pdb_file"]]).hexdigest(),
                           "full_assembly_sha256": None}
            for seed in SEEDS:
                record = lookup.get((profile, seed))
                output = {**case, "metal_profile": profile, "seed": seed, "status": "missing", "top1_score": None}
                if record is not None:
                    _require(record.get("case_id") == case_id and record.get("family_id") == "POL_CBM67"
                             and record.get("receptor_state") == case["state"]
                             and record.get("assembly_condition") == "complete_biological_dimer_" + profile
                             and record.get("information_regime") == "true-site-local"
                             and record.get("vina_version") == "AutoDock Vina v1.2.7", "POL run identity/version/regime differs")
                    _require(record.get("fingerprint") == expected_fp, "POL run frozen input fingerprint changed")
                    pose_name, log_name = f"{profile}_{seed}.pdbqt", f"{profile}_{seed}.log"
                    _require(record.get("pose_file") == pose_name and record.get("log_file") == log_name, "POL run artifact filenames changed")
                    _command(record, case_id, metadata, descriptor)
                    status = record.get("status")
                    _require(status in {"ok", "failed", "running"}, "Unknown POL run status")
                    output.update({"status": status, "command": record["command"], "input_fingerprint": expected_fp,
                                   "error": record.get("error")})
                    for field in ("wall_seconds", "cpu_seconds", "child_peak_rss_bytes"):
                        if field in record:
                            value = record[field]
                            _require(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value >= 0,
                                     "Invalid POL resource measurement")
                            output[field] = value
                    if status == "ok":
                        _require(record.get("exit_code") == 0, "Successful POL run has nonzero exit status")
                        raw = evidence.read(directory / pose_name, f"results/{case_id}/{pose_name}")
                        pose_hash = hashlib.sha256(raw).hexdigest()
                        _require(pose_hash == record.get("pose_sha256"), "Raw POL pose checksum changed")
                        poses = _raw_poses(raw, names, elements, molecule, atom_map, graphs[case["ccd_id"]]["canonical_smiles"])
                        _require(len(poses) == record.get("returned_poses") and 1 <= len(poses) <= 20, "POL returned pose denominator changed")
                        log = evidence.read(directory / log_name, f"results/{case_id}/{log_name}")
                        log_audit = _log(log, record, metadata, poses)
                        output.update({"top1_score": poses[0]["vina_score"], "pose_sha256": pose_hash,
                                       "log_sha256": hashlib.sha256(log).hexdigest(), "returned_poses": len(poses)})
                        output.update(log_audit)
                        scores.append({"case_id": case_id, "metal_profile": profile, "seed": seed, "top1_score": poses[0]["vina_score"]})
                outputs.append(output)
    counts = Counter(row["status"] for row in outputs)
    eligible = counts["ok"] == 96
    aggregation = summarize_pol_scores(scores) if eligible else None
    summary = {"schema_version": 1, "cohort_id": COHORT_ID,
               "status": "evaluated_descriptive" if eligible else "ineligible_incomplete_frozen_run_denominator",
               "information_regime": config["information_regime"], "configuration_sha256": config_hash,
               "report_implementation_sha256": sha256(Path(__file__)),
               "denominator": {"expected_runs": 96, "successful_runs": counts["ok"], "failed_runs": counts["failed"],
                               "running_runs": counts["running"], "missing_runs": counts["missing"],
                               "cases": 16, "model_states": 4, "experimental_species": 2,
                               "independent_families": 1, "experimental_epimer_contrasts": 1,
                               "primary_runs": 12, "sensitivity_runs": 84, "metrics_eligible": eligible},
               "metrics_eligible": eligible, "aggregation": aggregation, "runs": outputs,
               "verification_geometry_roundoff_tolerance_A": 1e-9,
               "resources": {"wall_seconds_sum": math.fsum(r.get("wall_seconds", 0) for r in outputs),
                             "cpu_seconds_sum": math.fsum(r.get("cpu_seconds", 0) for r in outputs),
                             "returned_poses": sum(r.get("returned_poses", 0) for r in outputs),
                             "logged_modes": sum(r.get("logged_modes", 0) for r in outputs),
                             "energy_window_omitted_modes": sum(r.get("energy_window_omitted_modes", 0) for r in outputs),
                             "wall_measurements_present": sum("wall_seconds" in r for r in outputs),
                             "cpu_measurements_present": sum("cpu_seconds" in r for r in outputs)},
               "receptor_audits": {state: native[state][1] for state in ("holo", "apo")}, "alignment": panel["alignment"],
               "assay": source["assay"], "experimental_anomer_fractions": None,
               "source_hash_manifest": {"verified_artifacts": dict(sorted(evidence.sources.items())),
                                        "article_sha256_declared": source["primary_source"]["article_pdf_sha256"],
                                        "publication_bytes_locally_reverified": False},
               "limitations": ["One natural epimer contrast in one additional protein family; no precision or reliability certification.",
                               "No AUROC/AUPRC is estimated from one pair and no score-to-dissociation calibration is performed.",
                               "Glucose NB is absence of detectable binding in the source ITC titration, not a numerical Kd bound.",
                               "Four model states represent two free-sugar species; equal anomer weighting is a heuristic and solution populations/furanose/open-chain fractions are unknown.",
                               "Known-site boxes use prior melibiose structural information; no blind-site claim follows.",
                               "Native glycosylated apo versus recombinant holo is a disclosed construct/context approximation; removed covalent glycans and Ser218 alternate selection are recorded.",
                               "Fixed +2 Vina calcium is not an assay-validated model of coordination or polarization; both profiles omit explicit waters.",
                               "Generated input conformers provide chemical identity only; no experimental pose truth or RMSD endpoint is available.",
                               "Exact assay-batch sequence, residual assay salt, ITC detection threshold and source Table SIII uncertainties remain unresolved."]}
    return summary, render_pol_markdown(summary)


def _table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
                     + ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def render_pol_markdown(summary):
    lines = ["# POL: an independent natural epimer contrast", "",
             "This frozen known-site analysis compares measured D-galactose binding with actually tested D-glucose showing no detectable ITC binding. NB is categorical; no numerical Kd bound is assigned. This is one contrast in the POL/CBM67 family.", ""]
    lines += ["The two sugars differ only at C4 and share formula C6H12O6, 12 heavy atoms and neutral model charge.", ""]
    denominator = summary["denominator"]
    lines += [f"Raw verification: **{denominator['successful_runs']}/96 successful runs**; {denominator['failed_runs']} failed, {denominator['running_runs']} running and {denominator['missing_runs']} missing.", ""]
    if not summary["metrics_eligible"]:
        lines += ["All aggregates are ineligible until the complete frozen 96-run denominator passes verification.", ""]
    else:
        aggregate = summary["aggregation"]
        primary = aggregate["primary"]
        gal = primary["species"]["D_galactose"]["equal_anomer_mean_kcal_mol"]
        glc = primary["species"]["D_glucose"]["equal_anomer_mean_kcal_mol"]
        outcome = "D-galactose scores lower" if gal < glc else "D-glucose scores lower" if glc < gal else "The aggregate scores tie"
        lines += [f"**Primary N-site, holo, calcium-retained condition: {outcome}.** Equal-anomer means are {gal:.4f} kcal/mol for D-galactose and {glc:.4f} kcal/mol for D-glucose. The glucose-minus-galactose difference is {glc-gal:.4f} kcal/mol.", "",
                  "Seeds are averaged within each modeled anomer, then the two anomers receive equal weight. This is a computational heuristic, not measured population weighting. Furanose, open-chain and solution anomer fractions are unknown. Sites, receptor states, calcium profiles, anomers and seeds are sensitivity conditions rather than independent biological replicates.", "",
                  "## Separate conditions", ""]
        sensitivities = [r for r in aggregate["conditions"] if r is not primary
                         and (r["state"], r["site"], r["metal_profile"]) != ("holo", "N", "Ca_retained")]
        if glc != gal:
            reversals = sum(r["glucose_minus_galactose_kcal_mol"] * (glc - gal) < 0 for r in sensitivities)
            direction_note = f"The ranking direction reverses in **{reversals} of {len(sensitivities)} prespecified sensitivity conditions**. These conditions are descriptive checks of this one contrast, not independent biological replications."
            insert_at = lines.index("## Separate conditions")
            lines[insert_at:insert_at] = [direction_note, ""]
        anomer_means = {species: {r["anomer"]: r["seed_mean_kcal_mol"] for r in primary["species"][species]["anomers"]}
                        for species in ("D_galactose", "D_glucose")}
        differences = [anomer_means["D_glucose"][a] - anomer_means["D_galactose"][a] for a in ("alpha", "beta")]
        if differences[0] * differences[1] < 0:
            insert_at = lines.index("## Separate conditions")
            lines[insert_at:insert_at] = ["The matched alpha and beta model pairs also favor different sugars in the primary condition. The equal-weight primary ranking is sensitive to this unmeasured solution-mixture choice.", ""]
        rows = []
        for condition in aggregate["conditions"]:
            g, c = condition["species"]["D_galactose"], condition["species"]["D_glucose"]
            rows.append([condition["state"], condition["site"], condition["metal_profile"],
                         f"{g['equal_anomer_mean_kcal_mol']:.4f}", f"{c['equal_anomer_mean_kcal_mol']:.4f}",
                         f"{condition['glucose_minus_galactose_kcal_mol']:.4f}",
                         f"{g['minimum_anomer_seed_mean_kcal_mol']:.4f} / {g['maximum_anomer_seed_mean_kcal_mol']:.4f}",
                         f"{c['minimum_anomer_seed_mean_kcal_mol']:.4f} / {c['maximum_anomer_seed_mean_kcal_mol']:.4f}"])
        lines += [_table(["State", "Site", "Calcium profile", "Gal mean", "Glc mean", "Glc-Gal", "Gal min/max anomer", "Glc min/max anomer"], rows), "",
                  "All energies are kcal/mol; lower scores are favorable. `dry` removes calcium; both profiles omit explicit waters.", "",
                  "## Every anomer and seed", ""]
        rows = []
        for condition in aggregate["conditions"]:
            for species in ("D_galactose", "D_glucose"):
                for anomer in condition["species"][species]["anomers"]:
                    rows.append([condition["state"], condition["site"], condition["metal_profile"], species,
                                 anomer["anomer"], ", ".join(f"{v:.3f}" for v in anomer["seed_top1_scores_kcal_mol"]),
                                 f"{anomer['seed_mean_kcal_mol']:.4f}"])
        lines += [_table(["State", "Site", "Profile", "Species", "Anomer", "Seeds 104729 / 130363 / 155921", "Mean"], rows), ""]
    resources = summary["resources"]
    lines += ["## Verification and scope", "",
              f"The audit checks frozen source/configuration/code snapshots, every preparation archive, exact CCD/SDF/PDBQT stereochemistry and atom maps, native dimer/calcium identities, each command, all raw pose hashes, every returned pose's heavy elements and stereochemistry, and raw log score tables. Recorded totals are {resources['returned_poses']} returned poses, {resources['wall_seconds_sum']:.1f} summed run wall-seconds and {resources['cpu_seconds_sum']:.1f} CPU-seconds.", "",
              f"The logs contain {resources['logged_modes']} modes; {resources['energy_window_omitted_modes']} are outside the saved energy window. Vina prints all search-result modes before its PDBQT writer applies the frozen 10 kcal/mol window. The audit requires the saved poses to be the corresponding rank prefix, allowing only documented print rounding. [Vina 1.2.7 source](https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp#L557-L613).", "",
              "The holo A/B model preserves 5,034 observed protein heavy atoms; the independently solved apo dimer preserves 5,058 after coherent highest-occupancy Ser218B selection. Calcium-retained profiles add four native ions. A 2,704-atom protein-backbone alignment transfers centers into the apo frame. Native apo glycosylation, missing N-terminal residues and assay/construct correspondence limitations remain disclosed.", "",
              "No AUROC/AUPRC or calibrated affinity is reported from one epimer pair. Fixed +2 Vina calcium approximates the native ions and does not establish accurate coordination energetics. Input coordinates are independently generated and supply no experimental pose truth.", "",
              "Assay source: [Perduca et al. 2020, Figure 4 and ITC methods](https://doi.org/10.1093/glycob/cwaa006). [Source qualification](INDEPENDENT_SPECIFICITY_CANDIDATES.md), [frozen experiment](../benchmarks/pol_experiment.json), and [machine-readable evidence](../results/pol_specificity/summary.json).", ""]
    return "\n".join(lines)


def write_pol_report(root, *, prepared_root=None, results_root=None, summary_output=None, markdown_output=None, check=False):
    root = Path(root).resolve()
    summary, markdown = generate_pol_report(root, prepared_root=prepared_root, results_root=results_root)
    json_path = Path(summary_output) if summary_output else root / "results/pol_specificity/summary.json"
    markdown_path = Path(markdown_output) if markdown_output else root / "docs/POL_SPECIFICITY_RESULTS.md"
    text = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if check:
        _require(json_path.read_text() == text and markdown_path.read_text() == markdown, "POL derived report is stale")
    else:
        write_text_atomic(json_path, text)
        write_text_atomic(markdown_path, markdown)
    return json_path, markdown_path
