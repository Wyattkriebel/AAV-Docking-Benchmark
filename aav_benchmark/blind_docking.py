"""Reference-independent predicted-box inputs and paired Vina execution.

This module never reads native ligand coordinate files or evaluation tables.
The oracle control takes a center explicitly supplied by its separate caller.
Archived ligand conformers were independently embedded: their old translation
is removed before they are placed at each new search center.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import shutil
import tempfile

import gemmi
import numpy as np

from .aav import PROTEIN_RESIDUES, _pdb_atoms, prepare_receptor, select_neighborhood
from .provenance import sha256, write_json, write_text_atomic
from .site_geometry import distinct_candidates, read_receptor
from .site_prediction_verification import restore_site_inputs


_LIGAND_DIRECTORIES = {"aav9": "aav9", "154L": "154l", "1GNY": "1GNY", "1KJL": "1KJL",
                       "1QFO": "1QFO", "2J1T": "2J1T", "1PWB": "1PWB", "2VUZ": "2VUZ"}
_PREDICTORS = {"fpocket": "fpocket_v2/execution.json", "p2rank": "p2rank/result.json"}


def _finite_vector(value, label, *, positive=False):
    xyz = np.asarray(value, dtype=float)
    if xyz.shape != (3,) or not np.isfinite(xyz).all() or (positive and np.any(xyz <= 0)):
        raise ValueError(f"{label} must be a finite {'positive ' if positive else ''}three-vector")
    return xyz


def _fresh(path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Preserve existing experiment inputs; use a fresh directory: {path}")
    path.mkdir(parents=True)
    return path


def _archived_file(directory, manifest, filename):
    matches = [entry for entry in manifest["files"] if entry["filename"] == filename]
    if len(matches) != 1:
        raise ValueError(f"Missing or ambiguous archived input: {filename}")
    entry = matches[0]
    if any(Path(entry[k]).name != entry[k] or entry[k] in {".", ".."} for k in ("filename", "archive")):
        raise ValueError("Archived input names must be plain filenames")
    archive = directory / entry["archive"]
    if archive.is_symlink() or sha256(archive) != entry["archive_sha256"]:
        raise ValueError("Archived ligand input bytes changed")
    content = gzip.decompress(archive.read_bytes())
    import hashlib
    if hashlib.sha256(content).hexdigest() != entry["sha256"]:
        raise ValueError("Restored ligand input bytes changed")
    return content, {"filename": filename, "sha256": entry["sha256"],
                     "archive_sha256": entry["archive_sha256"]}


def ligand_atoms(text):
    """Read independent input PDBQT atoms without accepting docking MODEL blocks."""
    if any(line.startswith(("MODEL", "REMARK VINA RESULT:")) for line in text.splitlines()):
        raise ValueError("A returned docking pose cannot be used as an independent ligand input")
    atoms = []
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip()
        xyz = _finite_vector([float(line[30:38]), float(line[38:46]), float(line[46:54])], "Ligand coordinate")
        atom_type = line[77:].strip()
        if not atom_type:
            raise ValueError("Ligand input atom has no Vina type")
        atoms.append({"name": name, "coordinates": xyz.tolist(), "hydrogen": atom_type in {"H", "HD", "HS"}})
    names = [a["name"] for a in atoms if not a["hydrogen"]]
    if not names or len(set(names)) != len(names):
        raise ValueError("Independent ligand must have complete unique heavy-atom names")
    return atoms


def ligand_source(root, complex_id):
    """Load only archived independent PDBQT and preparation provenance.

    The old reference JSON (which contains native coordinates) is deliberately
    neither opened nor restored. Atom names come from the input PDBQT itself.
    """
    root = Path(root)
    if complex_id not in _LIGAND_DIRECTORIES:
        raise ValueError("Unknown frozen ligand graph case")
    directory = root / "results/development" / _LIGAND_DIRECTORIES[complex_id] / "inputs"
    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"No independently validated ligand preparation is archived for {complex_id}")
    manifest = json.loads(manifest_path.read_text())
    raw, prep_receipt = _archived_file(directory, manifest, "preparation.json")
    preparation = json.loads(raw)
    chemistry = preparation.get("ligand", preparation.get("chemistry", {}))
    independent = (chemistry.get("native_coordinates_used_for_embedding") is False
                   or (chemistry.get("embedding") == "ETKDGv3; one thread; no reference coordinates"
                       and chemistry.get("stereochemistry") == "CCD ideal named-neighbor handedness checked against reference; reference coordinates discarded before embedding"))
    if preparation.get("status", "prepared") != "prepared" or not independent:
        raise ValueError("Ligand preparation does not establish independent embedding")
    filename = preparation.get("ligand_file", "galactose.pdbqt")
    raw, ligand_receipt = _archived_file(directory, manifest, filename)
    text = raw.decode("utf-8")
    atoms = ligand_atoms(text)
    xyz = np.asarray([a["coordinates"] for a in atoms if not a["hydrogen"]])
    radius = float(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1).max())
    side = int(max(24, math.ceil(2 * radius + 16)))
    if side > 42:
        raise ValueError("Independent conformer requires a search side beyond the frozen 42 A review limit")
    return {"complex_id": complex_id, "pdbqt": text, "names": [a["name"] for a in atoms if not a["hydrogen"]],
            "box_size": [side] * 3, "heavy_atom_centroid_radius_A": radius,
            "box_rule": "cube side ceil(max(24,2*max heavy-atom distance from independent-conformer centroid+16)) A; maximum42 A",
            "native_coordinates_used_for_embedding": False, "source_receipts": [prep_receipt, ligand_receipt],
            "archive_manifest": str(manifest_path.relative_to(root)), "archive_manifest_sha256": sha256(manifest_path),
            "chemistry": chemistry}


def translate_ligand(text, center):
    """Replace only absolute translation; preserve graph, torsions and orientation."""
    center = _finite_vector(center, "Search center")
    atoms = ligand_atoms(text)
    centroid = np.mean([a["coordinates"] for a in atoms if not a["hydrogen"]], axis=0)
    shift = center - centroid
    lines = []
    index = 0
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            point = np.asarray(atoms[index]["coordinates"]) + shift
            index += 1
            fields = [f"{x:8.3f}" for x in point]
            if any(len(field) != 8 for field in fields):
                raise ValueError("Translated PDBQT coordinate exceeds fixed-width representation")
            line = line[:30] + "".join(fields) + line[54:]
        lines.append(line)
    result = "\n".join(lines) + "\n"
    after = ligand_atoms(result)
    if [(a["name"], a["hydrogen"]) for a in atoms] != [(a["name"], a["hydrogen"]) for a in after]:
        raise AssertionError("Translation changed ligand atom identities")
    residual = np.asarray([a["coordinates"] for a in after]) - np.asarray([a["coordinates"] for a in atoms]) - shift
    if float(np.abs(residual).max()) > 0.000501:
        raise ValueError("Translation moved internal ligand coordinates beyond PDBQT rounding")
    return result, {"source_heavy_centroid_A": centroid.tolist(), "translation_A": shift.tolist(),
                    "old_translation_removed": True, "native_coordinate_access": False,
                    "max_coordinate_rounding_A": float(np.abs(residual).max())}


def frozen_candidates(root, condition_id, method, *, limit=5):
    """Read fixed predictor ranks and receptor-only symmetry, never site scoring."""
    root = Path(root)
    if method not in _PREDICTORS or limit != 5:
        raise ValueError("Use a frozen predictor and the declared five-box maximum")
    config_path = root / "benchmarks/site_experiment.json"
    config = json.loads(config_path.read_text())
    freeze = json.loads((root / "results/site_discovery/freeze_receipt.json").read_text())
    if sha256(config_path) != freeze["configuration_sha256"]:
        raise ValueError("Site configuration no longer matches its frozen receipt")
    cases = [case for case in config["cases"] if case["condition_id"] == condition_id]
    if len(cases) != 1:
        raise ValueError("Unknown site condition")
    case = cases[0]
    directory = root / "results/site_discovery" / condition_id
    path = directory / _PREDICTORS[method]
    record = json.loads(path.read_text())
    operations = []
    sources = {str(path.relative_to(root)): sha256(path)}
    if case["complex_id"] == "aav9":
        certificate = directory / "source_symmetry_amended.json"
        receipt = json.loads((root / "results/site_discovery/geometry_amendment_receipt.json").read_text())
        if sha256(certificate) != receipt["certificates"][condition_id]:
            raise ValueError("Source-only symmetry certificate changed")
        audit = json.loads(certificate.read_text())
        operations = [(np.asarray(o["rotation"]), np.asarray(o["translation"])) for o in audit["transforms"]]
        sources[str(certificate.relative_to(root))] = sha256(certificate)
    status = "ok" if record["status"] in {"ok", "completed"} else "predictor_failed"
    predictions = record.get("predictions", []) if status == "ok" else []
    candidates, mapping = distinct_candidates(predictions, operations, config["candidate_merge_distance_A"])
    return {"condition_id": condition_id, "complex_id": case["complex_id"], "state": case["state"], "method": method,
            "status": status if status != "ok" or candidates else "completed_no_candidates",
            "predictor_status": record["status"], "raw_candidate_count": len(predictions),
            "distinct_candidate_count": len(candidates), "candidates": candidates[:limit], "candidate_mapping": mapping,
            "prediction_input_sha256": record.get("superseded_input_sha256", record["input_sha256"]),
            "source_sha256": sources, "native_coordinate_access": False}


def _protein_text(receptor_cif, center, size, cutoff):
    keys, coordinates, elements = read_receptor(receptor_cif)
    if any(key[3] not in PROTEIN_RESIDUES for key in keys):
        raise ValueError("Docking input is not uniformly protein-only")
    atoms = dict(zip(keys, coordinates))
    selected, nearest = select_neighborhood(atoms, center, size, cutoff=cutoff)
    structure = gemmi.read_structure(str(receptor_cif))
    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    for chain in structure[0]:
        if chain.name in selected:
            if len(chain.name) != 1:
                raise ValueError("An explicit chain map is needed before PDB serialization")
            out[0].add_chain(chain.clone())
    selected_keys = {key for key in keys if key[0] in selected}
    if len(selected_keys) > 99999:
        raise ValueError("Whole-chain neighborhood exceeds PDB atom serial capacity")
    text = out.make_pdb_string()
    actual = {key: atom for key, atom in _pdb_atoms(text).items() if not atom["hydrogen"]}
    if set(actual) != selected_keys:
        raise ValueError("PDB serialization changed whole-chain atom identities")
    delta = max(float(np.linalg.norm(np.asarray(actual[k]["coordinates"]) - atoms[k])) for k in actual)
    if delta > 0.0009:
        raise ValueError("PDB serialization changed protein coordinates beyond source rounding")
    return text, {"rule": "all complete source protein chains with any heavy atom within6 A of any point of the search box",
                  "selected_source_chains": selected, "source_atom_count": len(keys), "selected_heavy_atoms": len(selected_keys),
                  "source_chain_count": len({key[0] for key in keys}), "closest_omitted_atom_to_box_A": nearest,
                  "interaction_cutoff_A": cutoff, "max_PDB_serialization_displacement_A": delta,
                  "source_heavy_atoms_unchanged": True, "native_coordinate_access": False}


def _prepare_numbered_protein(text, stem, definition, root):
    """Expose only independently curated sequence gaps; retain source output IDs."""
    mapping_name = definition.get("receptor_numbering_map") if definition else None
    if not mapping_name:
        return prepare_receptor(text, stem)
    path = Path(root) / mapping_name
    if sha256(path) != definition["receptor_numbering_map_sha256"]:
        raise ValueError("Curated protein sequence numbering map changed")
    mapping = json.loads(path.read_text())
    forward, reverse = {}, {}
    for row in mapping["residues"]:
        if row["prepared_resname"] not in PROTEIN_RESIDUES:
            continue  # The site arm already removed every metal uniformly.
        source = (row["prepared_chain"], row["prepared_residue_number"], "", row["prepared_resname"])
        author = (row["deposited_author_chain"], row["deposited_author_number"], row["deposited_insertion_code"], row["deposited_resname"])
        if source in forward or author in reverse or source[3] != author[3]:
            raise ValueError("Protein-only numbering map is not a chemical bijection")
        forward[source], reverse[author] = author, source
    observed = {key[:4] for key in _pdb_atoms(text)}
    if not observed.issubset(forward):
        raise ValueError("Selected protein contains residues absent from source numbering map")
    def replace(content, replacements):
        lines = []
        for line in content.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())
                if key not in replacements:
                    raise ValueError("Protein atom has no sequence mapping")
                ch, num, ins, res = replacements[key]
                line = line[:17] + f"{res:3s}" + line[20:21] + f"{ch:1s}{num:4d}{ins or ' ':1s}" + line[27:]
            lines.append(line)
        return "\n".join(lines) + "\n"
    log = prepare_receptor(replace(text, forward), stem)
    stem.with_suffix(".pdbqt").write_text(replace(stem.with_suffix(".pdbqt").read_text(), reverse))
    stem.with_suffix(".pdb").write_text(text)
    log.update(source_numbering_map=mapping_name, source_numbering_map_sha256=sha256(path),
               identity_frame="preparation audit uses deposited sequence IDs; output restored to released source IDs",
               removed_metals_not_reintroduced=True)
    write_json(stem.with_name(stem.name + "_atom_changes.json"), log)
    return log


def _cached_protein(text, stem, definition, root, cache_dir):
    """Reuse an exact audited whole-chain input, never a near-coordinate match."""
    if cache_dir is None:
        return _prepare_numbered_protein(text, stem, definition, root), None
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key_data = {"protein_PDB_sha256": hashlib.sha256(text.encode()).hexdigest(),
                "blind_preparation_code_sha256": sha256(Path(__file__)),
                "protein_preparation_code_sha256": sha256(Path(__file__).with_name("aav.py")),
                "meeko_version": importlib.metadata.version("meeko"),
                "numbering_map_sha256": definition.get("receptor_numbering_map_sha256") if definition else None}
    key = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
    entry = cache_dir / key
    names = ("neighborhood.pdb", "neighborhood.pdbqt", "neighborhood_atom_changes.json")
    if not entry.exists():
        temporary = Path(tempfile.mkdtemp(prefix=".protein-cache-", dir=cache_dir))
        try:
            _prepare_numbered_protein(text, temporary / "neighborhood", definition, root)
            write_json(temporary / "manifest.json", {"key": key, "inputs": key_data,
                       "files": {name: sha256(temporary / name) for name in names}})
            try:
                temporary.rename(entry)
            except OSError:
                if not entry.is_dir():
                    raise
                # Another worker produced the same exact input first. Its
                # complete manifest and all bytes are checked below.
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
    manifest_path = entry / "manifest.json"
    if entry.is_symlink() or manifest_path.is_symlink():
        raise ValueError("Audited receptor cache must not use symlinks")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("key") != key or manifest.get("inputs") != key_data or set(manifest.get("files", {})) != set(names):
        raise ValueError("Audited receptor cache identity differs")
    for name, digest in manifest["files"].items():
        source = entry / name
        if source.is_symlink() or sha256(source) != digest:
            raise ValueError("Audited receptor cache bytes changed")
        destination = stem.with_suffix(".pdb" if name.endswith(".pdb") else ".pdbqt") if name != names[2] else stem.with_name(stem.name + "_atom_changes.json")
        shutil.copyfile(source, destination)
    return json.loads(stem.with_name(stem.name + "_atom_changes.json").read_text()), {
        "key": key, "input_fingerprint": key_data, "manifest_sha256": sha256(manifest_path),
        "files": manifest["files"], "policy": "exact audited full-chain PDB input plus preparation code/version; no geometry approximation"}


def prepare_box(receptor_cif, ligand, output_dir, *, center, condition_id, complex_id,
                method, candidate_class_rank, state="holo", configuration_path,
                candidate_original_rank=None, definition=None, root=None, receptor_cache_dir=None):
    """Prepare one independently centered box using a uniform dry protein model."""
    receptor_cif, config_path = Path(receptor_cif), Path(configuration_path)
    center = _finite_vector(center, "Search center")
    size = _finite_vector(ligand["box_size"], "Search size", positive=True)
    if (method not in {*_PREDICTORS, "oracle"} or isinstance(candidate_class_rank, bool)
            or not isinstance(candidate_class_rank, int) or candidate_class_rank < 1 or candidate_class_rank > 5):
        raise ValueError("Unknown method or candidate class outside the fixed budget")
    output = _fresh(output_dir)
    attempt = {"condition_id": condition_id, "complex_id": complex_id, "method": method,
               "candidate_class_rank": candidate_class_rank, "candidate_original_rank": candidate_original_rank,
               "status": "preparing", "configuration_sha256": sha256(config_path)}
    write_json(output / "preparation_attempt.json", attempt)
    try:
        protein, neighborhood = _protein_text(receptor_cif, center, size, cutoff=6.0)
        log, cache = _cached_protein(protein, output / "neighborhood", definition, root, receptor_cache_dir)
        moved, translation = translate_ligand(ligand["pdbqt"], center)
        write_text_atomic(output / "ligand.pdbqt", moved)
        # This is an atom-map-only runner contract, never a native coordinate reference.
        write_json(output / "reference.json", {"names": ligand["names"], "pose_reference_eligible": False,
                                                "purpose": "input heavy-atom identities only; atomic truth is evaluated separately"})
        prep = {**attempt, "case_id": complex_id, "status": "prepared", "cohort_role": "development_only",
                "family_id": complex_id, "information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                "execution_information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                "box_center": center.tolist(), "box_size": size.tolist(), "ligand_file": "ligand.pdbqt",
                "receptors": [{"id": "neighborhood", "state": state, "assembly_condition": "complete_box_neighborhood",
                               "pdb_file": "neighborhood.pdb", "pdbqt_file": "neighborhood.pdbqt"}],
                "source_receptor_sha256": sha256(receptor_cif), "neighborhood": neighborhood,
                "ligand_source": {k: v for k, v in ligand.items() if k != "pdbqt"},
                "ligand_translation": translation, "ligand_sha256": sha256(output / "ligand.pdbqt"),
                "uniform_model": "protein-only, dry, rigid; all metals and glycans removed before candidate selection",
                "native_ligand_coordinate_access": method == "oracle", "native_coordinates_in_runner_reference": False,
                "preparation_code_sha256": sha256(Path(__file__)), "receptor_preparation": log}
        if cache is not None:
            prep["receptor_cache"] = cache
        write_json(output / "preparation.json", prep)
        write_json(output / "preparation_attempt.json", {**attempt, "status": "ok"})
        return prep
    except (ValueError, RuntimeError, OSError) as error:
        write_json(output / "preparation_attempt.json", {**attempt, "status": "failed", "error_type": type(error).__name__, "reason": str(error)})
        raise


def _condition(root, condition_id):
    config = json.loads((Path(root) / "benchmarks/site_experiment.json").read_text())
    rows = [row for row in config["cases"] if row["condition_id"] == condition_id]
    if len(rows) != 1:
        raise ValueError("Unknown frozen site condition")
    return rows[0]


def prepare_condition(root, condition_id, method, output_dir, configuration_path, *, oracle_center=None, receptor_cache_dir=None):
    """Prepare the fixed available candidate boxes, retaining every failed condition.

    An oracle center is accepted only for the explicitly labeled oracle arm.
    No native source, evaluation record, or ligand atom-coordinate reference is
    consulted to generate predicted boxes or their neighborhood boundaries.
    """
    root, config_path = Path(root), Path(configuration_path)
    case = _condition(root, condition_id)
    output = _fresh(output_dir)
    manifest = {"schema_version": 1, "condition_id": condition_id, "complex_id": case["complex_id"],
                "state": case["state"], "method": method, "configuration_sha256": sha256(config_path),
                "eligible_condition": True, "max_candidate_boxes": 1 if method == "oracle" else 5,
                "jobs": [], "status": "preparing", "native_coordinate_access": method == "oracle"}
    try:
        if method == "oracle":
            center = _finite_vector(oracle_center, "Explicit oracle control center")
            candidates = [{"rank": 1, "original_rank": 1, "center_A": center.tolist()}]
            manifest["source_center_policy"] = "separate caller supplies native-informed center for labeled oracle control"
        else:
            if oracle_center is not None:
                raise ValueError("Native-informed centers are forbidden in the predicted arm")
            selection = frozen_candidates(root, condition_id, method)
            manifest["candidate_selection"] = selection
            candidates = selection["candidates"]
            if selection["status"] != "ok":
                manifest["status"] = selection["status"]
                write_json(output / "condition.json", manifest)
                return manifest
        ligand = ligand_source(root, case["complex_id"])
        manifest["ligand_source"] = {k: v for k, v in ligand.items() if k != "pdbqt"}
        restore_site_inputs(root, condition_id, "original", output / "source_receptor")
        receptor = output / "source_receptor/receptor.cif"
        if method != "oracle" and sha256(receptor) != selection["prediction_input_sha256"]:
            raise ValueError("The docking receptor is not the frozen predictor receptor")
        definition = None
        if case["complex_id"] != "aav9":
            external = json.loads((root / "benchmarks/external_cohort.json").read_text())
            definition = next(row for row in external["cases"] if row["case_id"] == case["complex_id"])
        for candidate in candidates:
            name = f"box{candidate['rank']:02d}"
            job = {"id": name, "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"],
                   "center_A": candidate["center_A"], "prepared_directory": name}
            try:
                prepare_box(receptor, ligand, output / name, center=candidate["center_A"],
                            condition_id=condition_id, complex_id=case["complex_id"], method=method,
                            candidate_class_rank=candidate["rank"], candidate_original_rank=candidate["original_rank"],
                            state=case["state"], configuration_path=config_path, definition=definition, root=root,
                            receptor_cache_dir=receptor_cache_dir)
                job["status"] = "prepared"
                job["preparation_sha256"] = sha256(output / name / "preparation.json")
            except (ValueError, RuntimeError, OSError) as error:
                job.update(status="preparation_failed", reason=str(error), error_type=type(error).__name__)
            manifest["jobs"].append(job)
            write_json(output / "condition.json", manifest)
        manifest["status"] = "prepared" if all(j["status"] == "prepared" for j in manifest["jobs"]) else "partial_preparation_failure"
    except (ValueError, RuntimeError, OSError) as error:
        manifest.update(status="preparation_failed", reason=str(error), error_type=type(error).__name__)
    write_json(output / "condition.json", manifest)
    return manifest


def prepare_predicted_condition(root, condition_id, method, output_dir, configuration_path, *, receptor_cache_dir=None):
    return prepare_condition(root, condition_id, method, output_dir, configuration_path, receptor_cache_dir=receptor_cache_dir)


def planning_inventory(root):
    """Source-only inventory for a plan frozen before the new docking outputs."""
    root = Path(root)
    cases = json.loads((root / "benchmarks/site_experiment.json").read_text())["cases"]
    rows = []
    for case in cases:
        try:
            ligand = ligand_source(root, case["complex_id"])
            ligand_info = {key: value for key, value in ligand.items() if key not in {"pdbqt", "chemistry"}}
            ligand_status = "ok"
        except ValueError as error:
            ligand_info, ligand_status = {"reason": str(error)}, "preparation_failed"
        for method in (*_PREDICTORS, "oracle"):
            selection = None if method == "oracle" else frozen_candidates(root, case["condition_id"], method)
            rows.append({"condition_id": case["condition_id"], "complex_id": case["complex_id"], "state": case["state"],
                         "method": method, "eligible_condition": True, "ligand_status": ligand_status,
                         "ligand_source": ligand_info, "candidate_selection": selection,
                         "available_box_count": 1 if method == "oracle" else len(selection["candidates"]),
                         "native_coordinate_access": False})
    return {"schema_version": 1, "purpose": "source-only design inventory; no new docking or atomic accuracy evaluation",
            "eligible_complexes": len({c["complex_id"] for c in cases}), "eligible_conditions": len(cases),
            "eligible_method_conditions": len(rows), "rows": rows,
            "candidate_and_ligand_code_sha256": sha256(Path(__file__))}


def prepare_oracle_condition(root, condition_id, center, output_dir, configuration_path, *, receptor_cache_dir=None):
    return prepare_condition(root, condition_id, "oracle", output_dir, configuration_path, oracle_center=center,
                             receptor_cache_dir=receptor_cache_dir)


def run_condition(prepared_dir, output_dir, executable, configuration_path, *, resume=False):
    """Execute prepared boxes with the existing strict engine; preserve failures."""
    from .engine import dock
    prepared, config_path = Path(prepared_dir), Path(configuration_path)
    condition = json.loads((prepared / "condition.json").read_text())
    if condition["configuration_sha256"] != sha256(config_path):
        raise ValueError("Frozen configuration differs from preparation")
    output = Path(output_dir)
    if not resume:
        _fresh(output)
    else:
        output.mkdir(parents=True, exist_ok=True)
    runs = []
    for job in condition["jobs"]:
        if job["status"] != "prepared":
            runs.append({**job, "execution_status": "preparation_failed", "runs": []})
            continue
        directory = prepared / job["prepared_directory"]
        if sha256(directory / "preparation.json") != job["preparation_sha256"]:
            raise ValueError("Prepared candidate changed before execution")
        result = dock(directory, output / job["id"], executable, config_path, resume=resume)
        runs.append({**job, "execution_status": "ok" if all(row["status"] == "ok" for row in result) else "execution_failed", "runs": result})
    receipt = {"schema_version": 1, "condition_id": condition["condition_id"], "complex_id": condition["complex_id"],
               "method": condition["method"], "preparation_status": condition["status"],
               "condition_manifest_sha256": sha256(prepared / "condition.json"),
               "configuration_sha256": sha256(config_path), "eligible_condition": True,
               "jobs": runs, "status": "ok" if runs and all(j["execution_status"] == "ok" for j in runs) else "workflow_failed",
               "note": "Missing/failed predictor or chemistry remains in eligible denominator; no inferred biological nonbinding"}
    write_json(output / "condition_runs.json", receipt)
    return receipt
