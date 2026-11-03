"""Portable evaluation of the frozen predicted-box/oracle pose comparison.

Atomic truth enters only after the fixed workflow ledger has been checked. The
default refuses incomplete execution; a pending-only inspection never reports
atomic accuracy. Qualified full models and modeled cores retain separate
denominators. No ligand-only superposition is performed.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
from statistics import mean
import tempfile

import numpy as np
from scipy.spatial import cKDTree

from .aav import PROTEIN_RESIDUES
from .engine import read_poses
from .external import pdb_heavy_atoms
from .metrics import contact_metrics, contact_pairs, fit_receptor
from .provenance import sha256, write_json, write_text_atomic
from .site_geometry import read_receptor
from .site_prediction_verification import restore_site_inputs
from .site_reproduction import compare_site_values


_DERIVED_FLOAT_FIELDS = frozenset({
    "global_RMSD_A", "local_interface_RMSD_A", "top1_RMSD_A", "top5_min_RMSD_A",
    "best_returned_RMSD_A", "backbone_RMSD_A", "rotation", "translation",
    "known_translation_RMSD_A",
    "full_parent_min_protein_distance_A", "contact_recall", "contact_precision",
    "interface_residue_recall", "interface_residue_precision", "interface_residue_jaccard",
    "max_heavy_bond_length_change_from_input_A", "max_ring_pair_distance_change_from_input_A",
    "distance_A", "input_distance_A",
})


class Evidence:
    """Checksums of the actual portable files and decompressed input bytes."""
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.sources = {}
        self._bundles = {}

    def read(self, relative, expected=None):
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Evidence paths must remain repository-relative")
        full = self.root / path
        if full.is_symlink():
            raise ValueError("Evidence must be an ordinary file")
        raw = full.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError(f"Evidence checksum changed: {path}")
        self.sources[path.as_posix()] = {"sha256": digest, "bytes": len(raw)}
        return raw

    def json(self, relative, expected=None):
        return json.loads(self.read(relative, expected))

    def bundle(self, directory):
        directory = Path(directory)
        key = directory.as_posix()
        if key not in self._bundles:
            manifest = self.json(directory / "inputs/manifest.json")
            files = {}
            for entry in manifest["files"]:
                name, archive = entry["filename"], entry["archive"]
                if (any(Path(x).name != x or x in {".", ".."} for x in (name, archive))
                        or name in files):
                    raise ValueError("Invalid or duplicate archived input filename")
                raw = self.read(directory / "inputs" / archive, entry["archive_sha256"])
                content = gzip.decompress(raw)
                if (hashlib.sha256(content).hexdigest() != entry["sha256"]
                        or ("original_bytes" in entry and len(content) != entry["original_bytes"])):
                    raise ValueError("Decompressed prepared input changed")
                files[name] = content
            self._bundles[key] = files
        return self._bundles[key]


def _key(value):
    chain, number, insertion, residue, atom = value
    return str(chain), int(number), str(insertion).strip(), str(residue), str(atom)


def qualified_atom_map(record, reference):
    """Map the exact frozen atomic scope to the full parent ligand names."""
    names = reference["names"]
    coordinates = np.asarray(reference["coordinates"], dtype=float)
    atoms = reference["atom_map"]
    if (len(names) != len(set(names)) or coordinates.shape != (len(names), 3)
            or not np.isfinite(coordinates).all() or len(atoms) != len(names)):
        raise ValueError("Parent ligand reference has missing or duplicate atom identities")
    by_source = {}
    for atom in atoms:
        source = _key(atom["source_key"])
        if source in by_source or atom["name"] not in names:
            raise ValueError("Parent source atom map is not a complete bijection")
        by_source[source] = names.index(atom["name"])
    selected, native, keys = [], [], []
    for atom in record["scored_atom_map"]:
        source = (atom["prepared_chain"], int(atom["prepared_resnum"]), "", atom["prepared_resname"], atom["atom_name"])
        if source not in by_source or source in keys:
            raise ValueError("A qualified atomic-reference atom is absent or duplicated in the full input graph")
        index = by_source[source]
        xyz = np.asarray(atom["xyz_A"], dtype=float)
        if xyz.shape != (3,) or not np.isfinite(xyz).all() or np.linalg.norm(coordinates[index] - xyz) > 0.001:
            raise ValueError("Qualified and archived native atom coordinates disagree")
        selected.append(index)
        native.append(xyz.tolist())
        keys.append(source)
    if (not selected or len(selected) != record["scored_atom_count"]
            or len(set(selected)) != len(selected)):
        raise ValueError("Qualified scored atom denominator is incomplete")
    return {"indices": selected, "names": [names[i] for i in selected], "coordinates": native,
            "source_atom_keys": [list(key) for key in keys], "parent_atom_count": len(names)}


def pooled_poses(boxes, seed, site_budget):
    """Pool fixed returned candidates using score and prespecified identity ties."""
    rows = []
    for box in boxes:
        if box["candidate_class_rank"] > site_budget:
            continue
        for pose in box["poses_by_seed"].get(seed, []):
            rows.append({**pose, "candidate_class_rank": box["candidate_class_rank"],
                         "candidate_original_rank": box["candidate_original_rank"], "box_id": box["id"]})
    rows.sort(key=lambda x: (x["vina_score"], x["candidate_class_rank"], x["rank"]))
    return [{**row, "pooled_rank": i} for i, row in enumerate(rows, 1)]


def accuracy_summary(rows, field):
    """Zero recovery for failed eligible pools; no denominator shrinking."""
    values = [float(row[field]) for row in rows]
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Invalid pose-distance metric")
    result = {"returned_pose_count": len(rows), "top1_RMSD_A": values[0] if values else None,
              "top5_min_RMSD_A": min(values[:5]) if values else None,
              "best_returned_RMSD_A": min(values) if values else None,
              "best_returned_pooled_rank": int(np.argmin(values)) + 1 if values else None}
    for threshold in (1, 2, 3):
        result[f"top1_success_{threshold}A"] = bool(values and values[0] <= threshold)
        result[f"top5_success_{threshold}A"] = any(v <= threshold for v in values[:5])
        result[f"best_returned_success_{threshold}A"] = any(v <= threshold for v in values)
    return result


def native_metric_controls(native, protein):
    native = np.asarray(native, dtype=float)
    protein = np.asarray(protein, dtype=float)
    if native.ndim != 2 or native.shape[1:] != (3,) or not len(native):
        raise ValueError("Native control requires the complete qualified atomic map")
    contacts = contact_pairs(native, protein)
    exact = float(np.sqrt(np.mean(np.sum((native - native) ** 2, axis=1))))
    shifted = native + [3, 4, 0]
    shift = float(np.sqrt(np.mean(np.sum((shifted - native) ** 2, axis=1))))
    metrics = contact_metrics(native, protein, contacts)
    if exact != 0 or abs(shift - 5) > 1e-12 or metrics["contact_recall"] != 1 or metrics["contact_precision"] != 1:
        raise ValueError("Qualified native-coordinate controls failed")
    return {"native_RMSD_A": exact, "native_contact_recall": metrics["contact_recall"],
            "native_contact_precision": metrics["contact_precision"], "native_contact_count": len(contacts),
            "known_translation_A": [3, 4, 0], "known_translation_RMSD_A": shift,
            "alignment": "identity-only analytic controls; never ligand-fitted; docking endpoint uses receptor fits"}


def _validate_command(row, prep, config, condition, method, box):
    command = row.get("command", [])
    if not command or len(command) % 2 != 1 or any(not isinstance(v, str) for v in command):
        raise ValueError("Missing or malformed actual Vina command")
    options = dict(zip(command[1::2], command[2::2]))
    if len(options) * 2 != len(command) - 1:
        raise ValueError("Duplicate actual Vina command flag")
    expected = {"--scoring": "vina", "--seed": row["seed"], "--cpu": config["cpu"],
                "--exhaustiveness": config["exhaustiveness"], "--num_modes": config["num_modes"],
                "--min_rmsd": config["min_rmsd_A"], "--energy_range": config["energy_range_kcal_mol"]}
    for axis, center, size in zip("xyz", prep["box_center"], prep["box_size"]):
        expected["--center_" + axis], expected["--size_" + axis] = center, size
    for flag, value in expected.items():
        actual = options.get(flag)
        if (actual is None or (isinstance(value, str) and actual != value)
                or (not isinstance(value, str) and (not math.isfinite(float(actual)) or float(actual) != float(value)))):
            raise ValueError(f"Actual Vina command changed frozen option {flag}")
    suffixes = {"--receptor": Path("data/processed/blind_pose") / condition / method / box / "neighborhood.pdbqt",
                "--ligand": Path("data/processed/blind_pose") / condition / method / box / "ligand.pdbqt",
                "--out": Path("results/blind_pose") / condition / method / box / row["pose_file"]}
    for flag, suffix in suffixes.items():
        if flag not in options or Path(options[flag]).parts[-len(suffix.parts):] != suffix.parts:
            raise ValueError("Actual Vina command refers to a different condition/box artifact")
    if set(options) != set(expected) | set(suffixes):
        raise ValueError("Unexpected extra Vina invocation option")


def _pose_log(content, poses):
    matches = []
    for line in content.decode().splitlines():
        match = re.match(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s*$", line)
        if match:
            matches.append((int(match.group(1)), float(match.group(2))))
    if len(matches) != len(poses) or any(rank != pose["rank"] or abs(score - pose["vina_score"]) > 0.002
                                          for (rank, score), pose in zip(matches, poses)):
        raise ValueError("Raw Vina log score table disagrees with retained poses")


def _all_execution(evidence, config, freeze, inventory, oracle):
    """Validate every workflow first; return no atomic endpoint for a partial run."""
    frozen_rows = {(r["condition_id"], r["method"]): r for r in inventory["rows"]}
    if len(frozen_rows) != 27 or len(config["cases"]) != 9 or config["methods"] != ["fpocket", "p2rank", "oracle"]:
        raise ValueError("Frozen workflow denominator differs from the declared 27 conditions")
    vina_lock = evidence.json("benchmarks/vina.lock.json")
    executable_hashes = {asset["sha256"] for asset in vina_lock["assets"].values()}
    conditions, pending = [], []
    total_boxes = total_terminal = total_success = 0
    for case in config["cases"]:
        cid = case["condition_id"]
        for method in config["methods"]:
            directory = Path("results/blind_pose") / cid / method
            condition_path = directory / "condition_preparation.json"
            condition_raw = evidence.read(condition_path)
            condition = json.loads(condition_raw)
            if (condition.get("condition_id") != cid or condition.get("complex_id") != case["complex_id"]
                    or condition.get("method") != method or condition.get("configuration_sha256") != freeze["configuration_sha256"]
                    or condition.get("eligible_condition") is not True):
                raise ValueError("Condition metadata differs from the frozen scientific denominator")
            expected = frozen_rows[(cid, method)]
            expected_boxes = expected["available_box_count"] if expected["ligand_status"] == "ok" else 0
            if (len(condition["jobs"]) != expected_boxes
                    or [job["candidate_class_rank"] for job in condition["jobs"]] != list(range(1, expected_boxes + 1))):
                raise ValueError("Per-condition candidate-box denominator changed")
            if method != "oracle" and condition.get("candidate_selection") != expected["candidate_selection"]:
                raise ValueError("Frozen predictor candidate identities/order changed")
            if method == "oracle" and condition["jobs"] and condition["jobs"][0]["center_A"] != oracle["controls"][cid]["center_A"]:
                raise ValueError("Oracle center differs from its separate frozen control")
            result = {"condition_id": cid, "complex_id": case["complex_id"], "state": case["state"], "method": method,
                      "preparation_status": condition["status"], "boxes": []}
            for job in condition["jobs"]:
                box_id = job["id"]
                if box_id != f"box{job['candidate_class_rank']:02d}" or job["prepared_directory"] != box_id:
                    raise ValueError("Invalid candidate box identity")
                box = {"id": box_id, "candidate_class_rank": job["candidate_class_rank"],
                       "candidate_original_rank": job["candidate_original_rank"], "poses_by_seed": {}, "runs": []}
                if job["status"] != "prepared":
                    box["status"] = "preparation_failed"
                    result["boxes"].append(box)
                    continue
                total_boxes += 1
                relative = directory / box_id
                files = evidence.bundle(relative)
                preparation_raw = files["preparation.json"]
                if hashlib.sha256(preparation_raw).hexdigest() != job["preparation_sha256"]:
                    raise ValueError("Prepared box bytes differ from condition manifest")
                prep = json.loads(preparation_raw)
                reference = json.loads(files["reference.json"])
                regime = "true-site-local" if method == "oracle" else "predicted-site-local"
                if (prep["condition_id"] != cid or prep["method"] != method or prep["complex_id"] != case["complex_id"]
                        or prep["box_center"] != job["center_A"] or prep["configuration_sha256"] != freeze["configuration_sha256"]
                        or prep["box_size"] != expected["ligand_source"]["box_size"]
                        or prep["execution_information_regime"] != regime or reference.get("pose_reference_eligible") is not False
                        or set(reference) != {"names", "pose_reference_eligible", "purpose"}):
                    raise ValueError("Prepared box leaked or changed its scientific input contract")
                if prep["neighborhood"]["interaction_cutoff_A"] != 6:
                    raise ValueError("The original frozen six-Angstrom neighborhood rule changed")
                box["preparation"] = prep
                box["names"] = reference["names"]
                runs_path = relative / "runs.json"
                if not (evidence.root / runs_path).is_file():
                    pending.append(f"{cid}/{method}/{box_id}:missing")
                    result["boxes"].append(box)
                    continue
                evidence.read(relative / "configuration.json", freeze["configuration_sha256"])
                rows = evidence.json(runs_path)
                identities = [row["seed"] for row in rows]
                if len(set(identities)) != len(rows) or set(identities) - set(config["seeds"]):
                    raise ValueError("Duplicate or unknown receptor-by-seed identity")
                if set(identities) != set(config["seeds"]):
                    pending.append(f"{cid}/{method}/{box_id}:missing_seeds")
                for row in rows:
                    box["runs"].append(row)
                    if row["status"] in {"running", "pending"}:
                        pending.append(f"{cid}/{method}/{box_id}/{row['seed']}:running")
                        continue
                    if row["status"] not in {"ok", "failed"}:
                        raise ValueError("Unknown terminal docking state")
                    total_terminal += 1
                    fingerprint = row["fingerprint"]
                    expected_fp = {"config_sha256": freeze["configuration_sha256"], "preparation_sha256": job["preparation_sha256"],
                                   "reference_sha256": hashlib.sha256(files["reference.json"]).hexdigest(),
                                   "ligand_sha256": hashlib.sha256(files["ligand.pdbqt"]).hexdigest(),
                                   "receptor_sha256": hashlib.sha256(files["neighborhood.pdbqt"]).hexdigest(),
                                   "receptor_pdb_sha256": hashlib.sha256(files["neighborhood.pdb"]).hexdigest(), "full_assembly_sha256": None}
                    if (any(fingerprint.get(k) != value for k, value in expected_fp.items())
                            or fingerprint.get("executable_sha256") not in executable_hashes
                            or row.get("case_id") != case["complex_id"] or row.get("receptor") != "neighborhood"
                            or row.get("information_regime") != regime or row.get("vina_version") != "AutoDock Vina v1.2.7"):
                        raise ValueError("Run fingerprints or input-information label changed")
                    _validate_command(row, prep, config, cid, method, box_id)
                    for name in (row["pose_file"], row["log_file"]):
                        if Path(name).name != name or name in {".", ".."}:
                            raise ValueError("Invalid raw result artifact filename")
                    log = evidence.read(relative / row["log_file"])
                    if row["status"] == "ok":
                        if row.get("exit_code") != 0:
                            raise ValueError("Successful run has nonzero exit status")
                        raw_pose = evidence.read(relative / row["pose_file"], row["pose_sha256"])
                        with tempfile.TemporaryDirectory() as temporary:
                            path = Path(temporary) / "poses.pdbqt"
                            path.write_bytes(raw_pose)
                            poses = read_poses(path, reference["names"])
                        if len(poses) != row["returned_poses"] or not 1 <= len(poses) <= config["num_modes"]:
                            raise ValueError("Returned pose denominator changed")
                        _pose_log(log, poses)
                        box["poses_by_seed"][row["seed"]] = poses
                        total_success += 1
                result["boxes"].append(box)
            if len({box["candidate_class_rank"] for box in result["boxes"]}) != len(result["boxes"]):
                raise ValueError("Duplicate candidate class in a condition")
            conditions.append(result)
    if len(conditions) != 27 or total_boxes != 52:
        raise ValueError("The frozen 27-condition / 52-box workflow denominator changed")
    return conditions, sorted(pending), {"eligible_conditions": 27, "prepared_boxes": total_boxes,
        "scheduled_seed_runs": 156, "terminal_seed_runs": total_terminal, "successful_seed_runs": total_success,
        "failed_seed_runs": total_terminal - total_success}


def _old_inputs(evidence, case_id):
    directory = "154l" if case_id == "154L" else case_id
    relative = Path("results/development") / directory
    if not (evidence.root / relative / "inputs/manifest.json").exists():
        return None
    return evidence.bundle(relative)


def _check_atomic_sources(evidence, ledger):
    evidence.read(ledger["policy_archive"], ledger["policy_sha256"])
    evidence.read(ledger["manual_adjudication_archive"], ledger["manual_adjudication_sha256"])
    for case in ledger["cases"]:
        for source in case["sources"].values():
            evidence.read(source["archive"], source["sha256"])
        for key in ("released_ligand_source", "released_receptor_source"):
            source = case[key]
            compressed = evidence.read(source["archive"], source["compressed_sha256"])
            if hashlib.sha256(gzip.decompress(compressed)).hexdigest() != source["uncompressed_sha256"]:
                raise ValueError("Native source decompression checksum changed")


def _native_protein(evidence, record):
    source = record["released_receptor_source"]
    text = gzip.decompress(evidence.read(source["archive"], source["compressed_sha256"])).decode()
    atoms = [atom for atom in pdb_heavy_atoms(text) if atom["key"][3] in PROTEIN_RESIDUES]
    keys = [_key(atom["key"]) for atom in atoms]
    xyz = np.asarray([atom["coordinates"] for atom in atoms], dtype=float)
    if not keys or len(set(keys)) != len(keys) or not np.isfinite(xyz).all():
        raise ValueError("Native protein reference is empty, duplicated or nonfinite")
    return keys, xyz


def _protein_alignment(keys, xyz, native_keys, native_xyz, native_core):
    source, target = dict(zip(keys, xyz)), dict(zip(native_keys, native_xyz))
    if set(source) != set(target):
        raise ValueError("Frozen full input protein and native protein atom identities differ")
    native_contacts = contact_pairs(native_core, native_xyz)
    residues = {native_keys[j][:4] for _, j in native_contacts}
    fits = {}
    for label in ("global", "local_interface"):
        selected = [key for key in native_keys if key[-1] in {"N", "CA", "C", "O"}
                    and (label == "global" or key[:4] in residues)]
        a, b = np.asarray([source[k] for k in selected]), np.asarray([target[k] for k in selected])
        if len(selected) < 3 or np.linalg.matrix_rank(a - a.mean(axis=0), tol=1e-8) < 2:
            raise ValueError("Receptor-only backbone fit is underdetermined")
        rotation, translation = fit_receptor(a, b)
        fits[label] = {"rotation": rotation, "translation": translation,
                       "metadata": {"matched_backbone_atom_count": len(selected),
                                    "backbone_RMSD_A": float(np.sqrt(np.mean(np.sum((a @ rotation + translation - b) ** 2, axis=1)))),
                                    "rotation": rotation.tolist(), "translation": translation.tolist(),
                                    "selection": "all corresponding protein N/CA/C/O" if label == "global" else "N/CA/C/O of native qualified-core contact residues; evaluation only"}}
    return fits, native_contacts, residues


class ParentGeometry:
    """Fixed parent graph diagnostics, independent from the RMSD success event."""
    def __init__(self, sdf, reference, mapping):
        from rdkit import Chem
        supplier = Chem.ForwardSDMolSupplier(io.BytesIO(sdf), removeHs=False)
        molecules = list(supplier)
        if len(molecules) != 1 or molecules[0] is None:
            raise ValueError("Archived input SDF does not define one chemical graph")
        original = molecules[0]
        original_indices = [atom.GetIdx() for atom in original.GetAtoms() if atom.GetAtomicNum() != 1]
        heavy = Chem.RemoveHs(original)
        if heavy.GetNumAtoms() != len(reference["names"]) or set(original_indices) != {a["rdkit_index"] for a in reference["atom_map"]}:
            raise ValueError("Archived graph and input heavy-atom map differ")
        original_to_heavy = {index: i for i, index in enumerate(original_indices)}
        atom_by_name = {atom["name"]: atom for atom in reference["atom_map"]}
        self.name_to_heavy = [original_to_heavy[atom_by_name[name]["rdkit_index"]] for name in reference["names"]]
        self.heavy_to_name = np.argsort(self.name_to_heavy)
        self.molecule = heavy
        self.expected_smiles = Chem.MolToSmiles(heavy, isomericSmiles=True)
        self.input_xyz = np.asarray(heavy.GetConformer().GetPositions())
        self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in heavy.GetBonds()]
        self.rings = [list(ring) for ring in heavy.GetRingInfo().AtomRings()]
        by_heavy = {original_to_heavy[atom["rdkit_index"]]: _key(atom["source_key"]) for atom in reference["atom_map"]}
        self.links = [(a, b) for a, b in self.bonds if by_heavy[a][:4] != by_heavy[b][:4]]
        permutations = heavy.GetSubstructMatches(heavy, uniquify=False, useChirality=True, maxMatches=4096)
        if len(permutations) >= 4096:
            raise ValueError("Chemical automorphism enumeration exceeded its explicit complete-enumeration limit")
        scored_heavy = [self.name_to_heavy[i] for i in mapping["indices"]]
        unique = set()
        self.scored_maps = []
        for permutation in permutations:
            full_indices = tuple(int(self.heavy_to_name[permutation[i]]) for i in scored_heavy)
            if full_indices not in unique:
                unique.add(full_indices)
                self.scored_maps.append(list(full_indices))
        if not self.scored_maps:
            raise ValueError("No graph/stereochemistry-valid scored atom map")

    def diagnostics(self, full_pose):
        from rdkit import Chem
        full_pose = np.asarray(full_pose, dtype=float)
        if full_pose.shape != (len(self.name_to_heavy), 3) or not np.isfinite(full_pose).all():
            raise ValueError("Full-parent pose atom loss or nonfinite coordinate")
        xyz = full_pose[self.heavy_to_name]
        candidate = Chem.Mol(self.molecule)
        candidate.RemoveAllConformers()
        conformer = Chem.Conformer(len(xyz))
        for i, point in enumerate(xyz):
            conformer.SetAtomPosition(i, point.tolist())
        candidate.AddConformer(conformer)
        Chem.RemoveStereochemistry(candidate)
        Chem.AssignAtomChiralTagsFromStructure(candidate, confId=0, replaceExistingTags=True)
        Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
        ring_delta = 0.0
        for ring in self.rings:
            actual = np.linalg.norm(xyz[ring][:, None, :] - xyz[ring][None, :, :], axis=2)
            source = np.linalg.norm(self.input_xyz[ring][:, None, :] - self.input_xyz[ring][None, :, :], axis=2)
            ring_delta = max(ring_delta, float(np.abs(actual - source).max()))
        bond_delta = max((abs(float(np.linalg.norm(xyz[a] - xyz[b])) - float(np.linalg.norm(self.input_xyz[a] - self.input_xyz[b])))
                          for a, b in self.bonds), default=0.0)
        return {"fixed_parent_graph_and_heavy_map_preserved": True,
                "stereochemistry_from_pose_3D_matches_input": Chem.MolToSmiles(candidate, isomericSmiles=True) == self.expected_smiles,
                "max_heavy_bond_length_change_from_input_A": bond_delta,
                "max_ring_pair_distance_change_from_input_A": ring_delta,
                "glycosidic_links": [{"heavy_indices": [a, b], "distance_A": float(np.linalg.norm(xyz[a] - xyz[b])),
                                       "input_distance_A": float(np.linalg.norm(self.input_xyz[a] - self.input_xyz[b]))} for a, b in self.links],
                "ring_pucker_validity_certified": False, "glycosidic_torsion_energy_validity_certified": False,
                "scope": "continuous invariant/3D-stereo diagnostics; no fitted validity gate and no complete physical-validity certificate"}


def _score_pool(rows, native, source_protein, native_protein, native_keys, fits, native_contacts, graph,
                source_to_native_indices=None):
    native = np.asarray(native)
    tree = cKDTree(source_protein)
    target_residues = [key[:4] for key in native_keys]
    if source_to_native_indices is None:
        source_to_native_indices = list(range(len(source_protein)))
    evaluated = []
    for pose in rows:
        full = np.asarray(pose["coordinates"])
        geometry = graph.diagnostics(full)
        distances, _ = tree.query(full)
        severe_pairs = sum(len(group) for group in tree.query_ball_point(full, 1.5))
        result = {key: pose[key] for key in ["rank", "vina_score", "candidate_class_rank", "candidate_original_rank", "box_id", "pooled_rank"]}
        result.update(geometry_diagnostics=geometry, full_parent_min_protein_distance_A=float(distances.min()),
                      full_parent_protein_pairs_below_1p5A=severe_pairs)
        for label, fit in fits.items():
            transformed = full @ fit["rotation"] + fit["translation"]
            options = [transformed[index] for index in graph.scored_maps]
            errors = [float(np.sqrt(np.mean(np.sum((candidate - native) ** 2, axis=1)))) for candidate in options]
            best = int(np.argmin(errors))
            result[label + "_RMSD_A"] = errors[best]
            result[label + "_chemical_atom_map_index"] = best
            # A proper common receptor transform preserves distances. Compute
            # contacts in the input frame to avoid transform-roundoff changing
            # a strict 4 A boundary, retaining the same chemical atom mapping.
            input_core = full[graph.scored_maps[best]]
            result[label + "_contacts"] = contact_metrics(input_core, source_protein, native_contacts,
                receptor_residues=target_residues, receptor_native_indices=source_to_native_indices,
                native_receptor_count=len(native_protein))
        evaluated.append(result)
    return evaluated


def aggregate_conditions(rows, scope, method, site_budget, alignment):
    selected = [row for row in rows if row["atomic_scope"] == scope and row["method"] == method]
    fields = [f"{rank}_success_{threshold}A" for rank in ("top1", "top5", "best_returned") for threshold in (1, 2, 3)]
    values = [row["site_budgets"][str(site_budget)][alignment] for row in selected]
    return {"eligible_complexes": len(selected), "complex_ids": [r["complex_id"] for r in selected],
            **{field: mean(v[field] for v in values) if values else None for field in fields},
            "note": "three seeds averaged within each complex first; equal complex weighting; no independent-family reliability interval"}


def generate_qualified_pose_report(root, *, allow_pending=False):
    root = Path(root).resolve()
    evidence = Evidence(root)
    config = evidence.json("benchmarks/blind_pose_experiment.json")
    freeze = evidence.json("results/blind_pose/freeze_receipt.json")
    if (evidence.sources["benchmarks/blind_pose_experiment.json"]["sha256"] != freeze["configuration_sha256"]
            or config.get("frozen_before_new_docking") is not True or freeze.get("frozen_before_new_docking") is not True
            or freeze.get("new_pose_outcomes_consulted") is not False):
        raise ValueError("Experiment configuration does not match its pre-execution freeze")
    ledger = evidence.json(config["atomic_reference_ledger"], freeze["atomic_reference_ledger_sha256"])
    if config["atomic_reference_ledger_sha256"] != freeze["atomic_reference_ledger_sha256"]:
        raise ValueError("Atomic eligibility changed after freezing")
    inventory = evidence.json("results/blind_pose/planning_inventory.json", freeze["planning_inventory_sha256"])
    oracle = evidence.json("results/blind_pose/oracle_controls.json", freeze["oracle_controls_sha256"])
    _check_atomic_sources(evidence, ledger)
    conditions, pending, denominator = _all_execution(evidence, config, freeze, inventory, oracle)
    if pending:
        if not allow_pending:
            raise ValueError(f"Frozen execution remains incomplete; no atomic accuracy report generated ({len(pending)} pending box/seed records)")
        return {"schema_version": 1, "status": "pending_incomplete_workflow", "metrics_eligible": False,
                "denominator": denominator, "pending": pending}, "# Qualified pose report pending\n\nThe frozen run denominator is incomplete. No atomic accuracy aggregate is reported.\n"
    if denominator["terminal_seed_runs"] != 156:
        raise ValueError("The full frozen seed denominator is not terminal")
    records = {row["case_id"]: row for row in ledger["cases"]}
    qualified = [row for row in ledger["cases"] if row["atomic_pose_eligible"]]
    if (len(qualified) != 6 or sum(row["full_ligand_eligible"] for row in qualified) != 2
            or sum(row["modeled_core_eligible"] for row in qualified) != 4):
        raise ValueError("Frozen full/core atomic eligibility denominator changed")
    contexts = {}
    with tempfile.TemporaryDirectory(prefix="qualified-pose-reference-") as temporary:
        temporary = Path(temporary)
        for record in qualified:
            cid = record["case_id"]
            native_keys, native_protein = _native_protein(evidence, record)
            native = np.asarray([atom["xyz_A"] for atom in record["scored_atom_map"]])
            controls = native_metric_controls(native, native_protein)
            restore_site_inputs(root, cid, "original", temporary / cid)
            keys, protein, _ = read_receptor(temporary / cid / "receptor.cif")
            input_manifest_path = Path("results/site_discovery") / cid / "inputs/manifest.json"
            input_manifest = evidence.json(input_manifest_path)
            for entry in input_manifest["files"]:
                if entry["revision"] == "original":
                    evidence.read(input_manifest_path.parent / entry["archive"], entry["archive_sha256"])
            fits, contacts, residues = _protein_alignment(keys, protein, native_keys, native_protein, native)
            native_index = {key: i for i, key in enumerate(native_keys)}
            old = _old_inputs(evidence, cid)
            mapping = graph = None
            if old is not None:
                reference = json.loads(old["reference.json"])
                mapping = qualified_atom_map(record, reference)
                graph = ParentGeometry(old["ligand.sdf"], reference, mapping)
            contexts[cid] = {"record": record, "native": native, "protein": protein, "native_protein": native_protein,
                             "native_keys": native_keys, "fits": fits, "contacts": contacts, "mapping": mapping,
                             "graph": graph, "controls": controls, "interface_residue_count": len(residues),
                             "source_to_native_indices": [native_index[key] for key in keys]}
    rows = []
    for condition in conditions:
        cid = condition["complex_id"]
        record = records.get(cid)
        eligible = bool(record and record["atomic_pose_eligible"])
        output = {key: condition[key] for key in ("condition_id", "complex_id", "state", "method", "preparation_status")}
        output.update(atomic_pose_eligible=eligible, atomic_scope=record["qualification"] if eligible else "ineligible_atomic_reference",
                      eligibility_reason=record["decision_reason"] if record else "AAV9 remains site-only; no qualified atomic reference in this frozen study",
                      prepared_box_count=len(condition["boxes"]), returned_pose_count=sum(len(p) for box in condition["boxes"] for p in box["poses_by_seed"].values()))
        if not eligible:
            output["atomic_accuracy"] = None
            rows.append(output)
            continue
        context = contexts[cid]
        if context["graph"] is None and any(box["poses_by_seed"] for box in condition["boxes"]):
            raise ValueError("Unsupported full parent graph unexpectedly produced evaluated poses")
        output.update(scored_atom_count=record["scored_atom_count"], parent_ligand_atom_count=record["released_ligand_heavy_atom_count"],
                      scored_prepared_residue_numbers=record["scored_prepared_residue_numbers"],
                      native_controls=context["controls"], native_interface_residue_count=context["interface_residue_count"],
                      atom_mapping=context["mapping"], chemical_scored_atom_map_count=len(context["graph"].scored_maps) if context["graph"] else None,
                      receptor_alignments={label: fit["metadata"] for label, fit in context["fits"].items()}, site_budgets={})
        for budget in config["candidate_budgets"]:
            evaluations = []
            for seed in config["seeds"]:
                pool = pooled_poses(condition["boxes"], seed, budget)
                scored = _score_pool(pool, context["native"], context["protein"], context["native_protein"],
                                     context["native_keys"], context["fits"], context["contacts"], context["graph"],
                                     context["source_to_native_indices"]) if pool else []
                evaluations.append({"seed": seed, "poses": scored,
                    **{label: accuracy_summary(scored, label + "_RMSD_A") for label in ("global", "local_interface")}})
            averaged = {}
            for label in ("global", "local_interface"):
                fields = [field for field in evaluations[0][label] if "success_" in field]
                averaged[label] = {field: mean(float(row[label][field]) for row in evaluations) for field in fields}
                averaged[label]["seeds_with_returned_pose"] = sum(bool(row[label]["returned_pose_count"]) for row in evaluations)
            output["site_budgets"][str(budget)] = {**averaged, "seed_evaluations": evaluations}
        rows.append(output)
    aggregates = {scope: {method: {str(budget): {alignment: aggregate_conditions(rows, scope, method, budget, alignment)
                     for alignment in ("global", "local_interface")} for budget in (1, 5)}
                     for method in config["methods"]} for scope in ("full_ligand", "modeled_core")}
    paired = []
    for scope in ("full_ligand", "modeled_core"):
        for method in ("fpocket", "p2rank"):
            cases = [row for row in rows if row["atomic_scope"] == scope and row["method"] == method]
            oracle_cases = {row["complex_id"]: row for row in rows if row["atomic_scope"] == scope and row["method"] == "oracle"}
            for budget in (1, 5):
                for alignment in ("global", "local_interface"):
                    for rank in ("top1", "top5", "best_returned"):
                        field = rank + "_success_2A"
                        contrasts = [{"complex_id": row["complex_id"],
                                      "oracle_minus_predicted": oracle_cases[row["complex_id"]]["site_budgets"]["1"][alignment][field] - row["site_budgets"][str(budget)][alignment][field]} for row in cases]
                        paired.append({"atomic_scope": scope, "method": method, "predicted_site_budget": budget,
                                       "oracle_site_budget": 1, "alignment": alignment, "pose_rank_budget": rank,
                                       "complexes": contrasts, "mean_oracle_minus_predicted": mean(c["oracle_minus_predicted"] for c in contrasts) if contrasts else None,
                                       "uncertainty": "descriptive paired development contrast; no population confidence interval"})
    summary = {"schema_version": 1, "status": "complete_development_evaluation", "metrics_eligible": True,
               "denominator": {**denominator, "atomic_eligible_complexes": 6, "full_model_complexes": 2,
                               "modeled_core_complexes": 4, "atomic_eligible_method_conditions": 18,
                               "atomic_reference_ineligible_method_conditions": 9},
               "conditions": rows, "aggregates": aggregates, "paired_contrasts": paired, "source_sha256": evidence.sources,
               "evaluator_sha256": sha256(Path(__file__)), "metrics_sha256": sha256(Path(__file__).with_name("metrics.py")),
               "comparison_code_sha256": sha256(Path(__file__).with_name("site_reproduction.py")),
               "limitations": ["Existing development systems, not untouched families or a confirmatory precision claim.",
                   "Core-only recovery is not complete-ligand recovery; full parent graphs were never truncated to fit the reference.",
                   "RMSD successes are geometric endpoints; ring-pucker and glycosidic torsion energetic validity remain uncertified.",
                   "All docking uses the frozen dry protein model and six-Angstrom whole-chain neighborhood selection; no general eight-Angstrom Vina coverage claim.",
                   "Top5 site search uses up to five boxes versus one oracle box; total compute is not matched.",
                   "No retained near-native pose does not prove absence among internal unretained sampling candidates."]}
    return summary, render_qualified_pose_markdown(summary)


def render_qualified_pose_markdown(summary):
    lines = ["# Qualified predicted-site and oracle-site pose results", "",
             "This development comparison separates complete modeled ligands from independently qualified modeled cores. No core result is called whole-ligand recovery.", "",
             f"The fixed workflow contains {summary['denominator']['eligible_conditions']} conditions and {summary['denominator']['scheduled_seed_runs']} scheduled seed runs. Atomic accuracy is eligible for six external complexes: two complete models and four modeled cores. AAV9 and 154L remain atomic-reference ineligible. Unsupported 1QFO chemistry remains zero workflow recovery within the eligible core denominator.", "",
             "Each cell below is the mean of three seed-specific success indicators per complex, then the equal-weight mean across complexes. The primary cutoff is receptor-aligned RMSD ≤2 Å. Site budget means searched candidate boxes; pose rank means energy-ranked returned poses within that site budget.", ""]
    for scope, title in [("full_ligand", "Complete modeled ligands"), ("modeled_core", "Qualified modeled cores")]:
        lines.extend([f"## {title}", "", "| Method | Site budget | Global Top1 | Global Top5 | Global best returned | Local Top1 | Local Top5 | Local best returned |", "|---|---:|---:|---:|---:|---:|---:|---:|"])
        for method in ("fpocket", "p2rank", "oracle"):
            for budget in ((1,) if method == "oracle" else (1, 5)):
                value = summary["aggregates"][scope][method][str(budget)]
                numbers = [value[label][field] for label in ("global", "local_interface") for field in ("top1_success_2A", "top5_success_2A", "best_returned_success_2A")]
                cells = ["undefined" if x is None else f"{x:.2%}" for x in numbers]
                lines.append(f"| {method} | {budget if method != 'oracle' else '1 (oracle)'} | " + " | ".join(cells) + " |")
        lines.append("")
    lines.extend(["## Per-complex results", "", "The following uses the five-site predicted budget and the one-box oracle control. Each percentage averages all three frozen seeds, retaining failed eligible runs as zero success.", "",
                  "| Complex | Qualified scope | Method | Global Top1 / Top5 / best returned | Local Top1 / Top5 / best returned |", "|---|---|---|---|---|"])
    for row in summary["conditions"]:
        if not row["atomic_pose_eligible"]:
            continue
        value = row["site_budgets"]["1" if row["method"] == "oracle" else "5"]
        cells = [" / ".join(f"{value[label][field]:.2%}" for field in ("top1_success_2A", "top5_success_2A", "best_returned_success_2A")) for label in ("global", "local_interface")]
        lines.append(f"| {row['complex_id']} | {row['atomic_scope']} ({row['scored_atom_count']} atoms) | {row['method']} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.extend(["## Interpretation and verification", "",
        "No independent ligand fit is used. Global alignment uses all corresponding protein backbone N/CA/C/O atoms. Local alignment uses backbone atoms of native qualified-core contact residues, solely during evaluation. Chemically valid full-parent automorphisms preserve source graph/stereochemistry; the same chosen map supplies RMSD and contacts. Native controls give zero RMSD, full native contacts and 5 Å error after the known [3,4,0] Å translation.", "",
        "Per-pose 3D stereochemistry, bond-length changes, ring-distance changes, linkage distances and full-parent/protein severe overlaps are recorded separately. These diagnostics do not certify ring pucker, glycosidic torsion energies or overall physical validity. Primary RMSD success is not silently converted into a pose-validity claim.", "",
        "All failures retain their eligible workflow denominator. Three seeds, capsid copies and repeated boxes are not independent biological observations. Full-model and core strata are small development samples; no family-level reliability interval or universal docking conclusion is earned.", "",
        "Inputs restore from checked-in preparation archives, while native coordinates enter only this evaluator after execution completeness is verified. Frozen configuration, reference ledger, original structural sources, input/pose hashes, commands and raw Vina log tables are checked. The summary retains every pooled candidate identity and source checksum.", "",
        "Machine-readable results: [qualified_pose_summary.json](../results/blind_pose/qualified_pose_summary.json). Atomic eligibility: [atomic_reference_ledger.json](../benchmarks/atomic_reference_ledger.json).", ""])
    return "\n".join(lines)


def write_qualified_pose_report(root, *, summary_output=None, markdown_output=None, check=False, allow_pending=False,
                                numeric_tolerance=0.0):
    if (isinstance(numeric_tolerance, bool) or not math.isfinite(numeric_tolerance)
            or not 0 <= numeric_tolerance <= 1e-8):
        raise ValueError("Numeric tolerance must be finite, explicit and at most1e-8")
    root = Path(root).resolve()
    summary, markdown = generate_qualified_pose_report(root, allow_pending=allow_pending)
    summary_path = Path(summary_output) if summary_output else root / "results/blind_pose/qualified_pose_summary.json"
    markdown_path = Path(markdown_output) if markdown_output else root / "docs/QUALIFIED_BLIND_POSE_RESULTS.md"
    for path, content in [(summary_path, json.dumps(summary, sort_keys=True, indent=2, allow_nan=False) + "\n"), (markdown_path, markdown)]:
        if check:
            if not path.is_file():
                raise ValueError(f"Qualified pose report is stale: {path}")
            if path.read_text() != content:
                if path != summary_path or numeric_tolerance == 0:
                    raise ValueError(f"Qualified pose report is stale: {path}")
                comparison = compare_site_values(json.loads(path.read_text()), summary,
                    numeric_tolerance=numeric_tolerance, derived_fields=_DERIVED_FLOAT_FIELDS)
                changed = {k: v for k, v in comparison["max_numeric_difference_by_field"].items() if v}
                print(f"Qualified pose JSON derived-float comparison: maximum difference {comparison['max_numeric_difference']:.17g}; "
                      f"explicit tolerance {numeric_tolerance:g}; differences by field {changed}")
        else:
            write_text_atomic(path, content)
    return summary_path, markdown_path
