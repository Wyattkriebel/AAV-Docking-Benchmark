"""Prospective, portable evaluation of the fixed six structural test cases.

An explicit pinned evaluation plan names every source archive and condition.
Every scheduled run is checked before any reference-based endpoint is computed.
Preparation and execution failures remain eligible; full ligands and modeled
cores stay separate. These purposively selected cases do not certify unseen
families, unseen training data, or population reliability.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from pathlib import Path
from statistics import mean
import tarfile
import tempfile

import numpy as np
from scipy.spatial import cKDTree

from .engine import read_poses
from .metrics import contact_metrics, contact_pairs, fit_receptor
from .provenance import sha256, write_text_atomic
from .qualified_pose import Evidence, ParentGeometry, accuracy_summary, native_metric_controls, pooled_poses, _pose_log
from .site_geometry import read_receptor
from .statistics import paired_family_bootstrap


CASE_IDS = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
METHODS = ("fpocket", "p2rank", "oracle")
FAMILIES = dict(zip(CASE_IDS, ("cyanovirin_N", "BambL_RSL_propeller", "CpMnBP1_SBP_II",
                             "norovirus_P_domain", "RbmC_beta_prism", "GH116")))
TERMINAL_CONDITION = {"prepared", "preparation_failed", "predictor_failed", "completed_no_candidates", "unsupported"}


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _key(value):
    chain, number, insertion, residue, atom = value
    return str(chain), int(number), str(insertion).strip(), str(residue), str(atom)


def _pin(evidence, record):
    return evidence.read(record["path"], record["sha256"])


def archive_files(evidence, record):
    """Read a pinned portable source tar without extracting arbitrary paths."""
    manifest = json.loads(_pin(evidence, record))
    archive = Path(record["path"]).parent / manifest["archive"]
    if Path(manifest["archive"]).is_absolute() or ".." in Path(manifest["archive"]).parts:
        raise ValueError("Source archive must remain beside its manifest")
    raw = evidence.read(archive, manifest["archive_sha256"])
    expected = {}
    for entry in manifest["files"]:
        name = entry["path"]
        if Path(name).is_absolute() or ".." in Path(name).parts or name in expected:
            raise ValueError("Unsafe or duplicate source archive member")
        expected[name] = entry
    actual = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        for member in archive.getmembers():
            if member.isdir():
                continue
            if not member.isfile() or member.name not in expected or member.name in actual:
                raise ValueError("Source archive has unexpected or duplicate members")
            content = archive.extractfile(member).read()
            entry = expected[member.name]
            if _digest(content) != entry["sha256"] or len(content) != entry["bytes"]:
                raise ValueError("Source archive member changed")
            actual[member.name] = content
    if set(actual) != set(expected):
        raise ValueError("Source archive is incomplete")
    return actual


def check_reference_sources(evidence, record, source_manifest_record=None):
    ledger = json.loads(_pin(evidence, record))
    if ledger.get("source_only") is not True or ledger.get("prediction_outcomes_consulted") is not False:
        raise ValueError("Reference ledger is not a source-only freeze")
    if ledger["cohort_denominator"] != 12:
        raise ValueError("Original twelve screened candidates must remain recorded")
    for key, value in ledger.items():
        if key.endswith("_archive"):
            hash_key = key.removesuffix("_archive") + "_sha256"
            if hash_key in ledger:
                evidence.read(value, ledger[hash_key])
    if source_manifest_record is not None and source_manifest_record["path"] != ledger["source_record"]:
        raise ValueError("Reference source-manifest identity changed")
    source_manifest = (json.loads(_pin(evidence, source_manifest_record)) if source_manifest_record is not None
                       else evidence.json(ledger["source_record"]))
    if source_manifest["ledger_sha256"] != record["sha256"]:
        raise ValueError("Reference source inventory and ledger disagree")
    for entry in source_manifest["archived_files"]:
        raw = evidence.read(entry["archive"], entry["sha256"])
        if len(raw) != entry["bytes"]:
            raise ValueError("Pinned reference archive size changed")
    cases = {c["case_id"]: c for c in ledger["cases"]}
    if len(cases) != 12 or any(cid not in cases for cid in CASE_IDS):
        raise ValueError("Fixed structural reference denominator changed")
    selected = {cid: cases[cid] for cid in CASE_IDS}
    if (sum(c["qualification"] == "full_ligand" for c in selected.values()) != 5
            or selected["4R9F"]["qualification"] != "modeled_core" or selected["4R9F"]["scored_atom_count"] != 11):
        raise ValueError("Frozen five-full/one-core qualification changed")
    return selected


def canonical_mapping(case, reference):
    """Require the complete parent graph and exactly one fixed canonical mask."""
    names, atoms = reference["names"], reference["canonical_atom_map"]
    if (len(names) != len(set(names)) or len(atoms) != len(names)
            or len(names) != case["inventory_entry"]["expected_heavy_atom_count"]):
        raise ValueError("Full-parent atom loss or duplicate ligand names")
    by_key = {}
    for atom in atoms:
        key = atom["canonical_atom_key"]
        if key in by_key or atom["name"] not in names:
            raise ValueError("Canonical parent atom map is not a bijection")
        by_key[key] = names.index(atom["name"])
    if "residues" in case["inventory_entry"]:
        expected_keys = {f"{r['branch_num']}:{r['ccd_id']}:{name}" for r in case["inventory_entry"]["residues"]
                         for name in r["expected_heavy_atoms_after_leaving_atom_removal"]}
        if set(by_key) != expected_keys:
            raise ValueError("Complete parent canonical identities differ from the pinned chemical graph")
    if {a["name"] for a in atoms} != set(names):
        raise ValueError("Canonical map omits full-parent identities")
    mask = [a["canonical_atom_key"] for a in case["scored_atom_map"]]
    if len(mask) != case["scored_atom_count"] or len(set(mask)) != len(mask) or not set(mask) <= set(by_key):
        raise ValueError("Qualified canonical mask is incomplete")
    copies = []
    for copy in case["qualified_copies"]:
        cmap = copy["scored_atom_map"]
        if ([a["canonical_atom_key"] for a in cmap] != mask or copy["scored_atom_count"] != len(mask)
                or len(cmap) != len(mask)):
            raise ValueError("Qualified copies changed the fixed canonical scored mask")
        xyz = np.asarray([a["xyz_A"] for a in cmap], float)
        if xyz.shape != (len(mask), 3) or not np.isfinite(xyz).all():
            raise ValueError("Qualified copy has incomplete native coordinates")
        copies.append({"copy_id": copy["copy_id"], "coordinates": xyz})
    if not copies or len({c["copy_id"] for c in copies}) != len(copies):
        raise ValueError("Qualified assembly copies are empty or duplicated")
    return {"indices": [by_key[k] for k in mask], "keys": mask, "copies": copies,
            "parent_atom_count": len(names), "names": [names[by_key[k]] for k in mask]}


def protein_fits(input_keys, input_xyz, protein_map, native, *, allow_missing_observed=False):
    """Fit only observed protein backbones; generated atoms are not truth."""
    target = {_key(a["output_key"]): np.asarray(a["assembly_xyz_A"], float) for a in protein_map["atoms"]}
    if len(target) != len(protein_map["atoms"]) or not np.isfinite(np.asarray(list(target.values()))).all():
        raise ValueError("Observed protein map has duplicate identities")
    source = dict(zip(input_keys, input_xyz))
    if len(source) != len(input_keys) or (not allow_missing_observed and not set(target) <= set(source)):
        raise ValueError("Input protein lost an observed source atom")
    native_keys = list(target)
    native_xyz = np.asarray([target[k] for k in native_keys])
    contacts = contact_pairs(native, native_xyz)
    residues = {native_keys[j][:4] for _, j in contacts}
    fits = {}
    for label in ("global", "local_interface"):
        keys = [k for k in native_keys if k in source and k[-1] in {"N", "CA", "C", "O"}
                and (label == "global" or k[:4] in residues)]
        a, b = np.asarray([source[k] for k in keys]), np.asarray([target[k] for k in keys])
        # Identical holo frames have an exact identity fit. Avoid injecting SVD
        # roundoff into strict distance/contact boundaries in that common case.
        if len(a) >= 3 and np.array_equal(a, b) and np.linalg.matrix_rank(a - a.mean(axis=0)) >= 2:
            rotation, translation = np.eye(3), np.zeros(3)
        else:
            rotation, translation = fit_receptor(a, b)
        fits[label] = {"rotation": rotation, "translation": translation,
                       "metadata": {"matched_backbone_atoms": len(keys), "backbone_RMSD_A": float(np.sqrt(np.mean(np.sum((a @ rotation + translation - b) ** 2, axis=1)))),
                                    "rotation": rotation.tolist(), "translation": translation.tolist(),
                                    "selection": "all observed corresponding N/CA/C/O" if label == "global" else "observed N/CA/C/O of qualified-native <4A interface residues; evaluation only"}}
    universe = {k: i for i, k in enumerate(native_keys)}
    for key in input_keys:
        universe.setdefault(key, len(universe))
    residues_all = [None] * len(universe)
    for key, index in universe.items():
        residues_all[index] = key[:4]
    return {"fits": fits, "contacts": contacts, "native_keys": native_keys, "native_protein": native_xyz,
            "input_to_universe": [universe[k] for k in input_keys], "universe_residues": residues_all,
            "generated_protein_atom_count": len(universe) - len(native_keys),
            "missing_observed_protein_atom_count": len(set(target) - set(source))}


def score_pool(rows, mapping, input_keys, protein, protein_map, graph, *, allow_missing_observed=False):
    """Select qualified-copy correspondence within each fixed metric, never ranks."""
    contexts = [(copy, protein_fits(input_keys, protein, protein_map, copy["coordinates"], allow_missing_observed=allow_missing_observed)) for copy in mapping["copies"]]
    allowed_maps = [m for m in graph.scored_maps if set(m) == set(mapping["indices"])]
    if not allowed_maps:
        raise ValueError("No chemical automorphism preserves the fixed canonical mask")
    tree = cKDTree(protein)
    evaluated = []
    for pose in rows:
        full = np.asarray(pose["coordinates"], float)
        geometry = graph.diagnostics(full)
        distances, _ = tree.query(full)
        row = {k: pose[k] for k in ("rank", "vina_score", "candidate_class_rank", "candidate_original_rank", "box_id", "pooled_rank")}
        row.update(geometry_diagnostics=geometry, full_parent_min_protein_distance_A=float(distances.min()),
                   full_parent_protein_pairs_below_1p5A=sum(len(g) for g in tree.query_ball_point(full, 1.5)))
        global_fit = contexts[0][1]["fits"]["global"]
        world_pose = full @ global_fit["rotation"] + global_fit["translation"]
        observed_assembly_tree = cKDTree(contexts[0][1]["native_protein"])
        observed_distances, _ = observed_assembly_tree.query(world_pose)
        row["full_observed_assembly_min_protein_distance_A"] = float(observed_distances.min())
        row["full_observed_assembly_protein_pairs_below_1p5A"] = sum(len(g) for g in observed_assembly_tree.query_ball_point(world_pose, 1.5))
        for label in ("global", "local_interface"):
            options = []
            for copy, context in contexts:
                fit = context["fits"][label]
                transformed = full @ fit["rotation"] + fit["translation"]
                for index, atom_map in enumerate(allowed_maps):
                    rmsd = float(np.sqrt(np.mean(np.sum((transformed[atom_map] - copy["coordinates"]) ** 2, axis=1))))
                    options.append((rmsd, copy["copy_id"], index, context, atom_map))
            rmsd, copy_id, map_index, context, atom_map = min(options, key=lambda x: x[:3])
            row[label + "_RMSD_A"], row[label + "_qualified_copy_id"] = rmsd, copy_id
            row[label + "_chemical_atom_map_index"] = map_index
            # A shared proper protein transform preserves pair distances. Use
            # raw input-frame pairs to avoid roundoff at the strict4A boundary.
            row[label + "_contacts"] = contact_metrics(full[atom_map], protein, context["contacts"],
                receptor_residues=context["universe_residues"], receptor_native_indices=context["input_to_universe"],
                native_receptor_count=len(context["universe_residues"]))
        evaluated.append(row)
    return evaluated, [{"copy_id": copy["copy_id"], "alignments": {k: v["metadata"] for k, v in context["fits"].items()},
                        "native_controls": native_metric_controls(copy["coordinates"], context["native_protein"]),
                        "generated_atoms_excluded_from_reference": context["generated_protein_atom_count"],
                        "missing_observed_input_protein_atoms": context["missing_observed_protein_atom_count"]} for copy, context in contexts]


def site_summary(candidates, copies, *, box_sizes=None, rotation=None, translation=None):
    """DCA and coverage against any qualified copy; exact fixed rank order."""
    rotation = np.eye(3) if rotation is None else np.asarray(rotation)
    translation = np.zeros(3) if translation is None else np.asarray(translation)
    rows = []
    for candidate in candidates:
        center = np.asarray(candidate["center_A"]) @ rotation + translation
        distances = {c["copy_id"]: float(np.linalg.norm(c["coordinates"] - center, axis=1).min()) for c in copies}
        rows.append({"candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"],
                     "center_to_native_atom_A": min(distances.values()), "qualified_copy_distances_A": distances})
    output = {"definition": "predicted center to nearest qualified native heavy atom; global protein fit only", "candidates": rows,
              "qualified_copy_count": len(copies), "primary_threshold_A": 4, "sensitivity_thresholds_A": [3, 5]}
    for budget in (1, 5):
        selected = [r for r in rows if r["candidate_class_rank"] <= budget]
        values = {"available_candidates": len(selected), "minimum_DCA_A": min((r["center_to_native_atom_A"] for r in selected), default=None)}
        for threshold in (3, 4, 5):
            covered = sum(any(r["qualified_copy_distances_A"][c["copy_id"]] <= threshold for r in selected) for c in copies)
            values[f"success_{threshold}A"] = bool(covered)
            values[f"qualified_copy_coverage_{threshold}A"] = covered / len(copies)
        if box_sizes is not None:
            coverage = []
            for copy in copies:
                fractions = []
                for candidate in candidates:
                    if candidate["rank"] > budget:
                        continue
                    # Cubes are defined in input coordinates, so invert only the
                    # common global protein transform for the coverage check.
                    native_input = (copy["coordinates"] - translation) @ rotation.T
                    inside = np.all(np.abs(native_input - np.asarray(candidate["center_A"])) <= np.asarray(box_sizes) / 2, axis=1)
                    fractions.append(float(inside.mean()))
                coverage.append(max(fractions, default=0.0))
            values["qualified_copy_max_box_atom_coverage"] = coverage
            values["any_qualified_copy_entirely_in_box"] = any(v == 1 for v in coverage)
        output[str(budget)] = values
    return output


def _command(row, prep, config, preparation_path, execution_root):
    command = row.get("command", [])
    if not command or len(command) % 2 != 1 or any(not isinstance(x, str) for x in command):
        raise ValueError("Missing actual Vina invocation")
    options = dict(zip(command[1::2], command[2::2]))
    if len(options) * 2 != len(command) - 1:
        raise ValueError("Duplicate Vina command flags")
    expected = {"--scoring": "vina", "--seed": row["seed"], "--cpu": config["cpu"], "--exhaustiveness": config["exhaustiveness"],
                "--num_modes": config["num_modes"], "--min_rmsd": config["min_rmsd_A"], "--energy_range": config["energy_range_kcal_mol"]}
    for axis, center, side in zip("xyz", prep["box_center"], prep["box_size"]):
        expected["--center_" + axis], expected["--size_" + axis] = center, side
    for flag, value in expected.items():
        actual = options.get(flag)
        if actual is None or (isinstance(value, str) and actual != value) or (not isinstance(value, str) and (not math.isfinite(float(actual)) or float(actual) != float(value))):
            raise ValueError("Actual invocation changed a frozen Vina setting")
    suffixes = {"--receptor": preparation_path / "neighborhood.pdbqt", "--ligand": preparation_path / "ligand.pdbqt",
                "--out": execution_root / row["pose_file"]}
    if set(options) != set(expected) | set(suffixes):
        raise ValueError("Unexpected Vina invocation flags")
    for flag, path in suffixes.items():
        if Path(options[flag]).parts[-len(path.parts):] != path.parts:
            raise ValueError("Actual invocation refers to a different condition/box")


def read_workflows(evidence, profile, config, sources, *, prediction_sources=None, prediction_config_sha256=None):
    """Validate every condition before any reference-based accuracy is available."""
    declared, methods = config.get("primary_case_ids", config.get("case_ids")), config["methods"]
    cases = profile.get("case_ids", list(CASE_IDS))
    is_subset_control = profile["role"] == "assembly_first_chain_sensitivity"
    if (not cases or len(cases) != len(set(cases)) or not set(cases) <= set(CASE_IDS)
            or (not is_subset_control and cases != list(CASE_IDS)) or not set(cases) <= set(declared)
            or methods != list(METHODS)):
        raise ValueError("Primary profiles retain all six cases; explicit assembly controls retain their declared subset and all three methods")
    prediction_sources = sources if prediction_sources is None else prediction_sources
    prediction_config_sha256 = profile["configuration"]["sha256"] if prediction_config_sha256 is None else prediction_config_sha256
    expected = {(cid, method) for cid in cases for method in METHODS}
    descriptors = {(x["case_id"], x["method"]): x for x in profile["conditions"]}
    if set(descriptors) != expected or len(profile["conditions"]) != len(expected):
        raise ValueError("Explicit method-condition denominator must retain every declared identity")
    vina = evidence.json("benchmarks/vina.lock.json")
    executables = {v["sha256"] for v in vina["assets"].values()}
    conditions, pending = [], []
    counts = {"eligible_conditions": len(expected), "prepared_boxes": 0, "scheduled_seed_runs": 0,
              "terminal_seed_runs": 0, "successful_seed_runs": 0, "preparation_failed_boxes": 0,
              "conditions_without_prepared_boxes": 0, "predictor_failed_conditions": 0,
              "measured_terminal_cpu_seconds": 0.0, "measured_terminal_wall_seconds": 0.0,
              "previous_attempts_with_unknown_cost": 0}
    config_sha = profile["configuration"]["sha256"]
    for cid, method in sorted(expected):
        descriptor = descriptors[cid, method]
        path = Path(descriptor["preparation"])
        if not (evidence.root / path).is_file():
            pending.append(f"{cid}/{method}:missing_condition")
            continue
        condition = evidence.json(path)
        if (condition.get("case_id", condition.get("complex_id")) != cid or condition["method"] != method
                or condition["configuration_sha256"] != config_sha
                or condition.get("fixed_denominator_member", condition.get("eligible_condition")) is not True):
            raise ValueError("Condition changed its frozen scientific identity")
        if condition["status"] not in TERMINAL_CONDITION:
            pending.append(f"{cid}/{method}:condition_{condition['status']}")
        selection = None
        if method != "oracle":
            selection_path = Path(descriptor["prediction_selection"])
            if not (evidence.root / selection_path).is_file():
                pending.append(f"{cid}/{method}:missing_prediction_selection")
            else:
                selection = evidence.json(selection_path)
                record = selection["prediction_record"]
                predictor = evidence.json(selection_path.parent / record["filename"], record["sha256"])
                if selection["predictor_status"] != predictor["status"]:
                    raise ValueError("Selection changed the raw predictor execution status")
                if (selection["case_id"] != cid or selection["method"] != method or selection.get("native_coordinate_access") is not False
                        or selection["prediction_input_sha256"] != _digest(prediction_sources[f"{cid}/receptor.cif"])
                        or selection["configuration_sha256"] != prediction_config_sha256):
                    raise ValueError("Prediction selection changed its case or source-only input")
                if condition.get("candidate_selection") != selection:
                    raise ValueError("Preparation changed the selected predictor candidates")
                if selection["predictor_status"] not in {"ok", "completed"}:
                    counts["predictor_failed_conditions"] += 1
                if "primary_candidate_selection" in descriptor:
                    original = json.loads(_pin(evidence, descriptor["primary_candidate_selection"]))
                    fields = ("candidates", "candidate_mapping", "prediction_record", "prediction_input_sha256")
                    if any(original[k] != selection[k] for k in fields):
                        raise ValueError("Modeled sensitivity changed original centers/ranks/predictions")
        jobs = condition["jobs"]
        ranks = [j["candidate_class_rank"] for j in jobs]
        if ranks != list(range(1, len(jobs) + 1)) or len(jobs) > (1 if method == "oracle" else 5):
            raise ValueError("Candidate classes are missing, reordered or duplicated")
        if selection is not None:
            if len(jobs) != len(selection["candidates"]):
                raise ValueError("A failed preparation disappeared from the candidate denominator")
            for job, candidate in zip(jobs, selection["candidates"]):
                if job["center_A"] != candidate["center_A"] or job["candidate_original_rank"] != candidate["original_rank"]:
                    raise ValueError("Prepared candidate differs from its frozen predictor rank/center")
        output = {"case_id": cid, "method": method, "preparation_status": condition["status"], "selection": selection, "boxes": []}
        for job in jobs:
            box_id = job["id"]
            if box_id != f"box{job['candidate_class_rank']:02d}" or job["prepared_directory"] != box_id:
                raise ValueError("Unsafe candidate-box directory")
            box = {"id": box_id, "candidate_class_rank": job["candidate_class_rank"], "candidate_original_rank": job["candidate_original_rank"],
                   "center_A": job["center_A"], "status": job["status"], "poses_by_seed": {}}
            output["boxes"].append(box)
            if job["status"] != "prepared":
                if job["status"] not in {"preparation_failed", "unsupported"}:
                    pending.append(f"{cid}/{method}/{box_id}:preparation_{job['status']}")
                counts["preparation_failed_boxes"] += 1
                attempt = path.parent / box_id / "preparation_attempt.json"
                evidence.read(attempt)
                continue
            counts["prepared_boxes"] += 1
            counts["scheduled_seed_runs"] += len(config["seeds"])
            directory, runs_dir = path.parent / box_id, Path(descriptor["execution_root"]) / box_id
            raw_prep = evidence.read(directory / "preparation.json", job["preparation_sha256"])
            prep, reference = json.loads(raw_prep), evidence.json(directory / "reference.json")
            chemistry = json.loads(sources[f"{cid}/chemistry.json"])
            omitted = prep["neighborhood"].get("closest_omitted_atom_to_box_A")
            if (prep["case_id"] != cid or prep["box_center"] != job["center_A"] or prep["configuration_sha256"] != config_sha
                    or reference.get("pose_reference_eligible") is not False or "coordinates" in reference
                    or prep["neighborhood"]["interaction_cutoff_A"] != 8 or prep["box_size"] != chemistry["box_size"]
                    or (omitted is not None and omitted <= 8)):
                raise ValueError("Prepared box changed its source-only scientific contract")
            source_map = json.loads(sources[f"{cid}/ligand_atom_map.json"])["atoms"]
            if reference["canonical_atom_map"] != source_map or reference["names"] != [a["name"] for a in source_map]:
                raise ValueError("Prepared reference changed the complete canonical parent map")
            expected_fp = {"config_sha256": config_sha, "preparation_sha256": _digest(raw_prep),
                "reference_sha256": _digest(evidence.read(directory / "reference.json")),
                "ligand_sha256": _digest(evidence.read(directory / "ligand.pdbqt")),
                "receptor_sha256": _digest(evidence.read(directory / "neighborhood.pdbqt")),
                "receptor_pdb_sha256": _digest(evidence.read(directory / "neighborhood.pdb")), "full_assembly_sha256": None}
            if not (evidence.root / runs_dir / "runs.json").is_file():
                pending.append(f"{cid}/{method}/{box_id}:missing_runs")
                continue
            evidence.read(runs_dir / "configuration.json", config_sha)
            rows = evidence.json(runs_dir / "runs.json")
            seeds = [r["seed"] for r in rows]
            if len(seeds) != len(set(seeds)) or set(seeds) - set(config["seeds"]):
                raise ValueError("Unknown or duplicate frozen seed")
            if set(seeds) != set(config["seeds"]):
                pending.append(f"{cid}/{method}/{box_id}:missing_seeds")
            box["runs"] = rows
            for row in rows:
                if any(row["fingerprint"].get(k) != v for k, v in expected_fp.items()) or row["fingerprint"].get("executable_sha256") not in executables:
                    raise ValueError("Raw run input/configuration/executable fingerprints changed")
                regime = "true-site-local" if method == "oracle" else "predicted-site-local"
                if row.get("case_id") != cid or row.get("receptor") != "neighborhood" or row.get("information_regime") != regime or row.get("vina_version") != "AutoDock Vina v1.2.7":
                    raise ValueError("Raw run identity/information regime changed")
                _command(row, prep, config, directory, runs_dir)
                if row["status"] in {"running", "pending"}:
                    pending.append(f"{cid}/{method}/{box_id}/{row['seed']}:running")
                    continue
                if row["status"] not in {"ok", "failed"}:
                    raise ValueError("Unknown terminal run state")
                counts["terminal_seed_runs"] += 1
                for field in ("cpu_seconds", "wall_seconds"):
                    value = row.get(field)
                    if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                        raise ValueError("Terminal execution lacks measured finite resource use")
                    counts["measured_terminal_" + field] += value
                counts["previous_attempts_with_unknown_cost"] += sum(
                    any(field not in attempt for field in ("cpu_seconds", "wall_seconds")) for attempt in row.get("previous_attempts", []))
                for name in (row["pose_file"], row["log_file"]):
                    if Path(name).name != name or name in {".", ".."}:
                        raise ValueError("Unsafe raw result filename")
                log = evidence.read(runs_dir / row["log_file"])
                if row["status"] == "ok":
                    if row.get("exit_code") != 0:
                        raise ValueError("Successful raw run has nonzero exit status")
                    raw = evidence.read(runs_dir / row["pose_file"], row["pose_sha256"])
                    with tempfile.TemporaryDirectory() as tmp:
                        temporary = Path(tmp) / "pose.pdbqt"; temporary.write_bytes(raw)
                        poses = read_poses(temporary, reference["names"])
                    if len(poses) != row["returned_poses"] or not 1 <= len(poses) <= config["num_modes"]:
                        raise ValueError("Raw pose count changed")
                    _pose_log(log, poses)
                    box["poses_by_seed"][row["seed"]] = poses
                    counts["successful_seed_runs"] += 1
        if not any(j["status"] == "prepared" for j in jobs):
            counts["conditions_without_prepared_boxes"] += 1
        conditions.append(output)
    counts["failed_seed_runs"] = counts["terminal_seed_runs"] - counts["successful_seed_runs"]
    counts["cases_with_prepared_boxes"] = len({r["case_id"] for r in conditions if any(b["status"] == "prepared" for b in r["boxes"])})
    counts["fixed_case_denominator"] = len(cases)
    counts["pending_records"] = len(pending)
    return conditions, sorted(pending), counts


def _source_context(source, cid, case, *, evaluation_source=None):
    atoms = json.loads(source[f"{cid}/ligand_atom_map.json"])["atoms"]
    reference = {"names": [a["name"] for a in atoms], "canonical_atom_map": atoms}
    mapping = canonical_mapping(case, reference)
    converted = {"names": reference["names"], "atom_map": [{**a, "source_key": ["L", int(a["branch_num"]), "", a["ccd_id"], a["atom_name"]]} for a in atoms]}
    graph = ParentGeometry(source[f"{cid}/ligand.sdf"], converted, mapping)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "receptor.cif"; path.write_bytes(source[f"{cid}/receptor.cif"])
        keys, xyz, _ = read_receptor(path)
    observed = source if evaluation_source is None else evaluation_source
    return mapping, graph, keys, xyz, json.loads(observed[f"{cid}/protein_atom_map.json"])


def conditional_statistics(rows, case_to_family, seeds):
    """Per-case and family-cluster descriptions, conditional on this fixed sample."""
    output = []
    for scope in ("full_ligand", "modeled_core"):
        selected = [r for r in rows if r["atomic_scope"] == scope]
        if not selected:
            continue
        for budget in (1, 5):
            for alignment in ("global", "local_interface"):
                for rank in ("top1", "top5", "best_returned"):
                    for threshold in (1, 2, 3):
                        metric = f"{rank}_success_{threshold}A"
                        records = [{"family_id": case_to_family[r["case_id"]], "complex_id": r["case_id"], "method": r["method"],
                                    "seed": s["seed"], "status": s["site_budgets"][str(budget)]["status"], "success": s["site_budgets"][str(budget)][alignment][metric]}
                                   for r in selected for s in r["seeds"]]
                        for method in ("fpocket", "p2rank"):
                            value = paired_family_bootstrap(records, "oracle", method, n_resamples=10000, seed=20261103, expected_seeds=seeds)
                            value.update(atomic_scope=scope, candidate_budget=budget, alignment=alignment, pose_rank=rank, threshold_A=threshold,
                                interpretation="Conditional descriptive resampling of purposively selected source-family clusters; no population, unseen-family or training-unseen reliability claim")
                            output.append(value)
    return output


def macro_summaries(rows, case_to_family):
    """Equal seeds within case and equal cases within source-family cluster."""
    summaries = []
    fields = [f"{rank}_success_{threshold}A" for rank in ("top1", "top5", "best_returned") for threshold in (1, 2, 3)]
    for scope in ("full_ligand", "modeled_core"):
        for method in METHODS:
            selected = [r for r in rows if r["atomic_scope"] == scope and r["method"] == method]
            if not selected:
                continue
            for budget in (1, 5):
                for alignment in ("global", "local_interface"):
                    per_case = [{"case_id": r["case_id"], "family_id": case_to_family[r["case_id"]],
                                 **{field: mean(s["site_budgets"][str(budget)][alignment][field] for s in r["seeds"]) for field in fields}}
                                for r in selected]
                    per_family = [{"family_id": family, "case_ids": [r["case_id"] for r in per_case if r["family_id"] == family],
                                   **{field: mean(r[field] for r in per_case if r["family_id"] == family) for field in fields}}
                                  for family in sorted({r["family_id"] for r in per_case})]
                    summaries.append({"atomic_scope": scope, "method": method, "candidate_budget": budget, "alignment": alignment,
                        "eligible_cases": len(per_case), "source_family_clusters": len(per_family), "per_case": per_case, "per_family": per_family,
                        "case_macro": {field: mean(r[field] for r in per_case) for field in fields},
                        "family_macro": {field: mean(r[field] for r in per_family) for field in fields}})
    return summaries


def site_statistics(rows, case_to_family):
    comparisons = []
    for scope in ("full_ligand", "modeled_core"):
        selected = [r for r in rows if r["atomic_scope"] == scope and r["method"] != "oracle"]
        if not selected:
            continue
        for budget in (1, 5):
            for threshold in (3, 4, 5):
                records = [{"family_id": case_to_family[r["case_id"]], "complex_id": r["case_id"], "method": r["method"], "seed": 0,
                            "status": "ok" if r["site"][str(budget)]["available_candidates"] else "failed",
                            "success": r["site"][str(budget)][f"success_{threshold}A"]} for r in selected]
                result = paired_family_bootstrap(records, "fpocket", "p2rank", n_resamples=10000, seed=20261103, expected_seeds=[0])
                result.update(atomic_scope=scope, candidate_budget=budget, threshold_A=threshold,
                    interpretation="Conditional source-family cluster resampling of the fixed purposive structural sample; one site workflow per case, no pseudo-replicated docking seeds")
                comparisons.append(result)
    return comparisons


def profile_sensitivity_statistics(profiles, case_to_family):
    """Paired modeled-minus-strict contrast, always separate from primary arms."""
    primary = [p for p in profiles if p["role"] == "strict_primary"]
    modeled = [p for p in profiles if p["role"] == "modeled_sidechain_sensitivity"]
    if not modeled:
        return []
    if len(primary) != 1:
        raise ValueError("Modeled sensitivity needs exactly one strict primary profile")
    results = []
    baseline = primary[0]
    for sensitivity in modeled:
        for method in METHODS:
            for scope in ("full_ligand", "modeled_core"):
                for budget in (1, 5):
                    for alignment in ("global", "local_interface"):
                        records = []
                        for profile in (baseline, sensitivity):
                            for row in profile["conditions"]:
                                if row["method"] != method or row["atomic_scope"] != scope:
                                    continue
                                for seed in row["seeds"]:
                                    value = seed["site_budgets"][str(budget)]
                                    records.append({"family_id": case_to_family[row["case_id"]], "complex_id": row["case_id"], "method": profile["id"],
                                        "seed": seed["seed"], "status": value["status"], "success": value[alignment]["top1_success_2A"]})
                        if not records:
                            continue
                        seeds = sorted({r["seed"] for r in records})
                        result = paired_family_bootstrap(records, sensitivity["id"], baseline["id"], n_resamples=10000, seed=20261103, expected_seeds=seeds)
                        result.update(method=method, atomic_scope=scope, candidate_budget=budget, alignment=alignment,
                            pose_rank="top1", threshold_A=2, interpretation="Separate modeled-sidechain-minus-strict sensitivity; same original candidate centers and source-qualified reference; conditional fixed-cohort interval")
                        results.append(result)
    return results


def generate_report(root, plan_path, freeze_path, *, allow_pending=False):
    evidence = Evidence(root)
    plan = evidence.json(plan_path)
    freeze = evidence.json(freeze_path)
    if (freeze.get("evaluation_plan_sha256") != _digest(evidence.read(plan_path))
            or freeze.get("frozen_before_test_pose_accuracy") is not True or freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Evaluation requires its explicit pre-accuracy plan freeze")
    if plan["case_ids"] != list(CASE_IDS) or set(plan["case_to_family"]) != set(CASE_IDS):
        raise ValueError("Evaluation plan changed the fixed six-case/family denominator")
    if plan["statistics"] != {"bootstrap_replicates": 10000, "bootstrap_seed": 20261103, "confidence": 0.95, "interpretation": "conditional_purposive_sample"}:
        raise ValueError("Predeclared conditional family bootstrap changed")
    cases = check_reference_sources(evidence, plan["source_ledger"], plan["reference_source_manifest"])
    _pin(evidence, plan["family_audit"])
    prediction_sources = prediction_config_sha = None
    if "site_prediction" in plan:
        prediction = plan["site_prediction"]
        json.loads(_pin(evidence, prediction["configuration"]))
        prediction_freeze = json.loads(_pin(evidence, prediction["freeze_receipt"]))
        prediction_config_sha = prediction["configuration"]["sha256"]
        if (prediction_freeze.get("configuration_sha256") != prediction_config_sha
                or prediction_freeze.get("frozen_before_new_predictions") is not True):
            raise ValueError("Shared site predictions lack their pre-prediction freeze")
        prediction_sources = archive_files(evidence, prediction["source_archive_manifest"])
    loaded, pending = [], []
    if len({p["id"] for p in plan["profiles"]}) != len(plan["profiles"]) or not plan["profiles"]:
        raise ValueError("Profile identities must be unique and nonempty")
    for profile in plan["profiles"]:
        config = json.loads(_pin(evidence, profile["configuration"]))
        execution_freeze = json.loads(_pin(evidence, profile["freeze_receipt"]))
        if execution_freeze.get("configuration_sha256") != profile["configuration"]["sha256"] or execution_freeze.get("frozen_before_new_docking") is not True:
            raise ValueError("Scientific profile lacks its pre-execution configuration freeze")
        sources = archive_files(evidence, profile["source_archive_manifest"])
        source_receipt = json.loads(_pin(evidence, profile["source_preparation_receipt"]))
        if source_receipt.get("source_ledger_sha256", plan["source_ledger"]["sha256"]) != plan["source_ledger"]["sha256"]:
            raise ValueError("Profile source preparation changed the qualified reference ledger")
        if profile["role"] == "modeled_sidechain_sensitivity" and any(c["method"] != "oracle" and "primary_candidate_selection" not in c for c in profile["conditions"]):
            raise ValueError("Modeled sensitivity must pin all original primary candidate selections")
        conditions, unfinished, counts = read_workflows(evidence, profile, config, sources,
            prediction_sources=prediction_sources, prediction_config_sha256=prediction_config_sha)
        pending.extend(f"{profile['id']}/{s}" for s in unfinished)
        loaded.append((profile, config, sources, conditions, counts))
    strict_sources = [s for p, _, s, _, _ in loaded if p["role"] == "strict_primary"]
    for profile, _, sources, _, _ in loaded:
        if profile["role"] == "modeled_sidechain_sensitivity":
            if len(strict_sources) != 1:
                raise ValueError("Modeled sensitivity must pair with exactly one strict source profile")
            for cid in CASE_IDS:
                for name in ("ligand_atom_map.json", "ligand.sdf", "ligand.pdbqt", "protein_atom_map.json"):
                    if sources[f"{cid}/{name}"] != strict_sources[0][f"{cid}/{name}"]:
                        raise ValueError("Modeled sensitivity changed independent ligand or observed experimental protein map")
    # No native coordinates have been fitted and no site/pose metric computed.
    if pending:
        if not allow_pending:
            raise ValueError(f"Frozen structural execution remains incomplete ({len(pending)} pending records); no accuracy report")
        return {"schema_version": 1, "status": "pending_incomplete_workflow", "metrics_eligible": False, "pending": pending,
                "profiles": [{"id": p["id"], "counts": counts} for p, _, _, _, counts in loaded], "source_sha256_inventory": evidence.sources}
    profiles = []
    for profile, config, sources, conditions, counts in loaded:
        if counts["terminal_seed_runs"] != counts["scheduled_seed_runs"]:
            raise ValueError("The full scheduled seed denominator is not terminal")
        control = profile["role"] == "assembly_first_chain_sensitivity"
        evaluation_source = None
        if control:
            if "evaluation_source_archive_manifest" not in profile:
                raise ValueError("First-chain assembly sensitivity must pin the original full observed-assembly reference source")
            evaluation_source = archive_files(evidence, profile["evaluation_source_archive_manifest"])
        contexts = {cid: _source_context(sources, cid, cases[cid], evaluation_source=evaluation_source) for cid in profile.get("case_ids", CASE_IDS)}
        rows = []
        for condition in conditions:
            cid = condition["case_id"]
            mapping, graph, keys, protein, protein_map = contexts[cid]
            seeds, controls = [], None
            for seed in config["seeds"]:
                budgets = {}
                for budget in (1, 5):
                    pool = pooled_poses(condition["boxes"], seed, budget)
                    scored, controls = score_pool(pool, mapping, keys, protein, protein_map, graph, allow_missing_observed=control)
                    budgets[str(budget)] = {"status": "ok" if scored else "failed" if any(b["status"] == "prepared" and b["candidate_class_rank"] <= budget for b in condition["boxes"]) else "unsupported",
                                            "global": accuracy_summary(scored, "global_RMSD_A"),
                                            "local_interface": accuracy_summary(scored, "local_interface_RMSD_A"), "poses": scored}
                seeds.append({"seed": seed, "site_budgets": budgets})
            first_context = protein_fits(keys, protein, protein_map, mapping["copies"][0]["coordinates"], allow_missing_observed=control)
            fit = first_context["fits"]["global"]
            candidates = condition["selection"]["candidates"] if condition["selection"] is not None else [
                {"rank": box["candidate_class_rank"], "original_rank": box["candidate_original_rank"], "center_A": box["center_A"]} for box in condition["boxes"]]
            chemistry = json.loads(sources[f"{cid}/chemistry.json"])
            site = site_summary(candidates, mapping["copies"], box_sizes=chemistry["box_size"], rotation=fit["rotation"], translation=fit["translation"])
            rows.append({"case_id": cid, "method": condition["method"], "atomic_scope": cases[cid]["qualification"],
                         "scored_atoms": cases[cid]["scored_atom_count"], "parent_atoms": mapping["parent_atom_count"],
                         "preparation_status": condition["preparation_status"], "site": site, "native_controls": controls, "seeds": seeds})
        profiles.append({"id": profile["id"], "role": profile["role"], "counts": counts, "conditions": rows,
                         "scope_summaries": macro_summaries(rows, plan["case_to_family"]),
                         "paired_family_statistics": conditional_statistics(rows, plan["case_to_family"], config["seeds"]),
                         "paired_site_family_statistics": site_statistics(rows, plan["case_to_family"])})
    return {"schema_version": 1, "status": "complete", "metrics_eligible": True, "screened_cases": 12, "structural_cases": 6,
            "full_ligand_cases": 5, "modeled_core_cases": 1, "profiles": profiles, "case_to_family": plan["case_to_family"],
            "modeled_sidechain_sensitivity_contrasts": profile_sensitivity_statistics(profiles, plan["case_to_family"]),
            "source_sha256_inventory": evidence.sources,
            "code_sha256": {str(Path(__file__).name): sha256(__file__), "qualified_pose.py": sha256(Path(__file__).with_name("qualified_pose.py")),
                             "metrics.py": sha256(Path(__file__).with_name("metrics.py")), "statistics.py": sha256(Path(__file__).with_name("statistics.py"))},
            "limitations": ["Purposive structural test conditional on source availability and qualification; no population or model-training-unseen claim.",
                "Four strict source protein parameterization failures remain eligible; modeled sidechain completion is a separate sensitivity.",
                "Top5 candidate boxes use greater total compute than one oracle box; comparisons are matched per box, not total CPU.",
                "Native interfaces enter only evaluation protein fits; no independent ligand fit.",
                "Generated protein atoms never become experimental reference truth. Ring/stereo diagnostics do not certify complete physical validity."]}


def render_markdown(report):
    if not report["metrics_eligible"]:
        return "# Structural test pending\n\nThe frozen execution denominator is incomplete. No site or pose accuracy is reported.\n"
    text = ["# Fixed structural test", "", "Five full-ligand cases and the 11-atom 4R9F core are reported separately. All six cases remain eligible in every method/profile. These results do not certify model-unseen families or population reliability.", ""]
    for profile in report["profiles"]:
        text.extend([f"## {profile['id']} ({profile['role']})", "", "| Scope | Method | Cases | Site Top1 ≤4Å | Site Top5 ≤4Å | Pose Top1 ≤2Å, one box | Pose Top5 ≤2Å, five boxes |", "|---|---|---:|---:|---:|---:|---:|"])
        for scope in ("full_ligand", "modeled_core"):
            for method in METHODS:
                rows = [r for r in profile["conditions"] if r["atomic_scope"] == scope and r["method"] == method]
                if not rows:
                    continue
                a = mean(r["site"]["1"]["success_4A"] for r in rows)
                b = mean(r["site"]["5"]["success_4A"] for r in rows)
                c = mean(mean(s["site_budgets"]["1"]["global"]["top1_success_2A"] for s in r["seeds"]) for r in rows)
                d = mean(mean(s["site_budgets"]["5"]["global"]["top5_success_2A"] for s in r["seeds"]) for r in rows)
                text.append(f"| {scope} | {method} | {len(rows)} | {a:.3f} | {b:.3f} | {c:.3f} | {d:.3f} |")
        text.extend(["", "The JSON retains every case, all global/local fits, score ranks, best-returned endpoints, 1/2/3Å pose thresholds, 3/4/5Å site thresholds, copy and box coverage, and conditional paired family bootstrap intervals.", ""])
    text.extend(report["limitations"])
    return "\n".join(text) + "\n"


def write_report(root, plan_path, freeze_path, summary_output, markdown_output, *, check=False, allow_pending=False):
    report = generate_report(root, plan_path, freeze_path, allow_pending=allow_pending)
    summary = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    markdown = render_markdown(report)
    for path, content in ((Path(root) / summary_output, summary), (Path(root) / markdown_output, markdown)):
        if check:
            if path.read_text() != content:
                raise ValueError(f"Structural report replay changed: {path}")
        else:
            write_text_atomic(path, content)
    return report


def example_plan():
    """Coordinator fills/pins these prospective paths before any test accuracy."""
    conditions = [{"case_id": cid, "method": method,
                   "preparation": f"results/structural_pose_ideal/{cid}/{method}/prepared/condition.json",
                   "execution_root": f"results/structural_pose_ideal/{cid}/{method}/runs",
                   **({"prediction_selection": f"results/structural_test_v2/{cid}/{method}/predictor/candidate_selection.json"} if method != "oracle" else {})}
                  for cid in CASE_IDS for method in METHODS]
    pin = lambda path: {"path": path, "sha256": "REPLACE_WITH_SOURCE_ONLY_SHA256"}
    return {"schema_version": 1, "case_ids": list(CASE_IDS), "case_to_family": FAMILIES,
        "source_ledger": pin("benchmarks/heldout_reference_ledger.json"), "family_audit": pin("benchmarks/family_overlap_audit.json"),
        "reference_source_manifest": pin("data/archive/heldout_reference/source_manifest.json"),
        "site_prediction": {"configuration": pin("benchmarks/structural_site_v2_experiment.json"),
                            "freeze_receipt": pin("results/structural_test_v2/freeze_receipt.json"),
                            "source_archive_manifest": pin("results/structural_test_v2/source_archive_manifest.json")},
        "statistics": {"bootstrap_replicates": 10000, "bootstrap_seed": 20261103, "confidence": 0.95, "interpretation": "conditional_purposive_sample"},
        "profiles": [{"id": "strict_CCD_ideal", "role": "strict_primary", "configuration": pin("benchmarks/structural_pose_ideal_experiment.json"),
                      "freeze_receipt": pin("results/structural_pose_ideal/freeze_receipt.json"),
                      "source_archive_manifest": pin("results/structural_pose_ideal/source_archive_manifest.json"),
                      "source_preparation_receipt": pin("results/structural_pose_ideal/source_preparation.json"), "conditions": conditions},
                     {"id": "modeled_sidechain_CCD_ideal", "role": "modeled_sidechain_sensitivity", "configuration": pin("benchmarks/structural_pose_ideal_repair_experiment.json"),
                      "freeze_receipt": pin("results/structural_pose_ideal_repair/freeze_receipt.json"),
                      "source_archive_manifest": pin("results/structural_pose_ideal_repair/source_archive_manifest.json"),
                      "source_preparation_receipt": pin("results/structural_pose_ideal_repair/source_preparation.json"),
                      "conditions": [{**c, "preparation": c["preparation"].replace("structural_pose_ideal/", "structural_pose_ideal_repair/"),
                                      "execution_root": c["execution_root"].replace("structural_pose_ideal/", "structural_pose_ideal_repair/"),
                                      **({"primary_candidate_selection": pin(c["prediction_selection"])} if c["method"] != "oracle" else {})} for c in conditions]}],
        "freeze_receipt_schema": {"evaluation_plan_sha256": "SHA256_OF_COMPLETED_PLAN", "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False}}
