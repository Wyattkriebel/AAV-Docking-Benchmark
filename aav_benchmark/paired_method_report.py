"""Prospective paired modern-Vina/Vina-Carb report adapter.

This extension reads native search artifacts; it never starts search or scoring.
Both engines must share exact prepared inputs, site box and seed identities.
Native method scores rank poses only within their own method. Experimental
coordinates enter only after all fixed eligible outcomes are terminal.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
import tempfile

import numpy as np

from .additional_method import export_vina_carb_ligand, expected_chi_links, read_carb_poses, verify_chi_recognition
from .engine import read_poses
from .heldout_pose_report import (CASE_IDS, METHODS, _source_context, archive_files, canonical_mapping,
    protein_fits, score_pool, site_summary, site_statistics)
from .provenance import sha256, write_text_atomic
from . import public_source_replay
from .prospective_vina_carb import (PARAMETERS, PROFILES, SEEDS, normalize_atom_map,
    render_configuration, validate_frozen_plan)
from .qualified_pose import Evidence, accuracy_summary, _pose_log
from .statistics import paired_family_bootstrap


SEARCH_METHODS = ("modern_vina_1_2_7", "vina_carb_1_0")
TERMINAL_VC = {"ok", "failed", "prepared_input_failure", "unsupported_preflight"}
CODE_SOURCES = ("aav_benchmark/paired_method_report.py", "aav_benchmark/heldout_pose_report.py",
    "aav_benchmark/qualified_pose.py", "aav_benchmark/additional_method.py", "aav_benchmark/metrics.py",
    "aav_benchmark/statistics.py", "aav_benchmark/provenance.py", "aav_benchmark/public_source_replay.py",
    "aav_benchmark/heldout_reference.py")


def _pin(evidence, record):
    return evidence.read(record["path"], record["sha256"])


def check_public_reference_sources(evidence, ledger_record, source_manifest_record, distribution_record):
    """Use exactly the public replay's11-paper waiver; every numeric pin stays strict.

    Optional present papers are verified by the unchanged public replay, but
    are not numerical dependencies of this report. Their availability is not
    serialized, so public and full local trees produce identical metric inputs.
    """
    if (ledger_record["path"] != "benchmarks/heldout_reference_ledger.json"
            or ledger_record["sha256"] != public_source_replay.LEDGER_SHA256
            or source_manifest_record["path"] != str(public_source_replay.SOURCE_MANIFEST)
            or source_manifest_record["sha256"] != public_source_replay.SOURCE_MANIFEST_SHA256
            or distribution_record["path"] != str(public_source_replay.PUBLIC_MANIFEST)):
        raise ValueError("Public reference/distribution identities differ from the frozen source replay")
    ledger = json.loads(_pin(evidence, ledger_record))
    original = json.loads(_pin(evidence, source_manifest_record))
    distribution = json.loads(_pin(evidence, distribution_record))
    replay = public_source_replay.verify_public_source_references(evidence.root)
    if (replay.get("all_consumed_numerical_sources_verified") is not True
            or replay.get("retained_archive_records_verified") != 163
            or replay.get("metadata_only_publication_count") != 11):
        raise ValueError("Public numerical replay did not verify every retained source")
    omitted = {record["archive"] for record in distribution["metadata_only_publication_snapshots"]}
    for entry in original["archived_files"]:
        if entry["archive"] in omitted:
            continue  # Already scoped and strictly checked when present by public replay.
        raw = evidence.read(entry["archive"], entry["sha256"])
        if len(raw) != entry["bytes"]:
            raise ValueError("Retained public numerical archive size changed")
    for key, value in ledger.items():
        if key.endswith("_archive"):
            hash_key = key.removesuffix("_archive") + "_sha256"
            if hash_key in ledger:
                evidence.read(value, ledger[hash_key])
    by_id = {case["case_id"]: case for case in ledger["cases"]}
    cases = {cid: by_id[cid] for cid in CASE_IDS}
    if (len(by_id) != 12 or ledger["cohort_denominator"] != 12
            or sum(c["qualification"] == "full_ligand" for c in cases.values()) != 5
            or cases["4R9F"]["qualification"] != "modeled_core" or cases["4R9F"]["scored_atom_count"] != 11):
        raise ValueError("Frozen five-full/one-core structural reference denominator changed")
    summary = {"distribution_id": replay["distribution_id"], "public_manifest_sha256": distribution_record["sha256"],
               "retained_archive_records_verified": 163, "metadata_only_publication_count": 11,
               "original_source_manifest_records": 174, "all_consumed_numerical_sources_verified": True,
               "original_source_coordinate_frames_and_masks_unchanged": True,
               "publication_judgments_independently_reproven": False}
    return cases, summary


def portable_archive(evidence, record, *, archive_path=None):
    """Accept an explicit frozen archive path when the V2 manifest omits it."""
    manifest = json.loads(_pin(evidence, record))
    if "archive" in manifest:
        if archive_path is not None and Path(record["path"]).parent / manifest["archive"] != Path(archive_path):
            raise ValueError("Configuration and manifest disagree on the source archive identity")
        return archive_files(evidence, record)
    if archive_path is None:
        raise ValueError("A source manifest without archive name requires its explicit frozen configuration path")
    raw = evidence.read(archive_path, manifest["archive_sha256"])
    expected = {entry["path"]: entry for entry in manifest["files"]}
    if len(expected) != len(manifest["files"]) or any(Path(p).is_absolute() or ".." in Path(p).parts for p in expected):
        raise ValueError("Invalid source archive member identities")
    result = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        for member in archive.getmembers():
            if member.isdir():
                continue
            if not member.isfile() or member.name not in expected or member.name in result:
                raise ValueError("Unexpected source archive member")
            content = archive.extractfile(member).read()
            if hashlib.sha256(content).hexdigest() != expected[member.name]["sha256"] or len(content) != expected[member.name]["bytes"]:
                raise ValueError("Source archive member changed")
            result[member.name] = content
    if set(result) != set(expected):
        raise ValueError("Source archive is incomplete")
    return result


def _historical_result_pin(evidence, record):
    """Read a pinned original result after repository relocation, preserving it."""
    path = Path(record["path"])
    if path.is_absolute():
        indices = [i for i, part in enumerate(path.parts) if part == "results"]
        if len(indices) != 1:
            raise ValueError("Historical selection path has no unambiguous portable results suffix")
        path = Path(*path.parts[indices[0]:])
    return evidence.read(path, record["sha256"]), path


def common_selection(evidence, row, prediction_config_sha, prediction_sources):
    """Verify full V2 selection while executing only its first class candidate."""
    condition = json.loads(_pin(evidence, row["source_condition"]))
    if (condition.get("case_id", condition.get("complex_id")) != row["case_id"] or condition["method"] != row["method"]
            or condition.get("profile") != row["receptor_profile"] or condition.get("fixed_denominator_member") is not True
            or condition.get("executed_primary_candidate_budget") != 1):
        raise ValueError("Paired source condition changed receptor/site identity or one-box denominator")
    if row["method"] == "oracle":
        return condition, [{"rank": j["candidate_class_rank"], "original_rank": j["candidate_original_rank"], "center_A": j["center_A"]} for j in condition["jobs"]]
    selection = condition["candidate_selection"]
    raw, path = _historical_result_pin(evidence, condition["primary_candidate_selection"])
    if json.loads(raw) != selection:
        raise ValueError("Prepared condition changed the full original V2 candidate selection")
    if (selection["case_id"] != row["case_id"] or selection["method"] != row["method"]
            or selection.get("native_coordinate_access") is not False or selection["configuration_sha256"] != prediction_config_sha
            or selection["prediction_input_sha256"] != hashlib.sha256(prediction_sources[f"{row['case_id']}/receptor.cif"]).hexdigest()):
        raise ValueError("Original V2 prediction input/configuration changed")
    record = selection["prediction_record"]
    predictor = evidence.json(path.parent / record["filename"], record["sha256"])
    if predictor["status"] != selection["predictor_status"]:
        raise ValueError("Original prediction status changed after source selection")
    candidates = selection["candidates"]
    if ([c["rank"] for c in candidates] != list(range(1, len(candidates) + 1)) or len(candidates) > 5
            or condition.get("available_prediction_candidates") != len(candidates)
            or condition.get("scheduled_candidate_jobs") != min(1, len(candidates))
            or len(condition["jobs"]) != min(1, len(candidates))):
        raise ValueError("Full site ranking or one-box execution budget changed")
    if candidates:
        job, candidate = condition["jobs"][0], candidates[0]
        if job["candidate_class_rank"] != 1 or job["candidate_original_rank"] != candidate["original_rank"] or job["center_A"] != candidate["center_A"]:
            raise ValueError("The executed box is not the original top1 candidate class")
    return condition, candidates


def canonical_carb_poses(path, export, names):
    """Restore Lxxx identities and reference ordering by invariant serial map."""
    names = list(names)
    if len(names) != len(set(names)) or not names:
        raise ValueError("Canonical parent names must be complete and unique")
    result = []
    for pose in read_carb_poses(path, export):
        actual = pose["names"]
        if len(actual) != len(set(actual)) or set(actual) != set(names):
            raise ValueError("Vina-Carb pose lost or duplicated a complete parent atom")
        xyz = np.asarray(pose["coordinates"], float)
        result.append({"rank": pose["rank"], "method_score": pose["vina_carb_score"],
                       "coordinates": xyz[[actual.index(name) for name in names]]})
    return result


def verify_carb_score_table(content, poses):
    """Check the original six-column VC table, retaining native source ranks.

    The original main.cpp prints affinity, CHI and affinity-minus-CHI before
    the two RMSD columns. The modern Vina four-column table is a different
    format. This validates raw method scores, never experimental accuracy.
    """
    text = content.decode()
    if not re.search(r"mode\s*\|\s*affinity\s*\|\s*chi\s*\|\s*affinity\s*\|", text, re.I):
        raise ValueError("Original Vina-Carb six-column score-table header is missing")
    number = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    pattern = re.compile(r"^\s*(\d+)\s+" + r"\s+".join([number] * 5) + r"\s*$")
    rows = []
    for line in text.splitlines():
        match = pattern.match(line)
        if match:
            rank, *values = match.groups()
            values = [float(v) for v in values]
            if not all(math.isfinite(v) for v in values) or values[-2] < 0 or values[-1] < 0:
                raise ValueError("Vina-Carb raw score table contains invalid numeric fields")
            rows.append((int(rank), values[0]))
    if len(rows) != len(poses) or any(rank != pose["rank"] or abs(score - pose["method_score"]) > .002
                                       for (rank, score), pose in zip(rows, poses)):
        raise ValueError("Original Vina-Carb raw score table disagrees with retained poses")
    return {"matched_models": len(rows), "table_columns": 6, "raw_score_tolerance_kcal_mol": .002}


def _artifacts(evidence, directory, record):
    for name, digest in record.get("artifacts", {}).items():
        if Path(name).name != name or name in {".", "..", "execution.json"}:
            raise ValueError("Unsafe Vina-Carb artifact identity")
        evidence.read(directory / name, digest)


def read_vina_carb(evidence, directory, row, seed, experiment, experiment_pin, receipt_pin,
                   names, *, executable_sha256):
    """Validate one terminal native VC record or return its pending identity."""
    path = directory / "execution.json"
    if not (evidence.root / path).is_file():
        return None, "missing_execution"
    record = evidence.json(path)
    if (record.get("condition") != row or record.get("seed") != seed
            or record.get("configuration_sha256") != experiment_pin["sha256"]
            or record.get("freeze_receipt_sha256") != receipt_pin["sha256"]
            or record.get("driver_source_sha256") != experiment["adapter_source"]["sha256"]
            or record.get("atomic_accuracy_computed") is not False
            or record.get("experimental_ligand_coordinates_read") is not False):
        raise ValueError("Vina-Carb execution changed its frozen condition/input/seed identity")
    if record["status"] in {"running", "preparing"}:
        return None, record["status"]
    if record["status"] not in TERMINAL_VC:
        raise ValueError("Unknown Vina-Carb terminal state")
    _artifacts(evidence, directory, record)
    if not isinstance(record.get("full_search_started"), bool):
        raise ValueError("Vina-Carb record must distinguish unexecuted and executed searches")
    if record["full_search_started"]:
        if record.get("parameters") != PARAMETERS or record.get("executable_sha256") != executable_sha256:
            raise ValueError("Vina-Carb attempt changed frozen effort/executable")
        if evidence.read(directory / "configuration.txt") != render_configuration(row, seed).encode():
            raise ValueError("Actual Vina-Carb configuration differs from the frozen paired box/seed")
        command = record.get("command", [])
        if len(command) != 3 or command[1:] != ["--config", "configuration.txt"]:
            raise ValueError("Actual Vina-Carb invocation changed the native full-search mode")
        evidence.read(directory / "receptor.pdbqt", row["receptor"]["sha256"])
        if record.get("receptor_sha256") != row["receptor"]["sha256"]:
            raise ValueError("Vina-Carb attempt used a different paired receptor")
    poses = []
    if row["status"] != "prepared":
        if record["status"] != "prepared_input_failure" or record["full_search_started"] is not False:
            raise ValueError("Failed common preparation must remain an explicit unexecuted VC failure")
    elif record["status"] == "prepared_input_failure":
        raise ValueError("A prepared common input was replaced by an unsupported-source claim")
    if record["status"] == "ok":
        needed = {"poses.pdbqt", "ligand.pdbqt", "ligand.pdbqt.map.json", "receptor.pdbqt", "configuration.txt",
                  "VC_log.txt", "search.log", "stdout.log", "stderr.log"}
        if not needed <= set(record.get("artifacts", {})):
            raise ValueError("Successful Vina-Carb execution omitted its raw evidence")
        if (record.get("full_search_started") is not True or record.get("exit_code") != 0
                or record.get("parameters") != PARAMETERS or record.get("executable_sha256") != executable_sha256):
            raise ValueError("Successful Vina-Carb search changed frozen effort/executable")
        export = evidence.json(directory / "ligand.pdbqt.map.json")
        exported = evidence.read(directory / "ligand.pdbqt", record["exported_ligand_sha256"])
        source_ligand = _pin(evidence, row["ligand"])
        atom_map = normalize_atom_map(json.loads(_pin(evidence, row["canonical_atom_map"])))
        with tempfile.TemporaryDirectory() as tmp:
            source, destination = Path(tmp) / "source.pdbqt", Path(tmp) / "export.pdbqt"
            source.write_bytes(source_ligand)
            replay = export_vina_carb_ligand(source, atom_map, destination)
            if export != replay or exported != destination.read_bytes():
                raise ValueError("Vina-Carb alias export changed coordinates, charges, types, torsions or canonical identity")
            pose_path = Path(tmp) / "poses.pdbqt"
            pose_path.write_bytes(evidence.read(directory / "poses.pdbqt"))
            poses = canonical_carb_poses(pose_path, export, names)
            log_path = Path(tmp) / "VC_log.txt"
            log_path.write_bytes(evidence.read(directory / "VC_log.txt"))
            recognition = verify_chi_recognition(log_path, expected_chi_links(export))
        if (not export.get("coordinates_charges_types_and_torsion_tree_preserved")
                or export.get("native_coordinates_accessed") is not False
                or export.get("source_sha256") != row["ligand"]["sha256"]
                or len(poses) != record.get("returned_poses") or not 1 <= len(poses) <= PARAMETERS["num_modes"]
                or recognition["recognized_links"] != len(export["expected_glycosidic_links"])):
            raise ValueError("Vina-Carb full-parent/linkage/returned-pose denominator changed")
        preflight = record.get("preflight", {})
        if (preflight.get("glycosidic_links") != len(export["expected_glycosidic_links"])
                or preflight.get("fully_supported_links") != len(export["expected_glycosidic_links"])):
            raise ValueError("Vina-Carb successful search lacks complete frozen-input CHI preflight support")
        if len(record.get("generated_names", [])) != len(names) or set(record.get("generated_names", [])) != set(names):
            raise ValueError("Vina-Carb execution did not restore the canonical parent names")
        verify_carb_score_table(evidence.read(directory / "search.log"), poses)
    for field in ("cpu_seconds", "wall_seconds"):
        value = record.get(field)
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
            raise ValueError("Terminal Vina-Carb record lacks measured or explicit unexecuted resource use")
    return {"status": record["status"], "poses": poses, "record": record}, None


def _modern_command(run, prep, config, prepared, output):
    command = run.get("command", [])
    if not command or len(command) % 2 != 1 or any(not isinstance(v, str) for v in command):
        raise ValueError("Missing actual modern-Vina invocation")
    options = dict(zip(command[1::2], command[2::2]))
    if len(options) * 2 != len(command) - 1:
        raise ValueError("Duplicate modern-Vina command flags")
    expected = {"--scoring": "vina", "--seed": run["seed"], "--cpu": config["cpu"], "--exhaustiveness": config["exhaustiveness"],
                "--num_modes": config["num_modes"], "--energy_range": config["energy_range_kcal_mol"], "--min_rmsd": config["min_rmsd_A"]}
    for axis, center, size in zip("xyz", prep["box_center"], prep["box_size"]):
        expected["--center_" + axis], expected["--size_" + axis] = center, size
    for flag, value in expected.items():
        actual = options.get(flag)
        if actual is None or (isinstance(value, str) and actual != value) or (not isinstance(value, str) and (not math.isfinite(float(actual)) or float(actual) != float(value))):
            raise ValueError("Modern-Vina command changed a frozen paired search setting")
    suffixes = {"--ligand": prepared / prep["ligand_file"], "--receptor": prepared / "neighborhood.pdbqt", "--out": output / run["pose_file"]}
    if set(options) != set(expected) | set(suffixes):
        raise ValueError("Modern-Vina command changed the search mode")
    for flag, path in suffixes.items():
        if Path(options[flag]).parts[-len(path.parts):] != path.parts:
            raise ValueError("Modern-Vina command used a different paired input/output")


def read_modern(evidence, row, config, config_pin, names, *, input_manifest_sha256):
    """Read all three modern-Vina rows once; failed sources are explicit zeros."""
    condition = json.loads(_pin(evidence, row["source_condition"]))
    if condition.get("case_id", condition.get("complex_id")) != row["case_id"] or condition["method"] != row["method"]:
        raise ValueError("Common preparation failure/source condition changed identity")
    receipt_path = Path(row["modern_runs"]).parent.parent / "execution.json"
    if not (evidence.root / receipt_path).is_file():
        return {}, ["missing_condition_execution"]
    receipt = evidence.json(receipt_path)
    if (receipt.get("condition_id") != row["condition_id"] or receipt.get("configuration_sha256") != config_pin["sha256"]
            or receipt.get("frozen_input_manifest_sha256") != input_manifest_sha256 or receipt.get("source_condition") != row["source_condition"]):
        raise ValueError("Modern-Vina condition receipt changed its paired input/configuration identity")
    if receipt["status"] in {"preparing", "running", "pending"}:
        return {}, ["condition_execution_" + receipt["status"]]
    if row["status"] != "prepared":
        if any(job["status"] == "prepared" for job in condition["jobs"]):
            raise ValueError("Failed common input hides a prepared modern-Vina box")
        if (receipt.get("status") != "prepared_input_failure" or receipt.get("native_search_started") is not False
                or receipt.get("eligible_seed_slots") != SEEDS or receipt.get("no_runs_or_poses_invented") is not True):
            raise ValueError("Failed modern-Vina preparation lacks its explicit terminal three-seed receipt")
        return {seed: {"status": "prepared_input_failure", "poses": [], "record": receipt} for seed in SEEDS}, []
    prepared, manifest = Path(row["modern_vina_prepared_directory"]), Path(row["modern_runs"])
    prep_raw = _pin(evidence, row["preparation"])
    prep, reference = json.loads(prep_raw), evidence.json(prepared / "reference.json")
    if (reference.get("pose_reference_eligible") is not False or "coordinates" in reference
            or reference["names"] != names or prep["box_center"] != row["box"]["center"] or prep["box_size"] != row["box"]["size"]
            or prep["configuration_sha256"] != config_pin["sha256"]):
        raise ValueError("Modern-Vina preparation changed paired geometry or leaked atomic truth")
    jobs = [job for job in condition["jobs"] if job["candidate_class_rank"] == row["candidate_class_rank"]]
    if (len(jobs) != 1 or jobs[0]["status"] != "prepared" or jobs[0]["preparation_sha256"] != row["preparation"]["sha256"]
            or jobs[0]["center_A"] != row["box"]["center"]):
        raise ValueError("Paired prepared box differs from the frozen selected candidate")
    if not (evidence.root / manifest).is_file():
        return {}, ["missing_runs"]
    if (receipt.get("status") != "complete" or receipt.get("native_search_started") is not True
            or receipt.get("runs_path") != row["modern_runs"] or receipt.get("terminal_seed_runs") != len(SEEDS)):
        raise ValueError("Prepared modern-Vina condition lacks its complete three-seed execution receipt")
    evidence.read(manifest, receipt["runs_sha256"])
    evidence.read(manifest.parent / "configuration.json", config_pin["sha256"])
    runs = evidence.json(manifest)
    identities = [r["seed"] for r in runs]
    if len(identities) != len(set(identities)) or set(identities) - set(SEEDS):
        raise ValueError("Unknown or duplicate paired modern-Vina seed")
    if set(identities) != set(SEEDS) or any(r["status"] not in {"ok", "failed"} for r in runs):
        raise ValueError("Complete modern-Vina condition receipt lacks three actual terminal seed rows")
    if receipt.get("seed_statuses") != {str(r["seed"]): r["status"] for r in runs}:
        raise ValueError("Modern-Vina condition receipt changed its actual terminal seed outcomes")
    pending = [f"missing_seed:{seed}" for seed in SEEDS if seed not in identities]
    vina = evidence.json("benchmarks/vina.lock.json")
    executables = {v["sha256"] for v in vina["assets"].values()}
    fingerprints = {"config_sha256": config_pin["sha256"], "preparation_sha256": row["preparation"]["sha256"],
        "ligand_sha256": row["ligand"]["sha256"], "receptor_sha256": row["receptor"]["sha256"],
        "receptor_pdb_sha256": sha256(evidence.root / prepared / "neighborhood.pdb"),
        "reference_sha256": sha256(evidence.root / prepared / "reference.json"), "full_assembly_sha256": None}
    evidence.read(prepared / "neighborhood.pdb", fingerprints["receptor_pdb_sha256"])
    result = {}
    for run in runs:
        if (any(run["fingerprint"].get(k) != value for k, value in fingerprints.items())
                or run["fingerprint"].get("executable_sha256") not in executables or run.get("case_id") != row["case_id"]
                or run.get("receptor") != "neighborhood" or run.get("vina_version") != "AutoDock Vina v1.2.7"
                or run.get("information_regime") != ("true-site-local" if row["method"] == "oracle" else "predicted-site-local")):
            raise ValueError("Modern-Vina outcome does not use the exact paired input fingerprints")
        _modern_command(run, prep, config, prepared, manifest.parent)
        if run["status"] in {"running", "pending"}:
            pending.append(f"running_seed:{run['seed']}")
            continue
        if run["status"] not in {"ok", "failed"}:
            raise ValueError("Unknown terminal modern-Vina status")
        for name in (run["pose_file"], run["log_file"]):
            if Path(name).name != name or name in {".", ".."}:
                raise ValueError("Unsafe modern-Vina raw artifact")
        log = evidence.read(manifest.parent / run["log_file"])
        poses = []
        if run["status"] == "ok":
            if run.get("exit_code") != 0:
                raise ValueError("Successful modern-Vina outcome has nonzero exit status")
            raw = evidence.read(manifest.parent / run["pose_file"], run["pose_sha256"])
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "poses.pdbqt"; path.write_bytes(raw)
                parsed = read_poses(path, names)
            if len(parsed) != run["returned_poses"] or not 1 <= len(parsed) <= config["num_modes"]:
                raise ValueError("Modern-Vina complete-parent pose denominator changed")
            _pose_log(log, parsed)
            poses = [{"rank": p["rank"], "method_score": p["vina_score"], "coordinates": p["coordinates"]} for p in parsed]
        for field in ("cpu_seconds", "wall_seconds"):
            value = run.get(field)
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
                raise ValueError("Terminal modern-Vina record lacks measured resource use")
        result[run["seed"]] = {"status": run["status"], "poses": poses, "record": run}
    return result, pending


def read_execution_pair(evidence, row, modern_config, experiment, experiment_pin, receipt_pin, *, executable_sha256):
    """Pure computational-evidence verifier for one three-seed paired condition.

    No reference ledger, experimental ligand coordinates, fit, RMSD or native
    accuracy helper is read/called. Common-pool builders may reuse this reader.
    The caller separately validates the full experiment/freeze/36-row manifest.
    """
    names = ([a["name"] for a in normalize_atom_map(json.loads(_pin(evidence, row["canonical_atom_map"])))]
             if row["status"] == "prepared" else [])
    modern, pending = read_modern(evidence, row, modern_config, experiment["paired_modern_vina_policy"], names,
        input_manifest_sha256=experiment["prepared_inputs"]["sha256"])
    pending = ["modern/" + p for p in pending]
    if set(row["vina_carb_outputs"]) != {str(s) for s in SEEDS}:
        raise ValueError("Explicit Vina-Carb seed-output paths changed the frozen denominator")
    seeds = []
    for seed in SEEDS:
        carb, unfinished = read_vina_carb(evidence, Path(row["vina_carb_outputs"][str(seed)]), row, seed, experiment,
            experiment_pin, receipt_pin, names, executable_sha256=executable_sha256)
        if unfinished:
            pending.append(f"vina_carb/{seed}:{unfinished}")
        seeds.append({"seed": seed, "search_methods": {SEARCH_METHODS[0]: modern.get(seed), SEARCH_METHODS[1]: carb}})
    return {"condition": row, "canonical_names": names, "seeds": seeds, "pending": pending,
            "atomic_accuracy_computed": False, "experimental_reference_access": False}


def paired_statistics(conditions, case_to_family):
    results = []
    for profile in PROFILES:
        for site in METHODS:
            for scope in ("full_ligand", "modeled_core"):
                selected = [r for r in conditions if r["receptor_profile"] == profile and r["site_condition"] == site and r["atomic_scope"] == scope]
                for alignment in ("global", "local_interface"):
                    for rank in ("top1", "top5", "best_returned"):
                        for threshold in (1, 2, 3):
                            records = [{"family_id": case_to_family[r["case_id"]], "complex_id": r["case_id"], "method": method,
                                "seed": seed["seed"], "status": "ok" if value["poses"] else "unsupported" if value.get("status") in {"prepared_input_failure", "unsupported_preflight"} else "failed",
                                "success": value[alignment][f"{rank}_success_{threshold}A"]}
                                for r in selected for seed in r["seeds"] for method, value in seed["search_methods"].items()]
                            if not records:
                                continue
                            summary = paired_family_bootstrap(records, SEARCH_METHODS[1], SEARCH_METHODS[0], n_resamples=10000, seed=20261103, expected_seeds=SEEDS)
                            summary.update(receptor_profile=profile, site_condition=site, atomic_scope=scope, alignment=alignment,
                                pose_rank=rank, threshold_A=threshold, executed_candidate_budget=1,
                                interpretation="Conditional source-family-cluster resampling of a purposive structural sample; VC-minus-modern recovery, never a score difference; no model-unseen/population claim")
                            results.append(summary)
    return results


def generate_paired_report(root, plan_path, freeze_path, *, allow_pending=False):
    evidence = Evidence(root)
    plan, freeze = evidence.json(plan_path), evidence.json(freeze_path)
    if (freeze.get("evaluation_plan_sha256") != sha256(evidence.root / plan_path)
            or freeze.get("frozen_before_test_pose_accuracy") is not True or freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Paired method evaluation needs its explicit pre-accuracy freeze")
    if (plan["case_ids"] != list(CASE_IDS) or plan["seeds"] != SEEDS or plan["executed_candidate_budget"] != 1
            or plan["search_methods"] != list(SEARCH_METHODS) or plan.get("common_pool_rescoring_enabled") is not False):
        raise ValueError("Frozen paired cohort/rank design changed; common-pool scoring needs a separate prospective extension")
    if set(plan["metric_code_sha256"]) != set(CODE_SOURCES):
        raise ValueError("Every shared metric/parser implementation needs its frozen source pin")
    for path, digest in plan["metric_code_sha256"].items():
        evidence.read(path, digest)
    reference_plan = json.loads(_pin(evidence, plan["reference_evaluation_plan"]))
    reference_freeze = json.loads(_pin(evidence, plan["reference_evaluation_freeze"]))
    if (reference_freeze.get("evaluation_plan_sha256") != plan["reference_evaluation_plan"]["sha256"]
            or reference_freeze.get("frozen_before_test_pose_accuracy") is not True or reference_freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Shared reference/metric plan changed its pre-accuracy freeze")
    cases, reference_validation = check_public_reference_sources(evidence, reference_plan["source_ledger"],
        reference_plan["reference_source_manifest"], plan["public_source_distribution"])
    _pin(evidence, reference_plan["family_audit"])
    prediction = reference_plan["site_prediction"]
    prediction_config = json.loads(_pin(evidence, prediction["configuration"]))
    prediction_freeze = json.loads(_pin(evidence, prediction["freeze_receipt"]))
    if (prediction_freeze.get("configuration_sha256") != prediction["configuration"]["sha256"]
            or prediction_freeze.get("frozen_before_new_predictions") is not True):
        raise ValueError("Shared site input/ranking plan lacks its pre-prediction freeze")
    prediction_archive = prediction_config["source_preparation_archive"]
    if (prediction_archive["manifest"] != prediction["source_archive_manifest"]["path"]
            or prediction_archive["manifest_sha256"] != prediction["source_archive_manifest"]["sha256"]):
        raise ValueError("Shared prediction configuration changed its portable source-manifest pin")
    prediction_sources = portable_archive(evidence, prediction["source_archive_manifest"], archive_path=prediction_archive["path"])
    experiment = json.loads(_pin(evidence, plan["vina_carb_experiment"]))
    _pin(evidence, plan["vina_carb_freeze"])
    # Pure validation only: this function cannot invoke search/native scoring.
    validated, rows = validate_frozen_plan(evidence.root, evidence.root / plan["vina_carb_experiment"]["path"], evidence.root / plan["vina_carb_freeze"]["path"])
    if validated != experiment or experiment["secondary_top5_enabled"] is not False:
        raise ValueError("This paired report requires the reviewed top1-box-only native search plan")
    if plan["vina_carb_executable_sha256"] != experiment["native_binary_sha256"]:
        raise ValueError("Evaluation and native-execution freezes disagree on the Vina-Carb binary")
    for key in ("method_lock", "initialization_policy", "initialization_algorithm", "prepared_inputs", "adapter_source", "paired_modern_vina_policy", "native_build_receipt", "native_preflight_receipt"):
        _pin(evidence, experiment[key])
    for path, digest in experiment["helper_sources"].items():
        evidence.read(path, digest)
    modern_config = json.loads(_pin(evidence, experiment["paired_modern_vina_policy"]))
    if modern_config["seeds"] != SEEDS or any(modern_config[k] != PARAMETERS[k] for k in ("cpu", "exhaustiveness", "num_modes", "energy_range_kcal_mol")):
        raise ValueError("Modern-Vina nominal per-box search effort differs from the paired native plan")
    if len(rows) != 36 or len({(r["case_id"], r["receptor_profile"], r["method"]) for r in rows}) != 36:
        raise ValueError("All six cases, both receptor profiles and all three site conditions remain eligible")
    by_id = {p["id"]: p for p in reference_plan["profiles"]}
    sources = {}
    for profile in PROFILES:
        profile_id = plan["source_profiles"][profile]
        if profile_id not in by_id:
            raise ValueError("Paired source profile is absent from the frozen reference plan")
        sources[profile] = archive_files(evidence, by_id[profile_id]["source_archive_manifest"])
    for cid in CASE_IDS:
        for name in ("ligand_atom_map.json", "ligand.sdf", "ligand.pdbqt", "protein_atom_map.json"):
            if sources["strict"][f"{cid}/{name}"] != sources["repair"][f"{cid}/{name}"]:
                raise ValueError("Repair sensitivity changed the fixed ligand or observed reference protein map")
    loaded, pending = [], []
    counts = {"eligible_conditions": 36, "eligible_seed_pairs": 108, "eligible_method_outcomes": 216,
              "terminal_method_outcomes": 0, "prepared_input_failed_conditions": 0, "successful_modern_runs": 0, "successful_vina_carb_runs": 0}
    for row in rows:
        cid, profile = row["case_id"], row["receptor_profile"]
        atoms = json.loads(sources[profile][f"{cid}/ligand_atom_map.json"])["atoms"]
        names = [a["name"] for a in atoms]
        canonical_mapping(cases[cid], {"names": names, "canonical_atom_map": atoms})
        # Source preparation and all paired ordinary input pins are inventoried.
        condition, candidates = common_selection(evidence, row, prediction["configuration"]["sha256"], prediction_sources)
        if row["status"] == "prepared":
            for key in ("ligand", "canonical_atom_map", "initializer_receipt", "initializer_sdf", "receptor", "preparation"):
                _pin(evidence, row[key])
            if _pin(evidence, row["initializer_sdf"]) != sources[profile][f"{cid}/ligand.sdf"]:
                raise ValueError("Paired initializer SDF differs from the pinned metric parent graph")
            if normalize_atom_map(json.loads(_pin(evidence, row["canonical_atom_map"]))) != normalize_atom_map({"atoms": atoms}):
                raise ValueError("Paired input changed canonical chemical atom identities")
        else:
            counts["prepared_input_failed_conditions"] += 1
        pair = read_execution_pair(evidence, row, modern_config, experiment, plan["vina_carb_experiment"], plan["vina_carb_freeze"],
                                   executable_sha256=plan["vina_carb_executable_sha256"])
        pending.extend(f"{row['condition_id']}/{p}" for p in pair["pending"])
        for seed in pair["seeds"]:
            for method, result in seed["search_methods"].items():
                if result is not None:
                    counts["terminal_method_outcomes"] += 1
                    if result["status"] == "ok":
                        counts["successful_modern_runs" if method == SEARCH_METHODS[0] else "successful_vina_carb_runs"] += 1
        loaded.append({"row": row, "seeds": pair["seeds"], "candidates": candidates, "condition": condition})
    if pending:
        if not allow_pending:
            raise ValueError(f"Paired native execution remains incomplete ({len(pending)} pending records); no accuracy report")
        return {"schema_version": 1, "status": "pending_incomplete_workflow", "metrics_eligible": False,
                "counts": counts, "pending": sorted(pending), "source_sha256_inventory": evidence.sources,
                "reference_validation": reference_validation}
    if counts["terminal_method_outcomes"] != counts["eligible_method_outcomes"]:
        raise ValueError("The complete fixed216-outcome paired denominator is not terminal")
    contexts = {(profile, cid): _source_context(sources[profile], cid, cases[cid]) for profile in PROFILES for cid in CASE_IDS}
    evaluated = []
    for item in loaded:
        row = item["row"]; cid, profile = row["case_id"], row["receptor_profile"]
        mapping, graph, keys, protein, protein_map = contexts[profile, cid]
        fit = protein_fits(keys, protein, protein_map, mapping["copies"][0]["coordinates"])["fits"]["global"]
        chemistry = json.loads(sources[profile][f"{cid}/chemistry.json"])
        site = site_summary(item["candidates"], mapping["copies"], box_sizes=chemistry["box_size"], rotation=fit["rotation"], translation=fit["translation"])
        results, controls = [], None
        for seed in item["seeds"]:
            methods = {}
            for method, result in seed["search_methods"].items():
                raw = [{"rank": p["rank"], "vina_score": p["method_score"], "coordinates": p["coordinates"],
                        "candidate_class_rank": row.get("candidate_class_rank", 1), "candidate_original_rank": row.get("candidate_original_rank", 1),
                        "box_id": "box01", "pooled_rank": p["rank"]} for p in result["poses"]]
                scored, controls = score_pool(raw, mapping, keys, protein, protein_map, graph)
                for pose in scored:
                    pose["native_method_score"] = {"method": method, "value": pose.pop("vina_score"), "cross_method_comparable": False}
                methods[method] = {"status": result["status"], "global": accuracy_summary(scored, "global_RMSD_A"),
                                   "local_interface": accuracy_summary(scored, "local_interface_RMSD_A"), "poses": scored,
                                   "resource_use": {k: v for k, v in (result["record"] or {}).items() if k in
                                      {"cpu_seconds", "wall_seconds", "preflight_cpu_seconds", "preflight_wall_seconds", "total_native_cpu_seconds", "total_native_wall_seconds"}}}
            results.append({"seed": seed["seed"], "search_methods": methods})
        evaluated.append({"condition_id": row["condition_id"], "case_id": cid, "receptor_profile": profile, "site_condition": row["method"],
            "atomic_scope": cases[cid]["qualification"], "scored_atoms": cases[cid]["scored_atom_count"], "parent_atoms": mapping["parent_atom_count"],
            "common_preparation_status": row["status"], "executed_candidate_budget": 1, "native_controls": controls, "seeds": results,
            "site": site, "site_workflow_shared_across_search_methods": True,
            "exact_common_input_sha256": {k: row[k]["sha256"] for k in ("ligand", "receptor", "preparation", "canonical_atom_map", "initializer_sdf") if k in row}})
    return {"schema_version": 1, "status": "complete", "metrics_eligible": True, "counts": counts, "conditions": evaluated,
            "paired_family_statistics": paired_statistics(evaluated, reference_plan["case_to_family"]),
            "shared_site_statistics": [{"receptor_profile": profile, "paired_source_family_intervals": site_statistics(
                [{**r, "method": r["site_condition"]} for r in evaluated if r["receptor_profile"] == profile], reference_plan["case_to_family"])} for profile in PROFILES],
            "source_sha256_inventory": evidence.sources, "metric_code_sha256": plan["metric_code_sha256"],
            "reference_validation": reference_validation,
            "limitations": ["Vina-Carb1.0 uses the older Vina1.1.2 engine; this is a paired implementation comparison, not an isolated causal CHI ablation.",
                "Each engine ranks its own returned poses by its own native score. Absolute scores are never subtracted, calibrated or pooled across methods.",
                "CPU allocation/exhaustiveness/output caps/box/seed are nominally matched; wall time and algorithmic search behavior are not identical.",
                "The original VC minimum output RMSD is hardcoded; modern Vina has its separately frozen min_rmsd setting.",
                "All six cases and preparation/unsupported-search failures remain eligible; five full ligands and the11-atom4R9Fcore are separate strata.",
                "Strict and modeled-sidechain profiles remain separate; generated protein atoms are excluded from experimental reference truth.",
                "Conditional family-cluster intervals describe this purposive sample; no population or model-training-unseen reliability claim.",
                "HADDOCK ambiguous restraints encode different guidance and require a separately declared guided-method profile; no equivalent-information assertion is made."]}


def render_markdown(report):
    if not report["metrics_eligible"]:
        return "# Paired method comparison pending\n\nThe fixed paired execution denominator is incomplete. No native pose accuracy is reported.\n"
    text = ["# Paired modern Vina and Vina-Carb", "", "Both engines use the same independently initialized ligand, receptor, box and seed. All preparation failures remain eligible. Scores rank poses within each engine and are never compared across engines.", "",
            "| Receptor profile | Atomic scope | Site condition | Search method | Cases | Top1 ≤2Å | Top5 ≤2Å | Best retained ≤2Å |",
            "|---|---|---|---|---:|---:|---:|---:|"]
    from statistics import mean
    for profile in PROFILES:
        for scope in ("full_ligand", "modeled_core"):
            for site in METHODS:
                rows = [r for r in report["conditions"] if r["receptor_profile"] == profile and r["atomic_scope"] == scope and r["site_condition"] == site]
                for method in SEARCH_METHODS:
                    metrics = [mean(mean(seed["search_methods"][method]["global"][f"{rank}_success_2A"] for seed in r["seeds"]) for r in rows) for rank in ("top1", "top5", "best_returned")]
                    text.append(f"| {profile} | {scope} | {site} | {method} | {len(rows)} | " + " | ".join(f"{v:.3f}" for v in metrics) + " |")
    text.extend(["", "Pose Top5 is within the single executed candidate box. The JSON retains global/local protein fits, 1/2/3Å endpoints, native controls, full-parent geometry/clashes/contacts, per-case outcomes and conditional paired family intervals.", ""])
    text.extend(report["limitations"])
    return "\n".join(text) + "\n"


def write_report(root, plan_path, freeze_path, summary_output, markdown_output, *, check=False, allow_pending=False):
    report = generate_paired_report(root, plan_path, freeze_path, allow_pending=allow_pending)
    for path, content in ((Path(root) / summary_output, json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"),
                          (Path(root) / markdown_output, render_markdown(report))):
        if check:
            if path.read_text() != content:
                raise ValueError(f"Paired method report replay changed: {path}")
        else:
            write_text_atomic(path, content)
    return report


def example_plan(root=None):
    root = Path(root) if root is not None else None
    def pin(path):
        return {"path": path, "sha256": sha256(root / path) if root is not None and (root / path).is_file() else "FILL_SOURCE_SHA256"}
    return {"schema_version": 1, "case_ids": list(CASE_IDS), "seeds": SEEDS, "executed_candidate_budget": 1,
        "search_methods": list(SEARCH_METHODS), "common_pool_rescoring_enabled": False,
        "reference_evaluation_plan": pin("benchmarks/structural_pose_evaluation.json"),
        "reference_evaluation_freeze": pin("results/structural_pose_v3/evaluation_freeze.json"),
        "public_source_distribution": pin("benchmarks/source_public_distribution.json"),
        "source_profiles": {"strict": "strict_CCD_ideal", "repair": "modeled_sidechain_CCD_ideal"},
        "vina_carb_experiment": pin("benchmarks/structural_test_vina_carb.json"),
        "vina_carb_freeze": pin("results/structural_test_vina_carb/freeze_receipt.json"),
        "vina_carb_executable_sha256": "FILL_PRE_EXECUTION_BINARY_SHA256",
        "metric_code_sha256": {path: pin(path)["sha256"] for path in CODE_SOURCES},
        "freeze_receipt_schema": {"evaluation_plan_sha256": "SHA256_OF_COMPLETED_PLAN", "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False}}
