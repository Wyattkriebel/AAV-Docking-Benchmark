#!/usr/bin/env python3
"""Replay immutable source and box inputs without predictors, docking or metrics."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.aav import _pdb_atoms
from aav_benchmark.blind_docking import translate_ligand
from aav_benchmark.heldout_docking import LEDGER_SHA256, neighborhood_protein
from aav_benchmark.heldout_pose_inputs import validate_candidate_reuse
from scripts.verify_serialization_archives import read_archive, safe_name

MANIFEST = "results/prepared_input_foundation/prepared_inputs_manifest.json"
MANIFEST_SHA256 = "359ac7ac03395d7fc1fb9e666450882bf611ace3f12cafad9a95b2a5ee27d311"
CASES = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
METHODS = ("fpocket", "p2rank", "oracle")
PROFILE_CASES = {"strict": CASES, "repair": CASES, "strict_monomer": ("3ZWE", "5HZA"), "repair_monomer": ("3ZWE", "5HZA")}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def pinned(root, record):
    name = safe_name(record["path"])
    path = root / name
    if path.is_symlink() or not path.is_file() or any(p.is_symlink() for p in path.parents if p != root):
        raise ValueError("Input evidence must be an ordinary pinned file")
    raw = path.read_bytes()
    if sha(raw) != record["sha256"] or ("bytes" in record and len(raw) != record["bytes"]):
        raise ValueError(f"Input evidence SHA256 or size changed: {name}")
    return raw


def historical_path(name, historical_root):
    path = Path(name)
    if path.is_absolute():
        try:
            path = path.relative_to(historical_root)
        except ValueError as error:
            raise ValueError("Historical input escaped the declared repository root") from error
    return safe_name(path.as_posix())


def verify_profiles(profiles):
    if set(profiles) != set(PROFILE_CASES):
        raise ValueError("Keep all four fixed source profiles")
    shared = ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json", "chemistry.json", "initializer_report.json")
    for profile, cases in PROFILE_CASES.items():
        source = profiles[profile]
        record = json.loads(source["profile_manifest.json"])
        if record["profile"] != profile or record["case_ids"] != list(cases):
            raise ValueError("Source profile changed its fixed case scope")
        case_records = {row["case_id"]: row for row in record["cases"]}
        if len(case_records) != len(record["cases"]) or set(case_records) != set(cases):
            raise ValueError("Source profile lost or duplicated a case receipt")
        bindings = {"canonical_atom_map_sha256": "ligand_atom_map.json", "chemistry_sha256": "chemistry.json",
            "initializer_report_sha256": "initializer_report.json", "ligand_sdf_sha256": "ligand.sdf",
            "ligand_sha256": "ligand.pdbqt", "receptor_PDB_sha256": "receptor.pdb", "receptor_sha256": "receptor.cif",
            "observed_protein_atom_map_sha256": "protein_atom_map.json", "modeled_atom_map_sha256": "modeled_atom_map.json",
            "full_experimental_assembly_sha256": "full_assembly_receptor.cif", "full_observed_protein_map_sha256": "full_protein_atom_map.json"}
        for cid in cases:
            for key, name in bindings.items():
                if key in case_records[cid] and case_records[cid][key] != sha(source[f"{cid}/{name}"]):
                    raise ValueError("Source profile case receipt differs from archived atom/chemical inputs")
            for name in shared:
                if source[f"{cid}/{name}"] != profiles["strict"][f"{cid}/{name}"]:
                    raise ValueError("Paired profiles changed shared independent ligand bytes")
            if profile == "repair" and source[f"{cid}/protein_atom_map.json"] != profiles["strict"][f"{cid}/protein_atom_map.json"]:
                raise ValueError("Repair changed the original observed protein map")
            if "monomer" in profile:
                full = profiles["strict"]
                if (source[f"{cid}/full_assembly_receptor.cif"] != full[f"{cid}/receptor.cif"]
                        or source[f"{cid}/full_protein_atom_map.json"] != full[f"{cid}/protein_atom_map.json"]):
                    raise ValueError("Monomer changed retained full observed assembly inputs")
                observed = json.loads(source[f"{cid}/protein_atom_map.json"])["atoms"]
                all_observed = json.loads(full[f"{cid}/protein_atom_map.json"])["atoms"]
                if observed != [a for a in all_observed if a["output_key"][0] == "A"]:
                    raise ValueError("Monomer map is not the exact source-chain A subset")


def verify_prediction(files, condition, profiles, manifest):
    cid, method, profile = condition["case_id"], condition["method"], condition["profile"]
    pin = condition["primary_candidate_selection"]
    name = historical_path(pin["path"], manifest["historical_repository_root"])
    if sha(files[name]) != pin["sha256"]:
        raise ValueError("Original candidate-selection bytes changed")
    selection = json.loads(files[name])
    if condition["candidate_selection"] != selection:
        raise ValueError("Prepared condition changed original candidate ranks or centers")
    primary = "strict_monomer" if "monomer" in profile else "strict"
    expected_input = sha(profiles[primary][f"{cid}/receptor.cif"])
    config = "benchmarks/structural_assembly_site_experiment.json" if "monomer" in profile else "benchmarks/structural_site_v2_experiment.json"
    config_sha = next(r["sha256"] for r in manifest["frozen_documents"] if r["path"] == config)
    candidates = validate_candidate_reuse(selection, cid, method, config_sha, expected_input)
    directory = Path(name).parent
    record_pin = selection["prediction_record"]
    if Path(record_pin["filename"]).name != record_pin["filename"]:
        raise ValueError("Predictor record moved outside its original input directory")
    raw = files[(directory / record_pin["filename"]).as_posix()]
    if sha(raw) != record_pin["sha256"]:
        raise ValueError("Original predictor record changed")
    record = json.loads(raw)
    if record["input_sha256"] != expected_input or record["status"] != selection["predictor_status"]:
        raise ValueError("Original predictor input or status changed")
    artifacts = record["artifacts"]
    entries = [{"path": p, **v} for p, v in artifacts.items()] if isinstance(artifacts, dict) else artifacts
    for artifact in entries:
        raw = files[(directory / safe_name(artifact["path"])).as_posix()]
        if sha(raw) != artifact["sha256"] or len(raw) != artifact.get("bytes", artifact.get("size_bytes")):
            raise ValueError("Original predictor raw rank evidence changed")
    return candidates


def verify_contents(files, profiles, manifest):
    """Replay source-only copying, rigid translations and whole-chain input selection."""
    verify_profiles(profiles)
    conditions, statuses, predictions = [], Counter(), set()
    for index_pin in manifest["indexes"]:
        raw = files[index_pin["path"]]
        if sha(raw) != index_pin["sha256"]:
            raise ValueError("Preparation index SHA256 changed")
        index = json.loads(raw)
        if index["accuracy_evaluated"] is not False or index["no_docking_started"] is not True or index["status"] != "complete":
            raise ValueError("Index must remain an immutable input-only snapshot")
        if len(index["conditions"]) != index_pin["conditions"]:
            raise ValueError("Preparation index lost a condition")
        for row in index["conditions"]:
            raw = files[row["condition_path"]]
            if sha(raw) != row["condition_sha256"]:
                raise ValueError("Prepared condition differs from its frozen index")
            condition = json.loads(raw)
            config_name = "benchmarks/structural_assembly_pose_v3_experiment.json" if "monomer" in condition["profile"] else "benchmarks/structural_pose_v3_experiment.json"
            config_sha = next(r["sha256"] for r in manifest["frozen_documents"] if r["path"] == config_name)
            if (condition["profile"] != row["receptor_profile"] or condition["case_id"] != row["case_id"]
                    or condition["method"] != row["method"] or condition["status"] != row["status"]
                    or condition["fixed_denominator_member"] is not True or condition["executed_primary_candidate_budget"] != 1
                    or condition["configuration_sha256"] != config_sha or condition["source_ledger_sha256"] != LEDGER_SHA256
                    or ("monomer" not in condition["profile"] and condition.get("source_profile_manifest_sha256")
                        != sha(profiles[condition["profile"]]["profile_manifest.json"]))):
                raise ValueError("Fixed input condition identity or denominator changed")
            conditions.append((row, condition))
    expected = {(p, cid, m) for p, cases in PROFILE_CASES.items() for cid in cases for m in METHODS}
    identities = [(c["profile"], c["case_id"], c["method"]) for _, c in conditions]
    if len(identities) != 48 or len(set(identities)) != 48 or set(identities) != expected:
        raise ValueError("Keep all36 full and12 monomer input conditions")
    prepared = failures = empty = 0
    with tempfile.TemporaryDirectory(prefix="prepared-input-byte-replay-") as temporary:
        temporary = Path(temporary)
        for row, condition in conditions:
            profile, cid, method = condition["profile"], condition["case_id"], condition["method"]
            source = profiles[profile]
            if condition["native_coordinate_access"] is not (method == "oracle"):
                raise ValueError("Predicted/oracle information labels changed")
            if method == "oracle":
                if "primary_candidate_selection" in condition or len(condition["jobs"]) != 1:
                    raise ValueError("Oracle input must remain a labeled single supplied box")
                candidates = [{"rank": 1, "original_rank": 1, "center_A": condition["jobs"][0]["center_A"]}]
            else:
                candidates = verify_prediction(files, condition, profiles, manifest)
                predictions.add(("monomer" if "monomer" in profile else "full", cid, method))
            if (condition["available_prediction_candidates"] != len(candidates) or len(condition["jobs"]) != min(1, len(candidates))
                    or condition["scheduled_candidate_jobs"] != len(condition["jobs"])
                    or row["prepared_jobs"] != sum(j["status"] == "prepared" for j in condition["jobs"])):
                raise ValueError("Retain all candidates and only the declared top-one prepared input")
            statuses[condition["status"]] += 1
            if not condition["jobs"]:
                if condition["status"] != "completed_no_candidates":
                    raise ValueError("Empty candidate workflow lost its failure mechanism")
                empty += 1
            for job in condition["jobs"]:
                candidate = candidates[0]
                if (job["prepared_directory"] != "box01" or job["id"] != "box01"
                        or job["center_A"] != candidate["center_A"] or job["candidate_class_rank"] != candidate["rank"]
                        or job["candidate_original_rank"] != candidate["original_rank"]):
                    raise ValueError("Prepared input moved or reranked its original candidate")
                directory = Path(row["condition_path"]).parent / job["prepared_directory"]
                if job["status"] != "prepared":
                    if condition["status"] != "preparation_failed" or not json.loads(files[(directory / "preparation_attempt.json").as_posix()]):
                        raise ValueError("Unprepared condition lost its failure receipt")
                    failures += 1
                    continue
                raw = files[(directory / "preparation.json").as_posix()]
                if sha(raw) != job["preparation_sha256"]:
                    raise ValueError("Box preparation SHA256 changed")
                prep = json.loads(raw)
                chemistry = json.loads(source[f"{cid}/chemistry.json"])
                mapping = json.loads(source[f"{cid}/ligand_atom_map.json"])
                reference = json.loads(files[(directory / "reference.json").as_posix()])
                if (prep["box_center"] != job["center_A"] or prep["box_size"] != chemistry["box_size"]
                        or prep["configuration_sha256"] != condition["configuration_sha256"]
                        or (prep["profile"], prep["case_id"], prep["method"]) != (profile, cid, method)
                        or ("monomer" not in profile and prep.get("source_profile_manifest_sha256") != sha(source["profile_manifest.json"]))
                        or ("monomer" not in profile and prep.get("initializer_report_sha256") != sha(source[f"{cid}/initializer_report.json"]))
                        or prep["information_regime"] != ("true-site-local" if method == "oracle" else "predicted-site-local")
                        or prep.get("modeled_protein_atoms_are_atomic_truth", prep.get("modeled_atoms_are_atomic_truth")) is not False
                        or prep["source_receptor_sha256"] != sha(source[f"{cid}/receptor.cif"])
                        or prep["independent_ligand_sha256"] != sha(source[f"{cid}/ligand.pdbqt"])
                        or prep["native_coordinates_in_runner_reference"] is not False or reference["pose_reference_eligible"] is not False
                        or set(reference) != {"canonical_atom_map", "names", "pose_reference_eligible", "purpose"}
                        or reference["canonical_atom_map"] != mapping["atoms"]
                        or reference["names"] != [a["name"] for a in mapping["atoms"]]
                        or json.loads(files[(directory / "ligand_atom_map.json").as_posix()]) != mapping):
                    raise ValueError("Prepared box changed independent source bytes or full chemical identities")
                moved, _ = translate_ligand(source[f"{cid}/ligand.pdbqt"].decode(), job["center_A"])
                if moved.encode() != files[(directory / "ligand.pdbqt").as_posix()]:
                    raise ValueError("Prepared ligand differs from its exact frozen rigid translation")
                receptor = temporary / f"{profile}_{cid}.cif"
                if not receptor.exists():
                    receptor.write_bytes(source[f"{cid}/receptor.cif"])
                pdb, audit = neighborhood_protein(receptor, prep["box_center"], prep["box_size"])
                original = files[(directory / "neighborhood.pdb").as_posix()]
                if pdb.encode() != original or audit["selected_chains"] != prep["neighborhood"]["selected_chains"] or audit["interaction_cutoff_A"] != 8:
                    raise ValueError("Prepared neighborhood differs from whole-chain8A input selection")
                before = {k: a["coordinates"] for k, a in _pdb_atoms(original.decode()).items() if not a["hydrogen"]}
                typed = files[(directory / "neighborhood.pdbqt").as_posix()]
                after = {k: a["coordinates"] for k, a in _pdb_atoms(typed.decode(), pdbqt=True).items() if not a["hydrogen"]}
                if before != after:
                    raise ValueError("Prepared PDBQT changed or lost receptor heavy atoms")
                prepared += 1
    if dict(statuses) != manifest["condition_status_counts"] or (prepared, failures, empty) != (28, 18, 2) or len(predictions) != 16:
        raise ValueError("Input preparation/failure/candidate inventory changed")
    return {"conditions": 48, "full_assembly_conditions": 36, "monomer_conditions": 12, "prepared_boxes": prepared,
            "preparation_failed_conditions": failures, "no_candidate_conditions": empty, "original_predictor_streams": len(predictions),
            "source_profiles": 4, "archived_input_files": len(files)}


def verify(root):
    root = Path(root).resolve()
    manifest = json.loads(pinned(root, {"path": MANIFEST, "sha256": MANIFEST_SHA256}))
    if manifest["source_only"] is not True or manifest["no_new_predictions_docking_or_accuracy"] is not True:
        raise ValueError("Input packaging metadata changed its scope")
    for path, digest in manifest["frozen_producers"].items():
        pinned(root, {"path": path, "sha256": digest})
    pinned(root, {"path": "benchmarks/heldout_reference_ledger.json", "sha256": LEDGER_SHA256})
    documents = {record["path"]: pinned(root, record) for record in manifest["frozen_documents"]}
    for config, receipt, barrier in (
        ("benchmarks/structural_pose_v3_experiment.json", "results/structural_pose_v3/freeze_receipt.json", "frozen_before_new_docking"),
        ("benchmarks/structural_assembly_pose_v3_experiment.json", "results/structural_assembly_pose_v3/freeze_receipt.json", "frozen_before_new_docking"),
        ("benchmarks/structural_site_v2_experiment.json", "results/structural_test_v2/freeze_receipt.json", "frozen_before_new_predictions"),
        ("benchmarks/structural_assembly_site_experiment.json", "results/structural_assembly_site/freeze_receipt.json", "frozen_before_monomer_predictions")):
        freeze = json.loads(documents[receipt])
        if freeze["configuration_sha256"] != sha(documents[config]) or freeze[barrier] is not True:
            raise ValueError("Input design differs from its original scientific freeze")
    for name in ("results/structural_repair/prepared_sources_v2_manifest.json", "results/structural_assembly_site/source_archive_manifest.json"):
        record = json.loads(documents[name])
        archive_name = record.get("archive", "source_preparation.tar.gz")
        read_archive(documents[(Path(name).parent / archive_name).as_posix()], record, expected_count=len(record["files"]))
    archive = Path(MANIFEST).parent / manifest["archive"]
    files = read_archive(pinned(root, {"path": archive.as_posix(), "sha256": manifest["archive_sha256"]}), manifest, expected_count=698)
    if any({"vina", "vina_carb", "run", "common_pose_union", "supplemental_candidates_20261103"} & set(Path(p).parts) for p in files):
        raise ValueError("Input archive contains a protected or native-output path")
    profiles = {}
    for profile, records in zip(PROFILE_CASES, manifest["source_archives"], strict=True):
        source_manifest = json.loads(pinned(root, records["manifest"]))
        profiles[profile] = read_archive(pinned(root, records["archive"]), source_manifest, expected_count=len(source_manifest["files"]))
        design_path = "benchmarks/structural_assembly_pose_v3_experiment.json" if "monomer" in profile else "benchmarks/structural_pose_v3_experiment.json"
        descriptor = json.loads(documents[design_path])["source_profiles"][profile]
        if (descriptor["manifest"]["sha256"] != sha(profiles[profile]["profile_manifest.json"])
                or descriptor["source_archive_manifest"] != records["manifest"]):
            raise ValueError("Source archive differs from its original frozen design descriptor")
    return verify_contents(files, profiles, manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        result = verify(args.root)
    except (ValueError, KeyError, OSError, RuntimeError) as error:
        parser.exit(1, f"Prepared input verification failed: {error}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
