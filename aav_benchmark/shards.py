"""Merge independently executed frozen-seed runs without changing the study plan.

A shard is an execution subset, never a smaller scientific denominator. Stop
all writers before merging. Analysis remains unavailable until every configured
receptor/seed identity has one terminal run record. Commands and resource
measurements are retained from the actual worker processes.
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import re
import shutil

from .engine import read_poses
from .provenance import sha256, write_json


def _filename(name):
    if not isinstance(name, str) or not name or Path(name).name != name or name in {".", ".."}:
        raise ValueError("Shard artifact filenames must be simple filenames")
    return name


def _preserving_copy(source, destination):
    """Copy changed bytes only after moving the old artifact to a unique name."""
    source, destination = Path(source), Path(destination)
    if destination.exists() and sha256(destination) == sha256(source):
        return None
    previous = None
    if destination.exists():
        attempt = 1
        while destination.with_name(destination.name + f".previous-{attempt}").exists():
            attempt += 1
        previous = destination.with_name(destination.name + f".previous-{attempt}")
        destination.rename(previous)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return previous.name if previous is not None else None


def _validate_command(row, definition, preparation, receptor, prepared, directory):
    command = row.get("command")
    if not isinstance(command, list) or not command or any(not isinstance(value, str) for value in command):
        raise ValueError("Shard must retain its actual invocation argument list")

    def argument(flag):
        if command.count(flag) != 1 or command.index(flag) + 1 >= len(command):
            raise ValueError(f"Missing/duplicate shard invocation option {flag}")
        return command[command.index(flag) + 1]

    numeric = {"--seed": row["seed"], "--cpu": definition["cpu"],
               "--exhaustiveness": definition["exhaustiveness"], "--num_modes": definition["num_modes"],
               "--min_rmsd": definition["min_rmsd_A"], "--energy_range": definition["energy_range_kcal_mol"]}
    for axis, center, size in zip("xyz", preparation["box_center"], preparation["box_size"]):
        numeric["--center_" + axis], numeric["--size_" + axis] = center, size
    for flag, expected in numeric.items():
        value = float(argument(flag))
        if not math.isfinite(value) or value != float(expected):
            raise ValueError(f"Shard invocation changed frozen setting {flag}")
    if argument("--scoring") != "vina":
        raise ValueError("Shard scoring function differs from frozen Vina scoring")
    paths = {"--receptor": prepared / receptor["pdbqt_file"],
             "--ligand": prepared / preparation.get("ligand_file", "galactose.pdbqt"),
             "--out": directory / _filename(row["pose_file"])}
    if any(Path(argument(flag)).resolve() != path.resolve() for flag, path in paths.items()):
        raise ValueError("Shard invocation input/output paths differ from its artifacts")


def merge_shards(prepared, output, config, shard_dirs, baseline_runs=None):
    """Validate and assemble an exact full receptor-by-seed run manifest.

    ``baseline_runs`` optionally contains the old output manifest after its
    orchestrator was stopped. Terminal records contribute to coverage; interrupted
    records become explicit previous attempts on the replacement run. Missing
    interruption timings remain missing and are never filled with zero.

    Each shard must retain the *full* frozen ``configuration.json``. Input,
    executable and pose fingerprints must agree, logs must exist, and completed
    PDBQT atom maps are independently parsed before any destination is changed.
    Duplicate terminal identities are rejected even if their artifacts match.
    """
    prepared, output, config = (Path(value).resolve() for value in (prepared, output, config))
    definition = json.loads(config.read_text())
    seeds = definition["seeds"]
    if not seeds or len(set(seeds)) != len(seeds) or any(isinstance(s, bool) or not isinstance(s, int) for s in seeds):
        raise ValueError("Frozen seeds must be unique integers")
    preparation_path, reference_path = prepared / "preparation.json", prepared / "reference.json"
    preparation = json.loads(preparation_path.read_text())
    reference = json.loads(reference_path.read_text())
    receptors = {row["id"]: row for row in preparation["receptors"]}
    if not receptors or len(receptors) != len(preparation["receptors"]):
        raise ValueError("Prepared receptors must have unique nonempty identities")
    expected = {(receptor, seed) for receptor in receptors for seed in seeds}
    frozen_hash = sha256(config)
    full_file = reference.get("full_assembly_file")
    common = {
        "config_sha256": frozen_hash,
        "preparation_sha256": sha256(preparation_path), "reference_sha256": sha256(reference_path),
        "ligand_sha256": sha256(prepared / preparation.get("ligand_file", "galactose.pdbqt")),
        "full_assembly_sha256": sha256(prepared / full_file) if full_file else None,
    }
    directories = [Path(path).resolve() for path in shard_dirs]
    if len(directories) != len(set(directories)):
        raise ValueError("Duplicate shard directories")
    if (output / "runs.json").exists() and output not in directories and baseline_runs is None:
        raise ValueError("Existing output runs must be included as a shard or supplied as baseline_runs")
    sources = []
    for directory in directories:
        manifest = directory / "runs.json"
        if sha256(directory / "configuration.json") != frozen_hash:
            raise ValueError("Shard configuration differs from the full frozen configuration")
        rows = json.loads(manifest.read_text())
        sources.append((directory, manifest, sha256(manifest), rows))
    if baseline_runs is not None:
        manifest = output / "runs.json"
        if not manifest.is_file():
            raise ValueError("Baseline rows require the original output runs.json")
        if sha256(output / "configuration.json") != frozen_hash:
            raise ValueError("Baseline configuration differs from frozen configuration")
        original = json.loads(manifest.read_text())
        if baseline_runs != original:
            raise ValueError("Baseline rows must exactly preserve the original manifest checkpoint")
        sources.insert(0, (output, manifest, sha256(manifest), baseline_runs))
    if not sources:
        raise ValueError("No shard or baseline evidence supplied")
    completed, interrupted, artifacts = {}, {}, []
    executable_hashes = set()
    for source_index, (directory, manifest, manifest_hash, rows) in enumerate(sources):
        if not isinstance(rows, list) or not rows:
            raise ValueError("Shard run manifests must be nonempty lists")
        local = set()
        for raw in rows:
            row = copy.deepcopy(raw)
            identity = (row["receptor"], row["seed"])
            if identity not in expected or identity in local:
                raise ValueError("Unknown or duplicate run identity within a shard")
            local.add(identity)
            receptor = receptors[identity[0]]
            required = {**common,
                        "receptor_sha256": sha256(prepared / receptor["pdbqt_file"]),
                        "receptor_pdb_sha256": sha256(prepared / receptor["pdb_file"])}
            fingerprint = row.get("fingerprint", {})
            if any(fingerprint.get(key) != value for key, value in required.items()):
                raise ValueError("Shard prepared inputs or frozen configuration changed")
            executable_hash = fingerprint.get("executable_sha256")
            if not isinstance(executable_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", executable_hash):
                raise ValueError("Missing or invalid executable fingerprint")
            executable_hashes.add(executable_hash)
            for key in ("case_id", "family_id"):
                if key in preparation and row.get(key) != preparation[key]:
                    raise ValueError("Shard case/family identity differs from preparation")
            _validate_command(row, definition, preparation, receptor, prepared, directory)
            row["shard_provenance"] = {
                "source_directory": str(directory), "source_run_manifest": str(manifest),
                "source_run_manifest_sha256": manifest_hash,
                "archived_run_manifest": f"shard_manifests/source-{source_index}-{manifest_hash}.json",
            }
            if row["status"] in {"running", "pending", "cancelled", "interrupted"}:
                if baseline_runs is None or directory != output:
                    raise ValueError("A shard still has nonterminal work; stop its writer before merging")
                row["status_before_merge"] = row["status"]
                row["status"] = "cancelled"
                row["interruption_note"] = "Caller supplied stopped baseline; replaced by a separate frozen-seed execution"
                row["missing_resource_measurements"] = [key for key in ("wall_seconds", "cpu_seconds", "child_peak_rss_bytes") if key not in row]
                interrupted.setdefault(identity, []).append(row)
                continue
            if row["status"] not in {"ok", "failed", "unsupported"}:
                raise ValueError("Unknown terminal shard status")
            if identity in completed:
                raise ValueError("Duplicate terminal receptor-by-seed identity across shards")
            log_name = _filename(row["log_file"])
            log_path = directory / log_name
            if not log_path.is_file():
                raise ValueError("Shard execution log is missing")
            log_hash = sha256(log_path)
            if row.get("log_sha256", log_hash) != log_hash:
                raise ValueError("Shard execution log checksum mismatch")
            row["log_sha256"] = log_hash
            artifacts.append((log_path, output / log_name, identity, "log_file"))
            if row["status"] == "ok":
                pose_name = _filename(row["pose_file"])
                pose_path = directory / pose_name
                if not pose_path.is_file() or sha256(pose_path) != row.get("pose_sha256"):
                    raise ValueError("Shard pose artifact checksum mismatch")
                if len(read_poses(pose_path, reference["names"])) != row.get("returned_poses"):
                    raise ValueError("Shard pose count differs from parsed complete atom maps")
                artifacts.append((pose_path, output / pose_name, identity, "pose_file"))
            completed[identity] = row
    if len(executable_hashes) != 1:
        raise ValueError("Shards used different executable fingerprints")
    if set(completed) != expected:
        raise ValueError("Incomplete full frozen receptor-by-seed denominator")
    destinations = [destination for _, destination, _, _ in artifacts]
    if len(destinations) != len(set(destinations)):
        raise ValueError("Different shard identities reuse an output artifact filename")
    output.mkdir(parents=True, exist_ok=True)
    if (output / "configuration.json").exists() and sha256(output / "configuration.json") != frozen_hash:
        raise ValueError("Existing output configuration differs from frozen configuration")
    # Preserve the exact source manifests before replacing an output manifest.
    for index, (_, manifest, manifest_hash, _) in enumerate(sources):
        destination = output / "shard_manifests" / f"source-{index}-{manifest_hash}.json"
        _preserving_copy(manifest, destination)
    for source, destination, identity, artifact_type in artifacts:
        previous = _preserving_copy(source, destination)
        if previous:
            completed[identity].setdefault("preserved_previous_artifacts", {})[artifact_type] = previous
            for attempt in interrupted.get(identity, []):
                attempt.setdefault("preserved_artifacts", {})[artifact_type] = previous
    if config != output / "configuration.json":
        _preserving_copy(config, output / "configuration.json")
    result = []
    for receptor in preparation["receptors"]:
        for seed in seeds:
            identity = (receptor["id"], seed)
            row = completed[identity]
            if interrupted.get(identity):
                row["previous_attempts"] = row.get("previous_attempts", []) + interrupted[identity]
            result.append(row)
    write_json(output / "runs.json", result)
    return result
