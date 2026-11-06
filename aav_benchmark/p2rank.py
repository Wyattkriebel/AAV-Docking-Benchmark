"""Frozen P2Rank predictions from an entire reference-free receptor structure.

The adapter never reads a ligand, reference, search box, or contact definition.
Unsuccessful execution and incomplete residue maps are explicit result records.
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time

from .provenance import sha256, write_json


LOCK_PATH = Path(__file__).resolve().parent.parent / "benchmarks" / "p2rank.lock.json"
DEFAULT_TIMEOUT_SECONDS = 600


def _table(path, required):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        if reader.fieldnames is None:
            raise ValueError("P2Rank CSV has no header")
        headers = [value.strip() for value in reader.fieldnames]
        if len(set(headers)) != len(headers) or not set(required).issubset(headers):
            raise ValueError("P2Rank CSV is missing required or has duplicate columns")
        reader.fieldnames = headers
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("P2Rank CSV has an incomplete or oversized row")
            rows.append({key: value.strip() for key, value in row.items()})
    return rows


def _finite(value, label):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"nonfinite P2Rank {label}")
    return number


def _label(label):
    match = re.fullmatch(r"(-?\d+)([^\s\d]?)", label)
    if match is None:
        raise ValueError(f"invalid P2Rank residue label {label!r}")
    return int(match[1]), match[2]


def _input_residues(receptor_path):
    """Index author IDs in PDB/mmCIF, without selecting chains or coordinates."""
    import gemmi

    structure = gemmi.read_structure(str(receptor_path))
    if len(structure) != 1:
        raise ValueError("P2Rank input must contain exactly one complete receptor model")
    residues = {}
    atom_count = 0
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if not all(math.isfinite(value) for value in (atom.pos.x, atom.pos.y, atom.pos.z)):
                    raise ValueError("input receptor has nonfinite coordinates")
                atom_count += 1
            if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                continue
            key = (chain.name, residue.seqid.num, residue.seqid.icode.strip())
            if key in residues:
                raise ValueError("input receptor has ambiguous duplicate author residue IDs")
            residues[key] = {"chain": key[0], "number": key[1], "insertion_code": key[2],
                             "resname": residue.name}
    if not residues or not atom_count:
        raise ValueError("input receptor has no amino-acid residues")
    return residues, atom_count


def parse_p2rank_predictions(predictions_path, residues_path, receptor_path):
    """Validate upstream CSVs and preserve the input's author residue identities.

    Residue assignment in ``*_residues.csv`` is exclusive, while pocket residue
    footprints in ``*_predictions.csv`` can overlap. The latter is the footprint
    returned for each pocket; the former validates that the receptor was read.
    """
    source, _ = _input_residues(receptor_path)
    residue_rows = _table(residues_path, ("chain", "residue_label", "residue_name", "pocket"))
    observed = {}
    for row in residue_rows:
        number, insertion = _label(row["residue_label"])
        key = (row["chain"], number, insertion)
        if key in observed or key not in source:
            raise ValueError("P2Rank residue map is duplicate or absent from input receptor")
        if source[key]["resname"] != row["residue_name"]:
            raise ValueError("P2Rank residue name differs from the input receptor")
        if int(row["pocket"]) < 0:
            raise ValueError("negative P2Rank residue pocket assignment")
        observed[key] = source[key]
    if set(observed) != set(source):
        raise ValueError("P2Rank did not report all input amino-acid residues")
    by_token = {f"{chain}_{number}{insertion}": value
                for (chain, number, insertion), value in observed.items()}
    if len(by_token) != len(observed):
        raise ValueError("ambiguous P2Rank residue token map")
    rows = _table(predictions_path, ("rank", "score", "probability", "center_x", "center_y",
                                     "center_z", "residue_ids", "surf_atom_ids", "sas_points"))
    predictions = []
    previous_score = math.inf
    for expected_rank, row in enumerate(rows, 1):
        rank = int(row["rank"])
        score = _finite(row["score"], "score")
        probability = _finite(row["probability"], "probability")
        if rank != expected_rank or score > previous_score or not 0 <= probability <= 1:
            raise ValueError("P2Rank pockets have invalid ranks, score order, or probabilities")
        previous_score = score
        tokens = row["residue_ids"].split()
        if not tokens or len(tokens) != len(set(tokens)) or any(token not in by_token for token in tokens):
            raise ValueError("P2Rank pocket residue footprint has an unmapped or duplicate ID")
        atom_ids = [int(value) for value in row["surf_atom_ids"].split()]
        sas_points = int(row["sas_points"])
        if sas_points < 0 or any(value < 0 for value in atom_ids):
            raise ValueError("invalid P2Rank surface point/atom metadata")
        predictions.append({"rank": rank, "center": [_finite(row["center_" + axis], "center")
                                                       for axis in "xyz"],
                            "score": score, "probability": probability,
                            "residues": [by_token[token] for token in tokens],
                            "atom_ids": atom_ids, "sas_point_count": sas_points})
    if any(int(row["pocket"]) > len(predictions) for row in residue_rows):
        raise ValueError("P2Rank residue assignments refer to an absent pocket")
    return predictions


def _installation(executable, java_home=None):
    lock = json.loads(LOCK_PATH.read_text())
    installation = executable.parent
    for relative, expected in lock["resource_sha256"].items():
        path = installation / relative
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"P2Rank installation checksum mismatch: {relative}")
    if executable.name != "prank" or sha256(executable) != lock["resource_sha256"]["prank"]:
        raise ValueError("P2Rank requires the pinned official 2.5.1 launcher")
    home = Path(java_home).resolve() if java_home else LOCK_PATH.parent.parent / lock["runtime"]["java_home"]
    for relative, expected in lock["runtime"]["resource_sha256"].items():
        path = home / relative
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"P2Rank portable Java runtime checksum mismatch: {relative}")
    environment = os.environ.copy()
    environment["JAVA_HOME"] = str(home)
    # Remove user-injected JVM options: the unmodified official launcher fixes
    # the maximum heap at 2048 MiB and the frozen model/config fix prediction.
    for key in ("JAVA_OPTS", "JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS", "_JAVA_OPTIONS"):
        environment.pop(key, None)
    version = subprocess.run([str(executable), "-v"], env=environment, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30, check=True).stdout
    if version.splitlines()[0] != "P2Rank 2.5.1":
        raise ValueError("P2Rank reported a different version")
    return lock, environment, version


def _terminate_group(process):
    """A timeout must terminate Java as well as its shell launcher."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run_p2rank(receptor_pdb, output_dir, executable, *, timeout_seconds=DEFAULT_TIMEOUT_SECONDS,
               java_home=None):
    """Run pinned P2Rank and archive logs, CSVs, hashes, timing, and failures.

    ``receptor_pdb`` also accepts a full mmCIF receptor. It must already be
    prepared without ligand-derived chain selection or cropping. Output must be
    fresh: existing output files are never treated as a completed prediction.
    ``timeout_seconds`` is exposed for operational cancellation/testing; cohort
    comparisons must use the frozen 600 second limit for every receptor.
    """
    receptor = Path(receptor_pdb).resolve()
    output = Path(output_dir).resolve()
    executable = Path(executable).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("P2Rank output directory must be empty to prevent stale predictions")
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("P2Rank timeout must be finite and positive")
    start = time.perf_counter()
    result = {"schema_version": 1, "method": "p2rank", "version": "2.5.1", "status": "failed",
              "input_path": str(receptor), "input_sha256": None, "executable_sha256": None,
              "command": [], "runtime_seconds": 0.0, "predictions": [], "artifacts": [],
              "reference_independent": True, "input_scope": "entire supplied prepared receptor",
              "residue_identity_system": "input author chain, number, insertion code, residue name",
              "lock_sha256": sha256(LOCK_PATH), "timeout_seconds": timeout_seconds}
    phase = "input_validation"
    try:
        result["input_sha256"] = sha256(receptor)
        residues, atoms = _input_residues(receptor)
        result.update({"input_amino_acid_residue_count": len(residues), "input_atom_count": atoms})
        phase = "installation_validation"
        result["executable_sha256"] = sha256(executable)
        lock, environment, version = _installation(executable, java_home)
        (output / "version.log").write_text(version)
        result["parameters"] = lock["parameters"]
        raw = output / "raw"
        command = [str(executable), "predict", "-f", str(receptor), "-o", str(raw),
                   "-threads", "1", "-visualizations", "0", "-seed", "42", "-model", "default"]
        result["command"] = command
        result["java_home"] = environment["JAVA_HOME"]
        phase = "execution"
        with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=environment,
                                       start_new_session=True)
            try:
                result["returncode"] = process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                _terminate_group(process)
                result["returncode"] = process.returncode
                result["status"] = "timeout"
                raise TimeoutError(f"P2Rank exceeded {timeout_seconds:g} seconds")
        if result["returncode"]:
            raise RuntimeError(f"P2Rank exited with code {result['returncode']}")
        phase = "output_validation"
        predictions = raw / (receptor.name + "_predictions.csv")
        residue_csv = raw / (receptor.name + "_residues.csv")
        result["predictions"] = parse_p2rank_predictions(predictions, residue_csv, receptor)
        if sha256(receptor) != result["input_sha256"]:
            raise ValueError("P2Rank input receptor changed during prediction")
        result["status"] = "ok"
        result["prediction_count"] = len(result["predictions"])
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        result["failure_category"] = "timeout" if result["status"] == "timeout" else phase
        result["reason"] = str(error)
        result["predictions"] = []
        (output / "adapter_error.log").write_text(f"{type(error).__name__}: {error}\n")
    finally:
        result["runtime_seconds"] = time.perf_counter() - start
        result["artifacts"] = [{"path": str(path.relative_to(output)), "sha256": sha256(path),
                                 "size_bytes": path.stat().st_size}
                                for path in sorted(output.rglob("*")) if path.is_file()]
        write_json(output / "result.json", result)
    return result


def _extract_distribution(archive, destination):
    """Extract a pinned one-directory distribution without replacing local work."""
    import shutil
    import tarfile
    import tempfile

    destination = Path(destination)
    if destination.exists():
        return
    with tempfile.TemporaryDirectory(prefix=".p2rank-install-", dir=destination.parent) as temporary:
        with tarfile.open(archive) as bundle:
            roots = {Path(member.name).parts[0] for member in bundle.getmembers()
                     if member.name and Path(member.name).parts}
            if len(roots) != 1:
                raise ValueError("P2Rank/runtime archive must contain exactly one top-level directory")
            bundle.extractall(temporary, filter="data")
        extracted = Path(temporary) / next(iter(roots))
        if not extracted.is_dir():
            raise ValueError("P2Rank/runtime archive did not contain its distribution directory")
        shutil.move(str(extracted), destination)


def install_p2rank(project_root=None):
    """Install the pinned official release and portable JRE under local tools/.

    Downloads are checked before extraction. Existing installations are verified
    and never overwritten or repaired. The frozen portable runtime is macOS
    arm64; a different platform needs its own reviewed runtime lock.
    """
    import platform
    from .provenance import fetch_verified

    if platform.system() != "Darwin" or platform.machine() not in {"arm64", "aarch64"}:
        raise ValueError("the pinned portable P2Rank runtime requires macOS arm64")
    root = Path(project_root).resolve() if project_root else LOCK_PATH.parent.parent
    lock = json.loads(LOCK_PATH.read_text())
    tools = root / "tools"
    tools.mkdir(parents=True, exist_ok=True)
    archive = fetch_verified(lock["archive"]["url"], tools / "p2rank_2.5.1.tar.gz", lock["archive"]["sha256"])
    runtime_archive = fetch_verified(lock["runtime"]["archive_url"],
                                     tools / "temurin-jre-21.0.13+9-macos-aarch64.tar.gz",
                                     lock["runtime"]["archive_sha256"])
    installation = root / lock["installed_directory"]
    java_home = root / lock["runtime"]["java_home"]
    _extract_distribution(archive, installation)
    # The macOS JRE archive contains <distribution>/Contents/Home.
    _extract_distribution(runtime_archive, java_home.parent.parent)
    executable = root / lock["executable"]
    _, _, version = _installation(executable, java_home)
    return {"status": "verified", "executable": str(executable), "java_home": str(java_home),
            "version_output": version, "lock_sha256": sha256(LOCK_PATH),
            "archive_sha256": sha256(archive), "runtime_archive_sha256": sha256(runtime_archive)}
