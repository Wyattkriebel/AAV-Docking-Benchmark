"""Reference-independent receptor preparation and audited fpocket execution.

Discovery functions have no native-ligand, site-center, or neighborhood input.
The full source receptor is stripped by chemical category, never by distance.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import resource
import shutil
import subprocess
import tarfile
import time

import gemmi
import numpy as np

from .aav import PROTEIN_RESIDUES, materialize_assembly
from .provenance import fetch_verified, sha256, write_json, write_text_atomic


DEFAULT_LOCK = Path(__file__).resolve().parents[1] / "benchmarks/site_methods.lock.json"


def _lock(path=DEFAULT_LOCK):
    return json.loads(Path(path).read_text())


def _identity(chain, residue, atom):
    return (chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, atom.name)


def receptor_atoms(path):
    """Read exact author atom identities and coordinates from a one-model input."""
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("site input must contain exactly one receptor model")
    atoms = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                key = _identity(chain, residue, atom)
                if key in atoms:
                    raise ValueError(f"duplicate receptor atom identity: {key}")
                xyz = np.array(list(atom.pos), dtype=float)
                if not np.isfinite(xyz).all():
                    raise ValueError(f"nonfinite receptor coordinates: {key}")
                atoms[key] = xyz
    if not atoms:
        raise ValueError("empty receptor input")
    return atoms


def prepare_site_receptor(source, output_path, *, assembly_id=None, lock_path=DEFAULT_LOCK):
    """Write the entire protein heavy-atom receptor, plus a source atom map.

    ``assembly_id='1'`` uses the existing reviewed AAV identity/60-mer audit.
    With ``None`` the whole deposited receptor model is retained; this is not a
    certification that a third-party benchmark's receptor is a full assembly.
    Use mmCIF for large assemblies. PDB serial or chain overflow is rejected.
    """
    source, output_path = Path(source), Path(output_path)
    if output_path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError("site receptor output must be PDB or mmCIF")
    manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
    map_path = output_path.with_suffix(output_path.suffix + ".atoms.tsv.gz")
    if any(path.exists() for path in (output_path, map_path, manifest_path)):
        raise FileExistsError("preserve existing site inputs; choose a fresh output path")
    structure = gemmi.read_structure(str(source))
    source_has_label_ids = structure.input_format != gemmi.CoorFormat.Pdb
    if len(structure) != 1:
        raise ValueError("curate an explicit single source model before site preparation")
    if assembly_id is None:
        model = structure[0].clone()
        assembly = {"assembly_id": None, "materializations": 0,
                    "scope": "entire deposited receptor model; biological completeness not certified"}
    elif assembly_id == "1":
        model, assembly = materialize_assembly(source)
    else:
        raise ValueError("only the reviewed AAV assembly 1 materialization is supported")
    source_model = model.clone()
    omitted = Counter()
    omitted_residues = Counter()
    before = {}
    map_rows = []
    chains = set()
    for ci in reversed(range(len(model))):
        chain = model[ci]
        for ri in reversed(range(len(chain))):
            residue = chain[ri]
            water = residue.is_water()
            metal_only = bool(len(residue)) and all(atom.element.is_metal for atom in residue)
            # Some benchmark PDBs encode calcium CA2 as ATOM. Classify whole
            # water/metal residues before interpreting the record type as protein.
            if water or metal_only or residue.het_flag != "A":
                if gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                    raise ValueError(f"nonstandard/HETATM amino acid requires curation: {chain.name}:{residue.seqid}:{residue.name}")
                omitted["water_atoms" if water else "nonprotein_atoms"] += len(residue)
                omitted["metal_atoms"] += sum(atom.element.is_metal for atom in residue)
                omitted_residues[residue.name] += 1
                del chain[ri]
                continue
            if residue.name not in PROTEIN_RESIDUES:
                raise ValueError(f"unsupported protein residue requires curation: {chain.name}:{residue.seqid}:{residue.name}")
            for ai in reversed(range(len(residue))):
                atom = residue[ai]
                if atom.element.is_hydrogen:
                    omitted["protein_hydrogens"] += 1
                    del residue[ai]
                    continue
                if atom.altloc not in ("\x00", " "):
                    raise ValueError("alternate conformers require independent curation before site preparation")
                if not math.isfinite(atom.occ) or atom.occ <= 0:
                    raise ValueError("nonpositive or nonfinite protein occupancy requires curation")
        if not any(len(residue) for residue in chain):
            del model[ci]
    for chain in model:
        chains.add(chain.name)
        for residue in chain:
            for atom in residue:
                key = _identity(chain, residue, atom)
                if key in before:
                    raise ValueError(f"duplicate protein atom identity: {key}")
                xyz = np.array(list(atom.pos), dtype=float)
                if not np.isfinite(xyz).all():
                    raise ValueError(f"nonfinite protein coordinates: {key}")
                before[key] = xyz
                map_rows.append([len(map_rows) + 1, *key, residue.subchain if source_has_label_ids else "",
                                 residue.label_seq if source_has_label_ids and residue.label_seq is not None else "",
                                 atom.serial, atom.element.name, *xyz.tolist()])
    if not before:
        raise ValueError("empty protein receptor after chemical-category stripping")
    if output_path.suffix.lower() == ".pdb":
        if len(before) > 99999 or any(len(c) != 1 for c in chains):
            raise ValueError("PDB identity/serial overflow; use mmCIF for the complete receptor")
        if any(not -999 <= key[1] <= 9999 for key in before):
            raise ValueError("PDB residue number overflow; use mmCIF")
        if any(np.any(xyz < -999.999) or np.any(xyz > 9999.999) for xyz in before.values()):
            raise ValueError("PDB coordinate overflow; use mmCIF")
    output = gemmi.Structure()
    output.name = "site_receptor"
    output.add_model(model)
    output.setup_entities()
    # BioJava requires numeric entity identifiers. Entity IDs are serialization
    # metadata, never author atom identities; original label IDs are in the map.
    for index, entity in enumerate(output.entities, 1):
        entity.name = str(index)
    for chain in output[0]:
        for index, residue in enumerate(chain, 1):
            if residue.label_seq is None:
                residue.label_seq = index
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".pdb":
        text = output.make_pdb_string()
    else:
        text = output.make_mmcif_document().as_string()
    write_text_atomic(output_path, text)
    after = receptor_atoms(output_path)
    if set(before) != set(after):
        output_path.unlink()
        raise ValueError("site receptor serialization changed protein atom identities")
    max_displacement = max(float(np.linalg.norm(before[key] - after[key])) for key in before)
    if max_displacement > 0.001:
        output_path.unlink()
        raise ValueError("site receptor serialization moved protein heavy atoms")
    # A compressed map preserves source serials plus author and label IDs without
    # bloating the committed manifest for a 247,860-atom AAV assembly.
    with map_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
            import io
            with io.TextIOWrapper(stream, newline="") as text_stream:
                writer = csv.writer(text_stream, delimiter="\t", lineterminator="\n")
                writer.writerow(["input_index", "chain", "number", "insertion_code", "resname", "atom_name",
                                 "source_label_asym_id", "source_label_seq_id", "source_serial", "element", "x", "y", "z"])
                writer.writerows(map_rows)
    manifest = {
        "schema_version": 1, "status": "prepared", "information_regime": "full_surface_site_agnostic",
        "source": str(source), "source_sha256": sha256(source),
        "receptor": str(output_path), "receptor_sha256": sha256(output_path),
        "method_lock_sha256": sha256(lock_path), "assembly": assembly,
        "protein_subunits": len(chains), "protein_heavy_atoms": len(before),
        "source_atoms_before_stripping": sum(len(r) for c in source_model for r in c),
        "omitted_atoms": dict(omitted), "omitted_residue_counts": dict(omitted_residues),
        "atom_map": str(map_path), "atom_map_sha256": sha256(map_path),
        "max_serialization_displacement_A": max_displacement,
        "policy": _lock(lock_path)["receptor_policy"], "gemmi_version": importlib.metadata.version("gemmi"),
        "native_ligand_access": False, "distance_based_selection": False,
        "serialization_metadata": "numeric entity IDs; sequential label_seq_id only where source has none; author IDs unchanged",
        "source_label_identifiers_present": source_has_label_ids,
    }
    write_json(manifest_path, manifest)
    return manifest


def install_fpocket(tools_dir, *, evidence_dir=None, lock_path=DEFAULT_LOCK):
    """Fetch and build unmodified official source, preserving compiler evidence.

    Only source, headers, platform-specific bundled molfile libraries and license
    files are extracted. Large test/data trees and unrelated platforms are not.
    This builder currently supports the documented Apple arm64 environment.
    """
    config = _lock(lock_path)["methods"]["fpocket"]
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise RuntimeError("the frozen builder supports macOS arm64; record a reviewed build recipe for this platform")
    tools_dir = Path(tools_dir).resolve()
    package = tools_dir / f"fpocket-{config['release']}.tar.gz"
    source_dir = tools_dir / f"fpocket-{config['commit']}"
    executable = source_dir / "bin/fpocket"
    evidence_dir = Path(evidence_dir).resolve() if evidence_dir else source_dir / "build_evidence"
    installed = source_dir / "installation.json"
    if installed.exists():
        existing = json.loads(installed.read_text())
        if existing["source_sha256"] != config["source_sha256"] or sha256(executable) != existing["executable_sha256"]:
            raise ValueError("installed fpocket differs from its recorded build")
        return existing
    if source_dir.exists():
        raise FileExistsError("unverified source/build directory exists; preserve it and choose a fresh tools directory")
    if evidence_dir != source_dir / "build_evidence":
        evidence_dir.mkdir(parents=True, exist_ok=True)
    fetch_verified(config["source_url"], package, config["source_sha256"])
    with tarfile.open(package, "r:gz") as archive:
        selected = []
        for member in archive.getmembers():
            parts = Path(member.name).parts
            if not parts or parts[0] != source_dir.name or ".." in parts or member.issym() or member.islnk():
                if member.issym() or member.islnk():
                    continue
                raise ValueError("unexpected source archive path")
            relative = parts[1:]
            keep = (not relative or len(relative) == 1 or relative[0] in {"src", "headers"}
                    or relative[:2] in (("plugins", "include"), ("plugins", "MACOSXARM64")))
            if keep:
                selected.append(member)
        archive.extractall(tools_dir, members=selected, filter="data")
    evidence_dir.mkdir(parents=True, exist_ok=True)
    for name, digest in config["license_files"].items():
        if sha256(source_dir / name) != digest:
            raise ValueError("upstream license content changed")
        shutil.copy2(source_dir / name, evidence_dir / name)
    recipe = config["build"]["macos_arm64"]
    command = ["make", "bin/fpocket", "-j2", *(f"{key}={recipe[key]}" for key in ("ARCH", "CFLAGS", "QCFLAGS"))]
    with (evidence_dir / "build.stdout.log").open("w") as stdout, (evidence_dir / "build.stderr.log").open("w") as stderr:
        result = subprocess.run(command, cwd=source_dir, stdout=stdout, stderr=stderr, check=False)
    if result.returncode or not executable.is_file():
        raise RuntimeError(f"fpocket compilation failed (exit {result.returncode}); retained {evidence_dir}")
    compiler = subprocess.run(["clang", "--version"], capture_output=True, text=True, check=False)
    write_text_atomic(evidence_dir / "compiler.txt", compiler.stdout + compiler.stderr)
    banner = subprocess.run([str(executable), "-h"], capture_output=True, text=True, check=False)
    write_text_atomic(evidence_dir / "version_banner.txt", banner.stdout + banner.stderr)
    manifest = {
        "schema_version": 1, "method": "fpocket", "release": config["release"], "commit": config["commit"],
        "source_url": config["source_url"], "source_sha256": config["source_sha256"],
        "executable": str(executable), "executable_sha256": sha256(executable),
        "build_command": command, "build_cwd": str(source_dir), "platform": platform.platform(),
        "algorithm_source_modifications": [], "method_lock_sha256": sha256(lock_path),
        "build_evidence": {p.name: sha256(p) for p in sorted(evidence_dir.iterdir()) if p.is_file()},
        "build_evidence_dir": str(evidence_dir),
    }
    write_json(installed, manifest)
    write_json(evidence_dir / "installation.json", manifest)
    return manifest


def _fpocket_atom_rows(path):
    """Read upstream's atom loop despite its unquoted prose before the loop.

    fpocket 4.2.3 pocketN_atm.cif includes invalid, unquoted descriptor prose in
    a preceding loop. Preserve that raw file; parse only its _atom_site loop.
    """
    lines = Path(path).read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() == "_atom_site.group_PDB"), None)
    if start is None:
        raise ValueError(f"missing fpocket atom loop: {path}")
    atom_loop = "data_pocket\nloop_\n" + "\n".join(lines[start:]) + "\n"
    table = gemmi.cif.read_string(atom_loop).sole_block().find_mmcif_category("_atom_site.")
    tags = list(table.tags)
    return [{tag.removeprefix("_atom_site."): str(row[i]) for i, tag in enumerate(tags)} for row in table]


def _row_atom(row):
    insertion = row.get("pdbx_PDB_ins_code", "?")
    insertion = "" if insertion in {"?", "."} else insertion
    key = (row.get("auth_asym_id", row["label_asym_id"]),
           int(row.get("auth_seq_id", row["label_seq_id"])), insertion,
           row["label_comp_id"], row["label_atom_id"])
    return key, np.array([float(row[f"Cartn_{axis}"]) for axis in "xyz"])


def audit_fpocket_receptor(output_path, receptor_path):
    """Prove upstream's full protein output contains every supplied atom."""
    expected = receptor_atoms(receptor_path)
    observed = {}
    for row in _fpocket_atom_rows(output_path):
        if row["group_PDB"] != "ATOM":
            continue  # upstream's extra HETATM/STP records are alpha spheres
        key, xyz = _row_atom(row)
        if key in observed:
            raise ValueError(f"duplicate fpocket full-receptor identity: {key}")
        observed[key] = xyz
    if set(expected) != set(observed):
        raise ValueError(f"fpocket full-receptor atom loss/addition: expected {len(expected)}, observed {len(observed)}")
    maximum = max(float(np.linalg.norm(expected[key] - observed[key])) for key in expected)
    if maximum > 0.001:
        raise ValueError("fpocket full-receptor output changed protein heavy-atom coordinates")
    return {"status": "passed", "expected_protein_heavy_atoms": len(expected),
            "observed_protein_heavy_atoms": len(observed), "max_displacement_A": maximum}


def parse_fpocket(output_dir, receptor_path):
    """Normalize upstream score-ranked centers and audited receptor footprints."""
    output_dir = Path(output_dir)
    receptor = receptor_atoms(receptor_path)
    pockets = []
    for vertices in sorted((output_dir / "pockets").glob("pocket*_vert.pqr"),
                           key=lambda p: int(re.fullmatch(r"pocket(\d+)_vert.pqr", p.name).group(1))):
        rank = int(re.fullmatch(r"pocket(\d+)_vert.pqr", vertices.name).group(1))
        atom_file = vertices.with_name(f"pocket{rank}_atm.cif")
        if not atom_file.is_file():
            raise ValueError(f"missing footprint for pocket {rank}")
        coordinates, score, count = [], None, None
        for line in vertices.read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elif "Pocket Score" in line:
                score = float(line.rsplit(":", 1)[1])
            elif "Number of V. Vertices" in line:
                count = int(line.rsplit(":", 1)[1])
        xyz = np.asarray(coordinates, dtype=float)
        if not len(coordinates) or xyz.shape[1:] != (3,) or not np.isfinite(xyz).all():
            raise ValueError(f"invalid alpha-sphere centers for pocket {rank}")
        if score is None or not math.isfinite(score) or count != len(coordinates):
            raise ValueError(f"incomplete alpha-sphere/score output for pocket {rank}")
        footprint = set()
        for row in _fpocket_atom_rows(atom_file):
            if row["group_PDB"] not in {"ATOM", "HETATM"}:
                continue
            key, position = _row_atom(row)
            if key not in receptor or not np.allclose(position, receptor[key], rtol=0, atol=0.001):
                raise ValueError(f"fpocket footprint lost source identity or coordinates: {key}")
            if key in footprint:
                raise ValueError(f"duplicate fpocket footprint identity: {key}")
            footprint.add(key)
        if not footprint:
            raise ValueError(f"empty receptor footprint for pocket {rank}")
        residues = sorted({key[:4] for key in footprint})
        pockets.append({
            "rank": rank, "center": xyz.mean(axis=0).tolist(), "center_A": xyz.mean(axis=0).tolist(), "score": score,
            "alpha_sphere_count": len(coordinates),
            "atom_ids": [list(key) for key in sorted(footprint)],
            "residues": [{"chain": c, "number": n, "insertion_code": i, "resname": r}
                         for c, n, i, r in residues],
            "residue_keys": [list(key) for key in residues],
            "raw_vertices": str(vertices.relative_to(output_dir)),
            "raw_footprint": str(atom_file.relative_to(output_dir)),
        })
    if pockets and [p["rank"] for p in pockets] != list(range(1, len(pockets) + 1)):
        raise ValueError("fpocket pocket ranks are incomplete")
    return pockets


def run_fpocket(receptor_path, output_dir, *, executable, timeout_seconds=3600, lock_path=DEFAULT_LOCK):
    """Execute the locked method once, retaining all evidence, even on failure."""
    receptor_path, output_dir, executable = Path(receptor_path), Path(output_dir), Path(executable).resolve()
    if output_dir.exists():
        raise FileExistsError("preserve prior site execution; choose a new output directory")
    config = _lock(lock_path)["methods"]["fpocket"]
    installation_path = executable.parent.parent / "installation.json"
    if not installation_path.is_file():
        raise ValueError("fpocket needs a pinned installation.json from install_fpocket")
    installation = json.loads(installation_path.read_text())
    binary_sha = sha256(executable)
    if any(installation.get(key) != config[key] for key in ("release", "commit", "source_sha256")) or installation.get("executable_sha256") != binary_sha:
        raise ValueError("fpocket executable/build provenance differs from the frozen method")
    # Do not permit raw ligand-containing inputs or earlier reference-selected
    # neighborhoods to masquerade as full-surface discovery inputs.
    prepared_path = receptor_path.with_suffix(receptor_path.suffix + ".manifest.json")
    prepared = json.loads(prepared_path.read_text())
    if (prepared.get("receptor_sha256") != sha256(receptor_path)
            or prepared.get("native_ligand_access") is not False
            or prepared.get("distance_based_selection") is not False):
        raise ValueError("site input differs from its reference-independent preparation manifest")
    for key in receptor_atoms(receptor_path):
        if key[3] not in PROTEIN_RESIDUES:
            raise ValueError("fpocket discovery baseline accepts only prepared protein atoms")
    output_dir.mkdir(parents=True)
    output_dir = output_dir.resolve()
    local_input = output_dir / ("receptor" + receptor_path.suffix.lower())
    shutil.copy2(receptor_path, local_input)
    shutil.copy2(prepared_path, output_dir / "receptor_preparation.json")
    command = [str(executable), "-f", local_input.name, *config["parameters"]]
    manifest = {
        "schema_version": 1, "method": "fpocket", "status": "running", "release": config["release"],
        "commit": config["commit"], "method_lock_sha256": sha256(lock_path),
        "input_sha256": sha256(local_input), "executable_sha256": binary_sha,
        "source_sha256": config["source_sha256"], "installation": installation,
        "command": command, "cwd": str(output_dir), "information_regime": "full_surface_site_agnostic",
        "center_definition": config["center_definition"], "ranking": config["ranking"],
        "randomness": config["randomness"], "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "predictions": [], "platform": platform.platform(), "logical_cpus": os.cpu_count(),
    }
    write_json(output_dir / "execution.json", manifest)
    started, usage_before = time.monotonic(), resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        with (output_dir / "stdout.log").open("w") as stdout, (output_dir / "stderr.log").open("w") as stderr:
            result = subprocess.run(command, cwd=output_dir, stdout=stdout, stderr=stderr,
                                    timeout=timeout_seconds, check=False)
        manifest["returncode"] = result.returncode
        if result.returncode:
            raise RuntimeError(f"fpocket exited with status {result.returncode}")
        raw_output = output_dir / "receptor_out"
        info = raw_output / "receptor_info.txt"
        # fpocket can return 0 after input/parser failure: output evidence is
        # mandatory. A documented no-pocket result is retained as zero sites.
        no_pockets = "no pockets found" in (output_dir / "stdout.log").read_text().splitlines()
        if not info.is_file() and not no_pockets:
            raise RuntimeError("fpocket returned zero but did not produce its descriptor file or explicit no-pockets result")
        if info.is_file():
            manifest["receptor_output_audit"] = audit_fpocket_receptor(raw_output / "receptor_out.cif", local_input)
        else:
            manifest["receptor_output_audit"] = {"status": "unavailable", "reason": "upstream writes no full receptor when no pockets are found"}
        predictions = parse_fpocket(raw_output, local_input) if info.is_file() else []
        described = len(re.findall(r"(?m)^Pocket\s+\d+\s*:", info.read_text())) if info.is_file() else 0
        if described != len(predictions):
            raise ValueError("fpocket descriptors and normalized pocket counts differ")
        manifest.update(status="completed", predictions=predictions, pocket_count=len(predictions),
                        no_pockets_reported=no_pockets)
        write_json(output_dir / "predictions.json", predictions)
    except Exception as error:
        # Malformed upstream artifacts and adapter failures belong in the same
        # retained attempt ledger as process failures; never leave "running".
        manifest.update(status="failed", error_type=type(error).__name__, error=str(error))
    finally:
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        manifest["runtime_seconds"] = time.monotonic() - started
        manifest["cpu_seconds"] = ((usage_after.ru_utime + usage_after.ru_stime)
                                   - (usage_before.ru_utime + usage_before.ru_stime))
        manifest["peak_memory_bytes"] = usage_after.ru_maxrss if platform.system() == "Darwin" else usage_after.ru_maxrss * 1024
        manifest["peak_memory_scope"] = "process-family cumulative child high-water mark; not isolated per-run RSS"
        manifest["ended_at_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["artifacts"] = {str(p.relative_to(output_dir)): {"sha256": sha256(p), "bytes": p.stat().st_size}
                                 for p in sorted(output_dir.rglob("*")) if p.is_file() and p.name != "execution.json"}
        write_json(output_dir / "execution.json", manifest)
    return manifest
