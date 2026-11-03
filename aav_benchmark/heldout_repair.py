"""Separate source-only sidechain completion sensitivity for the fixed six.

Primary strict inputs remain immutable.  The isolated PDBFixer/OpenMM worker
sees only protein atoms, never carbohydrate references or docking results.
Modeled atoms are preparation assumptions and cannot become atomic truth.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

CASE_IDS = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
STANDARD_RESIDUES = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
BACKBONE = frozenset(("N", "CA", "C", "O", "OXT"))
SEED = 20261024


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text())


def _write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def parse_pdb_heavy(text):
    """Read exact source identities; no implicit atom selection or renaming."""
    atoms = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 78 or line[16:17] not in (" ", ""):
            raise ValueError("Require explicit element fields and preselected whole-residue alternatives")
        element = line[76:78].strip()
        if not element or element in ("H", "D"):
            raise ValueError("Repair input/output must contain only explicit heavy atoms")
        key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip(), line[12:16].strip())
        xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        if key in atoms or not all(math.isfinite(x) for x in xyz):
            raise ValueError("Duplicate or nonfinite source atom")
        atoms[key] = {"coordinates": xyz, "element": element, "line": line}
    if not atoms:
        raise ValueError("Empty source protein")
    return atoms


def audit_repair(before, after, allowed_added_keys, max_observed_displacement_A=.001):
    """Fail on any changed observed identity, element, or heavy coordinate."""
    if isinstance(before, str): before = parse_pdb_heavy(before)
    if isinstance(after, str): after = parse_pdb_heavy(after)
    for records in (before, after):
        for key, row in records.items():
            xyz = row["coordinates"]
            if (len(xyz) != 3 or not all(math.isfinite(float(x)) for x in xyz)
                    or not row.get("element") or row["element"] in ("H", "D")):
                raise ValueError(f"Repair atom must have finite coordinates and explicit heavy element: {key}")
    allowed = {tuple(k) for k in allowed_added_keys}
    deleted = set(before) - set(after)
    added = set(after) - set(before)
    if deleted:
        raise ValueError(f"Repair deleted observed heavy atoms: {sorted(deleted)}")
    if added != allowed:
        raise ValueError(f"Repair added atoms outside exact allowlist or failed to add allowed atoms: unexpected={sorted(added-allowed)}, missing={sorted(allowed-added)}")
    maximum = 0.
    for key, row in before.items():
        if row["element"] != after[key]["element"]:
            raise ValueError(f"Repair changed observed heavy element: {key}")
        delta = math.sqrt(sum((a-b)**2 for a, b in zip(row["coordinates"], after[key]["coordinates"])))
        maximum = max(maximum, delta)
    if maximum > max_observed_displacement_A:
        raise ValueError(f"Repair moved observed heavy atoms by {maximum:g} A")
    for key in added:
        if key[3] not in STANDARD_RESIDUES or key[4] in BACKBONE:
            raise ValueError(f"Repair attempted a nonstandard or backbone atom: {key}")
    return {"observed_heavy_atoms": len(before), "added_heavy_atoms": len(added),
            "max_observed_displacement_A": maximum, "observed_heavy_identity_and_elements_unchanged": True,
            "added_atoms": [{"key": list(key), "element": after[key]["element"], "coordinates_A": after[key]["coordinates"]}
                            for key in sorted(added)], "modeled_atoms_are_atomic_truth": False}


def split_modeling_segments(text):
    """Expose observed sequence gaps without inserting loops or guessing bonds.

    Temporary segment chain names affect modeling topology only.  Every output
    chain/residue ID is restored before preparation; full source chains remain
    the prediction and whole-chain box-selection units.
    """
    atoms = parse_pdb_heavy(text)
    residues = defaultdict(list)
    for key, row in atoms.items(): residues[key[:4]].append((key, row))
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    chain_map, residue_map, rows, previous, current = {}, [], [], None, None
    for residue_key, entries in residues.items():
        chain, number, insertion, resname = residue_key
        names = {key[-1]: row for key, row in entries}
        split, reason = previous is None, "source_chain_start"
        if previous is not None:
            pkey, pnames = previous
            consecutive = chain == pkey[0] and not insertion and not pkey[2] and number == pkey[1] + 1
            distance = (math.sqrt(sum((a-b)**2 for a, b in zip(pnames["C"]["coordinates"], names["N"]["coordinates"])))
                        if "C" in pnames and "N" in names else None)
            split = not consecutive or distance is None or not 1.1 - 1e-9 <= distance <= 1.7 + 1e-9
            reason = "source_chain_or_number_gap" if not consecutive else "missing_or_disconnected_peptide"
        if split:
            if len(chain_map) >= len(alphabet):
                raise ValueError("Modeling segmentation exceeds explicit PDB chain capacity")
            current = alphabet[len(chain_map)]; chain_map[current] = chain
            if rows: rows.append("TER")
        residue_map.append({"modeling_chain": current, "modeling_resnum": number, "modeling_insertion": insertion,
                            "source_chain": chain, "source_resnum": number, "source_insertion": insertion,
                            "resname": resname, "new_segment": split, "segment_reason": reason if split else None})
        for _, row in entries: rows.append(row["line"][:21] + current + row["line"][22:])
        previous = (residue_key, names)
    rows.extend(("TER", "END"))
    return "\n".join(rows) + "\n", {"chains": chain_map, "residues": residue_map,
                                             "source_only": True, "missing_loops_inserted": False}


def verify_tools(root, lock_path):
    """Check exact installed code/resources as well as official distribution pins."""
    root = Path(root).resolve(); lock = _load(lock_path)
    if lock.get("schema_version") != 1 or lock.get("platform") != "Reference" or lock.get("seed") != SEED:
        raise ValueError("Repair tool/profile lock differs from the proposed deterministic protocol")
    for record in lock["sources"].values():
        path = root / record["archive_path"]
        if _sha(path) != record["sha256"]:
            raise ValueError(f"Repair distribution archive changed: {record['archive_path']}")
    for relative, expected in lock["installed_resource_sha256"].items():
        path = root / relative
        if Path(relative).is_absolute() or ".." in Path(relative).parts or _sha(path) != expected:
            raise ValueError(f"Installed repair tool resource changed: {relative}")
    return lock


def _source_receipt(prepared):
    prepared = Path(prepared)
    receipt = _load(prepared / "source_preparation.json")
    if receipt["case_ids"] != list(CASE_IDS):
        raise ValueError("Repair sensitivity must retain exactly the original fixed six cases")
    for case in receipt["cases"]:
        directory = prepared / case["case_id"]
        for relative, expected in case["artifact_sha256"].items():
            path = (directory / relative).resolve()
            if not path.is_relative_to(directory.resolve()) or _sha(path) != expected:
                raise ValueError(f"Strict primary source artifact changed: {case['case_id']}/{relative}")
    return receipt


def _worker_command(root, lock, source, output, mode):
    return [str(Path(root) / lock["python_executable"]), "-I", str(Path(__file__).resolve()), "worker",
            "--source", str(Path(source).resolve()), "--output", str(Path(output).resolve()),
            "--pdbfixer-source", str((Path(root) / lock["pdbfixer_source_directory"]).resolve()),
            "--mode", mode]


def inspect_sources(root, prepared, output, lock_path):
    """Find source-only missing sidechains, without generating a single atom."""
    root, prepared, output = Path(root), Path(prepared), Path(output)
    lock = verify_tools(root, lock_path); source = _source_receipt(prepared)
    if output.exists(): raise FileExistsError("Preserve previous feasibility inventory")
    output.mkdir(parents=True)
    rows = []
    for case_id in CASE_IDS:
        directory = output / case_id; directory.mkdir()
        command = _worker_command(root, lock, prepared / case_id / "receptor.pdb", directory, "inspect")
        with (directory / "worker.log").open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=120,
                                    env={**os.environ, "OPENMM_CPU_THREADS": "1"}, check=False)
        row = _load(directory / "worker_receipt.json") if (directory / "worker_receipt.json").is_file() else {"status": "failed", "reason": "worker did not write a receipt"}
        rows.append({"case_id": case_id, "worker_exit_code": result.returncode, "source_protein_sha256": _sha(prepared / case_id / "receptor.pdb"), **row})
    receipt = {"schema_version": 1, "purpose": "source-only feasibility; no new atoms or docking outcomes generated",
               "fixed_case_ids": list(CASE_IDS), "primary_source_preparation_sha256": _sha(prepared / "source_preparation.json"),
               "tools_lock_sha256": _sha(lock_path), "source_ledger_sha256": source["source_ledger_sha256"],
               "native_glycan_coordinate_access": False, "generated_atom_count": 0, "cases": rows,
               "code_sha256": _sha(Path(__file__))}
    _write(output / "feasibility.json", receipt)
    return receipt


def _verify_secondary_freeze(root, prepared, lock_path, config_path, freeze_path):
    lock = verify_tools(root, lock_path); _source_receipt(prepared)
    config, freeze = _load(config_path), _load(freeze_path)
    if (freeze.get("configuration_sha256") != _sha(config_path)
            or freeze.get("secondary_profile_frozen_before_new_outcomes") is not True
            or freeze.get("tools_lock_sha256") != _sha(lock_path)
            or freeze.get("primary_source_preparation_sha256") != _sha(Path(prepared) / "source_preparation.json")
            or freeze.get("repair_code_sha256") != _sha(Path(__file__))):
        raise ValueError("Secondary atom generation requires coordinator freeze before new outcomes")
    if config.get("case_ids") != list(CASE_IDS) or config.get("seed") != SEED or config.get("platform") != "Reference":
        raise ValueError("Frozen secondary profile differs from fixed cases/seed/Reference platform")
    return lock


def prepare_repaired_sources(root, prepared, output, lock_path, config_path, freeze_path):
    """Generate a separately frozen sensitivity; never alter strict inputs."""
    root, prepared, output = Path(root), Path(prepared), Path(output)
    lock = _verify_secondary_freeze(root, prepared, lock_path, config_path, freeze_path)
    if output.exists(): raise FileExistsError("Preserve previous modeled sensitivity inputs")
    output.mkdir(parents=True)
    rows = []
    for case_id in CASE_IDS:
        directory = output / case_id; directory.mkdir()
        command = _worker_command(root, lock, prepared / case_id / "receptor.pdb", directory, "model")
        with (directory / "worker.log").open("w") as log:
            execution = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=3600,
                                       env={**os.environ, "OPENMM_CPU_THREADS": "1"}, check=False)
        row = _load(directory / "worker_receipt.json") if (directory / "worker_receipt.json").is_file() else {"status": "failed", "reason": "worker did not write a receipt"}
        row.update(case_id=case_id, worker_exit_code=execution.returncode, fixed_denominator_member=True)
        if row["status"] == "prepared":
            row["receptor_parameterization"] = _parameterize_sensitivity(prepared / case_id, directory)
        rows.append(row)
        receipt = {"schema_version": 1, "profile": "standard_sidechain_completion_sensitivity",
                   "primary_strict_inputs_changed": False, "fixed_case_ids": list(CASE_IDS),
                   "configuration_sha256": _sha(config_path), "freeze_sha256": _sha(freeze_path),
                   "tools_lock_sha256": _sha(lock_path), "cases": rows,
                   "native_glycan_coordinate_access": False, "modeled_atoms_are_atomic_truth": False}
        _write(output / "repair_preparation.json", receipt)
    return receipt


def _parameterize_sensitivity(primary, output):
    """Reuse the strict protein audit with a separate modeled-input manifest."""
    import shutil
    import gemmi
    from .aav import prepare_receptor

    primary, output = Path(primary), Path(output)
    for filename in ("ligand.sdf", "ligand.pdbqt", "chemistry.json", "ligand_atom_map.json"):
        shutil.copy2(primary / filename, output / filename)
    structure = gemmi.read_pdb_string((output / "receptor.pdb").read_text())
    structure.make_mmcif_document().write_file(str(output / "receptor.cif"))
    try:
        preparation = prepare_receptor((output / "receptor.pdb").read_text(), output / "parameterized")
        result = {"status": "prepared", "audit": preparation, "pdbqt_sha256": _sha(output / "parameterized.pdbqt")}
    except (ValueError, RuntimeError, OSError) as error:
        result = {"status": "preparation_failed", "reason": str(error), "error_type": type(error).__name__}
    _write(output / "receptor_parameterization.json", result)
    return result


def _isolated_worker(source, output, pdbfixer_source, mode):
    """Isolated worker: import only pinned OpenMM/NumPy and official source."""
    import io
    import importlib.metadata
    sys.path.insert(0, str(Path(pdbfixer_source).resolve()))
    import openmm
    from openmm import unit
    from pdbfixer import PDBFixer
    import pdbfixer.pdbfixer as implementation

    def forbid_download(*args, **kwargs):
        raise RuntimeError("Only pinned local standard templates are authorized; network template retrieval forbidden")
    implementation.urlopen = forbid_download
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    original_text = Path(source).read_text(); original = parse_pdb_heavy(original_text)
    unsupported = sorted({key[3] for key in original} - STANDARD_RESIDUES)
    metadata = {"source_sha256": _sha(source), "seed": SEED, "platform": "Reference", "mode": mode,
                "versions": {p: importlib.metadata.version(p) for p in ("openmm", "numpy")},
                "pdbfixer_release": "1.12", "pdbfixer_commit": "94cfa4c0ca551cdc5f13320f9a658efd59f2b881",
                "native_glycan_coordinate_access": False, "observed_heavy_atom_count": len(original),
                "missing_loops_inserted": False, "nonstandard_residues_converted": False,
                "terminal_atoms_added": False, "modeled_atoms_are_atomic_truth": False}
    if unsupported:
        result = {**metadata, "status": "unsupported", "reason": "Nonstandard protein retained without mutation or unpinned template/charge assumptions", "nonstandard_residues": unsupported, "new_atom_count": 0}
        _write(output / "worker_receipt.json", result); return result
    segment_text, segment_map = split_modeling_segments(original_text)
    _write(output / "modeling_segment_map.json", segment_map)

    def make_fixer():
        fixer = PDBFixer(pdbfile=io.StringIO(segment_text), platform=openmm.Platform.getPlatformByName("Reference"))
        fixer.findMissingResidues(); fixer.missingResidues = {}
        fixer.findMissingAtoms()
        missing, prohibited = {}, []
        for residue, atoms in fixer.missingAtoms.items():
            allowed_atoms = []
            for atom in atoms:
                key = (segment_map["chains"][residue.chain.id], int(residue.id), residue.insertionCode.strip(), residue.name, atom.name)
                if atom.name in BACKBONE:
                    prohibited.append(list(key))
                else:
                    allowed_atoms.append(atom)
            if allowed_atoms: missing[residue] = allowed_atoms
        suppressed_terminal = [{"chain": segment_map["chains"][r.chain.id], "number": int(r.id), "resname": r.name, "atoms": list(names)}
                               for r, names in fixer.missingTerminals.items()]
        fixer.missingAtoms = missing; fixer.missingTerminals = {}
        expected = [(segment_map["chains"][r.chain.id], int(r.id), r.insertionCode.strip(), r.name, atom.name)
                    for r, atoms in missing.items() for atom in atoms]
        return fixer, expected, prohibited, suppressed_terminal

    fixer, expected, prohibited, terminals = make_fixer()
    metadata.update(allowed_added_atom_keys=[list(k) for k in expected], missing_sidechain_heavy_atom_count=len(expected),
                    prohibited_missing_backbone_atom_keys=prohibited, suppressed_terminal_additions=terminals)
    if mode == "inspect" or prohibited:
        result = {**metadata, "status": "unsupported_missing_backbone" if prohibited else "feasible",
                  "new_atom_count": 0, "atom_generation_performed": False}
        _write(output / "worker_receipt.json", result); return result
    results = []
    # Repeat from identical source with the fixed seed; compare modeled output.
    for repeat in range(2):
        if repeat: fixer, _, _, _ = make_fixer()
        fixer.addMissingAtoms(seed=SEED)
        positions = fixer.positions.value_in_unit(unit.angstrom)
        coordinates = {}
        for atom in fixer.topology.atoms():
            residue = atom.residue
            key = (segment_map["chains"][residue.chain.id], int(residue.id), residue.insertionCode.strip(), residue.name, atom.name)
            if key in coordinates: raise ValueError("PDBFixer output has duplicate restored source keys")
            coordinates[key] = {"coordinates": list(positions[atom.index]), "element": atom.element.symbol}
        audit = audit_repair(original, coordinates, expected)
        results.append((coordinates, audit))
    new_difference = max((math.sqrt(sum((a-b)**2 for a,b in zip(results[0][0][key]["coordinates"], results[1][0][key]["coordinates"]))) for key in expected), default=0.)
    if new_difference > .000001:
        raise ValueError(f"Repeated fixed-seed Reference modeling differs by {new_difference:g} A")
    repaired_text = _serialize_preserving_observed(original, results[0][0], expected)
    serialization_audit = audit_repair(original, parse_pdb_heavy(repaired_text), expected)
    (output / "receptor.pdb").write_text(repaired_text)
    result = {**metadata, "status": "prepared", "new_atom_count": len(expected), "atom_generation_performed": True,
              "pre_serialization_audit": results[0][1], "serialization_audit": serialization_audit,
              "determinism": {"repetitions": 2, "max_new_atom_displacement_A": new_difference, "tolerance_A": .000001},
              "receptor_sha256": _sha(output / "receptor.pdb")}
    _write(output / "worker_receipt.json", result)
    return result


def _serialize_preserving_observed(original, coordinates, expected):
    residues = defaultdict(list)
    for key, row in original.items(): residues[key[:4]].append((key, row))
    added = defaultdict(list)
    for key in expected: added[key[:4]].append(key)
    rows, serial, previous_chain = [], 0, None
    for residue_key, entries in residues.items():
        chain, number, insertion, resname = residue_key
        if previous_chain is not None and previous_chain != chain: rows.append("TER")
        for _, row in entries:
            serial += 1; rows.append(row["line"][:6] + f"{serial:5d}" + row["line"][11:])
        b_iso = sum(float(row["line"][60:66]) for _, row in entries) / len(entries)
        for key in sorted(added[residue_key]):
            serial += 1; atom = coordinates[key]; name = key[-1]; xyz = atom["coordinates"]
            rows.append(f"ATOM  {serial:5d} {name:>4s} {resname:3s} {chain:1s}{number:4d}{insertion or ' ':1s}   "
                        + "".join(f"{x:8.3f}" for x in xyz) + f"{1.:6.2f}{b_iso:6.2f}          {atom['element']:>2s}")
        previous_chain = chain
    rows.extend(("TER", "END")); return "\n".join(rows) + "\n"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("worker", choices=("worker",))
    parser.add_argument("--source", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--pdbfixer-source", required=True); parser.add_argument("--mode", choices=("inspect", "model"), required=True)
    args = parser.parse_args()
    try:
        _isolated_worker(args.source, args.output, args.pdbfixer_source, args.mode)
    except Exception as error:
        _write(Path(args.output) / "worker_receipt.json", {"status": "failed", "error_type": type(error).__name__, "reason": str(error)})
        raise
