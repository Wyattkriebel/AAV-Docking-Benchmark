"""Read-only replay of the frozen source-only atomic-reference qualification.

This module reads experimental sources and their provenance only. It never reads
predicted poses, docking scores or recovery outcomes. Numerical tolerance applies
only to normalized density samples; atom identity, coordinates, graph, chirality
signs and eligibility classifications must agree exactly. Manual source-quality
adjudications are hash-pinned inputs, not conclusions proved by this verifier.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

import gemmi
import numpy as np

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEDGER = Path("benchmarks/atomic_reference_ledger.json")
DEFAULT_DENSITY_TOLERANCE_SIGMA = 1e-5
CASE_IDS = ("154L", "1GNY", "1KJL", "1QFO", "2J1T", "1PWB", "2VUZ")
DENSITY_FIELDS = ("two_mFo_DFc_atom_sigma", "two_mFo_DFc_max_within_0p5A_sigma", "mFo_DFc_atom_sigma")


class AtomicReferenceError(ValueError):
    """Frozen provenance or source replay failed."""


def _load(path: Path) -> Any:
    return json.loads(path.read_text())


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _path(root: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def _same(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AtomicReferenceError(f"{label}: source replay differs from frozen reference")


def _check_pin(root: Path, source: dict, checked: set[str]) -> Path:
    path = _path(root, source["archive"])
    if not path.is_file() or _sha(path) != source["sha256"]:
        raise AtomicReferenceError(f"checksum mismatch: {source['archive']}")
    if "bytes" in source and path.stat().st_size != source["bytes"]:
        raise AtomicReferenceError(f"byte count mismatch: {source['archive']}")
    checked.add(source["archive"])
    return path


def _xyz(category: dict, index: int) -> np.ndarray:
    return np.array([float(category["Cartn_" + axis][index]) for axis in "xyz"])


def _read_grid(path: Path, kind: str, sample_rate: float):
    block = gemmi.as_refln_blocks(gemmi.cif.read_file(str(path)))[0]
    amplitude, phase = (("pdbx_FWT", "pdbx_PHWT") if kind == "2fo-fc"
                        else ("pdbx_DELFWT", "pdbx_DELPHWT"))
    grid = block.transform_f_phi_to_map(amplitude, phase, sample_rate=sample_rate)
    values = np.asarray(grid)
    mean, sigma = float(values.mean()), float(values.std())
    if not np.isfinite(sigma) or sigma <= 0:
        raise AtomicReferenceError(f"Invalid full-cell density sigma: {path}")
    values[:] = (values - mean) / sigma
    return grid, dict(shape=list(grid.shape), full_cell_mean_before_normalization=mean,
                      full_cell_sigma=sigma, unit_cell=list(grid.unit_cell.parameters),
                      spacegroup=grid.spacegroup.hm)


def _sample(grid, position: np.ndarray) -> float:
    value = float(grid.interpolate_value(gemmi.Position(*map(float, position))))
    if not np.isfinite(value):
        raise AtomicReferenceError("Nonfinite interpolated source density")
    return value


def _ccd(path: Path):
    block = gemmi.cif.read_file(str(path))[0]
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    ideal = {name: np.array([float(atoms["pdbx_model_Cartn_" + axis + "_ideal"][i])
                            for axis in "xyz"])
             for i, name in enumerate(atoms["atom_id"])}
    edges = [(bonds["atom_id_1"][i], bonds["atom_id_2"][i], bonds["value_order"][i])
             for i in range(len(bonds["atom_id_1"]))
             if not (bonds["atom_id_1"][i].startswith("H") or bonds["atom_id_2"][i].startswith("H"))]
    return atoms, ideal, edges


def _chirality(atoms, ideal, edges, coordinates, deposited, links, category):
    results = []
    for i, name in enumerate(atoms["atom_id"]):
        if atoms["pdbx_stereo_config"][i] not in ["R", "S"] or name not in coordinates:
            continue
        neighbors = sorted([b if a == name else a for a, b, _ in edges if a == name or b == name])
        actual = {key: coordinates[key] for key in neighbors if key in coordinates}
        substitutions = {}
        for key in neighbors:
            if key in actual:
                continue
            candidates = []
            for link in links:
                for this, other in [(link["left"], link["right"]), (link["right"], link["left"])]:
                    if this == [deposited["author_chain"], deposited["author_resnum"], name]:
                        hits = [j for j in range(len(category["id"]))
                                if [category["auth_asym_id"][j], category["auth_seq_id"][j],
                                    category["label_atom_id"][j]] == other]
                        if len(hits) == 1 and category["type_symbol"][hits[0]] == "O":
                            candidates.append((other, _xyz(category, hits[0])))
            if len(candidates) == 1 and key.startswith("O"):
                actual[key] = candidates[0][1]
                substitutions[key] = candidates[0][0]
        if len(neighbors) not in [3, 4] or set(neighbors) - set(actual):
            results.append(dict(atom=name, status="unresolved", neighbors=neighbors))
            continue

        def determinant(positions, center):
            vectors = [positions[key] - center for key in neighbors]
            matrix = vectors if len(vectors) == 3 else [vectors[k] - vectors[3] for k in range(3)]
            return float(np.linalg.det(np.array(matrix)))

        ideal_volume = determinant(ideal, ideal[name])
        observed_volume = determinant(actual, coordinates[name])
        results.append(dict(atom=name, ccd_config=atoms["pdbx_stereo_config"][i],
                            neighbors=neighbors, substitutions=substitutions,
                            ideal_signed_volume=ideal_volume, observed_signed_volume=observed_volume,
                            status="consistent" if ideal_volume * observed_volume > 0 else "inverted"))
    return results


def _residue_failures(validation, atoms, flags, stereo, policy):
    required = policy["required_per_scored_residue"]
    failures = []
    for key, operation, limit in [("rscc", "min", required["RSCC_min"]),
                                   ("rsr", "max", required["RSR_max"]),
                                   ("avgoccu", "min", required["occupancy_min"])]:
        if key not in validation or (validation[key] < limit if operation == "min" else validation[key] > limit):
            failures.append(f"{key}_{operation}_gate")
    supported = sum(a["two_mFo_DFc_max_within_0p5A_sigma"] >= 1 for a in atoms)
    fraction = supported / len(atoms)
    negative = min(a["mFo_DFc_atom_sigma"] for a in atoms)
    if fraction < required["fraction_heavy_atoms_with_local_2mFo_DFc_at_least_1sigma_min"]:
        failures.append("2mFo_DFc_coverage_below_90_percent")
    if negative <= required["atom_center_mFo_DFc_min_sigma_exclusive"]:
        failures.append("negative_difference_density_at_atom")
    if any(flag["kind"] == "chiral-outlier" for flag in flags):
        failures.append("wwPDB_chirality_outlier")
    if any(flag["kind"] in ["mog-bond-outlier", "mog-angle-outlier"]
           and abs(float(flag["Zscore"])) >= required["absolute_mogul_bond_or_angle_Zscore_max_exclusive"]
           for flag in flags):
        failures.append("absolute_bond_or_angle_Z_at_least_5")
    if any(item["status"] != "consistent" for item in stereo):
        failures.append("CCD_signed_volume_stereochemistry_unresolved_or_inverted")
    if any(a["released_coordinate_delta_A"] > .000001 for a in atoms):
        failures.append("released_deposited_coordinate_mismatch")
    return failures, dict(supported_atoms=supported, total_atoms=len(atoms),
                          fraction_supported=fraction, min_atom_difference_sigma=negative)


def replay_case(case: dict, root: Path, policy: dict) -> dict:
    """Recompute all source residues for one pinned case; do not write files.

    The public verifier checks the source hashes before calling this function.
    This lower-level replay is also useful for controlled source perturbations.
    """
    category = gemmi.cif.read_file(str(_path(root, case["sources"]["coordinate"]["archive"])))[0]
    atom_site = category.get_mmcif_category("_atom_site.")
    required = policy["required_per_scored_residue"]
    maps = {kind: _read_grid(_path(root, case["sources"][source_key]["archive"]), kind,
                             required["map_fft_sample_rate"])
            for kind, source_key in [("2fo-fc", "2fofc_coefficients"), ("fo-fc", "fofc_coefficients")]}
    xml = ET.fromstring(gzip.open(_path(root, case["sources"]["validation_xml"]["archive"]), "rt").read())
    entry = xml.find("Entry").attrib
    connection = category.get_mmcif_category("_struct_conn.")
    links = []
    for i, kind in enumerate(connection.get("conn_type_id", [])):
        if kind == "covale":
            links.append(dict(left=[connection["ptnr1_auth_asym_id"][i], connection["ptnr1_auth_seq_id"][i], connection["ptnr1_label_atom_id"][i]],
                              right=[connection["ptnr2_auth_asym_id"][i], connection["ptnr2_auth_seq_id"][i], connection["ptnr2_label_atom_id"][i]],
                              distance_A=float(connection["pdbx_dist_value"][i])))
    released = gemmi.read_pdb_string(gzip.open(_path(root, case["released_ligand_source"]["archive"]), "rt").read())
    released_atoms = {}
    for chain in released[0]:
        for residue in chain:
            for atom in residue:
                if atom.element.name in ["H", "D"]:
                    continue
                key = (chain.name, residue.seqid.num, residue.name, atom.name)
                if key in released_atoms:
                    raise AtomicReferenceError(f"Duplicate released atom {case['case_id']}:{key}")
                released_atoms[key] = atom
    step = required["local_density_sampling_grid_A"]
    radius = required["local_density_search_radius_A"]
    nstep = round(radius / step)
    offsets = [np.array(offset) * step for offset in itertools.product(range(-nstep, nstep + 1), repeat=3)
               if sum(value * value for value in offset) * step * step <= radius * radius]
    residues, used = [], set()
    for frozen in case["residue_qualification"]:
        deposited = frozen["deposited"]
        prepared_chain = next((a["prepared_chain"] for a in case["scored_atom_map"]), "B")
        indices = [i for i in range(len(atom_site["id"]))
                   if atom_site["auth_asym_id"][i] == deposited["author_chain"]
                   and atom_site["auth_seq_id"][i] == deposited["author_resnum"]
                   and atom_site["type_symbol"][i] not in ["H", "D"]]
        coordinates = {atom_site["label_atom_id"][i]: _xyz(atom_site, i) for i in indices}
        if len(coordinates) != len(indices):
            raise AtomicReferenceError(f"Duplicate deposited atom {case['case_id']}:{deposited}")
        validation_record = next(x for x in xml.findall("ModelledSubgroup")
                                 if x.attrib.get("chain") == deposited["author_chain"]
                                 and x.attrib.get("resnum") == deposited["author_resnum"]
                                 and x.attrib.get("resname") == deposited["resname"])
        ccd_atoms, ideal, edges = _ccd(root / "data/archive/external" / f"{deposited['resname']}.cif.gz")
        names = set(coordinates)
        full = {name for i, name in enumerate(ccd_atoms["atom_id"]) if ccd_atoms["type_symbol"][i] not in ["H", "D"]}
        stereo = _chirality(ccd_atoms, ideal, edges, coordinates, deposited, links, atom_site)
        atoms = []
        for i in indices:
            name, position = atom_site["label_atom_id"][i], _xyz(atom_site, i)
            key = (prepared_chain, frozen["prepared_resnum"], frozen["prepared_resname"], name)
            if key not in released_atoms:
                raise AtomicReferenceError(f"Missing released atom {case['case_id']}:{key}")
            used.add(key)
            atom = released_atoms[key]
            density = dict(two_mFo_DFc_atom_sigma=_sample(maps["2fo-fc"][0], position),
                           two_mFo_DFc_max_within_0p5A_sigma=max(_sample(maps["2fo-fc"][0], position + offset) for offset in offsets),
                           mFo_DFc_atom_sigma=_sample(maps["fo-fc"][0], position))
            atoms.append(dict(prepared_chain=prepared_chain, prepared_resnum=frozen["prepared_resnum"],
                              prepared_resname=frozen["prepared_resname"], atom_name=name,
                              element=atom_site["type_symbol"][i], deposited_author_chain=deposited["author_chain"],
                              deposited_author_resnum=deposited["author_resnum"], deposited_label_chain=atom_site["label_asym_id"][i],
                              deposited_label_resnum=atom_site["label_seq_id"][i], deposited_atom_site_id=atom_site["id"][i],
                              xyz_A=position.tolist(), occupancy=float(atom_site["occupancy"][i]),
                              B_iso_A2=float(atom_site["B_iso_or_equiv"][i]),
                              released_coordinate_delta_A=float(np.linalg.norm(position - np.array([atom.pos.x, atom.pos.y, atom.pos.z]))), **density))
        validation = {key: float(validation_record.attrib[key]) for key in ["rscc", "rsr", "avgoccu", "owab", "EDIAm", "OPIA"] if key in validation_record.attrib}
        flags = [dict(kind=x.tag, **x.attrib) for x in validation_record]
        failures, local = _residue_failures(validation, atoms, flags, stereo, policy)
        residues.append(dict(prepared_resnum=frozen["prepared_resnum"], prepared_resname=frozen["prepared_resname"],
                             deposited=deposited, heavy_atom_count=len(atoms), missing_CCD_heavy_atoms=sorted(full - names),
                             CCD_missing_atom_note="Only named glycosidic leaving oxygen is allowed absent, and its external replacement is checked in stereo/links.",
                             validation=validation, validation_flags=flags, stereochemistry=stereo, local_density=local,
                             source_gate_pass=not failures, source_gate_failures=failures, atom_map=atoms,
                             CCD_internal_heavy_bonds=[list(edge) for edge in edges if edge[0] in names and edge[1] in names]))
    _same(used, set(released_atoms), f"{case['case_id']} full released atom inventory")
    return dict(case_id=case["case_id"], maps={key: value[1] for key, value in maps.items()},
                entry_validation={key: entry.get(key) for key in ["PDB-resolution", "PDB-R", "EDS_R", "DCC_R", "PDB-Rfree", "DCC_Rfree", "XMLcreationDate", "PDB-revision-date"]},
                covalent_links=links, residues=residues, sphere_sample_count=len(offsets))


def _compare_replay(actual: Any, expected: Any, label: str, tolerance: float, deltas: list[float]) -> None:
    """Compare normalized map floats only with tolerance; chirality is a sign test."""
    if isinstance(expected, dict):
        for key, value in expected.items():
            if key == "reference_atom_key":
                wanted = f"{actual['prepared_chain']}:{actual['prepared_resnum']}:{actual['prepared_resname']}:{actual['atom_name']}"
                _same(wanted, value, label + ".reference_atom_key")
            elif key in ["ideal_signed_volume", "observed_signed_volume"]:
                _same(float(np.sign(actual[key])), float(np.sign(value)), label + "." + key + " sign")
            elif key in DENSITY_FIELDS or key == "min_atom_difference_sigma":
                if key == "two_mFo_DFc_max_within_0p5A_sigma":
                    _same(actual[key] >= 1, value >= 1, label + ".atom_1sigma_coverage_classification")
                elif key == "mFo_DFc_atom_sigma":
                    _same(actual[key] <= -3, value <= -3, label + ".atom_negative_3sigma_classification")
                difference = abs(actual[key] - value)
                deltas.append(difference)
                if not np.isfinite(difference) or difference > tolerance:
                    raise AtomicReferenceError(f"{label}.{key}: density differs by {difference:g} sigma (tolerance {tolerance:g})")
            elif key in ["full_cell_mean_before_normalization", "full_cell_sigma"]:
                # Raw amplitudes are pinned; normalization scalars are informative.
                # Actual normalized atom samples and every threshold classification are checked.
                if not np.isfinite(actual[key]):
                    raise AtomicReferenceError(f"{label}.{key}: nonfinite normalization")
            elif key not in ["additional_manual_block", "scored"]:
                if key not in actual:
                    raise AtomicReferenceError(f"{label}: missing replay field {key}")
                _compare_replay(actual[key], value, label + "." + key, tolerance, deltas)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AtomicReferenceError(f"{label}: array length mismatch")
        for i, (item, frozen) in enumerate(zip(actual, expected)):
            _compare_replay(item, frozen, f"{label}[{i}]", tolerance, deltas)
    else:
        _same(actual, expected, label)


def _qualified_component(replayed: dict, adjudication: dict):
    passed = {r["prepared_resnum"] for r in replayed["residues"] if r["source_gate_pass"]
              and str(r["prepared_resnum"]) not in adjudication["additional_residue_blocks"]}
    residue_map = {(r["deposited"]["author_chain"], r["deposited"]["author_resnum"]): r["prepared_resnum"] for r in replayed["residues"]}
    edges = [(residue_map.get(tuple(link["left"][:2])), residue_map.get(tuple(link["right"][:2]))) for link in replayed["covalent_links"]]
    counts = {r["prepared_resnum"]: r["heavy_atom_count"] for r in replayed["residues"]}
    pending, components = set(passed), []
    while pending:
        component, queue = {min(pending)}, [min(pending)]
        while queue:
            atom = queue.pop()
            for left, right in edges:
                neighbor = right if left == atom else left if right == atom else None
                if neighbor in passed and neighbor not in component:
                    component.add(neighbor)
                    queue.append(neighbor)
        pending -= component
        components.append(sorted(component))
    components.sort(key=lambda c: (-sum(counts[n] for n in c), tuple(c)))
    selected = components[0] if components else []
    full = len(passed) == len(replayed["residues"]) and adjudication["full_ligand_complete"]
    return selected, "full_ligand" if full else "modeled_core" if selected else "rejected"


def verify_atomic_references(root: Path | str = DEFAULT_ROOT, *, ledger_path: Path | str = DEFAULT_LEDGER,
                             density_tolerance_sigma: float = DEFAULT_DENSITY_TOLERANCE_SIGMA) -> dict:
    """Hash-check and replay all seven frozen references, returning a review report."""
    root = Path(root)
    if not np.isfinite(density_tolerance_sigma) or density_tolerance_sigma < 0:
        raise ValueError("Density tolerance must be finite and nonnegative")
    ledger_file = _path(root, ledger_path)
    ledger = _load(ledger_file)
    _same(tuple(c["case_id"] for c in ledger["cases"]), CASE_IDS, "Frozen seven-case denominator/order")
    _same(ledger["cohort_denominator"], 7, "Cohort denominator")
    _same(ledger["source_only"], True, "Ledger source-only flag")
    _same(ledger["prediction_outcomes_consulted"], False, "Ledger outcome-blinding flag")
    checked: set[str] = set()
    for prefix in ["policy", "manual_adjudication"]:
        _check_pin(root, dict(archive=ledger[prefix + "_archive"], sha256=ledger[prefix + "_sha256"]), checked)
    policy = _load(_path(root, ledger["policy_archive"]))
    manual = _load(_path(root, ledger["manual_adjudication_archive"]))
    _same(policy, ledger["required_policy"], "Ledger policy")
    manifest = _load(_path(root, ledger["source_record"]))
    _same(manifest["ledger_sha256"], _sha(ledger_file), "Frozen ledger checksum")
    _same(manifest["source_only"], True, "Manifest source-only flag")
    _same(manifest["prediction_outcomes_consulted"], False, "Manifest outcome-blinding flag")
    _same(manifest["frozen_policy"], dict(archive=ledger["policy_archive"], sha256=ledger["policy_sha256"]), "Manifest policy pin")
    _same(manifest["manual_adjudications"], dict(archive=ledger["manual_adjudication_archive"], sha256=ledger["manual_adjudication_sha256"]), "Manifest manual adjudication pin")
    pins = {item["archive"]: item for item in manifest["archived_files"]}
    for item in pins.values():
        _check_pin(root, item, checked)
    reports, deltas = [], []
    for case in ledger["cases"]:
        label = case["case_id"]
        for source in case["sources"].values():
            _check_pin(root, source, checked)
        for name in ["released_ligand_source", "released_receptor_source"]:
            source = case[name]
            _check_pin(root, dict(archive=source["archive"], sha256=source["compressed_sha256"]), checked)
            _same(hashlib.sha256(gzip.decompress(_path(root, source["archive"]).read_bytes())).hexdigest(), source["uncompressed_sha256"], label + " " + name + " decompressed checksum")
        diagnostic_path = case["full_atom_diagnostic_archive"]
        required_paths = [diagnostic_path] + [f"data/archive/external/{r['deposited']['resname']}.cif.gz" for r in case["residue_qualification"]]
        if any(path not in pins for path in required_paths):
            raise AtomicReferenceError(f"{label}: missing source diagnostic or CCD checksum in manifest")
        diagnostics = _load(_path(root, diagnostic_path))
        _same(diagnostics["policy"], policy, label + " diagnostic policy")
        _same(tuple(c["case_id"] for c in diagnostics["cases"]), CASE_IDS, label + " diagnostic denominator")
        for flag, required in [("source_only", True), ("prediction_outcomes_consulted", False)]:
            if flag in diagnostics:
                _same(diagnostics[flag], required, label + " diagnostic " + flag)
        frozen = next(c for c in diagnostics["cases"] if c["case_id"] == label)
        actual = replay_case(case, root, policy)
        _same(actual.pop("sphere_sample_count"), diagnostics["local_sphere_sample_count"], label + " sphere sampling inventory")
        _compare_replay(actual, frozen, label + " full diagnostics", density_tolerance_sigma, deltas)
        _compare_replay(actual["residues"], case["residue_qualification"], label + " residue ledger", density_tolerance_sigma, deltas)
        _compare_replay(actual["maps"], case["map_evidence"], label + " map ledger", density_tolerance_sigma, deltas)
        _same(actual["entry_validation"], case["entry_validation"], label + " validation entry ledger")
        _same(actual["covalent_links"], case["covalent_links"], label + " covalent link ledger")
        adjudication = manual["decisions"][label]
        selected, scope = _qualified_component(actual, adjudication)
        _same(case["decision_reason"], adjudication["reason"], label + " manual decision reason")
        for residue in case["residue_qualification"]:
            _same(residue["additional_manual_block"], adjudication["additional_residue_blocks"].get(str(residue["prepared_resnum"])), label + " pinned manual residue block")
            _same(residue["scored"], residue["prepared_resnum"] in selected, label + " per-residue scored decision")
        _same(case["excluded_prepared_residue_numbers"], sorted({r["prepared_resnum"] for r in actual["residues"]} - set(selected)), label + " excluded whole residues")
        _same(selected, case["scored_prepared_residue_numbers"], label + " deterministic selected core")
        _same(scope, case["qualification"], label + " qualification scope")
        _same(bool(selected), case["atomic_pose_eligible"], label + " atomic eligibility")
        _same(scope == "full_ligand", case["full_ligand_eligible"], label + " full-ligand eligibility")
        _same(scope == "modeled_core", case["modeled_core_eligible"], label + " core eligibility")
        expected_atoms = [a for r in actual["residues"] if r["prepared_resnum"] in selected for a in r["atom_map"]]
        _compare_replay(expected_atoms, case["scored_atom_map"], label + " scored atom map", density_tolerance_sigma, deltas)
        _same(len(expected_atoms), case["scored_atom_count"], label + " scored atom count")
        total = sum(r["heavy_atom_count"] for r in actual["residues"])
        _same(total, case["released_ligand_heavy_atom_count"], label + " full released atom count")
        reports.append(dict(case_id=label, qualification=scope, scored_residues=selected,
                            scored_atoms=len(expected_atoms), released_atoms=total,
                            residue_count=len(actual["residues"]), status="verified"))
    _same(ledger["qualification_counts"], {scope: sum(c["qualification"] == scope for c in reports) for scope in ["full_ligand", "modeled_core", "rejected"]}, "Qualification counts")
    return dict(schema_version=1, source_only=True, prediction_outcomes_consulted=False, status="verified",
                ledger_sha256=_sha(ledger_file), policy_sha256=ledger["policy_sha256"],
                manual_adjudication_sha256=ledger["manual_adjudication_sha256"],
                manual_adjudications="Pinned source-quality judgments retained; not independently automated or proved.",
                density_tolerance_sigma=density_tolerance_sigma,
                tolerance_scope="Only normalized Fourier-map samples. Exact atom identities, coordinates, inventories, validation and cutoff decisions. Chirality determinant signs, named neighbors and substitutions are exact; magnitudes are informative.",
                maximum_absolute_density_difference_sigma=max(deltas, default=0), checked_archive_count=len(checked),
                software=dict(gemmi=gemmi.__version__, numpy=np.__version__), cases=reports)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--density-tolerance-sigma", type=float, default=DEFAULT_DENSITY_TOLERANCE_SIGMA)
    args = parser.parse_args(argv)
    try:
        result = verify_atomic_references(args.root, ledger_path=args.ledger,
                                           density_tolerance_sigma=args.density_tolerance_sigma)
    except (ValueError, KeyError, StopIteration, RuntimeError, OSError) as exc:
        print(json.dumps(dict(status="failed", error=str(exc))))
        return 1
    print(json.dumps(result, indent=2))
    return 0
