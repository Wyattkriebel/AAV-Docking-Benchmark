"""Offline replay of archived site inputs, surface nulls and native-site metrics.

Predictor executables are not installed or executed. The raw-prediction verifier
reparses their archived outputs. The surface null is resampled from protein-only
inputs, and the separate evaluator reads native ligand coordinates afterwards.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from .provenance import sha256
from .site_experiment import evaluate_experiment
from .site_geometry import capsid_symmetry, read_receptor
from .site_prediction_verification import restore_site_inputs, verify_site_predictions
from .surface import build_surface_geometry, sample_surface_patches


_METRIC_FLOAT_FIELDS = frozenset({
    "CA_rmsd_A", "rotation", "translation", "evaluated_center_A",
    "center_to_native_atom_A", "top1_distance_A", "top5_min_distance_A",
    "best_returned_distance_A", "recall", "precision",
    "top1_common_interface_recall", "top1_method_interface_recall",
    "top5_common_union_interface_recall", "top5_method_union_interface_recall",
})
_SURFACE_FLOAT_FIELDS = frozenset({
    "center_A", "footprint_max_atom_distance_A",
})
_SURFACE_AREA_FIELD = "estimated_probe_exposed_surface_area_A2"
_SYMMETRY_FLOAT_FIELDS = frozenset({
    "rotation", "translation", "fit_CA_rmsd_A", "chain_proposal_whole_capsid_max_A",
    "whole_capsid_max_CA_displacement_A", "whole_capsid_refinement_CA_rmsd_A",
})


def _validate_tolerances(numeric_tolerance, surface_area_tolerance):
    if (isinstance(numeric_tolerance, bool) or not math.isfinite(numeric_tolerance)
            or not 0 <= numeric_tolerance <= 1e-8):
        raise ValueError("Numeric tolerance must be finite and between 0 and 1e-8")
    if (isinstance(surface_area_tolerance, bool) or not math.isfinite(surface_area_tolerance)
            or not 0 <= surface_area_tolerance <= 1e-6):
        raise ValueError("Surface area tolerance must be finite and between 0 and 1e-6 square Angstroms")


def compare_site_values(expected, actual, *, numeric_tolerance=0.0, surface_area_tolerance=0.0,
                        derived_fields=_METRIC_FLOAT_FIELDS):
    """Compare JSON values, limiting tolerance to explicitly named derived fields.

    Shapes, types, discrete decisions, ranks, counts, residue/atom identities,
    copied coordinates, scores and provenance always require exact equality.
    Matrix/list entries inherit the field name. Nonfinite values are rejected.
    Surface area has a separate explicit tolerance in square Angstroms, applied
    only to the derived accumulated area field, never to distances or centers.
    """
    _validate_tolerances(numeric_tolerance, surface_area_tolerance)
    maximum = 0.0
    differences_by_field = {}

    def compare(left, right, field=None, path="root"):
        nonlocal maximum
        if type(left) is not type(right):
            raise ValueError(f"Site replay type differs at {path}")
        if isinstance(left, float):
            if not math.isfinite(left) or not math.isfinite(right):
                raise ValueError(f"Nonfinite site replay value at {path}")
            difference = abs(left - right)
            permitted = (surface_area_tolerance if field == _SURFACE_AREA_FIELD
                         else numeric_tolerance if field in derived_fields else 0.0)
            if difference > permitted:
                raise ValueError(f"Site replay floating value differs at {path}: {difference:g}")
            maximum = max(maximum, difference)
            field_name = field if field is not None else "<root>"
            differences_by_field[field_name] = max(differences_by_field.get(field_name, 0.0), difference)
        elif isinstance(left, dict):
            if left.keys() != right.keys():
                raise ValueError(f"Site replay fields differ at {path}")
            for key in left:
                compare(left[key], right[key], key, f"{path}.{key}")
        elif isinstance(left, list):
            if len(left) != len(right):
                raise ValueError(f"Site replay row count differs at {path}")
            for index, (a, b) in enumerate(zip(left, right)):
                compare(a, b, field, f"{path}[{index}]")
        elif left != right:
            raise ValueError(f"Site replay identity, count, decision or copied value differs at {path}")

    compare(expected, actual)
    return {"exact_recomputation": maximum == 0.0, "max_numeric_difference": maximum,
            "max_numeric_difference_by_field": dict(sorted(differences_by_field.items()))}


def _field_maxima(checks):
    maxima = {}
    for check in checks:
        for field, difference in check["max_numeric_difference_by_field"].items():
            maxima[field] = max(maxima.get(field, 0.0), difference)
    return dict(sorted(maxima.items()))


def _safe_relative(value):
    if not isinstance(value, str) or "\\" in value:
        raise ValueError("Site replay source path must be portable")
    path = Path(value)
    if path.is_absolute() or not path.parts or any(part in {".", ".."} for part in path.parts):
        raise ValueError("Site replay source path must stay inside the repository")
    return path


def _check_sources(root, sources):
    if not isinstance(sources, dict) or not sources:
        raise ValueError("Site evaluation has no source checksum inventory")
    for relative, digest in sources.items():
        path = root / _safe_relative(relative)
        if path.is_symlink() or not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Site evaluation source checksum mismatch: {relative}")


def restore_truth_sources(root, destination):
    """Restore only required native truth bytes from repository gzip archives.

    No existing ``data/raw`` tree and no network fallback are consulted. A missing
    archive or checksum mismatch therefore fails a clean-checkout verification.
    """
    root, destination = Path(root), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Native truth restoration requires a fresh destination")
    configuration = json.loads((root / "benchmarks/site_experiment.json").read_text())
    source_lock = json.loads((root / "benchmarks/sources.lock.json").read_text())
    external = json.loads((root / "benchmarks/external_cohort.json").read_text())
    holo = [source for source in source_lock["sources"] if source["filename"] == "7MTZ.cif"]
    if len(holo) != 1:
        raise ValueError("Native AAV truth source is missing or ambiguous")
    inventory = [(Path("7MTZ.cif"), holo[0])]
    for case in configuration["cases"]:
        if case["complex_id"] == "aav9":
            continue
        matches = [entry for entry in external["cases"] if entry["case_id"] == case["complex_id"]]
        if len(matches) != 1:
            raise ValueError("Native external truth source is missing or ambiguous")
        for suffix in ("r", "l"):
            filename = f"{case['complex_id']}_{suffix}_b.pdb"
            sources = [s for s in matches[0]["sources"] if s["filename"] == filename]
            if len(sources) != 1:
                raise ValueError("Native external truth file is missing or ambiguous")
            inventory.append((Path("external") / filename, sources[0]))
    destination.mkdir(parents=True)
    receipt = []
    for relative, source in inventory:
        archive = root / _safe_relative(source["archive"])
        if archive.is_symlink() or not archive.is_file():
            raise ValueError(f"Missing offline native truth archive: {source['archive']}")
        digest = hashlib.sha256()
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(archive, "rb") as compressed, target.open("xb") as output:
            for chunk in iter(lambda: compressed.read(1024 * 1024), b""):
                digest.update(chunk)
                output.write(chunk)
        if digest.hexdigest() != source["sha256"]:
            raise ValueError(f"Native truth archive checksum mismatch: {relative}")
        receipt.append({"filename": relative.as_posix(), "sha256": digest.hexdigest(),
                        "archive": source["archive"], "archive_sha256": sha256(archive)})
    return receipt


def _verify_surface(root, name, receptor, configuration, tolerance, surface_area_tolerance=0.0):
    directory = root / "results/site_discovery" / name
    record_path = directory / "random_surface_v2/result.json"
    record = json.loads(record_path.read_text())
    baseline = configuration["baseline"]
    keys, coordinates, elements = read_receptor(receptor)
    if (record.get("status") != "ok" or record.get("method") != "random_surface"
            or record.get("condition_id") != name or record.get("input_sha256") != sha256(receptor)
            or record.get("configuration_sha256") != sha256(root / "benchmarks/site_experiment.json")
            or record.get("points_per_atom") != baseline["points_per_atom"]
            or record.get("probe_radius_A") != baseline["probe_radius_A"]
            or record.get("native_information_used") is not False
            or record.get("surface_code_sha256") != sha256(Path(__file__).with_name("surface.py"))):
        raise ValueError("Surface baseline input, status, configuration or implementation differs")
    trials = record["replicates"]
    if (len(trials) != baseline["replicates"]
            or [trial["parameters"]["seed"] for trial in trials] != baseline["seeds"]):
        raise ValueError("Surface baseline seed/replicate denominator differs")
    transforms, symmetry_check = [], None
    symmetry_path = directory / "random_surface_v2/source_symmetry.json"
    if name.startswith("aav9_"):
        symmetry = json.loads(symmetry_path.read_text())
        amendment = json.loads((root / "benchmarks/site_geometry_amendment.json").read_text())
        _, recomputed = capsid_symmetry(
            keys, coordinates, tolerance_A=amendment["effective_operator_max_CA_displacement_A"],
            chain_fit_tolerance_A=amendment["chain_fit_max_CA_rmsd_A"])
        # The holo null was generated under the original stricter guard; both
        # guards accept the identical source-derived operators, not new fits.
        comparable = {key: symmetry[key] for key in recomputed if key not in {"tolerance_A", "chain_fit_tolerance_A"}}
        symmetry_check = compare_site_values(
            comparable, {key: recomputed[key] for key in comparable}, numeric_tolerance=tolerance,
            derived_fields=_SYMMETRY_FLOAT_FIELDS)
        transforms = [(np.asarray(t["rotation"]), np.asarray(t["translation"])) for t in symmetry["transforms"]]
        if (symmetry["operator_count"] != 60 or record.get("symmetry_operator_count") != 60
                or symmetry["tolerance_A"] > amendment["effective_operator_max_CA_displacement_A"]):
            raise ValueError("Surface baseline source-symmetry count or guard differs")
    elif record.get("symmetry_operator_count") != 1:
        raise ValueError("Noncapsid surface baseline unexpectedly uses source symmetry")
    geometry = build_surface_geometry(
        coordinates, elements, [tuple(key[:4]) for key in keys],
        points_per_atom=baseline["points_per_atom"], probe_radius_A=baseline["probe_radius_A"])
    if record.get("exposed_surface_points") != geometry.exposed_point_count:
        raise ValueError("Surface baseline exposed-point denominator differs")
    area_check = compare_site_values(
        {"estimated_probe_exposed_surface_area_A2": record["estimated_surface_area_A2"]},
        {"estimated_probe_exposed_surface_area_A2": geometry.estimated_surface_area_A2},
        numeric_tolerance=tolerance, surface_area_tolerance=surface_area_tolerance,
        derived_fields=_SURFACE_FLOAT_FIELDS)
    checks = []
    for trial in trials:
        recomputed = sample_surface_patches(
            geometry, seed=trial["parameters"]["seed"], max_sites=max(configuration["candidate_budgets"]),
            patch_radius_A=baseline["patch_radius_A"], min_center_separation_A=baseline["minimum_center_separation_A"],
            symmetry_transforms=transforms)
        # JSON serialization converts tuple residue IDs to lists in the archive.
        recomputed = json.loads(json.dumps(recomputed, allow_nan=False))
        checks.append(compare_site_values(trial, recomputed, numeric_tolerance=tolerance,
                                          surface_area_tolerance=surface_area_tolerance,
                                          derived_fields=_SURFACE_FLOAT_FIELDS))
    return {"condition_id": name, "replicates_resampled": len(checks),
            "record_sha256": sha256(record_path), "source_symmetry_recomputed": symmetry_check,
            "surface_area_recomputation": area_check,
            "all_replicates_exact": all(check["exact_recomputation"] for check in checks),
            "max_numeric_difference": max(check["max_numeric_difference"] for check in [area_check, *checks]),
            "max_numeric_difference_by_field": _field_maxima([area_check, *checks])}


def verify_site_experiment(root, *, numeric_tolerance=0.0, surface_area_tolerance=0.0):
    """Verify complete native-tool-free replay in temporary working directories."""
    _validate_tolerances(numeric_tolerance, surface_area_tolerance)
    root = Path(root).resolve()
    original = root / "results/site_discovery"
    expected_path = original / "evaluation.json"
    expected = json.loads(expected_path.read_text())
    _check_sources(root, expected["source_sha256"])
    if (expected["evaluator_sha256"] != sha256(Path(__file__).with_name("site_experiment.py"))
            or expected["geometry_code_sha256"] != sha256(Path(__file__).with_name("site_geometry.py"))):
        raise ValueError("Site evaluator/geometry implementation differs from the deposited result")
    configuration = json.loads((root / "benchmarks/site_experiment.json").read_text())
    prediction_check = verify_site_predictions(root)
    surface_checks, table_checks = [], []
    with tempfile.TemporaryDirectory(prefix="aav-site-replay-") as directory:
        temporary = Path(directory)
        inputs, outputs, raw = temporary / "inputs", temporary / "outputs", temporary / "raw"
        truth_sources = restore_truth_sources(root, raw)
        for case in configuration["cases"]:
            name = case["condition_id"]
            restore_site_inputs(root, name, "original", inputs / name)
            surface_checks.append(_verify_surface(root, name, inputs / name / "receptor.cif",
                                                   configuration, numeric_tolerance, surface_area_tolerance))
        actual = evaluate_experiment(root, input_root=inputs, output_root=outputs, raw_root=raw)
        generated = set()
        for case in configuration["cases"]:
            name = case["condition_id"]
            reference_relative = f"results/site_discovery/{name}/evaluation_reference.json"
            generated.add(reference_relative)
            reference_check = compare_site_values(
                json.loads((original / name / "evaluation_reference.json").read_text()),
                json.loads((outputs / name / "evaluation_reference.json").read_text()),
                numeric_tolerance=numeric_tolerance)
            table_checks.append({"file": reference_relative, **reference_check})
            relative = f"results/site_discovery/{name}/evaluation_details.json"
            generated.add(relative)
            left = json.loads((original / name / "evaluation_details.json").read_text())
            right = json.loads((outputs / name / "evaluation_details.json").read_text())
            if (left["reference_sha256"] != sha256(original / name / "evaluation_reference.json")
                    or right["reference_sha256"] != sha256(outputs / name / "evaluation_reference.json")):
                raise ValueError("Site details reference digest does not match its own reference artifact")
            # Both reference artifacts have just been compared. Hashes may
            # differ only because an allowed derived float serialized differently.
            right["reference_sha256"] = left["reference_sha256"]
            table_checks.append({"file": relative, **compare_site_values(
                left, right, numeric_tolerance=numeric_tolerance)})
        if actual["source_sha256"].keys() != expected["source_sha256"].keys():
            raise ValueError("Site replay source inventory differs")
        for relative, digest in actual["source_sha256"].items():
            path = outputs / Path(relative).relative_to("results/site_discovery") if relative in generated else root / relative
            if sha256(path) != digest:
                raise ValueError("Site replay source digest differs from its own artifact")
            if relative not in generated and digest != expected["source_sha256"][relative]:
                raise ValueError("Site replay changed a copied input/source digest")
        actual["source_sha256"] = expected["source_sha256"]
        table_checks.append({"file": "results/site_discovery/evaluation.json", **compare_site_values(
            expected, actual, numeric_tolerance=numeric_tolerance)})
    return {"schema_version": 2, "status": "verified", "scientific_predictors_executed": False,
            "network_used": False, "existing_raw_or_prepared_inputs_used": False,
            "numeric_tolerance": numeric_tolerance,
            "surface_area_tolerance_A2": surface_area_tolerance,
            "surface_area_tolerance_field": _SURFACE_AREA_FIELD,
            "comparison_policy": "Explicit derived geometry floats only; a separate square-Angstrom tolerance applies exclusively to accumulated estimated_probe_exposed_surface_area_A2. Copied coordinates, scores, provenance, ranks, counts, decisions and shapes are exact. Original and replay artifact hashes checked independently before semantic comparison.",
            "eligible_conditions": len(configuration["cases"]), "native_truth_sources": truth_sources,
            "predictor_archives": {key: prediction_check[key] for key in (
                "verified_corrected_method_conditions", "verified_initial_CIF_reader_failures",
                "corrected_completed_conditions", "corrected_workflow_failures")},
            "surface_nulls": surface_checks, "tables": table_checks,
            "surface_max_numeric_difference_by_field": _field_maxima(surface_checks),
            "table_max_numeric_difference_by_field": _field_maxima(table_checks),
            "all_tables_exact": all(check["exact_recomputation"] for check in table_checks),
            "source_evaluation_sha256": sha256(expected_path), "verifier_sha256": sha256(Path(__file__))}
