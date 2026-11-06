"""Frozen site-localization evaluation, separate from predictor execution."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil
from statistics import mean

import numpy as np

from .metrics import fit_receptor
from .provenance import sha256, write_json
from .site_geometry import read_receptor, distinct_candidates, evaluate_candidates
from .site_truth import load_site_truth


def archive_inputs(root):
    root = Path(root)
    config = json.loads((root / "benchmarks/site_experiment.json").read_text())
    for case in config["cases"]:
        name = case["condition_id"]
        output = root / "results/site_discovery" / name / "inputs"
        output.mkdir(parents=True, exist_ok=True)
        entries = []
        for revision, directory in (("original", "site_inputs"), ("fpocket_compatibility", "site_inputs_v2")):
            source = root / "data/processed" / directory / name
            for p in sorted(source.iterdir()):
                if not p.is_file():
                    continue
                filename = revision + "_" + p.name + ".gz"
                with (output / filename).open("wb") as stream:
                    with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
                        with p.open("rb") as original:
                            shutil.copyfileobj(original, zipped)
                entries.append({"revision": revision, "filename": p.name, "archive": filename,
                                "sha256": sha256(p), "archive_sha256": sha256(output / filename)})
        write_json(output / "manifest.json", {"condition_id": name, "files": entries})


def _case_inputs(root, case, input_root=None):
    base = Path(input_root) if input_root is not None else root / "data/processed/site_inputs"
    path = base / case["condition_id"] / "receptor.cif"
    manifest_path = path.with_suffix(".cif.manifest.json")
    manifest = json.loads(manifest_path.read_text())
    if manifest["source_sha256"] != case["source_sha256"] or manifest["receptor_sha256"] != sha256(path):
        raise ValueError("Site receptor/source changed after freezing")
    if manifest.get("native_ligand_access") is not False or manifest.get("distance_based_selection") is not False:
        raise ValueError("Site prediction inputs are not demonstrably reference independent")
    keys, xyz, elements = read_receptor(path)
    return path, keys, xyz


def evaluate_experiment(root, *, input_root=None, raw_root=None, output_root=None,
                        output_path=None, require_complete=True):
    """Recompute endpoints without executing predictors or modifying their records.

    ``output_root`` is the site-result directory for regenerated evaluation files;
    it can be outside the repository for independent, non-destructive replay.
    Archived input and native-source restorations may likewise live elsewhere.
    """
    root = Path(root).resolve()
    output_root = Path(output_root) if output_root is not None else root / "results/site_discovery"
    config_path = root / "benchmarks/site_experiment.json"
    config = json.loads(config_path.read_text())
    receipt = json.loads((root / "results/site_discovery/freeze_receipt.json").read_text())
    if sha256(config_path) != receipt["configuration_sha256"]:
        raise ValueError("Frozen site experiment configuration changed")
    if (receipt.get("parameters_frozen_before_scientific_predictions") is not True
            or receipt["reference_ledger_sha256"] != config["reference_ledger_sha256"]
            or receipt["method_lock_sha256"] != config["method_locks"]
            or receipt["input_source_sha256"] != {c["condition_id"]: c["source_sha256"] for c in config["cases"]}):
        raise ValueError("Site freeze receipt and configuration differ")
    if sha256(root / config["reference_ledger"]) != config["reference_ledger_sha256"]:
        raise ValueError("Site reference eligibility changed after freezing")
    for path, digest in config["method_locks"].items():
        if sha256(root / path) != digest:
            raise ValueError("Frozen site method settings changed")
    amendment = json.loads((root / "benchmarks/site_geometry_amendment.json").read_text())
    amendment_receipt_path = root / "results/site_discovery/geometry_amendment_receipt.json"
    amendment_receipt = json.loads(amendment_receipt_path.read_text())
    if sha256(root / "benchmarks/site_geometry_amendment.json") != amendment_receipt["amendment_sha256"]:
        raise ValueError("Source-only geometry amendment changed")
    if (amendment["original_configuration_sha256"] != receipt["configuration_sha256"]
            or amendment.get("decision_before_any_site_metric_evaluation") is not True
            or amendment_receipt.get("completed_before_site_endpoint_scoring") is not True):
        raise ValueError("Geometry amendment does not bind the frozen source-only decision")
    tolerance = amendment["effective_operator_max_CA_displacement_A"]
    loaded = {c["condition_id"]: _case_inputs(root, c, input_root) for c in config["cases"]}
    truths, rows, sources = {}, [], {}
    def track(path, logical_path=None):
        sources[str(logical_path or Path(path).relative_to(root))] = sha256(path)
    for path in (config_path, root / config["reference_ledger"], root / "benchmarks/site_geometry_amendment.json"):
        track(path)
    track(amendment_receipt_path)
    track(root / "results/site_discovery/freeze_receipt.json")
    for path in config["method_locks"]:
        track(root / path)
    pending = []
    for case in config["cases"]:
        name, complex_id = case["condition_id"], case["complex_id"]
        receptor, keys, xyz = loaded[name]
        if complex_id not in truths:
            truths[complex_id] = load_site_truth(
                complex_id, raw_dir=raw_root or root / "data/raw",
                ledger_path=root / config["reference_ledger"])
            # Provenance records use repository-relative names, independent of
            # the location where the verified native bytes were restored.
            truths[complex_id]["site_reference_ledger"] = config["reference_ledger"]
            for source in truths[complex_id]["source_files"]:
                source["path"] = "data/raw/" + ("external/" if complex_id != "aav9" else "") + Path(source["path"]).name
        truth = truths[complex_id]
        reference_path = output_root / name / "evaluation_reference.json"
        write_json(reference_path, truth)
        track(reference_path, Path("results/site_discovery") / name / reference_path.name)
        operations = []
        if complex_id == "aav9":
            certificate = root / "results/site_discovery" / name / "source_symmetry_amended.json"
            if sha256(certificate) != amendment_receipt["certificates"][name]:
                raise ValueError("Source-only symmetry certificate changed")
            audit = json.loads(certificate.read_text())
            if (audit.get("amendment_sha256") != amendment_receipt["amendment_sha256"]
                    or audit.get("operator_count") != 60 or audit.get("tolerance_A") != tolerance
                    or len(audit["transforms"]) != 60 or len(audit["operators"]) != 60
                    or max(r["whole_capsid_max_CA_displacement_A"] for r in audit["operators"]) > tolerance):
                raise ValueError("Source-only capsid symmetry certification failed")
            operations = [(np.asarray(o["rotation"]), np.asarray(o["translation"])) for o in audit["transforms"]]
            track(certificate)
        rotation, translation = np.eye(3), np.zeros(3)
        alignment = {"scope": "identity: paired released bound coordinate frame"}
        if case["state"] == "apo":
            _, holo_keys, holo_xyz = loaded["aav9_holo"]
            apo_ca = {k: v for k, v in zip(keys, xyz) if k[-1] == "CA"}
            holo_ca = {k: v for k, v in zip(holo_keys, holo_xyz) if k[-1] == "CA"}
            matching = sorted(apo_ca.keys() & holo_ca.keys())
            a, b = np.asarray([apo_ca[k] for k in matching]), np.asarray([holo_ca[k] for k in matching])
            rotation, translation = fit_receptor(a, b)
            alignment = {"scope": "evaluation-only whole-receptor CA alignment; no ligand fit", "matched_CA_count": len(matching),
                         "CA_rmsd_A": float(np.sqrt(np.mean(np.sum((a @ rotation + translation - b) ** 2, axis=1)))),
                         "rotation": rotation.tolist(), "translation": translation.tolist()}
        directory = root / "results/site_discovery" / name
        details = {"condition_id": name, "complex_id": complex_id, "receptor_state": case["state"],
                   "reference_sha256": sha256(reference_path), "alignment": alignment, "methods": {}}
        for method in ("fpocket", "p2rank", "random_surface"):
            result_path = directory / ({"fpocket": "fpocket_v2/execution.json", "p2rank": "p2rank/result.json",
                                        "random_surface": "random_surface_v2/result.json"}[method])
            if not result_path.exists():
                pending.append(f"{name}/{method}")
                continue
            result = json.loads(result_path.read_text())
            track(result_path)
            if result["status"] in ("running", "pending"):
                pending.append(f"{name}/{method}")
                continue
            status = "ok" if result["status"] in ("ok", "completed") else result["status"]
            # The compatibility retry changes only an mmCIF comment delimiter.
            expected_input = sha256(receptor)
            if method == "fpocket":
                if result.get("superseded_input_sha256") != expected_input:
                    raise ValueError("fpocket compatibility input lacks frozen receptor equivalence")
            elif result["input_sha256"] != expected_input:
                raise ValueError("Prediction input differs from the frozen whole receptor")
            trials = result["replicates"] if method == "random_surface" else [result]
            if method == "random_surface":
                if (len(trials) != config["baseline"]["replicates"]
                        or [t["parameters"]["seed"] for t in trials] != config["baseline"]["seeds"]
                        or result.get("configuration_sha256") != receipt["configuration_sha256"]
                        or result.get("condition_id") != name or result.get("native_information_used") is not False
                        or result.get("points_per_atom") != config["baseline"]["points_per_atom"]
                        or result.get("probe_radius_A") != config["baseline"]["probe_radius_A"]):
                    raise ValueError("Surface baseline replicate denominator changed")
            evaluated = []
            sensitivity = {str(distance): [] for distance in config["merge_sensitivity_A"]}
            for trial in trials:
                predictions = trial.get("sites", trial.get("predictions", [])) if status == "ok" else []
                unique, mapping = distinct_candidates(predictions, operations, config["candidate_merge_distance_A"])
                def evaluate(unique, mapping):
                    metrics = evaluate_candidates(unique, truth, keys, xyz, rotation=rotation, translation=translation,
                                                  footprint_radius_A=config["common_footprint_radius_A"], symmetry_transforms=operations,
                                                  symmetry_tolerance_A=tolerance)
                    metrics.update(raw_candidate_count=len(predictions), source_candidate_mapping=mapping,
                                   first_matching_class_4A=next((r["rank"] for r in metrics["candidates"] if r["center_to_native_atom_A"] <= 4), None),
                                   best_returned_distance_A=min((r["center_to_native_atom_A"] for r in metrics["candidates"]), default=None))
                    return metrics
                metrics = evaluate(unique, mapping)
                evaluated.append(metrics)
                for distance in config["merge_sensitivity_A"]:
                    other_unique, other_mapping = distinct_candidates(predictions, operations, distance)
                    # Identical class mappings imply identical representatives,
                    # avoiding repeated large-capsid footprint calculations.
                    other = metrics if other_mapping == mapping else evaluate(other_unique, other_mapping)
                    sensitivity[str(distance)].append(other)
            fields = [f"top{k}_success_{t}A" for k in (1, 5) for t in (3, 4, 5)] + [
                "top1_common_interface_recall", "top1_method_interface_recall", "top5_common_union_interface_recall", "top5_method_union_interface_recall"]
            def aggregate(evaluations):
                return {**{field: mean(float(e[field]) for e in evaluations) for field in fields},
                        "raw_candidate_count": mean(e["raw_candidate_count"] for e in evaluations),
                        "distinct_candidate_count": mean(e["returned_distinct_candidates"] for e in evaluations)}
            row = {"condition_id": name, "complex_id": complex_id, "state": case["state"], "method": method,
                   "status": status, "eligible_conditions": 1, "baseline_replicates": len(evaluated) if method == "random_surface" else None,
                   **aggregate(evaluated),
                   "merge_sensitivity": {distance: aggregate(values) for distance, values in sensitivity.items()},
                   "top1_distance_A": evaluated[0]["top1_distance_A"] if method != "random_surface" else None,
                   "top5_min_distance_A": evaluated[0]["top5_min_distance_A"] if method != "random_surface" else None,
                   "first_matching_class_4A": evaluated[0]["first_matching_class_4A"] if method != "random_surface" else None}
            rows.append(row)
            details["methods"][method] = {"status": status, "source_record_sha256": sha256(result_path), "evaluations": evaluated,
                                          "merge_sensitivity": {distance: {"evaluations": values} for distance, values in sensitivity.items()}}
        details_path = output_root / name / "evaluation_details.json"
        write_json(details_path, details)
        track(details_path, Path("results/site_discovery") / name / details_path.name)
    if require_complete and pending:
        raise ValueError(f"Site experiment still has unfinished eligible work: {pending}")
    if not pending and len(rows) != len(config["cases"]) * 3:
        raise ValueError("Site method-by-condition denominator is incomplete")
    result = {"schema_version": 1, "purpose": "published-site localization and coverage, not atomic pose or binding specificity",
              "complete": not pending, "pending": pending, "eligible_complexes": 8, "eligible_conditions": len(config["cases"]),
              "eligible_method_conditions": len(config["cases"]) * 3, "rows": rows,
              "source_sha256": sources, "evaluator_sha256": sha256(Path(__file__)),
              "geometry_code_sha256": sha256(Path(__file__).with_name("site_geometry.py")),
              "limitations": ["development references; no independent-family or held-out estimate", "100 surface draws are Monte Carlo null replicates, not biological observations",
                              "AAV apo/holo are paired conditions of one capsid system", "fpocket compatibility failures and source-only geometry amendments preserved separately"]}
    write_json(output_path or output_root / "evaluation.json", result)
    return result
