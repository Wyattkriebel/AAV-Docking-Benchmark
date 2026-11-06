"""Recompute the frozen Gal9C assay description from raw, audited docking files.

No structural pose metric is computed: reference.json contains independently
generated input coordinates and supplies atom identities only. Missing or failed
runs remain in the fixed 12-compound by 3-seed denominator; ranking metrics are
ineligible until all 36 runs have independently passed the evidence checks.
"""
from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import mean
import tempfile

from .engine import read_poses
from .provenance import sha256, write_text_atomic
from .reproduction import compare_tables
from .specificity import _load_panel, summarize_specificity_scores


_SEEDS = [104729, 130363, 155921]
_RECEPTOR = "holo_complete"


class _Evidence:
    def __init__(self):
        self.sources = {}

    def read(self, path, key):
        path = Path(path)
        if not path.is_file() and path.with_name(path.name + ".gz").is_file():
            path = path.with_name(path.name + ".gz")
        raw = path.read_bytes()
        content = gzip.decompress(raw) if path.suffix == ".gz" else raw
        self.sources[key] = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
        if path.suffix == ".gz":
            self.sources[key].update({"gzip_sha256": hashlib.sha256(raw).hexdigest(), "gzip_bytes": len(raw)})
        return content

    def json(self, path, key):
        return json.loads(self.read(path, key))

    def digest(self, path, key):
        self.read(path, key)
        return self.sources[key]["sha256"]


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _exists(path):
    return path.is_file() or path.with_name(path.name + ".gz").is_file()


def _file(directory, name):
    relative = Path(name)
    _require(isinstance(name, str) and not relative.is_absolute() and ".." not in relative.parts,
             "Artifact name must stay inside its case directory")
    return directory / relative


def _command(row, case, config, receptor, ligand, pose):
    command = row.get("command")
    _require(isinstance(command, list) and len(command) > 1 and all(isinstance(x, str) for x in command),
             f"Missing actual command: {case}/{row.get('seed')}")
    _require(Path(command[0]).name == "vina-1.2.7", "Unexpected Vina executable command")
    _require(len(command[1:]) % 2 == 0, "Malformed command option/value pairs")
    options = dict(zip(command[1::2], command[2::2]))
    _require(len(options) * 2 == len(command) - 1, "Duplicate actual command option")
    fixed = {"--scoring": "vina", "--cpu": "2", "--seed": str(row["seed"]),
             "--exhaustiveness": "16", "--num_modes": "20",
             "--min_rmsd": str(config["min_rmsd_A"]),
             "--energy_range": str(config["energy_range_kcal_mol"])}
    for axis, center, size in zip("xyz", config["box_center"], config["box_size"]):
        fixed["--center_" + axis], fixed["--size_" + axis] = str(center), str(size)
    _require(set(options) == set(fixed) | {"--receptor", "--ligand", "--out"},
             "Unexpected or missing actual docking command option")
    for key, value in fixed.items():
        _require(options[key] == value, f"Actual command differs from frozen {key}")
    for flag, filename in (("--receptor", receptor), ("--ligand", ligand), ("--out", pose)):
        # Archived files may move. Preserve the actual command, certify its case
        # and filename, then verify relocated bytes through the run fingerprints.
        _require(Path(options[flag]).parts[-2:] == (case, filename),
                 f"Actual command {flag} has the wrong case or artifact")


def _log(content, row, config, poses):
    text = content.decode()
    for needle in ("AutoDock Vina v1.2.7", "Scoring function : vina", "Exhaustiveness: 16", "CPU: 2",
                   f"Performing docking (random seed: {row['seed']})"):
        _require(needle in text, f"Raw log does not certify {needle}")
    for label, values in (("Grid center", config["box_center"]), ("Grid size", config["box_size"])):
        found = re.search(label + r"\s*:\s*X\s+(\S+)\s+Y\s+(\S+)\s+Z\s+(\S+)", text)
        _require(found is not None, f"Raw log is missing {label}")
        _require(all(math.isclose(float(actual), expected, abs_tol=0.0001, rel_tol=0)
                     for actual, expected in zip(found.groups(), values)), f"Raw log changed {label}")
    energies = [(int(rank), float(score)) for rank, score in re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+)\s+\S+\s+\S+\s*$", text, re.MULTILINE)]
    _require(energies == [(p["rank"], p["vina_score"]) for p in poses],
             "Raw log score table and raw pose ranks/scores differ")


def _prepared(evidence, prepared_root, row, panel, config, config_hash, audit_hash):
    case = row["case_id"]
    directory = prepared_root / case
    if not _exists(directory / "preparation.json") and (directory / "inputs/manifest.json").is_file():
        directory = directory / "inputs"
        manifest = evidence.json(directory / "manifest.json", f"prepared/{case}/archive_manifest.json")
        seen = set()
        for entry in manifest["files"]:
            name, archive = entry["filename"], entry["archive"]
            _require(Path(name).name == name and Path(archive).name == archive and name not in seen,
                     "Invalid or duplicate preparation archive filename")
            seen.add(name)
            compressed = directory / archive
            _require(sha256(compressed) == entry["archive_sha256"], "Preparation archive bytes changed")
            content = gzip.decompress(compressed.read_bytes())
            _require(hashlib.sha256(content).hexdigest() == entry["sha256"]
                     and len(content) == entry["original_bytes"], "Preparation archive contents changed")
    def read(name):
        return evidence.json(directory / name, f"prepared/{case}/{name}")
    def digest(name):
        return evidence.digest(_file(directory, name), f"prepared/{case}/{name}")
    preparation = read("preparation.json")
    reference = read("reference.json")
    graph = read("source_graph.json")
    receptor_preparation = read("receptor_preparation.json")
    _require(preparation.get("case_id") == case and preparation.get("ligand_id") == row["ligand_id"]
             and preparation.get("assay_row") == row, f"Prepared assay identity differs: {case}")
    _require(preparation.get("status") == "prepared" and preparation.get("pose_reference_eligible") is False
             and reference.get("pose_reference_eligible") is False,
             f"Generated input cannot become pose ground truth: {case}")
    _require(preparation.get("frozen_analysis") == panel["frozen_analysis"], f"Frozen analysis changed: {case}")
    _require(preparation.get("box_center") == config["box_center"] and preparation.get("box_size") == config["box_size"],
             f"Prepared box differs from frozen experiment: {case}")
    panel_hash = config["panel_sha256"]
    fingerprint = preparation["input_fingerprint"]
    _require(fingerprint.get("panel_sha256") == panel_hash and graph.get("panel_input_sha256") == panel_hash,
             f"Prepared source panel fingerprint differs: {case}")
    _require(graph.get("source_row") == row and graph.get("source") == panel["source"]
             and graph.get("atom_map") == reference.get("atom_map"), f"Prepared source graph/map differs: {case}")
    names, atom_map = reference.get("names"), reference.get("atom_map")
    _require(isinstance(names, list) and len(names) == row["expected_heavy_atoms"] and len(set(names)) == len(names)
             and isinstance(atom_map, list) and [a.get("name") for a in atom_map] == names,
             f"Incomplete or duplicate prepared heavy atom map: {case}")
    _require([a.get("source_smiles_atom_index") for a in atom_map] == list(range(len(names)))
             and [a.get("rdkit_index") for a in atom_map] == list(range(len(names))),
             f"Prepared atom map is not a complete source graph map: {case}")
    chemistry = preparation["chemistry"]
    _require(chemistry.get("heavy_atoms") == len(names) and chemistry.get("formal_charge") == row["formal_charge"]
             and chemistry.get("source_isomeric_smiles") == row["isomeric_smiles"]
             and chemistry.get("native_coordinates_used_for_embedding") is False
             and chemistry.get("native_pose_reference_available") is False
             and chemistry.get("stereochemistry_verified_from_3d") is True
             and chemistry.get("embedding_seed") == config["embedding_seed"],
             f"Prepared chemistry/source-conformer provenance differs: {case}")
    _require(digest("receptor_preparation.json") == fingerprint["receptor_preparation_sha256"],
             f"Receptor preparation hash differs: {case}")
    provenance = receptor_preparation["source_provenance"]
    _require(provenance.get("source_cif_sha256") == config["source_receptor_sha256"]
             and provenance.get("audit_sha256") == audit_hash, f"Source receptor provenance differs: {case}")
    _require(len(preparation.get("receptors", [])) == 1
             and preparation["receptors"] == receptor_preparation.get("receptors"), f"Receptor set differs: {case}")
    receptor = preparation["receptors"][0]
    _require(receptor.get("id") == _RECEPTOR, f"Unexpected receptor condition: {case}")
    _require(receptor_preparation.get("source_receptor_sha256") ==
             fingerprint["receptor_files_sha256"].get(receptor["pdb_file"]),
             f"Clean receptor PDB source hash differs: {case}")
    for filename, expected in fingerprint["receptor_files_sha256"].items():
        _require(digest(filename) == expected, f"Prepared receptor artifact hash differs: {case}/{filename}")
    ligand_file = preparation["ligand_file"]
    ligand_hash = digest(ligand_file)
    _require(ligand_hash == preparation["ligand_sha256"]
             and digest("ligand.sdf") == preparation["ligand_sdf_sha256"], f"Prepared ligand hash differs: {case}")
    ligand_text = evidence.read(_file(directory, ligand_file), f"prepared/{case}/{ligand_file}").decode()
    ligand_names = [line[12:16].strip() for line in ligand_text.splitlines()
                    if line.startswith(("ATOM  ", "HETATM")) and line[77:].strip() not in {"H", "HD", "HS"}]
    _require(Counter(ligand_names) == Counter(names), f"Prepared ligand PDBQT/map mismatch: {case}")
    expected = {"config_sha256": config_hash, "ligand_sha256": ligand_hash,
                "preparation_sha256": digest("preparation.json"), "reference_sha256": digest("reference.json"),
                "receptor_sha256": digest(receptor["pdbqt_file"]),
                "receptor_pdb_sha256": digest(receptor["pdb_file"]), "full_assembly_sha256": None}
    return names, receptor, ligand_file, expected


def generate_specificity_report(root, *, prepared_root=None, results_root=None):
    """Return deterministic ``(summary, markdown)`` from checked raw evidence.

    Prepared and result roots may contain ordinary files or per-file ``.gz``
    archives. Explicit paths support archival relocation without changing the
    recorded historical commands. Missing/failed runs produce an ineligible
    report; unexpected runs, duplicate identities, tampered hashes, invalid atom
    maps, or commands outside the frozen protocol raise ``ValueError``.
    """
    root = Path(root).resolve()
    results_root = Path(results_root).resolve() if results_root else root / "results/gal9_specificity"
    prepared_root = Path(prepared_root).resolve() if prepared_root else (
        results_root if (results_root / "gal9c_1/inputs/manifest.json").is_file() else root / "data/processed/gal9_panel")
    evidence = _Evidence()
    def read(relative):
        return evidence.json(root / relative, relative)
    config = read("benchmarks/gal9_experiment.json")
    config_hash = evidence.sources["benchmarks/gal9_experiment.json"]["sha256"]
    panel_path = root / "benchmarks/specificity_resolution.json"
    read("benchmarks/specificity_resolution.json")
    panel_hash = evidence.sources["benchmarks/specificity_resolution.json"]["sha256"]
    panel, rows, pairs, _ = _load_panel(panel_path)
    audit = read("benchmarks/gal9_receptor_audit.json")
    audit_hash = evidence.sources["benchmarks/gal9_receptor_audit.json"]["sha256"]
    freeze = evidence.json(results_root / "freeze_receipt.json", "results/freeze_receipt.json")
    _require(config.get("frozen_before_new_docking") is True and config.get("seeds") == _SEEDS
             and config.get("cases") == [r["case_id"] for r in rows] and len(rows) == 12,
             "Frozen experiment must schedule exactly the 12 source-qualified compounds by three unique seeds")
    _require(config.get("cpu") == 2 and config.get("exhaustiveness") == 16 and config.get("num_modes") == 20
             and config.get("classifier_threshold") is None, "Frozen protocol changed or classifier threshold introduced")
    _require(config.get("panel_sha256") == panel_hash and config.get("receptor_audit_sha256") == audit_hash,
             "Frozen experiment panel/audit hashes differ")
    for key, expected in (("experiment_sha256", config_hash), ("panel_sha256", panel_hash),
                          ("receptor_audit_sha256", audit_hash), ("source_receptor_sha256", config["source_receptor_sha256"])):
        _require(freeze.get(key) == expected, f"Pre-run freeze receipt changed: {key}")
    _require(freeze.get("frozen_before_new_docking") is True, "Missing pre-run freeze declaration")
    source_hash = evidence.digest(root / "data/archive/specificity/3NV3.cif", "source/3NV3.cif")
    _require(source_hash == config["source_receptor_sha256"]
             and source_hash == audit["sources"]["checksums"]["3NV3.cif"]["sha256"], "Source receptor bytes changed")
    _require(audit["decision"].get("supports_known_site_exploratory_CRD_docking") is True,
             "Receptor audit does not support this experiment")
    lock = read("benchmarks/vina.lock.json")
    executable_hashes = {entry["sha256"] for entry in lock["assets"].values()}
    unexpected = [p.name for p in results_root.iterdir() if p.is_dir() and _exists(p / "runs.json")
                  and p.name not in config["cases"]]
    _require(not unexpected, f"Unexpected specificity case directories: {unexpected}")
    run_rows, compound_rows = [], []
    shared_receptor = None
    shared_executable = None
    for row in rows:
        case = row["case_id"]
        names, receptor, ligand_file, expected = _prepared(
            evidence, prepared_root, row, panel, config, config_hash, audit_hash)
        shared = (expected["receptor_sha256"], expected["receptor_pdb_sha256"])
        _require(shared_receptor in (None, shared), "Ligands did not share identical receptor coordinates")
        shared_receptor = shared
        directory = results_root / case
        records = []
        if _exists(directory / "runs.json"):
            stored_config = evidence.json(directory / "configuration.json", f"results/{case}/configuration.json")
            _require(stored_config == config
                     and evidence.sources[f"results/{case}/configuration.json"]["sha256"] == config_hash,
                     f"Case frozen configuration differs: {case}")
            records = evidence.json(directory / "runs.json", f"results/{case}/runs.json")
        _require(isinstance(records, list), f"Run manifest must be a list: {case}")
        identities = [(r.get("receptor"), r.get("seed")) for r in records]
        _require(len(identities) == len(set(identities)), f"Duplicate run identities: {case}")
        _require(set(identities) <= {(_RECEPTOR, seed) for seed in _SEEDS}, f"Unexpected run identity: {case}")
        lookup = {r["seed"]: r for r in records}
        scores = []
        case_outputs = []
        for seed in _SEEDS:
            record = lookup.get(seed)
            output = {"case_id": case, "ligand_id": row["ligand_id"], "receptor": _RECEPTOR,
                      "seed": seed, "status": "missing", "top1_score_kcal_mol": None}
            if record is not None:
                _require(record.get("case_id") == case and record.get("vina_version") == "AutoDock Vina v1.2.7",
                         f"Run case/version differs: {case}/{seed}")
                fingerprint = record.get("fingerprint", {})
                _require(set(fingerprint) == set(expected) | {"executable_sha256"}
                         and all(fingerprint.get(key) == value for key, value in expected.items()),
                         f"Run frozen input fingerprint differs: {case}/{seed}")
                _require(fingerprint.get("executable_sha256") in executable_hashes,
                         f"Run executable is not a frozen official Vina asset: {case}/{seed}")
                _require(shared_executable in (None, fingerprint["executable_sha256"]),
                         "Specificity runs did not use one identical executable")
                shared_executable = fingerprint["executable_sha256"]
                pose_name, log_name = f"{_RECEPTOR}_{seed}.pdbqt", f"{_RECEPTOR}_{seed}.log"
                _require(record.get("pose_file") == pose_name and record.get("log_file") == log_name,
                         f"Run artifact names differ: {case}/{seed}")
                _command(record, case, config, receptor["pdbqt_file"], ligand_file, pose_name)
                status = record.get("status")
                _require(status in {"ok", "failed", "running"}, f"Unknown run status: {case}/{seed}")
                output.update({"status": status, "input_fingerprint": fingerprint,
                               "command": record["command"], "error": record.get("error")})
                if status == "ok":
                    _require(record.get("exit_code") == 0, f"Successful run has nonzero exit: {case}/{seed}")
                    content = evidence.read(directory / pose_name, f"results/{case}/{pose_name}")
                    pose_hash = evidence.sources[f"results/{case}/{pose_name}"]["sha256"]
                    _require(pose_hash == record.get("pose_sha256"), f"Raw pose hash differs: {case}/{seed}")
                    with tempfile.TemporaryDirectory() as temporary:
                        snapshot = Path(temporary) / pose_name
                        snapshot.write_bytes(content)
                        poses = read_poses(snapshot, names)
                    _require(len(poses) == record.get("returned_poses") and 1 <= len(poses) <= config["num_modes"],
                             f"Returned pose count differs: {case}/{seed}")
                    log = evidence.read(directory / log_name, f"results/{case}/{log_name}")
                    _log(log, record, config, poses)
                    score = poses[0]["vina_score"]
                    scores.append(score)
                    output.update({"top1_score_kcal_mol": score, "returned_pose_count": len(poses),
                                   "pose_sha256": pose_hash,
                                   "log_sha256": evidence.sources[f"results/{case}/{log_name}"]["sha256"]})
            run_rows.append(output)
            case_outputs.append(output)
        complete = len(scores) == len(_SEEDS)
        compound_rows.append({"case_id": case, "ligand_id": row["ligand_id"], "pair_id": row["pair_id"],
                              "stereochemical_series": row["stereochemical_series"],
                              "aryl_substituent": row["aryl_substituent"], "binary_label": row["binary_label"],
                              "source_measurement": row["source_measurement"], "kd_uM": row["kd_uM"],
                              "kd_sd_uM": row["kd_sd_uM"], "numeric_kd_lower_bound_uM": None,
                              "heavy_atoms": row["expected_heavy_atoms"], "formal_charge": row["formal_charge"],
                              "successful_seed_count": len(scores), "expected_seed_count": 3,
                              "complete": complete,
                              "seed_top1_scores_kcal_mol": [r["top1_score_kcal_mol"] for r in case_outputs],
                              "mean_top1_score_kcal_mol": mean(scores) if complete else None})
    status_counts = Counter(r["status"] for r in run_rows)
    complete = status_counts["ok"] == 36
    metrics = summarize_specificity_scores(panel_path, {r["ligand_id"]: r["mean_top1_score_kcal_mol"]
                                                     for r in compound_rows}) if complete else None
    pair_counts = None
    if metrics:
        concordances = [r["concordance"] for r in metrics["per_pair"]]
        pair_counts = {"measured_binder_lower_energy": concordances.count(1.),
                       "exact_aggregate_ties": concordances.count(.5),
                       "conditional_NB_lower_energy": concordances.count(0.),
                       "pairs": 6, "tie_aware_concordance": mean(concordances)}
    summary = {"schema_version": 1, "panel_id": panel["panel_id"],
               "status": "evaluated_descriptive" if complete else "ineligible_incomplete_frozen_run_denominator",
               "information_regime": "known site from independent prior Gal9C glycan structure; shared across alternatives",
               "aggregation": config["score_aggregation"], "score_direction": "lower Vina energy predicts measured binding",
               "denominator": {"expected_ligands": 12, "measured_binders": 6, "conditional_NB_at_2mM": 6,
                               "matched_pairs": 6, "expected_seeds": _SEEDS, "expected_runs": 36,
                               "successful_runs": status_counts["ok"], "failed_runs": status_counts["failed"],
                               "running_runs": status_counts["running"], "missing_runs": status_counts["missing"],
                               "complete_ligands": sum(r["complete"] for r in compound_rows),
                               "metrics_eligible": complete},
               "compounds": compound_rows, "runs": run_rows, "metrics": metrics,
               "matched_pair_counts": pair_counts, "classifier_threshold": None,
               "source_hash_manifest": {"verified_artifacts": dict(sorted(evidence.sources.items())),
                                        "published_sources_declared_checksums": panel["source"]["source_file_checksums"],
                                        "publication_bytes_locally_reverified": False},
               "limitations": ["One target domain, assay context and synthetic glycoside series; exploratory development.",
                               "Gal9 and Gal3 share the galectin fold; this is not independent receptor-family generalization.",
                               "NB at highest tested 2mM is categorical; no exact or lower-bound KD is inferred.",
                               "No energy-to-KD calibration, fitted decision threshold or classifier performance claim.",
                               "Generated ligand input coordinates provide atom identity only, not structural pose ground truth.",
                               panel["target"]["construct_caveat"],
                               audit["decision"]["meaning"]]}
    return summary, render_specificity_markdown(summary, panel)


def _table(headers, rows):
    def cell(value):
        return str(value).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
                     + ["| " + " | ".join(cell(v) for v in row) + " |" for row in rows])


def render_specificity_markdown(summary, panel):
    """Render only recomputed statistics; never substitute NB with numeric KD."""
    d = summary["denominator"]
    text = ["# Gal9C matched-epimer specificity results", "",
            "This fixed analysis compares six compounds with finite measured affinities and six matched epimers reported as NB at the highest tested concentration, 2 mM. NB remains a categorical assay outcome; it is not a numerical KD bound.", "",
            f"The fixed denominator is 12 compounds × 3 seeds = 36 runs, covering six matched pairs. Verified successful runs: {d['successful_runs']}/36; failed: {d['failed_runs']}; running: {d['running_runs']}; missing: {d['missing_runs']}. Complete compounds: {d['complete_ligands']}/12.", "",
            "AutoDock Vina 1.2.7 uses CPU 2, exhaustiveness 16, up to 20 modes and the same 24 Å box for every compound. Seeds are 104729, 130363 and 155921. Each compound's score is the mean of its three seed-specific top-ranked energies. Lower energy predicts measured binding. The site comes from the independent prior 3NV3 glycan structure; this is a known-site experiment.", ""]
    if summary["metrics"]:
        metrics = summary["metrics"]
        raw = metrics["classification"]["contexts"][0]["raw"]
        counts = summary["matched_pair_counts"]
        quant = metrics["positive_only_quantitative_association"]["contexts"][0]
        spearman = f"{quant['spearman_rho']:.4f}" if quant["spearman_rho"] is not None else "undefined (constant ranks)"
        text.extend([f"Descriptive AUROC is {raw['auroc']:.4f}; step-wise average precision is {raw['auprc_average_precision']:.4f} at measured-binder prevalence 0.5000. No score threshold or classifier was fitted.", "",
                     f"Across six source-matched epimer pairs, measured binders have lower mean energy in {counts['measured_binder_lower_energy']}, conditional NB epimers have lower mean energy in {counts['conditional_NB_lower_energy']}, and {counts['exact_aggregate_ties']} are exact aggregate ties. Ties receive half credit; concordance is {counts['tie_aware_concordance']:.4f}.", "",
                     f"Among the six finite-affinity positives only, Spearman correlation of −mean Vina energy with measured pKD is {spearman}. NB compounds contribute no quantitative affinity value, correlation or censored numeric contrast.", "",
                     "| Descriptor baseline | AUROC | Average precision |", "| --- | --- | --- |"])
        for name, baseline in metrics["baselines"].items():
            value = baseline["contexts"][0]["raw"]
            text.append(f"| {name} | {value['auroc']:.4f} | {value['auprc_average_precision']:.4f} |")
        text.extend(["", "Heavy-atom-count AUROC is also 0.5000 in the reverse orientation because each matched pair has equal size. Formal charges are all zero.", ""])
    else:
        text.extend(["Ranking, matched-pair and quantitative metrics are ineligible until all 36 frozen runs succeed and pass verification. Missing or failed runs do not shrink the denominator.", ""])
    rows = []
    for compound in summary["compounds"]:
        measurement = (f"{compound['kd_uM']} ± {compound['kd_sd_uM']} µM"
                       if compound["binary_label"] else "NB at 2 mM")
        scores = ", ".join(f"{s:.3f}" if s is not None else "unavailable"
                           for s in compound["seed_top1_scores_kcal_mol"])
        aggregate = compound["mean_top1_score_kcal_mol"]
        rows.append([compound["ligand_id"], compound["stereochemical_series"], compound["aryl_substituent"],
                     measurement, scores, f"{aggregate:.4f}" if aggregate is not None else "ineligible",
                     f"{compound['successful_seed_count']}/3"])
    text.extend([_table(["Compound", "Series", "Aryl group", "Source assay result", "Seed top-1 energies (kcal/mol)",
                         "Mean energy (kcal/mol)", "Successful seeds"], rows), "",
                 "Finite affinity errors above are source-reported SD from 4–8 single/double-point measurements. Seed energies follow the frozen seed order; incomplete compounds have no mean.", ""])
    if summary["metrics"]:
        rows = [[p["binder_id"] + "/" + p["conditional_nonbinder_id"],
                 f"{p['binder_minus_conditional_nonbinder_kcal_mol']:.4f}", f"{p['concordance']:.1f}"]
                for p in summary["metrics"]["per_pair"]]
        text.extend([_table(["Measured/NB epimer pair", "Measured minus NB mean energy (kcal/mol)", "Concordance"], rows), ""])
    text.extend(["The 3NV3 declared 138-residue structural core exactly matches the corresponding assay sequence, including the annotated binding region. The model contains 136 observed residues: 13 N-terminal assay residues lie outside the declared structural core and two additional declared residues are unmodeled. Fusion/batch provenance remains approximate. One independently generated conformer is used per compound. No atomic pose truth for these alternatives is established and no input-coordinate RMSD is reported.", "",
                 "The frozen configuration retains unused fields inherited from the earlier structural diagnostic. The [field interpretation](../benchmarks/gal9_field_interpretation.json) identifies them without changing its frozen bytes. The explicit Gal9 panel, score aggregation and specificity metrics govern this analysis; pose thresholds are inactive.", "",
                 "These results describe one Gal9C assay context and one synthetic glycoside series. Gal9 and Gal3 share a fold; this provides no independent receptor-family generalization, calibrated affinity prediction or reliability certification.", "",
                 "All aggregates are recomputed from raw pose scores after checking the frozen experiment, source receptor, preparation/reference/ligand/receptor hashes, command options, raw pose hashes, full heavy-atom maps and raw log score tables. The JSON summary lists every expected run and a source hash manifest. Published assay-source checksums are preserved as declarations; publication bytes are not redistributed or claimed locally reverified.", "",
                 f"Assay source: [Mahanti et al., Table 1 and supporting methods](https://doi.org/{panel['source']['doi']}). Receptor qualification: [Gal9C receptor audit](../benchmarks/gal9_receptor_audit.json). Frozen protocol: [gal9_experiment.json](../benchmarks/gal9_experiment.json). Machine-readable evidence: [summary.json](../results/gal9_specificity/summary.json).", ""])
    return "\n".join(text)


def _first_difference(left, right, path="root"):
    """Describe a difference for diagnostics; comparison policy lives elsewhere."""
    if type(left) is not type(right):
        return f"{path}: types {type(left).__name__}/{type(right).__name__}"
    if isinstance(left, dict):
        if left.keys() != right.keys():
            return f"{path}: object fields differ"
        for key in sorted(left):
            difference = _first_difference(left[key], right[key], f"{path}.{key}")
            if difference:
                return difference
    elif isinstance(left, list):
        if len(left) != len(right):
            return f"{path}: lengths {len(left)}/{len(right)}"
        for index, (a, b) in enumerate(zip(left, right)):
            difference = _first_difference(a, b, f"{path}[{index}]")
            if difference:
                return difference
    elif left != right:
        if isinstance(left, float):
            return f"{path}: {left!r}/{right!r}; absolute difference {abs(left - right):.17g}"
        return f"{path}: {str(left)[:80]!r}/{str(right)[:80]!r}"
    return None


def write_specificity_report(root, *, prepared_root=None, results_root=None,
                             summary_output=None, markdown_output=None, check=False,
                             numeric_tolerance=0.0):
    """Write/check JSON and Markdown; optional tolerance covers derived floats only.

    Raw scores, measured affinities, identities, counts and input hashes remain
    exact. The generated Markdown is always checked byte for byte.
    """
    _require(not isinstance(numeric_tolerance, bool) and math.isfinite(numeric_tolerance)
             and 0 <= numeric_tolerance <= 1e-8, "Numeric tolerance must be finite and between 0 and 1e-8")
    root = Path(root).resolve()
    summary, markdown = generate_specificity_report(root, prepared_root=prepared_root, results_root=results_root)
    summary_path = Path(summary_output) if summary_output else root / "results/gal9_specificity/summary.json"
    markdown_path = Path(markdown_output) if markdown_output else root / "docs/GAL9_SPECIFICITY_RESULTS.md"
    outputs = [(summary_path, json.dumps(summary, sort_keys=True, indent=2, allow_nan=False) + "\n"),
               (markdown_path, markdown)]
    for path, content in outputs:
        if check:
            _require(path.is_file(), f"Specificity report is stale: {path}")
            if path == summary_path:
                with tempfile.TemporaryDirectory() as temporary:
                    regenerated = Path(temporary) / "summary.json"
                    regenerated.write_text(content)
                    try:
                        comparison = compare_tables(path, regenerated, numeric_tolerance)
                    except ValueError as error:
                        try:
                            difference = _first_difference(json.loads(path.read_text()), summary)
                        except (ValueError, TypeError):
                            difference = "checked-in JSON cannot be parsed"
                        raise ValueError(f"Specificity report is stale: {path}; {error}; first difference: {difference}") from error
                if not comparison["exact_recomputation"]:
                    print(f"Specificity JSON derived-float comparison: maximum difference "
                          f"{comparison['max_numeric_difference']:.17g}; explicit tolerance {numeric_tolerance:g}")
            else:
                _require(path.read_text() == content, f"Specificity report is stale: {path}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            write_text_atomic(path, content)
    return summary_path, markdown_path
