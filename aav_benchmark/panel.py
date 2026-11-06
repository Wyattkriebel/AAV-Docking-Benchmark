"""Descriptive assay ranking; generated input coordinates are never pose truth."""

import json
from pathlib import Path
from statistics import mean

from .engine import read_poses
from .provenance import sha256, write_json
from .statistics import quantitative_specificity_by_context


def summarize_panel(prepared_root, results_root, panel_path):
    prepared_root, results_root = Path(prepared_root), Path(results_root)
    panel = json.loads(Path(panel_path).read_text())
    summaries, correlations = [], []
    for anomer in ("alpha", "beta"):
        measured, size_baseline = [], []
        for row in panel["rows"]:
            case_id = f"gal3_{row['ligand_id']}_{anomer}"
            prepared, output = prepared_root / case_id, results_root / case_id
            metadata = json.loads((prepared / "preparation.json").read_text())
            reference = json.loads((prepared / "reference.json").read_text())
            config_path = output / "configuration.json"
            config = json.loads(config_path.read_text())
            runs = json.loads((output / "runs.json").read_text())
            expected = {(r["id"], seed) for r in metadata["receptors"] for seed in config["seeds"]}
            actual = [(r["receptor"], r["seed"]) for r in runs]
            if len(set(actual)) != len(actual) or set(actual) != expected:
                raise ValueError("Incomplete or duplicate assay receptor-by-seed denominator")
            scores = []
            for run in runs:
                receptor = next(r for r in metadata["receptors"] if r["id"] == run["receptor"])
                expected_hashes = {"config_sha256": sha256(config_path),
                                   "preparation_sha256": sha256(prepared / "preparation.json"),
                                   "reference_sha256": sha256(prepared / "reference.json"),
                                   "ligand_sha256": sha256(prepared / metadata["ligand_file"]),
                                   "receptor_sha256": sha256(prepared / receptor["pdbqt_file"]),
                                   "receptor_pdb_sha256": sha256(prepared / receptor["pdb_file"])}
                if any(run["fingerprint"].get(k) != v for k, v in expected_hashes.items()):
                    raise ValueError("Assay scoring inputs changed since docking")
                if run["status"] != "ok":
                    raise ValueError(f"Assay run failed: {case_id}/{run['seed']}; ranking requires all frozen runs")
                pose_path = output / run["pose_file"]
                if sha256(pose_path) != run["pose_sha256"]:
                    raise ValueError("Assay pose output checksum mismatch")
                scores.append(read_poses(pose_path, reference["names"])[0]["vina_score"])
            entry = {"case_id": case_id, "ligand_id": row["ligand_id"], "anomer": anomer,
                     "top1_scores_kcal_mol": scores, "seed_mean_top1_score": mean(scores),
                     "kd_uM": row["kd_uM"], "kd_sem_uM": row["kd_sem_uM"],
                     "heavy_atoms": row["expected_heavy_atoms"], "eligible_runs": len(expected), "completed_runs": len(scores)}
            summaries.append(entry)
            common = {"context_id": panel["panel_id"], "target_id": panel["target"]["uniprot_id"],
                      "ligand_id": row["ligand_id"], "kd": row["kd_uM"], "kd_unit": "uM"}
            measured.append({**common, "score": -mean(scores)})
            size_baseline.append({**common, "score": row["expected_heavy_atoms"]})
        correlations.append({"anomer": anomer, "vina": quantitative_specificity_by_context(measured),
                             "heavy_atom_baseline": quantitative_specificity_by_context(size_baseline)})
    result = {"panel_id": panel["panel_id"], "role": "exploratory_development", "ligands": summaries,
              "correlations": correlations, "binary_classification": "ineligible: all three compounds are binders",
              "limitations": ["one target; three ligand identities; no held-out evaluation", panel["target"]["construct_mismatch"],
                              "unknown solution anomer weights; anomers are sensitivity conditions, not extra observations",
                              "size alone perfectly ranks this panel; matching is not available at n=3",
                              "experimental PBS pH/ionic composition unreported; rigid dry receptor and modeled protonation"]}
    write_json(results_root / "gal3_panel_summary.json", result)
    return result
