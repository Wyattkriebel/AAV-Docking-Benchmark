"""Portable prospective report for paired first-chain assembly controls.

Predicted full-versus-monomer arms are an end-to-end workflow sensitivity:
their source-only predictor centers may change. Oracle centers are held fixed.
Strict-versus-repaired comparisons within each assembly use identical ranks and
centers. All scoring retains the full observed experimental assembly as truth.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from statistics import mean

from .heldout_pose_report import CASE_IDS as FULL_CASE_IDS, _source_context, protein_fits, score_pool, site_summary
from .paired_method_report import (CODE_SOURCES as SHARED_CODE_SOURCES, check_public_reference_sources,
    common_selection, portable_archive, read_modern)
from .provenance import sha256, write_text_atomic
from .qualified_pose import Evidence, accuracy_summary
from .statistics import paired_family_bootstrap


CASE_IDS = ("3ZWE", "5HZA")
METHODS = ("fpocket", "p2rank", "oracle")
SEEDS = (104729, 130363, 155921)
PROFILE_PAIRS = {"strict_monomer": "strict", "repair_monomer": "repair"}
CENTER_POLICY = {
    "predicted_full_vs_monomer": "end_to_end_workflow_allow_predictor_centers_to_change",
    "oracle_full_vs_monomer": "same_full_native_parent_center",
    "strict_vs_repair_within_assembly": "exact_same_candidate_center_and_rank",
}
CODE_SOURCES = ("aav_benchmark/assembly_control_report.py", *SHARED_CODE_SOURCES)


def _pin(evidence, record):
    return evidence.read(record["path"], record["sha256"])


def _json(raw):
    return json.loads(raw)


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


class ControlReceiptView:
    """Read-only field alias for the frozen control driver's index binding.

    The shared modern reader uses ``frozen_input_manifest_sha256``. The control
    driver records the same source-index binding as ``prepared_conditions_sha256``.
    Only the designated receipt is normalized in memory; Evidence inventories
    the actual original bytes and no file is changed.
    """
    def __init__(self, evidence, receipt_path, prepared_conditions_sha256):
        self.evidence = evidence
        self.receipt_path = Path(receipt_path)
        self.expected = prepared_conditions_sha256
        self.root, self.sources = evidence.root, evidence.sources

    def read(self, relative, expected=None):
        return self.evidence.read(relative, expected)

    def json(self, relative, expected=None):
        value = self.evidence.json(relative, expected)
        if Path(relative) == self.receipt_path:
            if value.get("prepared_conditions_sha256") != self.expected:
                raise ValueError("Control execution receipt changed its frozen prepared-condition index")
            if "frozen_input_manifest_sha256" in value:
                raise ValueError("Control execution receipt has an unexpected duplicate index-binding field")
            value = {**value, "frozen_input_manifest_sha256": self.expected}
        return value


def condition_row(evidence, descriptor, condition_pin, configuration_sha256):
    """Build the existing modern-reader contract from one pinned preparation."""
    condition = _json(_pin(evidence, condition_pin))
    if (condition["case_id"] != descriptor["case_id"] or condition["method"] != descriptor["method"]
            or condition["profile"] != descriptor["profile"] or condition["configuration_sha256"] != configuration_sha256
            or condition.get("fixed_denominator_member") is not True
            or condition.get("executed_primary_candidate_budget") != 1):
        raise ValueError("Assembly condition changed its frozen case/profile/site/budget")
    jobs = condition["jobs"]
    if len(jobs) > 1 or any(j["candidate_class_rank"] != 1 for j in jobs):
        raise ValueError("Assembly control must retain its single top1 candidate class")
    row = {"condition_id": descriptor["condition_id"], "case_id": descriptor["case_id"],
           "method": descriptor["method"], "receptor_profile": descriptor["profile"],
           "source_condition": condition_pin, "modern_runs": descriptor["modern_runs"]}
    if not jobs:
        row.update(status="predictor_unavailable", failure_reason=condition["status"])
    elif jobs[0]["status"] != "prepared":
        if jobs[0]["status"] not in {"preparation_failed", "unsupported"}:
            raise ValueError("Assembly preparation is nonterminal")
        row.update(status="preparation_failed", failure_reason=jobs[0].get("reason", jobs[0]["status"]))
        evidence.read(Path(condition_pin["path"]).parent / "box01/preparation_attempt.json")
    else:
        job = jobs[0]
        directory = Path(condition_pin["path"]).parent / job["prepared_directory"]
        if job["prepared_directory"] != "box01" or job["id"] != "box01":
            raise ValueError("Unsafe assembly prepared-box identity")
        raw = evidence.read(directory / "preparation.json", job["preparation_sha256"])
        prep = _json(raw)
        row.update(status="prepared", candidate_class_rank=1, candidate_original_rank=job["candidate_original_rank"],
            modern_vina_prepared_directory=directory.as_posix(), box={"center": prep["box_center"], "size": prep["box_size"]},
            preparation={"path": (directory / "preparation.json").as_posix(), "sha256": _hash(raw)})
        for field, filename in (("ligand", "ligand.pdbqt"), ("receptor", "neighborhood.pdbqt"), ("canonical_atom_map", "ligand_atom_map.json")):
            content = evidence.read(directory / filename)
            row[field] = {"path": (directory / filename).as_posix(), "sha256": _hash(content)}
        if (prep["case_id"] != descriptor["case_id"] or prep["profile"] != descriptor["profile"]
                or prep["box_center"] != job["center_A"] or prep["configuration_sha256"] != configuration_sha256
                or prep.get("native_coordinates_in_runner_reference") is not False):
            raise ValueError("Assembly box changed identity, center or source-only runner contract")
    return row, condition


def validate_center_policy(conditions, policy):
    """Check controlled pairs without falsely equating end-to-end predictions."""
    if policy != CENTER_POLICY:
        raise ValueError("Assembly pairing-center policy differs from its predeclared semantics")
    by_key = {(c["profile"], c["case_id"], c["method"]): c for c in conditions}
    if len(by_key) != len(conditions):
        raise ValueError("Duplicate assembly comparison condition")
    def identity(item):
        return [(c["rank"], c["original_rank"], c["center_A"]) for c in item["candidates"]]
    for cid in CASE_IDS:
        for method in METHODS:
            for left, right in (("strict", "repair"), ("strict_monomer", "repair_monomer")):
                if identity(by_key[left, cid, method]) != identity(by_key[right, cid, method]):
                    raise ValueError("Repair sensitivity changed original candidate centers or ranks")
            if method == "oracle":
                if identity(by_key["strict", cid, method]) != identity(by_key["strict_monomer", cid, method]):
                    raise ValueError("Oracle assembly contrast changed its common native-informed center")
    return [{"case_id": cid, "method": method,
             "interpretation": "controlled_receptor_only_conditional_contrast" if method == "oracle" else "end_to_end_site_input_and_assembly_workflow_sensitivity",
             "full_candidates": by_key["strict", cid, method]["candidates"],
             "monomer_candidates": by_key["strict_monomer", cid, method]["candidates"],
             "exact_candidate_center_and_rank_match": identity(by_key["strict", cid, method]) == identity(by_key["strict_monomer", cid, method])}
            for cid in CASE_IDS for method in METHODS]


def validate_source_pairs(sources, observed):
    """Require fixed ligand bytes and the complete observed assembly as truth."""
    for cid in CASE_IDS:
        for name in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json"):
            baseline = sources["strict"][f"{cid}/{name}"]
            if any(source[f"{cid}/{name}"] != baseline for source in sources.values()):
                raise ValueError("Assembly comparison changed independent parent chemistry/initializer")
        target = _json(observed[f"{cid}/protein_atom_map.json"])
        target_atoms = {tuple(a["output_key"]): a for a in target["atoms"]}
        for profile in ("strict_monomer", "repair_monomer"):
            source = sources[profile]
            if (source[f"{cid}/full_protein_atom_map.json"] != observed[f"{cid}/protein_atom_map.json"]
                    or source[f"{cid}/full_assembly_receptor.cif"] != observed[f"{cid}/receptor.cif"]):
                raise ValueError("Monomer comparison replaced or truncated the full observed assembly truth")
            subset = _json(source[f"{cid}/protein_atom_map.json"])
            subset_atoms = {tuple(a["output_key"]): a for a in subset["atoms"]}
            expected = {key for key in target_atoms if key[0] == "A"}
            if set(subset_atoms) != expected or any(subset_atoms[k] != target_atoms[k] for k in expected):
                raise ValueError("First-chain observed atom map is not the exact fixed source-A subset")


def _load_prediction_sources(evidence, configuration, freeze_record, *, monomer=False):
    config = _json(_pin(evidence, configuration))
    receipt = _json(_pin(evidence, freeze_record))
    key = "frozen_before_monomer_predictions" if monomer else "frozen_before_new_predictions"
    if receipt.get("configuration_sha256") != configuration["sha256"] or receipt.get(key) is not True:
        raise ValueError("Assembly comparison site predictions lack a pre-prediction freeze")
    if monomer and (receipt.get("fixed_chain_rule_declared_before_working_full_predictions") is not True
                    or receipt.get("new_test_site_or_pose_accuracy_consulted") is not False):
        raise ValueError("First-chain selection was not fixed before working site predictions/accuracy")
    archive = config["source_preparation_archive"]
    sources = portable_archive(evidence, {"path": archive["manifest"], "sha256": archive["manifest_sha256"]}, archive_path=archive["path"])
    return sources


def comparison_statistics(conditions, case_to_family):
    output = []
    for mono, full in PROFILE_PAIRS.items():
        for method in METHODS:
            for alignment in ("global", "local_interface"):
                for rank in ("top1", "top5", "best_returned"):
                    for threshold in (1, 2, 3):
                        selected = [c for c in conditions if c["profile"] in {mono, full} and c["method"] == method]
                        records = [{"family_id": case_to_family[c["case_id"]], "complex_id": c["case_id"], "method": c["profile"],
                                    "seed": s["seed"], "status": "ok" if s["poses"] else "failed",
                                    "success": s[alignment][f"{rank}_success_{threshold}A"]} for c in selected for s in c["seeds"]]
                        value = paired_family_bootstrap(records, mono, full, n_resamples=10000, seed=20261103, expected_seeds=list(SEEDS))
                        value.update(site_condition=method, alignment=alignment, pose_rank=rank, threshold_A=threshold,
                            interpretation="Conditional two-source-family assembly-control sensitivity; no independent biological samples or population reliability claim",
                            pairing_center_policy=CENTER_POLICY["oracle_full_vs_monomer" if method == "oracle" else "predicted_full_vs_monomer"])
                        output.append(value)
    return output


def generate_report(root, plan_path, freeze_path, *, allow_pending=False):
    evidence = Evidence(root)
    plan, freeze = evidence.json(plan_path), evidence.json(freeze_path)
    if (freeze.get("evaluation_plan_sha256") != sha256(evidence.root / plan_path)
            or freeze.get("frozen_before_test_pose_accuracy") is not True or freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Assembly report requires its explicit pre-accuracy evaluation freeze")
    if (plan["case_ids"] != list(CASE_IDS) or plan["methods"] != list(METHODS) or plan["seeds"] != list(SEEDS)
            or plan["profile_pairs"] != PROFILE_PAIRS or plan["pairing_center_policy"] != CENTER_POLICY):
        raise ValueError("Fixed assembly-control case/profile/center design changed")
    if set(plan["metric_code_sha256"]) != set(CODE_SOURCES):
        raise ValueError("Assembly report metric/parser source pins are incomplete")
    for path, digest in plan["metric_code_sha256"].items():
        evidence.read(path, digest)
    paired = _json(_pin(evidence, plan["paired_evaluation_plan"]))
    paired_freeze = _json(_pin(evidence, plan["paired_evaluation_freeze"]))
    if (paired_freeze.get("evaluation_plan_sha256") != plan["paired_evaluation_plan"]["sha256"]
            or paired_freeze.get("frozen_before_test_pose_accuracy") is not True
            or paired_freeze.get("new_test_pose_accuracy_consulted") is not False
            or paired["reference_evaluation_plan"] != plan["reference_evaluation_plan"]
            or paired["reference_evaluation_freeze"] != plan["reference_evaluation_freeze"]
            or paired["public_source_distribution"] != plan["public_source_distribution"]
            or paired["source_profiles"] != plan["full_source_profiles"]
            or any(paired["metric_code_sha256"].get(path) != plan["metric_code_sha256"][path] for path in SHARED_CODE_SOURCES)):
        raise ValueError("Assembly comparator reference/source validation differs from its frozen paired evaluation")
    reference = _json(_pin(evidence, plan["reference_evaluation_plan"]))
    reference_freeze = _json(_pin(evidence, plan["reference_evaluation_freeze"]))
    if (reference_freeze.get("evaluation_plan_sha256") != plan["reference_evaluation_plan"]["sha256"]
            or reference_freeze.get("frozen_before_test_pose_accuracy") is not True or reference_freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Shared reference plan changed after its pre-accuracy freeze")
    cases, source_validation = check_public_reference_sources(evidence, reference["source_ledger"],
        reference["reference_source_manifest"], plan["public_source_distribution"])
    _pin(evidence, reference["family_audit"])
    assembly_config = _json(_pin(evidence, plan["assembly_configuration"]))
    assembly_freeze = _json(_pin(evidence, plan["assembly_freeze"]))
    full_config = _json(_pin(evidence, plan["full_configuration"]))
    full_freeze = _json(_pin(evidence, plan["full_freeze"]))
    for config, receipt, pin in ((assembly_config, assembly_freeze, plan["assembly_configuration"]),
                                 (full_config, full_freeze, plan["full_configuration"])):
        if receipt.get("configuration_sha256") != pin["sha256"] or receipt.get("frozen_before_new_docking") is not True:
            raise ValueError("Assembly/full scientific configurations require pre-execution freezes")
    for field in ("cpu", "exhaustiveness", "num_modes", "energy_range_kcal_mol", "min_rmsd_A", "seeds"):
        if assembly_config[field] != full_config[field]:
            raise ValueError("Assembly contrast changed matched Vina search settings")
    if (assembly_config["case_ids"] != list(CASE_IDS) or assembly_config["receptor_profiles"] != list(PROFILE_PAIRS)
            or assembly_config["primary_candidate_budget"] != 1 or assembly_config["seeds"] != list(SEEDS)):
        raise ValueError("Native assembly configuration changed fixed controls")
    observed = portable_archive(evidence, assembly_config["evaluation_source_archive_manifest"])
    sources = {}
    for profile in PROFILE_PAIRS:
        sources[profile] = portable_archive(evidence, assembly_config["source_profiles"][profile]["source_archive_manifest"])
        profile_manifest = sources[profile]["profile_manifest.json"]
        if (_hash(profile_manifest) != assembly_config["source_profiles"][profile]["manifest"]["sha256"]
                or _json(profile_manifest)["profile"] != profile or _json(profile_manifest)["case_ids"] != list(CASE_IDS)):
            raise ValueError("Archived assembly source profile differs from its frozen native input manifest")
    source_profiles = {p["id"]: p for p in reference["profiles"]}
    for full, profile_id in plan["full_source_profiles"].items():
        sources[full] = portable_archive(evidence, source_profiles[profile_id]["source_archive_manifest"])
    if set(sources) != {"strict", "repair", "strict_monomer", "repair_monomer"}:
        raise ValueError("Assembly report source profiles are missing or duplicated")
    validate_source_pairs(sources, observed)
    monomer_predictions = _load_prediction_sources(evidence, assembly_config["prediction_configuration"], assembly_config["prediction_freeze"], monomer=True)
    full_predictions = _load_prediction_sources(evidence, reference["site_prediction"]["configuration"], reference["site_prediction"]["freeze_receipt"])
    assembly_index = _json(_pin(evidence, plan["assembly_preparation_index"]))
    if (assembly_index["status"] != "complete" or assembly_index["condition_count"] != 12
            or assembly_index.get("expected_conditions") != 12 or assembly_index.get("accuracy_evaluated") is not False):
        raise ValueError("Assembly source preparation index is not a complete source-only twelve-condition inventory")
    index = {(r["receptor_profile"], r["case_id"], r["method"]): r for r in assembly_index["conditions"]}
    expected = {(p, cid, m) for p in PROFILE_PAIRS for cid in CASE_IDS for m in METHODS}
    descriptors = {(r["profile"], r["case_id"], r["method"]): r for r in plan["assembly_conditions"]}
    if (len(index) != len(assembly_index["conditions"]) or len(index) != 12 or set(index) != expected
            or len(descriptors) != len(plan["assembly_conditions"]) or len(descriptors) != 12 or set(descriptors) != expected):
        raise ValueError("Keep all twelve assembly-control method conditions in the denominator")
    full_inputs = _json(_pin(evidence, plan["full_prepared_inputs"]))
    full_rows = {(r["receptor_profile"], r["case_id"], r["method"]): r for r in full_inputs["conditions"]}
    full_expected = {(p, cid, method) for p in ("strict", "repair") for cid in FULL_CASE_IDS for method in METHODS}
    if (len(full_rows) != len(full_inputs["conditions"]) or len(full_rows) != 36 or set(full_rows) != full_expected
            or full_inputs.get("source_only_index") is not True or full_inputs.get("docking_or_accuracy_outcomes_consulted") is not False
            or full_inputs.get("primary_candidate_budget") != 1 or full_inputs.get("secondary_top5_enabled") is not False
            or full_inputs["modern_vina_configuration"] != plan["full_configuration"]):
        raise ValueError("Full comparator inputs changed their frozen source-only36-condition inventory")
    loaded, pending = [], []
    counts = {"control_conditions": 12, "comparator_conditions": 12, "eligible_seed_outcomes": 72,
              "terminal_seed_outcomes": 0, "successful_seed_outcomes": 0, "prepared_input_failure_seed_outcomes": 0}
    for profile, cid, method in sorted(expected):
        entry, descriptor = index[profile, cid, method], descriptors[profile, cid, method]
        condition_pin = {"path": entry["condition_path"], "sha256": entry["condition_sha256"]}
        row, condition = condition_row(evidence, descriptor, condition_pin, plan["assembly_configuration"]["sha256"])
        if descriptor["preparation"] != condition_pin:
            raise ValueError("Assembly descriptor/source preparation index disagree")
        if entry["status"] != condition["status"] or entry["prepared_jobs"] != sum(j["status"] == "prepared" for j in condition["jobs"]):
            raise ValueError("Assembly index changed actual preparation support/failure counts")
        if row["status"] == "prepared":
            prep = _json(_pin(evidence, row["preparation"]))
            if (prep["source_receptor_sha256"] != _hash(sources[profile][f"{cid}/receptor.cif"])
                    or prep["independent_ligand_sha256"] != _hash(sources[profile][f"{cid}/ligand.pdbqt"])
                    or prep["full_assembly_evaluation_archive_manifest"] != assembly_config["evaluation_source_archive_manifest"]
                    or row["box"]["size"] != _json(sources[profile][f"{cid}/chemistry.json"])["box_size"]):
                raise ValueError("Control prepared box changed its source protein, full parent or full-assembly truth binding")
        _, candidates = common_selection(evidence, row, assembly_config["prediction_configuration"]["sha256"], monomer_predictions)
        loaded.append({"profile": profile, "case_id": cid, "method": method, "row": row, "candidates": candidates,
                       "config": assembly_config, "config_pin": plan["assembly_configuration"], "input_index_sha256": plan["assembly_preparation_index"]["sha256"]})
        full_profile = PROFILE_PAIRS[profile]
        full_row = full_rows[full_profile, cid, method]
        _, full_candidates = common_selection(evidence, full_row, reference["site_prediction"]["configuration"]["sha256"], full_predictions)
        loaded.append({"profile": full_profile, "case_id": cid, "method": method, "row": full_row, "candidates": full_candidates,
                       "config": full_config, "config_pin": plan["full_configuration"], "input_index_sha256": plan["full_prepared_inputs"]["sha256"]})
    center_comparisons = validate_center_policy(loaded, plan["pairing_center_policy"])
    for item in loaded:
        cid, profile, row = item["case_id"], item["profile"], item["row"]
        atoms = _json(sources[profile][f"{cid}/ligand_atom_map.json"])["atoms"]
        reader_evidence = evidence
        if profile in PROFILE_PAIRS:
            receipt_path = Path(row["modern_runs"]).parent.parent / "execution.json"
            reader_evidence = ControlReceiptView(evidence, receipt_path, item["input_index_sha256"])
        outcomes, unfinished = read_modern(reader_evidence, row, item["config"], item["config_pin"], [a["name"] for a in atoms],
            input_manifest_sha256=item["input_index_sha256"])
        pending.extend(f"{profile}/{cid}/{item['method']}/{p}" for p in unfinished)
        item["outcomes"] = outcomes
        for result in outcomes.values():
            counts["terminal_seed_outcomes"] += 1
            counts["successful_seed_outcomes"] += result["status"] == "ok"
            counts["prepared_input_failure_seed_outcomes"] += result["status"] == "prepared_input_failure"
    if pending:
        if not allow_pending:
            raise ValueError(f"Assembly/control comparator execution remains incomplete ({len(pending)} pending records); no accuracy report")
        return {"schema_version": 1, "status": "pending_incomplete_workflow", "metrics_eligible": False, "counts": counts,
                "pending": sorted(pending), "reference_validation": source_validation, "source_sha256_inventory": evidence.sources,
                "control_receipt_binding_field": "prepared_conditions_sha256; read-only in-memory alias, original receipt bytes preserved"}
    if counts["terminal_seed_outcomes"] != 72:
        raise ValueError("All72 assembly/control-comparator seed outcomes must be terminal")
    contexts = {(profile, cid): _source_context(source, cid, cases[cid], evaluation_source=observed)
                for profile, source in sources.items() for cid in CASE_IDS}
    conditions = []
    for item in loaded:
        profile, cid = item["profile"], item["case_id"]
        mapping, graph, keys, protein, protein_map = contexts[profile, cid]
        allow_missing = profile in PROFILE_PAIRS
        seeds, controls = [], None
        for seed in SEEDS:
            result = item["outcomes"][seed]
            raw = [{"rank": p["rank"], "vina_score": p["method_score"], "coordinates": p["coordinates"],
                    "candidate_class_rank": 1, "candidate_original_rank": item["row"].get("candidate_original_rank", 1),
                    "box_id": "box01", "pooled_rank": p["rank"]} for p in result["poses"]]
            scored, controls = score_pool(raw, mapping, keys, protein, protein_map, graph, allow_missing_observed=allow_missing)
            seeds.append({"seed": seed, "status": result["status"], "poses": scored,
                          "global": accuracy_summary(scored, "global_RMSD_A"), "local_interface": accuracy_summary(scored, "local_interface_RMSD_A")})
        fit = protein_fits(keys, protein, protein_map, mapping["copies"][0]["coordinates"], allow_missing_observed=allow_missing)["fits"]["global"]
        chemistry = _json(sources[profile][f"{cid}/chemistry.json"])
        site = site_summary(item["candidates"], mapping["copies"], box_sizes=chemistry["box_size"], rotation=fit["rotation"], translation=fit["translation"])
        conditions.append({"profile": profile, "case_id": cid, "method": item["method"], "atomic_scope": cases[cid]["qualification"],
            "scored_atoms": cases[cid]["scored_atom_count"], "parent_atoms": mapping["parent_atom_count"], "native_controls": controls,
            "source_preparation_status": item["row"]["status"], "site": site, "seeds": seeds})
    return {"schema_version": 1, "status": "complete", "metrics_eligible": True, "counts": counts, "conditions": conditions,
            "full_ligand_cases": 2, "modeled_core_cases": 0,
            "center_comparisons": center_comparisons, "pairing_center_policy": CENTER_POLICY,
            "paired_family_statistics": comparison_statistics(conditions, reference["case_to_family"]),
            "reference_validation": source_validation, "source_sha256_inventory": evidence.sources, "metric_code_sha256": plan["metric_code_sha256"],
            "control_receipt_binding_field": "prepared_conditions_sha256; read-only in-memory alias, original receipt bytes preserved",
            "limitations": ["Both controls are full-ligand qualified; no modeled-core recovery rate is inferred from this zero-core subset.",
                "Predicted full-versus-monomer comparisons change site input and may change centers; this is an end-to-end workflow assembly sensitivity.",
                "Oracle centers are held fixed; only those contrasts isolate receptor assembly conditional on the known site.",
                "Strict-versus-repaired comparisons within each assembly reuse exactly the same candidate centers and ranks.",
                "Omitted chains remain in full observed native contacts and clash truth. Generated atoms never become experimental truth.",
                "Two source families and incomplete monomer controls do not constitute new biological samples or population reliability evidence."]}


def render_markdown(report):
    if not report["metrics_eligible"]:
        return "# Assembly-control report pending\n\nAll planned control and comparator outcomes must be terminal. No site or pose accuracy is reported.\n"
    text = ["# Paired first-chain assembly controls", "", "The 3ZWE and 5HZA controls retain their full experimental assemblies as truth. Predicted arms measure an end-to-end workflow sensitivity; oracle centers are fixed.", "",
            "| Profile | Site condition | Full-ligand cases | Top1 pose ≤2Å | Top5 pose ≤2Å | Best retained ≤2Å |",
            "|---|---|---:|---:|---:|---:|"]
    for profile in ("strict", "strict_monomer", "repair", "repair_monomer"):
        for method in METHODS:
            selected = [r for r in report["conditions"] if r["profile"] == profile and r["method"] == method]
            values = [mean(mean(s["global"][f"{rank}_success_2A"] for s in r["seeds"]) for r in selected) for rank in ("top1", "top5", "best_returned")]
            text.append(f"| {profile} | {method} | {len(selected)} | " + " | ".join(f"{v:.3f}" for v in values) + " |")
    text.extend(["", "Top5 poses are within the single executed candidate box. The JSON preserves every case/seed, global/local alignments, 1/2/3Å endpoints, source-native controls, full-assembly contact/clash diagnostics and input-center comparisons.", ""])
    text.extend(report["limitations"])
    return "\n".join(text) + "\n"


def write_report(root, plan_path, freeze_path, summary_output, markdown_output, *, check=False, allow_pending=False):
    report = generate_report(root, plan_path, freeze_path, allow_pending=allow_pending)
    for path, content in ((Path(root) / summary_output, json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"),
                          (Path(root) / markdown_output, render_markdown(report))):
        if check:
            if path.read_text() != content:
                raise ValueError(f"Assembly-control report replay changed: {path}")
        else:
            write_text_atomic(path, content)
    return report


def example_plan(root=None):
    root = Path(root) if root is not None else None
    def pin(path):
        return {"path": path, "sha256": sha256(root / path) if root is not None and (root / path).is_file() else "FILL_SOURCE_SHA256"}
    conditions = [{"profile": profile, "case_id": cid, "method": method,
        "condition_id": f"{profile}/{cid}/{method}/box01",
        "preparation": pin(f"results/structural_assembly_pose_v3/{profile}/{cid}/{method}/prepared/condition.json"),
        "modern_runs": f"results/structural_assembly_pose_v3/{profile}/{cid}/{method}/vina/box01/runs.json"}
        for profile in PROFILE_PAIRS for cid in CASE_IDS for method in METHODS]
    return {"schema_version": 1, "case_ids": list(CASE_IDS), "methods": list(METHODS), "seeds": list(SEEDS),
        "profile_pairs": PROFILE_PAIRS, "pairing_center_policy": CENTER_POLICY,
        "reference_evaluation_plan": pin("benchmarks/structural_pose_reference_evaluation.json"),
        "reference_evaluation_freeze": pin("results/structural_pose_v3/reference_evaluation_freeze.json"),
        "paired_evaluation_plan": pin("benchmarks/structural_paired_method_evaluation_v2.json"),
        "paired_evaluation_freeze": pin("results/structural_pose_v3/paired_evaluation_v2_freeze.json"),
        "public_source_distribution": pin("benchmarks/source_public_distribution.json"),
        "assembly_configuration": pin("benchmarks/structural_assembly_pose_v3_experiment.json"),
        "assembly_freeze": pin("results/structural_assembly_pose_v3/freeze_receipt.json"),
        "assembly_preparation_index": pin("results/structural_assembly_pose_v3/prepared_conditions.json"),
        "full_configuration": pin("benchmarks/structural_pose_v3_experiment.json"), "full_freeze": pin("results/structural_pose_v3/freeze_receipt.json"),
        "full_prepared_inputs": pin("results/structural_pose_v3/vina_carb_prepared_inputs.json"),
        "full_source_profiles": {"strict": "strict_CCD_ideal", "repair": "modeled_sidechain_CCD_ideal"},
        "assembly_conditions": conditions, "metric_code_sha256": {path: pin(path)["sha256"] for path in CODE_SOURCES},
        "freeze_receipt_schema": {"evaluation_plan_sha256": "SHA256_OF_COMPLETED_PLAN", "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False}}
