"""Portable external metrics for the separately declared guided HADDOCK assay.

No native software is invoked or imported. Numerical models, full-precision
individual score metadata and frozen inputs are required locally; installation
hashes remain historical provenance. Every eligible source failure stays in the
denominator. Experimental accuracy is unavailable until all108 cells terminate.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import tarfile
import tempfile

import numpy as np

from .haddock_paired_cohort import model_inventory
from .heldout_pose_report import CASE_IDS, METHODS, _source_context, accuracy_summary, score_pool
from .paired_method_report import (CODE_SOURCES as SHARED_CODE_SOURCES, check_public_reference_sources,
                                  common_selection, portable_archive)
from .provenance import sha256, write_text_atomic
from .qualified_pose import Evidence
from .site_geometry import read_receptor
from .statistics import paired_family_bootstrap


SEEDS = (104729, 130363, 155921)
PROFILES = ("strict", "repair")
CODE_SOURCES = SHARED_CODE_SOURCES + ("aav_benchmark/haddock_pose_report.py", "aav_benchmark/haddock_paired_cohort.py", "aav_benchmark/site_geometry.py")


class OutputEvidence(Evidence):
    """Read-only ordinary/archive view; verified bytes may be materialized in temp."""
    def __init__(self, root):
        super().__init__(root)
        self.output_files = {}

    def read(self, relative, expected=None):
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Evidence paths must remain repository-relative")
        name = path.as_posix()
        if name not in self.output_files or (self.root / path).exists():
            raw = super().read(relative, expected)
            if name in self.output_files and raw != self.output_files[name]:
                raise ValueError("Ordinary and archived native output bytes disagree")
            return raw
        raw = self.output_files[name]
        digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError("Archived native output checksum changed")
        self.sources[name] = {"sha256": digest, "bytes": len(raw)}
        return raw

    def present(self, path):
        return (self.root / path).is_file() or Path(path).as_posix() in self.output_files

    def add_output_archive(self, archive_record, manifest_record, archive_root):
        manifest = json.loads(_pin(self, manifest_record))
        expected = {}
        archive_root = Path(archive_root)
        for entry in manifest["files"]:
            path = Path(entry["path"])
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("Unsafe archived native output identity")
            try:
                relative = path.relative_to(archive_root).as_posix()
            except ValueError as error:
                raise ValueError("Archived native output escaped the declared cohort root") from error
            if relative in expected:
                raise ValueError("Duplicate archived native output identity")
            expected[relative] = entry
        actual = set()
        with tarfile.open(fileobj=io.BytesIO(_pin(self, archive_record)), mode="r:gz") as archive:
            for member in archive.getmembers():
                if member.isdir():
                    continue
                if not member.isfile() or member.name not in expected or member.name in actual:
                    raise ValueError("Unexpected or duplicate archived native model/output")
                raw = archive.extractfile(member).read()
                entry = expected[member.name]
                _source_bytes(raw, entry)
                path = entry["path"]
                if path in self.output_files and self.output_files[path] != raw:
                    raise ValueError("Two native output archives disagree")
                self.output_files[path] = raw
                actual.add(member.name)
        if actual != set(expected):
            raise ValueError("Native output archive is incomplete")

    def matching_files(self, directory, pattern):
        directory = Path(directory)
        names = {p.relative_to(self.root) for p in (self.root / directory).glob(pattern) if p.is_file()}
        names.update(Path(p) for p in self.output_files if Path(p).parent == directory and Path(p).match(pattern))
        return sorted(names)


def _present(evidence, path):
    return evidence.present(path) if isinstance(evidence, OutputEvidence) else (evidence.root / path).is_file()


def _inventory(evidence, run_dir, job_id):
    archived = [p for p in getattr(evidence, "output_files", {}) if Path(p).is_relative_to(run_dir)]
    if not archived:
        return model_inventory(evidence.root / run_dir, job_id)
    names = set(archived)
    names.update(p.relative_to(evidence.root).as_posix() for p in (evidence.root / run_dir).rglob("*") if p.is_file())
    # The frozen inventory consumes top-level stage PDBs/io only. Unrelated
    # local debug files must not create extra public-report dependencies.
    relevant = []
    for name in names:
        relative = Path(name).relative_to(run_dir)
        if len(relative.parts) == 2 and relative.parts[0][0].isdigit() and "_" in relative.parts[0]:
            if relative.name == "io.json" or relative.name.endswith((".pdb", ".pdb.gz")):
                relevant.append(name)
    with tempfile.TemporaryDirectory(prefix="haddock-archived-models-") as temporary:
        temporary = Path(temporary)
        for name in sorted(relevant):
            destination = temporary / Path(name).relative_to(run_dir)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(evidence.read(name))
        return model_inventory(temporary, job_id)


def _pin(evidence, record):
    raw = evidence.read(record["path"], record["sha256"])
    if "bytes" in record and len(raw) != record["bytes"]:
        raise ValueError("Pinned native input/metadata size changed")
    return raw


def _archive(evidence, archive_record, manifest_record):
    """Read all declared source preparation members without installing HADDOCK."""
    manifest_raw = _pin(evidence, manifest_record)
    manifest = json.loads(manifest_raw)
    if manifest.get("source_only") is not True or manifest.get("no_cohort_poses") is not True:
        raise ValueError("Prepared-input archive unexpectedly claims cohort poses")
    raw = _pin(evidence, archive_record)
    expected = {e["path"]: e for e in manifest["files"]}
    if len(expected) != len(manifest["files"]) or any(Path(p).is_absolute() or ".." in Path(p).parts for p in expected):
        raise ValueError("Invalid prepared-input archive member identities")
    result = {}
    embedded_manifest = Path(manifest_record["path"]).name
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        for member in archive.getmembers():
            if member.isdir():
                continue
            if not member.isfile() or member.name in result:
                raise ValueError("Unexpected or duplicate prepared-input archive member")
            content = archive.extractfile(member).read()
            if member.name == embedded_manifest and member.name not in expected:
                if content != manifest_raw:
                    raise ValueError("Embedded prepared-input manifest differs from its external frozen pin")
                result[member.name] = content
                continue
            if member.name not in expected:
                raise ValueError("Unexpected prepared-input archive member")
            entry = expected[member.name]
            if hashlib.sha256(content).hexdigest() != entry["sha256"] or len(content) != entry["bytes"]:
                raise ValueError("Prepared-input archive member changed")
            result[member.name] = content
    observed = set(result)
    if embedded_manifest in observed and embedded_manifest not in expected:
        observed.remove(embedded_manifest)
    if observed != set(expected):
        raise ValueError("Prepared-input archive is incomplete")
    return result


def _prepared(evidence, record, archive, archive_root):
    path = Path(record["path"])
    try:
        relative = path.relative_to(archive_root).as_posix()
    except ValueError as error:
        raise ValueError("Prepared numerical input is outside its declared portable archive") from error
    if relative not in archive:
        raise ValueError("Prepared numerical input is absent from its portable archive")
    raw = archive[relative]
    _source_bytes(raw, record)
    if (evidence.root / path).exists():
        full = evidence.root / path
        if full.is_symlink():
            raise ValueError("Prepared input must be an ordinary file or its explicit archived bytes")
        ordinary = full.read_bytes()
        _source_bytes(ordinary, record)
        if ordinary != raw:
            raise ValueError("Ordinary and archived prepared inputs disagree")
    return raw


def _source_bytes(raw, historical_record):
    """Resolve a historical processed-tree pin to explicit archived source bytes."""
    if hashlib.sha256(raw).hexdigest() != historical_record["sha256"] or ("bytes" in historical_record and len(raw) != historical_record["bytes"]):
        raise ValueError("Portable source does not match its recorded historical input pin")


def parse_model(raw_pdb, cns_atom_map, parent_atom_map, expected_protein_keys):
    """Map complete CCD/CNS ligand atoms without fitting or reading atomic truth."""
    parent_keys = [a["canonical_atom_key"] for a in parent_atom_map]
    if not parent_keys or len(set(parent_keys)) != len(parent_keys) or len({a["name"] for a in parent_atom_map}) != len(parent_keys):
        raise ValueError("Full-parent canonical map is incomplete or duplicated")
    cns = {}
    by_canonical = {}
    for atom in cns_atom_map:
        canonical = atom["canonical_atom_key"]
        key = (str(atom["cns_chain"]), int(atom["cns_resnum"]), "", str(atom["cns_resname"]), str(atom["cns_atom_name"]))
        if key in cns or canonical in by_canonical:
            raise ValueError("CNS/canonical full-parent map is not a bijection")
        if key[0] != "L":
            raise ValueError("Declared ligand CNS segment must be L")
        cns[key], by_canonical[canonical] = atom, key
    if set(by_canonical) != set(parent_keys):
        raise ValueError("CNS ligand map lost a full-parent atom; scored-core truncation is forbidden")
    expected_protein_keys = [tuple(k) for k in expected_protein_keys]
    if len(set(expected_protein_keys)) != len(expected_protein_keys):
        raise ValueError("Input receptor heavy identities are duplicated")
    all_atoms, ligand, protein, elements = set(), {}, {}, {}
    model_count = 0
    total = 0
    for line in raw_pdb.decode().splitlines():
        if line.startswith("MODEL"):
            model_count += 1
            if model_count > 1:
                raise ValueError("An individual HADDOCK artifact contains multiple coordinate models")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain, segment = line[21:22].strip(), line[72:76].strip()
        if chain and segment and chain != segment:
            raise ValueError("PDB chain and CNS segment identities disagree")
        identity = chain or segment
        if len(identity) != 1:
            raise ValueError("Model lacks an unambiguous single-character chain/segment identity")
        if line[16:17].strip():
            raise ValueError("Model has unresolved alternate conformations")
        key = (identity, int(line[22:26]), line[26:27].strip(), line[17:20].strip(), line[12:16].strip())
        if key in all_atoms:
            raise ValueError("Model has duplicated atom identities")
        all_atoms.add(key)
        xyz = np.asarray([float(line[start:start + 8]) for start in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError("Model has nonfinite atom coordinates")
        element = line[76:78].strip().upper() or next((c.upper() for c in key[-1] if c.isalpha()), "")
        if not element:
            raise ValueError("Model atom element is absent")
        total += 1
        if element in {"H", "D"}:
            continue
        if identity == "L":
            if key not in cns or element != cns[key]["element"].upper():
                raise ValueError("Model ligand changed a complete canonical heavy identity/element")
            ligand[key] = xyz
        else:
            protein[key], elements[key] = xyz, element
    if set(ligand) != set(cns) or set(protein) != set(expected_protein_keys):
        raise ValueError("Model lost or added full-parent/protein heavy atoms")
    return {"ligand_coordinates": np.asarray([ligand[by_canonical[k]] for k in parent_keys]),
            "protein_keys": list(protein), "protein_coordinates": np.asarray(list(protein.values())),
            "protein_elements": [elements[k] for k in protein], "total_atom_count": total}


def _inventory_shape(stages, max_models):
    if isinstance(max_models, bool) or not isinstance(max_models, int) or max_models < 1:
        raise ValueError("Model cap must be a positive integer")
    names, ids, selected = set(), set(), {}
    for stage in stages:
        name = stage["stage"]
        if name in names or not any(name.endswith("_" + suffix) for suffix in ("topoaa", "rigidbody", "flexref")):
            raise ValueError("Duplicate or undeclared native stage; no clustering/native accuracy stage is allowed")
        names.add(name)
        models = stage["models"]
        if stage["actual_physical_model_files"] != len(models):
            raise ValueError("Native physical-model count differs from retained inventory")
        for model in models:
            filename = model["file_name"]
            if Path(filename).name != filename or not filename.removesuffix(".gz").endswith(".pdb"):
                raise ValueError("Unsafe native model artifact name")
            if model["model_id"] in ids or not model["model_id"].endswith("/" + name + "/" + filename.removesuffix(".gz")):
                raise ValueError("Duplicated or inconsistent immutable model identity")
            ids.add(model["model_id"])
        for suffix in ("rigidbody", "flexref"):
            if name.endswith("_" + suffix):
                if suffix in selected or len(models) > max_models:
                    raise ValueError("Duplicate native stage or physical-model budget exceeded")
                selected[suffix] = stage
    return selected


def verify_stage_inventory(stages, max_models=20):
    """Check individual score ranks at full precision, without native evaluation."""
    selected = _inventory_shape(stages, max_models)
    if set(selected) != {"rigidbody", "flexref"}:
        raise ValueError("Completed native workflow lacks one initial and one final stage")
    for stage in selected.values():
        models = stage["models"]
        if not models:
            raise ValueError("Completed native stage returned no physical models")
        for model in models:
            score = model.get("score")
            if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
                raise ValueError("A retained physical model lacks finite full-precision native score metadata")
        ranked = sorted(models, key=lambda row: (row["score"], row["original_model_number"] or 0, row["model_id"]))
        if stage["individual_score_order"] != [m["model_id"] for m in ranked]:
            raise ValueError("Individual native scores/model-identity tie order changed")
        if any(m.get("individual_score_rank") != rank for rank, m in enumerate(ranked, 1)):
            raise ValueError("Native individual score rank changed")
    if len(selected["flexref"]["models"]) > len(selected["rigidbody"]["models"]):
        raise ValueError("Refinement generated undeclared replacement/model slots")
    return copy.deepcopy(selected)


def _verify_lineage(evidence, run_dir, selected):
    initial, final = selected["rigidbody"], selected["flexref"]
    initial_by_name = {m["file_name"].removesuffix(".gz"): m for m in initial["models"]}
    metadata = evidence.json(run_dir / final["stage"] / "io.json", final["io_sha256"])
    inputs = metadata.get("input")
    if not isinstance(inputs, list) or any(not isinstance(value, dict) or not isinstance(value.get("file_name"), str) for value in inputs):
        raise ValueError("Final stage lacks explicit native input-model lineage")
    input_names = [Path(value["file_name"]).name.removesuffix(".gz") for value in inputs]
    if len(input_names) != len(initial_by_name) or len(set(input_names)) != len(input_names) or set(input_names) != set(initial_by_name):
        raise ValueError("Not every initial physical model was passed once to refinement")
    outputs = metadata.get("output", [])
    if not isinstance(outputs, list) or any(not isinstance(value, dict) or not isinstance(value.get("file_name"), str) for value in outputs):
        raise ValueError("Final stage lacks explicit native output-model lineage")
    output_names = [value["file_name"] for value in outputs]
    physical_names = {model["file_name"].removesuffix(".gz") for model in final["models"]}
    if len(output_names) != len(set(output_names)) or set(output_names) != physical_names:
        raise ValueError("Native output score metadata duplicated or omitted a retained physical model")
    output_by_name = {value["file_name"]: value for value in outputs}
    parents = []
    for model in final["models"]:
        name = model["file_name"].removesuffix(".gz")
        native = output_by_name.get(name)
        if not native or not isinstance(native.get("ori_name"), str):
            raise ValueError("Retained refined model lacks its native original-model identity")
        parent = Path(native["ori_name"]).name.removesuffix(".gz")
        if parent not in initial_by_name or parent in parents:
            raise ValueError("Refinement replaced, duplicated or invented an original model slot")
        if model.get("actual_CNS_seed") != initial_by_name[parent].get("actual_CNS_seed"):
            raise ValueError("Refined model lost its actual native seed lineage")
        parents.append(parent)
    return {"initial_models_passed_to_refinement": len(input_names), "refined_parent_models": parents,
            "refinement_failed_slots_not_replaced": len(parents) <= len(input_names)}


def read_job(evidence, job, condition, prepared_plan, prepared_pin, freeze_pin, *, recovery=None, execution_receipt=None):
    """Verify one declared native job from retained artifacts, never software."""
    directory = Path(job["directory"])
    execution_path = Path(execution_receipt["path"]) if execution_receipt is not None else directory / "execution.json"
    if execution_path not in {directory / "execution.json", directory / "wrapper_failure.json"}:
        raise ValueError("Native execution receipt is outside its declared job directory")
    if not _present(evidence, execution_path):
        return None, "missing_native_execution"
    execution = (json.loads(_pin(evidence, execution_receipt)) if execution_receipt is not None else evidence.json(execution_path))
    if execution.get("status") == "orchestration_failed":
        if (execution_path != directory / "wrapper_failure.json" or recovery is None or execution.get("terminal") is not True
                or execution.get("eligible_workflow_pool") is not True or execution.get("job_id") != job["job_id"]
                or execution.get("prepared_plan", {}).get("sha256") != prepared_pin["sha256"] or execution.get("config") != job["config"]
                or execution.get("process_return_code") is not None or execution.get("native_retry") is not False
                or execution.get("no_accuracy_evaluation") is not True):
            raise ValueError("Orchestration failure lacks its explicit frozen-job/unknown-return-code provenance")
        partial = _inventory(evidence, directory / "run", job["job_id"])
        if partial != execution.get("retained_partial_stages"):
            raise ValueError("Retained orchestration-failure partial artifacts changed")
        selected = _inventory_shape(partial, 20)
        for stage in partial:
            stage_dir = directory / "run" / stage["stage"]
            if stage.get("io_sha256") is not None:
                evidence.read(stage_dir / "io.json", stage["io_sha256"])
            for model in stage["models"]:
                evidence.read(stage_dir / model["file_name"], model["sha256"])
        return {"status": "orchestration_failed", "execution": execution, "final_stage": None, "lineage": None,
                "retained_initial_models": selected.get("rigidbody", {}).get("actual_physical_model_files", 0),
                "retained_refined_models": selected.get("flexref", {}).get("actual_physical_model_files", 0)}, None
    if (execution.get("job_id") != job["job_id"] or execution.get("terminal") is not True
            or execution.get("eligible_workflow_pool") is not True or execution.get("prepared_plan", {}).get("sha256") != prepared_pin["sha256"]
            or execution.get("root_freeze", {}).get("sha256") != freeze_pin["sha256"] or execution.get("helper_source") != prepared_plan["driver"]
            or execution.get("CNS_binary") != prepared_plan["CNS_binary"] or execution.get("config") != job["config"]
            or execution.get("ncores") != 1 or execution.get("no_accuracy_evaluation") is not True
            or execution.get("failed_model_slots_not_replaced") is not True or execution.get("clustering_performed") is not False
            or execution.get("requested_initial_model_slots") != 20):
        raise ValueError("Native execution changed its declared job/configuration/provenance")
    recovered = execution.get("status") == "complete_outputs_recovered"
    if execution.get("status") not in {"completed", "native_execution_failed", "complete_outputs_recovered"}:
        raise ValueError("Unknown terminal native HADDOCK status")
    if (not isinstance(execution.get("elapsed_seconds"), (int, float)) or not math.isfinite(execution["elapsed_seconds"])
            or execution["elapsed_seconds"] < 0 or (not recovered and not isinstance(execution.get("process_return_code"), int))):
        raise ValueError("Native execution lacks its actual exit/resource measurements")
    if recovered:
        if (recovery is None or execution.get("recovery_receipt") != recovery["adoption_pin"]
                or execution.get("job_id") != recovery["adoption"]["job_id"] or execution.get("process_return_code") is not None
                or execution.get("native_rerun") is not False or execution.get("process_CPU_seconds") is not None
                or execution.get("external_measured_wall_seconds") is not None
                or execution.get("elapsed_seconds") != recovery["adoption"]["native_self_reported_elapsed_seconds"]):
            raise ValueError("Recovered output acceptance must retain explicit unknown return-code/CPU telemetry")
    else:
        command = execution.get("command", [])
        if len(command) != 2 or command[1] != "run.cfg" or Path(command[0]).name != "haddock3":
            raise ValueError("Native execution changed its declared module workflow")
    run_dir = directory / "run"
    replay = _inventory(evidence, run_dir, job["job_id"])
    if replay != execution.get("stages"):
        raise ValueError("Retained physical models/io metadata differ from native execution inventory")
    selected = _inventory_shape(replay, 20)
    for stage in replay:
        stage_dir = run_dir / stage["stage"]
        if stage.get("io_sha256") is not None:
            evidence.read(stage_dir / "io.json", stage["io_sha256"])
        for model in stage["models"]:
            evidence.read(stage_dir / model["file_name"], model["sha256"])
        # Effective weights/configurations and actual CNS topologies are raw
        # source evidence; no offline native-score formula is substituted.
        for pattern in ("params.cfg", "*.psf"):
            paths = (evidence.matching_files(stage_dir, pattern) if isinstance(evidence, OutputEvidence)
                     else [path.relative_to(evidence.root) for path in sorted((evidence.root / stage_dir).glob(pattern))])
            for path in paths:
                evidence.read(path)
    for name in ("stdout.log", "stderr.log"):
        evidence.read(directory / name)
    initial_count = selected.get("rigidbody", {}).get("actual_physical_model_files", 0)
    final_count = selected.get("flexref", {}).get("actual_physical_model_files", 0)
    if (execution.get("actual_initial_models") != initial_count or execution.get("actual_refined_models") != final_count
            or execution.get("unfilled_initial_slots") != 20 - initial_count):
        raise ValueError("Native physical model budgets differ from execution receipt")
    if recovered and (initial_count != 20 or final_count != 20):
        raise ValueError("Only the separately approved complete20/20 orphan may use recovered-output acceptance")
    final = None
    lineage = None
    if execution["status"] in {"completed", "complete_outputs_recovered"}:
        if (not recovered and execution["process_return_code"] != 0) or execution.get("final_score_metadata_complete") is not True:
            raise ValueError("Completed native execution did not pass its actual final model gate")
        strict = verify_stage_inventory(replay)
        final = strict["flexref"]
        for stage in strict.values():
            evidence.read(run_dir / stage["stage"] / "params.cfg")
        topology_stages = [stage for stage in replay if stage["stage"].endswith("_topoaa")]
        if len(topology_stages) != 1:
            raise ValueError("Completed native workflow lacks its retained topology stage")
        for name in ("receptor_haddock.psf", "glycan_haddock.psf"):
            evidence.read(run_dir / topology_stages[0]["stage"] / name)
        lineage = _verify_lineage(evidence, run_dir, strict)
    return {"status": execution["status"], "execution": execution, "final_stage": final,
            "run_dir": run_dir, "lineage": lineage}, None


def _load_recovery(evidence, evaluation, prepared):
    records = evaluation.get("orchestration_recovery")
    if records is None:
        return None
    recovery = json.loads(_pin(evidence, records["plan"]))
    authorization = json.loads(_pin(evidence, records["freeze"]))
    adoption = json.loads(_pin(evidence, records["adoption_receipt"]))
    if "preparation_receipt" in records:
        archived = json.loads(_pin(evidence, records["preparation_receipt"]))
        if (archived["frozen_original_plan_sha256"] != evaluation["prepared_plan"]["sha256"]
                or archived["native_first_job_rerun"] is not False or archived["original_scientific_inputs_changed"] is not False):
            raise ValueError("Recovery raw archive changed the original experiment")
        if not isinstance(evidence, OutputEvidence):
            raise ValueError("Archived recovered outputs require the read-only report evidence view")
        evidence.add_output_archive(archived["raw_archive"], archived["raw_manifest"], Path(evaluation["prepared_plan"]["path"]).parent)
    if (authorization.get("status") != "frozen_before_HADDOCK_orchestration_resumption" or authorization.get("execution_authorized") is not True
            or authorization.get("recovery_plan", {}).get("sha256") != records["plan"]["sha256"]
            or recovery.get("original_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]
            or recovery.get("original_execution_freeze", {}).get("sha256") != evaluation["native_execution_freeze"]["sha256"]
            or recovery.get("frozen_helper") != prepared["driver"] or recovery.get("native_model_parameters_changed") is not False
            or recovery.get("no_native_rerun") is not True or recovery.get("no_accuracy_evaluation") is not True):
        raise ValueError("Orchestration recovery changed the original experiment or lacks a separate root freeze")
    _pin(evidence, recovery["recovery_wrapper"])
    proof = json.loads(_pin(evidence, recovery["orphan_completion_proof"]))
    if (adoption.get("status") != "verified_completed_outputs_adopted_without_rerun" or adoption.get("native_rerun") is not False
            or adoption.get("original_native_return_code") is not None or adoption.get("original_process_CPU_seconds") is not None
            or adoption.get("no_accuracy_evaluation") is not True or adoption.get("job_id") != recovery["orphan_job_id"]
            or adoption.get("completion_proof") != recovery["orphan_completion_proof"]
            or adoption.get("recovery_plan", {}).get("sha256") != records["plan"]["sha256"]
            or adoption.get("recovery_freeze", {}).get("sha256") != records["freeze"]["sha256"]
            or proof.get("job_id") != adoption["job_id"] or proof.get("no_native_rerun") is not True
            or proof.get("no_accuracy_evaluation") is not True or proof.get("process_return_code") is not None
            or proof.get("process_CPU_seconds") is not None or proof.get("external_measured_wall_seconds") is not None
            or proof.get("frozen_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]
            or proof.get("original_execution_freeze", {}).get("sha256") != evaluation["native_execution_freeze"]["sha256"]
            or proof.get("frozen_helper") != prepared["driver"]
            or proof.get("native_completion", {}).get("normal_finished_marker") is not True
            or proof.get("native_completion", {}).get("elapsed_seconds") != adoption.get("native_self_reported_elapsed_seconds")):
        raise ValueError("Recovered outputs lack their exact independently approved complete-budget proof")
    for record in proof["raw_artifacts"]:
        _pin(evidence, record)
    _pin(evidence, adoption["original_driver_failure"])
    return {"plan_pin": records["plan"], "freeze_pin": records["freeze"], "adoption_pin": records["adoption_receipt"],
            "adoption": adoption, "proof": proof, "plan": recovery}


def _load_source_protocol(evidence, evaluation):
    """Validate portable pre-execution source evidence, not local software bytes."""
    prepared = json.loads(_pin(evidence, evaluation["prepared_plan"]))
    receipt = json.loads(_pin(evidence, evaluation["preparation_receipt"]))
    freeze = json.loads(_pin(evidence, evaluation["native_execution_freeze"]))
    if (freeze.get("status") != "frozen_before_first_HADDOCK_cohort_pose" or freeze.get("execution_authorized") is not True
            or freeze.get("prepared_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]
            or receipt.get("prepared_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]
            or receipt.get("prepared_plan", {}).get("path") != evaluation["prepared_plan"]["path"]):
        raise ValueError("Portable native provenance does not bind the independently frozen prepared plan")
    if (prepared.get("source_only_preparation") is not True or prepared.get("native_cohort_docking_started") is not False
            or prepared.get("no_accuracy_evaluation") is not True or prepared.get("ncores") != 1
            or prepared.get("condition_denominator") != 36 or prepared.get("seed_cell_denominator") != 108
            or prepared.get("seeds") != list(SEEDS) or prepared.get("native_modules") != ["topoaa", "rigidbody", "flexref"]
            or receipt.get("clustering") is not False or receipt.get("native_reference_or_interface_restraints") is not False):
        raise ValueError("Frozen guided workflow/source/CPU/denominator criteria changed")
    for key in ("proposal", "driver", "assembly_protocol_policy", "source_manifest", "preflight_receipt", "release_lock"):
        _pin(evidence, prepared[key])
    software = json.loads(_pin(evidence, receipt["native_software_manifest"]))
    weights = json.loads(_pin(evidence, receipt["native_score_weights"]))
    if (software.get("installed_source_files") != prepared["runtime_sources"]
            or software.get("CNS_binary") != prepared["CNS_binary"]
            or software.get("release_distribution") != prepared["release_distribution"]):
        raise ValueError("Historical native software provenance differs from the execution freeze")
    # Deliberately never open runtime_sources/CNS_binary/release_distribution/
    # native_environment_freeze paths: these are historical installation pins.
    archive = _archive(evidence, receipt["archive"], receipt["archive_manifest"])
    if len(archive) != receipt["archived_files"]:
        raise ValueError("Portable prepared-input archive denominator changed")
    root = Path(evaluation["prepared_plan"]["path"]).parent
    conditions, jobs = prepared["conditions"], prepared["jobs"]
    expected_conditions = {(p, cid, method) for p in PROFILES for cid in CASE_IDS for method in METHODS}
    if (len(conditions) != 36 or {(c["profile"], c["case_id"], c["method"]) for c in conditions} != expected_conditions
            or len({c["condition_id"] for c in conditions}) != 36):
        raise ValueError("Every six-case strict/repair/site condition remains eligible")
    expected_jobs = {(p, cid, method, seed) for p, cid, method in expected_conditions for seed in SEEDS}
    if (len(jobs) != 108 or {(j["profile"], j["case_id"], j["method"], j["seed"]) for j in jobs} != expected_jobs
            or len({j["job_id"] for j in jobs}) != 108):
        raise ValueError("Every condition/seed cell remains eligible")
    ready = sum(job["status"] == "ready" for job in jobs)
    if ready != receipt["ready_native_jobs"] or sum(c["status"] == "ready" for c in conditions) != receipt["ready_conditions"]:
        raise ValueError("Source gate changed the frozen native-search coverage")
    for condition in conditions:
        if condition.get("eligible_workflow_condition") is not True:
            raise ValueError("Source-failed guided condition disappeared from eligibility")
        _pin(evidence, condition["source_condition"])
        if condition["status"] == "ready":
            for record in condition["input_files"]:
                _prepared(evidence, record, archive, root)
    for job in jobs:
        if job.get("eligible_workflow_pool") is not True or job.get("initial_model_slots") != 20 or job.get("maximum_refined_model_slots") != 20:
            raise ValueError("Native guided per-seed model budget/eligibility changed")
        identity = f"{job['profile']}/{job['case_id']}/{job['method']}/box01"
        matching = [c for c in conditions if c["condition_id"] == identity]
        if (len(matching) != 1 or job["condition_id"] != identity or job["job_id"] != identity + "/s" + str(job["seed"])
                or job["status"] != matching[0]["status"] or job.get("model_seed_range") != [job["seed"] + 1, job["seed"] + 20]):
            raise ValueError("Native seed job is detached from its fixed case/profile/site identity")
        if job["status"] == "ready":
            _prepared(evidence, job["config"], archive, root)
    return prepared, receipt, archive, root, {"native_software_metadata": receipt["native_software_manifest"],
        "recorded_installed_source_files": len(prepared["runtime_sources"]), "CNS_binary_sha256": prepared["CNS_binary"]["sha256"],
        "local_installation_or_binary_reverified": False, "provenance_scope": "Historical native execution freeze; portable metric replay requires numerical models/inputs, not a local HADDOCK installation", "native_score_weights": weights,
        "source_parameterization_gate": "Pinned historical CNS preflight judgment; native parameterization is not rerun by this external metric report",
        "score_policy": "Exact retained stage io.json scores and params.cfg remain authoritative; no offline energy formula or cross-method score comparison"}


def read_workflows(evidence, evaluation, prepared, *, recovery=None):
    """All108 progress records and ready-job artifacts precede any accuracy."""
    progress_path = Path(evaluation["execution_progress"])
    jobs = {j["job_id"]: j for j in prepared["jobs"]}
    conditions = {c["condition_id"]: c for c in prepared["conditions"]}
    if not _present(evidence, progress_path):
        return {}, [job + ":missing_execution_progress" for job in jobs], {"eligible_cells": 108, "terminal_cells": 0, "native_ready_cells": sum(j["status"] == "ready" for j in jobs.values())}
    progress = evidence.json(progress_path)
    if (progress.get("prepared_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]
            or progress.get("root_freeze", {}).get("sha256") != evaluation["native_execution_freeze"]["sha256"]
            or progress.get("seed_cell_denominator") != 108 or progress.get("no_accuracy_evaluation") is not True):
        raise ValueError("Native progress index changed its source/configuration denominator")
    if progress.get("schema_version") == 2:
        if (recovery is None or progress.get("recovery_plan", {}).get("sha256") != recovery["plan_pin"]["sha256"]
                or progress.get("execution_recovery_freeze", {}).get("sha256") != recovery["freeze_pin"]["sha256"]):
            raise ValueError("Recovered orchestration progress lacks its independently frozen source provenance")
    results = progress.get("jobs", [])
    by_id = {r["job_id"]: r for r in results}
    if len(by_id) != len(results) or set(by_id) - set(jobs) or progress.get("terminal_cells") != len(results):
        raise ValueError("Native progress index duplicated or invented an eligible cell")
    pending, checked = [], {}
    counts = {"eligible_cells": 108, "terminal_cells": 0, "native_ready_cells": sum(j["status"] == "ready" for j in jobs.values()),
              "completed_native_cells": 0, "recovered_output_cells": 0, "orchestration_failed_cells": 0,
              "failed_native_cells": 0, "source_failed_cells": 0, "retained_initial_models": 0, "retained_refined_models": 0}
    for job_id, job in jobs.items():
        if job_id not in by_id:
            pending.append(job_id + ":missing_terminal_cell")
            continue
        record = by_id[job_id]
        if record.get("terminal") is not True or record.get("eligible_workflow_pool") is not True:
            pending.append(job_id + ":not_terminal")
            continue
        if job["status"] != "ready":
            if record.get("status") != job["status"] or record.get("physical_models") != 0 or record.get("no_native_execution") is not True:
                raise ValueError("Unsupported source/assembly cell lost its explicit unexecuted failure mechanism")
            checked[job_id] = {"status": record["status"], "execution": record, "final_stage": None, "lineage": None}
            counts["source_failed_cells"] += 1
        else:
            result, unfinished = read_job(evidence, job, conditions[job["condition_id"]], prepared, evaluation["prepared_plan"], evaluation["native_execution_freeze"],
                recovery=recovery, execution_receipt=record.get("execution_receipt"))
            if unfinished:
                pending.append(job_id + ":" + unfinished)
                continue
            core_record = {key: value for key, value in record.items() if key not in {"execution_receipt", "resource_usage"}}
            if result["execution"] != core_record:
                raise ValueError("Native condition progress and individual execution receipts disagree")
            if "resource_usage" in record:
                resource = json.loads(_pin(evidence, record["resource_usage"]))
                if (resource.get("job_id") != job_id or resource.get("native_cpus_per_job") != 1
                        or resource.get("original_plan", {}).get("sha256") != evaluation["prepared_plan"]["sha256"]):
                    raise ValueError("Measured resource sidecar changed job/CPU/configuration identity")
                for field in ("CPU_seconds", "wall_seconds"):
                    value = resource.get("measured_usage", {}).get(field)
                    if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                        raise ValueError("Resource sidecar lacks finite actual measurements")
                result["resource_measurement"] = resource
            checked[job_id] = result
            count_field = {"completed": "completed_native_cells", "complete_outputs_recovered": "recovered_output_cells", "orchestration_failed": "orchestration_failed_cells", "native_execution_failed": "failed_native_cells"}[result["status"]]
            counts[count_field] += 1
            counts["retained_initial_models"] += record.get("actual_initial_models", result.get("retained_initial_models", 0))
            counts["retained_refined_models"] += record.get("actual_refined_models", result.get("retained_refined_models", 0))
        counts["terminal_cells"] += 1
    if not pending and (progress.get("all_terminal") is not True or counts["terminal_cells"] != 108):
        raise ValueError("Native execution does not satisfy the complete frozen workflow gate")
    return checked, sorted(pending), counts


def _condition_map(evidence, condition, source, archive, root):
    """Bind independent aliases/protein identities to portable preflight inputs."""
    cid = condition["case_id"]
    _source_bytes(source[f"{cid}/chemistry.json"], condition["chemistry"])
    if condition["status"] != "ready":
        return None
    records = {Path(r["path"]).name: r for r in condition["input_files"]}
    preparation = json.loads(_prepared(evidence, records["preflight_preparation.json"], archive, root))
    if (preparation.get("source_only") is not True or preparation.get("native_coordinates_read") is not False
            or preparation.get("no_docking_stages") is not True or preparation["full_parent_chemistry"].get("experimental_coordinates_read") is not False):
        raise ValueError("Prepared glycan/topology source contract accessed atomic reference truth")
    aliases = {"original_independent_sdf": "ligand.sdf", "observed_protein_atom_map": "protein_atom_map.json",
               "receptor_source": "receptor.pdb", "source_report": "initializer_report.json"}
    for key, name in aliases.items():
        _source_bytes(source[f"{cid}/{name}"], preparation[key])
    if "modeled_repair_atom_map_not_truth" in preparation:
        _source_bytes(source[f"{cid}/modeled_atom_map.json"], preparation["modeled_repair_atom_map_not_truth"])
    _source_bytes(_prepared(evidence, records["glycan_translated.sdf"], archive, root), preparation["input_sdf"])
    _source_bytes(_prepared(evidence, records["glycan.pdb"], archive, root), preparation["input_PDB"])
    air = json.loads(_prepared(evidence, records["box_air.json"], archive, root))
    if (air.get("native_contacts_used") is not False or air.get("residue_expansion") is not False
            or air.get("random_removal") is not False or air.get("assignment_count") != 1 or air["box"] != condition["box"]):
        raise ValueError("The declared coarse box-to-AIR guidance changed")
    parent = json.loads(source[f"{cid}/ligand_atom_map.json"])["atoms"]
    if ({a["canonical_atom_key"] for a in preparation["atom_map"]} != {a["canonical_atom_key"] for a in parent}
            or len(preparation["atom_map"]) != len(parent)):
        raise ValueError("Prepared CNS aliases changed the complete independent parent graph")
    return preparation["atom_map"]


def _statistics(rows, case_to_family):
    result = []
    for profile in PROFILES:
        for scope in ("full_ligand", "modeled_core"):
            selected = [r for r in rows if r["profile"] == profile and r["atomic_scope"] == scope]
            for alignment in ("global", "local_interface"):
                for rank in ("top1", "top5", "best_returned"):
                    for threshold in (1, 2, 3):
                        records = [{"family_id": case_to_family[r["case_id"]], "complex_id": r["case_id"], "method": r["method"],
                            "seed": seed["seed"], "status": "ok" if seed["poses"] else "failed",
                            "success": seed[alignment][f"{rank}_success_{threshold}A"]} for r in selected for seed in r["seeds"]]
                        for method in ("fpocket", "p2rank"):
                            summary = paired_family_bootstrap(records, "oracle", method, n_resamples=10000, seed=20261103, expected_seeds=SEEDS)
                            summary.update(profile=profile, atomic_scope=scope, alignment=alignment, pose_rank=rank, threshold_A=threshold,
                                information_regime="guided coarse box-to-AIR patch contact, not a hard Cartesian box", executed_candidate_budget=1,
                                interpretation="Conditional source-family resampling of this purposive structural sample; no population/model-training-unseen reliability claim")
                            result.append(summary)
    return result


def generate_report(root, plan_path, freeze_path, *, allow_pending=False):
    """Read source/terminal evidence first; external accuracy is strictly second."""
    evidence = OutputEvidence(root)
    evaluation, evaluation_freeze = evidence.json(plan_path), evidence.json(freeze_path)
    if (evaluation_freeze.get("evaluation_plan_sha256") != sha256(evidence.root / plan_path)
            or evaluation_freeze.get("frozen_before_test_pose_accuracy") is not True
            or evaluation_freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("External guided evaluation requires its independent pre-accuracy freeze")
    if (evaluation.get("method_id") != "haddock3_guided_box_air" or evaluation.get("primary_stage") != "flexref"
            or evaluation.get("case_ids") != list(CASE_IDS) or evaluation.get("seeds") != list(SEEDS)
            or evaluation.get("rank_budgets") != [1, 5] or evaluation.get("pose_thresholds_A") != [1, 2, 3]
            or evaluation.get("clustering_used_for_primary") is not False or evaluation.get("no_equivalent_hard_box_information_claim") is not True):
        raise ValueError("Prospective guided-method metric/information design changed")
    if set(evaluation["metric_code_sha256"]) != set(CODE_SOURCES):
        raise ValueError("Every external reader/shared metric implementation must be source-pinned")
    for path, digest in evaluation["metric_code_sha256"].items():
        evidence.read(path, digest)
    reference = json.loads(_pin(evidence, evaluation["reference_evaluation_plan"]))
    reference_freeze = json.loads(_pin(evidence, evaluation["reference_evaluation_freeze"]))
    if (reference_freeze.get("evaluation_plan_sha256") != evaluation["reference_evaluation_plan"]["sha256"]
            or reference_freeze.get("frozen_before_test_pose_accuracy") is not True or reference_freeze.get("new_test_pose_accuracy_consulted") is not False):
        raise ValueError("Shared qualification/metric plan is not bound to its pre-accuracy freeze")
    cases, source_validation = check_public_reference_sources(evidence, reference["source_ledger"], reference["reference_source_manifest"], evaluation["public_source_distribution"])
    _pin(evidence, reference["family_audit"])
    prediction = reference["site_prediction"]
    prediction_config = json.loads(_pin(evidence, prediction["configuration"]))
    prediction_freeze = json.loads(_pin(evidence, prediction["freeze_receipt"]))
    if (prediction_freeze.get("configuration_sha256") != prediction["configuration"]["sha256"]
            or prediction_freeze.get("frozen_before_new_predictions") is not True):
        raise ValueError("Coarse guided patches lack their original pre-prediction input freeze")
    prediction_sources = portable_archive(evidence, prediction["source_archive_manifest"], archive_path=prediction_config["source_preparation_archive"]["path"])
    profile_records = {p["id"]: p for p in reference["profiles"]}
    sources = {profile: portable_archive(evidence, profile_records[evaluation["source_profiles"][profile]]["source_archive_manifest"]) for profile in PROFILES}
    prepared, receipt, archive, archive_root, provenance = _load_source_protocol(evidence, evaluation)
    aliases = {}
    for condition in prepared["conditions"]:
        profile, cid = condition["profile"], condition["case_id"]
        source_condition, candidates = common_selection(evidence,
            {"source_condition": condition["source_condition"], "case_id": cid, "method": condition["method"], "receptor_profile": profile},
            prediction["configuration"]["sha256"], prediction_sources)
        expected_center = candidates[0]["center_A"] if candidates else None
        if ((condition["box"] is None) != (expected_center is None)
                or (condition["box"] is not None and condition["box"]["center"] != expected_center)):
            raise ValueError("Guided patch changed its predeclared original top1/oracle center")
        del source_condition, candidates  # Metadata-only identity check; no site accuracy.
        aliases[condition["condition_id"]] = _condition_map(evidence, condition, sources[profile], archive, archive_root)
    recovery = _load_recovery(evidence, evaluation, prepared)
    for output in evaluation.get("native_output_archives", []):
        evidence.add_output_archive(output["archive"], output["manifest"], Path(evaluation["prepared_plan"]["path"]).parent)
    checked, pending, counts = read_workflows(evidence, evaluation, prepared, recovery=recovery)
    if pending:
        if not allow_pending:
            raise ValueError(f"Guided native execution remains incomplete ({len(pending)} pending cells); no accuracy report")
        return {"schema_version": 1, "status": "pending_incomplete_workflow", "metrics_eligible": False, "counts": counts,
            "pending": pending, "reference_validation": source_validation, "native_provenance": provenance,
            "source_sha256_inventory": evidence.sources}
    if counts["terminal_cells"] != 108:
        raise ValueError("All108 eligible guided cells must be terminal before native accuracy")
    contexts = {(profile, cid): _source_context(sources[profile], cid, cases[cid]) for profile in PROFILES for cid in CASE_IDS}
    rows = []
    jobs = {(j["condition_id"], j["seed"]): j for j in prepared["jobs"]}
    for condition in prepared["conditions"]:
        profile, cid = condition["profile"], condition["case_id"]
        mapping, graph, input_keys, _, protein_map = contexts[profile, cid]
        parent = json.loads(sources[profile][f"{cid}/ligand_atom_map.json"])["atoms"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "receptor.cif"; path.write_bytes(sources[profile][f"{cid}/receptor.cif"])
            keys, _, elements = read_receptor(path)
        expected_elements = dict(zip(keys, elements))
        seeds = []
        for seed in SEEDS:
            job = jobs[condition["condition_id"], seed]
            result = checked[job["job_id"]]
            poses = []
            model_controls = []
            if result["final_stage"] is not None:
                final = result["final_stage"]
                models = {m["model_id"]: m for m in final["models"]}
                for rank, identity in enumerate(final["individual_score_order"], 1):
                    model = models[identity]
                    raw = evidence.read(result["run_dir"] / final["stage"] / model["file_name"], model["sha256"])
                    if model["file_name"].endswith(".gz"):
                        raw = gzip.decompress(raw)
                    coordinates = parse_model(raw, aliases[condition["condition_id"]], parent, input_keys)
                    if coordinates["total_atom_count"] != model["atom_count"] or any(expected_elements[k].upper() != element for k, element in zip(coordinates["protein_keys"], coordinates["protein_elements"])):
                        raise ValueError("Native model atom count/protein element identity changed")
                    # HADDOCK moves/refines protein atoms. Each ligand is fitted
                    # using its OWN corresponding protein, never the input PDB.
                    scored, controls = score_pool([{"rank": rank, "vina_score": model["score"], "coordinates": coordinates["ligand_coordinates"],
                        "candidate_class_rank": 1, "candidate_original_rank": 1, "box_id": "top1_guided_patch", "pooled_rank": rank}],
                        mapping, coordinates["protein_keys"], coordinates["protein_coordinates"], protein_map, graph)
                    value = scored[0]
                    value["native_haddock_score"] = value.pop("vina_score")
                    value.update(model_id=identity, original_model_number=model["original_model_number"], actual_CNS_seed=model["actual_CNS_seed"],
                                 native_energy_terms=model["energies"], score_scope="Individual native full-precision HADDOCK score; no cross-method calibration")
                    poses.append(value)
                    model_controls.append({"model_id": identity, "protein_fits_and_native_controls": controls})
            seeds.append({"seed": seed, "workflow_status": result["status"], "poses": poses, "model_reference_controls": model_controls,
                          "global": accuracy_summary(poses, "global_RMSD_A"), "local_interface": accuracy_summary(poses, "local_interface_RMSD_A"),
                          "native_model_lineage": result["lineage"], "physical_initial_models": result["execution"].get("actual_initial_models", result.get("retained_initial_models", 0)),
                          "physical_refined_models": result["execution"].get("actual_refined_models", result.get("retained_refined_models", 0)),
                          "elapsed_seconds": result["execution"].get("elapsed_seconds"), "process_return_code": result["execution"].get("process_return_code"),
                          "process_return_code_provenance": result["execution"].get("process_return_code_provenance"),
                          "elapsed_seconds_provenance": result["execution"].get("elapsed_seconds_provenance"),
                          "process_CPU_seconds": result["execution"].get("process_CPU_seconds"),
                          "measured_resource_sidecar": result.get("resource_measurement")})
        rows.append({"condition_id": condition["condition_id"], "case_id": cid, "profile": profile, "method": condition["method"],
                     "atomic_scope": cases[cid]["qualification"], "scored_atoms": cases[cid]["scored_atom_count"], "parent_atoms": len(parent),
                     "source_preparation_status": condition["status"], "information_regime": condition["information"], "seeds": seeds})
    return {"schema_version": 1, "status": "complete", "metrics_eligible": True, "method_id": "haddock3_guided_box_air", "counts": counts,
        "conditions": rows, "conditional_family_statistics": _statistics(rows, reference["case_to_family"]),
        "reference_validation": source_validation, "native_provenance": provenance, "source_sha256_inventory": evidence.sources,
        "metric_code_sha256": evaluation["metric_code_sha256"],
        "limitations": ["Separate bounded guided feasibility protocol: ambiguous patch-contact AIRs are not hard Cartesian Vina boxes or equivalent information encoding.",
            "All108 eligible source/method/seed cells remain in the denominator; unavailable source/assembly workflows are not biological nonbinding or a general software capability verdict.",
            "Five full-ligand references and the11-atom4R9F modeled core are distinct strata; full parent chemistry is never truncated to the scored core.",
            "Each model uses its own protein coordinates for proper global/local receptor fits. No independent ligand fit is performed.",
            "All physical initial/refined models and native scores are retained; no clustering, score rounding, best-RMSD selection or missing-slot replacement.",
            "Generated protein atoms are not reference truth; ring/stereo/contact/clash diagnostics do not certify complete physical validity.",
            "Native installation checks are historical metadata; archived numeric inputs/models/io are locally verified without CNS/HADDOCK software.",
            "Native weighted scores rank only within this method and stage; absolute scores are not compared with Vina or Vina-Carb.",
            "Family intervals are conditional descriptions of a purposive structural cohort, not population or training-unseen reliability estimates."]}


def example_plan(root=None):
    root = Path(root) if root is not None else None
    def pin(path):
        record = {"path": path, "sha256": sha256(root / path) if root is not None and (root / path).is_file() else "FILL_SOURCE_SHA256"}
        if root is not None and (root / path).is_file():
            record["bytes"] = (root / path).stat().st_size
        return record
    return {"schema_version": 1, "method_id": "haddock3_guided_box_air", "primary_stage": "flexref", "case_ids": list(CASE_IDS),
        "seeds": list(SEEDS), "rank_budgets": [1, 5], "pose_thresholds_A": [1, 2, 3], "clustering_used_for_primary": False,
        "no_equivalent_hard_box_information_claim": True, "prepared_plan": pin("results/haddock_paired_cohort/prepared_plan.json"),
        "preparation_receipt": pin("results/haddock_paired_cohort/preparation_receipt.json"),
        "native_execution_freeze": pin("results/haddock_paired_cohort/execution_freeze.json"),
        "execution_progress": "results/haddock_paired_cohort/execution_progress.json",
        "native_output_archives": [],
        "orchestration_recovery": {"plan": pin("results/haddock_paired_cohort/orchestration_recovery/recovery_plan.json"),
                                   "freeze": pin("results/haddock_paired_cohort/orchestration_recovery/recovery_freeze.json"),
                                   "adoption_receipt": pin("results/haddock_paired_cohort/orchestration_recovery/adoption_receipt.json"),
                                   "preparation_receipt": pin("results/haddock_paired_cohort/orchestration_recovery/recovery_preparation_receipt.json")},
        "reference_evaluation_plan": pin("benchmarks/structural_pose_reference_evaluation.json"),
        "reference_evaluation_freeze": pin("results/structural_pose_v3/reference_evaluation_freeze.json"),
        "public_source_distribution": pin("benchmarks/source_public_distribution.json"),
        "source_profiles": {"strict": "strict_CCD_ideal", "repair": "modeled_sidechain_CCD_ideal"},
        "metric_code_sha256": {path: pin(path)["sha256"] for path in CODE_SOURCES},
        "freeze_receipt_schema": {"evaluation_plan_sha256": "SHA256_OF_COMPLETED_PLAN", "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False}}


def render_markdown(report):
    if not report["metrics_eligible"]:
        return "# Guided HADDOCK report pending\n\nThe fixed108-cell workflow is incomplete. No native pose accuracy is reported.\n"
    from statistics import mean
    lines = ["# Guided HADDOCK structural test", "", "This is a separate ambiguous-patch protocol. Information encoding and native scores are not equivalent to a hard-box Vina search.", "",
             "| Profile | Scope | Site condition | Cases | Top1 ≤2Å | Top5 ≤2Å | Best retained ≤2Å |", "|---|---|---|---:|---:|---:|---:|"]
    for profile in PROFILES:
        for scope in ("full_ligand", "modeled_core"):
            for method in METHODS:
                rows = [r for r in report["conditions"] if r["profile"] == profile and r["atomic_scope"] == scope and r["method"] == method]
                values = [mean(mean(s["global"][f"{rank}_success_2A"] for s in r["seeds"]) for r in rows) for rank in ("top1", "top5", "best_returned")]
                lines.append(f"| {profile} | {scope} | {method} | {len(rows)} | " + " | ".join(f"{v:.3f}" for v in values) + " |")
    lines.extend(["", "The JSON retains every eligible failure, full-precision individual model identity/score, global/local fits, 1/2/3Å endpoints, lineage and source hashes.", ""] + report["limitations"])
    return "\n".join(lines) + "\n"


def write_report(root, plan_path, freeze_path, summary_output, markdown_output, *, check=False, allow_pending=False):
    report = generate_report(root, plan_path, freeze_path, allow_pending=allow_pending)
    for path, text in ((Path(root) / summary_output, json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"),
                       (Path(root) / markdown_output, render_markdown(report))):
        if check:
            if path.read_text() != text:
                raise ValueError(f"Guided HADDOCK report replay changed: {path}")
        else:
            write_text_atomic(path, text)
    return report
