"""Shared prospective CCD-ideal ligand and paired receptor-profile inputs.

Frozen V2 protein predictions supply centers and ranks.  The independently
validated CCD-ideal conformer is exported once per case and copied identically
to strict, modeled-sidechain and separate monomer source profiles.  This module
never reads a returned pose or computes reference-based accuracy.
"""
from __future__ import annotations

import copy
import gzip
import io
import json
import math
from pathlib import Path
import shutil
import tarfile

import numpy as np

from .aav import prepare_receptor
from .blind_docking import ligand_atoms, translate_ligand
from .heldout_docking import CASE_IDS, LEDGER_SHA256, neighborhood_protein
from .provenance import sha256, write_json, write_text_atomic

POLICY_PATH = "benchmarks/ideal_ring_conformer_policy.json"
METHODS = ("fpocket", "p2rank", "oracle")


def _load(path): return json.loads(Path(path).read_text())


def independent_box_size(coordinates):
    xyz = np.asarray(coordinates, float)
    if xyz.ndim != 2 or xyz.shape[1:] != (3,) or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("Independent ligand coordinates must be a nonempty finite Nx3 array")
    radius = float(np.linalg.norm(xyz-xyz.mean(axis=0), axis=1).max())
    side = math.ceil(max(24, 2*radius+16))
    if side > 42: raise ValueError("Independent full-parent conformer exceeds the frozen 42 A cube limit")
    return side, radius


def export_independent_ligand(sdf_path, report_path, original_atom_map_path, original_chemistry_path, output):
    """Export existing validated 3D unchanged, restoring original canonical names."""
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from .ideal_ring_conformers import build_parent_graph, graph_fingerprint, validate_conformer

    sdf_path, report_path, output = Path(sdf_path), Path(report_path), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any((output / name).exists() for name in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json", "chemistry.json")):
        raise FileExistsError("Preserve previous prospective ligand exports")
    report, old_map, old_chemistry = _load(report_path), _load(original_atom_map_path)["atoms"], _load(original_chemistry_path)
    root = Path(__file__).resolve().parents[1]; policy = _load(root / POLICY_PATH)
    if (report["policy_sha256"] != sha256(root / POLICY_PATH) or report["selected_sdf_sha256"] != sha256(sdf_path)
            or report["algorithm_source_sha256"] != sha256(root / "aav_benchmark/ideal_ring_conformers.py")
            or report["search"]["status"] != "selected" or report["chemical_graph"]["experimental_coordinates_read"] is not False
            or report["no_docking_or_accuracy_evaluation"] is not True):
        raise ValueError("CCD-ideal input is not bound to its prospective source-only policy/report")
    molecules = list(Chem.SDMolSupplier(str(sdf_path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None: raise ValueError("Expected one complete independent SDF conformer")
    molecule = molecules[0]
    parent, fresh_map, rings, stereo, provenance = build_parent_graph(report["schema"], root)
    mapping = report["atom_map"]
    if fresh_map != mapping or graph_fingerprint(parent) != graph_fingerprint(molecule):
        raise ValueError("SDF full-parent graph, charges or stereo differ from pinned CCD chemistry")
    old_by_key = {row["canonical_atom_key"]: row for row in old_map}
    if len(old_by_key) != len(old_map) or set(old_by_key) != {row["canonical_atom_key"] for row in mapping}:
        raise ValueError("New and frozen full-parent canonical maps differ")
    if len({row["name"] for row in old_map}) != len(old_map): raise ValueError("Frozen ligand names are ambiguous")
    for row in mapping:
        old = old_by_key[row["canonical_atom_key"]]
        for key in ("branch_num", "ccd_id", "atom_name", "element", "formal_charge", "ccd_stereo_config"):
            if str(old[key]) != str(row[key]): raise ValueError(f"New ligand changed canonical chemical field {key}")
    validation = validate_conformer(molecule, 0, mapping, rings, stereo, policy)
    fresh_stereo = Chem.Mol(molecule); Chem.AssignAtomChiralTagsFromStructure(fresh_stereo, replaceExistingTags=True)
    if not validation["valid"] or graph_fingerprint(fresh_stereo) != graph_fingerprint(parent):
        raise ValueError("Independent serialized SDF fails geometry or fresh-coordinate stereochemistry")
    new_map = [{**row, "name": old_by_key[row["canonical_atom_key"]]["name"]} for row in mapping]
    heavy_by_index = {row["rdkit_index"]: row for row in new_map}
    hydrogen_number = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            row = heavy_by_index[atom.GetIdx()]
            if atom.GetProp("canonical_atom_key") != row["canonical_atom_key"]: raise ValueError("SDF canonical atom property/order changed")
            name, component, number = row["name"], row["ccd_id"], int(row["branch_num"])
        else:
            hydrogen_number += 1; name, component, number = f"H{hydrogen_number}", "LIG", 1
        info = Chem.AtomPDBResidueInfo(); info.SetName(name.rjust(4)); info.SetResidueName(component)
        info.SetResidueNumber(number); info.SetChainId("L"); atom.SetMonomerInfo(info)
    setups = MoleculePreparation().prepare(molecule)
    if len(setups) != 1: raise ValueError("Expected one Meeko full-parent export")
    text, success, error = PDBQTWriterLegacy.write_string(setups[0])
    if not success: raise ValueError(error)
    atoms = ligand_atoms(text); heavy = {a["name"]: a for a in atoms if not a["hydrogen"]}
    if set(heavy) != {row["name"] for row in new_map}: raise ValueError("PDBQT export lost a full-parent heavy atom")
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and not math.isfinite(float(line[70:76])):
            raise ValueError("Meeko exported a nonfinite partial charge")
    original_xyz = molecule.GetConformer().GetPositions()
    maximum = max(float(np.linalg.norm(original_xyz[row["rdkit_index"]]-heavy[row["name"]]["coordinates"])) for row in new_map)
    if maximum > .0009: raise ValueError("PDBQT serialization changed independent heavy geometry")
    # Verify the rounded actual docking coordinates against the same named gates.
    rounded = Chem.Mol(molecule)
    for row in new_map: rounded.GetConformer().SetAtomPosition(row["rdkit_index"], heavy[row["name"]]["coordinates"])
    rounded_validation = validate_conformer(rounded, 0, mapping, rings, stereo, policy)
    if not rounded_validation["valid"]: raise ValueError("PDBQT rounding violates the frozen CCD-ideal geometry/stereo gates")
    side, radius = independent_box_size([heavy[row["name"]]["coordinates"] for row in new_map])
    shutil.copy2(sdf_path, output / "ligand.sdf"); shutil.copy2(report_path, output / "initializer_report.json")
    write_text_atomic(output / "ligand.pdbqt", text)
    write_json(output / "ligand_atom_map.json", {"atoms": new_map})
    chemistry = {**old_chemistry, "isomeric_smiles": provenance["isomeric_smiles"], "heavy_atoms": len(new_map),
                 "formal_charge": provenance["formal_charge"], "embedding": report["search"]["coordinates_source"],
                 "native_coordinates_used_for_embedding": False, "initializer_profile": report["profile_id"],
                 "initializer_policy": {"path": POLICY_PATH, "sha256": report["policy_sha256"]},
                 "initializer_report_sha256": sha256(report_path), "initializer_sdf_sha256": sha256(sdf_path),
                 "signed_chirality_tests": rounded_validation["stereochemistry"],
                 "box_size": [side]*3, "heavy_atom_centroid_radius_A": radius,
                 "box_rule": "cube side ceil(max(24,2R+16)) A, R from actual serialized independent heavy PDBQT; maximum42 A",
                 "serialization": {"sdf_sha256": sha256(output / "ligand.sdf"), "pdbqt_sha256": sha256(output / "ligand.pdbqt"),
                                   "max_heavy_serialization_displacement_A": maximum, "rounded_CCD_ideal_validation": rounded_validation},
                 "full_parent_canonical_key_and_name_bijection_preserved": True,
                 "no_new_embedding_or_minimization_during_export": True}
    write_json(output / "chemistry.json", chemistry)
    return chemistry


def _portable_archive(source, output, stem):
    source, output = Path(source), Path(output); output.mkdir(parents=True, exist_ok=True)
    archive_path, manifest_path = output / f"{stem}.tar.gz", output / f"{stem}_manifest.json"
    if archive_path.exists() or manifest_path.exists(): raise FileExistsError("Preserve previous prospective input archives")
    files = []
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, filename="", mtime=0, mode="wb") as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in sorted(source.rglob("*")):
                    if not path.is_file(): continue
                    content = path.read_bytes(); name = str(path.relative_to(source))
                    info = tarfile.TarInfo(name); info.size = len(content); info.mode = 0o644; info.mtime = 0
                    archive.addfile(info, io.BytesIO(content)); files.append({"path": name, "sha256": sha256(path), "bytes": len(content)})
    record = {"schema_version": 1, "archive": archive_path.name, "archive_sha256": sha256(archive_path), "files": files,
              "purpose": "prospective CCD-ideal ligand and fixed receptor profile inputs; no pose results"}
    write_json(manifest_path, record)
    return {"path": str(manifest_path), "sha256": sha256(manifest_path)}


def integrate_sources(root, strict, repair, monomer, conformers, output, archive_output):
    """Create fresh paired source profiles, changing only ligand initializer files."""
    root, strict, repair, monomer, conformers, output, archive_output = map(Path, (root, strict, repair, monomer, conformers, output, archive_output))
    if output.exists(): raise FileExistsError("Preserve all earlier source/conformer profiles")
    output.mkdir(parents=True)
    export = output / "shared_ligands"; export.mkdir()
    ligand_receipts = {}
    for case_id in CASE_IDS:
        destination = export / case_id
        ligand_receipts[case_id] = export_independent_ligand(conformers / case_id / "ligand.sdf", conformers / case_id / "report.json",
                                                           strict / case_id / "ligand_atom_map.json", strict / case_id / "chemistry.json", destination)
    profiles = {}
    for profile, original, cases in (("strict", strict, CASE_IDS), ("repair", repair, CASE_IDS), ("monomer", monomer, ("3ZWE", "5HZA"))):
        source = output / profile; source.mkdir(); rows = []
        for case_id in cases:
            prior = original / case_id; directory = source / case_id; shutil.copytree(prior, directory)
            old_files = {str(p.relative_to(prior)): sha256(p) for p in prior.rglob("*") if p.is_file()}
            for name in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json", "chemistry.json", "initializer_report.json"):
                shutil.copy2(export / case_id / name, directory / name)
            protected = [name for name in old_files if name not in ("ligand.sdf", "ligand.pdbqt", "ligand_atom_map.json", "chemistry.json")]
            if any(sha256(directory / name) != old_files[name] for name in protected):
                raise ValueError("CCD-ideal integration changed a nonligand source artifact")
            row = {"case_id": case_id, "profile": profile, "receptor_sha256": sha256(directory / "receptor.cif"),
                   "receptor_PDB_sha256": sha256(directory / "receptor.pdb"),
                   "original_prediction_input_sha256": sha256(strict / case_id / "receptor.cif"),
                   "original_prediction_configuration": "benchmarks/structural_site_v2_experiment.json",
                   "canonical_atom_map_sha256": sha256(directory / "ligand_atom_map.json"),
                   "ligand_sha256": sha256(directory / "ligand.pdbqt"), "ligand_sdf_sha256": sha256(directory / "ligand.sdf"),
                   "chemistry_sha256": sha256(directory / "chemistry.json"), "initializer_report_sha256": sha256(directory / "initializer_report.json"),
                   "box_size": ligand_receipts[case_id]["box_size"], "all_prior_nonligand_files_byte_identical": True,
                   "source_receptor_geometry_unchanged": True, "native_ligand_coordinate_access": False,
                   "atomic_reference_scope_changed": False, "prior_source_path": str(prior)}
            if profile == "repair":
                row["receptor_support"] = _load(prior / "source_preparation.json")["status"]
            elif profile == "strict":
                support = _load(prior / "source_preparation.json").get("receptor_parameterization", {}).get("status")
                row["receptor_support"] = "prepared" if support == "prepared" else "unsupported" if support == "unsupported" else "preparation_failed"
            else: row["receptor_support"] = "source_only_control_not_parameterized"
            rows.append(row)
        manifest = {"schema_version": 1, "profile": profile, "case_ids": list(cases), "source_ledger_sha256": LEDGER_SHA256,
                    "native_ligand_coordinate_access": False, "new_predictions_or_docking_results_read": False,
                    "protein_geometry_changed_by_initializer_integration": False,
                    "ligand_initializer_policy_sha256": sha256(root / POLICY_PATH), "cases": rows,
                    "evaluation_reference": "unchanged full biological-assembly qualified source copies and original observed protein maps; generated protein atoms never truth",
                    "center_policy": "same original full-assembly V2 predictor centers/ranks for strict/repair; monomer site controls have a separate prospective freeze",
                    "code_sha256": sha256(Path(__file__))}
        write_json(source / "profile_manifest.json", manifest)
        archive = _portable_archive(source, archive_output, f"{profile}_sources")
        profiles[profile] = {"source_directory": str(source), "manifest": {"path": str(source / "profile_manifest.json"), "sha256": sha256(source / "profile_manifest.json")},
                             "source_archive_manifest": archive, "case_ids": list(cases)}
    summary = {"schema_version": 1, "profiles": profiles, "initializer_policy_sha256": sha256(root / POLICY_PATH),
               "shared_ligands": {case: {"box_size": r["box_size"], "heavy_atoms": r["heavy_atoms"], "native_coordinates_used_for_embedding": False,
                                               "ligand_pdbqt_sha256": sha256(export / case / "ligand.pdbqt")} for case, r in ligand_receipts.items()},
               "no_scientific_prediction_docking_or_accuracy_evaluation": True, "input_adapter_sha256": sha256(Path(__file__))}
    write_json(archive_output / "source_integration.json", summary)
    return summary


def validate_candidate_reuse(selection, case_id, method, source_prediction_config_sha256, expected_input_sha256):
    if case_id not in CASE_IDS or method not in ("fpocket", "p2rank"):
        raise ValueError("Unknown frozen case/predictor")
    if (selection["case_id"] != case_id or selection["method"] != method
            or selection["configuration_sha256"] != source_prediction_config_sha256
            or selection["prediction_input_sha256"] != expected_input_sha256
            or selection.get("native_coordinate_access") is not False):
        raise ValueError("Candidate selection changed original prediction identity/input")
    candidates = selection["candidates"]
    if len(candidates) > 5 or [row["rank"] for row in candidates] != list(range(1, len(candidates)+1)):
        raise ValueError("Original candidate classes are not contiguous within the fixed five budget")
    original_ranks = [row["original_rank"] for row in candidates]
    if (original_ranks != sorted(set(original_ranks)) or any(isinstance(rank, bool) or not isinstance(rank, int) or rank < 1 for rank in original_ranks)
            or any(np.asarray(row["center_A"], float).shape != (3,) or not np.isfinite(np.asarray(row["center_A"], float)).all() for row in candidates)):
        raise ValueError("Original candidate rank/center identity is invalid")
    if selection["predictor_status"] not in ("ok", "completed") and candidates:
        raise ValueError("Failed predictor cannot supply candidate boxes")
    return candidates


def select_execution_candidates(candidates, budget):
    """Execute the single predeclared site while retaining all original ranks."""
    if isinstance(budget, bool) or not isinstance(budget, int) or budget != 1:
        raise ValueError("The matched prospective Vina/VC comparison executes exactly the top-one site budget")
    return candidates[:budget]


def _verify_pose_freeze(root, source, profile, config_path, freeze_path):
    root, source = Path(root), Path(source); config, freeze = _load(config_path), _load(freeze_path)
    if (freeze.get("configuration_sha256") != sha256(config_path) or freeze.get("frozen_before_new_docking") is not True
            or freeze.get("input_adapter_sha256") != sha256(Path(__file__))):
        raise ValueError("Prospective pose preparation requires coordinator freeze of exact config and input adapter")
    select_execution_candidates([], config.get("primary_candidate_budget"))
    descriptor = config["source_profiles"][profile]
    if descriptor["manifest"]["sha256"] != sha256(source / "profile_manifest.json"):
        raise ValueError("Frozen source profile manifest changed")
    manifest = _load(source / "profile_manifest.json")
    if manifest["profile"] != profile or manifest["case_ids"] != list(CASE_IDS):
        raise ValueError("Paired strict/repair pose profiles retain exactly six cases")
    for row in manifest["cases"]:
        directory = source / row["case_id"]
        for filename, field in (("receptor.cif", "receptor_sha256"), ("receptor.pdb", "receptor_PDB_sha256"),
                                ("ligand.pdbqt", "ligand_sha256"), ("ligand.sdf", "ligand_sdf_sha256"),
                                ("ligand_atom_map.json", "canonical_atom_map_sha256"), ("chemistry.json", "chemistry_sha256"),
                                ("initializer_report.json", "initializer_report_sha256")):
            if sha256(directory / filename) != row[field]: raise ValueError(f"Prospective source input changed: {row['case_id']}/{filename}")
    return config, manifest


def prepare_paired_condition(root, source, profile, case_id, method, prediction_directory, output, config_path, freeze_path, *, oracle_center=None):
    """Prepare identical centers/ranks and ligand geometry for paired receptors."""
    config, manifest = _verify_pose_freeze(root, source, profile, config_path, freeze_path)
    if profile not in ("strict", "repair") or case_id not in CASE_IDS or method not in METHODS:
        raise ValueError("Use a fixed full-assembly paired receptor profile and method")
    source, output = Path(source), Path(output); directory = source / case_id
    if output.exists(): raise FileExistsError("Preserve earlier candidate preparations")
    output.mkdir(parents=True)
    source_record = next(row for row in manifest["cases"] if row["case_id"] == case_id)
    condition = {"schema_version": 1, "case_id": case_id, "condition_id": case_id, "complex_id": case_id,
                 "profile": profile, "method": method, "configuration_sha256": sha256(config_path),
                 "fixed_denominator_member": True, "status": "preparing", "jobs": [],
                 "source_profile_manifest_sha256": sha256(source / "profile_manifest.json"),
                 "native_coordinate_access": method == "oracle", "source_ledger_sha256": LEDGER_SHA256,
                 "executed_primary_candidate_budget": config["primary_candidate_budget"],
                 "all_prediction_candidates_retained_in_selection": True}
    if method == "oracle":
        center = np.asarray(oracle_center, float)
        if center.shape != (3,) or not np.isfinite(center).all() or prediction_directory is not None:
            raise ValueError("Oracle requires only the explicit common native-informed full-parent center")
        candidates = [{"rank": 1, "original_rank": 1, "center_A": center.tolist()}]
    else:
        if oracle_center is not None: raise ValueError("Native center forbidden in predicted-site arm")
        prediction_directory = Path(prediction_directory); selection_path = prediction_directory / "candidate_selection.json"
        selection = _load(selection_path); prediction_config = config["prediction_configuration"]
        if sha256(Path(root) / prediction_config["path"]) != prediction_config["sha256"]: raise ValueError("Original prediction config changed")
        candidates = validate_candidate_reuse(selection, case_id, method, prediction_config["sha256"], source_record["original_prediction_input_sha256"])
        record = selection["prediction_record"]
        if Path(record["filename"]).name != record["filename"] or sha256(prediction_directory / record["filename"]) != record["sha256"]:
            raise ValueError("Original predictor record changed")
        condition.update(candidate_selection=selection, primary_candidate_selection={"path": str(selection_path), "sha256": sha256(selection_path)},
                         candidate_centers_and_ranks_exactly_reused=True)
    condition["available_prediction_candidates"] = len(candidates)
    candidates = select_execution_candidates(candidates, config["primary_candidate_budget"])
    condition["scheduled_candidate_jobs"] = len(candidates)
    chemistry, atom_map = _load(directory / "chemistry.json"), _load(directory / "ligand_atom_map.json")["atoms"]
    for candidate in candidates:
        box_id = f"box{candidate['rank']:02d}"; box = output / box_id; box.mkdir()
        job = {"id": box_id, "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"],
               "center_A": candidate["center_A"], "prepared_directory": box_id}
        try:
            if source_record["receptor_support"] == "unsupported": raise ValueError("Source receptor parameterization unsupported; MSE preserved without mutation or unvalidated charge/type model")
            protein, neighborhood = neighborhood_protein(directory / "receptor.cif", candidate["center_A"], chemistry["box_size"])
            preparation = prepare_receptor(protein, box / "neighborhood")
            ligand, translation = translate_ligand((directory / "ligand.pdbqt").read_text(), candidate["center_A"])
            write_text_atomic(box / "ligand.pdbqt", ligand); shutil.copy2(directory / "ligand_atom_map.json", box / "ligand_atom_map.json")
            write_json(box / "reference.json", {"names": [row["name"] for row in atom_map], "canonical_atom_map": atom_map,
                                                 "pose_reference_eligible": False, "purpose": "complete canonical heavy identities only; atomic truth evaluated separately"})
            receipt = {"case_id": case_id, "family_id": case_id, "profile": profile, "status": "prepared",
                       "method": method, "candidate_class_rank": candidate["rank"], "candidate_original_rank": candidate["original_rank"],
                       "cohort_role": "fixed_workflow_untouched_structural_test", "configuration_sha256": sha256(config_path),
                       "information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                       "execution_information_regime": "true-site-local" if method == "oracle" else "predicted-site-local",
                       "box_center": candidate["center_A"], "box_size": chemistry["box_size"], "ligand_file": "ligand.pdbqt",
                       "receptors": [{"id": "neighborhood", "state": "holo", "assembly_condition": "complete_8A_box_neighborhood", "pdb_file": "neighborhood.pdb", "pdbqt_file": "neighborhood.pdbqt"}],
                       "neighborhood": neighborhood, "ligand_translation": translation, "receptor_preparation": preparation,
                       "source_profile_manifest_sha256": condition["source_profile_manifest_sha256"],
                       "source_receptor_sha256": source_record["receptor_sha256"], "independent_ligand_sha256": source_record["ligand_sha256"],
                       "initializer_report_sha256": source_record["initializer_report_sha256"], "native_coordinates_in_runner_reference": False,
                       "model_family_heldout_claim": False, "modeled_protein_atoms_are_atomic_truth": False}
            write_json(box / "preparation.json", receipt)
            job.update(status="prepared", preparation_sha256=sha256(box / "preparation.json"))
        except (ValueError, RuntimeError, OSError) as error:
            job.update(status="unsupported" if source_record["receptor_support"] == "unsupported" else "preparation_failed", reason=str(error), error_type=type(error).__name__)
            write_json(box / "preparation_attempt.json", job)
        condition["jobs"].append(job); write_json(output / "condition.json", condition)
    condition["status"] = ("prepared" if candidates and all(job["status"] == "prepared" for job in condition["jobs"])
                           else selection["status"] if method != "oracle" and not candidates else "preparation_failed")
    write_json(output / "condition.json", condition)
    return condition
