#!/usr/bin/env python3
"""Read-only verification of saved chemical-ideal ring initialization evidence.

No embedding, force field, optimization, docking or native-coordinate fit runs.
This verifies the archived development procedure and products, not a rerun of
the 20-slot searches or a claim about equilibrium/bound ring populations.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark import ideal_ring_conformers as ir
from aav_benchmark import ring_geometry as rg

BASE = Path("results/ideal_ring_conformer_feasibility")
V3 = BASE / "ccd_ideal_pyranose_v3"
REVIEW = Path("results/common_pose_rescoring/geometry_diagnosis")
RECEIPT_SHA256 = "03997f897fbedeba967e5f9fbd17418e57b9a467350ea8b44fb36a272f720a7e"
AMENDMENT_SHA256 = "857d1dd83b1cf71a27fafc94134c101b2f4eb88919ae0417132121e4fe5e1757"
DERIVED_GEOMETRY_TOLERANCE = 1e-7


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def pin(root, row):
    path = Path(row["path"])
    require(not path.is_absolute() and ".." not in path.parts, "Unsafe receipt artifact path")
    target = root / path
    require(target.is_file() and sha(target) == row["sha256"] and target.stat().st_size == row["bytes"],
            "Ring initialization artifact checksum mismatch: " + str(path))
    return target


def compare_geometry(actual, expected, label="derived geometry"):
    """Roundoff only for recomputed geometry; gates are always evaluated exactly."""
    if isinstance(actual, dict) and isinstance(expected, dict):
        require(set(actual) == set(expected), label + " fields differ")
        for key in actual:
            compare_geometry(actual[key], expected[key], label + "/" + key)
    elif isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        require(len(actual) == len(expected), label + " length differs")
        for i, (left, right) in enumerate(zip(actual, expected)):
            compare_geometry(left, right, label + "/" + str(i))
    elif isinstance(actual, float) and isinstance(expected, float):
        require(math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual, expected, rel_tol=0, abs_tol=DERIVED_GEOMETRY_TOLERANCE),
                label + " differs")
    else:
        require(type(actual) is type(expected) and actual == expected, label + " differs")


def validate_search(search, policy, *, selected):
    attempts = search["attempts"]
    slots = policy["sampling"]["candidate_slots"]
    seed = policy["sampling"]["seed_base"]
    require(slots == search["candidate_slots"] == len(attempts) == 20
            and [a["slot"] for a in attempts] == list(range(20))
            and [a["seed"] for a in attempts] == list(range(seed, seed + 20)),
            "Fixed candidate-slot or seed denominator differs")
    require(all(a["status"] in {"valid", "invalid", "embedding_failed"} for a in attempts),
            "Unknown candidate status")
    valid = [a for a in attempts if a["status"] == "valid"]
    require(len(valid) == search["valid_candidate_count"], "Valid candidate count differs")
    for attempt in valid:
        require(attempt["embedding_return_code"] >= 0 and attempt["mmff_return_code"] == 0
                and attempt["validation"]["valid"] is True and not attempt["validation"]["failures"]
                and math.isfinite(attempt["unconstrained_mmff_energy_at_selected_geometry"]),
                "A recorded valid candidate failed its declared gates")
    if selected:
        require(valid and search["status"] == "selected", "No eligible selected candidate")
        best = min(valid, key=lambda a: (a["unconstrained_mmff_energy_at_selected_geometry"], a["slot"]))
        require(search["selected_slot"] == best["slot"]
                and search["selected_unconstrained_mmff_energy"] == best["unconstrained_mmff_energy_at_selected_geometry"],
                "Selected candidate is not the minimum archived valid MMFF energy")
    else:
        require(not valid and search["status"] == "no_valid_conformer" and search["selected_slot"] is None,
                "Development failure cannot become a selected conformer")
    return Counter(a["status"] for a in attempts)


def inspect_product(molecule, mapping, rings, stereo, policy, expected_smiles, expected_fingerprint):
    require(molecule is not None and molecule.GetNumConformers() == 1, "Missing saved conformer")
    require(ir.graph_fingerprint(molecule) == expected_fingerprint, "Saved complete named chemical graph differs")
    require(all(molecule.GetAtomWithIdx(row["rdkit_index"]).GetProp("canonical_atom_key")
                == row["canonical_atom_key"] for row in mapping), "Serialized atom map differs")
    from3d = Chem.Mol(molecule)
    Chem.AssignStereochemistryFrom3D(from3d, confId=0, replaceExistingTags=True)
    require(Chem.MolToSmiles(Chem.RemoveHs(from3d), isomericSmiles=True) == expected_smiles,
            "Freshly reassigned saved 3D stereochemistry differs")
    result = ir.validate_conformer(molecule, 0, mapping, rings, stereo, policy)
    require(result["valid"], "Saved conformer fails declared chemical/physical gates: " + str(result["failures"]))
    return result


def verify_review_amendment(root, receipt):
    """Reconcile one overwritten derived review without altering its old pin."""
    path = root / BASE / "review_pin_amendment.json"
    require(sha(path) == AMENDMENT_SHA256, "Review provenance amendment changed")
    amendment = load(path)
    require(amendment["amendment_id"] == "initializer_independent_review_provenance_v2"
            and amendment["original_receipt"]["sha256"] == RECEIPT_SHA256
            and amendment["historical_review_pin"] == receipt["validation"]["independent_review"]
            and amendment["historical_review_pin_reproduced"] is False
            and amendment["primary_artifact_count"] == 32
            and all(amendment[key] is False for key in [
                "primary_initializer_artifacts_changed", "initializer_policy_changed",
                "initializer_algorithm_changed", "new_conformer_search_performed_for_this_amendment",
                "new_native_or_docking_execution_performed"]),
            "Review amendment cannot change primary evidence or claim historical byte recovery")
    pin(root, amendment["original_receipt"])
    current = pin(root, amendment["current_derived_review"])
    versioned = pin(root, amendment["versioned_review_copy"])
    require(amendment["current_derived_review"]["path"] == amendment["historical_review_pin"]["path"]
            and amendment["current_derived_review"]["sha256"] != amendment["historical_review_pin"]["sha256"]
            and current.read_bytes() == versioned.read_bytes(),
            "Current versus versioned derived-review identity differs")
    pin(root, amendment["standalone_review_source"])
    pin(root, amendment["independent_expected_parent_graphs"])
    return load(versioned)


def verify_ring_initialization(root=ir.ROOT):
    root = Path(root).resolve()
    receipt_path = root / BASE / "profile_receipt.json"
    require(sha(receipt_path) == RECEIPT_SHA256, "Frozen initializer receipt changed")
    receipt = load(receipt_path)
    artifacts = receipt["artifact_files"]
    require(len(artifacts) == len({row["path"] for row in artifacts}) == 32,
            "Initializer artifact denominator differs")
    for row in artifacts:
        pin(root, row)
    for key in ["policy", "algorithm", "ring_math", "documentation", "tests", "golden_chemical_graphs"]:
        pin(root, receipt[key])
    require(sha(ir.__file__) == receipt["algorithm"]["sha256"]
            and sha(rg.__file__) == receipt["ring_math"]["sha256"], "Imported frozen initializer/math changed")
    review = verify_review_amendment(root, receipt)
    require(sha(root / REVIEW / "review_initializer_outputs.py") == review["review_source_sha256"]
            and sha(root / REVIEW / "expected_parent_graphs.json") == review["golden_expected_graphs_sha256"]
            == receipt["golden_chemical_graphs"]["sha256"], "Independent review source/graph pins differ")
    policy = load(root / ir.POLICY_PATH)
    golden = load(root / receipt["golden_chemical_graphs"]["path"])["cases"]
    cases = policy["cases"]
    require([row["case_id"] for row in receipt["cases"]] == cases
            and set(golden) == set(cases) and len(cases) == 7,
            "Seven-case initializer denominator differs")
    review_cases = {row["case_id"]: row for row in review["cases"]}
    require(len(review_cases) == len(review["cases"]) == 7 and set(review_cases) == set(cases),
            "Independent review case denominator differs")
    math_controls = rg.self_validation()
    require(math_controls["status"] == "passed", "Ring mathematics self-validation failed")
    total = Counter()
    rows = []
    for case in cases:
        report = load(root / V3 / case / "report.json")
        require(report["case_id"] == case and report["no_docking_or_accuracy_evaluation"] is True
                and report["policy_sha256"] == receipt["policy"]["sha256"]
                and report["algorithm_source_sha256"] == receipt["algorithm"]["sha256"]
                and report["ring_math_source_sha256"] == receipt["ring_math"]["sha256"],
                "Case source/profile scope differs")
        schema = ir.case_schema(case, root)
        require(schema == report["schema"], "Frozen chemical input schema differs")
        parent, mapping, rings, stereo, provenance = ir.build_parent_graph(schema, root)
        require(parent.GetNumConformers() == 0 and provenance["experimental_coordinates_read"] is False
                and provenance["isomeric_smiles"] == golden[case]["canonical_isomeric_smiles"],
                "Independent complete chemical graph differs")
        require(mapping == report["atom_map"], "Source-derived atom map differs")
        compare_geometry(rings, report["ring_templates"], case + " chemical-ideal ring templates")
        compare_geometry(stereo, report["stereo_tests"], case + " named stereocenter tests")
        sdf = root / V3 / case / "ligand.sdf"
        require(sha(sdf) == report["selected_sdf_sha256"] == review_cases[case]["saved_sdf_sha256"]
                and sha(root / V3 / case / "report.json") == review_cases[case]["report_sha256"],
                "Saved product or independent review hash differs")
        molecule = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
        validation = inspect_product(molecule, mapping, rings, stereo, policy,
            provenance["isomeric_smiles"], ir.graph_fingerprint(parent))
        compare_geometry(validation, report["serialized_sdf_validation"], case + " serialized validation")
        total.update(validate_search(report["search"], policy, selected=True))
        # A rotation about an interresidue bond changes linkage orientation while
        # preserving all intraring observables. Work only on temporary copies.
        indices = {row["canonical_atom_key"]: row["rdkit_index"] for row in mapping}
        a, b = [indices[key] for key in provenance["declared_links"][0]]
        n1 = next(n.GetIdx() for n in molecule.GetAtomWithIdx(a).GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != b)
        n2 = next(n.GetIdx() for n in molecule.GetAtomWithIdx(b).GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != a)
        quartet = [n1, a, b, n2]
        angle = rdMolTransforms.GetDihedralDeg(molecule.GetConformer(), *quartet)
        freedom = []
        for change in [5., 10., 20., 30.]:
            moved = Chem.Mol(molecule)
            rdMolTransforms.SetDihedralDeg(moved.GetConformer(), *quartet, angle + change)
            altered = ir.validate_conformer(moved, 0, mapping, rings, stereo, policy)
            max_Q_change = max(abs(rg.cremer_pople6(molecule.GetConformer().GetPositions()[r["indices"]])["Q_A"]
                                   - rg.cremer_pople6(moved.GetConformer().GetPositions()[r["indices"]])["Q_A"]) for r in rings)
            require(max_Q_change < 1e-12, "Linkage rotation changed ring puckering")
            freedom.append(dict(rotation_deg=change, all_gates_valid=altered["valid"],
                                failures=altered["failures"], max_ring_Q_change_A=max_Q_change))
        require(any(row["all_gates_valid"] for row in freedom), "No tested linkage rotation preserves validity")
        compare_geometry(freedom, review_cases[case]["temporary_copy_glycosidic_rotation_controls"], case + " linkage freedom")
        rows.append(dict(case_id=case, heavy_atoms=len(mapping), rings=len(rings), stereocenters=len(stereo),
                         valid_candidates=report["search"]["valid_candidate_count"], selected_slot=report["search"]["selected_slot"]))
    require(dict(total) == receipt["validation"]["candidate_status_counts"]
            and sum(total.values()) == 140, "Complete 140-slot status denominator differs")
    history_failures = 0
    for history in receipt["development_history"][:2]:
        old_policy = load(pin(root, history["policy"]))
        old_report = load(pin(root, history["case_report"]))
        require(old_report["no_docking_or_accuracy_evaluation"] is True,
                "Development history information condition changed")
        counts = validate_search(old_report["search"], old_policy, selected=False)
        require(counts == Counter(embedding_failed=20), "Initial development failure history differs")
        history_failures += 20
    for entry in receipt["same_seed_full20_slot_determinism_replay"]["artifacts"].values():
        original, repeated = pin(root, entry["original"]), pin(root, entry["replay"])
        require(entry["byte_identical"] is True and original.read_bytes() == repeated.read_bytes(),
                "Archived same-seed determinism replay differs")
    return dict(status="verified", profile_id=policy["profile_id"], receipt_sha256=RECEIPT_SHA256,
                artifact_records_verified=32, primary_artifact_records_verified=32,
                review_pin_amendment_sha256=AMENDMENT_SHA256,
                historical_independent_review_pin_reproduced=False,
                versioned_independent_review_verified=True,
                cases=rows, complete_heavy_atoms=sum(row["heavy_atoms"] for row in rows),
                rings=sum(row["rings"] for row in rows), stereocenters=sum(row["stereocenters"] for row in rows),
                fixed_v3_slots=140, slot_status_counts=dict(total), earlier_failed_slots=history_failures,
                archived_same_seed_replay_byte_identical=True, independently_reassigned_3D_stereochemistry=True,
                ring_math_self_validation=math_controls["status"], temporary_linkage_rotation_controls=True,
                derived_geometry_roundoff_tolerance=DERIVED_GEOMETRY_TOLERANCE,
                scope="Saved source-only development products; geometry, labels and candidate metadata verified. No embedding, MMFF evaluation, optimization, docking, native fit, population estimate or prospective preregistration is claimed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ir.ROOT)
    args = parser.parse_args()
    try:
        result = verify_ring_initialization(args.root)
    except (ValueError, OSError, KeyError, RuntimeError, TypeError, AssertionError) as error:
        parser.exit(1, f"Ring initialization verification failed: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
