"""Saved ring-product replay must not launch conformer search or optimization."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import re
import shutil

import pytest
from rdkit import Chem

from aav_benchmark import ideal_ring_conformers as ir

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ring_product_verifier", ROOT / "scripts/verify_ring_initialization.py")
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


def test_actual_saved_product_replay_never_embeds_optimizes_or_evaluates_forcefields(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Read-only verifier attempted native conformer/force-field execution")
    for name in ["EmbedMolecule", "MMFFGetMoleculeForceField", "MMFFGetMoleculeProperties"]:
        monkeypatch.setattr(ir.AllChem, name, forbidden)
    before = verifier.sha(ROOT / verifier.BASE / "profile_receipt.json")
    result = verifier.verify_ring_initialization(ROOT)
    assert result["artifact_records_verified"] == 32
    assert result["primary_artifact_records_verified"] == 32
    assert result["historical_independent_review_pin_reproduced"] is False
    assert result["versioned_independent_review_verified"] is True
    assert result["complete_heavy_atoms"] == 226
    assert result["rings"] == 21
    assert result["stereocenters"] == 100
    assert result["fixed_v3_slots"] == 140
    assert result["slot_status_counts"] == {"valid": 86, "invalid": 19, "embedding_failed": 35}
    assert result["earlier_failed_slots"] == 40
    assert result["archived_same_seed_replay_byte_identical"] is True
    assert result["temporary_linkage_rotation_controls"] is True
    assert verifier.sha(ROOT / verifier.BASE / "profile_receipt.json") == before


def test_a_recorded_selected_slot_cannot_bypass_minimum_valid_energy_or_complete_seed_budget():
    report = verifier.load(ROOT / verifier.V3 / "2RDK/report.json")
    policy = verifier.load(ROOT / ir.POLICY_PATH)
    search = deepcopy(report["search"])
    search["selected_slot"] = (search["selected_slot"] + 1) % 20
    with pytest.raises(ValueError, match="minimum archived valid"):
        verifier.validate_search(search, policy, selected=True)
    search = deepcopy(report["search"])
    search["attempts"].pop()
    with pytest.raises(ValueError, match="seed denominator"):
        verifier.validate_search(search, policy, selected=True)


def test_saved_graph_tags_cannot_hide_reflected_three_dimensional_stereochemistry():
    parent, mapping, rings, stereo, provenance = ir.build_parent_graph(ir.case_schema("2RDK"))
    molecule = Chem.SDMolSupplier(str(ROOT / verifier.V3 / "2RDK/ligand.sdf"), removeHs=False)[0]
    # Preserve every stored graph/chiral tag while reflecting the coordinates.
    conformer = molecule.GetConformer()
    for index, point in enumerate(conformer.GetPositions()):
        conformer.SetAtomPosition(index, [-point[0], point[1], point[2]])
    policy = verifier.load(ROOT / ir.POLICY_PATH)
    with pytest.raises(ValueError, match="reassigned saved 3D stereochemistry"):
        verifier.inspect_product(molecule, mapping, rings, stereo, policy,
                                 provenance["isomeric_smiles"], ir.graph_fingerprint(parent))


def test_missing_preserved_development_artifact_cannot_be_skipped(tmp_path):
    receipt = verifier.load(ROOT / verifier.BASE / "profile_receipt.json")
    record = receipt["artifact_files"][0]
    with pytest.raises(ValueError, match="artifact checksum mismatch"):
        verifier.pin(tmp_path, record)


@pytest.mark.parametrize("relative", [
    "ccd_ideal_pyranose_v3/2RDK/ligand.sdf",
    "ccd_ideal_pyranose_v3/5V6F/report.json",
    "profile_v1/1GNY/report.json",
])
def test_review_amendment_cannot_waive_a_modified_primary_artifact(tmp_path, relative):
    receipt = verifier.load(ROOT / verifier.BASE / "profile_receipt.json")
    for row in receipt["artifact_files"]:
        target = tmp_path / row["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(ROOT / row["path"])
    for name in ["profile_receipt.json", "review_pin_amendment.json", "independent_review_v2.json"]:
        shutil.copyfile(ROOT / verifier.BASE / name, tmp_path / verifier.BASE / name)
    target = tmp_path / verifier.BASE / relative
    target.unlink()
    target.write_bytes(b"synthetic changed primary evidence\n")
    with pytest.raises(ValueError, match="checksum mismatch: " + re.escape(str(verifier.BASE / relative))):
        verifier.verify_ring_initialization(tmp_path)


def test_geometry_roundoff_never_changes_gating_decisions_or_integer_denominators():
    verifier.compare_geometry({"Q": .6 + 1e-10}, {"Q": .6})
    with pytest.raises(ValueError, match="differs"):
        verifier.compare_geometry({"valid": False}, {"valid": True})
    with pytest.raises(ValueError, match="differs"):
        verifier.compare_geometry({"ring_count": 6}, {"ring_count": 5})


def test_receipt_rebinding_is_rejected_before_any_case_artifact_access(tmp_path):
    path = tmp_path / verifier.BASE / "profile_receipt.json"
    path.parent.mkdir(parents=True)
    receipt = verifier.load(ROOT / verifier.BASE / "profile_receipt.json")
    receipt["artifact_files"].pop()
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="Frozen initializer receipt changed"):
        verifier.verify_ring_initialization(tmp_path)


def test_review_amendment_cannot_claim_recovery_of_unavailable_historical_bytes(tmp_path):
    path = tmp_path / verifier.BASE / "review_pin_amendment.json"
    path.parent.mkdir(parents=True)
    receipt = verifier.load(ROOT / verifier.BASE / "profile_receipt.json")
    amendment = verifier.load(ROOT / verifier.BASE / "review_pin_amendment.json")
    amendment["historical_review_pin_reproduced"] = True
    path.write_text(json.dumps(amendment))
    with pytest.raises(ValueError, match="Review provenance amendment changed"):
        verifier.verify_review_amendment(tmp_path, receipt)
