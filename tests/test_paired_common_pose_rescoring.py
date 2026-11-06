"""The exact union retains source identities, hydrogens and failed pools."""
from pathlib import Path

import pytest

from aav_benchmark.additional_method import _atom, export_vina_carb_ligand, read_carb_poses
from aav_benchmark.paired_common_pose_rescoring import (
    restore_carb_pose_aliases, rank_union_pools, score_available_union, export_source_pose,
)
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.common_pose_rescoring import model_blocks
import hashlib

ROOT = Path(__file__).resolve().parents[1]


def tutorial_aliases(tmp_path):
    source = ROOT / "results/additional_method_reproduction/example_inputs/LeY-xray.pdbqt"
    native = ROOT / "results/additional_method_reproduction/example_exact/LeY-xray_vc_out.pdbqt"
    lines = [];mapping = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = _atom(line);name = f"L{atom['serial'] - 1:03d}"
            if atom["type"] not in {"H", "HD", "HS"}:
                mapping.append({"name": name, "source_key": ["A", atom["resnum"], "", atom["resname"], atom["name"]]})
            line = line[:12] + name + line[16:17] + "LIG" + line[20:22] + f"{1:4d}" + line[26:]
        lines.append(line)
    generic = tmp_path / "prepared.pdbqt";generic.write_text("\n".join(lines) + "\n")
    export = export_vina_carb_ligand(generic, mapping, tmp_path / "canonical_input.pdbqt")
    atoms = {atom["serial"]: atom for atom in export["atoms"]}
    rows = []
    for line in native.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = atoms[_atom(line)["serial"]]
            line = line[:12] + f" {atom['carb_name']:<3}" + line[16:17] + f"{atom['carb_resname']:>3}" + line[20:22] + f"{atom['carb_resnum']:4d}" + line[26:]
        rows.append(line)
    raw = tmp_path / "raw_carb_poses.pdbqt";raw.write_text("\n".join(rows) + "\n")
    return generic, raw, export


def test_canonical_alias_roundtrip_preserves_every_model_and_hydrogen_coordinate(tmp_path):
    prepared, raw, export = tutorial_aliases(tmp_path)
    assert len(read_carb_poses(raw, export)) == 10
    restored = tmp_path / "generic_poses.pdbqt"
    receipt = restore_carb_pose_aliases(prepared, raw, export, restored)
    before = [_atom(line) for line in raw.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    after = [_atom(line) for line in restored.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    assert len(before) == len(after) == receipt["copied_atom_records"]
    assert any(atom["type"] == "HD" for atom in before)
    for a, b in zip(before, after):
        assert all(a[key] == b[key] for key in ("serial", "coordinates", "type", "charge"))
    assert receipt["all_heavy_and_hydrogen_coordinates_preserved"] is True


def test_canonical_alias_decoder_rejects_changed_atom_typing_or_map(tmp_path):
    prepared, raw, export = tutorial_aliases(tmp_path)
    export["atoms"][0]["atom_type"] = "OA"
    with pytest.raises(ValueError):
        restore_carb_pose_aliases(prepared, raw, export, tmp_path / "bad.pdbqt")


def test_carb_source_generic_canonical_roundtrip_keeps_all_fixed_model_coordinates(tmp_path):
    prepared, raw, export = tutorial_aliases(tmp_path)
    map_path = tmp_path / "original_export.json";write_json(map_path, export)
    chemical_map = tmp_path / "chemical_map.json"
    write_json(chemical_map, {"atoms": [{"name": atom["generated_name"], "source_key": atom["source_key"]}
                                        for atom in export["atoms"] if not atom["hydrogen"]]})
    receptor = tmp_path / "receptor.pdbqt";receptor.write_text("same paired receptor")
    row = {"ligand": {"path": prepared.name}, "canonical_atom_map": {"path": chemical_map.name},
           "receptor": {"path": receptor.name}}
    source = {"source_algorithm": "vina_carb_1_0", "source_pose_file": raw.name, "model_rank": 1,
              "source_ligand_sha256": sha256(prepared), "source_pose_file_sha256": sha256(raw),
              "source_model_sha256": hashlib.sha256(model_blocks(raw.read_bytes())[1]).hexdigest(),
              "source_export_map": {"path": map_path.name, "sha256": sha256(map_path)},
              "source_canonical_atom_map_sha256": sha256(chemical_map), "source_receptor_sha256": sha256(receptor)}
    target = tmp_path / "fixed_pose.pdbqt"
    receipt = export_source_pose(source, row, tmp_path, target)
    before = [_atom(line.decode()) for line in model_blocks(raw.read_bytes())[1].splitlines() if line.startswith((b"ATOM  ", b"HETATM"))]
    after = [_atom(line) for line in target.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    assert len(before) == len(after)
    assert all(a["coordinates"] == b["coordinates"] and a["type"] == b["type"] and a["charge"] == b["charge"]
               for a, b in zip(before, after))
    assert len(receipt["expected_glycosidic_links"]) == 3
    assert receipt["canonical_alias_roundtrip"]["all_heavy_and_hydrogen_coordinates_preserved"] is True


def pool_fixture():
    a = {"pose_id": "pool/modern/p01", "pool_id": "pool", "complex_id": "CASE", "source_algorithm_rank": 0,
         "candidate_class_rank": 1, "model_rank": 1, "source_model_sha256": "identical_raw_coordinates"}
    b = {**a, "pose_id": "pool/carb/p01", "source_algorithm_rank": 1}
    pool = {"pool_id": "pool", "case_id": "CASE", "source_statuses": {"modern": "ok", "carb": "ok"},
            "eligible_workflow_pool": True, "original_source_orders": {"modern": [a["pose_id"]], "carb": [b["pose_id"]]}}
    index = {"poses": [a, b], "pools": [pool]}
    records = {a["pose_id"]: {"status": "ok", "source": a, "score": {"base_shared0": -8., "combined_shared0": -3., "torsion_normalization_scale": .5}},
               b["pose_id"]: {"status": "ok", "source": b, "score": {"base_shared0": -8., "combined_shared0": -8., "torsion_normalization_scale": .5}}}
    return index, records


def test_identical_source_poses_remain_distinct_and_both_coefficients_use_exact_same_union():
    index, records = pool_fixture();rows, conflicts = rank_union_pools(index, records)
    pool = rows[0]
    assert not conflicts and pool["status"] == "completed"
    assert len(pool["union_pose_ids"]) == 2
    assert pool["coefficient0_order"] == ["pool/modern/p01", "pool/carb/p01"]
    assert pool["coefficient1_order"] == ["pool/carb/p01", "pool/modern/p01"]
    assert set(pool["coefficient0_order"]) == set(pool["coefficient1_order"]) == set(pool["union_pose_ids"])
    assert pool["original_source_orders"] == index["pools"][0]["original_source_orders"]


def test_any_failed_source_or_pose_invalidates_complete_union_instead_of_dropping_it():
    index, records = pool_fixture();records["pool/carb/p01"]["status"] = "pose_scoring_failed"
    pool = rank_union_pools(index, records)[0][0]
    assert pool["status"] == "pose_scoring_failure_whole_union_required" and pool["coefficient1_order"] == []
    assert len(pool["union_pose_ids"]) == 2 and pool["eligible_workflow_pool"] is True
    index["pools"][0]["source_statuses"]["carb"] = "failed"
    pool = rank_union_pools(index, records)[0][0]
    assert pool["status"] == "source_workflow_failure_whole_pair_required" and pool["coefficient0_order"] == []


def test_normalization_conflict_blocks_the_case_pool():
    index, records = pool_fixture()
    records["pool/carb/p01"]["score"]["torsion_normalization_scale"] = .6
    rows, conflicts = rank_union_pools(index, records)
    assert conflicts == ["CASE"] and rows[0]["status"] == "normalization_conflict"
    assert rows[0]["coefficient0_order"] == rows[0]["coefficient1_order"] == []


def test_unfrozen_union_plan_cannot_call_the_native_helper(tmp_path, monkeypatch):
    from aav_benchmark import paired_common_pose_rescoring as module
    plan = tmp_path / "proposal.json";write_json(plan, {"frozen_before_fixed_pose_scoring": False})
    monkeypatch.setattr(module, "_invoke_fixed_pose", lambda *a, **k: pytest.fail("Native scoring reached before freeze"))
    with pytest.raises(ValueError, match="Separate root freeze"):
        score_available_union(root=tmp_path, plan_path=plan, receipt_path=None, installation=tmp_path / "unused")
