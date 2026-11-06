"""Synthetic controls for the OLD evaluator's named full-parent ordering seam."""
import ast
import copy
import gzip
import hashlib
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from aav_benchmark import qualified_pose_v2 as v2
from aav_benchmark.engine import read_poses
from aav_benchmark.provenance import sha256, write_json


def _pose(coordinates, *, rank=1, score=-5.):
    return {"rank": rank, "vina_score": score, "coordinates": coordinates}


def _box(names, coordinates, *, box_id="box01", candidate_rank=1):
    return {"id": box_id, "names": names, "candidate_class_rank": candidate_rank,
            "candidate_original_rank": candidate_rank + 2,
            "poses_by_seed": {7: [_pose(coordinates)]},
            "atom_order_reference": {"archive": "raw/ref.gz", "archive_sha256": "a" * 64,
                                     "reference_sha256": "b" * 64}}


def _synthetic_parent():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles("C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=13) == 0
    heavy = Chem.RemoveHs(molecule)
    xyz = np.round(np.asarray(heavy.GetConformer().GetPositions()), 3)
    names = [f"L{i:03d}" for i in range(len(xyz))]
    reference = {"names": names, "coordinates": xyz.tolist(), "atom_map": [
        {"name": name, "rdkit_index": i, "source_key": ["B", 1, "", "LIG", str(i)]}
        for i, name in enumerate(names)]}
    return names, xyz, reference, Chem.MolToMolBlock(molecule).encode() + b"\n$$$$\n"


def _parsed_box(tmp_path, names, xyz, order, *, candidate_rank=1, shift=None):
    raw_names = [names[i] for i in order]
    points = xyz if shift is None else xyz + shift
    lines = ["MODEL 1", "REMARK VINA RESULT: -5.0 0 0"]
    for serial, i in enumerate(order, 1):
        p = points[i]
        lines.append(f"ATOM  {serial:5d} {names[i]:>4s} LIG A   1    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00     0.000  C")
    lines.append("ENDMDL")
    path = tmp_path / f"box{candidate_rank:02d}.pdbqt"
    path.write_text("\n".join(lines) + "\n")
    box = _box(raw_names, None, box_id=f"box{candidate_rank:02d}", candidate_rank=candidate_rank)
    box["poses_by_seed"][7] = read_poses(path, raw_names)
    return box


@pytest.mark.parametrize("scope", ["full", "core"])
def test_raw_pdbqt_to_pool_to_full_or_core_scoring_is_permutation_invariant(tmp_path, scope):
    names, xyz, reference, sdf = _synthetic_parent()
    indices = list(range(len(names))) if scope == "full" else [7, 2, 0, 5]
    graph = v2.ParentGeometry(sdf, reference, {"indices": indices})
    native = xyz[indices]
    keys = [("A", 1, "", "ALA", atom) for atom in ("N", "CA", "C", "O")]
    protein = np.array([[0., 0, 3], [1., 0, 3], [0., 1, 3], [1., 1, 3]])
    fits, contacts, _ = v2._protein_alignment(keys, protein, keys, protein, native)
    ordered = _parsed_box(tmp_path, names, xyz, list(range(len(names))))
    shuffled = _parsed_box(tmp_path, names, xyz, list(reversed(range(len(names)))), candidate_rank=2)
    raw_before = shuffled["poses_by_seed"][7][0]["coordinates"].copy()
    identity = v2.canonicalize_pooled_poses(v2.pooled_poses([ordered], 7, 1), [ordered], names)
    permuted = v2.canonicalize_pooled_poses(v2.pooled_poses([shuffled], 7, 5), [shuffled], names)
    np.testing.assert_array_equal(identity[0]["coordinates"], xyz)
    np.testing.assert_array_equal(permuted[0]["coordinates"], xyz)
    np.testing.assert_array_equal(shuffled["poses_by_seed"][7][0]["coordinates"], raw_before)
    score = lambda rows: v2._score_pool(rows, native, protein, protein, keys, fits, contacts, graph)
    expected, actual = score(identity)[0], score(permuted)[0]
    for field in ("global_RMSD_A", "local_interface_RMSD_A"):
        assert expected[field] == actual[field]
        assert actual[field] < 1e-12
    for field in ("global_contacts", "local_interface_contacts", "geometry_diagnostics",
                  "full_parent_min_protein_distance_A", "full_parent_protein_pairs_below_1p5A"):
        assert expected[field] == actual[field]
    translated = _parsed_box(tmp_path, names, xyz, list(reversed(range(len(names)))), shift=[3, 4, 0])
    pool = v2.canonicalize_pooled_poses(v2.pooled_poses([translated], 7, 1), [translated], names)
    assert score(pool)[0]["global_RMSD_A"] == pytest.approx(5., abs=1e-12)
    assert len(permuted[0]["coordinates"]) == len(names) > 4


def test_joint_receptor_transform_preserves_named_score_without_ligand_fit():
    names = ["a", "b", "c"]
    xyz = np.array([[0., 0, 0], [2., 0, 0], [0., 1, 0]])
    protein = np.array([[0., 0, 2], [1., 0, 2], [0., 1, 2], [1., 1, 2]])
    rotation = np.array([[0., 1, 0], [-1., 0, 0], [0., 0, 1]])
    shift = np.array([8., -2, 5])
    source = protein @ rotation + shift
    keys = [("A", 1, "", "ALA", n) for n in ["N", "CA", "C", "O"]]
    fits, contacts, _ = v2._protein_alignment(keys, source, keys, protein, xyz)
    box = _box(["c", "a", "b"], (xyz @ rotation + shift)[[2, 0, 1]])
    class Graph:
        scored_maps = [[0, 1, 2]]
        def diagnostics(self, pose):
            return {"all_parent_atoms": len(pose)}
    pool = v2.canonicalize_pooled_poses(v2.pooled_poses([box], 7, 1), [box], names)
    result = v2._score_pool(pool, xyz, source, protein, keys, fits, contacts, Graph())[0]
    assert result["global_RMSD_A"] < 1e-12
    assert result["local_interface_RMSD_A"] < 1e-12
    assert result["global_contacts"]["contact_recall"] == 1


def test_two_box_orders_keep_ranking_and_never_double_permute():
    names = ["a", "b", "c"]
    xyz = np.arange(9, dtype=float).reshape(3, 3)
    boxes = [_box(["c", "a", "b"], xyz[[2, 0, 1]]),
             _box(["b", "c", "a"], xyz[[1, 2, 0]], box_id="box02", candidate_rank=2)]
    boxes[0]["poses_by_seed"][7].append(_pose(xyz[[2, 0, 1]], rank=2, score=-4))
    before = [p["coordinates"].copy() for b in boxes for p in b["poses_by_seed"][7]]
    first = v2.canonicalize_pooled_poses(v2.pooled_poses(boxes, 7, 1), boxes, names)
    raw_pool = v2.pooled_poses(boxes, 7, 5)
    second = v2.canonicalize_pooled_poses(raw_pool, boxes, names)
    third = v2.canonicalize_pooled_poses(v2.pooled_poses(boxes, 7, 5), boxes, names)
    assert [(r["box_id"], r["rank"], r["pooled_rank"]) for r in second] == [("box01", 1, 1), ("box02", 1, 2), ("box01", 2, 3)]
    for original, corrected in zip(raw_pool, second):
        assert {k: v for k, v in original.items() if k != "coordinates"} == {k: v for k, v in corrected.items() if k != "coordinates"}
    for row in first + second + third:
        np.testing.assert_array_equal(row["coordinates"], xyz)
    for original, current in zip(before, [p["coordinates"] for b in boxes for p in b["poses_by_seed"][7]]):
        np.testing.assert_array_equal(original, current)
    assert not np.shares_memory(second[0]["coordinates"], boxes[0]["poses_by_seed"][7][0]["coordinates"])


@pytest.mark.parametrize("actual,canonical", [
    (["a", "a"], ["a", "b"]), (["a", "b"], ["a", "a"]),
    (["a"], ["a", "b"]), (["a", "b", "c"], ["a", "b"]),
    (["a", "x"], ["a", "b"]), ("ab", ["a", "b"]),
    (["a", 2], ["a", "b"]), (["a", ""], ["a", "b"]),
    (["a", "b"], ("a", "b")), ([], []), (None, ["a", "b"]),
])
def test_invalid_full_parent_name_bijections_fail_closed(actual, canonical):
    box = _box(actual, [[0, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError):
        v2.canonicalize_pooled_poses(v2.pooled_poses([box], 7, 1), [box], canonical)


@pytest.mark.parametrize("coordinates", [
    [[0, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]], [0, 1],
    [[0, 0], [1, 0]], [[0, 0, 0], [float("nan"), 0, 0]],
    [[0, 0, 0], [float("inf"), 0, 0]], [["0", "0", "0"], ["1", "0", "0"]],
    [[True, False, False], [False, True, False]], [[0, 0, 0], [1, 0]], None,
])
def test_invalid_coordinate_types_shapes_lengths_fail_closed(coordinates):
    box = _box(["b", "a"], coordinates)
    with pytest.raises(ValueError):
        v2.canonicalize_pooled_poses(v2.pooled_poses([box], 7, 1), [box], ["a", "b"])


def test_unknown_and_duplicate_boxes_fail_closed_and_identity_returns_copy():
    box = _box(["a", "b"], np.array([[0., 0, 0], [1., 0, 0]]))
    rows = v2.pooled_poses([box], 7, 1)
    with pytest.raises(ValueError, match="unknown"):
        v2.canonicalize_pooled_poses([{**rows[0], "box_id": "unknown"}], [box], box["names"])
    with pytest.raises(ValueError, match="Duplicate"):
        v2.canonicalize_pooled_poses(rows, [box, box], box["names"])
    corrected = v2.canonicalize_pooled_poses(rows, [box], box["names"])[0]
    np.testing.assert_array_equal(corrected["coordinates"], rows[0]["coordinates"])
    assert not np.shares_memory(corrected["coordinates"], rows[0]["coordinates"])


def test_permutation_provenance_keeps_both_orders_and_exact_source_pins(tmp_path):
    raw = b'{"names": ["b", "a"]}'
    compressed = gzip.compress(raw, mtime=0)
    directory = tmp_path / "box/inputs"
    directory.mkdir(parents=True)
    (directory / "reference.json.gz").write_bytes(compressed)
    write_json(directory / "manifest.json", {"files": [{"filename": "reference.json", "archive": "reference.json.gz",
        "sha256": hashlib.sha256(raw).hexdigest(), "archive_sha256": hashlib.sha256(compressed).hexdigest()}]})
    source = v2._atom_order_reference(v2.Evidence(tmp_path), "box", raw)
    box = _box(["b", "a"], [[1, 0, 0], [0, 0, 0]])
    box["atom_order_reference"] = source
    canonical = {"archive": "legacy/reference.json.gz", "archive_sha256": "c" * 64, "reference_sha256": "d" * 64}
    record = v2._atom_order_provenance([box], ["a", "b"], canonical)[0]
    assert record["raw_box_names"] == ["b", "a"] and record["canonical_parent_names"] == ["a", "b"]
    assert record["canonical_to_raw_indices"] == [1, 0]
    assert record["raw_reference"] == source and record["canonical_reference"] == canonical
    with pytest.raises(ValueError, match="archived bytes"):
        v2._atom_order_reference(v2.Evidence(tmp_path), "box", raw + b" ")


def _pending_fixture(root):
    (root / "benchmarks").mkdir()
    (root / "results/blind_pose").mkdir(parents=True)
    write_json(root / "benchmarks/ledger.json", {})
    write_json(root / "benchmarks/blind_pose_experiment.json", {"frozen_before_new_docking": True,
        "atomic_reference_ledger": "benchmarks/ledger.json", "atomic_reference_ledger_sha256": sha256(root / "benchmarks/ledger.json")})
    for name in ("planning_inventory", "oracle_controls"):
        write_json(root / f"results/blind_pose/{name}.json", {})
    freeze = {"frozen_before_new_docking": True, "new_pose_outcomes_consulted": False,
        "configuration_sha256": sha256(root / "benchmarks/blind_pose_experiment.json"),
        "atomic_reference_ledger_sha256": sha256(root / "benchmarks/ledger.json"),
        "planning_inventory_sha256": sha256(root / "results/blind_pose/planning_inventory.json"),
        "oracle_controls_sha256": sha256(root / "results/blind_pose/oracle_controls.json")}
    write_json(root / "results/blind_pose/freeze_receipt.json", freeze)
    return freeze


@pytest.mark.parametrize("pending,terminal", [(["box:running"], 155), ([], 155)])
def test_complete_156_gate_precedes_amendment_native_context_and_permutation(tmp_path, pending, terminal):
    _pending_fixture(tmp_path)
    with patch.object(v2, "_check_atomic_sources"), patch.object(v2, "_all_execution", return_value=([], pending, {"terminal_seed_runs": terminal})), patch.object(v2, "_implementation_amendment", side_effect=AssertionError("amendment before terminal gate")), patch.object(v2, "_native_protein", side_effect=AssertionError("native context before terminal gate")), patch.object(v2, "canonicalize_pooled_poses", side_effect=AssertionError("permutation before terminal gate")):
        with pytest.raises(ValueError, match="incomplete|not terminal"):
            v2.generate_qualified_pose_report(tmp_path)
        if pending:
            summary, _ = v2.generate_qualified_pose_report(tmp_path, allow_pending=True)
            assert summary["metrics_eligible"] is False and "conditions" not in summary


def _amendment(freeze):
    return {"schema_version": 1, "implementation": "qualified_pose_v2", "evaluator_sha256": sha256(Path(v2.__file__)),
        "superseded_evaluator_sha256": v2._V1_EVALUATOR_SHA256, "correction": "bijective_atom_name_coordinate_permutation",
        "outcome_timing": "after_invalid_v1_and_corrected_diagnostic_outcomes", "invalid_v1_outcomes_consulted": True,
        "corrected_diagnostic_outcomes_consulted": True, "configuration_sha256": freeze["configuration_sha256"],
        "atomic_reference_ledger_sha256": freeze["atomic_reference_ledger_sha256"]}


def test_explicit_amendment_binds_evaluator_sources_and_honest_post_outcome_timing(tmp_path):
    freeze = {"configuration_sha256": "a" * 64, "atomic_reference_ledger_sha256": "b" * 64}
    path = tmp_path / "amendment.json"
    good = _amendment(freeze)
    write_json(path, good)
    checked = v2._implementation_amendment(v2.Evidence(tmp_path), "amendment.json", sha256(path), freeze)
    assert checked["record"] == good
    for field, value in [("evaluator_sha256", "c" * 64), ("outcome_timing", "before_outcomes"),
                         ("corrected_diagnostic_outcomes_consulted", False), ("invalid_v1_outcomes_consulted", 1),
                         ("configuration_sha256", "d" * 64), ("atomic_reference_ledger_sha256", "e" * 64)]:
        write_json(path, {**good, field: value})
        with pytest.raises(ValueError, match="post-outcome"):
            v2._implementation_amendment(v2.Evidence(tmp_path), "amendment.json", sha256(path), freeze)
    with pytest.raises(ValueError, match="explicit"):
        v2._implementation_amendment(v2.Evidence(tmp_path), None, None, freeze)
    with pytest.raises(ValueError, match="checksum"):
        v2._implementation_amendment(v2.Evidence(tmp_path), "amendment.json", "0" * 64, freeze)


def test_check_uses_frozen_path_comparator_and_protects_permutation_metadata(tmp_path):
    summary = {"conditions": [{"site_budgets": {"1": {"seed_evaluations": [{"poses": [{"global_RMSD_A": 2.123456789,
        "vina_score": -5., "rank": 1}]}]}}, "atom_order_correction": [{"canonical_to_raw_indices": [1, 0]}]}]}
    with patch.object(v2, "generate_qualified_pose_report", return_value=(summary, "V2 synthetic report\n")):
        paths = v2.write_qualified_pose_report(tmp_path)
        assert paths[0] == tmp_path / "results/development_pose_report_release/qualified_pose_v2_summary.json"
        assert paths[1].name == "QUALIFIED_BLIND_POSE_V2_RESULTS.md"
        rounded = copy.deepcopy(summary)
        rounded["conditions"][0]["site_budgets"]["1"]["seed_evaluations"][0]["poses"][0]["global_RMSD_A"] += 1e-12
        write_json(paths[0], rounded)
        with pytest.raises(ValueError, match="stale"):
            v2.write_qualified_pose_report(tmp_path, check=True)
        v2.write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)
        for field in ("permutation", "score", "source"):
            changed = copy.deepcopy(rounded)
            condition = changed["conditions"][0]
            if field == "permutation":
                condition["atom_order_correction"][0]["canonical_to_raw_indices"][0] += 1e-12
            elif field == "score":
                condition["site_budgets"]["1"]["seed_evaluations"][0]["poses"][0]["vina_score"] += 1e-12
            else:
                condition["atom_order_correction"][0]["raw_reference_sha256"] = "changed"
            write_json(paths[0], changed)
            with pytest.raises(ValueError):
                v2.write_qualified_pose_report(tmp_path, check=True, numeric_tolerance=1e-10)


def test_frozen_v1_and_scientific_endpoint_functions_remain_byte_exact():
    old_path = Path(v2.__file__).with_name("qualified_pose.py")
    assert sha256(old_path) == v2._V1_EVALUATOR_SHA256
    old = {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(old_path.read_text()).body if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
    new = {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(Path(v2.__file__).read_text()).body if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
    for name in ("pooled_poses", "ParentGeometry", "_score_pool", "qualified_atom_map", "native_metric_controls",
                 "accuracy_summary", "aggregate_conditions", "_protein_alignment", "_validate_command", "_pose_log"):
        assert old[name] == new[name], name


def test_geometry_diagnostic_baseline_is_explicitly_archived_sdf_not_box_pdbqt(tmp_path):
    sdf = b"synthetic independent SDF byte source\n"
    compressed = gzip.compress(sdf, mtime=0)
    directory = tmp_path / "legacy/inputs"
    directory.mkdir(parents=True)
    (directory / "ligand.sdf.gz").write_bytes(compressed)
    write_json(directory / "manifest.json", {"files": [{"filename": "ligand.sdf", "archive": "ligand.sdf.gz",
        "sha256": hashlib.sha256(sdf).hexdigest(), "archive_sha256": hashlib.sha256(compressed).hexdigest()}]})
    baseline = v2._geometry_baseline(v2.Evidence(tmp_path), "legacy", sdf)
    assert baseline["source"] == "archived_independent_legacy_SDF"
    assert baseline["exact_executed_box_PDBQT_baseline"] is False
    assert baseline["archive"] == "legacy/inputs/ligand.sdf.gz"
    assert baseline["uncompressed_sha256"] == hashlib.sha256(sdf).hexdigest()
    assert baseline["archive_sha256"] == hashlib.sha256(compressed).hexdigest()
    box = _box(["a", "b"], [[0, 0, 0], [1, 0, 0]])
    assert v2._atom_order_provenance([box], box["names"], {}, baseline)[0]["geometry_diagnostic_baseline"] == baseline
    with pytest.raises(ValueError, match="archived SDF"):
        v2._geometry_baseline(v2.Evidence(tmp_path), "legacy", sdf + b"changed")
