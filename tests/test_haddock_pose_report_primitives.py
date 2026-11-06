"""Synthetic HADDOCK report controls; no native calls or actual pose outcomes."""
from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from aav_benchmark.haddock_pose_report import parse_model, verify_stage_inventory


def pdb_atom(serial, name, residue, number, xyz, *, chain="A", segment="", element="C", alternate="", insertion="", record="ATOM"):
    return (f"{record:<6s}{serial:5d} {name:>4s}{alternate or ' '}{residue:>3s} {chain or ' '}{number:4d}{insertion or ' '}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.:6.2f}{0.:6.2f}      {segment:<4s}{element:>2s}  ")


def model_fixture():
    parent = [
        {"canonical_atom_key": "1:MAN:C1", "name": "L000"},
        {"canonical_atom_key": "1:MAN:O5", "name": "L001"},
        {"canonical_atom_key": "2:MAN:C1", "name": "L002"},
    ]
    cns = [
        {"canonical_atom_key": "2:MAN:C1", "cns_chain": "L", "cns_resnum": "2", "cns_resname": "MAN", "cns_atom_name": "C1", "element": "C"},
        {"canonical_atom_key": "1:MAN:C1", "cns_chain": "L", "cns_resnum": "1", "cns_resname": "MAN", "cns_atom_name": "C1", "element": "C"},
        {"canonical_atom_key": "1:MAN:O5", "cns_chain": "L", "cns_resnum": "1", "cns_resname": "MAN", "cns_atom_name": "O5", "element": "O"},
    ]
    protein = [("A", 1, "", "ALA", name) for name in ("N", "CA", "C", "O")]
    protein_xyz = np.asarray([[10., 0, 0], [11., 0, 0], [10., 1., 0], [10., 0, 1.]])
    lines = [pdb_atom(i + 1, key[-1], key[-2], key[1], xyz, chain=key[0], segment=key[0], element="N" if key[-1] == "N" else "O" if key[-1] == "O" else "C")
             for i, (key, xyz) in enumerate(zip(protein, protein_xyz))]
    # Deliberately different file, CNS-map and canonical-parent orders.
    lines.extend([
        pdb_atom(5, "O5", "MAN", 1, [2., 0, 0], chain="L", segment="L", element="O", record="HETATM"),
        pdb_atom(6, "C1", "MAN", 2, [3., 0, 0], chain="L", segment="L", record="HETATM"),
        pdb_atom(7, "C1", "MAN", 1, [1., 0, 0], chain="L", segment="L", record="HETATM"),
        pdb_atom(8, "H1", "MAN", 1, [1.1, 0, 0], chain="L", segment="L", element="H", record="HETATM"),
        pdb_atom(9, "D", "ALA", 1, [10., 0, .1], chain="A", segment="A", element="D"),
    ])
    return lines, cns, parent, protein, protein_xyz


def parse(lines, cns, parent, protein):
    return parse_model(("\n".join(lines) + "\nTER\nEND\n").encode(), cns, parent, protein)


def test_canonical_ccd_mapping_reorders_full_parent_and_retains_model_protein_coordinates():
    lines, cns, parent, protein, xyz = model_fixture()
    result = parse(lines, cns, parent, protein)
    np.testing.assert_array_equal(result["ligand_coordinates"], [[1., 0, 0], [2., 0, 0], [3., 0, 0]])
    observed = dict(zip(result["protein_keys"], result["protein_coordinates"]))
    assert set(observed) == set(protein)
    for key, point in zip(protein, xyz):
        np.testing.assert_array_equal(observed[key], point)
    assert result["total_atom_count"] == 9  # H and D counted, excluded from heavy endpoints.
    assert result["protein_coordinates"].shape == (4, 3)
    assert result["ligand_coordinates"].shape == (3, 3)
    shifted = lines.copy()
    for i, point in enumerate(xyz):
        line = shifted[i]
        moved = point + [20., -30., 40.]
        shifted[i] = line[:30] + "".join(f"{v:8.3f}" for v in moved) + line[54:]
    own_model = parse(shifted, cns, parent, protein)
    actual = dict(zip(own_model["protein_keys"], own_model["protein_coordinates"]))
    for key, point in zip(protein, xyz):
        np.testing.assert_array_equal(actual[key], point + [20., -30., 40.])
    # The parser neither fits the ligand nor substitutes the initial receptor.
    np.testing.assert_array_equal(own_model["ligand_coordinates"], result["ligand_coordinates"])


@pytest.mark.parametrize("kind", ["missing_ligand", "extra_ligand", "missing_protein", "extra_protein", "duplicate_heavy", "duplicate_hydrogen"])
def test_atom_loss_extra_heavy_and_all_atom_duplicates_fail(kind):
    lines, cns, parent, protein, _ = model_fixture()
    if kind == "missing_ligand":
        lines.pop(5)
    elif kind == "extra_ligand":
        lines.append(pdb_atom(10, "O9", "MAN", 1, [0, 0, 0], chain="L", segment="L", element="O", record="HETATM"))
    elif kind == "missing_protein":
        lines.pop(0)
    elif kind == "extra_protein":
        lines.append(pdb_atom(10, "CB", "ALA", 1, [0, 0, 0], chain="A", segment="A"))
    elif kind == "duplicate_heavy":
        lines.append(lines[0][:6] + f"{10:5d}" + lines[0][11:])
    else:
        lines.append(lines[7][:6] + f"{10:5d}" + lines[7][11:])
    with pytest.raises(ValueError):
        parse(lines, cns, parent, protein)


@pytest.mark.parametrize("index", [0, 4, 7])
def test_alternates_are_rejected_even_for_hydrogen(index):
    lines, cns, parent, protein, _ = model_fixture()
    lines[index] = lines[index][:16] + "A" + lines[index][17:]
    with pytest.raises(ValueError):
        parse(lines, cns, parent, protein)


@pytest.mark.parametrize("index,value", [(0, float("nan")), (4, float("inf")), (7, float("nan"))])
def test_nonfinite_coordinates_are_rejected_before_hydrogen_filtering(index, value):
    lines, cns, parent, protein, _ = model_fixture()
    lines[index] = lines[index][:30] + f"{value:8.3f}" + lines[index][38:]
    with pytest.raises(ValueError):
        parse(lines, cns, parent, protein)


def test_chain_segment_agreement_and_single_character_fallback():
    lines, cns, parent, protein, _ = model_fixture()
    # Explicit segment identity is usable when the PDB chain column is blank.
    for index in (0, 4, 5, 6):
        lines[index] = lines[index][:21] + " " + lines[index][22:]
    result = parse(lines, cns, parent, protein)
    assert set(result["protein_keys"]) == set(protein)
    np.testing.assert_array_equal(result["ligand_coordinates"][:, 0], [1., 2., 3.])
    conflict = model_fixture()[0]
    conflict[0] = conflict[0][:72] + "B   " + conflict[0][76:]
    with pytest.raises(ValueError):
        parse(conflict, cns, parent, protein)
    long_segment = model_fixture()[0]
    long_segment[0] = long_segment[0][:21] + " " + long_segment[0][22:72] + "LONG" + long_segment[0][76:]
    with pytest.raises(ValueError):
        parse(long_segment, cns, parent, protein)


def test_full_parent_mapping_and_element_identity_cannot_be_relaxed_to_core():
    lines, cns, parent, protein, _ = model_fixture()
    absent = copy.deepcopy(cns[:-1])
    with pytest.raises(ValueError):
        parse(lines, absent, parent, protein)
    duplicated = copy.deepcopy(cns)
    duplicated[1]["canonical_atom_key"] = duplicated[0]["canonical_atom_key"]
    with pytest.raises(ValueError):
        parse(lines, duplicated, parent, protein)
    changed_element = lines.copy()
    changed_element[6] = changed_element[6][:76] + " N" + changed_element[6][78:]
    with pytest.raises(ValueError):
        parse(changed_element, cns, parent, protein)


def model(stage, number, score, *, job="strict/synthetic/fpocket/box01/s104729"):
    name = f"{stage.split('_', 1)[1]}_{number}.pdb" if number is not None else "flexref_unknown.pdb"
    return {"model_id": f"{job}/{stage}/{name}", "original_model_number": number, "file_name": name,
            "sha256": "a" * 64, "atom_count": 9, "score": score, "energies": {"vdw": score},
            "actual_CNS_seed": 104729 + (number or 0), "cluster_id": None, "cluster_rank": None,
            "within_cluster_rank": None}


def stage(name, models):
    rows = copy.deepcopy(models)
    ranked = sorted([row for row in rows if row["score"] is not None and math.isfinite(row["score"])],
                    key=lambda row: (row["score"], row["original_model_number"] or 0, row["model_id"]))
    for rank, row in enumerate(ranked, 1):
        row["individual_score_rank"] = rank
    return {"stage": name, "actual_physical_model_files": len(rows), "models": rows, "io_sha256": "b" * 64,
            "individual_score_order": [row["model_id"] for row in ranked]}


def stages_fixture():
    topo = stage("0_topoaa", [model("0_topoaa", 1, None), model("0_topoaa", 2, None)])
    initial = stage("1_rigidbody", [model("1_rigidbody", 1, -1.1), model("1_rigidbody", 2, -1.2), model("1_rigidbody", 3, -1.3)])
    final = stage("2_flexref", [model("2_flexref", 10, -2.123456789), model("2_flexref", 2, -2.123456788), model("2_flexref", 1, -.1)])
    return [topo, initial, final]


def chosen_final(value):
    if isinstance(value, dict):
        if str(value.get("stage", "")).endswith("_flexref"):
            return value
        for nested in value.values():
            found = chosen_final(nested)
            if found is not None:
                return found
    elif isinstance(value, list):
        for nested in value:
            found = chosen_final(nested)
            if found is not None:
                return found
    return None


def test_stage_selection_retains_all_individual_models_and_full_precision_without_mutation():
    stages = stages_fixture()
    before = copy.deepcopy(stages)
    result = verify_stage_inventory(stages)
    assert isinstance(result, dict) and stages == before
    final = chosen_final(result)
    assert final is not None and final["actual_physical_model_files"] == len(final["models"]) == 3
    # A rounding-to-six-decimals tie would produce the wrong first model.
    assert final["individual_score_order"][0].endswith("flexref_10.pdb")
    assert next(m["score"] for m in final["models"] if m["original_model_number"] == 10) == -2.123456789
    tied = copy.deepcopy(stages)
    tied[-1] = stage("2_flexref", [model("2_flexref", 10, -2.), model("2_flexref", 2, -2.), model("2_flexref", None, -2.)])
    tie_final = chosen_final(verify_stage_inventory(tied))
    assert [path.rsplit("/", 1)[-1] for path in tie_final["individual_score_order"]] == ["flexref_unknown.pdb", "flexref_2.pdb", "flexref_10.pdb"]


@pytest.mark.parametrize("kind", ["missing_initial", "missing_final", "two_initial", "two_final", "count", "over_cap", "missing_score", "nan_score", "infinite_score", "wrong_order", "wrong_rank", "duplicate_id"])
def test_invalid_or_incomplete_stage_inventory_fails(kind):
    stages = stages_fixture()
    if kind == "missing_initial":
        stages.pop(1)
    elif kind == "missing_final":
        stages.pop()
    elif kind == "two_initial":
        extra = copy.deepcopy(stages[1]); extra["stage"] = "11_rigidbody"; stages.append(extra)
    elif kind == "two_final":
        extra = copy.deepcopy(stages[2]); extra["stage"] = "12_flexref"; stages.append(extra)
    elif kind == "count":
        stages[-1]["actual_physical_model_files"] = 2
    elif kind == "over_cap":
        stages[-1] = stage("2_flexref", [model("2_flexref", i, -float(i)) for i in range(1, 22)])
    elif kind == "missing_score":
        stages[-1]["models"][0]["score"] = None
    elif kind == "nan_score":
        stages[-1]["models"][0]["score"] = float("nan")
    elif kind == "infinite_score":
        stages[-1]["models"][0]["score"] = float("inf")
    elif kind == "wrong_order":
        stages[-1]["individual_score_order"] = stages[-1]["individual_score_order"][::-1]
    elif kind == "wrong_rank":
        stages[-1]["models"][0]["individual_score_rank"] = 3
    else:
        duplicate = copy.deepcopy(stages[-1]["models"][0])
        duplicate["file_name"] += ".gz"
        duplicate["sha256"] = "c" * 64
        stages[-1]["models"].append(duplicate)
        stages[-1]["actual_physical_model_files"] += 1
        stages[-1]["individual_score_order"].append(duplicate["model_id"])
    with pytest.raises(ValueError):
        verify_stage_inventory(stages)
