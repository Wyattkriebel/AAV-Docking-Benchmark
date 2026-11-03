"""Source preservation gates for the separately declared repair sensitivity.

All fixtures are synthetic.  These tests do not import OpenMM or run modeling.
"""
from copy import deepcopy

import numpy as np
import pytest

from aav_benchmark.heldout_repair import (
    audit_repair, parse_pdb_heavy, split_modeling_segments,
)


def _atom(serial, atom_name, residue_number, xyz, *, chain="A", resname="GLY",
          element="C", insertion="", alternate=""):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {atom_name:>4s}{alternate or ' ':1s}{resname:3s} "
            f"{chain:1s}{residue_number:4d}{insertion or ' ':1s}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}\n")


def _residue(start_serial, number, x, *, chain="A", resname="GLY"):
    return "".join(_atom(start_serial + i, name, number, (x + dx, 0., 0.),
                         chain=chain, resname=resname, element=element)
                   for i, (name, dx, element) in enumerate((
                       ("N", 0., "N"), ("CA", 1., "C"),
                       ("C", 2., "C"), ("O", 2.5, "O"))))


@pytest.fixture
def observed():
    return parse_pdb_heavy(_residue(1, 7, 0., resname="LYS"))


def _allowed_row(coordinates=(2., 1., 0.), element="C"):
    return {"coordinates": list(coordinates), "element": element, "line": ""}


def test_parser_preserves_source_author_ids_elements_and_coordinates():
    text = (_atom(1, "N", -2, (1.125, -2.25, 3.5), chain="Z", element="N", insertion="A")
            + _atom(2, "SE", 19, (4.125, 5.25, -6.5), chain="b", resname="MSE", element="Se"))
    parsed = parse_pdb_heavy(text)
    assert set(parsed) == {("Z", -2, "A", "GLY", "N"), ("b", 19, "", "MSE", "SE")}
    assert parsed[("Z", -2, "A", "GLY", "N")]["coordinates"] == [1.125, -2.25, 3.5]
    assert parsed[("b", 19, "", "MSE", "SE")]["element"] == "Se"
    assert parsed[("b", 19, "", "MSE", "SE")]["line"][12:16].strip() == "SE"


def test_parser_rejects_ambiguous_duplicate_heavy_atom_identity():
    row = _atom(1, "CA", 7, (0., 0., 0.))
    with pytest.raises(ValueError):
        parse_pdb_heavy(row + _atom(2, "CA", 7, (1., 0., 0.)))


def test_audit_accepts_only_allowlisted_additions_without_moving_source(observed):
    added_key = ("A", 7, "", "LYS", "CB")
    after = deepcopy(observed)
    after[added_key] = _allowed_row()
    audit = audit_repair(observed, after, {added_key})
    assert audit["observed_heavy_atoms"] == len(observed)
    assert audit["added_heavy_atoms"] == 1
    assert audit["max_observed_displacement_A"] == 0.
    assert len(audit["added_atoms"]) == 1
    assert after[added_key]["coordinates"] == [2., 1., 0.]


def test_audit_accepts_unchanged_protein_when_no_additions_are_expected(observed):
    audit = audit_repair(observed, deepcopy(observed), set())
    assert audit["added_heavy_atoms"] == 0
    assert audit["added_atoms"] == []


def test_audit_rejects_missing_expected_additions(observed):
    with pytest.raises(ValueError):
        audit_repair(observed, deepcopy(observed), {("A", 7, "", "LYS", "CB")})


def test_audit_rejects_deleting_an_observed_heavy_atom(observed):
    after = deepcopy(observed)
    del after[("A", 7, "", "LYS", "O")]
    with pytest.raises(ValueError):
        audit_repair(observed, after, set())


def test_audit_rejects_changing_observed_element(observed):
    after = deepcopy(observed)
    after[("A", 7, "", "LYS", "C")]["element"] = "Se"
    with pytest.raises(ValueError):
        audit_repair(observed, after, set())


def test_audit_rejects_an_unlisted_atom_even_with_other_permitted_atoms(observed):
    after = deepcopy(observed)
    after[("A", 7, "", "LYS", "CB")] = _allowed_row()
    after[("A", 7, "", "LYS", "CG")] = _allowed_row((2., 2., 0.))
    with pytest.raises(ValueError):
        audit_repair(observed, after, {("A", 7, "", "LYS", "CB")})


@pytest.mark.parametrize("added_key", [
    ("A", 7, "", "LYS", "OXT"),
    ("A", 7, "", "LYS", "CA"),
    ("A", 8, "", "MSE", "SE"),
])
def test_audit_allowlist_cannot_authorize_backbone_or_nonstandard_completion(observed, added_key):
    before = deepcopy(observed)
    # Make CA a missing observed backbone atom for this fixture, so it is an
    # addition rather than a modified existing atom.
    before.pop(added_key, None)
    after = deepcopy(before)
    after[added_key] = _allowed_row((2., 1., 0.))
    with pytest.raises(ValueError):
        audit_repair(before, after, {added_key})


@pytest.mark.parametrize("element", ["H", "D"])
def test_audit_rejects_hydrogen_even_if_allowlisted(observed, element):
    after = deepcopy(observed)
    key = ("A", 7, "", "LYS", "H1")
    after[key] = _allowed_row((0., 1., 0.), element=element)
    with pytest.raises(ValueError):
        audit_repair(observed, after, {key})


def test_audit_accepts_displacement_within_euclidean_tolerance(observed):
    after = deepcopy(observed)
    after[("A", 7, "", "LYS", "N")]["coordinates"] = [.0006, .0006, 0.]
    audit = audit_repair(observed, after, set(), max_observed_displacement_A=.001)
    assert audit["max_observed_displacement_A"] == pytest.approx(np.sqrt(2) * .0006)


def test_audit_rejects_combined_displacement_over_euclidean_tolerance(observed):
    after = deepcopy(observed)
    after[("A", 7, "", "LYS", "N")]["coordinates"] = [.0008, .0008, 0.]
    with pytest.raises(ValueError):
        audit_repair(observed, after, set(), max_observed_displacement_A=.001)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_audit_rejects_nonfinite_source_or_added_coordinates(observed, bad):
    after = deepcopy(observed)
    after[("A", 7, "", "LYS", "N")]["coordinates"][0] = bad
    with pytest.raises(ValueError):
        audit_repair(observed, after, set())
    key = ("A", 7, "", "LYS", "CB")
    after = deepcopy(observed)
    after[key] = _allowed_row((bad, 1., 0.))
    with pytest.raises(ValueError):
        audit_repair(observed, after, {key})


def _chain_residues(text):
    groups = {}
    for key in parse_pdb_heavy(text):
        groups.setdefault(key[0], set()).add(key[1:4])
    return groups


def test_segments_split_number_gaps_and_broken_peptides_without_source_loss():
    # 10->11 is a 1.3 A peptide. 11->13 has a numbering gap despite
    # a 1.3 A distance. 13->14 is consecutive but has a 3.0 A C-N distance.
    text = (_residue(1, 10, 0.) + _residue(5, 11, 3.3)
            + _residue(9, 13, 6.6) + _residue(13, 14, 11.6)
            + _residue(17, 1, 30., chain="B") + _residue(21, 2, 33.3, chain="B"))
    rewritten, segment_map = split_modeling_segments(text)
    groups = _chain_residues(rewritten)
    assert sorted(sorted(number for number, _, _ in residues) for residues in groups.values()) == [
        [1, 2], [10, 11], [13], [14]]
    assert len(segment_map["chains"]) == 4
    assert len(segment_map["residues"]) == 6
    before = parse_pdb_heavy(text)
    restored = {}
    for key, row in parse_pdb_heavy(rewritten).items():
        source_key = (segment_map["chains"][key[0]], *key[1:])
        assert source_key not in restored
        restored[source_key] = row
    audit = audit_repair(before, restored, set())
    assert audit["observed_heavy_atoms"] == len(before)
    assert audit["added_heavy_atoms"] == 0
    assert audit["max_observed_displacement_A"] == 0.


@pytest.mark.parametrize("cn_distance", [.9, 1.8])
def test_segments_reject_both_short_and_long_peptide_connections(cn_distance):
    text = _residue(1, 1, 0.) + _residue(5, 2, 2. + cn_distance)
    rewritten, _ = split_modeling_segments(text)
    assert len(_chain_residues(rewritten)) == 2


@pytest.mark.parametrize("cn_distance", [1.1, 1.3, 1.7])
def test_segments_keep_connected_residues_in_the_declared_peptide_range(cn_distance):
    text = _residue(1, 1, 0.) + _residue(5, 2, 2. + cn_distance)
    rewritten, _ = split_modeling_segments(text)
    assert len(_chain_residues(rewritten)) == 1
