import json
from pathlib import Path

import numpy as np
import pytest

from aav_benchmark.heldout_docking import (
    CASE_IDS, INTERACTION_CUTOFF_A, select_residue_altloc, check_signed_chirality,
    source_cases, build_source_ligand, prepare_source_protein, neighborhood_protein,
    certified_protein_operators, verify_freeze,
)

ROOT = Path(__file__).resolve().parents[1]


def test_whole_residue_alt_mean_and_tie():
    rows = [{"label_atom_id": "N", "label_alt_id": False, "occupancy": "1"},
            {"label_atom_id": "CB", "label_alt_id": "A", "occupancy": ".4"},
            {"label_atom_id": "CG", "label_alt_id": "A", "occupancy": ".8"},
            {"label_atom_id": "CB", "label_alt_id": "B", "occupancy": ".6"},
            {"label_atom_id": "CG", "label_alt_id": "B", "occupancy": ".6"}]
    selected, audit = select_residue_altloc(rows)
    assert audit["selected_alt_id"] == "A"
    assert {r["label_alt_id"] for r in selected} == {False, "A"}
    rows[-1]["occupancy"] = ".9"
    assert select_residue_altloc(rows)[1]["selected_alt_id"] == "B"


def test_signed_volume_rejects_reflection():
    atom_map = [{"canonical_atom_key": name} for name in ("center", "a", "b", "c")]
    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    tests = [{"center": "center", "neighbors": ["a", "b", "c"], "ideal_signed_volume": 1.}]
    assert check_signed_chirality(atom_map, xyz, tests)[0]["status"] == "consistent"
    xyz[:, 0] *= -1
    with pytest.raises(ValueError, match="chirality"):
        check_signed_chirality(atom_map, xyz, tests)


def test_fixed_parent_chemistry_is_not_scored_core():
    cases = source_cases(ROOT)
    assert tuple(c["case_id"] for c in cases) == CASE_IDS
    case = next(c for c in cases if c["case_id"] == "4R9F")
    molecule, atom_map, chemistry = build_source_ligand(ROOT, case)
    assert len(atom_map) == case["inventory_entry"]["expected_heavy_atom_count"] == 23
    assert len(atom_map) > case["scored_atom_count"]
    assert chemistry["native_coordinates_used_for_embedding"] is False
    assert all(t["status"] == "consistent" for t in chemistry["signed_chirality_tests"])
    assert molecule.GetNumHeavyAtoms() == 23


def test_protein_source_selection_retains_modification_and_copy_map(tmp_path):
    case = next(c for c in source_cases(ROOT) if c["case_id"] == "4R9F")
    result = prepare_source_protein(ROOT, case, tmp_path)
    assert result["protein_modifications_preserved"] == ["MSE"]
    assert result["distance_based_selection"] is False
    mapping = json.loads((tmp_path / "protein_atom_map.json").read_text())
    assert any(r["element"] == "Se" for r in mapping["atoms"])
    assert len(mapping["atoms"]) == result["protein_heavy_atom_count"]
    assert all("source_author_resnum" in r and "operator_name" in r for r in mapping["atoms"])
    operators, audit = certified_protein_operators(tmp_path / "receptor.cif", mapping["chains"])
    assert operators and audit["orbit_deduplication"] == "no_orbit_dedupe"


def test_freeze_requires_both_prediction_and_docking_barriers(tmp_path):
    config = tmp_path / "config.json"; config.write_text("{}")
    freeze = tmp_path / "freeze.json"; freeze.write_text("{}")
    with pytest.raises(ValueError, match="pre-prediction"):
        verify_freeze(ROOT, config, freeze)


def test_neighborhood_retains_chain_between_six_and_eight(tmp_path):
    import gemmi
    structure = gemmi.Structure(); structure.add_model(gemmi.Model("1"))
    for name, x in (("A", 0.), ("B", 17.), ("C", 19.)):
        chain = gemmi.Chain(name); residue = gemmi.Residue(); residue.name = "GLY"
        residue.seqid = gemmi.SeqId(1, " "); residue.het_flag = "A"
        atom = gemmi.Atom(); atom.name = "CA"; atom.element = gemmi.Element("C"); atom.pos = gemmi.Position(x, 0, 0)
        residue.add_atom(atom); chain.add_residue(residue); structure[0].add_chain(chain)
    path = tmp_path / "protein.cif"; structure.make_mmcif_document().write_file(str(path))
    _, audit = neighborhood_protein(path, [0, 0, 0], [20, 20, 20])
    assert INTERACTION_CUTOFF_A == audit["interaction_cutoff_A"] == 8
    assert audit["selected_chains"] == ["A", "B"]
    assert audit["closest_omitted_atom_to_box_A"] == 9.
