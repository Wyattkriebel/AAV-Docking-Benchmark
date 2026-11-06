"""Source-only monomer controls; no native predictor, modeling, or docking."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tarfile

import gemmi
import pytest

from aav_benchmark.heldout_assembly_control import (
    CASE_IDS, LIGAND_FILES, REQUIRED_FILES, build_assembly_control,
    materialize_control_case, subset_monomer_input,
)
from aav_benchmark.heldout_serialization import atom_semantics, compare_atom_semantics, finalize_cif_metadata

ROOT = Path(__file__).resolve().parents[1]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def multichain_cif(tmp_path):
    document = gemmi.cif.Document()
    category = {field: [] for field in (
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
        "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
        "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
        "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num",
    )}
    # Z appears first and is much closer to the origin; neither ordering nor
    # spatial proximity is permitted to override lexical chain A selection.
    for serial, chain, number, insertion, resname, atom_name, element, x in (
        (1, "Z", "1", "?", "GLY", "CA", "C", "0.000"),
        (2, "A", "-2", "B", "MSE", "SE", "Se", "91.123450"),
        (3, "A", "17", "?", "GLY", "CA", "C", "95.765400"),
        (4, "B", "4", "?", "GLY", "CA", "C", "3.000"),
    ):
        values = ("HETATM" if resname == "MSE" else "ATOM", str(serial), element, atom_name, False,
                  resname, chain, False, False, insertion, x, "2.123450", "-3.500000",
                  "0.5000", "19.123400", number, chain, "1")
        for field, value in zip(category, values):
            category[field].append(value)
    document.add_new_block("source_multichain").set_mmcif_category("_atom_site.", category)
    source = tmp_path / "source.cif"
    source.write_text(document.as_string())
    return source


def test_subset_selects_lexical_author_chain_and_preserves_mse_and_every_atom_field(multichain_cif, tmp_path):
    before_bytes = multichain_cif.read_bytes()
    output = tmp_path / "monomer.cif"
    proof = subset_monomer_input(multichain_cif, output)
    assert proof["selected_author_chain"] == "A"
    assert proof["source_author_chains"] == ["A", "B", "Z"]
    assert proof["source_assembly_atom_count"] == 4
    assert proof["selected_atom_count"] == 2
    assert proof["omitted_author_chains"] == ["B", "Z"]
    assert proof["native_ligand_coordinate_access"] is False
    assert proof["distance_based_selection"] is False
    assert proof["physically_incomplete_assembly_control"] is True
    before = atom_semantics(multichain_cif)
    expected = {key: row for key, row in before.items() if key[1] == "A"}
    assert compare_atom_semantics(expected, atom_semantics(output))["all_nonlabel_atom_fields_exactly_unchanged"] is True
    assert multichain_cif.read_bytes() == before_bytes
    assert output.read_bytes().endswith(b"\n#\n")
    atoms = gemmi.cif.read_file(str(output)).sole_block().get_mmcif_category("_atom_site.")
    assert atoms["type_symbol"] == ["Se", "C"]
    assert atoms["occupancy"] == ["0.5000", "0.5000"]
    assert atoms["B_iso_or_equiv"] == ["19.123400", "19.123400"]
    assert atoms["auth_seq_id"] == ["-2", "17"]
    assert atoms["pdbx_PDB_ins_code"][0] == "B"


def test_subset_cannot_overwrite_existing_input(multichain_cif, tmp_path):
    output = tmp_path / "retained.cif"
    output.write_text("preserve previous control\n")
    with pytest.raises(FileExistsError):
        subset_monomer_input(multichain_cif, output)
    assert output.read_text() == "preserve previous control\n"


def test_subset_rejects_already_monomeric_source(multichain_cif, tmp_path):
    first = tmp_path / "first.cif"
    subset_monomer_input(multichain_cif, first)
    with pytest.raises(ValueError, match="multiple"):
        subset_monomer_input(first, tmp_path / "second.cif")


def _actual_case(case_id, temporary):
    directory = temporary / case_id
    directory.mkdir()
    with tarfile.open(ROOT / "results/structural_test_v2/source_preparation.tar.gz") as bundle:
        for name in REQUIRED_FILES:
            (directory / name).write_bytes(bundle.extractfile(f"{case_id}/{name}").read())
    return directory


@pytest.mark.parametrize("case_id", CASE_IDS)
def test_actual_controls_retain_exact_subset_and_unchanged_full_source_maps_and_ligands(case_id, tmp_path):
    primary = _actual_case(case_id, tmp_path)
    hashes = {name: _sha(primary / name) for name in REQUIRED_FILES}
    output = tmp_path / "control"
    receipt = materialize_control_case(case_id, primary, output)
    assert receipt["status"] == "prepared"
    assert receipt["subset_proof"]["selected_author_chain"] == "A"
    assert receipt["independent_ligand_bytes_and_map_unchanged"] is True
    for name in LIGAND_FILES:
        assert _sha(output / name) == hashes[name]
    assert _sha(output / "full_assembly_receptor.cif") == hashes["receptor.cif"]
    assert _sha(output / "full_protein_atom_map.json") == hashes["protein_atom_map.json"]
    full = json.loads((primary / "protein_atom_map.json").read_text())
    subset = json.loads((output / "protein_atom_map.json").read_text())
    assert subset["atoms"] == [row for row in full["atoms"] if row["output_key"][0] == "A"]
    assert subset["chains"] == [row for row in full["chains"] if row["output_chain"] == "A"]
    expected = {key: row for key, row in atom_semantics(primary / "receptor.cif").items() if key[1] == "A"}
    assert len(expected) < len(atom_semantics(primary / "receptor.cif"))
    assert compare_atom_semantics(expected, atom_semantics(output / "receptor.cif"))["max_coordinate_displacement_A"] == 0.
    assert receipt["subset_proof"]["selected_atom_count"] == len(subset["atoms"]) == len(expected)
    assert all(_sha(primary / name) == expected_hash for name, expected_hash in hashes.items())
    manifest = json.loads((output / "receptor.cif.manifest.json").read_text())
    assert manifest["prediction_scope_does_not_change_atomic_truth"] is True
    assert manifest["full_assembly_evaluation_receptor"] == "full_assembly_receptor.cif"


def test_fixed_case_scope_cannot_expand_from_input_availability(tmp_path):
    with pytest.raises(ValueError, match="predeclared"):
        materialize_control_case("2RDK", tmp_path, tmp_path / "new")


def test_control_build_proposal_preserves_full_assembly_evaluation_scope(tmp_path):
    primary = tmp_path / "primary"
    primary.mkdir()
    cases = []
    for case_id in CASE_IDS:
        source = _actual_case(case_id, tmp_path)
        directory = primary / case_id
        directory.mkdir()
        for name in REQUIRED_FILES:
            (directory / name).write_bytes((source / name).read_bytes())
        cases.append({"case_id": case_id, "artifact_sha256": {name: _sha(directory / name) for name in REQUIRED_FILES}})
    source_receipt = {"case_ids": ["2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U"], "serialization_revision": 2,
                      "source_ledger_sha256": "source_only_test_pin", "cases": cases}
    (primary / "source_preparation.json").write_text(json.dumps(source_receipt))
    output = tmp_path / "controls"
    result = build_assembly_control(primary, output)
    assert result["case_ids"] == list(CASE_IDS)
    assert result["fixed_control_denominator"] == 2
    assert result["prediction_outcomes_consulted"] is False
    proposal = json.loads((output / "proposal.json").read_text())
    assert proposal["status"] == "proposal_requires_separate_coordinator_freeze_before_new_predictions"
    assert proposal["prediction"]["p2rank_seed"] == 42
    assert proposal["candidate_budgets"] == [1, 5]
    assert proposal["receptor_policy"]["interaction_cutoff_A"] == 8.
    assert proposal["evaluation"]["clash_denominator"] == "unchanged full biological assembly, including omitted protein chains"
    assert "full-assembly source-qualified copies" in proposal["evaluation"]["atomic_reference"]
    assert "not biological nonbinding" in proposal["evaluation"]["known_sites"]
    assert proposal["maximum_method_conditions"] == 6
    assert proposal["maximum_box_conditions"] == 22
    assert proposal["maximum_seed_runs"] == 66
    with pytest.raises(FileExistsError):
        build_assembly_control(primary, output)
