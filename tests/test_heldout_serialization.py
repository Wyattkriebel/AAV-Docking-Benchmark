"""Metadata-only reader compatibility; no predictor or docking execution."""
from copy import deepcopy
import hashlib
from pathlib import Path
import tarfile

import gemmi
import pytest

from aav_benchmark.heldout_serialization import (
    atom_semantics, compare_atom_semantics, finalize_cif_metadata,
)

ROOT = Path(__file__).resolve().parents[1]
CASE_IDS = ("2RDK", "3ZWE", "4R9F", "5HZA", "5V6F", "8I5U")
LABEL_FIELDS = {"label_asym_id", "label_entity_id", "label_seq_id"}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _atom_category(path):
    return gemmi.cif.read_file(str(path)).sole_block().get_mmcif_category("_atom_site.")


def _write_category(path, category):
    doc = gemmi.cif.Document()
    block = doc.add_new_block("synthetic_metadata")
    block.set_mmcif_category("_atom_site.", category)
    path.write_text(doc.as_string())


@pytest.fixture
def unfinalized_cif(tmp_path):
    # Includes nonconsecutive author numbers, an insertion, MSE/Se, a HETATM
    # polymer record, and decimal precision that must survive unchanged.
    residues = [
        ("A", "-2", "A", "GLY", [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]),
        ("A", "9", "?", "MSE", [("N", "N"), ("CA", "C"), ("SE", "Se")]),
        ("A", "10", "?", "GLU", [("N", "N"), ("CA", "C")]),
        ("B", "41", "?", "GLY", [("N", "N"), ("CA", "C")]),
        ("B", "43", "?", "SER", [("N", "N"), ("CA", "C"), ("OG", "O")]),
    ]
    category = {field: [] for field in (
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
        "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
        "B_iso_or_equiv", "pdbx_formal_charge", "auth_atom_id", "auth_comp_id",
        "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num",
    )}
    for chain, number, insertion, resname, atoms in residues:
        for name, element in atoms:
            serial = len(category["id"]) + 1
            row = {
                "group_PDB": "HETATM" if resname == "MSE" else "ATOM",
                "id": str(serial), "type_symbol": element, "label_atom_id": name,
                "label_alt_id": ".", "label_comp_id": resname,
                "label_asym_id": "old" + chain, "label_entity_id": ".", "label_seq_id": ".",
                "pdbx_PDB_ins_code": insertion, "Cartn_x": f"{serial}.123450",
                "Cartn_y": f"-{serial}.765400", "Cartn_z": "0.000010",
                "occupancy": "0.5000" if serial % 2 else "1.0000",
                "B_iso_or_equiv": f"{serial}.123400", "pdbx_formal_charge": ".",
                "auth_atom_id": name, "auth_comp_id": resname,
                "auth_seq_id": number, "auth_asym_id": chain, "pdbx_PDB_model_num": "1",
            }
            for key in category:
                category[key].append(row[key])
    source = tmp_path / "prior.cif"
    _write_category(source, category)
    return source


def _assert_raw_atoms_unchanged(before, after):
    assert set(before) == set(after)
    for field in before:
        if field not in LABEL_FIELDS:
            assert after[field] == before[field], field


def _assert_entity_references_consistent(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = block.get_mmcif_category("_atom_site.")
    entities = block.get_mmcif_category("_entity.")
    asym = block.get_mmcif_category("_struct_asym.")
    poly = block.get_mmcif_category("_entity_poly.")
    sequence = block.get_mmcif_category("_entity_poly_seq.")
    entity_ids = set(entities["id"])
    assert entity_ids and all(str(value).isdigit() and int(value) > 0 for value in entity_ids)
    assert set(entities["type"]) == {"polymer"}
    assert set(poly["entity_id"]) == entity_ids
    assert set(poly["type"]) == {"polypeptide(L)"}
    assert all(value not in (None, False, "", "?", ".") for value in poly["pdbx_seq_one_letter_code"])
    assert set(asym["id"]) == set(atoms["label_asym_id"])
    asym_to_entity = dict(zip(asym["id"], asym["entity_id"]))
    for label_asym, entity_id, label_seq in zip(atoms["label_asym_id"], atoms["label_entity_id"], atoms["label_seq_id"]):
        assert entity_id in entity_ids
        assert asym_to_entity[label_asym] == entity_id
        assert str(label_seq).isdigit() and int(label_seq) > 0
    expected_residues = {}
    observed_sequence = {}
    for i in range(len(atoms["id"])):
        author = tuple(str(atoms[field][i]) for field in ("auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code")) + (
            str(atoms.get("auth_comp_id", atoms["label_comp_id"])[i]),)
        label = (atoms["label_entity_id"][i], atoms["label_seq_id"][i])
        assert expected_residues.setdefault(author, label) == label
        assert observed_sequence.setdefault(label, atoms["label_comp_id"][i]) == atoms["label_comp_id"][i]
    declared_sequence = dict(zip(zip(sequence["entity_id"], sequence["num"]), sequence["mon_id"]))
    assert declared_sequence == observed_sequence
    for entity_id in entity_ids:
        numbers = sorted(int(number) for entity, number in declared_sequence if entity == entity_id)
        assert numbers == list(range(1, len(numbers) + 1))


def test_finalization_preserves_all_physical_atom_fields_and_prior_bytes(unfinalized_cif, tmp_path):
    before_bytes = unfinalized_cif.read_bytes()
    before = _atom_category(unfinalized_cif)
    output = tmp_path / "v2.cif"
    receipt = finalize_cif_metadata(unfinalized_cif, output)
    assert unfinalized_cif.read_bytes() == before_bytes
    _assert_raw_atoms_unchanged(before, _atom_category(output))
    assert receipt["atom_count"] == len(before["id"])
    assert receipt["source_receptor_sha256"] == _sha(unfinalized_cif)
    assert receipt["receptor_sha256"] == _sha(output)
    assert receipt["source_receptor_sha256"] != receipt["receptor_sha256"]
    for field in ("author_atom_identities_unchanged", "elements_unchanged", "coordinates_exactly_unchanged", "metadata_scope_only"):
        assert receipt[field] is True
    assert receipt["max_coordinate_displacement_A"] == 0.
    assert receipt["serialization_sequence"] == "observed residue names only; not full unresolved construct"
    assert receipt["entity_count"] == 2
    assert {field for field, count in receipt["atom_field_change_counts"].items() if count} <= LABEL_FIELDS


def test_finalization_writes_complete_numeric_polymer_metadata_and_terminal_delimiter(unfinalized_cif, tmp_path):
    output = tmp_path / "v2.cif"
    finalize_cif_metadata(unfinalized_cif, output)
    assert output.read_bytes().endswith(b"\n#\n")
    _assert_entity_references_consistent(output)
    atoms = _atom_category(output)
    assert "MSE" in atoms["auth_comp_id"]
    assert "Se" in atoms["type_symbol"]
    assert ("-2", "A") in set(zip(atoms["auth_seq_id"], atoms["pdbx_PDB_ins_code"]))
    assert {"9", "43"} <= set(atoms["auth_seq_id"])


def test_finalization_preserves_existing_output(unfinalized_cif, tmp_path):
    output = tmp_path / "existing.cif"
    output.write_text("retain this previous input\n")
    with pytest.raises(FileExistsError):
        finalize_cif_metadata(unfinalized_cif, output)
    assert output.read_text() == "retain this previous input\n"


@pytest.mark.parametrize("field,new_value", [
    ("Cartn_x", "0.125"), ("type_symbol", "Se"),
    ("occupancy", "0.2500"), ("B_iso_or_equiv", "98.7654"),
    ("auth_seq_id", "99"), ("auth_asym_id", "Z"),
])
def test_semantic_comparison_rejects_physical_atom_or_author_changes(unfinalized_cif, tmp_path, field, new_value):
    category = deepcopy(_atom_category(unfinalized_cif))
    category[field][0] = new_value
    changed = tmp_path / "changed.cif"
    _write_category(changed, category)
    with pytest.raises(ValueError):
        compare_atom_semantics(atom_semantics(unfinalized_cif), atom_semantics(changed))


def test_semantic_comparison_accepts_metadata_only_changes(unfinalized_cif, tmp_path):
    output = tmp_path / "v2.cif"
    finalize_cif_metadata(unfinalized_cif, output)
    summary = compare_atom_semantics(atom_semantics(unfinalized_cif), atom_semantics(output))
    assert summary["coordinates_exactly_unchanged"] is True
    assert summary["max_coordinate_displacement_A"] == 0.


@pytest.mark.parametrize("case_id", CASE_IDS)
def test_actual_six_source_inputs_retain_all_atom_values(case_id, tmp_path):
    archive = ROOT / "results/structural_test/source_preparation.tar.gz"
    if archive.is_file():
        with tarfile.open(archive) as bundle:
            source_bytes = bundle.extractfile(f"{case_id}/receptor.cif").read()
    else:
        source_bytes = (ROOT / "data/processed/heldout_structural_preflight" / case_id / "receptor.cif").read_bytes()
    source = tmp_path / f"{case_id}_prior.cif"
    source.write_bytes(source_bytes)
    output = tmp_path / f"{case_id}_v2.cif"
    before = _atom_category(source)
    receipt = finalize_cif_metadata(source, output)
    assert source.read_bytes() == source_bytes
    _assert_raw_atoms_unchanged(before, _atom_category(output))
    _assert_entity_references_consistent(output)
    assert receipt["atom_count"] == len(before["id"])
    assert receipt["max_coordinate_displacement_A"] == 0.
    assert output.read_bytes().endswith(b"\n#\n")
