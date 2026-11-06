"""Carbohydrate score activation must survive PDBQT identity export."""
import gzip
import json
from pathlib import Path
import tarfile

import numpy as np
import pytest

from aav_benchmark.additional_method import (
    _atom, _extract, expected_chi_links, export_vina_carb_ligand,
    export_vina_carb_pose, read_carb_poses, validate_chi_probe,
    verify_chi_recognition, verify_released_example, verify_fixed_pose_validation,
)

ROOT = Path(__file__).resolve().parents[1]


def prepared_lactosamine(tmp_path):
    archive = ROOT / "results/development/1KJL/inputs"
    reference = json.loads(gzip.decompress((archive / "reference.json.gz").read_bytes()))
    source = tmp_path / "ligand.pdbqt"
    source.write_bytes(gzip.decompress((archive / "ligand.pdbqt.gz").read_bytes()))
    return source, reference["atom_map"]


def atoms(path):
    return [_atom(line) for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]


def test_export_restores_source_ring_and_residue_ids_without_moving_or_retyping_atoms(tmp_path):
    source, mapping = prepared_lactosamine(tmp_path)
    target = tmp_path / "vc.pdbqt"
    receipt = export_vina_carb_ligand(source, mapping, target)
    before, after = atoms(source), atoms(target)
    assert receipt["heavy_atoms"] == 26
    assert len(receipt["expected_glycosidic_links"]) == 1
    assert expected_chi_links(receipt) == [[("NAG", 1), ("GAL", 2)]]
    assert receipt["native_coordinates_accessed"] is False
    for original, exported in zip(before, after):
        assert {k: original[k] for k in ("serial", "coordinates", "charge", "type")} == {
            k: exported[k] for k in ("serial", "coordinates", "charge", "type")}
    assert [line for line in source.read_text().splitlines() if line.startswith(("BRANCH", "ENDBRANCH", "TORSDOF"))] == [
        line for line in target.read_text().splitlines() if line.startswith(("BRANCH", "ENDBRANCH", "TORSDOF"))]
    for line in target.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "C1":
            assert line[13:16].strip() == "C1"  # Actual upstream parser slice.


def test_export_refuses_atom_loss_unsupported_donor_and_changed_atom_serial_order(tmp_path):
    source, mapping = prepared_lactosamine(tmp_path)
    with pytest.raises(ValueError, match="map|universes"):
        export_vina_carb_ligand(source, mapping[:-1], tmp_path / "missing.pdbqt")
    altered = json.loads(json.dumps(mapping))
    for row in altered:
        if row["source_key"][3:] == ["GAL", "C1"]:
            row["source_key"][-1] = "C9"
    with pytest.raises(ValueError, match="donor"):
        export_vina_carb_ligand(source, altered, tmp_path / "unsupported.pdbqt")
    lines = source.read_text().splitlines()
    i = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM")))
    lines[i] = lines[i][:6] + "  999" + lines[i][11:]
    source.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="serial"):
        export_vina_carb_ligand(source, mapping, tmp_path / "serial.pdbqt")
    assert not (tmp_path / "serial.pdbqt").exists()


def test_pose_export_preserves_fixed_source_chemistry_but_uses_every_retained_coordinate(tmp_path):
    source, mapping = prepared_lactosamine(tmp_path)
    moved = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            original = _atom(line)
            xyz = np.asarray(original["coordinates"]) + [5., -2., 3.]
            line = line[:30] + "".join(f"{x:8.3f}" for x in xyz) + line[54:]
        moved.append(line)
    poses = tmp_path / "poses.pdbqt"
    poses.write_text("MODEL 1\nREMARK VINA RESULT: -7.0 0 0\n" + "\n".join(moved) + "\nENDMDL\n")
    destination = tmp_path / "pose_for_vc.pdbqt"
    receipt = export_vina_carb_pose(source, mapping, poses, 1, destination)
    assert receipt["model_rank"] == 1
    assert receipt["native_coordinates_accessed"] is False
    assert np.asarray([a["coordinates"] for a in atoms(destination)]) == pytest.approx(
        np.asarray([a["coordinates"] for a in atoms(source)]) + [5., -2., 3.])
    text = poses.read_text().replace("MODEL 1", "MODEL 2")
    poses.write_text(text)
    with pytest.raises(ValueError, match="contiguous"):
        export_vina_carb_pose(source, mapping, poses, 2, tmp_path / "bad.pdbqt")


def test_serial_pose_parser_restores_unique_names_and_rejects_identity_loss(tmp_path):
    source, mapping = prepared_lactosamine(tmp_path)
    exported = tmp_path / "carb.pdbqt"
    receipt = export_vina_carb_ligand(source, mapping, exported)
    pose = tmp_path / "pose.pdbqt"
    pose.write_text("MODEL 1\nREMARK VINA RESULT: -7.0 0 0\n" + exported.read_text() + "ENDMDL\n")
    result = read_carb_poses(pose, receipt)
    assert len(result) == 1 and len(result[0]["names"]) == 26
    assert len(set(result[0]["names"])) == 26
    lines = pose.read_text().splitlines()
    lines = [line for line in lines if not (line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) == 1)]
    pose.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="Incomplete"):
        read_carb_poses(pose, receipt)


def test_activation_is_not_inferred_from_exit_status_or_empty_log(tmp_path):
    log = tmp_path / "VC_log.txt"
    log.write_text("CHI energy penalties applied to NAG 1 and GAL 2 linkage.\n")
    assert verify_chi_recognition(log, [[("NAG", 1), ("GAL", 2)]])["recognized_links"] == 1
    for text in ("", "CHI energy penalties NOT applied to phi torsion in linkage b/w NAG 1 and GAL 2.\n"):
        log.write_text(text)
        with pytest.raises(ValueError, match="activate"):
            verify_chi_recognition(log, [[("NAG", 1), ("GAL", 2)]])


def test_fixed_coordinate_coefficient_control_rejects_inactive_or_non_linear_terms():
    good = {"initial_chi_kcal_mol": 0, "coefficient_0": 0, "coefficient_1": 10.5,
            "coefficient_2": 21., "glycosidic_links": 3, "fully_supported_links": 3, "torsion_perturbed": 3}
    assert validate_chi_probe(good)["status"] == "verified"
    for changed in ({"coefficient_1": 0}, {"coefficient_2": 22.}, {"fully_supported_links": 2},
                    {"coefficient_1": float("nan")}, {"coefficient_0": False}):
        with pytest.raises(ValueError):
            validate_chi_probe({**good, **changed})


def test_source_archive_cannot_escape_installation_directory(tmp_path):
    archive = tmp_path / "malicious.tar"
    import io
    with tarfile.open(archive, "w") as stream:
        item = tarfile.TarInfo("../outside");item.size = 1
        stream.addfile(item, io.BytesIO(b"x"))
    with pytest.raises(ValueError, match="Unsafe"):
        _extract(archive, tmp_path / "destination")
    assert not (tmp_path / "outside").exists()


def test_retained_tutorial_and_active_original_chi_objects_verify_without_native_tools():
    result = verify_released_example(ROOT)
    assert result["poses"] == 10 and result["heavy_atoms_per_pose"] == 46
    assert result["chi_recognition"]["recognized_links"] == 3
    assert result["actual_random_seed"] == 0


def test_original_scoring_objects_preserve_tutorial_rank_with_shared_reference():
    result = verify_fixed_pose_validation(ROOT)
    assert result["source_fullsearch_rank_order_preserved"] is True
    assert result["common_reference_offset_admits_all_printed_score_intervals"] is True
    assert result["offset_range_kcal_mol"][1] - result["offset_range_kcal_mol"][0] < .1


def test_export_never_overwrites_a_previous_method_input(tmp_path):
    source, mapping = prepared_lactosamine(tmp_path)
    destination = tmp_path / "prior.pdbqt"
    destination.write_bytes(b"prior evidence")
    with pytest.raises(FileExistsError):
        export_vina_carb_ligand(source, mapping, destination)
    assert destination.read_bytes() == b"prior evidence"
