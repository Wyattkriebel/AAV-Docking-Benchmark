"""Input quarantine, whole-chain coverage, immutable preparation and arm labels."""
import copy
import gzip
import json
from pathlib import Path
from unittest.mock import patch

import gemmi
import numpy as np
import pytest

from aav_benchmark.blind_docking import (
    _cached_protein, _protein_text, frozen_candidates, ligand_atoms, ligand_source,
    planning_inventory, prepare_box, prepare_condition, translate_ligand,
)
from aav_benchmark.engine import dock
from aav_benchmark.provenance import sha256, write_json


ROOT = Path(__file__).resolve().parents[1]


def _atom(name, point, atom_type="C", serial=1):
    return f"ATOM  {serial:5d} {name:>4s} LIG A   1    {point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00  0.00     0.000 {atom_type:>2s}"


def _protein(path):
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model("1"))
    for chain_name, points in [("A", [(0, 0, 0), (80, 0, 0)]), ("B", [(16, 0, 0)]), ("C", [(50, 0, 0)])]:
        chain = gemmi.Chain(chain_name)
        for number, xyz in enumerate(points, 1):
            residue = gemmi.Residue()
            residue.name = "ALA"
            residue.seqid = gemmi.SeqId(number, " ")
            residue.het_flag = "A"
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(*xyz)
            residue.add_atom(atom)
            chain.add_residue(residue)
        structure[0].add_chain(chain)
    structure.make_mmcif_document().write_file(str(path))


def test_all_supported_archived_graphs_are_independent_and_no_native_reference_is_opened():
    original_read = Path.read_bytes
    original_text = Path.read_text
    def guarded(path):
        if path.name.startswith("reference.json") or path.name.endswith("_l_b.pdb"):
            raise AssertionError("Native ligand coordinate reference accessed during blind preparation")
        return original_read(path)
    def guarded_text(path, *args, **kwargs):
        if path.name.startswith("reference.json") or path.name.endswith("_l_b.pdb"):
            raise AssertionError("Native ligand coordinate reference accessed during blind preparation")
        return original_text(path, *args, **kwargs)
    with patch.object(Path, "read_bytes", guarded), patch.object(Path, "read_text", guarded_text):
        sources = {case: ligand_source(ROOT, case) for case in ["aav9", "154L", "1GNY", "1KJL", "2J1T", "1PWB", "2VUZ"]}
    assert {case: value["box_size"][0] for case, value in sources.items()} == {
        "aav9": 24, "154L": 34, "1GNY": 41, "1KJL": 28, "2J1T": 32, "1PWB": 27, "2VUZ": 40}
    for source in sources.values():
        assert source["native_coordinates_used_for_embedding"] is False
        assert source["box_size"][0] / 2 - source["heavy_atom_centroid_radius_A"] >= 8
        assert {x["filename"] for x in source["source_receipts"]}.issubset({"preparation.json", "ligand.pdbqt", "galactose.pdbqt"})
    with pytest.raises(ValueError, match="No independently validated ligand"):
        ligand_source(ROOT, "1QFO")


def test_translation_removes_old_site_center_preserving_internal_geometry_and_names():
    text = "ROOT\n" + _atom("C1", [100, 20, -10]) + "\n" + _atom("O1", [102, 20, -10], "OA", 2) + "\nENDROOT\nTORSDOF 0\n"
    moved, audit = translate_ligand(text, [-5, 8, 12])
    xyz = np.asarray([a["coordinates"] for a in ligand_atoms(moved)])
    np.testing.assert_array_equal(xyz.mean(axis=0), [-5, 8, 12])
    np.testing.assert_array_equal(xyz[1] - xyz[0], [2, 0, 0])
    assert audit["old_translation_removed"] and audit["native_coordinate_access"] is False
    assert "ROOT" in moved and "TORSDOF 0" in moved
    with pytest.raises(ValueError, match="returned docking pose"):
        ligand_atoms("MODEL 1\n" + text)
    with pytest.raises(ValueError, match="three-vector"):
        translate_ligand(text, [0, float("nan"), 0])


def test_neighborhood_retains_complete_chains_and_audits_every_omitted_atom(tmp_path):
    path = tmp_path / "protein.cif"
    _protein(path)
    text, audit = _protein_text(path, [0, 0, 0], [20, 20, 20], cutoff=6)
    assert audit["selected_source_chains"] == ["A", "B"]
    assert audit["selected_heavy_atoms"] == 3
    assert audit["closest_omitted_atom_to_box_A"] == 40
    # Atom A:2 at x=80 remains, even though only its chain's first atom is near.
    assert any(float(line[30:38]) == 80 for line in text.splitlines() if line.startswith("ATOM"))


def test_predictor_failures_zero_candidates_and_native_center_quarantine(tmp_path):
    for condition, method, status in [("aav9_holo", "fpocket", "predictor_failed"), ("1QFO", "p2rank", "completed_no_candidates")]:
        result = frozen_candidates(ROOT, condition, method)
        assert result["status"] == status and result["candidates"] == []
    config = tmp_path / "config.json"
    config.write_text("{}")
    result = prepare_condition(ROOT, "154L", "p2rank", tmp_path / "forbidden", config, oracle_center=[1, 2, 3])
    assert result["status"] == "preparation_failed" and "forbidden" in result["reason"]
    with pytest.raises(FileExistsError):
        prepare_condition(ROOT, "154L", "p2rank", tmp_path / "forbidden", config)
    inventory = planning_inventory(ROOT)
    assert inventory["eligible_method_conditions"] == 27
    assert sum(r["available_box_count"] for r in inventory["rows"] if r["ligand_status"] == "ok") == 52


def test_receptor_cache_reuses_only_exact_inputs_and_rejects_changed_artifacts(tmp_path):
    calls = []
    def preparation(text, stem, definition, root):
        calls.append(text)
        stem.with_suffix(".pdb").write_text(text)
        stem.with_suffix(".pdbqt").write_text("typed " + text)
        write_json(stem.with_name(stem.name + "_atom_changes.json"), {"status": "ok"})
        return {"status": "ok"}
    cache = tmp_path / "cache"
    with patch("aav_benchmark.blind_docking._prepare_numbered_protein", preparation):
        for name in ("one", "two"):
            directory = tmp_path / name
            directory.mkdir()
            log, receipt = _cached_protein("identical input", directory / "receptor", None, None, cache)
            assert log["status"] == "ok" and receipt["key"]
        assert calls == ["identical input"]
        entry = next(cache.iterdir())
        (entry / "neighborhood.pdbqt").write_text("changed")
        (tmp_path / "three").mkdir()
        with pytest.raises(ValueError, match="cache bytes changed"):
            _cached_protein("identical input", tmp_path / "three/receptor", None, None, cache)


def test_box_runner_reference_never_contains_native_coordinate_truth(tmp_path):
    protein = tmp_path / "protein.cif"
    _protein(protein)
    config = tmp_path / "config.json"
    config.write_text("{}")
    ligand = {"pdbqt": _atom("C1", [100, 50, -30]) + "\n", "names": ["C1"], "box_size": [24, 24, 24]}
    def receptor(text, stem):
        stem.with_suffix(".pdb").write_text(text)
        stem.with_suffix(".pdbqt").write_text(text)
        return {"status": "ok"}
    with patch("aav_benchmark.blind_docking.prepare_receptor", receptor):
        prep = prepare_box(protein, ligand, tmp_path / "box", center=[0, 0, 0], condition_id="example",
                           complex_id="example", method="p2rank", candidate_class_rank=1, configuration_path=config)
    reference = json.loads((tmp_path / "box/reference.json").read_text())
    assert reference["pose_reference_eligible"] is False
    assert set(reference) == {"names", "pose_reference_eligible", "purpose"}
    assert prep["execution_information_regime"] == "predicted-site-local"
    assert prep["native_ligand_coordinate_access"] is False


@pytest.mark.parametrize("regime,expected", [(None, "true-site-local"), ("predicted-site-local", "predicted-site-local")])
def test_engine_preserves_legacy_label_and_explicitly_labels_predicted_input(tmp_path, monkeypatch, regime, expected):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    executable = tmp_path / "vina"
    executable.write_text("frozen executable fixture")
    (prepared / "protein.pdb").write_text("receptor fixture")
    (prepared / "protein.pdbqt").write_text("typed receptor fixture")
    (prepared / "ligand.pdbqt").write_text(_atom("C1", [0, 0, 0]))
    write_json(prepared / "reference.json", {"names": ["C1"], "pose_reference_eligible": False})
    prep = {"box_center": [0, 0, 0], "box_size": [24, 24, 24], "ligand_file": "ligand.pdbqt",
            "information_regime": "legacy descriptive label must not replace the execution contract",
            "receptors": [{"id": "box", "state": "holo", "assembly_condition": "complete", "pdb_file": "protein.pdb", "pdbqt_file": "protein.pdbqt"}]}
    if regime:
        prep["execution_information_regime"] = regime
    write_json(prepared / "preparation.json", prep)
    config = tmp_path / "config.json"
    write_json(config, {"seeds": [1], "cpu": 2, "exhaustiveness": 16, "num_modes": 20, "min_rmsd_A": .5, "energy_range_kcal_mol": 10})
    monkeypatch.setattr("aav_benchmark.engine.subprocess.check_output", lambda *a, **k: "AutoDock Vina v1.2.7")
    monkeypatch.setattr("aav_benchmark.engine.platform.platform", lambda: "fixture platform")
    def invoke(command, **kwargs):
        output = Path(command[command.index("--out") + 1])
        output.write_text("MODEL 1\nREMARK VINA RESULT: -4.2 0 0\n" + _atom("C1", [0, 0, 0]) + "\nENDMDL\n")
        class Completed:
            returncode = 0
        return Completed()
    monkeypatch.setattr("aav_benchmark.engine.subprocess.run", invoke)
    result = dock(prepared, tmp_path / "result", executable, config)
    assert result[0]["information_regime"] == expected
