#!/usr/bin/env python3
"""Independently read saved initializer products; never embed or optimize."""
import json
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.ring_geometry import cremer_pople6
from aav_benchmark.ideal_ring_conformers import validate_conformer


def review():
    base = ROOT / "results/ideal_ring_conformer_feasibility/ccd_ideal_pyranose_v3"
    policy = json.loads((base / "policy.json").read_text())
    golden_path = Path(__file__).parent / "expected_parent_graphs.json"
    golden = json.loads(golden_path.read_text())["cases"]
    rows = []
    for case, expected in golden.items():
        directory = base / case
        report = json.loads((directory / "report.json").read_text())
        sdf = directory / "ligand.sdf"
        molecule = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
        assert molecule is not None and sha256(sdf) == report["selected_sdf_sha256"]
        expected_graph = expected["canonical_isomeric_smiles"]
        assert Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True) == expected_graph
        from3d = Chem.Mol(molecule)
        Chem.AssignStereochemistryFrom3D(from3d, confId=0, replaceExistingTags=True)
        assert Chem.MolToSmiles(Chem.RemoveHs(from3d), isomericSmiles=True) == expected_graph
        mapping = report["atom_map"]
        keys = {atom["canonical_atom_key"]: atom["rdkit_index"] for atom in mapping}
        assert all(molecule.GetAtomWithIdx(atom["rdkit_index"]).GetProp("canonical_atom_key") == atom["canonical_atom_key"] for atom in mapping)
        validation = validate_conformer(molecule, 0, mapping, report["ring_templates"], report["stereo_tests"], policy)
        assert validation["valid"]
        left, right = report["chemical_graph"]["declared_links"][0]
        a, b = keys[left], keys[right]
        n1 = next(n.GetIdx() for n in molecule.GetAtomWithIdx(a).GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != b)
        n2 = next(n.GetIdx() for n in molecule.GetAtomWithIdx(b).GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != a)
        quartet = [n1, a, b, n2]
        original = rdMolTransforms.GetDihedralDeg(molecule.GetConformer(), *quartet)
        freedom = []
        for change in [5., 10., 20., 30.]:
            temporary = Chem.Mol(molecule)
            rdMolTransforms.SetDihedralDeg(temporary.GetConformer(), *quartet, original + change)
            altered = validate_conformer(temporary, 0, mapping, report["ring_templates"], report["stereo_tests"], policy)
            delta = max(abs(cremer_pople6(molecule.GetConformer().GetPositions()[ring["indices"]])["Q_A"]
                            - cremer_pople6(temporary.GetConformer().GetPositions()[ring["indices"]])["Q_A"]) for ring in report["ring_templates"])
            assert delta < 1e-12
            freedom.append({"rotation_deg": change, "all_gates_valid": altered["valid"], "failures": altered["failures"], "max_ring_Q_change_A": delta})
        assert any(row["all_gates_valid"] for row in freedom)
        assert sha256(sdf) == report["selected_sdf_sha256"]
        rows.append({"case_id": case, "saved_sdf_sha256": sha256(sdf), "report_sha256": sha256(directory / "report.json"),
                     "expected_full_parent_isomeric_graph_matches": True, "independently_reassigned_3D_stereo_matches_expected_graph": True,
                     "serialized_canonical_keys_complete": True, "all_saved_product_gates_pass": True,
                     "ring_count": len(validation["rings"]), "stereocenter_count": len(validation["stereochemistry"]),
                     "max_named_pucker_angle_to_CCD_deg": max(row["comparison"]["pucker_vector_angle_to_CCD_ideal_deg"] for row in validation["rings"]),
                     "max_chair_pole_distance_deg": max(row["pucker"]["chair_pole_distance_deg"] for row in validation["rings"]),
                     "minimum_abs_chiral_volume_A3": min(abs(row["observed_signed_volume"]) for row in validation["stereochemistry"]),
                     "minimum_nonbonded_heavy_vdw_sum_ratio": validation["min_nonbonded_heavy_vdw_sum_ratio"],
                     "temporary_copy_glycosidic_rotation_controls": freedom, "source_sdf_unchanged": True})
    result = {"schema_version": 1, "status": "independently_reviewed_all_saved_products", "cases": rows,
              "policy_sha256": sha256(base / "policy.json"), "algorithm_snapshot_sha256": sha256(base / "algorithm.py"),
              "ring_geometry_snapshot_sha256": sha256(base / "ring_geometry.py"), "golden_expected_graphs_sha256": sha256(golden_path),
              "review_source_sha256": sha256(Path(__file__)), "no_embedding_optimization_docking_or_native_scoring_executed": True,
              "interpretation": "Saved-product chemistry and geometry checks; no free-energy, population, affinity or experimental accuracy claim"}
    write_json(Path(__file__).parent / "initializer_v3_all7_independent_review.json", result)
    return result


if __name__ == "__main__":
    result = review()
    print(f"Independently checked {len(result['cases'])} saved products; no embedding, optimization or native scoring")
