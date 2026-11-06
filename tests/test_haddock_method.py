"""Read-only tests of the actual released HADDOCK example and corruption controls.

No test executes CNS, changes scientific inputs, or reads benchmark docking runs.
"""
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import gemmi
import numpy as np

from aav_benchmark.haddock_method import (
    HaddockAuditError, ROOT, audit_air, glycan_chirality, pdb_atoms,
    pdb_inventory, preparation_correspondence, read_cns_topology,
    verify_lock, workflow_inventory,
)

EXAMPLE = ROOT / "results/haddock_reproduction/official_glycan_test"
RUN = EXAMPLE / "run1-test"
INPUT = EXAMPLE / "data/1LMQ_l_u.pdb"
PREPARED = RUN / "00_topoaa/1LMQ_l_u_haddock.pdb.gz"
TOPOLOGY = RUN / "00_topoaa/1LMQ_l_u_haddock.psf.gz"
CCDS = {name: ROOT / f"data/archive/external/{name}.cif.gz" for name in ["NDG", "NAG"]}
LOCK = Path("benchmarks/haddock.lock.json")


def _text(path):
    return gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()


def _replace_atom_coordinates(lines, *, residue, atom_name, xyz):
    replaced = 0
    result = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and line[22:26].strip() == residue and line[12:16].strip() == atom_name:
            line = line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]
            replaced += 1
        result.append(line)
    if replaced != 1:
        raise AssertionError("Control must change exactly one existing atom")
    return result


class ActualHaddockChemicalAudit(unittest.TestCase):
    def test_actual_input_and_preparation_expose_rebuilding(self):
        before, after = pdb_inventory(INPUT), pdb_inventory(PREPARED)
        self.assertEqual((before["residue_count"], before["heavy_atom_count"]), (4, 56))
        self.assertEqual((after["residue_count"], after["heavy_atom_count"]), (4, 57))
        observed = preparation_correspondence(INPUT, PREPARED)
        self.assertEqual(observed["exact_name_shared_heavy_atoms"], 44)
        self.assertEqual(observed["shared_max_displacement_A"], 0.0)
        self.assertEqual(len(observed["input_only_atoms"]), 12)
        self.assertEqual(len(observed["prepared_only_atoms"]), 13)
        self.assertIn(["1", "NDG", "O1"], observed["prepared_only_atoms"])
        shift = next(row for row in observed["semantic_alias_coordinate_comparisons"]
                     if row["input_atom"] == ["3", "NAG", "O2N"])
        self.assertEqual(shift["standard_semantic_atom"], ["3", "NAG", "O7"])
        self.assertTrue(shift["output_present"])
        # A rename preserves coordinates; this measured displacement disproves it.
        self.assertGreater(shift["coordinate_displacement_A"], 4.5)
        self.assertLess(shift["coordinate_displacement_A"], 4.6)

    def test_actual_topology_matches_pdb_has_four_rings_and_neutral_charge(self):
        topology = read_cns_topology(TOPOLOGY)
        atoms = pdb_atoms(PREPARED)
        pdb_keys = {(a["chain"], a["resnum"], a["resname"], a["atom_name"]) for a in atoms}
        topology_keys = {(a["segment"], a["resnum"], a["resname"], a["atom_name"])
                         for a in topology["atoms"].values()}
        self.assertEqual(topology_keys, pdb_keys)
        self.assertEqual(topology["atom_count"], len(atoms))
        self.assertEqual(topology["atom_count"], 111)
        self.assertEqual(topology["sentinel_atom_rows_ignored"], 1)
        self.assertEqual(topology["sentinel_bond_rows_ignored"], 1)
        self.assertAlmostEqual(topology["total_charge"], 0.0, places=12)
        adjacency = {key: set() for key in topology["atoms"]}
        for left, right in topology["bonds"]:
            adjacency[left].add(right)
            adjacency[right].add(left)
        seen, pending = set(), [next(iter(adjacency))]
        while pending:
            atom = pending.pop()
            if atom not in seen:
                seen.add(atom)
                pending.extend(adjacency[atom] - seen)
        self.assertEqual(seen, set(adjacency))
        self.assertEqual(topology["bond_count"] - topology["atom_count"] + 1, 4)

    def test_actual_ccd_handedness_uses_three_real_glycosidic_links(self):
        report = glycan_chirality(PREPARED, TOPOLOGY, CCDS)
        self.assertEqual(report["status"], "consistent")
        self.assertEqual(report["tested_centers"], 20)
        self.assertTrue(all(c["ideal_signed_volume"] * c["observed_signed_volume"] > 0
                            for c in report["centers"]))
        substitutions = {tuple(c["residue"]): c["substitutions"]["O1"]
                         for c in report["centers"] if c["atom"] == "C1" and c["substitutions"]}
        self.assertEqual(substitutions, {
            ("B", "2", "NAG"): ["B", "1", "NDG", "O4"],
            ("B", "3", "NAG"): ["B", "2", "NAG", "O4"],
            ("B", "4", "NAG"): ["B", "3", "NAG", "O4"],
        })
        self.assertEqual(len(report["actual_interresidue_CNS_bonds"]), 3)

    def test_local_neighbor_reflection_detects_wrong_handed_c1(self):
        by_name = {a["atom_name"]: np.array(a["xyz_A"]) for a in pdb_atoms(PREPARED)
                   if a["resnum"] == "1"}
        center = by_name["C1"]
        normal = np.cross(by_name["C2"] - center, by_name["O5"] - center)
        normal /= np.linalg.norm(normal)
        vector = by_name["O1"] - center
        reflected = center + vector - 2 * np.dot(vector, normal) * normal
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wrong_hand.pdb"
            changed = _replace_atom_coordinates(_text(PREPARED).splitlines(keepends=True),
                                                residue="1", atom_name="O1", xyz=reflected)
            path.write_text("".join(changed))
            report = glycan_chirality(path, TOPOLOGY, CCDS)
        tested = next(c for c in report["centers"] if c["residue"] == ["B", "1", "NDG"] and c["atom"] == "C1")
        self.assertEqual(tested["status"], "inverted")
        self.assertLess(tested["ideal_signed_volume"] * tested["observed_signed_volume"], 0)
        self.assertEqual(report["status"], "failed")

    def test_missing_prepared_chiral_center_fails(self):
        lines = [line for line in _text(PREPARED).splitlines(keepends=True)
                 if not (line.startswith(("ATOM  ", "HETATM")) and line[22:26].strip() == "1"
                         and line[12:16].strip() == "C1")]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lost_center.pdb"
            path.write_text("".join(lines))
            with self.assertRaisesRegex(HaddockAuditError, "inventory differs|Missing chiral center"):
                glycan_chirality(path, TOPOLOGY, CCDS)

    def test_missing_rebuilt_nonchiral_atom_fails_topology_inventory(self):
        # O7 is not itself a chiral center or a named neighbor in the C1-C5
        # handedness checks. A chirality-only audit would silently miss its loss.
        lines = [line for line in _text(PREPARED).splitlines(keepends=True)
                 if not (line.startswith(("ATOM  ", "HETATM")) and line[22:26].strip() == "3"
                         and line[12:16].strip() == "O7")]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lost_nonchiral.pdb"
            path.write_text("".join(lines))
            with self.assertRaisesRegex(HaddockAuditError, "inventory differs"):
                glycan_chirality(path, TOPOLOGY, CCDS)

    def test_duplicate_pdb_atom_identity_fails(self):
        text = _text(PREPARED)
        first = next(line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.pdb"
            path.write_text(text + first + "\n")
            with self.assertRaisesRegex(HaddockAuditError, "Duplicate PDB atom"):
                pdb_inventory(path)

    def test_unknown_positive_topology_bond_endpoint_fails(self):
        document = gemmi.cif.read_string(_text(TOPOLOGY))
        document[0].find_values("_cns_mtf_bond.id[2]")[0] = "999999"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dangling.psf"
            path.write_text(document.as_string())
            with self.assertRaisesRegex(HaddockAuditError, "absent atom"):
                read_cns_topology(path)

    def test_air_inventory_retains_known_site_information_condition(self):
        report = audit_air(EXAMPLE / "data/ambig.tbl")
        self.assertEqual(report["assignment_count"], 15)
        self.assertEqual(report["glycan_residues"], [1, 2, 3, 4])
        self.assertEqual(report["assigned_protein_residues"],
                         [35, 46, 52, 57, 58, 59, 62, 63, 75, 101, 102, 103, 107, 108, 109])
        self.assertIn("not blind", report["information_condition"])

    def test_actual_workflow_model_counts_follow_io_without_new_cns_run(self):
        stages = {row["stage"]: row for row in workflow_inventory(RUN)["stages"]}
        lock = json.loads((ROOT / LOCK).read_text())
        initial = json.loads((RUN / "01_rigidbody/io.json").read_text())
        refinement = json.loads((RUN / "07_flexref/io.json").read_text())
        selected = json.loads((RUN / "11_seletopclusts/io.json").read_text())
        self.assertEqual(len(stages), 14)
        self.assertTrue(all(row["io_present"] and row["params_present"] for row in stages.values()))
        self.assertEqual(stages["01_rigidbody"]["physical_pdb_files"], lock["official_example"]["sampling"])
        self.assertEqual(stages["01_rigidbody"]["physical_pdb_files"], len(initial["output"]))
        self.assertEqual(stages["07_flexref"]["physical_pdb_files"], len(refinement["input"]))
        self.assertEqual(len(refinement["output"]), len(refinement["input"]))
        self.assertEqual(len(selected["output"]), len(refinement["output"]))
        for stage in ["01_rigidbody", "07_flexref"]:
            for model in stages[stage]["model_inventory"]:
                self.assertEqual(model["chain_heavy_atoms"]["B"], 57)


@unittest.skipUnless(
    (ROOT / "tools/haddock3-2026.10.1/archives/haddock3-2026.10.1.tar.gz").is_file()
    and (ROOT / "tools/haddock3-2026.10.1/env/lib/python3.11/site-packages/haddock/cns/bin/arm64-darwin.bin").is_file(),
    "requires the pinned local HADDOCK archives and installation; chemical artifact tests do not execute CNS",
)
class HaddockReleasePinControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.lock = json.loads((ROOT / LOCK).read_text())
        # Actual release artifacts are read through symlinks. Any perturbation
        # replaces the local link before writing, never modifying source files.
        def visit(value):
            if isinstance(value, dict):
                if {"path", "sha256", "bytes"} <= set(value):
                    target = self.root / value["path"]
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if not target.exists():
                        target.symlink_to(ROOT / value["path"])
                for item in value.values():
                    visit(item)
            elif isinstance(value, list):
                for item in value:
                    visit(item)
        visit(self.lock)
        self._write_lock()

    def tearDown(self):
        self.temporary.cleanup()

    def _write_lock(self):
        path = self.root / LOCK
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.lock))

    def test_actual_pinned_release_verifies(self):
        report = verify_lock(self.root)
        self.assertEqual(report["status"], "verified")
        self.assertGreaterEqual(report["checked_paths"], 14)
        self.assertEqual(report["cns_sha256"], self.lock["cns_binary"]["sha256"])

    def test_lost_official_input_atom_fails_source_checksum(self):
        source = self.lock["official_example"]["input_files"][0]
        path = self.root / source["path"]
        original = path.read_text().splitlines(keepends=True)
        index = next(i for i, line in enumerate(original) if line.startswith(("ATOM  ", "HETATM")))
        del original[index]
        path.unlink()
        path.write_text("".join(original))
        with self.assertRaisesRegex(HaddockAuditError, "checksum mismatch"):
            verify_lock(self.root)

    def test_deleting_binary_pin_field_cannot_silently_skip_integrity_check(self):
        del self.lock["cns_binary"]["bytes"]
        self._write_lock()
        with self.assertRaises(HaddockAuditError):
            verify_lock(self.root)

    def test_known_site_example_cannot_be_relabeled_blind(self):
        self.lock["official_example"]["information_condition"] = "blind"
        self._write_lock()
        with self.assertRaisesRegex(HaddockAuditError, "condition changed|relabeled blind"):
            verify_lock(self.root)


if __name__ == "__main__":
    unittest.main()
