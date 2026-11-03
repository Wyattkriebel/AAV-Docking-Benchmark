"""Actual experimental-source replay and controlled corruption; no docking run."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from aav_benchmark.atomic_reference import (
    AtomicReferenceError, DEFAULT_ROOT, replay_case, verify_atomic_references,
)

ROOT = DEFAULT_ROOT
LEDGER = Path("benchmarks/atomic_reference_ledger.json")
MANIFEST = Path("benchmarks/sources/atomic_reference_sources.json")


class AtomicReferenceReplayControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        for folder in ["data/archive/atomic_reference", "data/archive/external"]:
            shutil.copytree(ROOT / folder, self.root / folder)
        for path in [LEDGER, MANIFEST]:
            (self.root / path).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / path, self.root / path)
        self.ledger = json.loads((self.root / LEDGER).read_text())

    def tearDown(self):
        self.temporary.cleanup()

    def _write_ledger(self):
        (self.root / LEDGER).write_text(json.dumps(self.ledger))
        manifest = json.loads((self.root / MANIFEST).read_text())
        manifest["ledger_sha256"] = hashlib.sha256((self.root / LEDGER).read_bytes()).hexdigest()
        (self.root / MANIFEST).write_text(json.dumps(manifest))

    def _repin_mutated_ligand(self, case_id, payload):
        """Bypass only byte-integrity gates to exercise independent atom replay."""
        case = next(c for c in self.ledger["cases"] if c["case_id"] == case_id)
        source = case["released_ligand_source"]
        path = self.root / source["archive"]
        path.write_bytes(gzip.compress(payload, mtime=0))
        source["compressed_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        source["uncompressed_sha256"] = hashlib.sha256(payload).hexdigest()
        manifest = json.loads((self.root / MANIFEST).read_text())
        item = next(item for item in manifest["archived_files"] if item["archive"] == source["archive"])
        item["sha256"] = source["compressed_sha256"]
        item["bytes"] = path.stat().st_size
        (self.root / MANIFEST).write_text(json.dumps(manifest))
        self._write_ledger()

    def test_actual_all_seven_sources_replay_and_keep_denominator(self):
        report = verify_atomic_references(self.root)
        self.assertEqual(report["status"], "verified")
        self.assertEqual(len(report["cases"]), 7)
        self.assertEqual(sum(c["released_atoms"] for c in report["cases"]), 317)
        self.assertEqual(sum(c["scored_atoms"] for c in report["cases"]), 179)
        self.assertEqual(sum(c["residue_count"] for c in report["cases"]), 26)
        self.assertEqual([c["qualification"] for c in report["cases"]],
                         ["rejected", "full_ligand", "modeled_core", "modeled_core",
                          "modeled_core", "full_ligand", "modeled_core"])
        self.assertLessEqual(report["maximum_absolute_density_difference_sigma"], 1e-5)
        self.assertFalse(report["prediction_outcomes_consulted"])

    def test_lost_released_atom_fails_even_with_updated_source_checksum(self):
        path = self.root / "data/archive/external/1PWB_l_b.pdb.gz"
        lines = gzip.decompress(path.read_bytes()).splitlines(keepends=True)
        index = next(i for i, line in enumerate(lines) if line.startswith((b"HETATM", b"ATOM  ")))
        del lines[index]
        self._repin_mutated_ligand("1PWB", b"".join(lines))
        with self.assertRaisesRegex(AtomicReferenceError, "Missing released atom"):
            verify_atomic_references(self.root)

    def test_moved_released_atom_fails_even_with_updated_source_checksums(self):
        path = self.root / "data/archive/external/1PWB_l_b.pdb.gz"
        lines = gzip.decompress(path.read_bytes()).decode().splitlines(keepends=True)
        index = next(i for i, line in enumerate(lines) if line.startswith(("HETATM", "ATOM  ")))
        line = lines[index]
        lines[index] = line[:30] + f"{float(line[30:38]) + .01:8.3f}" + line[38:]
        self._repin_mutated_ligand("1PWB", "".join(lines).encode())
        with self.assertRaisesRegex(AtomicReferenceError, "source_gate_pass|source_gate_failures|released_coordinate_delta"):
            verify_atomic_references(self.root)

    def test_map_byte_tampering_fails_strict_checksum_before_fourier_parse(self):
        path = self.root / "data/archive/atomic_reference/1gny_validation_2fo-fc_map_coef.cif.gz"
        payload = bytearray(path.read_bytes())
        payload[len(payload) // 2] ^= 1
        path.write_bytes(payload)
        with self.assertRaisesRegex(AtomicReferenceError, "checksum mismatch:.*1gny_validation_2fo-fc"):
            verify_atomic_references(self.root)

    def test_scored_atom_loss_cannot_be_hidden_by_changing_count(self):
        case = next(c for c in self.ledger["cases"] if c["case_id"] == "1GNY")
        case["scored_atom_map"].pop()
        case["scored_atom_count"] -= 1
        self._write_ledger()
        with self.assertRaisesRegex(AtomicReferenceError, "scored atom map: array length mismatch"):
            verify_atomic_references(self.root)

    def test_density_tolerance_never_applies_to_native_atom_coordinates(self):
        case = next(c for c in self.ledger["cases"] if c["case_id"] == "1GNY")
        case["scored_atom_map"][0]["xyz_A"][0] += .001
        self._write_ledger()
        with self.assertRaisesRegex(AtomicReferenceError, "xyz_A"):
            verify_atomic_references(self.root, density_tolerance_sigma=100)

    def test_manual_adjudication_bytes_are_pinned_and_not_recomputed_as_truth(self):
        path = self.root / self.ledger["manual_adjudication_archive"]
        data = json.loads(path.read_text())
        data["decisions"]["2J1T"]["additional_residue_blocks"] = {}
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(AtomicReferenceError, "checksum mismatch:.*manual_adjudications"):
            verify_atomic_references(self.root)

    def test_actual_1qfo_chirality_defect_replays_from_named_neighbors(self):
        case = next(c for c in self.ledger["cases"] if c["case_id"] == "1QFO")
        actual = replay_case(case, self.root, self.ledger["required_policy"])
        glucose = actual["residues"][0]
        center = next(s for s in glucose["stereochemistry"] if s["atom"] == "C1")
        self.assertEqual(center["status"], "inverted")
        self.assertLess(center["ideal_signed_volume"] * center["observed_signed_volume"], 0)
        self.assertIn("wwPDB_chirality_outlier", glucose["source_gate_failures"])
        self.assertIn("CCD_signed_volume_stereochemistry_unresolved_or_inverted", glucose["source_gate_failures"])
        self.assertTrue(all(r["source_gate_pass"] for r in actual["residues"][1:]))


if __name__ == "__main__":
    unittest.main()
