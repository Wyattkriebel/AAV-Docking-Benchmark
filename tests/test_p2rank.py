"""Reference-free P2Rank parsing, identity, and failure-denominator controls."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import gemmi

from aav_benchmark.p2rank import parse_p2rank_predictions, run_p2rank
from aav_benchmark.provenance import sha256


PREDICTIONS_HEADER = "rank,score,probability,center_x,center_y,center_z,residue_ids,surf_atom_ids,sas_points\n"
RESIDUES_HEADER = "chain,residue_label,residue_name,pocket\n"


class P2RankControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.receptor = self.root / "receptor.cif"
        structure = gemmi.Structure()
        model = gemmi.Model("1")
        for chain_name, number, insertion, name in [("AA", -3, "", "ALA"), ("B", 12, "A", "GLY")]:
            chain = gemmi.Chain(chain_name)
            residue = gemmi.Residue()
            residue.name, residue.seqid, residue.het_flag = name, gemmi.SeqId(number, insertion or " "), "A"
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = "CA", gemmi.Element("C"), gemmi.Position(1, 2, 3)
            residue.add_atom(atom)
            chain.add_residue(residue)
            model.add_chain(chain)
        structure.add_model(model)
        structure.make_mmcif_document().write_file(str(self.receptor))
        self.predictions = self.root / "predictions.csv"
        self.residues = self.root / "residues.csv"
        self.predictions.write_text(PREDICTIONS_HEADER + "1,2.5,0.6,1,2,3,AA_-3 B_12A,1 2,4\n")
        self.residues.write_text(RESIDUES_HEADER + "AA,-3,ALA,1\nB,12A,GLY,1\n")

    def parse(self):
        return parse_p2rank_predictions(self.predictions, self.residues, self.receptor)

    def test_signed_author_numbers_long_chains_and_insertion_codes_preserved(self):
        source_hash = sha256(self.receptor)
        rows = self.parse()
        self.assertEqual(rows[0]["residues"], [
            {"chain": "AA", "number": -3, "insertion_code": "", "resname": "ALA"},
            {"chain": "B", "number": 12, "insertion_code": "A", "resname": "GLY"}])
        self.assertEqual(rows[0]["center"], [1.0, 2.0, 3.0])
        self.assertEqual(sha256(self.receptor), source_hash)

    def test_zero_predictions_is_a_valid_executed_result(self):
        self.predictions.write_text(PREDICTIONS_HEADER)
        self.residues.write_text(RESIDUES_HEADER + "AA,-3,ALA,0\nB,12A,GLY,0\n")
        self.assertEqual(self.parse(), [])

    def test_footprints_may_overlap_between_pockets(self):
        self.predictions.write_text(PREDICTIONS_HEADER + "1,2.5,0.6,1,2,3,AA_-3 B_12A,1 2,4\n"
                                    "2,1.5,0.2,4,5,6,B_12A,2,2\n")
        self.assertEqual(len(self.parse()), 2)

    def test_unmapped_duplicated_and_incomplete_identity_outputs_fail(self):
        valid_predictions = self.predictions.read_text()
        valid_residues = self.residues.read_text()
        variants = [
            (valid_predictions.replace("B_12A", "B_12"), valid_residues),
            (valid_predictions.replace("B_12A", "AA_-3"), valid_residues),
            (valid_predictions, RESIDUES_HEADER + "AA,-3,ALA,1\n"),
            (valid_predictions, valid_residues + "B,12A,GLY,1\n"),
            (valid_predictions, valid_residues.replace("GLY", "ALA")),
        ]
        for predictions, residues in variants:
            with self.subTest(predictions=predictions, residues=residues):
                self.predictions.write_text(predictions)
                self.residues.write_text(residues)
                with self.assertRaises(ValueError):
                    self.parse()

    def test_invalid_ranks_probability_or_nonfinite_coordinates_fail(self):
        for row in ["2,2.5,0.6,1,2,3,AA_-3,1,4", "1,2.5,1.6,1,2,3,AA_-3,1,4",
                    "1,2.5,0.6,nan,2,3,AA_-3,1,4", "1,inf,0.6,1,2,3,AA_-3,1,4"]:
            with self.subTest(row=row):
                self.predictions.write_text(PREDICTIONS_HEADER + row + "\n")
                with self.assertRaises(ValueError):
                    self.parse()

    def test_missing_executable_is_an_archived_failure(self):
        output = self.root / "failed"
        result = run_p2rank(self.receptor, output, self.root / "absent")
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["failure_category"], "installation_validation")
        self.assertEqual(result["predictions"], [])
        self.assertEqual(result["input_sha256"], sha256(self.receptor))
        self.assertGreater(result["runtime_seconds"], 0)
        self.assertEqual(json.loads((output / "result.json").read_text()), result)
        self.assertTrue(any(row["path"] == "adapter_error.log" for row in result["artifacts"]))

    def test_existing_output_is_rejected_before_stale_csvs_can_be_read(self):
        output = self.root / "stale"
        output.mkdir()
        (output / "stdout.log").write_text("stale")
        with self.assertRaisesRegex(ValueError, "empty"):
            run_p2rank(self.receptor, output, self.root / "absent")

    def test_timeout_retains_denominator_and_terminates_process_group(self):
        executable = self.root / "prank"
        executable.write_text("#!/bin/sh\nsleep 30 &\nwait\n")
        executable.chmod(0o755)
        lock = {"parameters": {"model": "default"}}
        output = self.root / "timeout"
        environment = {**os.environ, "JAVA_HOME": str(self.root)}
        with patch("aav_benchmark.p2rank._installation", return_value=(lock, environment, "fixture")):
            result = run_p2rank(self.receptor, output, executable, timeout_seconds=0.05)
        self.assertEqual(result["status"], "timeout")
        self.assertEqual(result["failure_category"], "timeout")
        self.assertEqual(result["predictions"], [])
        self.assertLess(result["runtime_seconds"], 10)
        self.assertNotEqual(result["returncode"], 0)
        self.assertTrue((output / "stdout.log").is_file())
        self.assertTrue((output / "stderr.log").is_file())


if __name__ == "__main__":
    unittest.main()
