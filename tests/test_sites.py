"""Controls for reference-independent site inputs and fpocket evidence parsing."""
import gzip
import json
import os
from pathlib import Path
import platform
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import gemmi
import numpy as np

from aav_benchmark.provenance import sha256
from aav_benchmark.sites import (
    DEFAULT_LOCK, audit_fpocket_receptor, parse_fpocket, prepare_site_receptor, receptor_atoms,
    revise_site_receptor_serialization, run_fpocket,
)


def atom(serial, name="CA", chain="A", number=1, xyz=(0, 0, 0), element="C", *,
         record="ATOM", resname="ALA", icode="", alt=""):
    return (f"{record:6s}{serial:5d} {name:>4s}{alt:1s}{resname:3s} {chain:1s}{number:4d}{icode:1s}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n")


def pocket_fixture(directory, *, moved=False, sphere_count=2):
    """Minimal upstream-shaped fixture; never presented as scientific execution."""
    pockets = Path(directory) / "pockets"
    pockets.mkdir(parents=True)
    pqr = f"HEADER 0 - Pocket Score : 0.5000\nHEADER 2 - Number of V. Vertices : {sphere_count}\n"
    pqr += atom(1, "O", "", 1, (1, 2, 3), "O", resname="STP")
    pqr += atom(2, "O", "", 1, (3, 4, 5), "O", resname="STP")
    (pockets / "pocket1_vert.pqr").write_text(pqr)
    # Include fpocket's malformed unquoted descriptor prose before the atom
    # loop, and a real insertion code; parsing must preserve the raw artifact.
    cif = """data_pocket1_atm
loop_
_struct.pdbx_descriptor
This is invalid unquoted upstream descriptor prose.
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
"""
    cif += f"ATOM 1 C CA . ALA A 42 B {0.1 if moved else 0.0} 0.0 0.0 1 ? 42 A\n"
    (pockets / "pocket1_atm.cif").write_text(cif)
    (Path(directory) / "receptor_info.txt").write_text("Pocket 1 :\n Score : 0.500\n")


class SitePreparationControls(unittest.TestCase):
    def test_entire_protein_retained_while_all_nonprotein_removed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, output = directory / "source.pdb", directory / "receptor.cif"
            source.write_text(
                atom(1, xyz=(1, 2, 3)) + atom(2, chain="B", xyz=(150, 0, 0))
                + atom(3, name="H", element="H")
                + atom(4, name="C1", record="HETATM", resname="GAL", number=3)
                + atom(5, name="CA", element="Ca", record="HETATM", resname="CA", number=4)
                + atom(6, name="O", element="O", record="HETATM", resname="HOH", number=5))
            manifest = prepare_site_receptor(source, output)
            self.assertEqual(set(receptor_atoms(output)), {("A", 1, "", "ALA", "CA"), ("B", 1, "", "ALA", "CA")})
            self.assertEqual(manifest["protein_subunits"], 2)
            self.assertEqual(manifest["omitted_atoms"]["metal_atoms"], 1)
            self.assertEqual(manifest["omitted_atoms"]["water_atoms"], 1)
            self.assertEqual(manifest["omitted_atoms"]["nonprotein_atoms"], 2)
            self.assertEqual(manifest["omitted_atoms"]["protein_hydrogens"], 1)
            self.assertEqual(manifest["max_serialization_displacement_A"], 0)
            self.assertFalse(manifest["native_ligand_access"])
            self.assertFalse(manifest["distance_based_selection"])
            with gzip.open(manifest["atom_map"], "rt") as stream:
                self.assertIn("source_label_asym_id", stream.readline())
                self.assertEqual(len(stream.readlines()), 2)
            block = gemmi.cif.read(str(output)).sole_block()
            self.assertTrue(all(str(x).isdigit() for x in block.find_values("_entity.id")))
            self.assertTrue(output.read_text().endswith("\n#\n"))
            self.assertEqual(manifest["serialization_revision"], 2)

    def test_ligand_coordinates_cannot_select_or_change_receptor(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            hashes = []
            for index, ligand_position in enumerate(((0, 0, 0), (400, 500, 600))):
                source = directory / f"source{index}.pdb"
                source.write_text(atom(1, xyz=(1, 2, 3)) + atom(2, chain="B", xyz=(99, 0, 0))
                                  + atom(3, "C1", number=2, record="HETATM", resname="GAL", xyz=ligand_position))
                output = directory / f"receptor{index}.cif"
                hashes.append(prepare_site_receptor(source, output)["receptor_sha256"])
            self.assertEqual(hashes[0], hashes[1])

    def test_atom_encoded_calcium_is_removed_as_metal(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, output = directory / "source.pdb", directory / "receptor.cif"
            source.write_text(atom(1) + atom(2, "CA", number=2, element="Ca", resname="CA2"))
            manifest = prepare_site_receptor(source, output)
            self.assertEqual(manifest["protein_heavy_atoms"], 1)
            self.assertEqual(manifest["omitted_atoms"]["metal_atoms"], 1)
            self.assertEqual(manifest["omitted_residue_counts"], {"CA2": 1})

    def test_altloc_and_duplicate_identity_fail_without_selection(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(atom(1, alt="A"))
            with self.assertRaisesRegex(ValueError, "alternate"):
                prepare_site_receptor(source, directory / "receptor.cif")
            source.write_text(atom(1) + atom(2))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                prepare_site_receptor(source, directory / "receptor.cif")

    def test_cannot_overwrite_frozen_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, output = directory / "source.pdb", directory / "receptor.cif"
            source.write_text(atom(1))
            prepare_site_receptor(source, output)
            with self.assertRaises(FileExistsError):
                prepare_site_receptor(source, output)

    def test_pdb_chain_overflow_fails_and_mmcif_preserves_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            structure = gemmi.read_pdb_string(atom(1))
            structure[0][0].name = "LONG"
            source = directory / "source.cif"
            structure.make_mmcif_document().write_file(str(source))
            with self.assertRaisesRegex(ValueError, "overflow"):
                prepare_site_receptor(source, directory / "receptor.pdb")
            prepare_site_receptor(source, directory / "receptor.cif")
            self.assertIn(("LONG", 1, "", "ALA", "CA"), receptor_atoms(directory / "receptor.cif"))

    def test_serialization_revision_keeps_first_bytes_and_exact_atom_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, old = directory / "source.pdb", directory / "receptor_v1.cif"
            source.write_text(atom(1, number=42, icode="B", xyz=(1.234, 2.345, 3.456)))
            original_manifest = prepare_site_receptor(source, old)
            # Reproduce the pre-fix Gemmi serialization, without replacing any
            # actual scientific input: this is a synthetic revision fixture.
            old.write_text(old.read_text().removesuffix("#\n"))
            old_manifest_path = old.with_suffix(".cif.manifest.json")
            original_manifest["receptor_sha256"] = sha256(old)
            original_manifest.pop("serialization_revision")
            old_manifest_path.write_text(json.dumps(original_manifest))
            old_bytes, old_manifest_bytes = old.read_bytes(), old_manifest_path.read_bytes()
            new = directory / "receptor_v2.cif"
            manifest = revise_site_receptor_serialization(old, new)
            self.assertEqual(old.read_bytes(), old_bytes)
            self.assertEqual(old_manifest_path.read_bytes(), old_manifest_bytes)
            self.assertEqual(manifest["atom_map_sha256"], original_manifest["atom_map_sha256"])
            for key, xyz in receptor_atoms(old).items():
                np.testing.assert_array_equal(xyz, receptor_atoms(new)[key])
            receipt = json.loads(new.with_suffix(".cif.revision.json").read_text())
            self.assertEqual(receipt["prior_receptor_sha256"], sha256(old))
            self.assertEqual(receipt["receptor_sha256"], sha256(new))
            self.assertFalse(receipt["method_parameters_changed"])
            with self.assertRaises(FileExistsError):
                revise_site_receptor_serialization(old, new)

    @unittest.skipUnless(os.environ.get("AAV_FPOCKET_EXECUTABLE"), "set AAV_FPOCKET_EXECUTABLE for the actual pinned-tool smoke test")
    def test_real_fpocket_adapter_reads_100_atom_prepared_cif(self):
        """Actual executable integration; synthetic geometry, no site endpoint."""
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, receptor = directory / "synthetic100.pdb", directory / "receptor.cif"
            xyz = np.random.default_rng(20261030).uniform(0, 20, (100, 3))
            source.write_text("".join(atom(i, number=i, xyz=p) for i, p in enumerate(xyz, 1)))
            prepare_site_receptor(source, receptor)
            result = run_fpocket(receptor, directory / "execution", executable=os.environ["AAV_FPOCKET_EXECUTABLE"])
            self.assertEqual(result["status"], "completed", result.get("error"))
            self.assertEqual(result["receptor_output_audit"]["observed_protein_heavy_atoms"], 100)
            self.assertEqual(result["receptor_output_audit"]["max_displacement_A"], 0)


class FpocketOutputControls(unittest.TestCase):
    def setUp(self):
        # platform.platform() may invoke uname through subprocess on its first
        # call; resolve that independently of the runner's subprocess doubles.
        platform.platform()

    def test_alpha_sphere_barycenter_and_original_insertion_id(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            receptor = directory / "receptor.pdb"
            receptor.write_text(atom(1, number=42, icode="B"))
            pocket_fixture(directory / "raw")
            prediction, = parse_fpocket(directory / "raw", receptor)
            np.testing.assert_allclose(prediction["center_A"], [2, 3, 4])
            self.assertEqual(prediction["residue_keys"], [["A", 42, "B", "ALA"]])
            self.assertEqual(prediction["atom_ids"], [["A", 42, "B", "ALA", "CA"]])
            self.assertEqual(prediction["rank"], 1)

    def test_shifted_footprint_and_missing_spheres_are_failures(self):
        for moved, count in ((True, 2), (False, 3)):
            with tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                receptor = directory / "receptor.pdb"
                receptor.write_text(atom(1, number=42, icode="B"))
                pocket_fixture(directory / "raw", moved=moved, sphere_count=count)
                with self.assertRaises(ValueError):
                    parse_fpocket(directory / "raw", receptor)

    def test_full_receptor_audit_rejects_silent_missing_subunit(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            receptor = directory / "receptor.pdb"
            receptor.write_text(atom(1, number=42, icode="B") + atom(2, chain="Z", xyz=(100, 0, 0)))
            pocket_fixture(directory / "raw")
            with self.assertRaisesRegex(ValueError, "loss/addition"):
                audit_fpocket_receptor(directory / "raw/pockets/pocket1_atm.cif", receptor)

    def _fake_install(self, directory):
        executable = directory / "software/bin/fpocket"
        executable.parent.mkdir(parents=True)
        executable.write_text("test double; not a scientific method execution\n")
        lock = json.loads(DEFAULT_LOCK.read_text())["methods"]["fpocket"]
        installation = {k: lock[k] for k in ("release", "commit", "source_sha256")}
        installation["executable_sha256"] = sha256(executable)
        (executable.parent.parent / "installation.json").write_text(json.dumps(installation))
        return executable

    def test_runner_retains_nonzero_exit_and_logs_in_failure_ledger(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(atom(1))
            receptor = directory / "receptor.cif"
            prepare_site_receptor(source, receptor)
            executable = self._fake_install(directory)

            def failed(command, **kwargs):
                kwargs["stderr"].write("controlled subprocess failure\n")
                return SimpleNamespace(returncode=7)

            with patch("aav_benchmark.sites.subprocess.run", failed):
                result = run_fpocket(receptor, directory / "execution", executable=executable)
            self.assertEqual(result["status"], "failed")
            self.assertEqual(result["returncode"], 7)
            self.assertEqual(result["predictions"], [])
            self.assertIn("stderr.log", result["artifacts"])
            self.assertEqual(result["input_sha256"], sha256(receptor))
            self.assertIn("-m", result["command"])
            self.assertNotIn("-r", result["command"])
            self.assertNotIn("-P", result["command"])

    def test_zero_exit_without_outputs_is_not_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(atom(1))
            receptor = directory / "receptor.cif"
            prepare_site_receptor(source, receptor)
            executable = self._fake_install(directory)
            with patch("aav_benchmark.sites.subprocess.run", return_value=SimpleNamespace(returncode=0)):
                result = run_fpocket(receptor, directory / "execution", executable=executable)
            self.assertEqual(result["status"], "failed")
            self.assertIn("did not produce", result["error"])

    def test_explicit_no_pockets_is_a_completed_zero_candidate_result(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(atom(1))
            receptor = directory / "receptor.cif"
            prepare_site_receptor(source, receptor)
            executable = self._fake_install(directory)

            def no_pockets(command, **kwargs):
                kwargs["stdout"].write("no pockets found\n")
                return SimpleNamespace(returncode=0)

            with patch("aav_benchmark.sites.subprocess.run", no_pockets):
                result = run_fpocket(receptor, directory / "execution", executable=executable)
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["pocket_count"], 0)

    def test_mutated_input_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.pdb"
            source.write_text(atom(1))
            receptor = directory / "receptor.cif"
            prepare_site_receptor(source, receptor)
            receptor.write_text(receptor.read_text() + "# altered\n")
            executable = self._fake_install(directory)
            with self.assertRaisesRegex(ValueError, "differs"):
                run_fpocket(receptor, directory / "execution", executable=executable)
            self.assertFalse((directory / "execution").exists())


if __name__ == "__main__":
    unittest.main()
