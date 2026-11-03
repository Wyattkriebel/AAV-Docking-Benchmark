"""Source-only twelve-case audit controls, with no docking or outcome access."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import shutil
import unittest

import gemmi

from aav_benchmark.atomic_reference import AtomicReferenceError, DEFAULT_ROOT
from aav_benchmark.heldout_reference import audit_all, audit_candidate, derive_copy_entries, verify_heldout_references
from aav_benchmark.public_source_replay import (
    OMITTED_PUBLICATION_NAMES, PAPER_DIRECTORY, source_bond_distance_replay, verify_public_source_references,
)


class HeldoutReferenceAuditControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = DEFAULT_ROOT
        cls.inputs = json.loads((cls.root / "data/archive/heldout_reference/audit_inputs.json").read_text())
        cls.entries = {e["pdb_id"]: e for e in cls.inputs["candidates"]}
        cls.policy = json.loads((cls.root / cls.inputs["policy_path"]).read_text())
        cls.actual = audit_all(root=cls.root)
        cls.diagnostics = {c["case_id"]: c for c in cls.actual["cases"]}

    def test_actual_all_twelve_keep_every_native_atom_row_and_aglycone(self):
        self.assertEqual(len(self.actual["cases"]), 12)
        self.assertEqual(sum(c["native_heavy_atom_row_count"] for c in self.actual["cases"]), 388)
        self.assertEqual(sum(c["native_unique_heavy_atom_count"] for c in self.actual["cases"]), 377)
        self.assertTrue(self.actual["source_only"])
        self.assertFalse(self.actual["prediction_outcomes_consulted"])
        self.assertTrue(all(c["local_sphere_sample_count"] == 33 for c in self.actual["cases"]))
        self.assertTrue(all(len(c["maps"]) == 2 for c in self.actual["cases"]))
        self.assertTrue(all("released_coordinate_delta_A" not in a for c in self.actual["cases"] for a in c["native_atom_map"]))

    def _mutated_coordinate_entry(self, directory, case_id, operation):
        entry = deepcopy(self.entries[case_id])
        source = entry["sources"]["coordinates"]
        doc = gemmi.cif.read_file(str(self.root / source["archive_path"]))
        category = doc[0].get_mmcif_category("_atom_site.")
        index = next(i for i in range(len(category["id"])) if category["label_asym_id"][i] == entry["selected_label_asym_id"])
        if operation == "remove":
            category = {key: [value for i, value in enumerate(values) if i != index] for key, values in category.items()}
        elif operation == "move":
            category["Cartn_x"][index] = f"{float(category['Cartn_x'][index]) + .01:.3f}"
        doc[0].set_mmcif_category("_atom_site.", category)
        payload = doc.as_string().encode()
        path = Path(directory) / "mutated.cif.gz"
        path.write_bytes(gzip.compress(payload, mtime=0))
        source.update(archive_path=str(path), archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      archive_bytes=path.stat().st_size, sha256=hashlib.sha256(payload).hexdigest(), bytes=len(payload))
        return entry

    def test_lost_native_atom_fails_inventory_identity_even_after_source_checksum_update(self):
        with tempfile.TemporaryDirectory() as directory:
            entry = self._mutated_coordinate_entry(directory, "5V6F", "remove")
            with self.assertRaisesRegex(AtomicReferenceError, "native parent coordinate identity"):
                audit_candidate(entry, self.root, self.policy)

    def test_moved_native_atom_fails_inventory_identity_even_after_source_checksum_update(self):
        with tempfile.TemporaryDirectory() as directory:
            entry = self._mutated_coordinate_entry(directory, "5V6F", "move")
            with self.assertRaisesRegex(AtomicReferenceError, "native parent coordinate identity"):
                audit_candidate(entry, self.root, self.policy)

    def test_map_byte_tampering_fails_before_fourier_parse(self):
        entry = deepcopy(self.entries["5V6F"])
        source = entry["audit_sources"]["2fofc_coefficients"]
        payload = bytearray((self.root / source["archive"]).read_bytes())
        payload[len(payload) // 2] ^= 1
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tampered.cif.gz"
            path.write_bytes(payload)
            source["archive"] = str(path)
            with self.assertRaisesRegex(AtomicReferenceError, "checksum mismatch"):
                audit_candidate(entry, self.root, self.policy)

    def test_2wnb_complete_parent_preserves_covalent_aromatic_aglycone(self):
        actual = self.diagnostics["2WNB"]
        self.assertEqual(actual["native_heavy_atom_row_count"], 35)
        aglycone = next(r for r in actual["residues"] if r["external_aglycone"])
        self.assertEqual(aglycone["ccd_id"], "OPO")
        self.assertEqual(aglycone["heavy_atom_row_count"], 10)
        self.assertFalse(aglycone["source_gate_pass"])
        self.assertTrue(any(a["aromatic"] for a in aglycone["atom_map"]))
        attachment = next(b for b in actual["native_bond_graph"] if b.get("struct_conn_id") == "covale7")
        self.assertEqual(attachment["struct_conn_distance_A"], 1.478)
        self.assertEqual(attachment["order"], "SING")

    def test_2wab_alternate_rows_are_retained_and_not_selected_by_fit(self):
        actual = self.diagnostics["2WAB"]
        self.assertEqual(actual["native_heavy_atom_row_count"], 67)
        terminal = actual["residues"][-1]
        self.assertEqual(terminal["alternate_ids"], ["A", "B"])
        self.assertEqual(terminal["heavy_atom_row_count"], 22)
        self.assertEqual(terminal["unique_heavy_atom_count"], 11)
        self.assertIn("multiple_native_alternates_no_unique_atomic_target", terminal["source_gate_failures"])
        self.assertFalse(actual["residues"][3]["source_gate_pass"])

    def test_2wab_alternate_link_records_connect_the_named_alternate_once(self):
        graph = self.diagnostics["2WAB"]["native_bond_graph"]
        links = [b for b in graph if b.get("struct_conn_id") in ["covale4", "covale5"]]
        self.assertEqual(len(links), 2)
        self.assertEqual([(b["struct_conn_id"], b["atom2"].split(":")[-1]) for b in links],
                         [("covale4", "A"), ("covale5", "B")])
        self.assertEqual(len({(b["atom1"], b["atom2"]) for b in links}), 2)


    def test_unknown_validation_is_not_replaced_with_unit_atom_occupancies(self):
        actual = self.diagnostics["3A2S"]
        self.assertEqual(len(actual["maps"]), 2)
        self.assertTrue(all(a["occupancy"] == 1 for a in actual["native_atom_map"]))
        self.assertTrue(all("avgoccu" not in r["validation"] for r in actual["residues"]))
        self.assertTrue(all("avgoccu_min_gate" in r["source_gate_failures"] for r in actual["residues"]))
        self.assertTrue(all(not r["source_gate_pass"] for r in actual["residues"]))

    def test_same_full_chemical_copies_keep_separate_biological_assembly_membership(self):
        copies = derive_copy_entries(self.entries["2RDK"], self.root)
        self.assertEqual([(c["selected_label_asym_id"], [m["assembly_id"] for m in c["assembly_memberships"]])
                          for c in copies["same_full_chemical_entries"]], [("C", ["1"]), ("D", ["2"])])
        copies = derive_copy_entries(self.entries["2G7C"], self.root)
        selected = [c["selected_label_asym_id"] for c in copies["same_full_chemical_entries"]
                    if any(m["assembly_id"] == "1" for m in c["assembly_memberships"])]
        self.assertEqual(selected, ["C", "D"])

    def test_partial_ligands_remain_known_sites_outside_full_chemical_copy_group(self):
        result = derive_copy_entries(self.entries["3ZWE"], self.root)
        self.assertEqual([c["selected_label_asym_id"] for c in result["same_full_chemical_entries"]], ["E", "F"])
        fragments = {c["label_asym_id"]: c for c in result["other_carbohydrate_sites"]}
        self.assertEqual(len(fragments["D"]["atom_map"]), 22)
        self.assertEqual(set(fragments), {"D", "G", "H", "I"})
        self.assertTrue(all(not c["scored"] for c in fragments.values()))
        self.assertTrue(all(not c["same_full_chemical_graph_as_primary"] for c in fragments.values()))


class HeldoutFrozenLedgerControls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        source_dir = DEFAULT_ROOT / "data/archive/heldout_reference"
        destination = self.root / "data/archive/heldout_reference"
        destination.mkdir(parents=True)
        # Reuse immutable real source bytes. Mutated ledger/manifest are private
        # regular files, so controls never write through source symlinks.
        for path in source_dir.iterdir():
            (destination / path.name).symlink_to(path)
        for name in ["heldout_candidates", "atomic_reference"]:
            (self.root / "data/archive" / name).symlink_to(DEFAULT_ROOT / "data/archive" / name)
        (self.root / "benchmarks").mkdir()
        (self.root / "benchmarks/heldout_candidate_inventory.json").symlink_to(DEFAULT_ROOT / "benchmarks/heldout_candidate_inventory.json")
        self.ledger_path = self.root / "benchmarks/heldout_reference_ledger.json"
        shutil.copyfile(DEFAULT_ROOT / "benchmarks/heldout_reference_ledger.json", self.ledger_path)
        self.ledger = json.loads(self.ledger_path.read_text())
        self.manifest_path = self.root / self.ledger["source_record"]
        self.manifest = json.loads(self.manifest_path.read_text())
        # These private mutation fixtures exercise numerical-mask/assembly
        # controls through the original strict helper. Publication bytes are
        # irrelevant to those controls and absent in a public checkout. This
        # synthetic fixture manifest is never published as scientific evidence.
        publication_paths = {PAPER_DIRECTORY + name for name in OMITTED_PUBLICATION_NAMES}
        self.manifest["archived_files"] = [row for row in self.manifest["archived_files"]
                                           if row["archive"] not in publication_paths]
        self.manifest_path.unlink()
        self.manifest_path.write_text(json.dumps(self.manifest))

    def tearDown(self):
        self.temporary.cleanup()

    def _repin_fixture_ledger(self):
        self.ledger_path.write_text(json.dumps(self.ledger))
        self.manifest["ledger_sha256"] = hashlib.sha256(self.ledger_path.read_bytes()).hexdigest()
        self.manifest_path.write_text(json.dumps(self.manifest))

    def test_actual_frozen_all_twelve_primary_masks_and_assembly_copies_replay(self):
        public_report = verify_public_source_references(DEFAULT_ROOT)
        self.assertEqual(public_report["retained_archive_records_verified"], 163)
        self.assertEqual(public_report["metadata_only_publication_count"], 11)
        report = public_report["numerical_reference_replay"]
        self.assertEqual(report["status"], "verified")
        self.assertEqual(report["case_count"], 12)
        self.assertEqual(report["source_atom_row_count"], 388)
        self.assertEqual(report["primary_scored_atom_count"], 298)
        self.assertEqual(report["qualified_assembly_copy_count"], 16)
        self.assertFalse(report["prediction_outcomes_consulted"])
        self.assertLessEqual(report["maximum_absolute_density_difference_sigma"], 1e-5)

    def test_frozen_primary_mask_atom_loss_fails_after_ledger_digest_is_rebound(self):
        case = next(case for case in self.ledger["cases"] if case["case_id"] == "5V6F")
        case["scored_atom_map"].pop()
        case["scored_atom_count"] -= 1
        self._repin_fixture_ledger()
        with source_bond_distance_replay(), self.assertRaisesRegex(AtomicReferenceError, "primary scored map: array length mismatch"):
            verify_heldout_references(self.root)

    def test_misplaced_assembly_copy_fails_even_with_large_density_tolerance(self):
        case = next(case for case in self.ledger["cases"] if case["case_id"] == "1I3H")
        case["qualified_copies"][1]["scored_atom_map"][0]["xyz_A"][0] += .01
        self._repin_fixture_ledger()
        with source_bond_distance_replay(), self.assertRaisesRegex(AtomicReferenceError, "transformed qualified copies.*xyz_A"):
            verify_heldout_references(self.root, density_tolerance_sigma=100)

    def test_equally_valid_copy_from_another_assembly_cannot_be_substituted(self):
        case = next(case for case in self.ledger["cases"] if case["case_id"] == "2RDK")
        copy = case["qualified_copies"][0]
        copy["source_copy_id"] = "2RDK:D"
        copy["assembly_id"] = "2"
        copy["copy_id"] = "2RDK:D@assembly2:1"
        self._repin_fixture_ledger()
        with source_bond_distance_replay(), self.assertRaisesRegex(AtomicReferenceError, "transformed qualified copies.*copy_id"):
            verify_heldout_references(self.root)


if __name__ == "__main__":
    unittest.main()
