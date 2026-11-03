import io
import json
from pathlib import Path
import tarfile
import tempfile

import pytest

from aav_benchmark.heldout_serialization import atom_semantics, compare_atom_semantics
from aav_benchmark.provenance import sha256

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("case_id,observed_count,modeled_count", [("3ZWE", 661, 0), ("5HZA", 2414, 14)])
def test_paired_monomer_archive_is_exact_modeled_subset_with_original_truth(case_id, observed_count, modeled_count):
    manifest_path = ROOT / "results/structural_pose_v3/repair_monomer_sources_manifest.json"
    assert manifest_path.is_file(), "The fixed public monomer source archive is required"
    manifest = json.loads(manifest_path.read_text()); archive_path = manifest_path.parent / manifest["archive"]
    assert sha256(archive_path) == manifest["archive_sha256"]
    with tarfile.open(archive_path) as archive:
        members = {member.name: archive.extractfile(member).read() for member in archive if member.isfile()}
    prefix = case_id + "/"
    observed = json.loads(members[prefix + "protein_atom_map.json"])
    modeled = json.loads(members[prefix + "modeled_atom_map.json"])
    assert len(observed["atoms"]) == observed_count
    assert len(modeled["atoms"]) == modeled_count
    assert all(row["atomic_truth_eligible"] is False and row["source_atom_site_id"] is None for row in modeled["atoms"])
    assert not ({tuple(row["output_key"]) for row in observed["atoms"]} & {tuple(row["output_key"]) for row in modeled["atoms"]})
    with tempfile.TemporaryDirectory() as tmp:
        full, subset, experimental = (Path(tmp) / name for name in ("full.cif", "subset.cif", "experimental.cif"))
        full.write_bytes(members[prefix + "full_modeled_assembly_receptor.cif"])
        subset.write_bytes(members[prefix + "receptor.cif"])
        experimental.write_bytes(members[prefix + "full_assembly_receptor.cif"])
        expected = {key: row for key, row in atom_semantics(full).items() if key[1] == "A"}
        proof = compare_atom_semantics(expected, atom_semantics(subset))
        assert proof["atom_count"] == observed_count + modeled_count
        assert proof["coordinates_exactly_unchanged"] is True
        # Full experimental truth retains both/trimer chains, without modeled atoms.
        assert len({key[1] for key in atom_semantics(experimental)}) > 1
        original_map = json.loads(members[prefix + "full_protein_atom_map.json"])
        assert len(original_map["atoms"]) == len(atom_semantics(experimental))
