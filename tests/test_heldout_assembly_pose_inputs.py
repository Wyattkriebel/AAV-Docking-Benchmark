import json
from pathlib import Path
import tarfile

import pytest

from aav_benchmark.heldout_assembly_pose_inputs import validate_control_sources
from aav_benchmark.provenance import sha256

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("profile,manifest_relative", [
    ("strict_monomer", "results/structural_assembly_pose_v3/strict_monomer_sources_manifest.json"),
    ("repair_monomer", "results/structural_pose_v3/repair_monomer_sources_manifest.json"),
])
def test_portable_paired_control_source_profiles_and_tamper_rejection(profile, manifest_relative, tmp_path):
    path = ROOT / manifest_relative
    assert path.is_file(), "The fixed public assembly source archive is required"
    manifest = json.loads(path.read_text()); archive = path.parent / manifest["archive"]
    assert sha256(archive) == manifest["archive_sha256"]
    with tarfile.open(archive) as bundle:
        for member in bundle:
            assert member.isfile() and not Path(member.name).is_absolute() and ".." not in Path(member.name).parts
            output = tmp_path / member.name; output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(bundle.extractfile(member).read())
    result = validate_control_sources(tmp_path, profile)
    assert result["case_ids"] == ["3ZWE", "5HZA"]
    assert {row["selected_author_chain"] for row in result["cases"]} == {"A"}
    # Omitted chains and their native sites remain represented by the full truth.
    assert all((tmp_path / row["case_id"] / "full_assembly_receptor.cif").is_file() for row in result["cases"])
    target = tmp_path / "3ZWE/ligand.pdbqt"; target.write_text(target.read_text() + "\n")
    with pytest.raises(ValueError, match="Assembly source changed"):
        validate_control_sources(tmp_path, profile)
