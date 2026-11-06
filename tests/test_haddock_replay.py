"""Archive-only HADDOCK replay and corruption controls; never execute CNS."""
import csv
import json
from pathlib import Path

import pytest

from aav_benchmark import haddock_method as method
from aav_benchmark.haddock_replay import BASE, _compare, _raw_final_metrics, verify_haddock_replay


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def archive_root(tmp_path):
    root = tmp_path / "portable"
    manifest = json.loads((ROOT / BASE / "artifact_manifest.json").read_text())
    paths = {row["path"] for row in manifest["artifacts"]}
    paths.update({str(BASE / "artifact_manifest.json"), "benchmarks/haddock.lock.json",
                  "aav_benchmark/haddock_method.py", "data/archive/external/NAG.cif.gz",
                  "data/archive/external/NDG.cif.gz"})
    for name in paths:
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(ROOT / name)
    assert not (root / "tools").exists()
    return root


def replace_private(path, payload):
    path.unlink()
    path.write_bytes(payload)


def test_actual_archive_only_replay_without_any_software_installation(archive_root):
    result = verify_haddock_replay(archive_root)
    assert result["status"] == "verified"
    assert result["artifact_records_verified"] == 423
    assert result["retained_lock_paths_verified"] == 9
    assert result["external_release_pin_metadata_records"] == 5
    assert result["release_archives_and_CNS_installation_reverified"] is False
    assert result["execution_performed"] is False
    assert result["prepared_plus_model_chirality_records"] == 45
    assert result["rigidbody_models"] == 20
    assert result["refined_models"] == result["final_model_records"] == 12
    assert result["final_clusters"] == 6
    assert result["actual_glycan_atoms"] == 111
    assert result["actual_glycan_heavy_atoms"] == 57
    assert result["glycan_physical_bonds"] == 114
    assert result["known_interface_AIRs"] == 15
    assert result["generated_model_archives"] == 5
    assert result["released_example_self_evaluation"]["score_ranked_first"]["dockq"] == "0.215"
    assert result["released_example_self_evaluation"]["best_sampled_DockQ"]["caprieval_rank"] == "11"


def test_full_release_verifier_does_not_silently_become_archive_only(archive_root):
    with pytest.raises(method.HaddockAuditError, match="release/input checksum mismatch.*tools/"):
        method.verify_lock(archive_root)


@pytest.mark.parametrize("relative", [
    "official_glycan_test/run1-test/00_topoaa/1LMQ_l_u_haddock.psf.gz",
    "official_glycan_test/run1-test/13_caprieval/capri_ss.tsv",
])
def test_changed_topology_or_raw_score_table_fails_retained_pin(archive_root, relative):
    replace_private(archive_root / BASE / relative, b"synthetic changed evidence\n")
    with pytest.raises(method.HaddockAuditError, match="artifact checksum mismatch"):
        verify_haddock_replay(archive_root)


def test_rebound_manifest_cannot_discard_missing_model(archive_root):
    path = archive_root / BASE / "artifact_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["artifacts"].pop()
    manifest["artifact_count"] -= 1
    replace_private(path, json.dumps(manifest).encode())
    with pytest.raises(method.HaddockAuditError, match="artifact manifest changed"):
        verify_haddock_replay(archive_root)


def test_external_binary_metadata_cannot_be_rebound_in_release_lock(archive_root):
    path = archive_root / "benchmarks/haddock.lock.json"
    lock = json.loads(path.read_text())
    lock["cns_binary"]["sha256"] = "0" * 64
    replace_private(path, json.dumps(lock).encode())
    with pytest.raises(method.HaddockAuditError, match="release/input lock changed"):
        verify_haddock_replay(archive_root)


def test_ccd_handedness_sources_are_pinned_separately(archive_root):
    replace_private(archive_root / "data/archive/external/NAG.cif.gz", b"changed CCD")
    with pytest.raises(method.HaddockAuditError, match="CCD source checksum"):
        verify_haddock_replay(archive_root)


@pytest.mark.parametrize("change", ["duplicate_rank", "nonfinite_score"])
def test_raw_table_parser_rejects_invalid_rank_denominator_and_nonfinite_metrics(tmp_path, change):
    source = ROOT / BASE / "official_glycan_test/run1-test/13_caprieval/capri_ss.tsv"
    with source.open() as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        rows, fields = list(reader), reader.fieldnames
    if change == "duplicate_rank":
        rows[-1]["caprieval_rank"] = rows[0]["caprieval_rank"]
    else:
        rows[-1]["score"] = "nan"
    path = tmp_path / "synthetic_corrupt_capri.tsv"
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(method.HaddockAuditError, match="denominator|Nonfinite"):
        _raw_final_metrics(path)


def test_relocation_and_roundoff_do_not_relax_raw_coordinates_or_scores():
    old = {"path": "/old/location/results/haddock_reproduction/model.pdb",
           "observed_signed_volume": 2.1, "xyz_A": [1.0, 2.0, 3.0], "score": "-117.275"}
    new = {**old, "path": "/new/location/results/haddock_reproduction/model.pdb",
           "observed_signed_volume": 2.1 + 1e-12}
    _compare(new, old, "relocated informative geometry")
    with pytest.raises(method.HaddockAuditError, match="xyz_A"):
        _compare({**new, "xyz_A": [1.0 + 1e-12, 2.0, 3.0]}, old, "raw coordinates")
    with pytest.raises(method.HaddockAuditError, match="score"):
        _compare({**new, "score": "-117.274"}, old, "raw score")
