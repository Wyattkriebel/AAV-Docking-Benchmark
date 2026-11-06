"""Archive controls and parser replay without installed scientific programs."""
import gzip
import io
import json
import os
from pathlib import Path
import shutil
import tarfile

import numpy as np
import pytest

from aav_benchmark.archive import bundle, restore
from aav_benchmark.p2rank import parse_p2rank_predictions
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.site_archive import pack_run, restore_run
from aav_benchmark.sites import parse_fpocket


def _record(run, *, method="fpocket", status=None, record_name=None):
    record_name = record_name or ("execution.json" if method == "fpocket" else "result.json")
    status = status or ("completed" if method == "fpocket" else "ok")
    files = [path for path in sorted(run.rglob("*")) if path.is_file()]
    if method == "fpocket":
        artifacts = {path.relative_to(run).as_posix(): {"sha256": sha256(path), "bytes": path.stat().st_size}
                     for path in files}
    else:
        artifacts = [{"path": path.relative_to(run).as_posix(), "sha256": sha256(path), "size_bytes": path.stat().st_size}
                     for path in files]
    record = {"schema_version": 1, "method": method, "status": status, "input_sha256": "1" * 64,
              "command": ["/never/execute/this/tool"], "input_path": "/never/read/this/input.cif",
              "artifacts": artifacts}
    write_json(run / record_name, record)
    return record_name


def _minimal_run(tmp_path, *, method="fpocket", status=None):
    run = tmp_path / "run"
    run.mkdir()
    (run / "stdout.log").write_bytes(b"raw stdout\n")
    (run / "stderr.log").write_bytes(b"")
    (run / "raw").mkdir()
    (run / "raw/output.bin").write_bytes(bytes(range(256)))
    return run, _record(run, method=method, status=status)


def _repackage(run, transform):
    """Rewrite an adversarial tar and update only the container checksum."""
    archive_path = run / "raw_artifacts.tar.gz"
    with tarfile.open(archive_path, "r:gz") as archive:
        members = [(member, archive.extractfile(member).read()) for member in archive]
    members = transform(members)
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                for member, data in members:
                    archive.addfile(member, io.BytesIO(data) if member.isfile() else None)
    manifest_path = run / "raw_archive.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update(archive_sha256=sha256(archive_path), archive_bytes=archive_path.stat().st_size)
    write_json(manifest_path, manifest)


@pytest.mark.parametrize("method", ["fpocket", "p2rank"])
def test_roundtrip_is_exact_portable_deterministic_and_preserves_originals(tmp_path, method):
    run, record_name = _minimal_run(tmp_path, method=method)
    originals = {path.relative_to(run).as_posix(): path.read_bytes() for path in run.rglob("*") if path.is_file()}
    manifest = pack_run(run, record_name)
    assert manifest["record_sha256"] == sha256(run / record_name)
    assert all(not Path(entry["path"]).is_absolute() for entry in manifest["entries"])
    assert pack_run(run, record_name) == manifest
    for name, data in originals.items():
        assert (run / name).read_bytes() == data
    independent = tmp_path / "independent"
    shutil.copytree(run, independent, ignore=shutil.ignore_patterns("raw_archive.json", "raw_artifacts.tar.gz"))
    for path in independent.rglob("*"):
        if path.is_file():
            os.utime(path, (123456, 123456))
            path.chmod(0o600)
    other_manifest = pack_run(independent, record_name)
    assert other_manifest == manifest
    assert (independent / "raw_artifacts.tar.gz").read_bytes() == (run / "raw_artifacts.tar.gz").read_bytes()
    # Reproduction needs only the committed archive and manifest, not the old
    # workstation's source files, external input/tool paths, or installed tools.
    portable = tmp_path / "portable"
    portable.mkdir()
    for name in ("raw_archive.json", "raw_artifacts.tar.gz"):
        shutil.copy2(run / name, portable / name)
    destination = tmp_path / "restored"
    assert restore_run(portable, destination) == manifest
    assert {path.relative_to(destination).as_posix(): path.read_bytes()
            for path in destination.rglob("*") if path.is_file()} == originals


def test_running_attempt_and_changed_or_missing_artifacts_cannot_be_packed(tmp_path):
    run, record_name = _minimal_run(tmp_path, status="running")
    with pytest.raises(ValueError, match="terminal"):
        pack_run(run, record_name)
    record = json.loads((run / record_name).read_text())
    record["status"] = "failed"
    write_json(run / record_name, record)
    (run / "stdout.log").write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum or size"):
        pack_run(run, record_name)
    (run / "stdout.log").unlink()
    with pytest.raises(FileNotFoundError):
        pack_run(run, record_name)
    assert not (run / "raw_artifacts.tar.gz").exists()


@pytest.mark.parametrize("path", ["../escape", "/absolute", "raw/../escape", "raw//file", "./file", "C:/file", "raw\\file"])
def test_record_artifact_path_traversal_is_rejected(tmp_path, path):
    run, record_name = _minimal_run(tmp_path)
    record = json.loads((run / record_name).read_text())
    record["artifacts"][path] = {"sha256": "0" * 64, "bytes": 0}
    write_json(run / record_name, record)
    with pytest.raises(ValueError, match="path"):
        pack_run(run, record_name)


@pytest.mark.parametrize("parent_link", [False, True])
def test_pack_rejects_symlink_file_or_directory_components(tmp_path, parent_link):
    run, record_name = _minimal_run(tmp_path)
    if parent_link:
        original = run / "raw"
        outside = tmp_path / "outside"
        original.rename(outside)
        original.symlink_to(outside, target_is_directory=True)
    else:
        (run / "stdout.log").unlink()
        outside = tmp_path / "outside.log"
        outside.write_bytes(b"raw stdout\n")
        (run / "stdout.log").symlink_to(outside)
    with pytest.raises(ValueError, match="symlink"):
        pack_run(run, record_name)


def test_restore_rejects_existing_destination_and_corrupted_container(tmp_path):
    run, record_name = _minimal_run(tmp_path)
    pack_run(run, record_name)
    destination = tmp_path / "existing"
    destination.mkdir()
    sentinel = destination / "keep"
    sentinel.write_text("original")
    with pytest.raises(FileExistsError, match="fresh"):
        restore_run(run, destination)
    assert sentinel.read_text() == "original"
    (run / "raw_artifacts.tar.gz").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum or size"):
        restore_run(run, tmp_path / "new")
    assert not (tmp_path / "new").exists()


@pytest.mark.parametrize("attack", ["traversal", "symlink", "hardlink", "duplicate", "unexpected", "missing", "bytes"])
def test_restore_checks_tar_members_even_if_container_hash_is_updated(tmp_path, attack):
    run, record_name = _minimal_run(tmp_path)
    pack_run(run, record_name)

    def transform(members):
        member, data = members[0]
        if attack == "traversal":
            member.name = "../escape"
        elif attack == "symlink":
            member.type, member.linkname, member.size = tarfile.SYMTYPE, "../escape", 0
        elif attack == "hardlink":
            member.type, member.linkname, member.size = tarfile.LNKTYPE, "../escape", 0
        elif attack == "duplicate":
            members.append((member, data))
        elif attack == "unexpected":
            extra = tarfile.TarInfo("unexpected")
            extra.size = 1
            members.append((extra, b"x"))
        elif attack == "missing":
            members = members[1:]
        elif attack == "bytes":
            members[0] = (member, bytes([data[0] ^ 1]) + data[1:])
        return members

    _repackage(run, transform)
    destination = tmp_path / "restored"
    with pytest.raises(ValueError):
        restore_run(run, destination)
    assert not destination.exists()
    assert not (tmp_path / "escape").exists()
    assert not list(tmp_path.glob(".site-restore.*"))


def test_archive_entry_cannot_escape_through_manifest_paths(tmp_path):
    run, record_name = _minimal_run(tmp_path)
    pack_run(run, record_name)
    path = run / "raw_archive.json"
    manifest = json.loads(path.read_text())
    manifest["entries"][0]["path"] = "../escape"
    write_json(path, manifest)
    with pytest.raises(ValueError, match="path"):
        restore_run(run, tmp_path / "restored")


def test_restore_verifies_the_original_record_hash(tmp_path):
    run, record_name = _minimal_run(tmp_path)
    pack_run(run, record_name)
    path = run / "raw_archive.json"
    manifest = json.loads(path.read_text())
    manifest["record_sha256"] = "0" * 64
    write_json(path, manifest)
    with pytest.raises(ValueError, match="record identity"):
        restore_run(run, tmp_path / "restored")


def _protein_pdb():
    return "ATOM      1   CA ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n"


def test_restored_p2rank_raw_parses_with_separately_restored_exact_input(tmp_path):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    receptor = prepared / "receptor.pdb"
    receptor.write_text(_protein_pdb())
    bundle(prepared, tmp_path / "input_bundle")
    run = tmp_path / "p2rank"
    raw = run / "raw"
    raw.mkdir(parents=True)
    (raw / "receptor.pdb_predictions.csv").write_text(
        "rank,score,probability,center_x,center_y,center_z,residue_ids,surf_atom_ids,sas_points\n"
        "1,2.0,0.5,1.0,2.0,3.0,A_1,1,10\n")
    (raw / "receptor.pdb_residues.csv").write_text("chain,residue_label,residue_name,pocket\nA,1,ALA,1\n")
    record_name = _record(run, method="p2rank")
    record = json.loads((run / record_name).read_text())
    record["input_sha256"] = sha256(receptor)
    write_json(run / record_name, record)
    pack_run(run, record_name)
    restored = tmp_path / "restored_run"
    restore_run(run, restored)
    restored_input = tmp_path / "restored_input"
    restore(tmp_path / "input_bundle", restored_input)
    assert sha256(restored_input / "receptor.pdb") == record["input_sha256"]
    predictions = parse_p2rank_predictions(restored / "raw/receptor.pdb_predictions.csv",
                                          restored / "raw/receptor.pdb_residues.csv", restored_input / "receptor.pdb")
    assert predictions[0]["residues"] == [{"chain": "A", "number": 1, "insertion_code": "", "resname": "ALA"}]


def test_restored_fpocket_raw_parses_and_retains_exact_embedded_input(tmp_path):
    run = tmp_path / "fpocket"
    pockets = run / "receptor_out/pockets"
    pockets.mkdir(parents=True)
    receptor = run / "receptor.pdb"
    receptor.write_text(_protein_pdb())
    (pockets / "pocket1_vert.pqr").write_text(
        "HEADER 0 - Pocket Score : 0.5\nHEADER 2 - Number of V. Vertices : 1\n"
        "ATOM      1    O STP     1       2.000   3.000   4.000    0.00     4.35\n")
    (pockets / "pocket1_atm.cif").write_text("""data_pocket
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 CA ALA A 1 ? 1.0 2.0 3.0 1 A
#
""")
    record_name = _record(run)
    record = json.loads((run / record_name).read_text())
    record["input_sha256"] = sha256(receptor)
    write_json(run / record_name, record)
    pack_run(run, record_name)
    restored = tmp_path / "restored_run"
    restore_run(run, restored)
    assert sha256(restored / "receptor.pdb") == record["input_sha256"]
    prediction, = parse_fpocket(restored / "receptor_out", restored / "receptor.pdb")
    np.testing.assert_array_equal(prediction["center_A"], [2, 3, 4])
    assert prediction["residue_keys"] == [["A", 1, "", "ALA"]]
