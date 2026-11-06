"""Deterministic, portable bundles of terminal site-method attempts.

The original execution/result record is archived verbatim. Only its declared
regular-file artifacts are packed; no command or path from that record is
executed during restoration. Scientific tools are not needed to verify a run.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tarfile
import tempfile

from .provenance import sha256, write_json


ARCHIVE_NAME = "raw_artifacts.tar.gz"
MANIFEST_NAME = "raw_archive.json"
TERMINAL_STATUSES = frozenset({"completed", "failed", "ok", "timeout", "cancelled", "canceled", "error"})


def _relative(value):
    """Require one canonical, portable relative filename without traversal."""
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value or ":" in value:
        raise ValueError("archive entry must have a portable relative path")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {".", ".."} for part in path.parts) or path.as_posix() != value:
        raise ValueError("archive entry path traversal or noncanonical path")
    if value in {ARCHIVE_NAME, MANIFEST_NAME}:
        raise ValueError("an archive cannot contain its own control files")
    return value


def _digest(value):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError("archive entry requires a lowercase SHA256 digest")
    return value


def _size(value):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError("archive entry requires a nonnegative integer size")
    return value


def _record_name(value):
    value = _relative(value)
    if "/" in value:
        raise ValueError("record_name must be a filename in the run directory")
    return value


def _regular_file(root, relative):
    """Reject symlinks in every component, including links out of the run."""
    current = Path(root)
    if not stat.S_ISDIR(current.lstat().st_mode):
        raise ValueError("run directory must be an actual directory, not a symlink")
    parts = PurePosixPath(relative).parts
    for index, part in enumerate(parts):
        current = current / part
        mode = current.lstat().st_mode
        if index == len(parts) - 1:
            if not stat.S_ISREG(mode):
                raise ValueError(f"artifact must be a regular file, not a symlink: {relative}")
        elif not stat.S_ISDIR(mode):
            raise ValueError(f"artifact parent must be an actual directory, not a symlink: {relative}")
    return current


def _artifact_entries(record, record_name):
    if record.get("status") not in TERMINAL_STATUSES:
        raise ValueError("only a terminal execution record can be archived")
    artifacts = record.get("artifacts")
    if isinstance(artifacts, dict):
        raw = [{"path": path, "sha256": value["sha256"], "bytes": value["bytes"]}
               for path, value in artifacts.items()]
    elif isinstance(artifacts, list):
        raw = [{"path": value["path"], "sha256": value["sha256"], "bytes": value["size_bytes"]}
               for value in artifacts]
    else:
        raise ValueError("record must contain an artifact dictionary or list")
    entries = {}
    for value in raw:
        name = _relative(value["path"])
        if name == record_name or name in entries:
            raise ValueError("duplicate or self-referential record artifact")
        entries[name] = {"path": name, "sha256": _digest(value["sha256"]), "bytes": _size(value["bytes"]),
                         "role": "artifact"}
    return entries


def _validate_entries(values):
    if not isinstance(values, list):
        raise ValueError("archive entry inventory must be a list")
    entries = {}
    for value in values:
        name = _relative(value["path"])
        if name in entries:
            raise ValueError("duplicate archive entry path")
        if value.get("role") not in {"record", "artifact"}:
            raise ValueError("archive entry has an unknown role")
        entries[name] = {"path": name, "sha256": _digest(value["sha256"]), "bytes": _size(value["bytes"]),
                         "role": value["role"]}
    for name in entries:
        if any(parent.as_posix() in entries for parent in PurePosixPath(name).parents if str(parent) != "."):
            raise ValueError("archive file is also a parent of another entry")
    return entries


def _verify_file(path, entry):
    if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
        raise ValueError(f"artifact checksum or size mismatch: {entry['path']}")


def _inspect_tar(archive_path, entries, destination=None):
    """Validate every member before copying its bytes; never use extractall."""
    observed = set()
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            name = _relative(member.name)
            if member.type != tarfile.REGTYPE or member.issym() or member.islnk() or member.sparse is not None:
                raise ValueError(f"archive contains a nonregular file or link: {name}")
            if name not in entries or name in observed:
                raise ValueError(f"archive contains an unexpected or duplicate entry: {name}")
            expected = entries[name]
            if member.size != expected["bytes"]:
                raise ValueError(f"archive member size mismatch: {name}")
            observed.add(name)
            digest, count = hashlib.sha256(), 0
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(f"archive member has no file bytes: {name}")
            output = None
            try:
                if destination is not None:
                    path = Path(destination).joinpath(*PurePosixPath(name).parts)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    output = path.open("xb")
                with source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        digest.update(chunk)
                        count += len(chunk)
                        if count > expected["bytes"]:
                            raise ValueError(f"archive member exceeds its declared size: {name}")
                        if output is not None:
                            output.write(chunk)
            finally:
                if output is not None:
                    output.close()
            if count != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
                raise ValueError(f"archive member checksum mismatch: {name}")
    if observed != set(entries):
        raise ValueError("archive is missing declared artifacts")


def pack_run(run_dir, record_name):
    """Verify and bundle every declared artifact plus the unchanged record.

    fpocket ``execution.json`` dictionaries and P2Rank ``result.json`` lists
    are supported. The record's external input path is metadata only: P2Rank's
    receptor must be retained in the separate input archive identified by its
    recorded SHA256. Existing bundles are verified and reused, never replaced.
    """
    run_dir, record_name = Path(run_dir), _record_name(record_name)
    record_path = _regular_file(run_dir, record_name)
    record_bytes = record_path.read_bytes()
    record = json.loads(record_bytes)
    entries = _artifact_entries(record, record_name)
    entries[record_name] = {"path": record_name, "sha256": hashlib.sha256(record_bytes).hexdigest(),
                            "bytes": len(record_bytes), "role": "record"}
    entries = _validate_entries([entries[name] for name in sorted(entries)])
    for name, entry in entries.items():
        _verify_file(_regular_file(run_dir, name), entry)
    archive_path, manifest_path = run_dir / ARCHIVE_NAME, run_dir / MANIFEST_NAME
    if archive_path.exists() or archive_path.is_symlink() or manifest_path.exists() or manifest_path.is_symlink():
        if not archive_path.exists() or not manifest_path.exists():
            raise FileExistsError("partial prior bundle exists; preserve it before choosing a fresh run directory")
        existing = _load_manifest(run_dir)
        if existing["record_name"] != record_name or existing["record_sha256"] != entries[record_name]["sha256"] or existing["entries"] != list(entries.values()):
            raise ValueError("existing bundle differs from the terminal record or artifact inventory")
        _inspect_tar(archive_path, entries)
        return existing
    descriptor, temporary = tempfile.mkstemp(prefix=".raw_artifacts.", suffix=".tmp", dir=run_dir)
    os.close(descriptor)
    temporary = Path(temporary)
    try:
        with temporary.open("wb") as stream:
            with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0, compresslevel=9) as compressed:
                with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                    for name, entry in entries.items():
                        info = tarfile.TarInfo(name)
                        info.size, info.mode, info.mtime = entry["bytes"], 0o644, 0
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        with _regular_file(run_dir, name).open("rb") as source:
                            archive.addfile(info, source)
        # Detect source mutations between validation and copying as well as any
        # serialization defect before publishing the archive's checksum.
        _inspect_tar(temporary, entries)
        manifest = {
            "schema_version": 1, "archive": ARCHIVE_NAME, "archive_sha256": sha256(temporary),
            "archive_bytes": temporary.stat().st_size, "format": "deterministic tar+gzip; regular files only",
            "record_name": record_name, "record_sha256": entries[record_name]["sha256"],
            "method": record.get("method"), "status": record["status"],
            "input_sha256": record.get("input_sha256"), "entries": list(entries.values()),
            "uncompressed_bytes": sum(entry["bytes"] for entry in entries.values()),
        }
        temporary.replace(archive_path)
        write_json(manifest_path, manifest)
        return manifest
    finally:
        temporary.unlink(missing_ok=True)


def _load_manifest(run_dir):
    manifest = json.loads(_regular_file(run_dir, MANIFEST_NAME).read_text())
    if manifest.get("schema_version") != 1 or manifest.get("archive") != ARCHIVE_NAME:
        raise ValueError("unsupported site archive manifest")
    _record_name(manifest["record_name"])
    _digest(manifest["record_sha256"])
    _digest(manifest["archive_sha256"])
    archive = _regular_file(run_dir, ARCHIVE_NAME)
    if archive.stat().st_size != _size(manifest["archive_bytes"]) or sha256(archive) != manifest["archive_sha256"]:
        raise ValueError("site archive checksum or size mismatch")
    entries = _validate_entries(manifest["entries"])
    records = [entry for entry in entries.values() if entry["role"] == "record"]
    if len(records) != 1 or records[0]["path"] != manifest["record_name"] or records[0]["sha256"] != manifest["record_sha256"]:
        raise ValueError("archive record identity differs from its manifest")
    if sum(entry["bytes"] for entry in entries.values()) != _size(manifest["uncompressed_bytes"]):
        raise ValueError("archive uncompressed size differs from its inventory")
    return manifest


def restore_run(run_dir, destination):
    """Verify a bundle and restore it to a fresh destination without any tools.

    Files are staged in a sibling temporary directory. An invalid archive never
    leaves partial artifacts in the requested destination or changes originals.
    Recorded absolute tool/input paths are never used for file extraction.
    """
    run_dir, destination = Path(run_dir), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("site archive restoration requires a fresh destination")
    manifest = _load_manifest(run_dir)
    entries = _validate_entries(manifest["entries"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".site-restore.", dir=destination.parent))
    try:
        _inspect_tar(run_dir / ARCHIVE_NAME, entries, temporary)
        record = json.loads((temporary / manifest["record_name"]).read_text())
        declared = _artifact_entries(record, manifest["record_name"])
        restored_artifacts = {name: entry for name, entry in entries.items() if entry["role"] == "artifact"}
        if declared != restored_artifacts or any(record.get(key) != manifest.get(key) for key in ("method", "status", "input_sha256")):
            raise ValueError("restored record differs from the archive artifact inventory or provenance")
        if destination.exists() or destination.is_symlink():
            raise FileExistsError("restoration destination appeared while validating the archive")
        temporary.rename(destination)
        return manifest
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
