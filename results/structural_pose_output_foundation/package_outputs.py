#!/usr/bin/env python3
"""Approved opaque structural-output packaging; standard library only."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tarfile
import tempfile
import unicodedata

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOTS = {"main": "results/structural_pose_v3", "assembly": "results/structural_assembly_pose_v3",
                "haddock": "results/haddock_paired_cohort"}
LIMITS = {"max_raw_bytes": 70_000_000, "max_files": 1000, "max_gzip_bytes": 100_000_000}
HISTORICAL = {
    "archive": {"path": "results/haddock_paired_cohort/orchestration_recovery/completed_first_job_raw.tar.gz",
                "sha256": "c15799f6e92d38cf98d4d4c509a367ce14633faeb6a434dd5d0156ba23895ed7", "bytes": 3_254_709},
    "manifest": {"path": "results/haddock_paired_cohort/orchestration_recovery/completed_first_job_manifest.json",
                 "sha256": "7bbae28624ea8556e6b903ef638598e996f7fdb6e4695d7e9e4aff49af043ac2", "bytes": 55_318},
    "member_count": 240, "target_prefix": SOURCE_ROOTS["haddock"],
}
READINESS_COUNTS = {"main": 216, "assembly": 36, "haddock": 108, "union_sources": 216, "pools": 108}
INSTALLED_COMPONENTS = {"tools", ".venv", "venv", "env", "site-packages", "dist-packages", "node_modules"}
CACHE_COMPONENTS = {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
EXECUTABLE_MAGICS = (b"\x7fELF", b"MZ", b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
                     b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca", b"\xca\xfe\xba\xbf",
                     b"\xbf\xba\xfe\xca", b"\x00asm")
CHUNK = 1024 * 1024
MUST_BE_ABSENT = ["results/haddock_paired_cohort/orchestration_recovery/coordinator.lock"]


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def snapshot_sha256(snapshot):
    return digest(encoded(snapshot))


def _json(raw):
    def unique(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("Duplicate JSON key")
            value[key] = item
        return value
    def constant(value):
        raise ValueError("Nonfinite JSON value: " + value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=constant)


def portable_path(name):
    if (not isinstance(name, str) or not name or "\\" in name or ":" in name
            or unicodedata.normalize("NFC", name) != name or any(ord(c) < 32 or ord(c) == 127 for c in name)):
        raise ValueError("Use a portable repository-relative path")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != name or name == ".":
        raise ValueError("Use a canonical repository-relative path")
    devices = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)), *(f"lpt{i}" for i in range(1, 10))}
    if any(p.rstrip(". ") != p or p.casefold() == ".git" or p.split(".")[0].casefold() in devices for p in path.parts):
        raise ValueError("Unsafe portable path component")
    return name


def _sha(value):
    if not isinstance(value, str) or re.fullmatch("[0-9a-f]{64}", value) is None:
        raise ValueError("Expected a lowercase SHA256")
    return value


def _integer(value, *, positive=False):
    if type(value) is not int or value < (1 if positive else 0):
        raise ValueError("Expected a nonnegative integer" if not positive else "Expected a positive integer")
    return value


def _pin(record):
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise ValueError("Pins require exact path/bytes/SHA256 fields")
    portable_path(record["path"]); _integer(record["bytes"]); _sha(record["sha256"])
    return record


def _inside_sources(name):
    return any(name == prefix or name.startswith(prefix + "/") for prefix in SOURCE_ROOTS.values())


def _namespace(names):
    nodes = {}
    for name, final_kind in names:
        portable_path(name)
        parts = PurePosixPath(name).parts
        for index in range(1, len(parts) + 1):
            prefix = PurePosixPath(*parts[:index]).as_posix()
            kind = final_kind if index == len(parts) else "directory"
            prior = nodes.get(prefix.casefold())
            if prior is not None and prior != (prefix, kind):
                raise ValueError("Case or file/directory path collision")
            nodes[prefix.casefold()] = (prefix, kind)


def _exclusion(name, rules):
    matches = [r for r in rules if name == r["path"] or
               (r.get("scope", "file") == "tree" and name.startswith(r["path"] + "/"))]
    if len(matches) > 1:
        raise ValueError("Overlapping exclusion rules")
    return matches[0] if matches else None


def _source_path_policy(name, exclusion=None):
    parts = [p.casefold() for p in PurePosixPath(name).parts]
    if any("supplemental" in p or "worksharing" in p for p in parts):
        raise ValueError("Protected supplemental/worksharing path")
    if any(p in INSTALLED_COMPONENTS for p in parts) or parts[:2] == ["data", "processed"]:
        raise ValueError("Installed software and processed trees are forbidden")
    if (any(p in CACHE_COMPONENTS for p in parts) or parts[-1].endswith((".pyc", ".pyo"))
            or parts[-1] == ".ds_store") and (exclusion is None or exclusion["kind"] != "cache"):
        raise ValueError("Cache files require a finite explicit cache exclusion")


@contextmanager
def ordinary_file(root, name):
    """Open through stable no-follow descriptors for every path component."""
    portable_path(name)
    if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
        raise ValueError("No-follow directory descriptors are required")
    descriptors = []
    try:
        fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW); descriptors.append(fd)
        parts = PurePosixPath(name).parts
        for part in parts[:-1]:
            fd = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd); descriptors.append(fd)
        fd = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_NONBLOCK", 0), dir_fd=fd); descriptors.append(fd)
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            raise ValueError("Source must be an ordinary file: " + name)
        with os.fdopen(os.dup(fd), "rb") as stream:
            yield stream
    except OSError as error:
        raise ValueError("Missing, linked or nonordinary source: " + name) from error
    finally:
        for fd in reversed(descriptors):
            os.close(fd)


def _executable(name, prefix):
    return prefix.startswith(EXECUTABLE_MAGICS) or PurePosixPath(name).suffix.casefold() in {".so", ".dylib", ".dll", ".exe", ".class"}


def _record(root, name, *, reject_executable=False, max_bytes=None):
    hasher = hashlib.sha256(); size = 0
    with ordinary_file(root, name) as stream:
        before = os.fstat(stream.fileno())
        if max_bytes is not None and before.st_size > max_bytes:
            raise ValueError("A single file exceeds the raw-byte chunk limit")
        prefix = stream.read(16)
        if reject_executable and _executable(name, prefix):
            raise ValueError("Installed executable bytes are forbidden: " + name)
        hasher.update(prefix); size += len(prefix)
        for chunk in iter(lambda: stream.read(CHUNK), b""):
            hasher.update(chunk); size += len(chunk)
        after = os.fstat(stream.fileno())
        if (size != before.st_size or size != after.st_size
                or before.st_mtime_ns != after.st_mtime_ns or before.st_ctime_ns != after.st_ctime_ns):
            raise ValueError("Source changed while hashing: " + name)
    return {"path": name, "bytes": size, "sha256": hasher.hexdigest()}


def _read(root, name, *, reject_executable=False):
    with ordinary_file(root, name) as stream:
        raw = stream.read()
    if reject_executable and _executable(name, raw[:16]):
        raise ValueError("Installed executable bytes are forbidden: " + name)
    return raw


def _checked(root, record):
    _pin(record); raw = _read(root, record["path"])
    if len(raw) != record["bytes"] or digest(raw) != record["sha256"]:
        raise ValueError("Pinned bytes changed: " + record["path"])
    return raw


def validate_plan(plan, *, allow_synthetic_fixture=False):
    keys = {"schema_version", "profile", "source_roots", "destination", "limits", "exclusions", "historical_first_job",
            "source_contracts", "source_protocol_commit", "must_be_absent"}
    if not isinstance(plan, dict) or set(plan) != keys or type(plan["schema_version"]) is not int or plan["schema_version"] != 1:
        raise ValueError("Unexpected opaque packaging plan schema")
    synthetic = plan["profile"] == "synthetic_test_fixture"
    if plan["profile"] != "real_frozen_structural_outputs" and not (synthetic and allow_synthetic_fixture):
        raise ValueError("Real frozen structural profile required")
    if plan["source_roots"] != SOURCE_ROOTS:
        raise ValueError("Exactly three fixed disjoint structural roots are required")
    portable_path(plan["destination"])
    if _inside_sources(plan["destination"]) or any(p.casefold() in INSTALLED_COMPONENTS for p in PurePosixPath(plan["destination"]).parts):
        raise ValueError("Destination must be outside source and installed trees")
    limits = plan["limits"]
    if not isinstance(limits, dict) or set(limits) != set(LIMITS):
        raise ValueError("All three packaging limits are required")
    for key, maximum in LIMITS.items():
        _integer(limits[key], positive=True)
        if limits[key] > maximum or (not synthetic and limits[key] != maximum):
            raise ValueError("Packaging limits exceed or change the reviewed policy")
    rules = plan["exclusions"]
    if not isinstance(rules, list):
        raise ValueError("Use an explicit finite exclusion list")
    for row in rules:
        if (not isinstance(row, dict) or not {"path", "kind", "reason"} <= set(row) or set(row) - {"path", "kind", "reason", "scope"}
                or row.get("scope", "file") not in {"file", "tree"} or row["kind"] not in {"derived", "cache"}
                or not isinstance(row["reason"], str) or not row["reason"].strip()):
            raise ValueError("Exclusions require exact path, file/tree scope, derived/cache kind and reason")
        portable_path(row["path"])
        if not _inside_sources(row["path"]) or any(c in row["path"] for c in "*?[") or row["path"] in SOURCE_ROOTS.values():
            raise ValueError("Exclusions must name an exact subordinate source path")
        _source_path_policy(row["path"], row)
    for index, row in enumerate(rules):
        if any(row["path"] == other["path"] or
               (other.get("scope", "file") == "tree" and row["path"].startswith(other["path"] + "/"))
               for j, other in enumerate(rules) if index != j):
            raise ValueError("Duplicate or overlapping exclusion rules")
    history = plan["historical_first_job"]
    if not isinstance(history, dict) or set(history) != set(HISTORICAL):
        raise ValueError("Historical first-job archive/manifest contract required")
    _pin(history["archive"]); _pin(history["manifest"]); _integer(history["member_count"], positive=True)
    if history["target_prefix"] != SOURCE_ROOTS["haddock"] or (not synthetic and history != HISTORICAL):
        raise ValueError("All240 original historical first-job pins are immutable")
    if any(not history[key]["path"].startswith(SOURCE_ROOTS["haddock"] + "/") for key in ["archive", "manifest"]):
        raise ValueError("Historical evidence must remain in the HADDOCK root")
    if not isinstance(plan["source_protocol_commit"], str) or re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", plan["source_protocol_commit"]) is None:
        raise ValueError("Use the exact published source/protocol commit")
    contracts = plan["source_contracts"]
    if not isinstance(contracts, list) or not contracts:
        raise ValueError("A bounded source/protocol pin list is required")
    for row in contracts:
        _pin(row); _source_path_policy(row["path"])
    if len({r["path"] for r in contracts}) != len(contracts):
        raise ValueError("Duplicate source/protocol contract")
    _namespace((r["path"], "file") for r in contracts)
    if plan["must_be_absent"] != MUST_BE_ABSENT:
        raise ValueError("The stopped coordinator lock must be explicitly absent")
    return plan


def _fixture_root(root):
    bases = {Path(tempfile.gettempdir()).resolve(), Path("/tmp").resolve()}
    if not any(root.is_relative_to(base) and root != base for base in bases):
        raise ValueError("Synthetic fixtures require a private TEMP root")
    if any((p / ".git").exists() or (p / ".git").is_symlink() or
           ((p / "HEAD").is_file() and (p / "config").is_file() and (p / "objects").is_dir()) for p in (root, *root.parents)):
        raise ValueError("Synthetic fixtures cannot use a Git repository")


def checked_root(root):
    """Reject supplied path aliases before canonicalization can hide them."""
    requested = Path(root).absolute()
    if any(path.is_symlink() for path in (requested, *requested.parents)):
        raise ValueError("Supplied source root traverses a symlink")
    if not requested.is_dir():
        raise ValueError("Source root must be an existing ordinary directory")
    return requested.resolve(strict=True)


def inventory(root, plan):
    """Snapshot all ordinary files, including every explicitly excluded byte."""
    root = checked_root(root); included, excluded, names = [], [], []
    for relative in SOURCE_ROOTS.values():
        source = root / relative
        for ancestor in (source, *source.parents):
            if ancestor.is_symlink():
                raise ValueError("Linked source directory")
            if ancestor == root:
                break
        if not source.is_dir():
            raise ValueError("Missing source root: " + relative)
        count = 0
        def walk_error(error):
            raise ValueError("Unreadable source directory") from error
        for directory, dirs, files in os.walk(source, topdown=True, followlinks=False, onerror=walk_error):
            directory_path = Path(directory); names.append((directory_path.relative_to(root).as_posix(), "directory"))
            for name in dirs:
                child = directory_path / name
                if child.is_symlink() or not stat.S_ISDIR(child.lstat().st_mode):
                    raise ValueError("Linked or nonordinary source directory")
                child_name = child.relative_to(root).as_posix(); portable_path(child_name)
                if any(p.casefold() in INSTALLED_COMPONENTS or "supplemental" in p.casefold() or "worksharing" in p.casefold()
                       for p in PurePosixPath(child_name).parts):
                    raise ValueError("Forbidden installed/protected source directory")
            for name in files:
                relative_name = (directory_path / name).relative_to(root).as_posix(); portable_path(relative_name)
                rule = _exclusion(relative_name, plan["exclusions"]); _source_path_policy(relative_name, rule)
                row = _record(root, relative_name, reject_executable=True,
                              max_bytes=plan["limits"]["max_raw_bytes"] if rule is None else None)
                names.append((relative_name, "file")); count += 1
                if rule is None:
                    included.append(row)
                else:
                    excluded.append({**row, "kind": rule["kind"], "reason": rule["reason"],
                                     "scope": rule.get("scope", "file"), "rule_path": rule["path"]})
        if not count:
            raise ValueError("Empty ordinary-file source root: " + relative)
    _namespace(names)
    if not included:
        raise ValueError("No included source files")
    return {"included": sorted(included, key=lambda r: r["path"]), "excluded": sorted(excluded, key=lambda r: r["path"])}


def partition(records, limits):
    """Greedy sorted whole-file chunks; only paths and sizes determine layout."""
    for key in LIMITS:
        _integer(limits[key], positive=True)
        if limits[key] > LIMITS[key]:
            raise ValueError("Chunk limit exceeds the reviewed maximum")
    rows = sorted(records, key=lambda r: r["path"])
    if len({r["path"] for r in rows}) != len(rows):
        raise ValueError("Duplicate chunk input path")
    _namespace((r["path"], "file") for r in rows)
    chunks, current, size = [], [], 0
    for row in rows:
        _integer(row["bytes"]); _sha(row["sha256"])
        if row["bytes"] > limits["max_raw_bytes"]:
            raise ValueError("A single file exceeds the raw-byte chunk limit")
        if current and (size + row["bytes"] > limits["max_raw_bytes"] or len(current) >= limits["max_files"]):
            chunks.append(current); current, size = [], 0
        current.append(row); size += row["bytes"]
    if current:
        chunks.append(current)
    if not chunks:
        raise ValueError("No chunk records")
    return chunks


def verify_historical(root, plan, snapshot):
    """Check original full-root paths and opaque historical archive bytes."""
    history = plan["historical_first_job"]
    included = {r["path"]: r for r in snapshot["included"]}
    for key in ["archive", "manifest"]:
        if included.get(history[key]["path"]) != history[key]:
            raise ValueError("Historical archive/manifest must be included unchanged")
    original = _json(_checked(root, history["manifest"])); rows = original.get("files")
    if not isinstance(rows, list) or len(rows) != history["member_count"]:
        raise ValueError("Historical manifest member count changed")
    expected = {}
    for row in rows:
        _pin(row)
        if not row["path"].startswith(history["target_prefix"] + "/") or row["path"] in expected:
            raise ValueError("Historical full-root path is missing or duplicate")
        if included.get(row["path"]) != row:
            raise ValueError("Historical path is absent, excluded or changed: " + row["path"])
        expected[row["path"]] = row
    _namespace((name, "file") for name in expected)
    seen = set()
    with tarfile.open(fileobj=io.BytesIO(_checked(root, history["archive"])), mode="r:gz") as archive:
        for member in archive:
            portable_path(member.name)
            name = history["target_prefix"] + "/" + member.name
            if (member.type not in (tarfile.REGTYPE, tarfile.AREGTYPE) or member.linkname
                    or getattr(member, "sparse", None) is not None or name in seen or name not in expected):
                raise ValueError("Unexpected or nonordinary historical archive member")
            payload = archive.extractfile(member).read()
            if member.size != expected[name]["bytes"] or len(payload) != member.size or digest(payload) != expected[name]["sha256"]:
                raise ValueError("Historical archive member bytes changed")
            seen.add(name)
    if seen != set(expected):
        raise ValueError("Historical archive is incomplete")
    return {"archive": history["archive"], "manifest": history["manifest"], "exact_original_paths_verified": len(seen)}


def _metadata_name(root, path):
    path = Path(path)
    if path.is_absolute():
        try:
            name = path.relative_to(root).as_posix()
        except ValueError as error:
            raise ValueError("Packaging metadata must belong to the source repository") from error
    else:
        name = path.as_posix()
    portable_path(name)
    if _inside_sources(name):
        raise ValueError("Approval, plan and readiness must be outside the three source trees")
    return name


def _ensure_absent(root, names):
    for name in names:
        path = root / name
        if any(p.is_symlink() for p in (path, *path.parents)) or path.exists():
            raise ValueError("Stopped-writer lock must remain absent: " + name)


def _authorization(root, plan_name, approval_name, *, allow_synthetic_fixture=False):
    plan_raw = _read(root, plan_name); plan = validate_plan(_json(plan_raw), allow_synthetic_fixture=allow_synthetic_fixture)
    if plan["profile"] == "synthetic_test_fixture":
        _fixture_root(root)
    approval_raw = _read(root, approval_name); approval = _json(approval_raw)
    plan_pin = {"path": plan_name, "bytes": len(plan_raw), "sha256": digest(plan_raw)}
    if (type(approval.get("schema_version")) is not int or approval["schema_version"] != 1
            or approval.get("status") != "approved_terminal_structural_raw_packaging"
            or approval.get("packaging_authorized") is not True or approval.get("packaging_plan") != plan_pin
            or approval.get("packager_sha256") != digest(Path(__file__).read_bytes())):
        raise ValueError("Explicit root packaging approval and exact packager/plan pins are required")
    validation_pin = approval.get("terminal_validation"); _pin(validation_pin)
    _metadata_name(root, validation_pin["path"])
    validation = _json(_checked(root, validation_pin))
    if (type(validation.get("schema_version")) is not int or validation["schema_version"] != 1
            or validation.get("status") != "all_structural_sources_terminal_metadata_only"
            or validation.get("accuracy_computed") is not False or validation.get("numeric_pose_or_energy_payloads_interpreted") is not False
            or validation.get("writers_stopped") is not True or validation.get("pending") != []):
        raise ValueError("External metadata-only stopped-writer terminal readiness is required")
    # The external status also attests unique union pose IDs equal score-record
    # keys. Matching counters alone do not prove that; no reader is called here.
    counts = validation.get("counts")
    if not isinstance(counts, dict) or set(counts) != set(READINESS_COUNTS) | {"retained", "recorded"}:
        raise ValueError("All fixed terminal readiness counts are required")
    for value in counts.values():
        _integer(value)
    if any(counts[k] != v for k, v in READINESS_COUNTS.items()) or counts["retained"] != counts["recorded"]:
        raise ValueError("Terminal readiness denominator or retained/recorded binding changed")
    _sha(validation.get("source_snapshot_sha256")); _ensure_absent(root, plan["must_be_absent"])
    for record in plan["source_contracts"]:
        raw = _checked(root, record)
        if _executable(record["path"], raw[:16]):
            raise ValueError("Source contracts cannot contain installed executable bytes")
    gate = {"packaging_plan": plan_pin,
            "root_approval": {"path": approval_name, "bytes": len(approval_raw), "sha256": digest(approval_raw)},
            "terminal_validation": validation_pin, "source_snapshot_sha256": validation["source_snapshot_sha256"],
            "source_contracts": plan["source_contracts"], "source_protocol_commit": plan["source_protocol_commit"],
            "must_be_absent": plan["must_be_absent"], "writers_stopped": True, "counts": counts}
    return plan, gate


def _write_chunk(root, stage, rows, number, total, plan, snapshot):
    name = f"structural_outputs_{number:04d}.tar.gz"; target = stage / name
    with target.open("xb") as output:
        with gzip.GzipFile(fileobj=output, filename="", mode="wb", mtime=0, compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                for row in rows:
                    raw = _read(root, row["path"], reject_executable=True)
                    if len(raw) != row["bytes"] or digest(raw) != row["sha256"]:
                        raise ValueError("Source changed during packaging: " + row["path"])
                    info = tarfile.TarInfo(row["path"])
                    info.size, info.mode, info.uid, info.gid, info.mtime = len(raw), 0o644, 0, 0, 0
                    info.uname = info.gname = ""
                    archive.addfile(info, io.BytesIO(raw))
    if target.stat().st_size >= plan["limits"]["max_gzip_bytes"]:
        raise ValueError("Compressed archive reaches the exclusive size limit")
    expected = {r["path"]: r for r in rows}; seen = set()
    with tarfile.open(target, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or member.name in seen or member.name not in expected:
                raise ValueError("Staged archive member set changed")
            raw = archive.extractfile(member).read(); row = expected[member.name]
            if member.size != row["bytes"] or len(raw) != row["bytes"] or digest(raw) != row["sha256"]:
                raise ValueError("Staged opaque archive bytes changed")
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError("Staged archive is incomplete")
    archive_pin = {"path": name, "bytes": target.stat().st_size, "sha256": digest(target.read_bytes())}
    manifest = {"schema_version": 1, "archive": name, "archive_sha256": archive_pin["sha256"], "target_prefix": "",
                "files": rows, "chunk_number": number, "total_chunks": total, "source_roots": SOURCE_ROOTS,
                "source_snapshot_sha256": snapshot_sha256(snapshot), "excluded_files": snapshot["excluded"],
                "exclusion_rules": plan["exclusions"], "raw_member_bytes": sum(r["bytes"] for r in rows), "numeric_payloads_interpreted": False}
    manifest_name = f"structural_outputs_{number:04d}_manifest.json"; raw = encoded(manifest)
    with (stage / manifest_name).open("xb") as output:
        output.write(raw)
    return {"archive": archive_pin, "manifest": {"path": manifest_name, "bytes": len(raw), "sha256": digest(raw)},
            "target_prefix": "", "files": len(rows), "uncompressed_bytes": manifest["raw_member_bytes"]}


def _identity(value):
    return value.st_dev, value.st_ino


def _destination_identity(destination, expected):
    value = destination.lstat()
    if not stat.S_ISDIR(value.st_mode) or _identity(value) != expected:
        raise ValueError("Claimed publication directory was replaced")


def _publish_file(source, directory_fd, created, expected):
    name = source.name
    fd = os.open(name, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o644, dir_fd=directory_fd)
    record = {"name": name, "identity": _identity(os.fstat(fd)), "bytes": 0, "sha256": digest(b"")}; created.append(record)
    hasher = hashlib.sha256()
    try:
        with source.open("rb") as stream:
            for chunk in iter(lambda: stream.read(CHUNK), b""):
                offset = 0
                while offset < len(chunk):
                    written = os.write(fd, chunk[offset:])
                    if written <= 0:
                        raise OSError("Publication write made no progress")
                    hasher.update(chunk[offset:offset + written]); offset += written
                    record["bytes"] += written; record["sha256"] = hasher.hexdigest()
        os.lseek(fd, 0, os.SEEK_SET)
        with os.fdopen(os.dup(fd), "rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if (actual != record["sha256"] or actual != expected["sha256"]
                or record["bytes"] != expected["bytes"]):
            raise ValueError("Published bytes differ from verified staging")
    finally:
        os.close(fd)


def _verify_published(directory_fd, created, expected):
    """Bind visible entries to the exact inodes and bytes published here."""
    if len(created) != len(expected) or {r["name"] for r in created} != set(expected):
        raise ValueError("Published artifact identity set changed")
    for row in created:
        pin = expected[row["name"]]
        try:
            entry = os.stat(row["name"], dir_fd=directory_fd, follow_symlinks=False)
            if (not stat.S_ISREG(entry.st_mode) or _identity(entry) != row["identity"]
                    or entry.st_size != row["bytes"] or entry.st_size != pin["bytes"]):
                raise ValueError("Published directory entry changed: " + row["name"])
            fd = os.open(row["name"], os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_fd)
            with os.fdopen(fd, "rb") as stream:
                opened = os.fstat(stream.fileno())
                if _identity(opened) != row["identity"]:
                    raise ValueError("Published entry changed while opening")
                actual = hashlib.file_digest(stream, "sha256").hexdigest()
                after = os.fstat(stream.fileno())
            visible = os.stat(row["name"], dir_fd=directory_fd, follow_symlinks=False)
            if (not stat.S_ISREG(visible.st_mode) or _identity(visible) != row["identity"]
                    or visible.st_size != pin["bytes"] or opened.st_size != after.st_size
                    or opened.st_mtime_ns != after.st_mtime_ns or opened.st_ctime_ns != after.st_ctime_ns
                    or visible.st_mtime_ns != after.st_mtime_ns or visible.st_ctime_ns != after.st_ctime_ns
                    or actual != row["sha256"] or actual != pin["sha256"]):
                raise ValueError("Published artifact bytes changed: " + row["name"])
        except OSError as error:
            raise ValueError("Missing or linked published entry: " + row["name"]) from error


def _cleanup(destination, directory_fd, identity, created):
    for row in reversed(created):
        try:
            current = os.stat(row["name"], dir_fd=directory_fd, follow_symlinks=False)
            if not stat.S_ISREG(current.st_mode) or _identity(current) != row["identity"] or current.st_size != row["bytes"]:
                continue
            fd = os.open(row["name"], os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_fd)
            with os.fdopen(fd, "rb") as stream:
                unchanged = hashlib.file_digest(stream, "sha256").hexdigest() == row["sha256"]
            if unchanged:
                os.unlink(row["name"], dir_fd=directory_fd)
        except OSError:
            pass
    try:
        _destination_identity(destination, identity); destination.rmdir()
    except (OSError, ValueError):
        pass


def package(root, plan_path, approval_path, *, allow_synthetic_fixture=False):
    """Package externally approved terminal sources; invoke no native/readers.

    Readiness describes all108 HADDOCK cells. No file-count inference is made:
    the69 source-failed cells need not have physical native-run receipts.
    """
    root = checked_root(root)
    plan_name, approval_name = _metadata_name(root, plan_path), _metadata_name(root, approval_path)
    plan, gate = _authorization(root, plan_name, approval_name, allow_synthetic_fixture=allow_synthetic_fixture)
    destination = root / plan["destination"]
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Packaging destination must be fresh")
    if not destination.parent.is_dir() or any(p.is_symlink() for p in (destination.parent, *destination.parent.parents)):
        raise ValueError("Destination parent must be an existing ordinary directory")
    before = inventory(root, plan)
    if snapshot_sha256(before) != gate["source_snapshot_sha256"]:
        raise ValueError("External readiness does not bind included and excluded source bytes")
    history = verify_historical(root, plan, before); chunks = partition(before["included"], plan["limits"])
    with tempfile.TemporaryDirectory(prefix=".structural-package-", dir=destination.parent) as temporary:
        stage = Path(temporary); outputs = []
        for index, rows in enumerate(chunks):
            outputs.append(_write_chunk(root, stage, rows, index + 1, len(chunks), plan, before))
            if index == (len(chunks) - 1) // 2 and inventory(root, plan) != before:
                raise ValueError("Complete source snapshot changed during packaging")
        if _authorization(root, plan_name, approval_name, allow_synthetic_fixture=allow_synthetic_fixture) != (plan, gate):
            raise ValueError("Packaging authorization changed")
        expected_outputs = {pin["path"]: pin for row in outputs for pin in [row["archive"], row["manifest"]]}
        for pin in expected_outputs.values():
            _checked(stage, pin)
        destination.mkdir(mode=0o700)
        identity = _identity(destination.lstat())
        try:
            directory_fd = os.open(destination, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        except OSError:
            try:
                _destination_identity(destination, identity); destination.rmdir()
            except (OSError, ValueError):
                pass
            raise
        if _identity(os.fstat(directory_fd)) != identity:
            os.close(directory_fd)
            raise ValueError("Claimed publication directory changed before opening")
        created = []
        try:
            for source in sorted(stage.iterdir()):
                _destination_identity(destination, identity)
                _publish_file(source, directory_fd, created, expected_outputs[source.name])
            if inventory(root, plan) != before:
                raise ValueError("Complete source snapshot changed after packaging")
            if _authorization(root, plan_name, approval_name, allow_synthetic_fixture=allow_synthetic_fixture) != (plan, gate):
                raise ValueError("Packaging authorization changed before receipt publication")
            _destination_identity(destination, identity)
            _verify_published(directory_fd, created, expected_outputs)
            receipt = {"schema_version": 1, "stage": "Post-execution opaque packaging; no scientific freeze replaced.",
                       "source_roots": SOURCE_ROOTS, "limits": plan["limits"], "source_snapshot": before,
                       "source_snapshot_sha256": snapshot_sha256(before), "exclusion_rules": plan["exclusions"], "authorization": gate,
                       "archives": outputs, "historical_first_job": history, "all_before_during_after_file_snapshots_equal": True,
                       "all_included_ordinary_files_preserved": True, "receipt_published_last": True,
                       "numeric_pose_or_energy_payloads_interpreted": False, "accuracy_computed": False,
                       "native_executables_used": False, "resource_measurements_invented": False,
                       "synthetic_test_fixture": plan["profile"] == "synthetic_test_fixture", "packager_sha256": digest(Path(__file__).read_bytes())}
            with (stage / "packaging_receipt.json").open("xb") as output:
                output.write(encoded(receipt))
            receipt_path = stage / "packaging_receipt.json"
            receipt_pin = {"path": receipt_path.name, "bytes": receipt_path.stat().st_size, "sha256": digest(receipt_path.read_bytes())}
            _publish_file(receipt_path, directory_fd, created, receipt_pin)
            _verify_published(directory_fd, created, {**expected_outputs, receipt_pin["path"]: receipt_pin})
            _destination_identity(destination, identity)
            return receipt
        except BaseException:
            _cleanup(destination, directory_fd, identity, created)
            raise
        finally:
            os.close(directory_fd)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--approval", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(package(args.root, args.plan, args.approval), indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
