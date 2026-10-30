"""Checksummed downloads and atomic, finite JSON records."""

from __future__ import annotations

import hashlib
import gzip
import json
import os
from pathlib import Path
import tempfile
import shutil
from urllib.request import Request, urlopen


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text_atomic(path, data):
    """Replace a text artifact only after all bytes have been written."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(data)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def write_json(path, value):
    """Replace a manifest only when a complete, valid JSON document exists."""
    data = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    write_text_atomic(path, data)


def fetch_verified(url, destination, expected_sha256):
    """Never accept a stale or mutated source, including an existing local file."""
    destination = Path(destination)
    if len(expected_sha256) != 64 or any(c not in "0123456789abcdef" for c in expected_sha256):
        raise ValueError("Expected an explicit lowercase SHA256 digest")
    if not url.startswith("https://"):
        raise ValueError("Source URLs must use HTTPS")
    if destination.exists():
        if sha256(destination) != expected_sha256:
            raise ValueError(f"Checksum mismatch: {destination}; preserve it for inspection")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(fd, "wb") as target:
            with urlopen(Request(url, headers={"User-Agent": "aav-docking-benchmark/0.1"}), timeout=120) as source:
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    target.write(chunk)
        if sha256(temporary) != expected_sha256:
            raise ValueError(f"Downloaded content does not match frozen checksum: {url}")
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return destination


def fetch_manifest(manifest_path, destination):
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    destination = Path(destination)
    for source in manifest["sources"]:
        relative = Path(source["filename"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Source filename must stay inside the raw directory")
        target = destination / relative
        archive_relative = Path(source["archive"]) if source.get("archive") else None
        if archive_relative and (archive_relative.is_absolute() or ".." in archive_relative.parts):
            raise ValueError("Archive path must stay inside the repository")
        archive = manifest_path.parent.parent / archive_relative if archive_relative else None
        if not target.exists() and archive and archive.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            fd, temporary = tempfile.mkstemp(prefix=target.name + ".", dir=target.parent)
            try:
                with os.fdopen(fd, "wb") as stream, gzip.open(archive, "rb") as compressed:
                    shutil.copyfileobj(compressed, stream)
                if sha256(temporary) != source["sha256"]:
                    raise ValueError(f"Archived input checksum mismatch: {archive}")
                os.replace(temporary, target)
            finally:
                Path(temporary).unlink(missing_ok=True)
        fetch_verified(source["url"], target, source["sha256"])
    return manifest
