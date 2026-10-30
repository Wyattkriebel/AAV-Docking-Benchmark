"""Portable, checksummed preparation bundles stored beside raw results."""

import gzip
import json
from pathlib import Path
import shutil
import tempfile

from .provenance import sha256, write_json


def bundle(prepared_dir, output_dir):
    prepared_dir, output_dir = Path(prepared_dir), Path(output_dir)
    target = output_dir / "inputs"
    target.mkdir(parents=True, exist_ok=True)
    entries = []
    for source in sorted(prepared_dir.iterdir()):
        if not source.is_file():
            continue
        compressed = target / (source.name + ".gz")
        with compressed.open("wb") as stream:
            with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
                with source.open("rb") as original:
                    shutil.copyfileobj(original, zipped)
        entries.append({"filename": source.name, "sha256": sha256(source),
                        "archive": compressed.name, "archive_sha256": sha256(compressed),
                        "original_bytes": source.stat().st_size})
    write_json(target / "manifest.json", {"schema_version": 1, "files": entries})
    return entries


def restore(output_dir, destination):
    target, destination = Path(output_dir) / "inputs", Path(destination)
    entries = json.loads((target / "manifest.json").read_text())["files"]
    destination.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        for key in ("filename", "archive"):
            if Path(entry[key]).name != entry[key]:
                raise ValueError("Archive filenames must not contain directories")
        source, output = target / entry["archive"], destination / entry["filename"]
        if sha256(source) != entry["archive_sha256"]:
            raise ValueError(f"Corrupt preparation archive: {source}")
        if output.exists():
            if sha256(output) != entry["sha256"]:
                raise ValueError(f"Refusing to overwrite changed input: {output}")
            continue
        with tempfile.TemporaryDirectory(dir=destination) as temporary:
            restored = Path(temporary) / "restored"
            with gzip.open(source, "rb") as compressed, restored.open("wb") as stream:
                shutil.copyfileobj(compressed, stream)
            if sha256(restored) != entry["sha256"]:
                raise ValueError(f"Restored input checksum mismatch: {output}")
            restored.replace(output)
    return entries
