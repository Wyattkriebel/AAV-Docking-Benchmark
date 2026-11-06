import hashlib
import errno
from unittest.mock import patch

import pytest

from aav_benchmark.provenance import fetch_verified, sha256, write_json, write_text_atomic


def test_existing_corrupt_source_is_preserved_and_rejected(tmp_path):
    path = tmp_path / "source.cif"
    path.write_bytes(b"\r\n")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        fetch_verified("https://example.org/reference", path, "a" * 64)
    assert path.read_bytes() == b"\r\n"


def test_verified_existing_source_needs_no_network(tmp_path):
    path = tmp_path / "source.cif"
    path.write_bytes(b"reference bytes")
    digest = hashlib.sha256(b"reference bytes").hexdigest()
    assert sha256(path) == digest
    assert fetch_verified("https://example.org/reference", path, digest) == path


def test_json_write_rejects_nan_without_replacing_prior_manifest(tmp_path):
    path = tmp_path / "run.json"
    write_json(path, {"status": "ok"})
    before = path.read_bytes()
    with pytest.raises(ValueError):
        write_json(path, {"rmsd": float("nan")})
    assert path.read_bytes() == before


def test_failed_report_replacement_preserves_original_and_removes_temporary(tmp_path):
    path = tmp_path / "report.md"
    path.write_text("prior complete report\n")
    with patch("aav_benchmark.provenance.os.replace", side_effect=OSError(errno.ENOSPC, "Disk full")):
        with pytest.raises(OSError):
            write_text_atomic(path, "new report\n")
    assert path.read_text() == "prior complete report\n"
    assert list(tmp_path.iterdir()) == [path]
