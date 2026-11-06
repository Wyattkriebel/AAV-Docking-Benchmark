import pytest

from aav_benchmark.archive import bundle, restore


def test_bundle_restores_exact_bytes_and_rejects_corruption(tmp_path):
    prepared, output, destination = (tmp_path / n for n in ("prepared", "output", "restored"))
    prepared.mkdir()
    (prepared / "reference.json").write_bytes(b'{"reference":"raw bytes"}\n')
    entries = bundle(prepared, output)
    restore(output, destination)
    assert (destination / "reference.json").read_bytes() == (prepared / "reference.json").read_bytes()
    (destination / "reference.json").write_text("changed")
    with pytest.raises(ValueError, match="overwrite"):
        restore(output, destination)
    archive = output / "inputs" / entries[0]["archive"]
    archive.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Corrupt"):
        restore(output, tmp_path / "new")
