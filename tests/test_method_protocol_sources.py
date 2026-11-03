"""Actual archive integrity/source aliases; no native execution or pose metrics."""
import json
from pathlib import Path
import subprocess
from unittest.mock import patch

import pytest

from aav_benchmark import method_protocol_sources as sources


@pytest.fixture(scope="module")
def public_sources(tmp_path_factory):
    root = tmp_path_factory.mktemp("method-public-sources")
    with patch.object(subprocess, "run", side_effect=AssertionError("process call")), \
         patch.object(subprocess, "Popen", side_effect=AssertionError("process call")):
        receipt = sources.prepare_source_fixture(root, groups=("prepared", "haddock", "development"))
    return root, receipt


def test_exact_archives_and_frozen_code_without_live_trees(public_sources):
    root, receipt = public_sources
    assert receipt["hydration"]["verified_archive_members"] == 2239
    assert receipt["copied_public_files"] == 735
    assert receipt["accuracy_computed"] is receipt["native_executables_used"] is False
    assert not (root / "tools").exists() and not (root / "data/processed").exists()
    assert not (root / "results/structural_pose_v3/paired_execution_progress.json").exists()
    assert not (root / "results/blind_pose/execution_progress.json").exists()
    assert not (root / "results/structural_pose_v3/common_pose_union/scores_manifest.json").exists()
    assert not list(root.glob("**/*supplemental*"))
    assert not list(root.glob("**/*worksharing*"))
    assert len(list((root / "results/haddock_paired_cohort").glob("*/*/*/s*/run"))) == 1
    assert sources.verify_publication_amendment()["historical_archived_test_pins_verified"] == 4


def test_source_tamper_is_fatal_and_restored(public_sources):
    root, _ = public_sources
    path = root / "aav_benchmark/prospective_vina_carb.py"
    original = path.read_bytes()
    try:
        path.write_bytes(original + b"\n# changed\n")
        with pytest.raises(ValueError, match="pin changed"):
            sources.verify_source_distribution(root)
    finally:
        path.write_bytes(original)


def test_aliases_reject_symlink_then_restore_only_twelve_sources(public_sources, tmp_path):
    root, _ = public_sources
    processed = root / "data/processed"
    processed.symlink_to(tmp_path, target_is_directory=True)
    try:
        with pytest.raises(ValueError, match="symlink"):
            sources.restore_preflight_receptor_aliases(root)
        assert list(tmp_path.iterdir()) == []
    finally:
        processed.unlink()
    with patch.object(subprocess, "run", side_effect=AssertionError("native call")):
        result = sources.restore_preflight_receptor_aliases(root)
    assert result["restored_historical_receptor_aliases"] == 12
    assert result["copied_from_ignored_live_tree"] is False
    assert len(list(processed.rglob("*.*"))) == 12
    descriptor = json.loads((root / sources.DISTRIBUTION).read_text())
    for row in descriptor["historical_preflight_receptor_aliases"]:
        assert sources._digest(root / row["destination"]) == row["sha256"]


def test_fixture_rejects_existing_or_git_roots_and_unknown_groups(tmp_path):
    occupied = tmp_path / "occupied"; occupied.mkdir(); (occupied / "keep").write_text("keep")
    with pytest.raises(ValueError, match="empty"):
        sources.prepare_source_fixture(occupied)
    repository = tmp_path / "repository"; repository.mkdir(); (repository / ".git").mkdir()
    with pytest.raises(ValueError, match="Git"):
        sources.prepare_source_fixture(repository / "child")
    bare = tmp_path / "bare"; bare.mkdir()
    (bare / "HEAD").write_text("ref: refs/heads/main\n")
    (bare / "config").write_text("[core]\n bare=true\n")
    (bare / "objects").mkdir()
    with pytest.raises(ValueError, match="Git"):
        sources.prepare_source_fixture(bare / "child")
    with pytest.raises(ValueError, match="Unknown"):
        sources.prepare_source_fixture(tmp_path / "unused", groups=("typo",))
    assert (occupied / "keep").read_text() == "keep"
    assert not (tmp_path / "unused").exists()


def test_public_source_symlink_ancestors_are_rejected(tmp_path):
    directory = tmp_path / "real"; directory.mkdir()
    source = directory / "source"; source.write_bytes(b"source bytes")
    (tmp_path / "alias").symlink_to(directory, target_is_directory=True)
    with pytest.raises(ValueError, match="regular public source"):
        sources._checked(tmp_path, {"path": "alias/source", "sha256": sources._digest(source),
                                    "bytes": source.stat().st_size})
