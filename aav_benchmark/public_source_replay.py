"""Replay the fixed twelve references without redistributing uncleared papers.

The original ledger, 174-record manifest and numerical verifier are immutable.
Exactly eleven publication/context snapshots may be absent in the public tree.
No missing-byte waiver is available to a numerical source consumer. Present
optional snapshots still have to match their original hashes.
"""
from __future__ import annotations

from contextvars import ContextVar
import hashlib
import inspect
import json
from pathlib import Path
from threading import RLock
from unittest.mock import patch

from . import heldout_reference as reference


PUBLIC_MANIFEST = Path("benchmarks/source_public_distribution.json")
SOURCE_MANIFEST = Path("data/archive/heldout_reference/source_manifest.json")
SOURCE_MANIFEST_SHA256 = "c18d983eb6b970664ad71a34c410847827c15a28db9f007936b430270ccd3e57"
LEDGER_SHA256 = "5e1142c49191ee9881efc1e918630abc0f57cffb5f2fb0041350cd6df0ac7c4e"
HELPER_SHA256 = "fa4ea12beb98d6eaee70fdb5387f8bcd90c03d6489cca008f1d4dc47f2d83e28"
PAPER_DIRECTORY = "data/archive/heldout_reference/papers/"
OMITTED_PUBLICATION_NAMES = frozenset({
    "1I3H.primary_abstract.html", "2G7C.primary_abstract.html", "2GUP.deposition_page.html",
    "2RDK.PMC2327273.html", "2WNB.publisher_response.html", "2WNB.supplement.pdf",
    "3A2S.PMC2872391.html", "3ZWE.supplement.pdf.response.html",
    "5HZA.PMC4836343.html", "5HZA.primary_density_figure.jpg",
    "8I5U.supplement.pdf.response.html",
})
RETAINED_LICENSED_NAMES = frozenset({
    "2WAB.PMC2661963.html", "3ZWE.PMC3281725.html", "4R9F.PMC4263893.html",
    "5V6F.PMC5825169.html", "8I5U.PMC10167657.html",
    "3ZWE.primary_density_figure.jpg", "4R9F.primary_density_figure.jpg",
    "5V6F.primary_density_figure.jpg",
})
_ACTIVE = ContextVar("public_source_replay_active", default=False)
_LOCK = RLock()


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text())


def validate_public_manifest(root):
    """Bind the finite omission allowlist to the unchanged original evidence."""
    root = Path(root).resolve()
    manifest = _load(root / PUBLIC_MANIFEST)
    _require(manifest.get("schema_version") == 1
             and manifest.get("distribution_id") == "fixed12_source_public_distribution_v1"
             and manifest.get("frozen_local_sources_unchanged") is True
             and manifest.get("download_required_for_numerical_replay") is False,
             "Unknown public-source distribution policy")
    bindings = (
        ("frozen_source_manifest", str(SOURCE_MANIFEST), SOURCE_MANIFEST_SHA256),
        ("frozen_reference_ledger", "benchmarks/heldout_reference_ledger.json", LEDGER_SHA256),
        ("frozen_replay_helper", "aav_benchmark/heldout_reference.py", HELPER_SHA256),
    )
    for key, path, expected in bindings:
        _require(manifest.get(key) == path and manifest.get(key + "_sha256") == expected
                 and _sha(root / path) == expected, "Frozen public-source binding changed: " + key)
    # The imported helper also must be the pinned implementation, even if the
    # caller supplies a separate archive root while running an installed package.
    _require(_sha(reference.__file__) == HELPER_SHA256, "Imported frozen reference helper changed")
    original = _load(root / SOURCE_MANIFEST)
    rows = original["archived_files"]
    pins = {row["archive"]: row for row in rows}
    _require(len(rows) == len(pins) == manifest.get("frozen_manifest_records") == 174,
             "Original source manifest denominator changed")
    publications = {row["archive"]: row for row in original["primary_paper_retrievals"]
                    if "archive" in row}
    omissions = manifest.get("metadata_only_publication_snapshots", [])
    expected = {PAPER_DIRECTORY + name for name in OMITTED_PUBLICATION_NAMES}
    _require(len(omissions) == len(expected) == 11
             and {row.get("archive") for row in omissions} == expected,
             "Public omission allowlist changed; numerical sources cannot be omitted")
    for row in omissions:
        path = row["archive"]
        _require(all(row.get(key) == pins[path][key] for key in ("sha256", "bytes"))
                 and row.get("source_url") == publications[path]["url"]
                 and isinstance(row.get("reason"), str) and bool(row["reason"].strip())
                 and row.get("role", "").startswith("Noncomputational publication/context snapshot;"),
                 "Public publication omission differs from original retrieval pin")
    _require(manifest.get("retained_archive_records") == len(pins) - len(omissions) == 163,
             "Public retained-source denominator changed")
    licensed = manifest.get("retained_licensed_publications", [])
    _require(len(licensed) == len(RETAINED_LICENSED_NAMES)
             and {row.get("archive") for row in licensed}
             == {PAPER_DIRECTORY + name for name in RETAINED_LICENSED_NAMES},
             "Retained publication attribution denominator changed")
    for row in licensed:
        path = row["archive"]
        _require(all(row.get(key) == pins[path][key] for key in ("sha256", "bytes"))
                 and row.get("source_url") == publications[path]["url"]
                 and bool(row.get("authors")) and bool(row.get("title"))
                 and bool(row.get("doi")) and bool(row.get("license_evidence"))
                 and row.get("license_url") in {
                     "https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000071",
                     "https://creativecommons.org/licenses/by/4.0/"},
                 "Retained publication attribution/license pin is incomplete")
    return manifest, pins


def verify_public_source_references(root=reference.DEFAULT_ROOT, *, density_tolerance_sigma=1e-5):
    """Run unchanged numerical replay with a narrowly scoped publication waiver.

    This does not create placeholder files. A missing optional snapshot may be
    acknowledged only at the original verifier's top-level manifest-pin loop.
    Its use by audit_candidate or another numerical consumer is an error. The
    hook is restored even after failure and does not waive checks in another
    thread/context. Every available optional snapshot is verified strictly.
    """
    root = Path(root).resolve()
    manifest, pins = validate_public_manifest(root)
    omitted = {row["archive"]: row for row in manifest["metadata_only_publication_snapshots"]}
    skipped, optional_present = set(), set()

    def public_pin(source_root, source):
        path = source.get("archive", source.get("archive_path"))
        if not _ACTIVE.get() or path not in omitted:
            return original_pin(source_root, source)
        caller = inspect.currentframe().f_back
        top_level_manifest_check = caller.f_code is reference.verify_heldout_references.__code__
        del caller
        _require(top_level_manifest_check and source == pins[path]
                 and Path(source_root).resolve() == root,
                 "Publication omission requested by a numerical source consumer")
        target = root / path
        if target.exists() or target.is_symlink():
            result = original_pin(source_root, source)
            optional_present.add(path)
            return result
        skipped.add(path)
        return target  # Never opened: original manifest loop discards this path.

    with _LOCK:
        original_pin = reference._pin
        token = _ACTIVE.set(True)
        try:
            with patch.object(reference, "_pin", public_pin):
                numerical = reference.verify_heldout_references(
                    root, density_tolerance_sigma=density_tolerance_sigma)
        finally:
            _ACTIVE.reset(token)
    _require(skipped | optional_present == set(omitted),
             "Frozen verifier did not visit every declared publication snapshot")
    _require(numerical.pop("checked_archive_count") == 174,
             "Frozen numerical helper manifest denominator changed")
    numerical.update(checked_archive_count=163, original_manifest_record_count=174,
                     metadata_only_publication_count=11,
                     optional_local_publication_byte_checks=len(optional_present))
    return dict(status="verified", distribution_id=manifest["distribution_id"],
                source_only=True, prediction_outcomes_consulted=False,
                public_manifest_sha256=_sha(root / PUBLIC_MANIFEST),
                retained_archive_records_verified=163, metadata_only_publication_count=11,
                optional_publication_snapshots_present=len(optional_present),
                optional_publication_snapshots_absent=len(skipped),
                all_consumed_numerical_sources_verified=True,
                publication_judgments_independently_reproven=False,
                numerical_reference_replay=numerical,
                omitted_publication_snapshot_paths=sorted(skipped),
                scope=manifest["replay_scope"])
