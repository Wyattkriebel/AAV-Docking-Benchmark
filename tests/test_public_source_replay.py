"""Public packaging may omit named papers, never computational evidence."""
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil

import pytest

from aav_benchmark import heldout_reference as reference
from aav_benchmark import atomic_reference as atomic
from aav_benchmark.public_source_replay import (
    OMITTED_PUBLICATION_NAMES, PAPER_DIRECTORY, PUBLIC_MANIFEST,
    source_bond_distance_replay, validate_public_manifest, verify_public_source_references,
)


REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def public_root(tmp_path):
    root = tmp_path / "public"
    manifest = json.loads((REPO / "data/archive/heldout_reference/source_manifest.json").read_text())
    omitted = {PAPER_DIRECTORY + name for name in OMITTED_PUBLICATION_NAMES}
    for row in manifest["archived_files"]:
        if row["archive"] not in omitted:
            target = root / row["archive"]
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(REPO / row["archive"])
    for name in ("benchmarks/heldout_candidate_inventory.json", "benchmarks/heldout_reference_ledger.json",
                 str(PUBLIC_MANIFEST), "data/archive/heldout_reference/source_manifest.json",
                 "data/archive/atomic_reference/qualification_policy.json", "aav_benchmark/heldout_reference.py"):
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_symlink():
            target.unlink()
        shutil.copyfile(REPO / name, target)
    return root


def edit_json(path, operation):
    data = json.loads(path.read_text())
    operation(data)
    if path.is_symlink():
        path.unlink()
    path.write_text(json.dumps(data))


def test_all_eleven_absent_publication_snapshots_preserve_exact_numerical_replay(public_root):
    result = verify_public_source_references(public_root)
    assert result["retained_archive_records_verified"] == 163
    assert result["metadata_only_publication_count"] == result["optional_publication_snapshots_absent"] == 11
    assert result["optional_publication_snapshots_present"] == 0
    assert result["all_consumed_numerical_sources_verified"] is True
    assert result["publication_judgments_independently_reproven"] is False
    numerical = result["numerical_reference_replay"]
    assert numerical["checked_archive_count"] == 163
    assert numerical["original_manifest_record_count"] == 174
    assert numerical["case_count"] == 12
    assert numerical["primary_scored_atom_count"] == 298
    assert numerical["qualified_assembly_copy_count"] == 16
    assert numerical["maximum_absolute_density_difference_sigma"] <= 1e-5
    roundoff = result["derived_bond_distance_roundoff"]
    assert roundoff["tolerance_A"] == 1e-12
    assert roundoff["compared_distance_fields"] > 0
    assert roundoff["maximum_absolute_difference_A"] <= 1e-12


def test_original_strict_verifier_still_rejects_missing_publication_bytes(public_root):
    with pytest.raises(reference.AtomicReferenceError, match="checksum mismatch.*papers"):
        reference.verify_heldout_references(public_root)


def test_numerical_map_cannot_join_the_publication_allowlist(public_root):
    path = public_root / PUBLIC_MANIFEST
    original = json.loads((public_root / "data/archive/heldout_reference/source_manifest.json").read_text())
    map_pin = next(row for row in original["archived_files"] if "5v6f_validation_2fo-fc" in row["archive"])
    edit_json(path, lambda data: data["metadata_only_publication_snapshots"].append(deepcopy(map_pin)))
    with pytest.raises(ValueError, match="omission allowlist"):
        verify_public_source_references(public_root)


def test_missing_map_still_fails_and_restores_the_original_pin_function(public_root):
    (public_root / "data/archive/heldout_reference/5v6f_validation_2fo-fc_map_coef.cif.gz").unlink()
    original_pin = reference._pin
    with pytest.raises(reference.AtomicReferenceError, match="checksum mismatch.*5v6f_validation_2fo-fc"):
        verify_public_source_references(public_root)
    assert reference._pin is original_pin


@pytest.mark.parametrize("name", ["5V6F.PMC5825169.html", "2RDK.PMC2327273.html"])
def test_present_retained_or_optional_article_must_match_its_original_bytes(public_root, name):
    path = public_root / PAPER_DIRECTORY / name
    if path.is_symlink():
        path.unlink()
    path.write_bytes(b"synthetic tampered publication bytes\n")
    with pytest.raises(reference.AtomicReferenceError, match="checksum mismatch"):
        verify_public_source_references(public_root)


def test_omission_hash_cannot_be_rebound_in_public_manifest(public_root):
    edit_json(public_root / PUBLIC_MANIFEST,
              lambda data: data["metadata_only_publication_snapshots"][0].update(sha256="0" * 64))
    with pytest.raises(ValueError, match="original retrieval pin"):
        validate_public_manifest(public_root)


def test_original_manifest_cannot_be_changed_to_drop_a_computational_pin(public_root):
    edit_json(public_root / "data/archive/heldout_reference/source_manifest.json",
              lambda data: data["archived_files"].pop())
    with pytest.raises(ValueError, match="Frozen public-source binding"):
        verify_public_source_references(public_root)


def test_a_numerical_consumer_cannot_receive_missing_publication_bytes(public_root, monkeypatch):
    manifest = json.loads((public_root / "data/archive/heldout_reference/source_manifest.json").read_text())
    paper = next(row for row in manifest["archived_files"]
                 if row["archive"] == PAPER_DIRECTORY + "2RDK.PMC2327273.html")
    def numerical_consumer(*args, **kwargs):
        reference._pin(public_root, paper)
    monkeypatch.setattr(reference, "audit_all", numerical_consumer)
    original_pin = reference._pin
    with pytest.raises(ValueError, match="numerical source consumer"):
        verify_public_source_references(public_root)
    assert reference._pin is original_pin


@pytest.mark.parametrize("label", [
    "heldout primary source diagnostics.cases[0].native_bond_graph[8].distance_A",
    "1I3H native bond graph[8].distance_A",
])
def test_observed_one_ulp_bond_distance_roundoff_is_explicitly_reported(label):
    actual, expected = 1.4351466823986994, 1.4351466823986996
    density_deltas = []
    with source_bond_distance_replay() as report:
        atomic._compare_replay(actual, expected, label, 1e-5, density_deltas)
    assert report["maximum_absolute_difference_A"] == abs(actual - expected)
    assert report["compared_distance_fields"] == report["nonzero_difference_fields"] == 1
    assert density_deltas == []


def test_native_bond_distance_profile_reaches_only_the_named_nested_scalar():
    actual = {"cases": [{"native_bond_graph": [{"atom1": "B:1:MAN:C1", "distance_A": 1.4351466823986994}]}]}
    expected = deepcopy(actual)
    expected["cases"][0]["native_bond_graph"][0]["distance_A"] = 1.4351466823986996
    with source_bond_distance_replay() as report:
        reference._compare_replay(actual, expected, "heldout primary source diagnostics", 1e-5, [])
    assert report["nonzero_difference_fields"] == 1


@pytest.mark.parametrize("label", [
    "heldout primary source diagnostics.cases[0].native_atom_map[8].xyz_A[0]",
    "heldout primary source diagnostics.cases[0].native_bond_graph[8].struct_conn_distance_A",
    "heldout primary source diagnostics.cases[0].other_bond_graph[8].distance_A",
    "heldout primary source diagnostics.cases[0].native_bond_graph_extra[8].distance_A",
    "heldout primary source diagnostics.cases[0].native_bond_graph[8].distance_A_extra",
])
def test_raw_coordinates_deposited_distances_and_other_fields_remain_exact(label):
    with source_bond_distance_replay(), pytest.raises(atomic.AtomicReferenceError, match="source replay differs"):
        atomic._compare_replay(1.4351466823986994, 1.4351466823986996, label, 100, [])


@pytest.mark.parametrize("field,actual,expected", [
    ("atom1", "B:1:MAN:C1", "B:1:MAN:C2"),
    ("order", "SING", "DOUB"),
    ("source_gate_pass", False, True),
    ("sha256", "0" * 64, "1" * 64),
])
def test_graph_identities_decisions_and_hashes_cannot_use_bond_tolerance(field, actual, expected):
    with source_bond_distance_replay(), pytest.raises(atomic.AtomicReferenceError, match="source replay differs"):
        atomic._compare_replay({field: actual}, {field: expected},
                              "heldout primary source diagnostics.cases[0].native_bond_graph[8]", 1e-5, [])


def test_larger_bond_difference_and_zero_tolerance_fail_and_restore_original_comparator():
    original = atomic._compare_replay
    label = "heldout primary source diagnostics.cases[0].native_bond_graph[8].distance_A"
    for tolerance, actual, expected in [(1e-12, 1.4 + 2e-12, 1.4),
                                        (0, 1.4351466823986994, 1.4351466823986996)]:
        with pytest.raises(atomic.AtomicReferenceError, match="derived bond distance differs"):
            with source_bond_distance_replay(tolerance):
                atomic._compare_replay(actual, expected, label, 1e-5, [])
        assert atomic._compare_replay is original
    with pytest.raises(atomic.AtomicReferenceError, match="source replay differs"):
        atomic._compare_replay(1.4351466823986994, 1.4351466823986996, label, 1e-5, [])


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True, "1e-12"])
def test_bond_tolerance_requires_an_explicit_finite_nonnegative_number(value):
    with pytest.raises(ValueError, match="finite and nonnegative"):
        with source_bond_distance_replay(value):
            pass


def test_bond_roundoff_context_does_not_relax_another_thread():
    label = "heldout primary source diagnostics.cases[0].native_bond_graph[8].distance_A"
    with source_bond_distance_replay(), ThreadPoolExecutor(max_workers=1) as worker:
        future = worker.submit(atomic._compare_replay, 1.4351466823986994, 1.4351466823986996,
                               label, 1e-5, [])
        with pytest.raises(atomic.AtomicReferenceError, match="source replay differs"):
            future.result()


def test_existing_density_threshold_classifications_remain_exact():
    key = "two_mFo_DFc_max_within_0p5A_sigma"
    with source_bond_distance_replay(), pytest.raises(atomic.AtomicReferenceError, match="coverage_classification"):
        atomic._compare_replay({key: .9999999999999999}, {key: 1.0}, "native density", 1e-5, [])
