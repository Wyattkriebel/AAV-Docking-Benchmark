"""Public archive-only replay controls; no binaries or native pose accuracy."""
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

import aav_benchmark.paired_method_report as paired
from aav_benchmark.heldout_pose_report import FAMILIES
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.qualified_pose import Evidence


ROOT = Path(__file__).resolve().parents[1]


def _copy_source_tree(destination):
    """Restore exact archived preparations; never inspect a live output tree."""
    from aav_benchmark.method_protocol_sources import prepare_source_fixture
    prepare_source_fixture(destination, source_root=ROOT, groups=("prepared",))
    distribution = json.loads((destination / "benchmarks/source_public_distribution.json").read_text())
    omissions = {record["archive"] for record in distribution["metadata_only_publication_snapshots"]}
    experiment = json.loads((destination / "benchmarks/structural_test_vina_carb.json").read_text())
    return omissions, experiment


def _evaluation_plan(root, experiment):
    def pin(path):
        return {"path": path, "sha256": sha256(root / path)}
    reference = {"source_ledger": pin("benchmarks/heldout_reference_ledger.json"),
        "reference_source_manifest": pin("data/archive/heldout_reference/source_manifest.json"),
        "family_audit": pin("benchmarks/family_overlap_audit.json"), "case_to_family": FAMILIES,
        "site_prediction": {"configuration": pin("benchmarks/structural_site_v2_experiment.json"),
            "freeze_receipt": pin("results/structural_test_v2/freeze_receipt.json"),
            "source_archive_manifest": pin("results/structural_test_v2/source_archive_manifest.json")},
        "profiles": [{"id": name, "source_archive_manifest": pin(f"results/structural_pose_v3/{profile}_sources_manifest.json")}
                     for profile, name in (("strict", "strict_CCD_ideal"), ("repair", "modeled_sidechain_CCD_ideal"))]}
    reference_path = "benchmarks/structural_pose_evaluation.json"
    write_json(root / reference_path, reference)
    reference_freeze_path = "results/structural_pose_v3/evaluation_freeze.json"
    write_json(root / reference_freeze_path, {"evaluation_plan_sha256": sha256(root / reference_path),
        "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    plan = paired.example_plan(root)
    plan["vina_carb_freeze"] = pin("results/structural_pose_v3/vina_carb_freeze_receipt.json")
    plan["vina_carb_executable_sha256"] = experiment["native_binary_sha256"]
    plan_path, freeze_path = "benchmarks/paired_evaluation.json", "results/paired_evaluation_freeze.json"
    write_json(root / plan_path, plan)
    write_json(root / freeze_path, {"evaluation_plan_sha256": sha256(root / plan_path),
        "frozen_before_test_pose_accuracy": True, "new_test_pose_accuracy_consulted": False})
    return plan_path, freeze_path, reference, plan


def test_archive_only_public_checkout_missing11_and_empty_tools_stays_strictly_pending(tmp_path, monkeypatch):
    omitted, experiment = _copy_source_tree(tmp_path)
    assert len(omitted) == 11 and all(not (tmp_path / p).exists() for p in omitted)
    assert not (tmp_path / "tools").exists() and not (tmp_path / "data/processed").exists()
    receipts = [experiment[k]["path"] for k in ("native_build_receipt", "native_preflight_receipt")]
    assert len(receipts) == 2 and all(p.startswith("results/") and (tmp_path / p).is_file() for p in receipts)
    plan_path, freeze_path, reference, plan = _evaluation_plan(tmp_path, experiment)
    def forbidden(*args, **kwargs):
        pytest.fail("Archive-only pending validation entered actual native accuracy")
    monkeypatch.setattr(paired, "_source_context", forbidden)
    monkeypatch.setattr(paired, "score_pool", forbidden)
    monkeypatch.setattr(paired, "site_summary", forbidden)
    report = paired.generate_paired_report(tmp_path, plan_path, freeze_path, allow_pending=True)
    assert report["status"] == "pending_incomplete_workflow" and not report["metrics_eligible"]
    assert report["counts"]["eligible_method_outcomes"] == 216 and report["counts"]["terminal_method_outcomes"] == 0
    assert report["reference_validation"]["retained_archive_records_verified"] == 163
    assert report["reference_validation"]["metadata_only_publication_count"] == 11
    assert not set(report["source_sha256_inventory"]) & omitted
    assert not any(path.startswith(("tools/", "data/processed/")) for path in report["source_sha256_inventory"])
    write_json(tmp_path / "public_dependency_inventory.json", report["source_sha256_inventory"])
    # Import and execute the copied package itself, rather than relying on an
    # installed source module elsewhere in the developer workspace.
    script = """
import json
from pathlib import Path
import aav_benchmark.paired_method_report as report
assert Path(report.__file__).resolve().is_relative_to(Path.cwd())
def forbidden(*args, **kwargs):
    raise AssertionError('Native accuracy reached in archive-only pending replay')
report._source_context = forbidden
report.score_pool = forbidden
report.site_summary = forbidden
result = report.generate_paired_report(Path.cwd(), 'benchmarks/paired_evaluation.json', 'results/paired_evaluation_freeze.json', allow_pending=True)
print(json.dumps({'status': result['status'], 'metrics_eligible': result['metrics_eligible'], 'terminal': result['counts']['terminal_method_outcomes'], 'tools_present': Path('tools').exists()}))
"""
    replay = json.loads(subprocess.check_output([sys.executable, "-c", script], cwd=tmp_path, text=True))
    assert replay == {"status": "pending_incomplete_workflow", "metrics_eligible": False, "terminal": 0, "tools_present": False}
    # Availability metadata is normalized independently of metric inputs. Use a
    # synthetic availability variant of an actual verified public replay; no
    # omitted/private article bytes are required or treated as valid evidence.
    verified = paired.public_source_replay.verify_public_source_references(tmp_path)
    available = {**verified, "optional_publication_snapshots_present": 1,
                 "optional_publication_snapshots_absent": 10}
    with monkeypatch.context() as metadata:
        metadata.setattr(paired.public_source_replay, "verify_public_source_references",
                         lambda *args, **kwargs: available)
        cases, validation = paired.check_public_reference_sources(Evidence(tmp_path), reference["source_ledger"],
            reference["reference_source_manifest"], plan["public_source_distribution"])
    assert set(cases) == set(paired.CASE_IDS) and validation == report["reference_validation"]
    # Actual optional bytes, when present, still require their original hash.
    one = sorted(omitted)[0]
    (tmp_path / one).parent.mkdir(parents=True, exist_ok=True)
    (tmp_path / one).write_bytes(b"synthetic invalid optional publication bytes\n")
    with pytest.raises(ValueError, match="checksum|changed|mismatch"):
        paired.check_public_reference_sources(Evidence(tmp_path), reference["source_ledger"],
            reference["reference_source_manifest"], plan["public_source_distribution"])
    (tmp_path / one).unlink()
    # The finite paper waiver never excuses a numerical archive, even though
    # this pending execution has no pose/native-accuracy calculation at all.
    numeric = tmp_path / "data/archive/heldout_candidates/2RDK.cif.gz"
    numeric.write_bytes(numeric.read_bytes() + b"numeric tamper")
    with pytest.raises(ValueError, match="checksum|changed|mismatch"):
        paired.check_public_reference_sources(Evidence(tmp_path), reference["source_ledger"],
            reference["reference_source_manifest"], plan["public_source_distribution"])
    shutil.copyfile(ROOT / numeric.relative_to(tmp_path), numeric)
    assert all(not (tmp_path / path).exists() for path in omitted)


def test_original_vc_six_column_table_is_not_modern_four_column_format():
    poses = [{"rank": 1, "method_score": -5.0}, {"rank": 2, "method_score": -4.0}]
    raw = (b"mode | affinity | chi | affinity | dist from best mode\n"
           b"1 -5.0 0.2 -5.2 0.000 0.000\n2 -4.0 0.3 -4.3 1.234 2.345\n")
    assert paired.verify_carb_score_table(raw, poses)["raw_score_tolerance_kcal_mol"] == .002
    with pytest.raises(ValueError, match="six-column"):
        paired.verify_carb_score_table(b"1 -5.0 0.000 0.000\n2 -4.0 1.234 2.345\n", poses)
    with pytest.raises(ValueError, match="disagrees"):
        paired.verify_carb_score_table(raw.replace(b"1 -5.0", b"1 -5.1"), poses)
