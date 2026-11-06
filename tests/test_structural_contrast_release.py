"""Synthetic release provenance and stdlib-only contrast-reader portability.

No live report outputs, poses, atomic references, or native executables are read.
The complete report headers below deliberately contain no pose metrics.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from aav_benchmark import structural_contrast_report as report


ROOT = Path(__file__).resolve().parents[1]
PLAN = "benchmarks/synthetic_structural_contrast.json"
FREEZE = "results/synthetic_structural_contrast_freeze.json"
SCRIPT = "scripts/regenerate_structural_contrast_report.py"


def _bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _write(root, relative, value):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value if isinstance(value, bytes) else _bytes(value))
    return path


def _read(root, relative):
    return json.loads((root / relative).read_bytes())


def _freeze(root, plan):
    _write(root, PLAN, plan)
    _write(root, FREEZE, {
        "aggregation_plan_sha256": _sha((root / PLAN).read_bytes()),
        "frozen_before_any_structural_accuracy": True,
        "new_structural_accuracy_consulted": False,
    })


@pytest.fixture
def released(tmp_path):
    """Copy only source and frozen metadata, then create a fake release packet."""
    for relative in (*report.CODE, "aav_benchmark/__init__.py", SCRIPT,
                     "PROTOCOL.md", "benchmarks/structural_pose_reference_evaluation.json",
                     report.GLOBAL, report.AUTHORIZATION):
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    plan = report.example_plan(tmp_path)
    plan["frozen_before_any_structural_accuracy"] = True
    _freeze(tmp_path, plan)
    contract = _read(tmp_path, report.GLOBAL)
    readers = []
    for role in report.ROLES:
        counts = deepcopy(report.COUNTS[role])
        readers.append({
            "reader": role, "status": "locally_ready_accuracy_blocked",
            "accuracy_computed": False, "source_context_not_entered": True,
            "counts": counts,
        })
        if role == "fixed_union":
            readers[-1]["pure_validation_only"] = True
        spec = contract["pure_union"] if role == "fixed_union" else contract["readers"][role]
        _write(tmp_path, report.DIRECTORY + "/" + role + ".json", {
            "schema_version": 1, "status": "complete", "metrics_eligible": True,
            "counts": counts, "metric_code_sha256": spec["metric_sources"],
            "synthetic_test_header_only": True,
        })
    _write(tmp_path, report.GENERATION, {
        "schema_version": 1, "status": "all_four_reports_generated",
        "all_four_ready": True, "accuracy_computed": True,
        "real_reports_generated": True, "output": report.DIRECTORY,
        "contract": plan["global_contract"], "probe_freeze": plan["global_authorization"],
        "readers": readers,
    })
    assert not (tmp_path / "tools").exists()
    assert not (tmp_path / "data/processed").exists()
    return tmp_path


def _bind(root):
    return report.bind_released_reports(root, PLAN, FREEZE)


@pytest.mark.parametrize("field", ["hypothesis_labels", "endpoints", "new_contrasts"])
def test_example_plan_cannot_mutate_frozen_choice_constants(released, field):
    plan = report.example_plan(released)
    baseline = report.example_plan(released)[field]
    if field == "hypothesis_labels": plan[field]["H3"] = "preparation"
    elif field == "endpoints": plan[field]["thresholds_A"].append(4)
    else: plan[field]["main_H1_guidance"] = 0
    assert report.example_plan(released)[field] == baseline
    _freeze(released, plan)
    with pytest.raises(ValueError): report.validate_plan(released, PLAN, FREEZE)


@pytest.mark.parametrize("purity", [None, False])
def test_fixed_union_requires_the_original_pure_validation_marker(released, purity):
    receipt = _read(released, report.GENERATION)
    row = next(r for r in receipt["readers"] if r["reader"] == "fixed_union")
    if purity is None: row.pop("pure_validation_only")
    else: row["pure_validation_only"] = purity
    _write(released, report.GENERATION, receipt)
    with pytest.raises(ValueError, match="pure validation"): _bind(released)


@pytest.mark.parametrize("role", [*report.ROLES, "generation_receipt"])
def test_bound_report_and_generation_byte_counts_cannot_be_removed(released, monkeypatch, role):
    _bind(released); _no_aggregation(monkeypatch)
    def mutate(binding):
        pin = binding["generation_receipt"] if role == "generation_receipt" else binding["reports"][role]
        pin.pop("bytes")
    digest = _rewrite_binding(released, mutate)
    with pytest.raises(ValueError, match="byte"):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=digest)


def test_extra_file_after_binding_cannot_enter_aggregation(released, monkeypatch):
    pin = _bind(released); _no_aggregation(monkeypatch)
    _write(released, report.DIRECTORY + "/extra.json", {})
    with pytest.raises(ValueError, match="four-file packet"):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=pin["sha256"])


@pytest.mark.parametrize("empty_pin", [{}, None, False, []])
def test_empty_generation_pin_cannot_trigger_posthoc_fallback(released, monkeypatch, empty_pin):
    _bind(released); _no_aggregation(monkeypatch)
    digest = _rewrite_binding(released, lambda binding: binding.update(generation_receipt=empty_pin))
    with pytest.raises(ValueError, match="receipt pin/path"):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=digest)


def _rewrite_binding(root, mutation):
    binding = _read(root, report.BINDING)
    mutation(binding)
    _write(root, report.BINDING, binding)
    return _sha((root / report.BINDING).read_bytes())


def _no_aggregation(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Invalid release provenance reached outcome aggregation")
    monkeypatch.setattr(report, "_conditions", forbidden)
    monkeypatch.setattr(report, "paired_family_bootstrap", forbidden)
    monkeypatch.setattr(report, "aggregate_runs", forbidden)


def test_four_report_binding_has_exact_roles_bytes_hashes_and_honest_timing(released):
    before = {role: (released / report.DIRECTORY / (role + ".json")).read_bytes()
              for role in report.ROLES}
    pin = _bind(released)
    raw = (released / report.BINDING).read_bytes()
    assert pin == {"path": report.BINDING, "sha256": _sha(raw), "bytes": len(raw)}
    binding = json.loads(raw)
    assert binding["status"] == "four_released_reports_bound"
    assert set(binding["reports"]) == set(report.ROLES)
    assert binding["new_pose_metrics_computed"] is False
    assert binding["native_execution"] is False
    assert binding["hash_binding_timing"] == (
        "post_generation; the frozen generation receipt itself does not carry report-file hashes")
    assert "reports" not in _read(released, report.GENERATION)
    for role, expected_raw in before.items():
        assert binding["reports"][role] == {
            "path": report.DIRECTORY + "/" + role + ".json",
            "sha256": _sha(expected_raw), "bytes": len(expected_raw),
        }
        assert (released / binding["reports"][role]["path"]).read_bytes() == expected_raw
    generation_raw = (released / report.GENERATION).read_bytes()
    assert binding["generation_receipt"] == {
        "path": report.GENERATION, "sha256": _sha(generation_raw), "bytes": len(generation_raw)}
    # Exact replay may reuse the binding, and never rewrites the producer reports.
    assert _bind(released) == pin
    assert (released / report.BINDING).read_bytes() == raw


@pytest.mark.parametrize("mutation", [
    lambda value: value.update(status="all_four_reports_ready_accuracy_blocked"),
    lambda value: value.update(all_four_ready=False),
    lambda value: value.update(accuracy_computed=False),
    lambda value: value.update(real_reports_generated=False),
    lambda value: value.update(output="results/some_other_reports"),
    lambda value: value["contract"].update(sha256="0" * 64),
    lambda value: value["probe_freeze"].update(path="results/wrong_freeze.json"),
    lambda value: value["readers"].pop(),
    lambda value: value["readers"][0].update(reader=value["readers"][1]["reader"]),
    lambda value: value["readers"][0].update(status="complete"),
    lambda value: value["readers"][0].update(accuracy_computed=True),
    lambda value: value["readers"][0].update(source_context_not_entered=False),
    lambda value: value["readers"][0]["counts"].update(eligible_conditions=35),
    lambda value: value["readers"][0]["counts"].update(eligible_conditions=True),
], ids=("not-generated", "not-all-ready", "accuracy-not-generated", "no-real-release",
         "wrong-output", "wrong-contract", "wrong-authorization", "missing-reader",
         "duplicate-role", "reader-no-block", "reader-used-accuracy", "reader-entered-context",
         "wrong-denominator", "boolean-denominator"))
def test_binding_rejects_unsuccessful_or_mismatched_generation_receipt(released, mutation):
    value = _read(released, report.GENERATION)
    mutation(value)
    _write(released, report.GENERATION, value)
    with pytest.raises(ValueError):
        _bind(released)
    assert not (released / report.BINDING).exists()


def test_binding_requires_generation_receipt(released):
    (released / report.GENERATION).unlink()
    with pytest.raises(FileNotFoundError):
        _bind(released)
    assert not (released / report.BINDING).exists()


@pytest.mark.parametrize("change", ("missing", "extra", "renamed", "directory", "symlink"))
def test_binding_requires_exact_four_ordinary_report_paths(released, change):
    path = released / report.DIRECTORY / "main_paired.json"
    if change == "missing":
        path.unlink()
    elif change == "extra":
        _write(released, report.DIRECTORY + "/extra.json", {})
    elif change == "renamed":
        path.rename(path.with_name("other_role.json"))
    elif change == "directory":
        path.unlink()
        path.mkdir()
    else:
        saved = released / "synthetic_header.json"
        path.rename(saved)
        path.symlink_to(saved)
    with pytest.raises((ValueError, IsADirectoryError)):
        _bind(released)
    assert not (released / report.BINDING).exists()


@pytest.mark.parametrize("mutation", [
    lambda value: value.update(status="pending_incomplete_workflow"),
    lambda value: value.update(metrics_eligible=False),
    lambda value: value["counts"].update(terminal_method_outcomes=215),
    lambda value: value["metric_code_sha256"].update({"aav_benchmark/statistics.py": "0" * 64}),
], ids=("pending", "ineligible", "denominator", "metric-code-provenance"))
def test_binding_rejects_invalid_report_headers(released, mutation):
    path = report.DIRECTORY + "/main_paired.json"
    value = _read(released, path)
    mutation(value)
    _write(released, path, value)
    with pytest.raises(ValueError):
        _bind(released)
    assert not (released / report.BINDING).exists()


def test_report_and_readiness_extra_counters_must_agree(released):
    receipt = _read(released, report.GENERATION)
    receipt["readers"][0]["counts"]["additional_terminal_count"] = 7
    _write(released, report.GENERATION, receipt)
    with pytest.raises(ValueError, match="counters differ"):
        _bind(released)


@pytest.mark.parametrize("layer", ("plan", "freeze", "generation", "report", "binding"))
def test_duplicate_json_keys_are_rejected_at_every_release_layer(released, layer, monkeypatch):
    if layer == "binding":
        _bind(released)
    paths = {"plan": PLAN, "freeze": FREEZE, "generation": report.GENERATION,
             "report": report.DIRECTORY + "/main_paired.json", "binding": report.BINDING}
    value = _read(released, paths[layer])
    raw = b'{"schema_version":1,"schema_version":1,' + _bytes(value).lstrip()[1:]
    _write(released, paths[layer], raw)
    if layer == "binding":
        _no_aggregation(monkeypatch)
        operation = lambda: report.generate_report(released, PLAN, FREEZE, binding_sha256=_sha(raw))
    else:
        operation = lambda: _bind(released)
    with pytest.raises(ValueError, match="Duplicate JSON"):
        operation()


@pytest.mark.parametrize("change", ("plan-unfrozen", "missing-freeze", "bad-plan-hash", "late-freeze", "accuracy-consulted"))
def test_separate_pre_accuracy_freeze_is_required(released, change):
    if change == "plan-unfrozen":
        plan = _read(released, PLAN)
        plan["frozen_before_any_structural_accuracy"] = False
        _freeze(released, plan)
    elif change == "missing-freeze":
        (released / FREEZE).unlink()
    else:
        freeze = _read(released, FREEZE)
        field, value = {
            "bad-plan-hash": ("aggregation_plan_sha256", "0" * 64),
            "late-freeze": ("frozen_before_any_structural_accuracy", False),
            "accuracy-consulted": ("new_structural_accuracy_consulted", True),
        }[change]
        freeze[field] = value
        _write(released, FREEZE, freeze)
    with pytest.raises((ValueError, FileNotFoundError)):
        _bind(released)
    assert not (released / report.BINDING).exists()


@pytest.mark.parametrize("consulted", (True, None, 0, "false"))
def test_plan_itself_must_declare_no_structural_accuracy_consulted(released, consulted):
    plan = _read(released, PLAN)
    plan["new_structural_accuracy_consulted"] = consulted
    # Refreshing the checksum cannot authorize a changed scientific declaration.
    _freeze(released, plan)
    with pytest.raises(ValueError, match="pre-accuracy declaration"):
        _bind(released)
    assert not (released / report.BINDING).exists()


@pytest.mark.parametrize("mutation", [
    lambda plan: plan["endpoints"]["atomic_strata"].reverse(),
    lambda plan: plan["endpoints"]["alignments"].pop(),
    lambda plan: plan["endpoints"]["rank_budgets"].append("top10"),
    lambda plan: plan["endpoints"].update(thresholds_A=[1, 2, 4]),
    lambda plan: plan["endpoints"].update(unfrozen_endpoint=True),
    lambda plan: plan["new_contrasts"].update(main_H1_guidance=287),
    lambda plan: plan["new_contrasts"].update(main_preparation_sensitivity=215),
    lambda plan: plan["new_contrasts"].update(haddock_preparation_sensitivity=107),
    lambda plan: plan["new_contrasts"].update(additional_contrast=1),
], ids=("strata-order", "missing-alignment", "extra-rank", "threshold", "extra-endpoint",
         "main-H1-count", "main-preparation-count", "haddock-preparation-count", "extra-contrast"))
def test_bounded_endpoint_and_contrast_declarations_cannot_be_refrozen_to_changes(released, mutation):
    plan = _read(released, PLAN)
    mutation(plan)
    _freeze(released, plan)
    with pytest.raises(ValueError, match="Bounded contrast choices"):
        _bind(released)
    assert not (released / report.BINDING).exists()


@pytest.mark.parametrize("relative", (*report.CODE, report.GLOBAL, report.AUTHORIZATION,
                                      "PROTOCOL.md", "benchmarks/structural_pose_reference_evaluation.json"))
def test_frozen_source_and_metadata_tampering_is_rejected(released, relative):
    path = released / relative
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        _bind(released)


@pytest.mark.parametrize("target", ("report", "generation", "binding"))
def test_generate_rejects_after_binding_tampering_before_aggregation(released, target, monkeypatch):
    pin = _bind(released)
    relative = {"report": report.DIRECTORY + "/main_paired.json",
                "generation": report.GENERATION, "binding": report.BINDING}[target]
    path = released / relative
    path.write_bytes(path.read_bytes() + b"\n")
    _no_aggregation(monkeypatch)
    with pytest.raises(ValueError, match="checksum"):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=pin["sha256"])


@pytest.mark.parametrize("mutation", [
    lambda value: value["reports"].pop("haddock"),
    lambda value: value["reports"].update(extra=deepcopy(value["reports"]["main_paired"])),
    lambda value: value["reports"]["main_paired"].update(path=report.DIRECTORY + "/haddock.json"),
    lambda value: value["reports"]["main_paired"].update(bytes=0),
    lambda value: value["generation_receipt"].update(path="results/wrong_receipt.json"),
    lambda value: value["global_contract"].update(sha256="0" * 64),
    lambda value: value["global_authorization"].update(sha256="0" * 64),
    lambda value: value.update(hash_binding_timing="producer-time hashes"),
], ids=("missing-role", "extra-role", "wrong-role-path", "wrong-byte-count", "wrong-receipt-path",
         "wrong-contract", "wrong-authorization", "dishonest-timing"))
def test_generate_rejects_inconsistent_rehashed_bindings(released, mutation, monkeypatch):
    _bind(released)
    digest = _rewrite_binding(released, mutation)
    _no_aggregation(monkeypatch)
    with pytest.raises(ValueError):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=digest)


@pytest.mark.parametrize("field", ("new_pose_metrics_computed", "native_execution"))
@pytest.mark.parametrize("value", (True, None, 0, "false", "missing"))
def test_rehashed_binding_must_declare_no_new_metrics_or_native_execution(released, field, value, monkeypatch):
    _bind(released)
    def mutation(binding):
        if value == "missing":
            binding.pop(field)
        else:
            binding[field] = value
    digest = _rewrite_binding(released, mutation)
    _no_aggregation(monkeypatch)
    with pytest.raises(ValueError, match="Post-generation report binding"):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=digest)


@pytest.mark.parametrize("digest", (None, "", "0" * 63, "0" * 64))
def test_generate_requires_explicit_matching_binding_hash(released, digest, monkeypatch):
    _bind(released)
    _no_aggregation(monkeypatch)
    with pytest.raises(ValueError):
        report.generate_report(released, PLAN, FREEZE, binding_sha256=digest)


def test_rebinding_refuses_to_replace_different_existing_binding(released):
    _bind(released)
    original = (released / report.BINDING).read_bytes()
    path = report.DIRECTORY + "/main_paired.json"
    value = _read(released, path)
    value["synthetic_header_revision"] = 2
    _write(released, path, value)
    with pytest.raises(ValueError, match="Refusing to replace"):
        _bind(released)
    assert (released / report.BINDING).read_bytes() == original


@pytest.mark.parametrize("relative", (PLAN, report.DIRECTORY + "/main_paired.json"))
def test_input_symlink_loops_are_safe_validation_failures(released, relative):
    path = released / relative
    path.unlink()
    path.symlink_to(path.name)
    with pytest.raises(ValueError, match="resolved safely|ordinary file"):
        _bind(released)
    assert not (released / report.BINDING).exists()


def test_fresh_output_is_idempotent_and_conflicting_output_is_preserved(released):
    relative = "results/synthetic_contrast_output.json"
    value = {"synthetic_indicator_only": True}
    report._write_new(released, relative, value)
    original = (released / relative).read_bytes()
    report._write_new(released, relative, value)
    assert (released / relative).read_bytes() == original
    with pytest.raises(ValueError, match="Refusing to replace"):
        report._write_new(released, relative, {"synthetic_indicator_only": False})
    assert (released / relative).read_bytes() == original


@pytest.mark.parametrize("target", ("absolute", "traversal", "packet", "packet-existing", "packet-new"))
def test_output_cannot_escape_root_or_enter_immutable_report_packet(released, target):
    outside = released.parent / (released.name + "_escaped.json")
    relative = {
        "absolute": outside,
        "traversal": "../" + outside.name,
        "packet": report.DIRECTORY,
        "packet-existing": report.DIRECTORY + "/main_paired.json",
        "packet-new": report.DIRECTORY + "/new_contrasts.json",
    }[target]
    before = {p.name: p.read_bytes() for p in (released / report.DIRECTORY).iterdir()}
    with pytest.raises(ValueError, match="repository-relative|immutable four-report packet"):
        report._write_new(released, relative, {"synthetic_only": True})
    assert not outside.exists()
    assert before == {p.name: p.read_bytes() for p in (released / report.DIRECTORY).iterdir()}


def test_output_rejects_symlinked_directory_without_writing_target(released):
    target = released / "ordinary_target"
    target.mkdir()
    link = released / "symlinked_output"
    link.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="escaped its root"):
        report._write_new(released, "symlinked_output/contrast.json", {"synthetic_only": True})
    assert list(target.iterdir()) == []


def _cli(root, *args):
    env = dict(os.environ, PYTHONPATH="", PYTHONNOUSERSITE="1")
    return subprocess.run([sys.executable, "-S", str(root / SCRIPT), *args,
                           "--root", str(root)], cwd=root, env=env,
                          text=True, capture_output=True, timeout=30, check=False)


def test_copied_package_cli_uses_only_stdlib_and_synthetic_release_metadata(released):
    assert {p.name for p in (released / "aav_benchmark").iterdir()} == {
        "__init__.py", "structural_contrast_report.py", "statistics.py"}
    assert not (released / "tools").exists()
    assert not (released / "data/processed").exists()
    before = {role: (released / report.DIRECTORY / (role + ".json")).read_bytes()
              for role in report.ROLES}
    example = _cli(released, "example")
    assert example.returncode == 0, example.stderr
    plan = json.loads(example.stdout)
    assert plan["frozen_before_any_structural_accuracy"] is False
    assert plan["new_structural_accuracy_consulted"] is False
    assert plan["source_code_sha256"] == _read(released, PLAN)["source_code_sha256"]
    bound = _cli(released, "bind", "--plan", PLAN, "--freeze", FREEZE)
    assert bound.returncode == 0, bound.stderr
    pin = json.loads(bound.stdout)
    assert pin["sha256"] == _sha((released / report.BINDING).read_bytes())
    binding = _read(released, report.BINDING)
    assert set(binding["reports"]) == set(report.ROLES)
    for role, original in before.items():
        assert (released / report.DIRECTORY / (role + ".json")).read_bytes() == original
    assert not (released / "tools").exists()
    assert not (released / "data/processed").exists()


@pytest.mark.parametrize("args", (
    ("bind",),
    ("aggregate", "--plan", PLAN, "--freeze", FREEZE),
    ("aggregate", "--plan", PLAN, "--freeze", FREEZE, "--binding-sha256", "0" * 64),
))
def test_cli_requires_explicit_stage_inputs_and_fresh_aggregation_output(released, args):
    result = _cli(released, *args)
    assert result.returncode == 2
    assert "required" in result.stderr or "requires" in result.stderr
    assert not (released / report.BINDING).exists()


def test_portable_cli_refuses_unfrozen_release_without_creating_binding(released):
    plan = _read(released, PLAN)
    plan["frozen_before_any_structural_accuracy"] = False
    _freeze(released, plan)
    result = _cli(released, "bind", "--plan", PLAN, "--freeze", FREEZE)
    assert result.returncode != 0
    assert "freeze before any structural accuracy" in result.stderr
    assert not (released / report.BINDING).exists()
