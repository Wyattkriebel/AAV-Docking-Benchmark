"""Read-only replay comparison for the frozen qualified development summary.

This separate path policy reuses the frozen strict JSON parser and exact-ratio
comparison core. It never imports the empirical evaluator, changes a global
policy, reads a reference or computes a scientific metric.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
from pathlib import Path

from . import report_replay as _core


POLICY_ID = "qualified_development_derived_paths_v1"
EVALUATOR_SHA256 = "c77333c8fbda7fe81f79fbe99211a65f3c91f51382daaa80f2fb89bd12dfc525"
CORE_SHA256 = "a6340ec4779fa546eb5673fecaca7b976415b78dbdfb1804b1b3c55c8f2bb048"
_INDEX = object()
_ALIGNMENTS = ("global", "local_interface")
_SUCCESS_FIELDS = tuple(f"{rank}_success_{threshold}A" for rank in ("top1", "top5", "best_returned") for threshold in (1, 2, 3))


def _rules():
    paths = []
    condition = ("conditions", _INDEX)
    for budget in ("1", "5"):
        prefix = condition + ("site_budgets", budget)
        seed = prefix + ("seed_evaluations", _INDEX)
        pose = seed + ("poses", _INDEX)
        for field in ("global_RMSD_A", "local_interface_RMSD_A", "full_parent_min_protein_distance_A"):
            paths.append(pose + (field,))
        for alignment in _ALIGNMENTS:
            for field in ("contact_recall", "contact_precision", "interface_residue_recall",
                          "interface_residue_precision", "interface_residue_jaccard"):
                paths.append(pose + (alignment + "_contacts", field))
            for field in ("top1_RMSD_A", "top5_min_RMSD_A", "best_returned_RMSD_A"):
                paths.append(seed + (alignment, field))
            for field in _SUCCESS_FIELDS:
                paths.append(prefix + (alignment, field))  # Fractional seed averages, not seed booleans.
        geometry = pose + ("geometry_diagnostics",)
        for field in ("max_heavy_bond_length_change_from_input_A", "max_ring_pair_distance_change_from_input_A"):
            paths.append(geometry + (field,))
        for field in ("distance_A", "input_distance_A"):
            paths.append(geometry + ("glycosidic_links", _INDEX, field))
    for alignment in _ALIGNMENTS:
        fit = condition + ("receptor_alignments", alignment)
        paths.extend((fit + ("backbone_RMSD_A",), fit + ("rotation", _INDEX, _INDEX),
                      fit + ("translation", _INDEX)))
    paths.append(condition + ("native_controls", "known_translation_RMSD_A"))
    for scope in ("full_ligand", "modeled_core"):
        for method in ("fpocket", "p2rank", "oracle"):
            for budget in ("1", "5"):
                for alignment in _ALIGNMENTS:
                    for field in _SUCCESS_FIELDS:
                        paths.append(("aggregates", scope, method, budget, alignment, field))
    paths.extend((("paired_contrasts", _INDEX, "mean_oracle_minus_predicted"),
                  ("paired_contrasts", _INDEX, "complexes", _INDEX, "oracle_minus_predicted")))
    by_final = defaultdict(list)
    for path in dict.fromkeys(paths):
        by_final[path[-1]].append(path)
    return dict(by_final)


_RULES_BY_FINAL = _rules()


def _rule(path):
    if not path:
        return None
    final = _INDEX if type(path[-1]) is int else path[-1]
    for candidate in _RULES_BY_FINAL.get(final, ()):
        if len(path) == len(candidate) and all(type(part) is int if token is _INDEX else token == part for token, part in zip(candidate, path)):
            return _core.json_pointer(tuple("[]" if x is _INDEX else x for x in candidate))
    return None


def _check_core():
    if hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest() != CORE_SHA256:
        raise _core.ReportReplayError("Frozen strict comparison core changed")


def _project(value, allowed, path=()):
    """Mask only approved float leaves; retain every other value and position."""
    if type(value) is float and _rule(path) is not None:
        allowed[path] = value
        return None
    if type(value) is dict:
        return {key: _project(item, allowed, path + (key,)) for key, item in value.items()}
    if type(value) is list:
        return [_project(item, allowed, path + (index,)) for index, item in enumerate(value)]
    return value


def _original_path(carrier_path, paths):
    parts = carrier_path.split("/")
    if len(parts) == 4 and parts[1] == "per_source_pose_metrics" and parts[2].isdigit() and parts[3] == "global_RMSD_A":
        index = int(parts[2])
        if index < len(paths):
            return paths[index]
    raise _core.ReportReplayError("Unexpected scalar-comparator diagnostic path")


def compare_development_report_values(expected, actual, *, numeric_tolerance=0.0, differences_limit=50):
    """Exact original structure plus separately compared approved scalar pairs.

    The internal carrier is solely an adapter to immutable comparison arithmetic;
    it does not create poses, RMSDs or metric data. The original scalar path sets
    and types must agree before the original protected tree is compared exactly.
    """
    _check_core()
    tolerance = _core._tolerance(numeric_tolerance)
    _core._validate_complete(expected); _core._validate_complete(actual)
    left, right = {}, {}
    try:
        left_tree, right_tree = _project(expected, left), _project(actual, right)
    except RecursionError as error:
        raise _core.ReportReplayError("Cyclic or excessively nested development JSON value") from error
    if left.keys() != right.keys():
        differing = next(path for path in (*left, *right) if (path in left) != (path in right))
        raise _core.ReportReplayError("Approved derived floating path or type changed", differing)
    protected = _core.compare_report_values(left_tree, right_tree)
    paths = list(left)
    carrier = lambda values: {"per_source_pose_metrics": [{"global_RMSD_A": values[path]} for path in paths]}
    try:
        result = _core.compare_report_values(carrier(left), carrier(right), numeric_tolerance=tolerance,
                                              differences_limit=differences_limit)
    except _core.ReportReplayError as error:
        # Configuration errors have no carrier path and already identify the issue.
        if not error.path:
            raise
        path = _original_path(error.path, paths)
        difference = Fraction(error.exact_difference) if error.exact_difference is not None else None
        raise _core.ReportReplayError(error.reason, path, difference=difference) from error
    for difference in result["displayed_differences"]:
        path = _original_path(difference["path"], paths)
        difference.update(path=_core.json_pointer(path), approved_path_rule=_rule(path))
    if result["max_difference_path"] is not None:
        result["max_difference_path"] = _core.json_pointer(_original_path(result["max_difference_path"], paths))
    result.update(policy_id=POLICY_ID, evaluator_schema_source_sha256=EVALUATOR_SHA256,
                  comparison_core_sha256=CORE_SHA256, compared_leaf_values=protected["compared_leaf_values"],
                  compared_float_values=protected["compared_float_values"] + len(paths),
                  approved_derived_float_values=len(paths), global_comparison_policy_modified=False)
    return result


def compare_development_report_files(expected_path, actual_path, *, numeric_tolerance=None, structural=False, differences_limit=50):
    """Compare two existing development JSON files; default is exact bytes."""
    if type(structural) is not bool:
        raise _core.ReportReplayError("Structural mode must be boolean")
    if numeric_tolerance is not None:
        numeric_tolerance = _core._tolerance(numeric_tolerance)
    expected_raw, actual_raw = Path(expected_path).read_bytes(), Path(actual_path).read_bytes()
    expected, actual = _core.load_report_json(expected_raw), _core.load_report_json(actual_raw)
    exact = expected_raw == actual_raw
    if not structural and numeric_tolerance is None and not exact:
        raise _core.ReportReplayError("Report bytes differ; structural/numeric comparison was not requested")
    result = compare_development_report_values(expected, actual, numeric_tolerance=numeric_tolerance or 0.0,
                                                differences_limit=differences_limit)
    if not structural and numeric_tolerance is None:
        result.update(status="exact_byte_match", comparison_mode="exact_bytes")
    result.update(exact_bytes=exact,
        expected_file=dict(path=str(expected_path), sha256=hashlib.sha256(expected_raw).hexdigest(), bytes=len(expected_raw)),
        actual_file=dict(path=str(actual_path), sha256=hashlib.sha256(actual_raw).hexdigest(), bytes=len(actual_raw)))
    return result
