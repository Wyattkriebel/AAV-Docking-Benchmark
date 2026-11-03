"""Compare saved report JSON without reading references or executing metrics.

Files are byte-exact by default. An explicit bounded absolute tolerance permits
only reviewed derived floating values at complete report paths. No key-name
heuristic, relative tolerance, list sorting or source-data waiver is used.
This comparison does not authenticate or recompute a report's scientific inputs.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


MAX_NUMERIC_TOLERANCE = 1e-8
POLICY_ID = "global_report_derived_paths_v1"
_INDEX = object()
_COPY_ID = object()
_SEARCH_METHODS = ("modern_vina_1_2_7", "vina_carb_1_0")
_RANK_METHODS = (
    "VinaCarb1.0_Vina1.1.2_core_shared_zero_base_coeff0",
    "VinaCarb1.0_Vina1.1.2_core_shared_zero_base_plus_CHI_coeff1",
)
_ALIGNMENTS = ("global", "local_interface")


class ReportReplayError(ValueError):
    """A report is invalid or differs outside the explicit replay policy."""

    def __init__(self, reason, path=(), *, difference=None):
        self.reason = reason
        self.path = json_pointer(path)
        self.exact_difference = str(difference) if difference is not None else None
        try:
            self.difference = float(difference) if difference is not None else None
        except OverflowError:
            self.difference = None
        if self.difference is not None and not math.isfinite(self.difference):
            self.difference = None
        super().__init__(f"{self.path or '/'}: {reason}")


def json_pointer(path):
    """Unambiguous JSON Pointer, retaining object keys and list positions."""
    return "".join("/" + str(part).replace("~", "~0").replace("/", "~1") for part in path)


@dataclass(frozen=True)
class _Rule:
    path: tuple

    @property
    def label(self):
        return json_pointer(tuple("[]" if x is _INDEX else "{copy_id}" if x is _COPY_ID else x for x in self.path))

    def matches(self, path):
        return len(path) == len(self.path) and all(
            (type(value) is int if token is _INDEX else type(value) is str if token is _COPY_ID else value == token)
            for token, value in zip(self.path, path)
        )


def _make_rules():
    paths = []
    condition = ("conditions", _INDEX)
    seed = condition + ("seeds", _INDEX)
    poses = [seed + ("poses", _INDEX), ("per_source_pose_metrics", _INDEX)]
    summaries = []
    for container, methods in (("search_methods", _SEARCH_METHODS), ("ranking_methods", _RANK_METHODS)):
        for method in methods:
            branch = seed + (container, method)
            poses.append(branch + ("poses", _INDEX))
            summaries.extend(branch + (alignment,) for alignment in _ALIGNMENTS)
    summaries.extend(seed + (alignment,) for alignment in _ALIGNMENTS)
    for prefix in poses:
        for field in ("global_RMSD_A", "local_interface_RMSD_A", "full_parent_min_protein_distance_A",
                      "full_observed_assembly_min_protein_distance_A"):
            paths.append(prefix + (field,))
        for alignment in _ALIGNMENTS:
            for field in ("contact_recall", "contact_precision", "interface_residue_recall",
                          "interface_residue_precision", "interface_residue_jaccard"):
                paths.append(prefix + (alignment + "_contacts", field))
        geometry = prefix + ("geometry_diagnostics",)
        for field in ("max_heavy_bond_length_change_from_input_A", "max_ring_pair_distance_change_from_input_A"):
            paths.append(geometry + (field,))
        for field in ("distance_A", "input_distance_A"):
            paths.append(geometry + ("glycosidic_links", _INDEX, field))
    for prefix in summaries:
        for field in ("top1_RMSD_A", "top5_min_RMSD_A", "best_returned_RMSD_A"):
            paths.append(prefix + (field,))
    controls = [condition + ("native_controls", _INDEX),
                seed + ("model_reference_controls", _INDEX, "protein_fits_and_native_controls", _INDEX)]
    for prefix in controls:
        for alignment in _ALIGNMENTS:
            fit = prefix + ("alignments", alignment)
            paths.extend((fit + ("backbone_RMSD_A",), fit + ("rotation", _INDEX, _INDEX),
                          fit + ("translation", _INDEX)))
        paths.append(prefix + ("native_controls", "known_translation_RMSD_A"))
    site = condition + ("site",)
    paths.extend((site + ("candidates", _INDEX, "center_to_native_atom_A"),
                  site + ("candidates", _INDEX, "qualified_copy_distances_A", _COPY_ID)))
    for budget in ("1", "5"):
        paths.append(site + (budget, "minimum_DCA_A"))
        paths.append(site + (budget, "qualified_copy_max_box_atom_coverage", _INDEX))
        for threshold in (3, 4, 5):
            paths.append(site + (budget, f"qualified_copy_coverage_{threshold}A"))
    statistics = [("paired_family_statistics", _INDEX), ("conditional_family_statistics", _INDEX),
                  ("shared_site_statistics", _INDEX, "paired_source_family_intervals", _INDEX)]
    for prefix in statistics:
        paths.extend((prefix + ("difference_a_minus_b",), prefix + ("confidence_interval", _INDEX)))
        for field in ("recovery_a", "recovery_b", "difference"):
            paths.append(prefix + ("per_family", _INDEX, field))
        aggregate = prefix + ("eligible_and_conditional_summary",)
        for level in ("per_complex", "per_family", "overall"):
            for field in ("recovery", "conditional_completed_recovery", "conditional_near_native_recovery",
                          "family_macro_recovery", "complex_macro_recovery"):
                paths.append(aggregate + (level, _INDEX, field))
    return tuple(_Rule(path) for path in dict.fromkeys(paths))


_RULES = _make_rules()
_RULES_BY_FINAL = defaultdict(list)
for _rule in _RULES:
    _RULES_BY_FINAL[_rule.path[-1]].append(_rule)


def _derived_rule(path):
    if not path:
        return None
    final = _INDEX if type(path[-1]) is int else path[-1]
    candidates = _RULES_BY_FINAL.get(final, ())
    if type(path[-1]) is str:
        candidates = (*candidates, *_RULES_BY_FINAL.get(_COPY_ID, ()))
    return next((rule.label for rule in candidates if rule.matches(path)), None)


def _tolerance(value):
    if type(value) not in (int, float) or not 0 <= value <= MAX_NUMERIC_TOLERANCE or not math.isfinite(value):
        raise ReportReplayError("Numeric tolerance must be finite, nonboolean and between 0 and 1e-8")
    return float(value)


def _validate(value, path=()):
    if type(value) is float:
        if not math.isfinite(value):
            raise ReportReplayError("Nonfinite JSON number", path)
    elif type(value) is dict:
        if any(type(key) is not str for key in value):
            raise ReportReplayError("JSON object keys must be strings", path)
        for key, item in value.items():
            _validate(item, path + (key,))
    elif type(value) is list:
        for index, item in enumerate(value):
            _validate(item, path + (index,))
    elif type(value) not in (str, int, bool, type(None)):
        raise ReportReplayError("Unsupported non-JSON value", path)


def _no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ReportReplayError("Duplicate JSON object key: " + key)
        result[key] = value
    return result


def _validate_complete(value):
    try:
        _validate(value)
    except RecursionError as error:
        raise ReportReplayError("Cyclic or excessively nested JSON value") from error


def load_report_json(raw):
    """Parse strict UTF-8 JSON, rejecting duplicate keys and nonfinite values."""
    def constant(value):
        raise ReportReplayError("Nonfinite JSON constant: " + value)
    try:
        if type(raw) is bytes:
            raw = raw.decode("utf-8")
        elif type(raw) is not str:
            raise ReportReplayError("Expected UTF-8 JSON bytes or text")
        value = json.loads(raw, object_pairs_hook=_no_duplicates, parse_constant=constant)
    except ReportReplayError:
        raise
    except (UnicodeError, ValueError, RecursionError) as error:
        raise ReportReplayError("Invalid UTF-8 JSON: " + str(error)) from error
    _validate_complete(value)
    return value


def compare_report_values(expected, actual, *, numeric_tolerance=0.0, differences_limit=50):
    """Compare exact types/shapes/order; permit only approved derived floats.

    Unknown floats, input coordinates, native energy scores, measurements,
    telemetry, constants, counts, identifiers and decisions remain exact.
    Every value is compared even after the diagnostic display limit is reached.
    Signed zero is exact unless it occurs at an explicitly tolerated metric path.
    """
    tolerance = _tolerance(numeric_tolerance)
    if type(differences_limit) is not int or differences_limit < 0:
        raise ReportReplayError("Differences display limit must be a nonnegative integer")
    _validate_complete(expected); _validate_complete(actual)
    exact_tolerance = Fraction.from_float(tolerance)
    maximum_difference = Fraction(0)
    state = dict(compared_leaf_values=0, compared_float_values=0, approved_derived_float_values=0,
                 tolerated_difference_count=0, max_numeric_difference=0.0, max_difference_path=None)
    differences = []

    def compare(left, right, path):
        nonlocal maximum_difference
        if type(left) is not type(right):
            raise ReportReplayError("JSON value types differ", path)
        if type(left) is dict:
            if left.keys() != right.keys():
                raise ReportReplayError("JSON object fields differ", path)
            for key in left:
                compare(left[key], right[key], path + (key,))
        elif type(left) is list:
            if len(left) != len(right):
                raise ReportReplayError("JSON list lengths differ", path)
            for index, (a, b) in enumerate(zip(left, right)):
                compare(a, b, path + (index,))
        else:
            state["compared_leaf_values"] += 1
            if type(left) is float:
                state["compared_float_values"] += 1
                rule = _derived_rule(path)
                if rule is not None:
                    state["approved_derived_float_values"] += 1
                if left.hex() == right.hex():
                    return
                # Compare the exact represented binary values. A rounded
                # subtraction can otherwise accept an over-bound difference
                # (for example tolerance versus a tiny negative statistic).
                difference = abs(Fraction.from_float(left) - Fraction.from_float(right))
                if tolerance == 0 or rule is None:
                    raise ReportReplayError("Exact floating value differs outside the requested derived-roundoff policy", path, difference=difference)
                if difference > exact_tolerance:
                    raise ReportReplayError("Derived floating difference exceeds the absolute tolerance", path, difference=difference)
                state["tolerated_difference_count"] += 1
                if state["max_difference_path"] is None or difference > maximum_difference:
                    maximum_difference = difference
                    state["max_numeric_difference"], state["max_difference_path"] = float(difference), json_pointer(path)
                if len(differences) < differences_limit:
                    differences.append(dict(path=json_pointer(path), expected=left, actual=right,
                                            absolute_difference=float(difference), absolute_difference_exact_ratio=str(difference), approved_path_rule=rule))
            elif left != right:
                raise ReportReplayError("Exact identity, count, decision or value differs", path)

    compare(expected, actual, ())
    return dict(schema_version=1, status="match_with_derived_roundoff" if differences or state["tolerated_difference_count"] else "exact_structural_match",
                comparison_mode="derived_float_tolerance" if tolerance else "exact_structure_and_types",
                policy_id=POLICY_ID, numeric_tolerance=tolerance, relative_tolerance=0.0, **state,
                max_numeric_difference_exact_ratio=str(maximum_difference), tolerance_exact_ratio=str(exact_tolerance),
                displayed_differences=differences, differences_display_truncated=state["tolerated_difference_count"] > len(differences),
                list_order_preserved=True, unknown_fields_exact=True, source_validation_performed=False,
                metrics_recomputed=False)


def compare_report_files(expected_path, actual_path, *, numeric_tolerance=None, structural=False, differences_limit=50):
    """Validate both JSONs; compare bytes unless structural mode is explicit."""
    if type(structural) is not bool:
        raise ReportReplayError("Structural mode must be boolean")
    if numeric_tolerance is not None:
        numeric_tolerance = _tolerance(numeric_tolerance)
    expected_raw, actual_raw = Path(expected_path).read_bytes(), Path(actual_path).read_bytes()
    expected, actual = load_report_json(expected_raw), load_report_json(actual_raw)
    exact = expected_raw == actual_raw
    if not structural and numeric_tolerance is None and not exact:
        raise ReportReplayError("Report bytes differ; structural/numeric comparison was not requested")
    result = compare_report_values(expected, actual, numeric_tolerance=numeric_tolerance or 0.0,
                                   differences_limit=differences_limit)
    if not structural and numeric_tolerance is None:
        result.update(status="exact_byte_match", comparison_mode="exact_bytes")
    result.update(exact_bytes=exact,
                  expected_file=dict(path=str(expected_path), sha256=hashlib.sha256(expected_raw).hexdigest(), bytes=len(expected_raw)),
                  actual_file=dict(path=str(actual_path), sha256=hashlib.sha256(actual_raw).hexdigest(), bytes=len(actual_raw)))
    return result
