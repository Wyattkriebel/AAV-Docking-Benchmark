"""Compare derived tables while retaining strict scientific identities/counts."""

import csv
import json
import math
from pathlib import Path

from .provenance import sha256

_DERIVED_FLOAT_FIELDS = frozenset({
    "rmsd_A", "top1_rmsd_A", "top5_min_rmsd_A", "best_returned_rmsd_A",
    "contact_recall", "contact_precision", "top1_contact_recall", "top1_contact_precision",
    "interface_residue_recall", "interface_residue_precision", "interface_residue_jaccard",
    "evaluation_min_distance_A", "seed_mean_top1_score", "spearman_rho", "kendall_tau_b", "pkd",
})


def compare_tables(expected_path, actual_path, numeric_tolerance=0.0):
    """Byte exact by default; explicit tolerance permits only derived float noise.

    Input/pose hashes are independently checked before table comparison. Boolean
    decisions, integer counts, atom/case identities and collection shapes must
    match exactly even when a cross-platform tolerance is requested.
    """
    if not math.isfinite(numeric_tolerance) or not 0 <= numeric_tolerance <= 1e-8:
        raise ValueError("Numeric tolerance must be finite and between 0 and 1e-8")
    expected_path, actual_path = Path(expected_path), Path(actual_path)
    expected_sha, actual_sha = sha256(expected_path), sha256(actual_path)
    if expected_sha == actual_sha:
        return {"exact_recomputation": True, "max_numeric_difference": 0.0, "recomputed_sha256": actual_sha}
    if numeric_tolerance == 0:
        raise ValueError(f"Recomputed table bytes differ: {expected_path.name}")
    differences = []

    def compare(a, b, field=None):
        if type(a) is not type(b):
            raise ValueError("Derived table types differ")
        if isinstance(a, float):
            if field not in _DERIVED_FLOAT_FIELDS:
                if a != b:
                    raise ValueError("Copied measurement, resource or provenance value differs")
                return
            if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > numeric_tolerance:
                raise ValueError("Derived floating value differs beyond tolerance")
            differences.append(abs(a - b))
        elif isinstance(a, dict):
            if a.keys() != b.keys():
                raise ValueError("Derived table fields differ")
            for key in a:
                compare(a[key], b[key], key)
        elif isinstance(a, list):
            if len(a) != len(b):
                raise ValueError("Derived table row counts differ")
            for left, right in zip(a, b):
                compare(left, right, field)
        elif a != b:
            raise ValueError("Derived table identity, count or decision differs")

    if expected_path.suffix == ".json":
        compare(json.loads(expected_path.read_text()), json.loads(actual_path.read_text()))
    elif expected_path.suffix == ".csv":
        def rows(path):
            with path.open(newline="") as stream:
                data = list(csv.reader(stream))
            if not data:
                raise ValueError("Empty derived CSV")
            if any(len(row) != len(data[0]) for row in data[1:]):
                raise ValueError("Malformed derived CSV row width")
            return data
        left, right = rows(expected_path), rows(actual_path)
        if len(left) != len(right) or left[0] != right[0]:
            raise ValueError("Derived CSV fields or row counts differ")
        for a, b in zip(left[1:], right[1:]):
            for column, avalue, bvalue in zip(left[0], a, b):
                if column in _DERIVED_FLOAT_FIELDS:
                    compare(float(avalue), float(bvalue), column)
                elif avalue != bvalue:
                    raise ValueError("Derived CSV identity, count or copied value differs")
    else:
        raise ValueError("Only JSON and CSV derived tables support numeric comparison")
    return {"exact_recomputation": False, "max_numeric_difference": max(differences, default=0.0),
            "recomputed_sha256": actual_sha}
