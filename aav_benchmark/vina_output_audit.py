"""Pure I/O checks for ordinary-verbosity AutoDock Vina 1.2.7 output.

Only MODEL ranks, score tokens, execution metadata and complete-file hashes are
read. Atom records, experimental references and pose accuracy are not parsed.
This module deliberately does not import or change the frozen metric readers.
"""
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, DecimalException, ROUND_HALF_EVEN, localcontext
import hashlib
import json
import os
from pathlib import Path
import re


DEFAULT_SCOPES = (
    "blind_pose", "blind_pose_neighborhood8", "structural_pose_v3",
    "pol_specificity", "gal9_specificity", "development", "shards",
    "clean_environment",
)
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
MODE_ROW = re.compile(r"^\s*(\d+)\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s*$")
MODEL = re.compile(r"^MODEL\s+(\d+)\s*$")
POSE_SCORE = re.compile(r"^REMARK VINA RESULT:\s*([+-]?\d+\.\d{3})\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s*$")
SOURCE_URL = "https://github.com/ccsb-scripps/AutoDock-Vina/blob/v1.2.7/src/lib/vina.cpp"


@dataclass(frozen=True)
class _Bounds:
    lower: Decimal
    upper: Decimal
    lower_inclusive: bool = True
    upper_inclusive: bool = True

    def __getitem__(self, index):
        return (self.lower, self.upper)[index]

    def __iter__(self):
        return iter((self.lower, self.upper))

    def nonempty(self):
        return self.lower < self.upper or (self.lower == self.upper and self.lower_inclusive and self.upper_inclusive)


def _intersection(a, b):
    lower, upper = max(a.lower, b.lower), min(a.upper, b.upper)
    return _Bounds(lower, upper,
                   (a.lower != lower or a.lower_inclusive) and (b.lower != lower or b.lower_inclusive),
                   (a.upper != upper or a.upper_inclusive) and (b.upper != upper or b.upper_inclusive))


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _decimal(token):
    _require(len(token) <= 400 and re.fullmatch(NUMBER, token), "Invalid finite score token")
    try:
        value = Decimal(token)
    except DecimalException as error:
        raise ValueError("Invalid finite score token") from error
    _require(value.is_finite() and (not value or -324 <= value.adjusted() <= 308),
             "Score token is outside the supported decimal exponent envelope")
    return value if value else Decimal(0)


def _text(value):
    text = format(value, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def _range_json(bounds):
    return {"lower": _text(bounds[0]), "upper": _text(bounds[1]),
            "lower_inclusive": bounds.lower_inclusive, "upper_inclusive": bounds.upper_inclusive}


def _log_bounds(token):
    """Nearest-even rounding bounds for four significant digits.

    General formatting suppresses trailing zeroes. At powers of ten the
    predecessor spacing is ten times smaller than the successor spacing.
    Nonzero tiny values use scientific notation, so a logged zero is exact.
    Endpoint inclusion is checked against nearest-even formatting, preventing
    incompatible adjacent score bins from sharing an artificial tie point.
    """
    value = _decimal(token)
    if not value:
        return _Bounds(value, value)
    magnitude = abs(value)
    exponent = magnitude.adjusted()
    unit = Decimal(1).scaleb(exponent - 3)
    _require(magnitude % unit == 0, "Log score has more than four significant digits")
    down = unit / 10 if magnitude == Decimal(1).scaleb(exponent) else unit
    positive = magnitude - down / 2, magnitude + unit / 2
    lower, upper = positive if value > 0 else (-positive[1], -positive[0])
    return _Bounds(lower, upper, Decimal(format(lower, ".4g")) == value, Decimal(format(upper, ".4g")) == value)


def _saved_bounds(token):
    _require(re.fullmatch(r"[+-]?\d+\.\d{3}", token), "PDBQT score must have three decimal places")
    value = _decimal(token)
    lower, upper = value - Decimal("0.0005"), value + Decimal("0.0005")
    return _Bounds(lower, upper, Decimal(format(lower, ".3f")) == value, Decimal(format(upper, ".3f")) == value)


def score_token_bounds(token, *, format_kind):
    """Expose exact decimal token intervals without importing scientific code."""
    with localcontext() as context:
        context.prec = 700
        context.rounding = ROUND_HALF_EVEN
        if format_kind == "log_four_significant_digits":
            bounds = _log_bounds(token)
        elif format_kind == "pdbqt_three_decimal_places":
            bounds = _saved_bounds(token)
        else:
            raise ValueError("Unsupported score format")
        return _range_json(bounds)


def _log_rows(raw):
    text = raw.decode("utf-8")
    _require(re.findall(r"^\s*AutoDock Vina\s+(\S+)\s*$", text, re.MULTILINE) == ["v1.2.7"],
             "Expected exactly one AutoDock Vina v1.2.7 banner")
    headers = [i for i, line in enumerate(text.splitlines())
               if re.fullmatch(r"\s*mode\s*\|\s*affinity\s*\|\s*dist from best mode\s*", line)]
    _require(len(headers) == 1, "Expected exactly one Vina score table")
    lines = text.splitlines()[headers[0] + 1:]
    _require(len(lines) >= 2 and re.fullmatch(r"\s*\|\s*\(kcal/mol\)\s*\|\s*rmsd l\.b\.\|\s*rmsd u\.b\.\s*", lines[0]),
             "Malformed Vina score table units header")
    _require(re.fullmatch(r"\s*-+\+-+\+-+\+-+\s*", lines[1]), "Malformed Vina score table separator")
    rows = []
    for line in lines[2:]:
        if not line.strip():
            continue
        match = MODE_ROW.fullmatch(line)
        _require(match is not None, "Malformed row or unexpected content after Vina score table")
        # The two distance tokens are syntactically checked, never interpreted
        # as an experimental comparison or used to calculate a metric.
        rows.append({"rank": int(match[1]), "score_token": match[2]})
    return rows


def _saved_rows(raw):
    rows = []
    current = None
    for line in raw.decode("utf-8").splitlines():
        if line.lstrip().startswith(("MODEL", "REMARK VINA RESULT", "ENDMDL")):
            _require(line == line.lstrip(), "Indented reserved MODEL/score/ENDMDL control record")
        if line.startswith("MODEL"):
            match = MODEL.fullmatch(line)
            _require(match is not None and current is None, "Malformed or nested saved MODEL")
            current = {"rank": int(match[1]), "scores": []}
        elif line.startswith("REMARK VINA RESULT"):
            match = POSE_SCORE.fullmatch(line)
            _require(match is not None and current is not None, "Malformed Vina result or score outside MODEL")
            current["scores"].append(match[1])
        elif line.startswith("ENDMDL"):
            _require(line.strip() == "ENDMDL" and current is not None and len(current["scores"]) == 1,
                     "Saved MODEL must contain exactly one score and ENDMDL")
            rows.append({"rank": current["rank"], "score_token": current["scores"][0]})
            current = None
    _require(current is None, "Unterminated saved MODEL")
    return rows


def validate_vina_1_2_7_output(log_raw, pose_raw, *, energy_range_kcal_mol=10,
                             expected_returned_poses=None, num_modes=20):
    """Validate a saved rank prefix, allowing only writer-compatible suffixes.

    A boundary-compatible omission is accepted but explicitly ambiguous. It is
    never counted as proven outside. The writer retains energies <= best+10
    and omits energies > best+10. Bounds also require one shared latent best
    energy, rather than independently choosing a different best for each row.
    """
    with localcontext() as context:
        context.prec = 700
        context.rounding = ROUND_HALF_EVEN
        window = _decimal(str(energy_range_kcal_mol))
        _require(window == 10, "This audit requires the explicit frozen 10 kcal/mol window")
        _require(type(num_modes) is int and 1 <= num_modes <= 20, "Invalid frozen maximum mode count")
        logs, saved = _log_rows(log_raw), _saved_rows(pose_raw)
        _require(1 <= len(saved) <= len(logs) <= num_modes, "Invalid logged/saved mode counts")
        if expected_returned_poses is not None:
            _require(type(expected_returned_poses) is int, "returned_poses metadata must be an integer")
            _require(expected_returned_poses == len(saved), "returned_poses metadata disagrees with saved MODEL count")
        for rows, label in ((logs, "log"), (saved, "saved")):
            _require([row["rank"] for row in rows] == list(range(1, len(rows) + 1)),
                     f"Noncontiguous or reordered {label} ranks")
            values = [_decimal(row["score_token"]) for row in rows]
            _require(all(a <= b for a, b in zip(values, values[1:])), f"Unordered {label} scores")
        intervals = [_log_bounds(row["score_token"]) for row in logs]
        comparisons = []
        for index, row in enumerate(saved):
            intervals[index] = _intersection(intervals[index], _saved_bounds(row["score_token"]))
            _require(intervals[index].nonempty(), "Saved/log score prefix has incompatible rounding intervals")
            comparisons.append({"rank": row["rank"], "log_score_token": logs[index]["score_token"],
                                "saved_score_token": row["score_token"], "shared_score_bounds": _range_json(intervals[index])})
        running_lower, running_inclusive = intervals[0].lower, intervals[0].lower_inclusive
        for interval in intervals:
            if interval.lower > running_lower:
                running_lower, running_inclusive = interval.lower, interval.lower_inclusive
            elif interval.lower == running_lower:
                running_inclusive = running_inclusive and interval.lower_inclusive
            _require(running_lower < interval.upper or
                     (running_lower == interval.upper and running_inclusive and interval.upper_inclusive),
                     "No score-ordered realization of score token intervals")
        best = intervals[0]
        feasible = best
        omitted = []
        saved_boundary = []
        for index, interval in enumerate(intervals):
            lower, upper = interval
            delta = _Bounds(lower - best[1], upper - best[0], interval.lower_inclusive and best.upper_inclusive,
                            interval.upper_inclusive and best.lower_inclusive)
            feasible = _intersection(feasible, _Bounds(Decimal("-Infinity"), upper, True, interval.upper_inclusive))
            if index < len(saved):
                _require(delta[0] < window or (delta[0] == window and delta.lower_inclusive),
                         "Saved model is proven outside the output energy window")
                feasible = _intersection(feasible, _Bounds(lower - window, Decimal("Infinity"), interval.lower_inclusive))
                if delta[1] > window:
                    saved_boundary.append(logs[index]["rank"])
            else:
                _require(delta[1] > window, "Omitted logged mode is inside the output energy window")
                feasible = _intersection(feasible, _Bounds(Decimal("-Infinity"), upper - window, True, False))
                omitted.append({"rank": logs[index]["rank"], "log_score_token": logs[index]["score_token"],
                                "delta_from_best_bounds_kcal_mol": _range_json(delta),
                                "window_classification": "proven_outside_window" if delta[0] > window or
                                (delta[0] == window and not delta.lower_inclusive)
                                else "boundary_ambiguous_writer_compatible"})
        _require(feasible.nonempty(),
                 "No single best score makes the complete saved prefix and omitted suffix writer-compatible")
        proven = sum(row["window_classification"] == "proven_outside_window" for row in omitted)
        return {
            "status": "valid_writer_compatible", "source_version": "AutoDock Vina v1.2.7",
            "energy_range_kcal_mol": "10", "logged_modes": len(logs), "saved_models": len(saved),
            "omitted_modes": len(omitted), "proven_outside_window_omissions": proven,
            "boundary_ambiguous_omissions": len(omitted) - proven, "saved_boundary_ambiguous_ranks": saved_boundary,
            "saved_ranks_are_complete_log_prefix": True, "score_prefix": comparisons, "omitted_suffix": omitted,
            "best_score_bounds_kcal_mol": _range_json(best),
            "shared_writer_compatible_best_bounds_kcal_mol": _range_json(feasible),
            "log": {"sha256": hashlib.sha256(log_raw).hexdigest(), "bytes": len(log_raw)},
            "poses": {"sha256": hashlib.sha256(pose_raw).hexdigest(), "bytes": len(pose_raw)},
            "format_convention": {"log": "defaultfloat_four_significant_digits_verbosity_1",
                                  "saved_score": "fixed_three_decimal_places", "source": SOURCE_URL,
                                  "rounding_endpoints": "nearest_even; ambiguous window boundaries are reported"},
            "experimental_reference_access": False, "accuracy_computed": False,
        }


def _read_stable(path):
    raw = path.read_bytes()
    _require(path.read_bytes() == raw, "Raw artifact changed during read")
    return raw


def _option(command, flag):
    _require(isinstance(command, list) and command.count(flag) == 1, f"Expected exactly one explicit {flag}")
    index = command.index(flag)
    _require(index + 1 < len(command), f"Missing value for {flag}")
    return str(command[index + 1])


def _artifact(root, directory, filename, suffix):
    _require(isinstance(filename, str), "Missing raw artifact filename")
    _require(Path(filename).name == filename and Path(filename).suffix == suffix, "Unexpected raw artifact filename or suffix")
    original = directory / filename
    _require(not original.is_symlink(), "Raw artifact is symlinked or outside its run directory")
    path = original.resolve()
    _require(path == original and path.is_relative_to(root) and path.parent == directory.resolve(),
             "Raw artifact is symlinked or outside its run directory")
    return path


def audit_main_vina_outputs(root, *, scopes=DEFAULT_SCOPES):
    """Read current main/control run lists; emit a reproducible pure-I/O snapshot.

    Historical recovery copies and supplemental assembly outputs are excluded.
    Live running rows remain explicit and are never interpreted as failures.
    This function writes nothing and never opens preparation/reference files.
    """
    root = Path(root).resolve()
    _require(bool(scopes) and len(set(scopes)) == len(scopes) and set(scopes) <= set(DEFAULT_SCOPES),
             "Only the named existing main/control scopes are permitted")
    records, summaries = [], {}
    for scope in scopes:
        counts = Counter()
        scope_root = root / "results" / scope
        try:
            _require(not scope_root.is_symlink() and scope_root.resolve() == scope_root,
                     "Symlinked audit scope is not permitted")
        except (ValueError, OSError, RuntimeError) as error:
            records.append({"scope": scope, "status": "invalid_scope", "error": str(error)})
            counts["invalid_records"] += 1
            summaries[scope] = dict(counts)
            continue
        if not scope_root.exists():
            counts["missing_scopes"] += 1
            summaries[scope] = dict(counts)
            continue
        manifests = []
        def walk_error(error):
            records.append({"scope": scope, "status": "invalid_directory", "error": str(error)})
            counts["invalid_records"] += 1
        for directory, directories, filenames in os.walk(scope_root, followlinks=False, onerror=walk_error):
            for name in list(directories):
                child = Path(directory) / name
                if name in {"_recovery", "vina_carb"}:
                    directories.remove(name)
                elif child.is_symlink():
                    directories.remove(name)
                    records.append({"scope": scope, "status": "invalid_directory",
                                    "error": "Symlinked audit directory is not permitted", "path": str(child.relative_to(root))})
                    counts["invalid_records"] += 1
            directories.sort()
            if "runs.json" in filenames:
                manifests.append(Path(directory) / "runs.json")
        for manifest in sorted(manifests):
            counts["run_manifests"] += 1
            try:
                _require(not manifest.is_symlink() and manifest.resolve() == manifest and manifest.resolve().is_relative_to(scope_root),
                         "Symlinked or outside-scope run manifest is not permitted")
                manifest_raw = _read_stable(manifest)
                rows = json.loads(manifest_raw)
                _require(isinstance(rows, list), "Expected raw run-list JSON")
            except (ValueError, OSError, UnicodeError, RuntimeError) as error:
                records.append({"scope": scope, "manifest": str(manifest.relative_to(root)),
                                "status": "invalid_manifest", "error": str(error)})
                counts["invalid_records"] += 1
                continue
            manifest_pin = {"path": str(manifest.relative_to(root)), "sha256": hashlib.sha256(manifest_raw).hexdigest(),
                            "bytes": len(manifest_raw)}
            for index, row in enumerate(rows):
                if not isinstance(row, dict):
                    records.append({"scope": scope, "manifest": manifest_pin, "record_index": index,
                                    "status": "invalid_record", "error": "Run-list record is not an object"})
                    counts["records"] += 1
                    counts["invalid_records"] += 1
                    continue
                record = {"scope": scope, "manifest": manifest_pin, "record_index": index,
                          "seed": row.get("seed"), "receptor": row.get("receptor"), "source_status": row.get("status")}
                counts["records"] += 1
                counts["source_status_" + str(row.get("status"))] += 1
                if row.get("status") != "ok":
                    record["status"] = "pending_source" if row.get("status") == "running" else "terminal_source_failure"
                    records.append(record)
                    continue
                if row.get("vina_version") != "AutoDock Vina v1.2.7":
                    record["status"] = "other_version_not_audited"
                    counts["other_version_records"] += 1
                    records.append(record)
                    continue
                try:
                    command = row["command"]
                    _require(row.get("exit_code") == 0, "Successful record has a nonzero native exit code")
                    _require(type(row.get("returned_poses")) is int and row["returned_poses"] >= 1,
                             "Successful record lacks a positive returned_poses count")
                    _require(_option(command, "--scoring") == "vina", "Expected the Vina scoring function")
                    if "--verbosity" in command:
                        _require(_option(command, "--verbosity") == "1", "Only ordinary verbosity 1 formatting is supported")
                    mode_token = _option(command, "--num_modes")
                    _require(mode_token.isdigit(), "Invalid --num_modes token")
                    pose_path = _artifact(root, manifest.parent, row.get("pose_file"), ".pdbqt")
                    log_path = _artifact(root, manifest.parent, row.get("log_file"), ".log")
                    pose_raw, log_raw = _read_stable(pose_path), _read_stable(log_path)
                    for label, path, raw in (("poses", pose_path, pose_raw), ("log", log_path, log_raw)):
                        actual = hashlib.sha256(raw).hexdigest()
                        record[label] = {"path": str(path.relative_to(root)), "sha256": actual, "bytes": len(raw)}
                        expected = row.get("pose_sha256" if label == "poses" else "log_sha256")
                        if expected is not None:
                            _require(actual == expected, f"Raw {label} metadata hash mismatch")
                    record["validation"] = validate_vina_1_2_7_output(
                        log_raw, pose_raw, energy_range_kcal_mol=_option(command, "--energy_range"),
                        expected_returned_poses=row.get("returned_poses"), num_modes=int(mode_token))
                    record["status"] = "valid_writer_compatible"
                    counts["audited_records"] += 1
                    for key in ("logged_modes", "saved_models", "omitted_modes", "proven_outside_window_omissions", "boundary_ambiguous_omissions"):
                        counts[key] += record["validation"][key]
                except (ValueError, OSError, UnicodeError, KeyError, TypeError, RuntimeError) as error:
                    record.update(status="invalid_output", error=str(error))
                    counts["invalid_records"] += 1
                records.append(record)
        summaries[scope] = dict(sorted(counts.items()))
    totals = Counter()
    for counts in summaries.values():
        totals.update(counts)
    return {"schema_version": 1, "audit": "existing_main_vina_1_2_7_log_saved_prefix",
            "status": "invalid_output_present" if totals["invalid_records"] else
            ("all_audited_outputs_writer_compatible" if totals["audited_records"] else "no_auditable_outputs_found"),
            "scopes": list(scopes), "excluded": ["historical _recovery copies", "supplemental assembly outputs", "Vina-Carb outputs"],
            "scope_counts": summaries, "totals": dict(sorted(totals.items())), "records": records,
            "experimental_reference_access": False, "accuracy_computed": False, "native_execution": False,
            "frozen_helpers_modified": False, "eligibility_changed": False}
