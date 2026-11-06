#!/usr/bin/env python3
"""Compare two existing report JSON files; no scientific calculation is run."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.report_replay import ReportReplayError, compare_report_files


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expected", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--structural", action="store_true", help="Exact JSON types/values and indexed list order; ignore whitespace/object key order")
    parser.add_argument("--numeric-tolerance", type=float, help="Explicit absolute derived-metric tolerance, 0 to 1e-8; all other paths exact")
    parser.add_argument("--differences-limit", type=int, default=50)
    args = parser.parse_args(argv)
    try:
        result = compare_report_files(args.expected, args.actual, numeric_tolerance=args.numeric_tolerance,
                                      structural=args.structural, differences_limit=args.differences_limit)
    except ReportReplayError as error:
        print(json.dumps(dict(status="comparison_rejected", path=error.path, reason=error.reason,
                              absolute_difference=error.difference, absolute_difference_exact_ratio=error.exact_difference), allow_nan=False))
        return 1
    except OSError as error:
        print(json.dumps(dict(status="input_error", reason=str(error))))
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
