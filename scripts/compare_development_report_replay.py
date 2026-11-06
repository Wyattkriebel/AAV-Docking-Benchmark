#!/usr/bin/env python3
"""Compare two existing qualified-development JSON files; run no metrics."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.development_report_replay import compare_development_report_files
from aav_benchmark.report_replay import ReportReplayError


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expected", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--structural", action="store_true")
    parser.add_argument("--numeric-tolerance", type=float, help="Explicit absolute tolerance for approved development metric paths only, at most1e-8")
    parser.add_argument("--differences-limit", type=int, default=50)
    args = parser.parse_args(argv)
    try:
        result = compare_development_report_files(args.expected, args.actual, numeric_tolerance=args.numeric_tolerance,
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
