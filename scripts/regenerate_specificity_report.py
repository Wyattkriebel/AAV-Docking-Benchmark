#!/usr/bin/env python3
"""Recompute the Gal9C specificity report from frozen inputs and raw run files."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.specificity_report import write_specificity_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--prepared-root", type=Path, help="prepared case directories; ordinary or per-file gzip archives")
    parser.add_argument("--results-root", type=Path, help="raw run case directories and pre-run freeze receipt")
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true", help="fail if either checked-in output is stale")
    parser.add_argument("--numeric-tolerance", type=float, default=0.0,
                        help="explicit derived-float tolerance for JSON checks; Markdown remains exact")
    args = parser.parse_args()
    try:
        paths = write_specificity_report(args.root, prepared_root=args.prepared_root, results_root=args.results_root,
                                         summary_output=args.summary_output, markdown_output=args.markdown_output,
                                         check=args.check, numeric_tolerance=args.numeric_tolerance)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"Specificity report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + ", ".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
