#!/usr/bin/env python3
"""Verify raw POL evidence and regenerate the descriptive fixed-panel report."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.pol_report import write_pol_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--prepared-root", type=Path)
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        paths = write_pol_report(args.root, prepared_root=args.prepared_root, results_root=args.results_root,
                                 summary_output=args.summary_output, markdown_output=args.markdown_output, check=args.check)
    except (ValueError, OSError, KeyError, RuntimeError, TypeError) as error:
        parser.exit(1, f"POL report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + ", ".join(map(str, paths)))


if __name__ == "__main__":
    main()
