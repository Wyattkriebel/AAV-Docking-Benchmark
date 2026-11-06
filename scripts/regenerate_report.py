#!/usr/bin/env python3
"""Regenerate the development report without rerunning or changing docking."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.report import write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true", help="fail if the generated report is stale")
    args = parser.parse_args()
    try:
        path = write_report(args.root, args.output, check=args.check)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, f"Report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + str(path))


if __name__ == "__main__":
    main()
