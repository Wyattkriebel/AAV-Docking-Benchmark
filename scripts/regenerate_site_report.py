#!/usr/bin/env python3
"""Validate and render the descriptive published-site localization report."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.site_report import write_site_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true", help="fail if the checked-in report is stale")
    args = parser.parse_args()
    try:
        output = write_site_report(args.root, markdown_output=args.markdown_output, check=args.check)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"Site report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + str(output))


if __name__ == "__main__":
    main()
