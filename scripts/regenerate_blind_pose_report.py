#!/usr/bin/env python3
"""Verify the complete frozen run denominator and regenerate qualified pose results."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.qualified_pose import write_qualified_pose_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--numeric-tolerance", type=float, default=0.0,
                        help="explicit tolerance for named derived JSON geometry; copied coordinates/scores/counts/decisions and Markdown remain exact")
    parser.add_argument("--allow-pending", action="store_true", help="emit pending workflow status only, never partial atomic accuracy")
    args = parser.parse_args(argv)
    try:
        paths = write_qualified_pose_report(args.root, summary_output=args.summary_output,
                                            markdown_output=args.markdown_output, check=args.check, allow_pending=args.allow_pending,
                                            numeric_tolerance=args.numeric_tolerance)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"Qualified pose report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + ", ".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
