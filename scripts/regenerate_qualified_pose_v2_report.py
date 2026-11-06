#!/usr/bin/env python3
"""Regenerate atom-order-corrected OLD results under a pinned post-outcome amendment."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.qualified_pose_v2 import write_qualified_pose_report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--implementation-amendment", type=Path, required=True,
                        help="repository-relative post-outcome implementation amendment JSON")
    parser.add_argument("--implementation-amendment-sha256", required=True)
    parser.add_argument("--numeric-tolerance", type=float, default=0.0,
                        help="explicit tolerance for named derived JSON geometry; copied coordinates/scores/counts/decisions and Markdown remain exact")
    parser.add_argument("--allow-pending", action="store_true", help="emit pending workflow status only, never partial atomic accuracy")
    args = parser.parse_args(argv)
    try:
        paths = write_qualified_pose_report(args.root, summary_output=args.summary_output,
                                            markdown_output=args.markdown_output, check=args.check, allow_pending=args.allow_pending,
                                            numeric_tolerance=args.numeric_tolerance,
                                            implementation_amendment=args.implementation_amendment,
                                            implementation_amendment_sha256=args.implementation_amendment_sha256)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"Qualified pose report validation failed: {error}\n")
    print(("Verified " if args.check else "Wrote ") + ", ".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
