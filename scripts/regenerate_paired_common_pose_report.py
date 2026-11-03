#!/usr/bin/env python3
"""Validate H2 readiness; the frozen global guard releases atomic reporting."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.paired_common_pose_report import example_plan, probe_workflow


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--example", action="store_true")
    args = parser.parse_args(argv)
    if args.example:
        result = example_plan(args.root)
    else:
        if args.plan is None or args.freeze is None:
            parser.error("--plan and --freeze are required for the validation-only probe")
        result = probe_workflow(args.root, args.plan, args.freeze)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
