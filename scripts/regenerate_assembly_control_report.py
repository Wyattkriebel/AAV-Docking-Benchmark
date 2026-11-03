#!/usr/bin/env python3
"""Replay the prospective first-chain assembly-control report plan."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.assembly_control_report import example_plan, write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--allow-pending", action="store_true")
    parser.add_argument("--example", action="store_true")
    args = parser.parse_args()
    if args.example:
        print(json.dumps(example_plan(args.root), indent=2)); return
    if any(getattr(args, field) is None for field in ("plan", "freeze", "summary_output", "markdown_output")):
        parser.error("--plan, --freeze, --summary-output and --markdown-output are required")
    result = write_report(args.root, args.plan, args.freeze, args.summary_output, args.markdown_output,
                          check=args.check, allow_pending=args.allow_pending)
    print(json.dumps({"status": result["status"], "metrics_eligible": result["metrics_eligible"]}))


if __name__ == "__main__":
    main()
