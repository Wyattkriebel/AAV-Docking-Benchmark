#!/usr/bin/env python3
"""Bind released structural reports or aggregate frozen indicator contrasts."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.structural_contrast_report import bind_released_reports, example_plan, generate_report, _write_new


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("example", "bind", "aggregate"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--binding-sha256")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.stage == "example":
        result = example_plan(args.root)
    else:
        if args.plan is None or args.freeze is None:
            parser.error("--plan and --freeze are required")
        if args.stage == "bind":
            result = bind_released_reports(args.root, args.plan, args.freeze)
        else:
            if args.binding_sha256 is None or args.output is None:
                parser.error("Aggregation requires --binding-sha256 and a fresh --output")
            result = generate_report(args.root, args.plan, args.freeze, binding_sha256=args.binding_sha256)
            _write_new(args.root.resolve(), args.output, result)
            result = {"status": result["status"], "main_H1_contrasts": len(result["main_H1_guidance_contrasts"]),
                      "preparation_contrasts": len(result["preparation_sensitivity_contrasts"])}
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
