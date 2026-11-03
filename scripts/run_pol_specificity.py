#!/usr/bin/env python3
"""Prepare/preflight the source-pinned POL panel; docking needs a frozen config."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.pol_specificity import (
    ARCHIVE, LEDGER, dock_pol_panel, experiment_template, preflight, prepare_pol_panel,
)
from aav_benchmark.provenance import write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("template", "preflight", "prepare", "dock"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-ledger", type=Path, default=LEDGER)
    parser.add_argument("--source-dir", type=Path, default=ARCHIVE)
    parser.add_argument("--prepared-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--vina", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--case-id", action="append", dest="case_ids")
    parser.add_argument("--seed", action="append", type=int, dest="seeds")
    args = parser.parse_args()
    kwargs = {"ledger_path": args.source_ledger, "source_dir": args.source_dir}
    try:
        if args.action == "template":
            if args.config.exists():
                raise ValueError("Refusing to replace an existing experiment proposal/freeze")
            result = experiment_template(**kwargs)
            write_json(args.config, result)
        elif args.action == "preflight":
            result = preflight(args.config, **kwargs)
            if args.output:
                if args.output.exists():
                    raise ValueError("Preflight output must be fresh")
                write_json(args.output, result)
        elif args.action == "prepare":
            if args.prepared_root is None:
                parser.error("prepare requires --prepared-root")
            result = prepare_pol_panel(args.config, args.prepared_root, **kwargs)
        else:
            if args.prepared_root is None or args.output is None or args.vina is None:
                parser.error("dock requires --prepared-root, --output and --vina")
            result = dock_pol_panel(args.config, args.prepared_root, args.output, args.vina,
                                    resume=args.resume, case_ids=args.case_ids, seeds=args.seeds, **kwargs)
    except (OSError, ValueError, KeyError, RuntimeError, TypeError) as error:
        parser.exit(1, f"POL validation failed: {error}\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
