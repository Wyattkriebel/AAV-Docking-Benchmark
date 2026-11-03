#!/usr/bin/env python3
"""Inspect or prepare the separately frozen standard-sidechain sensitivity."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.heldout_repair import inspect_sources, prepare_repaired_sources


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("command", choices=("inspect", "prepare"))
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tools-lock", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "inspect":
            result = inspect_sources(args.root, args.prepared, args.output, args.tools_lock)
        else:
            if args.config is None or args.freeze_receipt is None:
                raise ValueError("Atom generation requires an explicit secondary config and coordinator freeze receipt")
            result = prepare_repaired_sources(args.root, args.prepared, args.output, args.tools_lock, args.config, args.freeze_receipt)
        for row in result["cases"]: print(f"{row['case_id']}: {row['status']}")
        return 0
    except (ValueError, OSError, RuntimeError) as error:
        parser.exit(1, f"Heldout repair sensitivity failed: {error}\n")


if __name__ == "__main__": main()
