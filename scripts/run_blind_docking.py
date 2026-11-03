#!/usr/bin/env python3
"""Prepare and execute the fixed predicted-box/oracle development comparison."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.blind_docking import (
    planning_inventory, prepare_predicted_condition, prepare_oracle_condition, run_condition,
)
from aav_benchmark.provenance import sha256, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    commands = parser.add_subparsers(dest="command", required=True)
    inventory = commands.add_parser("inventory", help="source-only candidate/chemical-input inventory")
    inventory.add_argument("--output", type=Path)
    prepare = commands.add_parser("prepare", help="prepare one condition; existing directories are rejected")
    prepare.add_argument("--condition", required=True)
    prepare.add_argument("--method", choices=["fpocket", "p2rank", "oracle"], required=True)
    prepare.add_argument("--center", nargs=3, type=float, help="native-informed center; accepted only for the labeled oracle control")
    prepare.add_argument("--config", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--receptor-cache", type=Path)
    run = commands.add_parser("run", help="execute all prepared boxes, retaining failures")
    run.add_argument("--prepared", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--freeze-receipt", type=Path, required=True,
                     help="pre-run receipt containing configuration_sha256 and frozen_before_new_docking:true")
    run.add_argument("--vina", type=Path, required=True)
    run.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "inventory":
            result = planning_inventory(args.root)
            if args.output:
                write_json(args.output, result)
            print(f"Inventoried {result['eligible_method_conditions']} eligible method conditions")
        elif args.command == "prepare":
            if args.method == "oracle":
                if args.center is None:
                    raise ValueError("The oracle control requires an explicitly supplied center")
                result = prepare_oracle_condition(args.root, args.condition, args.center, args.output, args.config,
                                                  receptor_cache_dir=args.receptor_cache)
            else:
                if args.center is not None:
                    raise ValueError("A native-informed center is forbidden for predicted-site input")
                result = prepare_predicted_condition(args.root, args.condition, args.method, args.output, args.config,
                                                     receptor_cache_dir=args.receptor_cache)
            print(f"{result['condition_id']}/{result['method']}: {result['status']}; {len(result['jobs'])} retained candidate jobs")
        else:
            freeze = json.loads(args.freeze_receipt.read_text())
            if (freeze.get("configuration_sha256") != sha256(args.config)
                    or freeze.get("frozen_before_new_docking") is not True):
                raise ValueError("Pre-execution freeze receipt does not bind this configuration")
            lock = json.loads((args.root / "benchmarks/vina.lock.json").read_text())
            if sha256(args.vina) not in {asset["sha256"] for asset in lock["assets"].values()}:
                raise ValueError("Executable is not a frozen official Vina binary")
            result = run_condition(args.prepared, args.output, args.vina, args.config, resume=args.resume)
            print(f"{result['condition_id']}/{result['method']}: {result['status']}")
        return 0
    except (ValueError, OSError, RuntimeError, KeyError) as error:
        parser.exit(1, f"Blind-docking workflow failed: {error}\n")


if __name__ == "__main__":
    main()
