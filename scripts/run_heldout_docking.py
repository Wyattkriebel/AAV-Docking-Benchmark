#!/usr/bin/env python3
"""Prepare and execute the fixed six-case new structural docking experiment."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.heldout_docking import (
    planning_inventory, prepare_sources, finalize_preflight, native_ligand_control, run_predictor,
    prepare_condition, run_condition,
)
from aav_benchmark.provenance import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    commands = parser.add_subparsers(dest="command", required=True)
    inventory = commands.add_parser("inventory", help="source-only denominator and input pins")
    inventory.add_argument("--prepared", type=Path)
    inventory.add_argument("--output", type=Path)
    source = commands.add_parser("prepare-sources", help="mechanical protein and independent CCD ligand preparation only")
    source.add_argument("--output", type=Path, required=True)
    preflight = commands.add_parser("preflight", help="mechanical receptor parameterization, native stereo and source symmetry audits")
    preflight.add_argument("--prepared", type=Path, required=True)
    native = commands.add_parser("native-control", help="native-heavy stereo and serialization audit; no docking")
    native.add_argument("--prepared", type=Path, required=True)
    native.add_argument("--case", required=True)
    native.add_argument("--output", type=Path, required=True)
    predict = commands.add_parser("predict", help="execute a frozen full-assembly site predictor")
    predict.add_argument("--prepared", type=Path, required=True)
    predict.add_argument("--case", required=True)
    predict.add_argument("--method", choices=("fpocket", "p2rank"), required=True)
    predict.add_argument("--executable", type=Path, required=True)
    predict.add_argument("--output", type=Path, required=True)
    box = commands.add_parser("prepare", help="prepare fixed top-five or explicitly labeled oracle boxes")
    box.add_argument("--prepared", type=Path, required=True)
    box.add_argument("--case", required=True)
    box.add_argument("--method", choices=("fpocket", "p2rank", "oracle"), required=True)
    box.add_argument("--predictions", type=Path)
    box.add_argument("--oracle-center", nargs=3, type=float)
    box.add_argument("--output", type=Path, required=True)
    run = commands.add_parser("run", help="execute prepared boxes with the strict frozen engine")
    run.add_argument("--prepared", type=Path, required=True)
    run.add_argument("--source-prepared", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--vina", type=Path, required=True)
    run.add_argument("--resume", action="store_true")
    for command in (predict, box, run):
        command.add_argument("--config", type=Path, required=True)
        command.add_argument("--freeze-receipt", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "inventory":
            result = planning_inventory(args.root, args.prepared)
            if args.output:
                write_json(args.output, result)
        elif args.command == "prepare-sources":
            result = prepare_sources(args.root, args.output)
        elif args.command == "preflight":
            result = finalize_preflight(args.root, args.prepared)
        elif args.command == "native-control":
            result = native_ligand_control(args.root, args.case, args.prepared, args.output)
        elif args.command == "predict":
            result = run_predictor(args.root, args.prepared, args.case, args.method, args.output,
                                   args.config, args.freeze_receipt, executable=args.executable)
        elif args.command == "prepare":
            result = prepare_condition(args.root, args.prepared, args.case, args.method, args.output,
                                       args.config, args.freeze_receipt, predictions=args.predictions,
                                       oracle_center=args.oracle_center)
        else:
            result = run_condition(args.root, args.prepared, args.output, args.vina, args.config,
                                   args.freeze_receipt, args.source_prepared, resume=args.resume)
        print(f"{args.command}: {result.get('status', 'completed')}")
        if "cases" in result:
            for row in result["cases"]:
                print(f"  {row['case_id']}: {row.get('status', 'inventoried')}")
        return 0
    except (ValueError, OSError, RuntimeError, KeyError) as error:
        parser.exit(1, f"Heldout structural workflow failed: {error}\n")


if __name__ == "__main__":
    main()
