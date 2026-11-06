#!/usr/bin/env python3
"""Prepare one frozen prospective strict/repair CCD-ideal candidate condition."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.heldout_pose_inputs import prepare_paired_condition


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--profile", choices=("strict", "repair"), required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--method", choices=("fpocket", "p2rank", "oracle"), required=True)
    parser.add_argument("--predictions", type=Path)
    parser.add_argument("--oracle-center", type=float, nargs=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--freeze-receipt", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        receipt = prepare_paired_condition(args.root, args.sources, args.profile, args.case, args.method, args.predictions,
                                           args.output, args.config, args.freeze_receipt, oracle_center=args.oracle_center)
        print(f"{args.profile}/{args.case}/{args.method}: {receipt['status']}; {len(receipt['jobs'])} retained candidate jobs")
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"Prospective pose preparation failed: {error}\n")


if __name__ == "__main__": main()
