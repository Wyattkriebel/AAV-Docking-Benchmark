#!/usr/bin/env python3
"""Index or explicitly frozen-score the complete paired Vina/VC union pool."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.paired_common_pose_rescoring import index_union_pool, score_available_union, verify_union_scores
from aav_benchmark.provenance import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["index", "score", "verify"])
    parser.add_argument("--plan", type=Path, default=Path("benchmarks/structural_pose_v3_common_pool.proposal.json"))
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--installation", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--max-poses", type=int)
    parser.add_argument("--allow-pending", action="store_true")
    parser.add_argument("--verification-output", type=Path)
    args = parser.parse_args(argv)
    if args.stage in {"score", "verify"} and args.receipt is None:
        parser.error("A separate root fixed-scoring freeze receipt is required")
    if args.stage == "score" and args.installation is None:
        parser.error("Native scoring requires an installation")
    if args.stage == "index":
        result = index_union_pool(plan_path=args.plan, output_root=args.output_root)
        print(f"Union index: {result['retained_pose_count']} distinct poses; {result['terminal_source_method_slots']}/216 terminal source slots; no fixed scores or accuracy computed")
    elif args.stage == "score":
        result = score_available_union(plan_path=args.plan, receipt_path=args.receipt, installation=args.installation,
                                       output_root=args.output_root, max_poses=args.max_poses)
        print(f"Union scores: {result['recorded_pose_count']}/{result['indexed_pose_count']}; no accuracy computed")
    else:
        result = verify_union_scores(plan_path=args.plan, receipt_path=args.receipt, output_root=args.output_root,
                                     require_complete=not args.allow_pending)
        if args.verification_output:
            write_json(args.verification_output, result)
        print(f"Portable union verification: {result['recorded_pose_count']} score records; no native executable or accuracy computation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
