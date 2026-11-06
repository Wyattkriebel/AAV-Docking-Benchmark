#!/usr/bin/env python3
"""Index, score, or verify the frozen unchanged retained-pose comparison."""
import argparse
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.common_pose_rescoring import index_source_pool, rescore_available, verify_common_scores
from aav_benchmark.provenance import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["index", "score", "watch", "verify"])
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--installation", type=Path)
    parser.add_argument("--max-poses", type=int, help="execution shard only; scientific denominator remains unchanged")
    parser.add_argument("--allow-pending", action="store_true")
    parser.add_argument("--interval", type=int, default=30)
    parser.add_argument("--verification-output", type=Path)
    args = parser.parse_args(argv)
    if args.stage in {"score", "watch"} and args.installation is None:
        parser.error("Native scoring requires --installation")
    if args.interval < 1 or args.interval > 60:
        parser.error("Watch interval must be between 1 and 60 seconds")
    try:
        if args.stage == "index":
            result = index_source_pool(args.root, output_root=args.output_root)
            print(f"Indexed {result['retained_pose_count']} retained poses; source terminal slots "
                  f"{result['source_global_denominator']['terminal_seed_runs']}/156; no accuracy computed")
        elif args.stage == "verify":
            result = verify_common_scores(args.root, output_root=args.output_root, require_complete=not args.allow_pending)
            if args.verification_output:
                write_json(args.verification_output, result)
            print(f"Verified {result['recorded_pose_count']} raw score records and "
                  f"{result['pose_identity_exports_reconstructed']} coordinate/identity exports without native programs")
        else:
            while True:
                result = rescore_available(args.root, installation=args.installation,
                                           output_root=args.output_root, max_poses=args.max_poses)
                print(f"Common pool: {result['recorded_pose_count']}/{result['indexed_pose_count']} recorded; "
                      f"{result['successful_pose_scores']} successful; source complete: "
                      f"{result['source_global_complete_and_verified']}; no accuracy computed", flush=True)
                if args.stage == "score" or result["ready_for_separate_accuracy_evaluation"]:
                    break
                time.sleep(args.interval)
        return 0
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as error:
        parser.exit(1, f"Frozen common-pool workflow failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
