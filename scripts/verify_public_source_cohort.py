#!/usr/bin/env python3
"""Verify public fixed-cohort references and family sources, without article downloads."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.family_audit import verify as verify_family
from aav_benchmark.public_source_replay import verify_public_source_references


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--density-tolerance-sigma", type=float, default=1e-5)
    parser.add_argument("--recompute-alignments", action="store_true")
    args = parser.parse_args()
    try:
        references = verify_public_source_references(args.root, density_tolerance_sigma=args.density_tolerance_sigma)
        families = verify_family(args.root, args.recompute_alignments)
    except (ValueError, KeyError, RuntimeError, OSError, ImportError) as error:
        parser.exit(1, f"Public source verification failed: {error}\n")
    print(json.dumps(dict(status="verified", references=references, family_audit=families), indent=2))


if __name__ == "__main__":
    main()
