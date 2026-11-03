#!/usr/bin/env python3
"""Export an exact Git commit and hydrate only explicit, verified archive members."""
import argparse
import json
from pathlib import Path
import sys
import tarfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.replay_workspace import restore_workspace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = restore_workspace(args.repository, args.plan, args.destination)
    except (ValueError, OSError, KeyError, EOFError, tarfile.TarError) as error:
        parser.exit(1, f"Replay workspace restoration failed: {error}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
