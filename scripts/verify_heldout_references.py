#!/usr/bin/env python3
"""Replay source-only qualification of the twelve-candidate inventory."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.heldout_reference import main

if __name__ == "__main__":
    raise SystemExit(main())
