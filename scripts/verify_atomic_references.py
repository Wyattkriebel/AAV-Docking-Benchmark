#!/usr/bin/env python3
"""Replay frozen experimental atomic-reference sources without docking outcomes."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.atomic_reference import main


if __name__ == "__main__":
    raise SystemExit(main())
