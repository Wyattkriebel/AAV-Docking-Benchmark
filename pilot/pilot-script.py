"""Compatibility entry point for the repaired AAV9 diagnostic.

From the repository root use: python -m aav_benchmark replay aav9
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.cli import main

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python pilot/pilot-script.py prepare|dock|analyze [options]")
    raise SystemExit(main([sys.argv[1], "aav9", *sys.argv[2:]]))
