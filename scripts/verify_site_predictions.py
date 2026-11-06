#!/usr/bin/env python3
"""Verify archived site inputs and raw predictor outputs without native tools."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.site_prediction_verification import main


if __name__ == "__main__":
    raise SystemExit(main())
