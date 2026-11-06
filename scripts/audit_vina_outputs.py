#!/usr/bin/env python3
"""Print a read-only I/O audit of existing main/control Vina 1.2.7 outputs."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.vina_output_audit import DEFAULT_SCOPES, audit_main_vina_outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--scope", action="append", choices=DEFAULT_SCOPES,
                        help="Repeat to restrict the audit; default is all named existing main/control scopes")
    args = parser.parse_args(argv)
    report = audit_main_vina_outputs(args.root, scopes=tuple(args.scope) if args.scope else DEFAULT_SCOPES)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "all_audited_outputs_writer_compatible" else 1


if __name__ == "__main__":
    raise SystemExit(main())
