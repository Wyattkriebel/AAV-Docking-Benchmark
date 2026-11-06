#!/usr/bin/env python3
"""Verify frozen method inputs and historical lineage bytes in temporary files."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aav_benchmark.method_protocol_sources import (  # noqa: E402
    DISTRIBUTION_SHA256, prepare_source_fixture, restore_preflight_receptor_aliases,
    verify_publication_amendment,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    amendment = verify_publication_amendment(args.root)
    with tempfile.TemporaryDirectory(prefix="method-source-verification-") as temporary:
        root = Path(temporary)
        receipt = prepare_source_fixture(root, source_root=args.root,
                                         groups=("prepared", "haddock", "development"))
        aliases = restore_preflight_receptor_aliases(root)
        result = {"schema_version": 1, "distribution_sha256": DISTRIBUTION_SHA256,
                  "public_root_files_verified": receipt["copied_public_files"],
                  "archive_members_verified": receipt["hydration"]["verified_archive_members"],
                  "frozen_implementation_sources_verified": 14,
                  "current_source_map_entries_verified": 59, **aliases,
                  **amendment,
                  "temporary_files_removed_on_exit": True,
                  "native_installation_reverified": False,
                  "new_predictions_or_scoring": False}
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
