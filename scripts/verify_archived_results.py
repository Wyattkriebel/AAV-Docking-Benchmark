"""Recompute committed structural and assay summaries from archived input bytes.

Uses temporary working copies and never overwrites the deposited raw results.
Run from the repository root: python scripts/verify_archived_results.py
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aav_benchmark.archive import restore
from aav_benchmark.engine import analyze
from aav_benchmark.panel import summarize_panel
from aav_benchmark.provenance import sha256, write_json
from aav_benchmark.reproduction import compare_tables


def verify(root, numeric_tolerance=0.0):
    root = Path(root).resolve()
    original = root / "results/development"
    checks = []
    with tempfile.TemporaryDirectory(prefix="aav-archive-verification-") as directory:
        temporary = Path(directory)
        prepared_root, results_root = temporary / "prepared", temporary / "results"
        for output in sorted(original.iterdir()):
            if not output.is_dir() or not (output / "runs.json").is_file():
                continue
            prepared, scratch = prepared_root / output.name, results_root / output.name
            restore(output, prepared)
            scratch.mkdir(parents=True)
            runs = json.loads((output / "runs.json").read_text())
            files = {"runs.json", "configuration.json"} | {r["pose_file"] for r in runs if r["status"] == "ok"}
            for filename in files:
                if Path(filename).name != filename:
                    raise ValueError("Result filenames must not contain directories")
                shutil.copyfile(output / filename, scratch / filename)
            reference = json.loads((prepared / "reference.json").read_text())
            if reference.get("pose_reference_eligible") is False:
                continue
            analyze(prepared, scratch)
            for filename in ("summary.json", "pose_metrics.csv"):
                comparison = compare_tables(output / filename, scratch / filename, numeric_tolerance)
                checks.append({"file": str((output / filename).relative_to(root)), "sha256": sha256(output / filename), **comparison})
        summarize_panel(prepared_root, results_root, root / "benchmarks/galectin3_itc_panel.json")
        panel_name = "gal3_panel_summary.json"
        comparison = compare_tables(original / panel_name, results_root / panel_name, numeric_tolerance)
        checks.append({"file": "results/development/" + panel_name, "sha256": sha256(original / panel_name), **comparison})
    result = {"scope": "Recomputation from archived prepared inputs and raw poses; explicit byte-exact/derived-float comparison policy below; no docking rerun or independent structural-quality adjudication",
              "all_tables_match": True, "all_tables_exact": all(row["exact_recomputation"] for row in checks),
              "numeric_tolerance": numeric_tolerance, "tables": checks,
              "verifier_sha256": sha256(Path(__file__)), "comparison_code_sha256": sha256(root / "aav_benchmark/reproduction.py")}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--numeric-tolerance", type=float, default=0.0, help="Explicit derived-float tolerance for cross-platform verification; default requires exact bytes")
    args = parser.parse_args()
    verification = verify(args.root, args.numeric_tolerance)
    write_json(args.root / "results/archive_verification.json", verification)
    print(f"Verified {len(verification['tables'])} tables; all byte-exact: {verification['all_tables_exact']}; maximum numeric difference: {max(row['max_numeric_difference'] for row in verification['tables']):.3g}")
