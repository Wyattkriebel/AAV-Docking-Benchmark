"""Command line entry point; expensive work always has an explicit command."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

from .provenance import fetch_manifest, fetch_verified

ROOT = Path(__file__).resolve().parents[1]


def install_vina(destination, root=ROOT):
    system = {"Darwin": "mac", "Linux": "linux"}.get(platform.system())
    machine = {"arm64": "aarch64", "aarch64": "aarch64", "x86_64": "x86_64", "AMD64": "x86_64"}.get(platform.machine())
    key = f"{system}_{machine}"
    assets = json.loads((root / "benchmarks/vina.lock.json").read_text())["assets"]
    if key not in assets:
        raise ValueError(f"No frozen Vina executable for {key}; use a supported platform")
    asset = assets[key]
    fetch_verified(asset["url"], destination, asset["sha256"])
    Path(destination).chmod(0o755)
    return Path(destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Reproducible glycan docking diagnostics")
    parser.add_argument("--root", type=Path, default=ROOT, help="Repository root containing frozen benchmark ledgers")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("fetch", help="Fetch and verify frozen AAV source coordinates and CCD")
    tool = commands.add_parser("install-vina", help="Download the frozen official Vina 1.2.7 executable")
    tool.add_argument("--destination", type=Path)
    for name in ("prepare", "dock", "analyze", "replay", "archive", "restore", "merge"):
        command = commands.add_parser(name)
        command.add_argument("case", choices=["aav9", "154L", "1GNY", "1KJL", "1QFO", "2J1T", "2VUZ", "1PWB"])
        command.add_argument("--prepared", type=Path)
        command.add_argument("--output", type=Path)
        if name in ("dock", "replay"):
            command.add_argument("--vina", type=Path)
            command.add_argument("--resume", action="store_true")
        if name == "dock":
            command.add_argument("--only-seed", type=int, action="append", help="Execution shard; preserves the full frozen scientific denominator")
        if name == "merge":
            command.add_argument("--shard", type=Path, action="append", required=True)
            command.add_argument("--include-stopped-baseline", action="store_true", help="Preserve current runs.json after stopping its writer")
    status = commands.add_parser("status", help="Report evidence gates; fail if requested scientific closure is incomplete")
    status.add_argument("--require-complete", action="store_true")
    panel = commands.add_parser("panel", help="Prepare, dock or summarize the frozen exploratory Gal3 panel")
    panel.add_argument("stage", choices=["prepare", "dock", "analyze"])
    panel.add_argument("--resume", action="store_true")
    panel.add_argument("--prepared-root", type=Path)
    panel.add_argument("--results-root", type=Path)
    commands.add_parser("native-control", help="Score supplied AAV native coordinates; never count as docking recovery")
    commands.add_parser("report", help="Regenerate descriptive development tables from validated archived run outputs")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        if args.command == "fetch":
            fetch_manifest(root / "benchmarks/sources.lock.json", root / "data/raw")
        elif args.command == "install-vina":
            print(install_vina(args.destination or root / "tools/vina-1.2.7", root))
        elif args.command == "status":
            gates = json.loads((root / "benchmarks/closure_gates.json").read_text())
            print(json.dumps(gates, indent=2))
            if args.require_complete and not all(g["status"] == "passed" for g in gates["gates"]):
                return 2
        elif args.command == "panel":
            prepared_root = args.prepared_root or root / "data/processed"
            results_root = args.results_root or root / "results/development"
            if args.stage == "prepare":
                from .external import prepare_gal3_panel
                from .archive import restore
                if (root / "data/archive/external/inputs/manifest.json").is_file():
                    restore(root / "data/archive/external", root / "data/raw/external")
                prepare_gal3_panel(root / "data/raw/external", prepared_root / "1KJL", prepared_root)
            elif args.stage == "dock":
                from .engine import dock
                config = json.loads((root / "benchmarks/gal3_run.json").read_text())
                failed = False
                for case in config["cases"]:
                    runs = dock(prepared_root / case, results_root / case,
                                root / "tools/vina-1.2.7", root / "benchmarks/gal3_run.json", resume=args.resume)
                    failed = failed or any(r["status"] != "ok" for r in runs)
                if failed:
                    return 1
            else:
                from .panel import summarize_panel
                summarize_panel(prepared_root, results_root, root / "benchmarks/galectin3_itc_panel.json")
        elif args.command == "native-control":
            from .native_control import run_native_control
            run_native_control(root / "data/processed/aav9", root / "results/development/aav9/native_control",
                               root / "tools/vina-1.2.7", root / "results/development/aav9")
        elif args.command == "report":
            from .report import write_report
            print(write_report(root))
        else:
            from .engine import analyze, dock
            case_dir = args.case.lower() if args.case in ("aav9", "154L") else args.case
            prepared = args.prepared or root / "data/processed" / case_dir
            output = args.output or root / "results/development" / case_dir
            if args.command in ("archive", "restore"):
                from .archive import bundle, restore
                (bundle if args.command == "archive" else lambda p, o: restore(o, p))(prepared, output)
                return 0
            config_name = "development.json" if args.case in ("aav9", "154L") else "external_run.json"
            if args.command == "merge":
                from .shards import merge_shards
                baseline = json.loads((output / "runs.json").read_text()) if args.include_stopped_baseline else None
                merge_shards(prepared, output, root / "benchmarks" / config_name, args.shard, baseline)
                return 0
            if args.command in ("prepare", "replay"):
                if args.case == "aav9":
                    from .aav import prepare
                    fetch_manifest(root / "benchmarks/sources.lock.json", root / "data/raw")
                    prepare(root / "data/raw", prepared)
                else:
                    from .external import fetch_external, prepare
                    from .archive import restore
                    if (root / "data/archive/external/inputs/manifest.json").is_file():
                        restore(root / "data/archive/external", root / "data/raw/external")
                    fetch_external(root / "data/raw/external", args.case)
                    prepare(root / "data/raw/external", prepared, args.case)
            if args.command in ("dock", "replay"):
                vina = args.vina or root / "tools/vina-1.2.7"
                if not vina.exists():
                    install_vina(vina, root)
                runs = dock(prepared, output, vina, root / "benchmarks" / config_name, resume=args.resume,
                            seeds=args.only_seed if args.command == "dock" else None)
                if args.command == "dock" and any(r["status"] != "ok" for r in runs):
                    return 1
            if args.command in ("analyze", "replay"):
                summaries = analyze(prepared, output)
                print(json.dumps([{k: r.get(k) for k in ("case_id", "receptor", "seed", "status", "top1_rmsd_A", "best_returned_rmsd_A")} for r in summaries], indent=2))
                if any(r["status"] != "ok" for r in summaries):
                    return 1
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0
