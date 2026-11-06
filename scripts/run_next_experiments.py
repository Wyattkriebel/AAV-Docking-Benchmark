#!/usr/bin/env python3
"""Replay the frozen site and Gal9C experiments in a fresh, separate directory.

A replay retains the original pre-run plan; it is not a newly blinded study.
Committed evidence is never used as the output of a scientific rerun.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aav_benchmark.provenance import sha256, write_json


def _load(path):
    return json.loads(Path(path).read_text())


def _outside_repository(root, destination):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    if destination == root or root in destination.parents:
        raise ValueError("Choose a replay directory outside the source repository")
    return destination


def _fresh_workspace(root, destination, arm):
    root = Path(root).resolve()
    destination = _outside_repository(root, destination)
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise FileExistsError("Replay preparation requires a new or empty directory")
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copytree(root / "benchmarks", destination / "benchmarks")
    (destination / "docs").mkdir()
    write_json(destination / "replay.json", {
        "schema_version": 1, "arm": arm, "source_repository": str(root),
        "role": "replay of a previously frozen development experiment; not newly prospective or held out",
        "configuration_sha256": sha256(root / "benchmarks" / f"{'site' if arm == 'site' else 'gal9'}_experiment.json"),
    })
    return destination


def _workspace(root, destination, arm):
    destination = _outside_repository(root, destination)
    marker = _load(destination / "replay.json")
    config = destination / "benchmarks" / f"{'site' if arm == 'site' else 'gal9'}_experiment.json"
    if marker.get("arm") != arm or marker.get("configuration_sha256") != sha256(config):
        raise ValueError("Replay workspace arm or frozen configuration differs from its preparation receipt")
    plan = _load(config)
    for relative, digest in plan.get("method_locks", {}).items():
        if sha256(destination / relative) != digest:
            raise ValueError("Replay method lock differs from the frozen configuration")
    if arm == "gal9" and "panel_source" in plan:
        for relative, digest in ((plan["panel_source"], plan["panel_sha256"]),
                                 ("benchmarks/gal9_receptor_audit.json", plan["receptor_audit_sha256"])):
            if sha256(destination / relative) != digest:
                raise ValueError("Replay assay panel or receptor audit differs from the frozen configuration")
    return destination


def _restore_source(archive, target, digest):
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(archive, "rb") as stream, target.open("wb") as output:
        shutil.copyfileobj(stream, output)
    if sha256(target) != digest:
        raise ValueError(f"Archived source checksum differs: {archive}")


def prepare_site(root, destination, *, reprepare=False):
    from aav_benchmark.site_prediction_verification import restore_site_inputs
    from aav_benchmark.sites import prepare_site_receptor, receptor_atoms

    root = Path(root).resolve()
    work = _fresh_workspace(root, destination, "site")
    config = _load(work / "benchmarks/site_experiment.json")
    source_dir = work / "data/raw"
    for source in _load(root / "benchmarks/sources.lock.json")["sources"]:
        _restore_source(root / source["archive"], source_dir / source["filename"], source["sha256"])
    for definition in _load(root / "benchmarks/external_cohort.json")["cases"]:
        for source in definition["sources"]:
            _restore_source(root / source["archive"], source_dir / "external" / source["filename"], source["sha256"])
    results = work / "results/site_discovery"
    results.mkdir(parents=True)
    # These are historical plan/revision receipts, copied unchanged. A replay
    # does not manufacture a new claim that its outcomes were unseen.
    for filename in ("freeze_receipt.json", "geometry_amendment_receipt.json", "fpocket_serialization_revision.json"):
        shutil.copyfile(root / "results/site_discovery" / filename, results / filename)
    preparation_checks = []
    for case in config["cases"]:
        name = case["condition_id"]
        if sha256(work / case["source"]) != case["source_sha256"]:
            raise ValueError(f"Frozen source changed for {name}")
        original = work / "data/processed/site_inputs" / name
        corrected = work / "data/processed/site_inputs_v2" / name
        restore_site_inputs(root, name, "original", original)
        restore_site_inputs(root, name, "fpocket_compatibility", corrected)
        # Relocate only private sidecar paths used by the runner. The original
        # byte archives remain in the source repository; receptor/map/receipt
        # bytes and coordinate identities stay unchanged.
        path = corrected / "receptor.cif.manifest.json"
        manifest = _load(path)
        manifest["receptor"] = str(corrected / "receptor.cif")
        manifest["atom_map"] = str(corrected / "receptor.cif.atoms.tsv.gz")
        manifest["serialization_revision_receipt"] = str(corrected / "receptor.cif.revision.json")
        write_json(path, manifest)
        if case["complex_id"] == "aav9":
            condition_results = results / name
            condition_results.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / "results/site_discovery" / name / "source_symmetry_amended.json",
                            condition_results / "source_symmetry_amended.json")
        if reprepare:
            rebuilt = work / "data/reprepared_sites" / name / "receptor.cif"
            prepare_site_receptor(work / case["source"], rebuilt, assembly_id=case["assembly_id"],
                                  lock_path=work / "benchmarks/site_methods.lock.json")
            a, b = receptor_atoms(original / "receptor.cif"), receptor_atoms(rebuilt)
            if set(a) != set(b) or any(not np.array_equal(a[key], b[key]) for key in a):
                raise ValueError(f"Fresh whole-receptor preparation differs from frozen atom data for {name}")
            preparation_checks.append({"condition_id": name, "author_atom_identities_and_coordinates_exact": True,
                                       "rebuilt_receptor_sha256": sha256(rebuilt),
                                       "note": "current serialization revision 2; archived original revision 1 is retained for exact input replay"})
    write_json(work / "preparation_checks.json", {"reprepare_requested": reprepare, "checks": preparation_checks})
    return work


def _surface(work, case, receptor, output):
    from aav_benchmark.site_geometry import capsid_symmetry, read_receptor
    from aav_benchmark.surface import build_surface_geometry, sample_surface_patches
    import aav_benchmark.surface as surface

    if output.exists():
        raise FileExistsError("Preserve the prior surface-baseline execution; use a new replay workspace")
    output.mkdir(parents=True)
    config_path = work / "benchmarks/site_experiment.json"
    config = _load(config_path)
    amendment = _load(work / "benchmarks/site_geometry_amendment.json")
    keys, xyz, elements = read_receptor(receptor)
    operations = []
    if case["complex_id"] == "aav9":
        operations, audit = capsid_symmetry(keys, xyz,
            tolerance_A=amendment["effective_operator_max_CA_displacement_A"],
            chain_fit_tolerance_A=config["source_symmetry"]["fit_max_CA_rmsd_A"])
        # Recompute the source-only certificate, then reuse the already frozen
        # operators after checking them. This avoids changing the historical
        # certificate hash through platform-level linear-algebra roundoff.
        frozen = _load(output.parent / "source_symmetry_amended.json")
        prior_operations = [(np.asarray(o["rotation"]), np.asarray(o["translation"])) for o in frozen["transforms"]]
        if len(operations) != len(prior_operations) or any(
                not np.allclose(r, fr, atol=1e-10, rtol=0) or not np.allclose(t, ft, atol=1e-10, rtol=0)
                for (r, t), (fr, ft) in zip(operations, prior_operations)):
            raise ValueError("Recomputed receptor-only symmetry differs from the frozen geometry")
        operations = prior_operations
        write_json(output / "source_symmetry.json", audit)
    baseline = config["baseline"]
    start = time.perf_counter()
    geometry = build_surface_geometry(xyz, elements, [tuple(k[:4]) for k in keys],
        points_per_atom=baseline["points_per_atom"], probe_radius_A=baseline["probe_radius_A"])
    geometry_time = time.perf_counter() - start
    trials = [sample_surface_patches(geometry, seed=seed, max_sites=max(config["candidate_budgets"]),
        patch_radius_A=baseline["patch_radius_A"], min_center_separation_A=baseline["minimum_center_separation_A"],
        symmetry_transforms=operations) for seed in baseline["seeds"]]
    result = {"condition_id": case["condition_id"], "method": "random_surface", "status": "ok",
        "configuration_sha256": sha256(config_path), "input_sha256": sha256(receptor),
        "native_information_used": False, "points_per_atom": baseline["points_per_atom"],
        "probe_radius_A": baseline["probe_radius_A"], "symmetry_operator_count": len(operations) or 1,
        "exposed_surface_points": geometry.exposed_point_count,
        "estimated_surface_area_A2": geometry.estimated_surface_area_A2,
        "geometry_seconds": geometry_time, "total_seconds": time.perf_counter() - start,
        "surface_code_sha256": sha256(Path(surface.__file__)), "replicates": trials}
    write_json(output / "result.json", result)
    return result


def run_site(root, destination, methods, *, fpocket=None, p2rank=None, java_home=None):
    from aav_benchmark.sites import run_fpocket
    from aav_benchmark.p2rank import run_p2rank
    work = _workspace(root, destination, "site")
    config = _load(work / "benchmarks/site_experiment.json")
    root = Path(root).resolve()
    fpocket_lock = _load(work / "benchmarks/site_methods.lock.json")["methods"]["fpocket"]
    p2_lock = _load(work / "benchmarks/p2rank.lock.json")
    fpocket = fpocket or root / "tools" / f"fpocket-{fpocket_lock['commit']}" / "bin/fpocket"
    p2rank = p2rank or root / p2_lock["executable"]
    java_home = java_home or root / p2_lock["runtime"]["java_home"]
    outcomes = []
    for case in config["cases"]:
        name = case["condition_id"]
        original = work / "data/processed/site_inputs" / name / "receptor.cif"
        corrected = work / "data/processed/site_inputs_v2" / name / "receptor.cif"
        directory = work / "results/site_discovery" / name
        for method in methods:
            if method == "fpocket":
                result = run_fpocket(corrected, directory / "fpocket_v2", executable=fpocket,
                    lock_path=work / "benchmarks/site_methods.lock.json")
            elif method == "p2rank":
                result = run_p2rank(original, directory / "p2rank", p2rank, java_home=java_home)
            else:
                result = _surface(work, case, original, directory / "random_surface_v2")
            outcomes.append({"condition_id": name, "method": method, "status": result["status"]})
            print(f"{name}/{method}: {result['status']}", flush=True)
    return outcomes


def prepare_gal9(root, destination):
    from aav_benchmark.archive import restore
    root = Path(root).resolve()
    work = _fresh_workspace(root, destination, "gal9")
    config = _load(work / "benchmarks/gal9_experiment.json")
    output = work / "results/gal9_specificity"
    output.mkdir(parents=True)
    shutil.copyfile(root / "results/gal9_specificity/freeze_receipt.json", output / "freeze_receipt.json")
    target = work / "data/archive/specificity"
    target.mkdir(parents=True)
    shutil.copyfile(root / "data/archive/specificity/3NV3.cif.gz", target / "3NV3.cif.gz")
    for case in config["cases"]:
        restore(root / "results/gal9_specificity" / case, work / "data/processed/gal9_panel" / case)
    return work


def run_gal9(root, destination, *, vina=None):
    from aav_benchmark.cli import install_vina
    from aav_benchmark.engine import dock
    work = _workspace(root, destination, "gal9")
    vina = Path(vina) if vina else Path(root).resolve() / "tools/vina-1.2.7"
    if not vina.exists():
        install_vina(vina, Path(root).resolve())
    allowed = {entry["sha256"] for entry in _load(work / "benchmarks/vina.lock.json")["assets"].values()}
    if vina.name != "vina-1.2.7" or sha256(vina) not in allowed:
        raise ValueError("Use the frozen official executable with the canonical vina-1.2.7 filename")
    config = work / "benchmarks/gal9_experiment.json"
    failed = False
    for case in _load(config)["cases"]:
        output = work / "results/gal9_specificity" / case
        if output.exists():
            raise FileExistsError("Preserve the prior docking execution; use a new replay workspace")
        runs = dock(work / "data/processed/gal9_panel" / case, output, vina, config)
        failed = failed or any(row["status"] != "ok" for row in runs)
        print(f"{case}: {sum(row['status'] == 'ok' for row in runs)}/{len(runs)} completed", flush=True)
    return failed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="source repository with retained archives and tools")
    commands = parser.add_subparsers(dest="stage", required=True)
    install = commands.add_parser("install-site-tools", help="install locked fpocket/P2Rank on supported macOS arm64")
    for name in ("site-prepare", "site-run", "site-evaluate", "gal9-prepare", "gal9-run", "gal9-report"):
        command = commands.add_parser(name)
        command.add_argument("--work-dir", type=Path, required=True, help="separate replay workspace outside the source repository")
        if name == "site-prepare":
            command.add_argument("--reprepare", action="store_true", help="also rebuild whole source receptors and verify exact coordinate/identity equivalence")
        if name == "site-run":
            command.add_argument("--method", choices=("fpocket", "p2rank", "random_surface"), action="append")
            command.add_argument("--fpocket", type=Path)
            command.add_argument("--p2rank", type=Path)
            command.add_argument("--java-home", type=Path)
        if name == "gal9-run":
            command.add_argument("--vina", type=Path)
    args = parser.parse_args(argv)
    try:
        root = args.root.resolve()
        if args.stage == "install-site-tools":
            from aav_benchmark.sites import install_fpocket
            from aav_benchmark.p2rank import install_p2rank
            print(install_fpocket(root / "tools")["executable"])
            print(install_p2rank(root)["executable"])
        elif args.stage == "site-prepare":
            print(prepare_site(root, args.work_dir, reprepare=args.reprepare))
        elif args.stage == "site-run":
            methods = args.method or ["fpocket", "p2rank", "random_surface"]
            if len(methods) != len(set(methods)):
                raise ValueError("Select each method at most once")
            run_site(root, args.work_dir, methods, fpocket=args.fpocket, p2rank=args.p2rank, java_home=args.java_home)
        elif args.stage == "site-evaluate":
            from aav_benchmark.site_experiment import evaluate_experiment
            work = _workspace(root, args.work_dir, "site")
            evaluate_experiment(work)
            print(work / "results/site_discovery/evaluation.json")
        elif args.stage == "gal9-prepare":
            print(prepare_gal9(root, args.work_dir))
        elif args.stage == "gal9-run":
            return int(run_gal9(root, args.work_dir, vina=args.vina))
        elif args.stage == "gal9-report":
            from aav_benchmark.specificity_report import write_specificity_report
            work = _workspace(root, args.work_dir, "gal9")
            print(write_specificity_report(work))
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"Replay failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
