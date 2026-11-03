"""Portable read-only replay of the archived guided HADDOCK example.

This rechecks retained evidence, not an installation or a new CNS execution.
The original release verifier and scientific helper remain unchanged.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile

from . import haddock_method as method


BASE = Path("results/haddock_reproduction")
LOCK = Path("benchmarks/haddock.lock.json")
LOCK_SHA256 = "5eec702d354597aa56b41a2fda30479bd8b481d35ee9498bb2ec0201f8f62bcb"
MANIFEST_SHA256 = "773e73c3ade68a67b354210ab16c26cd58c9112d5f627125bc12bbf2f98f4986"
HELPER_SHA256 = "04af6d150fcdb7195d81184dd7b3d7bdeb4e1649d509c2d3867295cf85b1558f"
CCD_SHA256 = {"NAG": "34acfbdbf655e8ee704112cb7ede8fca8489db9e5ae6211dbe1cc1469da726bb",
              "NDG": "e89ed48f0ded032c687b5ecbf7c88d0366e467f75e2d386448565b020d72d615"}
GEOMETRY_ROUNDOFF = 1e-9


def _require(condition, message):
    if not condition:
        raise method.HaddockAuditError(message)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load(path):
    return json.loads(Path(path).read_text())


def _local(root, relative):
    path = Path(relative)
    _require(not path.is_absolute() and ".." not in path.parts and path.parts,
             "Archive path must be local and relative")
    return root / path


def _pin(root, record):
    _require(isinstance(record, dict) and {"path", "sha256", "bytes"} <= set(record),
             "Incomplete retained evidence pin")
    path = _local(root, record["path"])
    _require(path.is_file() and _sha(path) == record["sha256"]
             and path.stat().st_size == record["bytes"],
             "Retained HADDOCK artifact checksum mismatch: " + record["path"])
    return path


def _without_historical_paths(value):
    """Remove only documented path roots from comparisons of derived reports."""
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in {"path", "pdb", "topology", "run_dir"} and isinstance(item, str):
                marker = str(BASE) + "/"
                if marker in item:
                    item = marker + item.split(marker, 1)[1]
            result[key] = _without_historical_paths(item)
        return result
    if isinstance(value, list):
        return [_without_historical_paths(item) for item in value]
    return value


def _same(actual, expected, label="report", field=None):
    if isinstance(actual, dict) and isinstance(expected, dict):
        _require(set(actual) == set(expected), label + " fields differ")
        for key in actual:
            _same(actual[key], expected[key], label + "/" + key, key)
    elif isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        _require(len(actual) == len(expected), label + " denominator differs")
        for index, (left, right) in enumerate(zip(actual, expected)):
            _same(left, right, label + "/" + str(index), field)
    elif field in {"coordinate_displacement_A", "ideal_signed_volume", "observed_signed_volume"}:
        _require(isinstance(actual, (int, float)) and isinstance(expected, (int, float))
                 and math.isfinite(actual) and math.isfinite(expected)
                 and math.isclose(actual, expected, rel_tol=0, abs_tol=GEOMETRY_ROUNDOFF),
                 label + " derived geometry differs")
    else:
        _require(type(actual) is type(expected) and actual == expected, label + " differs")


def _compare(actual, expected, label):
    _same(_without_historical_paths(actual), _without_historical_paths(expected), label)


def _retained_lock_pins(value):
    if isinstance(value, dict):
        if {"path", "sha256", "bytes"} <= set(value):
            yield value
        for item in value.values():
            yield from _retained_lock_pins(item)
    elif isinstance(value, list):
        for item in value:
            yield from _retained_lock_pins(item)


def _raw_final_metrics(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t"))
    ranked = sorted(rows, key=lambda row: int(row["caprieval_rank"]))
    _require(len(rows) == 12 and [int(row["caprieval_rank"]) for row in ranked] == list(range(1, 13))
             and len({row["model"] for row in rows}) == 12,
             "Raw final CAPRI model/rank denominator differs")
    for row in rows:
        _require(all(math.isfinite(float(row[key])) for key in
                     ["score", "irmsd", "fnat", "lrmsd", "ilrmsd", "dockq"]),
                 "Nonfinite upstream CAPRI metric")
    _require(all(float(a["score"]) <= float(b["score"]) for a, b in zip(ranked, ranked[1:])),
             "Upstream CAPRI ranks are not score ranks")
    best = max(rows, key=lambda row: float(row["dockq"]))
    fields = ["model", "caprieval_rank", "score", "irmsd", "fnat", "lrmsd", "ilrmsd", "dockq", "cluster_id"]
    return dict(final_model_count=12, final_cluster_count=len({row["cluster_id"] for row in rows}),
                score_ranked_first={key: ranked[0][key] for key in fields},
                best_sampled_DockQ={key: best[key] for key in fields})


def verify_haddock_replay(root=method.ROOT):
    root = Path(root).resolve()
    manifest_path = root / BASE / "artifact_manifest.json"
    _require(_sha(manifest_path) == MANIFEST_SHA256, "Frozen HADDOCK artifact manifest changed")
    _require(_sha(root / LOCK) == LOCK_SHA256, "Frozen HADDOCK release/input lock changed")
    _require(_sha(root / "aav_benchmark/haddock_method.py") == HELPER_SHA256
             and _sha(method.__file__) == HELPER_SHA256, "Frozen HADDOCK scientific helper changed")
    manifest, lock = _load(manifest_path), _load(root / LOCK)
    records = manifest["artifacts"]
    _require(manifest["lock_sha256"] == LOCK_SHA256 and manifest["artifact_count"] == len(records) == 423
             and len({row["path"] for row in records}) == 423
             and sum(row["bytes"] for row in records) == manifest["artifact_bytes"] == 10137416,
             "Frozen HADDOCK artifact denominator differs")
    for record in records:
        _require(Path(record["path"]).is_relative_to(BASE), "Artifact outside owned reproduction archive")
        _pin(root, record)
    retained, external = [], []
    for record in _retained_lock_pins(lock):
        if Path(record["path"]).parts[0] == "tools":
            external.append(record["path"])
        else:
            _pin(root, record)
            retained.append(record["path"])
    _require(len(set(retained)) == 9 and len(set(external)) == 5,
             "Retained versus external release-pin denominator differs")
    example = root / BASE / "official_glycan_test"
    run = example / "run1-test"
    upstream = (example / "upstream.cfg").read_text()
    _require((example / "reproduce.cfg").read_text() == upstream.replace("ncores = 10", "ncores = 1"),
             "Scientific example configuration changed beyond one CPU")
    source_example = root / BASE / "source_review/examples/docking-protein-glycan"
    _require((source_example / "docking-protein-glycan-test.cfg").read_bytes()
             == (example / "upstream.cfg").read_bytes(), "Retained official upstream configuration differs")
    for name in ["1LMQ_l_u.pdb", "1LMQ_r_u.pdb", "ambig.tbl", "target.pdb"]:
        _require((source_example / "data" / name).read_bytes() == (example / "data" / name).read_bytes(),
                 "Official example input differs: " + name)
    ccds = {name: root / f"data/archive/external/{name}.cif.gz" for name in CCD_SHA256}
    for name, path in ccds.items():
        _require(_sha(path) == CCD_SHA256[name], "Experimental CCD source checksum changed")
    summary = _load(root / BASE / "reproduction_summary.json")
    _require(summary["lock_verification"]["lock_sha256"] == LOCK_SHA256
             and summary["status"] == "completed_with_documented_upstream_input_rebuilding"
             and "not blind" in summary["workflow_condition"], "Recorded reproduction scope changed")
    inventory = method.workflow_inventory(run)
    _compare(inventory, summary["workflow"], "All fourteen stage inventories")
    stages = {row["stage"]: row for row in inventory["stages"]}
    _require(len(stages) == 14 and all(row["io_present"] and row["params_present"] for row in stages.values()),
             "Missing workflow stage I/O or parameters")
    for stage, expected in {"01_rigidbody": 20, "05_seletopclusts": 12,
                            "07_flexref": 12, "11_seletopclusts": 12}.items():
        io = _load(run / stage / "io.json")
        _require(stages[stage]["physical_pdb_files"] == len(io["output"]) == expected,
                 "Physical versus I/O model denominator differs: " + stage)
    top = run / "00_topoaa"
    prepared = top / "1LMQ_l_u_haddock.pdb.gz"
    topology_path = top / "1LMQ_l_u_haddock.psf.gz"
    inventories = {
        "input_receptor": example / "data/1LMQ_r_u.pdb", "input_glycan": example / "data/1LMQ_l_u.pdb",
        "prepared_receptor": top / "1LMQ_r_u_haddock.pdb.gz", "prepared_glycan": prepared,
    }
    for key, path in inventories.items():
        _compare(method.pdb_inventory(path), summary[key], key)
    for key, before, after in [
        ("preparation", inventories["input_glycan"], prepared),
        ("receptor_preparation", inventories["input_receptor"], inventories["prepared_receptor"]),
    ]:
        _compare(method.preparation_correspondence(before, after), summary[key], key)
    topology = method.read_cns_topology(topology_path)
    _same(topology, summary["actual_cns_glycan_topology"], "Typed CNS topology")
    _require(topology["atom_count"] == 111 and topology["bond_count"] == 114
             and abs(topology["total_charge"]) < 1e-12, "Typed glycan topology/charge differs")
    prepared_chirality = method.glycan_chirality(prepared, topology_path, ccds)
    _compare(prepared_chirality, summary["prepared_chirality"], "Prepared glycan handedness")
    chirality = []
    for stage in ["00_topoaa", "01_rigidbody", "07_flexref", "11_seletopclusts"]:
        paths = [prepared] if stage == "00_topoaa" else sorted((run / stage).glob("*.pdb.gz"))
        for path in paths:
            report = method.glycan_chirality(path, topology_path, ccds)
            _require(report["status"] == "consistent" and report["tested_centers"] == 20
                     and report["typed_heavy_atom_count"] == report["model_heavy_atom_count"] == 57,
                     "Model heavy-atom or handedness check failed")
            chirality.append(dict(stage=stage, path=str(path.relative_to(root)), pdb_sha256=_sha(path),
                status=report["status"], tested_centers=report["tested_centers"],
                typed_heavy_atom_count=report["typed_heavy_atom_count"],
                failed_centers=[row for row in report["centers"] if row["status"] != "consistent"]))
    _same(chirality, summary["all_model_chirality"], "Prepared plus 44 model handedness records")
    _require(len(chirality) == 45, "Prepared plus model handedness denominator differs")
    _compare(method.audit_air(example / "data/ambig.tbl"), summary["air"], "Native-informed AIRs")
    _require(summary["air"]["assignment_count"] == 15, "Known-interface AIR denominator differs")
    final = _raw_final_metrics(run / "13_caprieval/capri_ss.tsv")
    for key, value in final.items():
        _same(value, summary["released_example_self_evaluation"][key], "Upstream raw CAPRI " + key)
    log = (root / BASE / "official_workflow.stdout.log").read_text()
    times = re.findall(r"This HADDOCK3 run took: (.+)", log)
    _require("Finished at" in log and times == [summary["reported_runtime"]],
             "Recorded normal workflow completion/runtime differs")
    generated_archives = sorted(run.glob("analysis/*/summary.tgz"))
    _require(len(generated_archives) == 5, "Generated model archive denominator differs")
    for path in generated_archives:
        with tarfile.open(path) as archive:
            members = archive.getmembers()
            _require(members and all(row.isfile() and Path(row.name).name == row.name
                                    and row.name.endswith(".pdb") for row in members),
                     "Generated summary archive contains unexpected members")
    return dict(status="verified", method="HADDOCK3", version="2026.10.1",
                artifact_records_verified=423, artifact_manifest_sha256=MANIFEST_SHA256,
                retained_lock_paths_verified=9, external_release_pin_metadata_records=5,
                release_archives_and_CNS_installation_reverified=False, execution_performed=False,
                stage_count=14, rigidbody_models=20, refined_models=12, final_model_records=12,
                final_clusters=final["final_cluster_count"], prepared_plus_model_chirality_records=45,
                generated_model_archives=5, actual_glycan_atoms=111, actual_glycan_heavy_atoms=57,
                glycan_physical_bonds=114, known_interface_AIRs=15,
                derived_geometry_roundoff_tolerance=GEOMETRY_ROUNDOFF,
                tolerance_scope="Only semantic-alias distances and informative signed volumes; source bytes, identities, coordinates, scores, counts and handedness decisions remain exact.",
                released_example_self_evaluation=final,
                scope="Archived supplied-interface software example only; raw upstream CAPRI values are not qualified blind accuracy or an independently rerun docking result.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=method.ROOT)
    args = parser.parse_args(argv)
    try:
        result = verify_haddock_replay(args.root)
    except (ValueError, OSError, KeyError, TypeError, tarfile.TarError) as error:
        parser.exit(1, f"HADDOCK archive replay failed: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
