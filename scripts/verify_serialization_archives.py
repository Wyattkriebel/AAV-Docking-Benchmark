#!/usr/bin/env python3
"""Verify the fixed metadata-only CIF amendment without native tools or network."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gemmi
from aav_benchmark.heldout_serialization import (
    CASE_IDS, SEQUENCE_SCOPE, atom_semantics, compare_atom_semantics, finalize_cif_metadata,
)

PINNED_FILES = {
    "aav_benchmark/heldout_serialization.py": "8e977d7576a60344b2793a8b805ab132f92e5c54494e6d3332af2b9ffba91d0d",
    "tests/test_heldout_serialization.py": "0e6da2f5c176e26a6493b539730ddeec30643df59a2d9e12bcc7b9c737615cc5",
    "results/structural_test/source_preparation.tar.gz": "88915e34c2ae8e5202831799dec3e53be1ad08af097c8cb0845d536cde29ef40",
    "results/structural_test/source_archive_manifest.json": "d658f6be95ba3c01ba0f6d354b69c280403168dc304170c0dc0fcde030f33077",
    "results/structural_test/source_preparation.json": "2444c9c6bdd953b5cdbd8050b3b7eb525c5e3b95d93e47f4138803c2a436631c",
    "results/structural_test_v2/source_preparation.tar.gz": "e913b9ca15acde3646c33307c188d19f9a44c8ee2439ed181fd17d876fa11116",
    "results/structural_test_v2/source_archive_manifest.json": "d35dc3b4433b234f9d2b95c8ee00b55bf450ec40d51c2680cf55927a82e21edd",
    "results/structural_test_v2/source_preparation.json": "645728a27e358ee42452cb765b12a48a4593a43d6377597dafb064592d2a149d",
    "results/structural_test_v2/serialization_amendment.json": "04e48516585be23de30de565d815d653088254c1b85bad73486090baa9a19556",
}
GEMMI_VERSION = "0.7.5"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def safe_name(name):
    if not isinstance(name, str) or not name or "\\" in name or "\0" in name:
        raise ValueError("Invalid archive member name")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != name or name == ".":
        raise ValueError("Archive member path must be canonical and relative")
    return name


def read_archive(raw, manifest, *, expected_count):
    """Return strictly validated member bytes; never extract onto the filesystem."""
    if digest(raw) != manifest["archive_sha256"]:
        raise ValueError("Preparation archive SHA256 changed")
    expected = {}
    for entry in manifest["files"]:
        name = safe_name(entry["path"])
        if (name in expected or type(entry["bytes"]) is not int or entry["bytes"] < 0
                or not isinstance(entry["sha256"], str) or len(entry["sha256"]) != 64
                or any(c not in "0123456789abcdef" for c in entry["sha256"])):
            raise ValueError("Invalid or duplicate archive manifest record")
        expected[name] = entry
    if len(expected) != expected_count:
        raise ValueError("Preparation archive member denominator changed")
    result = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        for member in archive:
            name = safe_name(member.name)
            if (member.type not in (tarfile.REGTYPE, tarfile.AREGTYPE) or member.linkname
                    or getattr(member, "sparse", None) is not None):
                raise ValueError("Only ordinary regular archive members are permitted")
            if name in result or name not in expected:
                raise ValueError("Duplicate or unexpected preparation archive member")
            entry = expected[name]
            if member.size != entry["bytes"]:
                raise ValueError("Preparation archive member size changed")
            content = archive.extractfile(member).read(entry["bytes"] + 1)
            if len(content) != entry["bytes"] or digest(content) != entry["sha256"]:
                raise ValueError("Preparation archive member SHA256 or size changed")
            result[name] = content
    if set(result) != set(expected):
        raise ValueError("Preparation archive is missing declared members")
    return result


def unique_rows(rows, key):
    result = {}
    for row in rows:
        value = row[key]
        if value in result:
            raise ValueError("Duplicate provenance record")
        result[value] = row
    return result


def verify_case_receipts(files, receipt):
    if receipt["case_ids"] != list(CASE_IDS):
        raise ValueError("Fixed six-case preparation cohort changed")
    cases = unique_rows(receipt["cases"], "case_id")
    if set(cases) != set(CASE_IDS):
        raise ValueError("Preparation receipt lost or added a case")
    for cid, case in cases.items():
        prefix = cid + "/"
        if json.loads(files[prefix + "source_preparation.json"]) != case:
            raise ValueError("Per-case and aggregate preparation receipts disagree")
        if (json.loads(files[prefix + "protein_preparation.json"]) != case["protein"]
                or json.loads(files[prefix + "receptor.cif.manifest.json"]) != case["protein"]
                or case["protein"]["receptor_sha256"] != digest(files[prefix + "receptor.cif"])
                or case["protein"]["atom_map_sha256"] != digest(files[prefix + "protein_atom_map.json"])):
            raise ValueError("Protein preparation receipts disagree with their case or archived atom inputs")
        expected = {name.removeprefix(prefix): digest(raw) for name, raw in files.items()
                    if name.startswith(prefix) and name != prefix + "source_preparation.json"}
        if case["artifact_sha256"] != expected:
            raise ValueError("Case artifact pins differ from its complete archived input set")
    return cases


def verify_metadata_pair(original, amended, original_receipt, amended_receipt, amendment, *, module_sha256):
    """Verify exact metadata changes, source decisions, and archived/fresh CIFs."""
    if gemmi.__version__ != GEMMI_VERSION:
        raise ValueError("Replay requires the repository-pinned gemmi version")
    if len(original) != 99 or len(amended) != 107:
        raise ValueError("Keep the original99 and amended107 preparation members")
    allowed = {"source_preparation.json"} | {f"{cid}/{name}" for cid in CASE_IDS for name in
        ("receptor.cif", "receptor.cif.manifest.json", "protein_preparation.json", "source_preparation.json")}
    added = {f"{cid}/serialization_v2.json" for cid in CASE_IDS} | {
        "configuration.proposal.json", "serialization_amendment.json"}
    if set(amended) != set(original) | added or set(original) & added:
        raise ValueError("Metadata amendment changed an unexpected source identity")
    changed = {p: {"path": p, "prior_sha256": digest(raw), "sha256": digest(amended[p])}
               for p, raw in original.items() if raw != amended[p]}
    unchanged = {p: {"path": p, "sha256": digest(raw)} for p, raw in original.items() if raw == amended[p]}
    if set(changed) != allowed:
        raise ValueError("Only the exact predeclared metadata files may change")
    if (unique_rows(amendment["changed_original_files"], "path") != changed
            or unique_rows(amendment["unchanged_original_files"], "path") != unchanged):
        raise ValueError("Amendment provenance does not enumerate the actual source changes")
    if (json.loads(original["source_preparation.json"]) != original_receipt
            or json.loads(amended["source_preparation.json"]) != amended_receipt
            or json.loads(amended["serialization_amendment.json"]) != amendment):
        raise ValueError("External and archived source/amendment receipts disagree")
    original_cases = verify_case_receipts(original, original_receipt)
    amended_cases = verify_case_receipts(amended, amended_receipt)
    added_receipt_fields = {"geometry_and_chemical_graphs_changed", "prior_scientific_inputs_and_failures_retained",
        "prior_source_preparation_sha256", "serialization_amendment_code_sha256", "serialization_revision"}
    if (set(amended_receipt) != set(original_receipt) | added_receipt_fields
            or {k: v for k, v in original_receipt.items() if k != "cases"}
            != {k: v for k, v in amended_receipt.items() if k != "cases" and k not in added_receipt_fields}
            or amended_receipt["serialization_revision"] != 2):
        raise ValueError("Aggregate source decisions or exact serialization receipt additions changed")
    flags = ("all_PDB_PDBQT_SDF_files_byte_identical", "all_author_atoms_elements_coordinates_exactly_unchanged",
             "native_ligand_graphs_and_controls_byte_identical", "prior_tree_exactly_unchanged")
    if (amendment["fixed_case_ids"] != list(CASE_IDS) or amendment["original_file_count"] != 99
            or amendment["amendment"] != "serialization_metadata_only" or amendment["native_ligand_coordinate_access"] is not False
            or any(amendment[k] is not True for k in flags)
            or amendment["prior_source_preparation_sha256"] != digest(original["source_preparation.json"])
            or amendment["source_preparation_v2_sha256"] != digest(amended["source_preparation.json"])
            or amendment["code_sha256"] != module_sha256
            or amended_receipt["serialization_amendment_code_sha256"] != module_sha256
            or amended_receipt["prior_source_preparation_sha256"] != digest(original["source_preparation.json"])
            or amended_receipt["prior_scientific_inputs_and_failures_retained"] is not True
            or amended_receipt["geometry_and_chemical_graphs_changed"] is not False):
        raise ValueError("Metadata-only amendment provenance or source pins changed")
    proofs = unique_rows(amendment["cases"], "case_id")
    if set(proofs) != set(CASE_IDS):
        raise ValueError("Amendment must retain exactly six atom-preservation proofs")
    protected = [p for p in original if PurePosixPath(p).suffix.lower() in (".pdb", ".pdbqt", ".sdf")]
    graphs = [f"{cid}/{name}" for cid in CASE_IDS for name in ("chemistry.json", "ligand_atom_map.json")]
    if any(original[p] != amended[p] for p in protected + graphs):
        raise ValueError("Metadata amendment changed PDB/PDBQT/SDF or a chemical source graph")
    atom_counts, label_counts = {}, {}
    with tempfile.TemporaryDirectory(prefix="cif-metadata-replay-") as temporary:
        temporary = Path(temporary)
        for cid in CASE_IDS:
            prior, actual, fresh = (temporary / f"{cid}_{suffix}.cif" for suffix in ("prior", "archived", "fresh"))
            prior.write_bytes(original[f"{cid}/receptor.cif"])
            actual.write_bytes(amended[f"{cid}/receptor.cif"])
            comparison = compare_atom_semantics(atom_semantics(prior), atom_semantics(actual))
            proof = json.loads(amended[f"{cid}/serialization_v2.json"])
            fresh_proof = finalize_cif_metadata(prior, fresh)
            if fresh.read_bytes() != amended[f"{cid}/receptor.cif"]:
                raise ValueError("Fresh finalization does not reproduce archived V2 CIF bytes")
            if ({k: v for k, v in proof.items() if k not in ("source_receptor", "receptor")}
                    != {k: v for k, v in fresh_proof.items() if k not in ("source_receptor", "receptor")}):
                raise ValueError("Archived atom-preservation proof differs from fresh replay")
            if (proof["serialization_sequence"] != SEQUENCE_SCOPE or proof["gemmi_version"] != GEMMI_VERSION
                    or proof["source_receptor_sha256"] != digest(original[f"{cid}/receptor.cif"])
                    or proof["receptor_sha256"] != digest(amended[f"{cid}/receptor.cif"])
                    or proofs[cid] != {"case_id": cid, **proof}
                    or amended_cases[cid]["serialization_revision"] != proof
                    or any(proof[k] != v for k, v in comparison.items())):
                raise ValueError("Per-case CIF proof changed its scope, atom values or input/output pins")
            old = {k: v for k, v in original_cases[cid].items() if k not in ("artifact_sha256", "protein")}
            new = {k: v for k, v in amended_cases[cid].items() if k not in ("artifact_sha256", "protein", "serialization_revision")}
            if old != new:
                raise ValueError("Metadata amendment changed a chemical/preparation source decision")
            expected_protein = {**original_cases[cid]["protein"],
                "receptor_sha256": digest(amended[f"{cid}/receptor.cif"]), "serialization_revision": 2,
                "superseded_input_sha256": digest(original[f"{cid}/receptor.cif"]),
                "serialization_revision_reason": proof["reason"],
                "serialization_revision_receipt": f"{cid}/serialization_v2.json",
                "serialization_revision_receipt_sha256": digest(amended[f"{cid}/serialization_v2.json"])}
            if amended_cases[cid]["protein"] != expected_protein:
                raise ValueError("Metadata amendment changed a protein source decision")
            atom_counts[cid], label_counts[cid] = comparison["atom_count"], comparison["atom_field_change_counts"]
    return {"archives": 2, "original_members": 99, "v2_members": 107, "cases": 6,
            "changed_metadata_members": len(changed), "unchanged_original_members": len(unchanged),
            "unchanged_PDB_PDBQT_SDF_members": len(protected), "unchanged_source_graph_members": len(graphs),
            "fresh_CIF_byte_matches": 6, "atom_counts": atom_counts, "label_field_change_counts": label_counts}


def verify(root):
    root = Path(root).resolve()
    files = {}
    for name, expected in PINNED_FILES.items():
        path = root / name
        if path.is_symlink() or not path.is_file() or any(p.is_symlink() for p in path.parents if p != root):
            raise ValueError(f"Pinned source must be an ordinary file: {name}")
        raw = path.read_bytes()
        if digest(raw) != expected:
            raise ValueError(f"Pinned source SHA256 changed: {name}")
        files[name] = raw
    archives, receipts = [], []
    for parent, count in (("results/structural_test", 99), ("results/structural_test_v2", 107)):
        manifest = json.loads(files[parent + "/source_archive_manifest.json"])
        if "archive" in manifest and manifest["archive"] != "source_preparation.tar.gz":
            raise ValueError("Preparation archive filename changed")
        archives.append(read_archive(files[parent + "/source_preparation.tar.gz"], manifest, expected_count=count))
        receipts.append(json.loads(files[parent + "/source_preparation.json"]))
    amendment = json.loads(files["results/structural_test_v2/serialization_amendment.json"])
    return verify_metadata_pair(*archives, *receipts, amendment,
        module_sha256=PINNED_FILES["aav_benchmark/heldout_serialization.py"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        result = verify(args.root)
    except (ValueError, KeyError, OSError, EOFError, tarfile.TarError, json.JSONDecodeError) as error:
        parser.exit(1, f"Serialization archive verification failed: {error}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
