"""Receptor-only symmetry classes and evaluation-only site diagnostics."""

from __future__ import annotations

from pathlib import Path
import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .metrics import fit_receptor


def read_receptor(path):
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ValueError("Site experiment requires exactly one receptor model")
    keys, xyz, elements = [], [], []
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.element.is_hydrogen:
                    raise ValueError("Site experiment inputs must already be heavy-atom protein only")
                keys.append((chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name, atom.name))
                xyz.append(list(atom.pos))
                elements.append(atom.element.name)
    coords = np.asarray(xyz, dtype=float)
    if not keys or len(set(keys)) != len(keys) or not np.isfinite(coords).all():
        raise ValueError("Empty, duplicated or nonfinite site receptor")
    return keys, coords, elements


def capsid_symmetry(keys, coordinates, tolerance_A=0.05, chain_fit_tolerance_A=0.05):
    """Certify whole-capsid operators from protein C-alpha coordinates alone.

    Every transformed C-alpha must coincide with the same residue identity in
    another subunit. Merely fitting two homologous chains is insufficient.
    """
    ca = [(k, x) for k, x in zip(keys, coordinates) if k[-1] == "CA"]
    chain_ids = list(dict.fromkeys(k[0] for k, _ in ca))
    if len(chain_ids) != 60:
        raise ValueError("The audited capsid symmetry procedure requires 60 subunits")
    chains = {c: {k[1:4]: x for k, x in ca if k[0] == c} for c in chain_ids}
    identities = sorted(chains[chain_ids[0]])
    if any(set(chains[c]) != set(identities) for c in chain_ids):
        raise ValueError("Capsid subunit C-alpha identities differ")
    source = np.asarray([chains[chain_ids[0]][k] for k in identities])
    all_ca = np.asarray([x for _, x in ca])
    tree = cKDTree(all_ca)
    operations, audits = [], []
    for c in chain_ids:
        target = np.asarray([chains[c][k] for k in identities])
        rotation, translation = fit_receptor(source, target)
        fit_error = float(np.sqrt(np.mean(np.sum((source @ rotation + translation - target) ** 2, axis=1))))
        distances, nearest = tree.query(all_ca @ rotation + translation)
        identity_match = all(ca[i][0][1:4] == ca[int(j)][0][1:4] for i, j in enumerate(nearest))
        if fit_error > chain_fit_tolerance_A or not identity_match or len(set(nearest.tolist())) != len(all_ca):
            raise ValueError("Chain fit cannot propose a bijective residue-preserving capsid correspondence")
        # Fit all source-only correspondences, rather than amplify a tiny
        # single-subunit angular error across the diameter of the capsid.
        proposed_max = float(distances.max())
        rotation, translation = fit_receptor(all_ca, all_ca[nearest])
        distances, nearest = tree.query(all_ca @ rotation + translation)
        identity_match = all(ca[i][0][1:4] == ca[int(j)][0][1:4] for i, j in enumerate(nearest))
        if float(distances.max()) > tolerance_A or not identity_match or len(set(nearest.tolist())) != len(all_ca):
            raise ValueError("Proposed operator is not an audited whole-receptor symmetry")
        operations.append((rotation, translation))
        audits.append({"target_chain": c, "fit_CA_rmsd_A": fit_error,
                       "chain_proposal_whole_capsid_max_A": proposed_max,
                       "whole_capsid_max_CA_displacement_A": float(distances.max()),
                       "whole_capsid_refinement_CA_rmsd_A": float(np.sqrt(np.mean(distances ** 2))),
                       "whole_capsid_atom_correspondence_bijective": True,
                       "all_residue_identities_match": identity_match})
    return operations, {"source": "protein C-alpha chain fits propose identities, whole-capsid fits refine operators; no native ligand", "reference_chain": chain_ids[0],
                        "operator_count": len(operations), "tolerance_A": tolerance_A, "operators": audits,
                        "chain_fit_tolerance_A": chain_fit_tolerance_A,
                        "transforms": [{"rotation": r.tolist(), "translation": t.tolist()} for r, t in operations]}


def distinct_candidates(predictions, symmetry_transforms=None, merge_distance_A=4.0):
    """Preserve score order, collapsing only prespecified spatial orbit matches."""
    if not np.isfinite(merge_distance_A) or merge_distance_A <= 0:
        raise ValueError("Candidate merge distance must be positive")
    operations = symmetry_transforms or [(np.eye(3), np.zeros(3))]
    rows = sorted(predictions, key=lambda p: p["rank"])
    if [p["rank"] for p in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Candidate ranks must be complete and unique")
    unique, mapping = [], []
    for row in rows:
        center = np.asarray(row.get("center_A", row.get("center")), dtype=float)
        if center.shape != (3,) or not np.isfinite(center).all():
            raise ValueError("Candidate center is not a finite 3-vector")
        orbit = np.asarray([center @ np.asarray(r) + np.asarray(t) for r, t in operations])
        assigned = None
        for index, representative in enumerate(unique):
            if float(np.linalg.norm(orbit - representative["center_A"], axis=1).min()) <= merge_distance_A:
                assigned = index + 1
                break
        if assigned is None:
            assigned = len(unique) + 1
            unique.append({**row, "center_A": center.tolist(), "original_rank": row["rank"], "rank": assigned})
        mapping.append({"original_rank": row["rank"], "candidate_class": assigned})
    return unique, mapping


def evaluate_candidates(predictions, truth, keys, coordinates, *, rotation=None, translation=None,
                        thresholds=(3.0, 4.0, 5.0), footprint_radius_A=8.0, symmetry_transforms=None,
                        symmetry_tolerance_A=0.05):
    """Use native information only here, after candidate centers/ranks are fixed."""
    if not truth.get("evaluation_only") or not truth["sites"]:
        raise ValueError("A qualified evaluation-only reference is required")
    rotation = np.eye(3) if rotation is None else np.asarray(rotation)
    translation = np.zeros(3) if translation is None else np.asarray(translation)
    if (rotation.shape != (3, 3) or translation.shape != (3,) or not np.isfinite(rotation).all()
            or not np.isfinite(translation).all() or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8)
            or not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8)):
        raise ValueError("Evaluation requires a proper finite receptor transform")
    if not np.isfinite(footprint_radius_A) or footprint_radius_A <= 0:
        raise ValueError("Footprint radius must be positive")
    predictions = sorted(predictions, key=lambda row: row["rank"])
    if any(isinstance(p["rank"], bool) or not isinstance(p["rank"], int) for p in predictions) or [p["rank"] for p in predictions] != list(range(1, len(predictions) + 1)):
        raise ValueError("Evaluation candidate ranks must be complete and unique")
    for threshold in thresholds:
        if not np.isfinite(threshold) or threshold <= 0:
            raise ValueError("Site thresholds must be positive")
    tree = cKDTree(coordinates)
    residue_ids = [tuple(k[:4]) for k in keys]
    residue_universe = set(residue_ids)
    evaluated = []
    for row in predictions:
        original = np.asarray(row["center_A"], dtype=float)
        if original.shape != (3,) or not np.isfinite(original).all():
            raise ValueError("Candidate center must be a finite 3-vector")
        center = original @ rotation + translation
        distances = [float(np.linalg.norm(np.asarray(site["coordinates"]) - center, axis=1).min()) for site in truth["sites"]]
        # One consistent native copy supplies both distance and interface target.
        nearest = int(np.argmin(distances))
        native = {tuple(k) for k in truth["sites"][nearest]["interface_residue_keys"]}
        common = {residue_ids[i] for i in tree.query_ball_point(original, footprint_radius_A)}
        supplied = row.get("residue_keys", row.get("footprint_residues"))
        if supplied is None:
            supplied = [(r["chain"], r["number"], r["insertion_code"], r["resname"]) for r in row.get("residues", [])]
        supplied = {tuple(k) for k in supplied}
        if not supplied.issubset(residue_universe):
            raise ValueError("Method footprint contains a residue absent from its receptor input")
        def overlap(values):
            return {"recall": len(values & native) / len(native) if native else None,
                    "precision": len(values & native) / len(values) if values else 0.0}
        evaluated.append({**row, "evaluated_center_A": center.tolist(), "closest_native_site_id": truth["sites"][nearest]["id"],
                          "center_to_native_atom_A": distances[nearest], "native_interface_residue_count": len(native),
                          "common_radius_footprint": overlap(common), "method_footprint": overlap(supplied),
                          "normalized_method_footprint_residue_keys": [list(k) for k in sorted(supplied)],
                          "common_radius_residue_keys": [list(k) for k in sorted(common)]})
    result = {"returned_distinct_candidates": len(evaluated), "candidates": evaluated,
              "common_footprint_radius_A": footprint_radius_A,
              "coverage_definition": "Protein residue if any input heavy atom lies within the common radius; method footprint separately retained",
              "top1_distance_A": evaluated[0]["center_to_native_atom_A"] if evaluated else None,
              "top5_min_distance_A": min(r["center_to_native_atom_A"] for r in evaluated[:5]) if evaluated else None,
              "top1_common_interface_recall": evaluated[0]["common_radius_footprint"]["recall"] if evaluated else 0.0,
              "top1_method_interface_recall": evaluated[0]["method_footprint"]["recall"] if evaluated else 0.0}
    for t in thresholds:
        result[f"top1_success_{t:g}A"] = bool(evaluated and evaluated[0]["center_to_native_atom_A"] <= t)
        result[f"top5_success_{t:g}A"] = any(r["center_to_native_atom_A"] <= t for r in evaluated[:5])
    # Spatial classes may have representatives in different capsid copies.
    # Expand their source-certified receptor orbits before evaluating a union
    # against ONE native-copy denominator, keeping the metric copy-invariant.
    operations = symmetry_transforms or [(np.eye(3), np.zeros(3))]
    common_union, method_union = set(), set()
    ca = [(tuple(k[:4]), coordinates[i]) for i, k in enumerate(keys) if k[-1] == "CA"]
    ca_ids = [k for k, _ in ca]
    ca_map = dict(ca)
    ca_tree = cKDTree([x for _, x in ca]) if ca else None
    for row in evaluated[:5]:
        point = np.asarray(row["center_A"])
        for r, t in operations:
            moved = point @ np.asarray(r) + np.asarray(t)
            common_union.update(residue_ids[i] for i in tree.query_ball_point(moved, footprint_radius_A))
            for key in map(tuple, row["normalized_method_footprint_residue_keys"]):
                if key not in ca_map:
                    if len(operations) == 1:
                        method_union.add(key)
                        continue
                    raise ValueError("Source-symmetry footprint mapping requires a C-alpha for each method residue")
                distance, target = ca_tree.query(np.asarray(ca_map[key]) @ np.asarray(r) + np.asarray(t))
                if distance > symmetry_tolerance_A or ca_ids[int(target)][1:] != key[1:]:
                    raise ValueError("Method footprint symmetry does not preserve residue identity")
                method_union.add(ca_ids[int(target)])
    for label, union in [("common", common_union), ("method", method_union)]:
        scores = []
        for site in truth["sites"]:
            native = {tuple(k) for k in site.get("interface_residue_keys", [])}
            scores.append((len(union & native) / len(native) if native else 0.0, site.get("id"), len(native)))
        best = max(scores, key=lambda item: item[0])
        result[f"top5_{label}_union_interface_recall"] = best[0]
        result[f"top5_{label}_union_native_site_id"] = best[1]
        result[f"top5_{label}_union_native_denominator"] = best[2]
    result["top5_union_definition"] = "Union of source-symmetry-expanded Top5 candidate-class footprints, compared separately against each single native copy; maximum recall reported"
    return result
