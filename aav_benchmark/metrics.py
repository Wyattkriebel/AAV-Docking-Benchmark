"""Strict, coordinate-frame-aware metrics for retained docking poses.

All coordinates are heavy atoms in Angstrom. The caller must establish chemical
identity, stereochemistry and atom ordering before invoking coordinate metrics.
An RMSD is a pose error: only a receptor fit may remove rigid-body motion.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np


def _coordinates(value: Any, name: str, *, allow_empty: bool = False) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if allow_empty and result.shape in ((0,), (0, 3)):
        return np.empty((0, 3), dtype=float)
    if result.ndim != 2 or result.shape[1] != 3 or not len(result):
        raise ValueError(f"{name} must be a nonempty N x 3 coordinate array")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} contains nonfinite coordinates")
    return result


def mapped_coordinates(
    names: Sequence[str], coordinates: Any, required_names: Sequence[str]
) -> np.ndarray:
    """Order atoms by an exact identity map, rejecting loss or duplicate names.

    Filter hydrogens explicitly before this function. Extra atoms are rejected
    as well as missing atoms; silently intersecting atom sets biases RMSD.
    """
    names = list(names)
    required_names = list(required_names)
    coordinates = _coordinates(coordinates, "coordinates")
    if len(names) != len(coordinates):
        raise ValueError("atom-name and coordinate counts differ")
    if len(names) != len(set(names)) or len(required_names) != len(set(required_names)):
        raise ValueError("duplicate atom identity in map")
    if not required_names:
        raise ValueError("required atom map is empty")
    missing = sorted(set(required_names) - set(names))
    extra = sorted(set(names) - set(required_names))
    if missing or extra:
        raise ValueError(f"atom identity mismatch: missing={missing}; extra={extra}")
    index = {name: i for i, name in enumerate(names)}
    return coordinates[[index[name] for name in required_names]].copy()


def fit_receptor(source: Any, target: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return a proper least-squares receptor fit using ``x @ rotation + shift``.

    Arrays must already have the same atom identity ordering. At least three
    noncollinear matched receptor atoms are required to determine a stable fit.
    Reflections are forbidden.
    """
    source = _coordinates(source, "source receptor")
    target = _coordinates(target, "target receptor")
    if source.shape != target.shape:
        raise ValueError("matched receptor coordinate counts differ")
    a = source - source.mean(axis=0)
    b = target - target.mean(axis=0)
    if len(source) < 3 or np.linalg.matrix_rank(a) < 2 or np.linalg.matrix_rank(b) < 2:
        raise ValueError("receptor fit requires three noncollinear matched atoms")
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if np.linalg.det(u @ vt) > 0 else -1.0
    rotation = u @ correction @ vt
    return rotation, target.mean(axis=0) - source.mean(axis=0) @ rotation


def _transform(value: Any) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError("symmetry transform must contain rotation and translation")
    rotation, translation = (np.asarray(v, dtype=float) for v in value)
    if rotation.shape != (3, 3) or translation.shape != (3,):
        raise ValueError("symmetry transform has invalid dimensions")
    if not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        raise ValueError("symmetry transform contains nonfinite values")
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6, rtol=0):
        raise ValueError("symmetry rotation is not orthonormal")
    if not np.isclose(np.linalg.det(rotation), 1, atol=1e-6, rtol=0):
        raise ValueError("symmetry rotation must be proper, not a reflection")
    return rotation, translation


def ligand_rmsd(
    pose: Any,
    reference: Any,
    *,
    receptor_source: Any = None,
    receptor_target: Any = None,
    permutations: Iterable[Sequence[int]] | None = None,
    symmetry_transforms: Iterable[tuple[Any, Any]] | None = None,
) -> float:
    """Calculate ligand RMSD after an optional *receptor-only* rigid fit.

    ``permutations`` must be graph- and stereochemistry-valid atom automorphisms
    established by the caller, never arbitrary spatial nearest-neighbor maps.
    ``symmetry_transforms`` must be independently established biological assembly
    operators in the reference frame. This function validates their rigid-body
    form, but cannot certify that an operator is a biological symmetry. Identity
    is always included. No ligand superposition or symmetry inference is done.
    """
    pose = _coordinates(pose, "pose")
    reference = _coordinates(reference, "reference")
    if pose.shape != reference.shape:
        raise ValueError("ligand atom counts differ; refusing partial-atom RMSD")
    if (receptor_source is None) != (receptor_target is None):
        raise ValueError("both receptor coordinate arrays are required for alignment")
    if receptor_source is not None:
        rotation, translation = fit_receptor(receptor_source, receptor_target)
        pose = pose @ rotation + translation
    maps = [np.arange(len(pose))]
    if permutations is not None:
        for permutation in permutations:
            p = np.asarray(permutation)
            if p.shape != (len(pose),) or p.dtype.kind not in "iu":
                raise ValueError("atom permutation must contain integer indices")
            if sorted(p.tolist()) != list(range(len(pose))):
                raise ValueError("atom permutation is not a complete bijection")
            maps.append(p)
    transforms = [(np.eye(3), np.zeros(3))]
    if symmetry_transforms is not None:
        transforms.extend(_transform(value) for value in symmetry_transforms)
    return min(
        float(np.sqrt(np.mean(np.sum((pose[p] @ r + t - reference) ** 2, axis=1))))
        for p in maps
        for r, t in transforms
    )


def contact_pairs(ligand: Any, receptor: Any, cutoff_A: float = 4.0) -> set[tuple[int, int]]:
    """Return heavy-atom index pairs at strictly less than the given distance.

    Atom identity/order must be consistent for every pose and its holo native
    contacts, including when the evaluated receptor is an aligned apo structure.
    """
    ligand = _coordinates(ligand, "ligand")
    receptor = _coordinates(receptor, "receptor", allow_empty=True)
    if not np.isfinite(cutoff_A) or cutoff_A <= 0:
        raise ValueError("contact cutoff must be finite and positive")
    distances_squared = np.sum((ligand[:, None, :] - receptor[None, :, :]) ** 2, axis=2)
    return set(map(tuple, np.argwhere(distances_squared < cutoff_A**2).tolist()))


def contact_metrics(
    pose: Any,
    receptor: Any,
    native_contacts: Iterable[tuple[int, int]],
    *,
    cutoff_A: float = 4.0,
    receptor_residues: Sequence[Any] | None = None,
    receptor_native_indices: Sequence[int] | None = None,
    native_receptor_count: int | None = None,
) -> dict[str, int | float | None]:
    """Compare pose contacts with a fixed, holo-derived native contact set.

    No native contacts means recall is undefined (``None``). No predicted
    contacts gives precision 0, retaining a failed-contact prediction visibly.
    Residue recall measures interface coverage; residue precision and Jaccard
    guard against large, nonspecific interfaces. For an apo receptor missing
    holo atoms, provide ``receptor_native_indices`` mapping each observed apo
    atom to the *full holo* atom universe, and ``native_receptor_count``. Native
    contacts remain in that full universe, so missing atoms reduce recall.
    ``receptor_residues``, if provided, indexes the full native universe.
    """
    ligand = _coordinates(pose, "pose")
    receptor = _coordinates(receptor, "receptor", allow_empty=True)
    if receptor_native_indices is None:
        if native_receptor_count is not None and native_receptor_count != len(receptor):
            raise ValueError("missing receptor atoms require an explicit native atom map")
        native_receptor_count = len(receptor)
        receptor_native_indices = list(range(len(receptor)))
    else:
        receptor_native_indices = list(receptor_native_indices)
        if not isinstance(native_receptor_count, (int, np.integer)) or native_receptor_count < 0:
            raise ValueError("apo atom mapping requires the full native receptor atom count")
        if len(receptor_native_indices) != len(receptor) or len(set(receptor_native_indices)) != len(receptor_native_indices):
            raise ValueError("apo-to-native receptor atom map has wrong count or duplicates")
        if any(not isinstance(i, (int, np.integer)) or not 0 <= i < native_receptor_count for i in receptor_native_indices):
            raise ValueError("apo-to-native receptor atom index outside the native universe")
    native = set()
    for value in native_contacts:
        if len(value) != 2 or any(not isinstance(v, (int, np.integer)) for v in value):
            raise ValueError("native contact indices must be integer pairs")
        i, j = value
        if not (0 <= i < len(ligand) and 0 <= j < native_receptor_count):
            raise ValueError("native contact index outside the fixed atom map")
        native.add((int(i), int(j)))
    predicted = {(i, receptor_native_indices[j]) for i, j in contact_pairs(ligand, receptor, cutoff_A)}
    shared = len(native & predicted)
    result: dict[str, int | float | None] = {
        "native_contact_count": len(native),
        "predicted_contact_count": len(predicted),
        "shared_contact_count": shared,
        "contact_recall": shared / len(native) if native else None,
        "contact_precision": shared / len(predicted) if predicted else 0.0,
    }
    if receptor_residues is not None:
        if len(receptor_residues) != native_receptor_count:
            raise ValueError("residue identity count differs from full native receptor atom count")
        native_residues = {receptor_residues[j] for _, j in native}
        predicted_residues = {receptor_residues[j] for _, j in predicted}
        intersection = len(native_residues & predicted_residues)
        union = len(native_residues | predicted_residues)
        result.update(
            native_interface_residue_count=len(native_residues),
            predicted_interface_residue_count=len(predicted_residues),
            interface_residue_recall=intersection / len(native_residues) if native_residues else None,
            interface_residue_precision=intersection / len(predicted_residues) if predicted_residues else 0.0,
            interface_residue_jaccard=intersection / union if union else None,
        )
    return result


def clash_metrics(
    ligand: Any,
    receptor: Any,
    ligand_radii: Sequence[float],
    receptor_radii: Sequence[float],
    *,
    overlap_tolerance_A: float = 0.4,
) -> dict[str, float | int | None]:
    """Report inter-molecular van der Waals overlaps using explicit atom radii.

    This is a geometric diagnostic, not a complete physical-validity test. The
    caller records the radius source and treats coordinated metals separately.
    A pair clashes if ``r_ligand + r_receptor - distance > tolerance``.
    """
    ligand = _coordinates(ligand, "ligand")
    receptor = _coordinates(receptor, "receptor", allow_empty=True)
    lr = np.asarray(ligand_radii, dtype=float)
    rr = np.asarray(receptor_radii, dtype=float)
    if lr.shape != (len(ligand),) or rr.shape != (len(receptor),):
        raise ValueError("radius and atom counts differ")
    if not np.isfinite(lr).all() or not np.isfinite(rr).all() or np.any(lr <= 0) or np.any(rr <= 0):
        raise ValueError("van der Waals radii must be finite and positive")
    if not np.isfinite(overlap_tolerance_A) or overlap_tolerance_A < 0:
        raise ValueError("overlap tolerance must be finite and nonnegative")
    if not len(receptor):
        return {"clash_pair_count": 0, "clashing_ligand_atom_count": 0, "max_overlap_A": None, "min_distance_A": None}
    distances = np.linalg.norm(ligand[:, None, :] - receptor[None, :, :], axis=2)
    overlaps = lr[:, None] + rr[None, :] - distances
    clashes = overlaps > overlap_tolerance_A
    return {
        "clash_pair_count": int(clashes.sum()),
        "clashing_ligand_atom_count": int(clashes.any(axis=1).sum()),
        "max_overlap_A": float(max(0.0, overlaps.max())),
        "min_distance_A": float(distances.min()),
    }


def summarize_poses(
    rows: Iterable[Mapping[str, Any]], thresholds: Sequence[float] = (1.0, 2.0, 3.0)
) -> dict[str, Any]:
    """Summarize retained poses by explicit score rank, including empty failures.

    Rank 1 must be the method's best-ranked pose, with contiguous unique ranks.
    Missing near-native retained poses cannot diagnose internal sampling failure.
    Conditional Top-1 success is undefined when no retained near-native exists.
    """
    rows = sorted(list(rows), key=lambda row: row["rank"])
    ranks = [row["rank"] for row in rows]
    if any(isinstance(rank, (bool, np.bool_)) or not isinstance(rank, (int, np.integer)) for rank in ranks) or ranks != list(range(1, len(rows) + 1)):
        raise ValueError("pose ranks must be unique contiguous integers starting at 1")
    rmsds = [float(row["rmsd_A"]) for row in rows]
    if any(not np.isfinite(value) or value < 0 for value in rmsds):
        raise ValueError("pose RMSDs must be finite and nonnegative")
    thresholds = tuple(float(threshold) for threshold in thresholds)
    if any(not np.isfinite(value) or value <= 0 for value in thresholds):
        raise ValueError("success thresholds must be finite and positive")
    result: dict[str, Any] = {
        "returned_poses": len(rows),
        "top1_rmsd_A": rmsds[0] if rows else None,
        "top5_min_rmsd_A": min(rmsds[:5]) if rows else None,
        "best_returned_rmsd_A": min(rmsds) if rows else None,
        "best_returned_rank": rmsds.index(min(rmsds)) + 1 if rows else None,
    }
    for threshold in thresholds:
        suffix = f"{threshold:g}A"
        retained = bool(rows and min(rmsds) <= threshold)
        top1 = bool(rows and rmsds[0] <= threshold)
        result[f"top1_success_{suffix}"] = top1
        result[f"top5_success_{suffix}"] = bool(rows and min(rmsds[:5]) <= threshold)
        result[f"retained_success_{suffix}"] = retained
        result[f"top1_success_given_retained_{suffix}"] = top1 if retained else None
    return result


def site_metrics(
    ranked_centers: Any,
    native_ligand: Any,
    *,
    equivalent_classes: Sequence[Any] | None = None,
    symmetry_transforms: Iterable[tuple[Any, Any]] | None = None,
    thresholds: Sequence[float] = (3.0, 4.0, 5.0),
) -> dict[str, Any]:
    """Measure ranked center recovery, merging caller-curated equivalent sites.

    The center endpoint alone does not measure whole-interface coverage; report
    residue coverage separately for extended glycans. The supplied centers are
    already ranked, and native coordinates must only be used at evaluation time.
    For equivalent sites, supply verified biological assembly operators to
    generate equivalent native coordinates, or include those coordinates in
    ``native_ligand``. Deduplicating prediction classes alone cannot make a
    distant equivalent site match one reference. Operators follow ``x @ R+t``.
    """
    centers = _coordinates(ranked_centers, "ranked centers", allow_empty=True)
    native = _coordinates(native_ligand, "native ligand")
    if symmetry_transforms is not None:
        operators = [_transform(value) for value in symmetry_transforms]
        native = np.vstack([native] + [native @ rotation + translation for rotation, translation in operators])
    if equivalent_classes is not None:
        if len(equivalent_classes) != len(centers):
            raise ValueError("site-class count differs from center count")
        seen = set()
        indices = []
        for index, site_class in enumerate(equivalent_classes):
            if site_class not in seen:
                seen.add(site_class)
                indices.append(index)
        centers = centers[indices]
    distances = np.linalg.norm(centers[:, None, :] - native[None, :, :], axis=2).min(axis=1) if len(centers) else np.array([])
    result: dict[str, Any] = {
        "distinct_returned_sites": len(centers),
        "center_to_native_min_distances_A": distances.tolist(),
        "top1_distance_A": float(distances[0]) if len(distances) else None,
        "top5_min_distance_A": float(distances[:5].min()) if len(distances) else None,
    }
    for threshold in thresholds:
        if not np.isfinite(threshold) or threshold <= 0:
            raise ValueError("site distance thresholds must be finite and positive")
        result[f"top1_success_{threshold:g}A"] = bool(len(distances) and distances[0] <= threshold)
        result[f"top5_success_{threshold:g}A"] = bool(len(distances) and distances[:5].min() <= threshold)
    return result
