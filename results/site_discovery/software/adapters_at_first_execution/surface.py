"""Reference-independent random patches on the receptor's probe-exposed surface.

The Shrake--Rupley construction samples a deterministic Fibonacci sphere on each
heavy atom at ``van der Waals radius + 1.4 A`` and removes points inside another
expanded atom. A surviving point is a solvent-probe *center*, and is used as the
patch center without moving it toward a ligand or an inferred pocket. Radii are
``gemmi.Element.vdw_r`` from the pinned Gemmi dependency, not covalent radii.

Sources: Shrake & Rupley, J. Mol. Biol. 79 (1973), 351--371,
https://doi.org/10.1016/0022-2836(73)90011-9; Gemmi's documented element table,
https://gemmi.readthedocs.io/en/latest/chemistry.html#elements.

This tests local probe exclusion only. Both capsid cavity and exterior surfaces
are retained. It does not establish a solvent path from bulk water, passage of
a ligand through a capsid, or accessibility to a finite ligand.
"""
from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass
import importlib.metadata
from typing import Any, Iterator

import gemmi
import numpy as np
from scipy.spatial import cKDTree


DEFAULT_SEED = 20261030
DEFAULT_POINTS_PER_ATOM = 48
RADIUS_SOURCE = "gemmi.Element.vdw_r (Angstrom); pinned Gemmi element table"
_DISTANCE_TOLERANCE_A2 = 1e-9


def _positive(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be finite and positive")
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _integer(value: Any, name: str, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def fibonacci_sphere(points_per_atom: int = DEFAULT_POINTS_PER_ATOM) -> np.ndarray:
    """Return a fixed, approximately equal-area unit-sphere quadrature.

    The half-step in z avoids placing a point exactly at either pole. The same
    directions and atom order are used for every run; only patch selection uses
    the frozen random seed. This is a finite-resolution surface approximation.
    """
    count = _integer(points_per_atom, "points_per_atom", 12)
    index = np.arange(count, dtype=float)
    z = 1 - 2 * (index + 0.5) / count
    angle = index * (np.pi * (3 - np.sqrt(5)))
    planar_radius = np.sqrt(np.maximum(0, 1 - z * z))
    return np.column_stack((planar_radius * np.cos(angle), planar_radius * np.sin(angle), z))


def element_radii_A(elements: Sequence[str]) -> np.ndarray:
    """Look up explicit heavy-element van der Waals radii in Gemmi.

    Unknown elements and hydrogen/deuterium fail rather than silently receiving
    a carbon radius. Callers must supply curated receptor heavy atoms only.
    """
    lookup: dict[str, float] = {}
    normalized = []
    for value in elements:
        if not isinstance(value, str) or not value.strip():
            raise ValueError("elements must contain explicit heavy-element symbols")
        symbol = value.strip().capitalize()
        if symbol not in lookup:
            element = gemmi.Element(symbol)
            if element.name != symbol or element.atomic_number <= 1 or element.is_hydrogen:
                raise ValueError(f"unknown or non-heavy receptor element: {value!r}")
            radius = float(element.vdw_r)
            if not np.isfinite(radius) or radius <= 0:
                raise ValueError(f"missing van der Waals radius for receptor element {value!r}")
            lookup[symbol] = radius
        normalized.append(symbol)
    return np.array([lookup[symbol] for symbol in normalized], dtype=float)


def _geometry(coordinates: Any, elements: Sequence[str], probe_radius_A: float,
              points_per_atom: int, batch_atoms: int):
    # Own the array so freezing cached geometry never freezes the caller's
    # coordinates, and caller edits cannot alter later null replicates.
    xyz = np.array(coordinates, dtype=float, copy=True)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("receptor coordinates must be a nonempty finite N x 3 array")
    if len(elements) != len(xyz):
        raise ValueError("element and receptor coordinate counts differ")
    probe = _positive(probe_radius_A, "probe_radius_A")
    directions = fibonacci_sphere(points_per_atom)
    batch_atoms = _integer(batch_atoms, "batch_atoms", 1)
    radii = element_radii_A(elements) + probe
    return xyz, radii, cKDTree(xyz), directions, batch_atoms


def _accessible_batches(xyz, radii, tree, directions, batch_atoms):
    """Yield atom slices and exposure masks, bounding temporary point storage."""
    point_count = len(directions)
    max_radius = float(np.max(radii))
    # Coincident identical spheres have one geometric surface, not two areas.
    # Keep the first owner deterministically; distinct nested spheres still
    # undergo the usual occlusion test.
    _, unique_indices = np.unique(np.column_stack((xyz, radii)), axis=0, return_index=True)
    first_sphere = np.zeros(len(xyz), dtype=bool)
    first_sphere[unique_indices] = True
    for start in range(0, len(xyz), batch_atoms):
        stop = min(start + batch_atoms, len(xyz))
        points = (xyz[start:stop, None, :] + radii[start:stop, None, None]
                  * directions[None, :, :]).reshape(-1, 3)
        # Query every possibly occluding sphere. Nearest-atom-only rejection is
        # incorrect when neighbors have different van der Waals radii.
        neighbors = tree.query_ball_point(points, max_radius, return_sorted=False)
        counts = np.fromiter((len(row) for row in neighbors), dtype=np.intp, count=len(points))
        point_indices = np.repeat(np.arange(len(points), dtype=np.intp), counts)
        atom_indices = np.fromiter((j for row in neighbors for j in row), dtype=np.intp,
                                   count=int(counts.sum()))
        delta = points[point_indices] - xyz[atom_indices]
        buried = np.einsum("ij,ij->i", delta, delta) < radii[atom_indices] ** 2 - _DISTANCE_TOLERANCE_A2
        exposed = np.ones(len(points), dtype=bool)
        exposed[point_indices[buried]] = False
        exposed = exposed.reshape(stop - start, point_count)
        exposed[~first_sphere[start:stop]] = False
        yield start, stop, points, exposed


def iter_solvent_accessible_points(
    coordinates: Any, elements: Sequence[str], *, probe_radius_A: float = 1.4,
    points_per_atom: int = DEFAULT_POINTS_PER_ATOM, batch_atoms: int = 256,
) -> Iterator[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Yield ``(probe_centers_A, source_atom_indices, area_weights_A2)`` batches.

    Each surviving point represents ``4*pi*(vdW+probe)**2/points_per_atom``
    square Angstroms. No full coordinate cloud or pairwise distance matrix is
    retained. Tangent spheres do not bury a point; overlaps do.
    """
    xyz, radii, tree, directions, batch_atoms = _geometry(
        coordinates, elements, probe_radius_A, points_per_atom, batch_atoms)
    atom_areas = 4 * np.pi * radii ** 2 / len(directions)
    for start, stop, points, exposed in _accessible_batches(xyz, radii, tree, directions, batch_atoms):
        retained = exposed.ravel()
        atom_indices = np.repeat(np.arange(start, stop), len(directions))[retained]
        yield points[retained], atom_indices, atom_areas[atom_indices]


def _transforms(values: Iterable[tuple[Any, Any]] | None):
    result = [(np.eye(3), np.zeros(3))]
    for value in values if values is not None else ():
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            raise ValueError("symmetry transform must contain rotation and translation")
        rotation, translation = (np.asarray(v, dtype=float) for v in value)
        if (rotation.shape != (3, 3) or translation.shape != (3,)
                or not np.isfinite(rotation).all() or not np.isfinite(translation).all()):
            raise ValueError("symmetry transform must contain finite 3 x 3 rotation and 3-vector translation")
        if (not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6, rtol=0)
                or not np.isclose(np.linalg.det(rotation), 1, atol=1e-6, rtol=0)):
            raise ValueError("symmetry rotation must be proper and orthonormal")
        if not any(np.allclose(rotation, r, atol=1e-12, rtol=0)
                   and np.allclose(translation, t, atol=1e-12, rtol=0) for r, t in result):
            result.append((rotation, translation))
    return result


def _exclude_patch(mask, eligible_counts, xyz, radii, tree, directions, centers,
                   separation, batch_atoms):
    """Remove only local points; never allocate all atom-sphere coordinates."""
    for center in centers:
        neighbors = np.asarray(tree.query_ball_point(center, separation + float(np.max(radii))), dtype=np.intp)
        for start in range(0, len(neighbors), batch_atoms):
            indices = neighbors[start:start + batch_atoms]
            indices = indices[eligible_counts[indices] > 0]
            if not len(indices):
                continue
            points = xyz[indices, None, :] + radii[indices, None, None] * directions[None, :, :]
            distance_squared = np.sum((points - center) ** 2, axis=2)
            mask[indices] &= distance_squared >= separation ** 2
            eligible_counts[indices] = mask[indices].sum(axis=1)


@dataclass(frozen=True)
class SurfaceGeometry:
    """Cached receptor-only surface quadrature, reusable across random seeds.

    Arrays are read-only. Sampling copies the boolean mask and counts; it never
    modifies this geometry, so all null replicates use the identical surface.
    """
    coordinates_A: np.ndarray
    expanded_radii_A: np.ndarray
    residue_ids: tuple[Hashable, ...]
    tree: cKDTree
    directions: np.ndarray
    exposed_mask: np.ndarray
    exposed_counts: np.ndarray
    atom_point_areas_A2: np.ndarray
    probe_radius_A: float
    batch_atoms: int
    exposed_point_count: int
    estimated_surface_area_A2: float


def build_surface_geometry(
    coordinates: Any, elements: Sequence[str], residue_ids: Sequence[Hashable], *,
    points_per_atom: int = DEFAULT_POINTS_PER_ATOM, probe_radius_A: float = 1.4,
    batch_atoms: int = 256,
) -> SurfaceGeometry:
    """Calculate receptor-only probe exposure once for all baseline replicates."""
    xyz, radii, tree, directions, batch_atoms = _geometry(
        coordinates, elements, probe_radius_A, points_per_atom, batch_atoms)
    if len(residue_ids) != len(xyz):
        raise ValueError("residue identity and receptor coordinate counts differ")
    if any(not isinstance(value, Hashable) or value is None for value in residue_ids):
        raise ValueError("residue identities must be non-null hashable values")
    mask = np.empty((len(xyz), len(directions)), dtype=bool)
    for start, stop, _, exposed in _accessible_batches(xyz, radii, tree, directions, batch_atoms):
        mask[start:stop] = exposed
    counts = mask.sum(axis=1)
    atom_point_areas = 4 * np.pi * radii ** 2 / len(directions)
    for array in (xyz, radii, directions, mask, counts, atom_point_areas):
        array.setflags(write=False)
    return SurfaceGeometry(xyz, radii, tuple(residue_ids), tree, directions, mask,
                           counts, atom_point_areas, float(probe_radius_A), batch_atoms,
                           int(counts.sum()), float(np.dot(counts, atom_point_areas)))


def sample_surface_patches(
    geometry: SurfaceGeometry, *, seed: int = DEFAULT_SEED, max_sites: int = 5,
    patch_radius_A: float = 8.0, min_center_separation_A: float | None = None,
    symmetry_transforms: Iterable[tuple[Any, Any]] | None = None,
) -> dict[str, Any]:
    """Draw distinct, area-weighted random patches from cached receptor geometry.

    See :func:`random_surface_patches` for the center, footprint, random ranking,
    and symmetry definitions. Repeated calls reuse the fixed exposure calculation
    and copy only the boolean eligibility mask and per-atom counts.
    """
    if not isinstance(geometry, SurfaceGeometry):
        raise ValueError("geometry must be built with build_surface_geometry")
    max_sites = _integer(max_sites, "max_sites", 1)
    seed = _integer(seed, "seed", 0)
    patch_radius = _positive(patch_radius_A, "patch_radius_A")
    separation = 2 * patch_radius if min_center_separation_A is None else _positive(
        min_center_separation_A, "min_center_separation_A")
    transforms = _transforms(symmetry_transforms)
    xyz, radii, tree = geometry.coordinates_A, geometry.expanded_radii_A, geometry.tree
    directions, batch_atoms = geometry.directions, geometry.batch_atoms
    residue_ids = geometry.residue_ids
    mask, eligible_counts = geometry.exposed_mask.copy(), geometry.exposed_counts.copy()
    atom_point_areas = geometry.atom_point_areas_A2
    rng = np.random.Generator(np.random.PCG64(seed))
    sites = []
    for rank in range(1, max_sites + 1):
        cumulative_area = np.cumsum(eligible_counts * atom_point_areas)
        if cumulative_area[-1] <= 0:
            break
        draw = min(rng.random() * cumulative_area[-1], np.nextafter(cumulative_area[-1], 0))
        atom_index = int(np.searchsorted(cumulative_area, draw, side="right"))
        point_options = np.flatnonzero(mask[atom_index])
        point_index = int(point_options[rng.integers(len(point_options))])
        center = xyz[atom_index] + radii[atom_index] * directions[point_index]
        footprint_indices = sorted(tree.query_ball_point(center, patch_radius))
        residues = list(dict.fromkeys(residue_ids[i] for i in footprint_indices))
        actual_radius = max((float(np.linalg.norm(xyz[i] - center)) for i in footprint_indices), default=0.0)
        sites.append({"rank": rank, "center_A": center.tolist(),
                      "footprint_residues": residues, "footprint_atom_count": len(footprint_indices),
                      "footprint_radius_A": patch_radius, "footprint_max_atom_distance_A": actual_radius,
                      "source_atom_index": atom_index, "source_sphere_point_index": point_index,
                      "ranking": "seeded area-weighted random draw order"})
        centers = [center @ rotation + translation for rotation, translation in transforms]
        _exclude_patch(mask, eligible_counts, xyz, radii, tree, directions, centers, separation, batch_atoms)
    return {
        "method": "random_solvent_accessible_surface_patches",
        "reference_independent": True,
        "center_definition": "surviving solvent-probe center on a receptor atom vdW+probe sphere",
        "footprint_definition": "residues with any receptor heavy atom within footprint_radius_A of center",
        "accessibility_definition": "local probe exclusion; bulk-solvent or ligand-path connectivity not tested",
        "surface_scope": "all probe-exposed receptor surfaces, including capsid cavity and exterior",
        "sampling_distribution": "sequential area-weighted draws from eligible surface quadrature points",
        "parameters": {"seed": seed, "random_generator": "numpy.random.PCG64",
                       "max_sites": max_sites, "patch_radius_A": patch_radius,
                       "min_center_separation_A": separation, "probe_radius_A": geometry.probe_radius_A,
                       "points_per_atom": len(directions), "sphere_quadrature": "deterministic Fibonacci",
                       "radius_source": RADIUS_SOURCE, "symmetry_transform_count": len(transforms),
                       "symmetry_certification": "caller-supplied receptor/source operators; rigid form checked only"},
        "sampling": {"receptor_heavy_atom_count": len(xyz),
                     "tested_surface_point_count": len(xyz) * len(directions),
                     "probe_exposed_point_count": geometry.exposed_point_count,
                     "estimated_probe_exposed_surface_area_A2": geometry.estimated_surface_area_A2,
                     "surface_mask_bytes": int(mask.nbytes), "batch_atoms": batch_atoms,
                     "remaining_eligible_point_count": int(eligible_counts.sum())},
        "versions": {name: importlib.metadata.version(name) for name in ("gemmi", "numpy", "scipy")},
        "requested_site_count": max_sites, "returned_site_count": len(sites),
        "budget_complete": len(sites) == max_sites, "sites": sites,
    }


def random_surface_patches(
    coordinates: Any, elements: Sequence[str], residue_ids: Sequence[Hashable], *,
    seed: int = DEFAULT_SEED, max_sites: int = 5, patch_radius_A: float = 8.0,
    points_per_atom: int = DEFAULT_POINTS_PER_ATOM, probe_radius_A: float = 1.4,
    min_center_separation_A: float | None = None, batch_atoms: int = 256,
    symmetry_transforms: Iterable[tuple[Any, Any]] | None = None,
) -> dict[str, Any]:
    """Draw ranked, distinct random surface patches without native-site input.

    Rank is the seeded random draw order, without a pocket or docking score.
    At each draw, every eligible surface quadrature point has probability
    proportional to its represented surface area. Points closer than
    ``min_center_separation_A`` to a chosen center are removed before the next
    draw. The default is twice ``patch_radius_A``, so spherical patch interiors
    cannot overlap. The footprint contains every receptor residue with at least
    one heavy atom at or within the supplied radius; this same definition can
    be used for another method's centers.

    Optional symmetry operators use row vectors, ``x @ rotation + translation``.
    They must be independently certified from receptor/source data by the caller;
    this function checks rigid-body form and never infers symmetry or a site.
    Exclusion applies around every transformed selected center.

    Memory is O(N_atoms * points_per_atom) *booleans*, plus the receptor tree
    and bounded point batches (about 12 MB for 248k atoms at 48 points). No
    ligand, reference coordinates, native residues, or native-derived box enter
    this API. A small receptor may exhaust eligible patches before the requested
    budget; ``budget_complete`` makes this explicit rather than duplicating sites.
    """
    geometry = build_surface_geometry(coordinates, elements, residue_ids,
                                      points_per_atom=points_per_atom,
                                      probe_radius_A=probe_radius_A, batch_atoms=batch_atoms)
    return sample_surface_patches(geometry, seed=seed, max_sites=max_sites,
                                  patch_radius_A=patch_radius_A,
                                  min_center_separation_A=min_center_separation_A,
                                  symmetry_transforms=symmetry_transforms)
