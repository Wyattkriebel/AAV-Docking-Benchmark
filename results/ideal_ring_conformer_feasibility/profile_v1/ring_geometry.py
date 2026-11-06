"""Pure six-atom Cremer–Pople geometry; no benchmark/reference-coordinate I/O.

Original convention: Cremer and Pople, JACS97(1975)1354–1358,
    DOI 10.1021/ja00839a011, equations 8–14 and 22.
"""
from __future__ import annotations

import numpy as np

RING_ATOM_ORDER = ("C1", "C2", "C3", "C4", "C5", "O5")
PAPER = "https://s3.smu.edu/dedman/catco/publications/pdf/JACS_97_1354_1975.pdf"


def cremer_pople6(coordinates):
    """Original 1975 equations 8–14 and 22, with cyclic j=0,...,5.

    Mean-plane normal is R_sin cross R_cos. The phase sine term carries
    the original minus sign. This CP plane differs from the SVD plane.
    """
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.shape != (6, 3) or not np.isfinite(xyz).all():
        raise ValueError("Six finite, cyclically ordered ring coordinates are required")
    centered = xyz - xyz.mean(axis=0)
    j = np.arange(6)
    phase = 2 * np.pi * j / 6
    normal = np.cross(np.sin(phase) @ centered, np.cos(phase) @ centered)
    norm = float(np.linalg.norm(normal))
    if norm < 1e-12:
        raise ValueError("Degenerate Cremer–Pople mean plane")
    normal /= norm
    z = centered @ normal
    q2cos = float(np.sqrt(1 / 3) * (np.cos(2 * phase) @ z))
    q2sin = float(-np.sqrt(1 / 3) * (np.sin(2 * phase) @ z))
    q2 = float(np.hypot(q2cos, q2sin))
    q3 = float(((-1.0) ** j @ z) / np.sqrt(6))
    amplitude = float(np.hypot(q2, q3))
    theta = float(np.degrees(np.arctan2(q2, q3))) if amplitude > 1e-12 else None
    phi = float(np.degrees(np.arctan2(q2sin, q2cos)) % 360) if q2 > 1e-12 else None
    if phi is not None and abs(phi - 360) < 1e-10:
        phi = 0.0
    svd_normal = np.linalg.svd(centered, full_matrices=False)[2][-1]
    svd_z = centered @ svd_normal
    bonds = np.linalg.norm(xyz - np.roll(xyz, -1, axis=0), axis=1)
    angles = []
    for atom in j:
        a, b = xyz[(atom - 1) % 6] - xyz[atom], xyz[(atom + 1) % 6] - xyz[atom]
        angles.append(float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))))
    reconstructed = np.sqrt(1 / 3) * (q2cos * np.cos(2 * phase) - q2sin * np.sin(2 * phase)) + q3 / np.sqrt(6) * (-1.0) ** j
    return {"Q_A": amplitude, "q2_A": q2, "q3_A": q3, "theta_deg": theta, "phi_deg": phi,
            "q2_cos_phi_A": q2cos, "q2_sin_phi_A": q2sin,
            "chair_pole_distance_deg": min(theta, 180 - theta) if theta is not None else None,
            "cp_plane_rms_A": float(np.sqrt(np.mean(z ** 2))), "cp_plane_max_abs_A": float(np.max(np.abs(z))),
            "least_squares_plane_rms_A": float(np.sqrt(np.mean(svd_z ** 2))),
            "ring_bond_min_A": float(min(bonds)), "ring_bond_max_A": float(max(bonds)),
            "ring_angle_min_deg": min(angles), "ring_angle_max_deg": max(angles),
            "mean_plane_z_A": z.tolist(), "ring_bond_lengths_A": bonds.tolist(), "ring_angles_deg": angles,
            "fourier_reconstruction_max_error_A": float(np.max(np.abs(reconstructed - z))),
            "cp_mean_plane_constraint_max_error_A": float(max(abs(z.sum()), abs(np.cos(phase) @ z), abs(np.sin(phase) @ z)))}


def compare_to_ideal(coordinates, ideal_coordinates):
    """Compare identically named cyclic rings; do not impose validity thresholds.

    The full puckering-vector angle retains the distinction between opposite
    chairs and avoids a phase-only threshold near a chair pole.
    """
    observed, ideal = cremer_pople6(coordinates), cremer_pople6(ideal_coordinates)
    a = np.array([observed[key] for key in ("q2_cos_phi_A", "q2_sin_phi_A", "q3_A")])
    b = np.array([ideal[key] for key in ("q2_cos_phi_A", "q2_sin_phi_A", "q3_A")])
    angle = float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))) if min(np.linalg.norm(a), np.linalg.norm(b)) > 1e-12 else None
    return {"pucker_vector_angle_to_CCD_ideal_deg": angle,
            "theta_difference_to_CCD_ideal_deg": abs(observed["theta_deg"] - ideal["theta_deg"]) if observed["theta_deg"] is not None and ideal["theta_deg"] is not None else None,
            "Q_difference_to_CCD_ideal_A": observed["Q_A"] - ideal["Q_A"],
            "max_abs_ring_bond_difference_to_CCD_ideal_A": float(np.max(np.abs(np.asarray(observed["ring_bond_lengths_A"]) - ideal["ring_bond_lengths_A"]))),
            "max_abs_ring_angle_difference_to_CCD_ideal_deg": float(np.max(np.abs(np.asarray(observed["ring_angles_deg"]) - ideal["ring_angles_deg"])))}


def analytic_ring(q2, q3, phi_deg=0):
    phase = np.arange(6) * np.pi / 3
    z = q2 / np.sqrt(3) * np.cos(2 * phase + np.radians(phi_deg)) + q3 / np.sqrt(6) * (-1.0) ** np.arange(6)
    # Clockwise ordering as seen from positive z gives the original n=+z.
    return np.column_stack([1.3 * np.sin(phase), 1.3 * np.cos(phase), z])


def self_validation():
    controls = []
    for label, q2, q3, phi, theta in (("ideal_chair_north", 0, .6, 0, 0), ("ideal_chair_south", 0, -.6, 0, 180),
                                     ("ideal_boat", .6, 0, 0, 90), ("ideal_twist_boat", .6, 0, 30, 90)):
        xyz = analytic_ring(q2, q3, phi)
        value = cremer_pople6(xyz)
        assert abs(value["Q_A"] - .6) < 1e-12 and abs(value["theta_deg"] - theta) < 1e-10
        if q2:
            assert abs((value["phi_deg"] - phi + 180) % 360 - 180) < 1e-10
        for shift in range(6):
            moved = cremer_pople6(np.roll(xyz, shift, axis=0))
            assert abs(moved["Q_A"] - value["Q_A"]) < 1e-12
            assert abs(moved["chair_pole_distance_deg"] - value["chair_pole_distance_deg"]) < 1e-10
        reversed_ring = cremer_pople6(xyz[::-1])
        assert abs(reversed_ring["chair_pole_distance_deg"] - value["chair_pole_distance_deg"]) < 1e-10
        rotation, _ = np.linalg.qr(np.array([[.3, .4, .8], [.2, -.7, .4], [.6, .2, -.4]]))
        rotation[:, 0] *= np.linalg.det(rotation)
        rigid = cremer_pople6(xyz @ rotation + [21.3, -10.4, 7.2])
        assert abs(rigid["Q_A"] - value["Q_A"]) < 1e-12
        assert abs(rigid["theta_deg"] - value["theta_deg"]) < 1e-10
        controls.append({"control": label, "Q_A": value["Q_A"], "theta_deg": value["theta_deg"], "phi_deg": value["phi_deg"]})
    # Published Table I, p1357. This test intentionally uses its O5,C1,...,C5
    # ordering, not the C1,...,O5 ordering used for the diagnostic table.
    sucrose = [[0, 1.3839, .1976], [1.1997, .7624, -.2106], [1.2356, -.7040, .2393],
               [.0110, -1.4564, -.2550], [-1.2300, -.7208, .2420], [-1.2164, .7350, -.2133]]
    value = cremer_pople6(sucrose)
    assert abs(value["Q_A"] - .556) < .001 and abs(value["theta_deg"] - 5.2) < .1
    assert abs(value["q2_A"] - .050) < .001 and abs(value["q3_A"] - .554) < .001
    assert abs(value["phi_deg"] - 183.7) < 1
    controls.append({"control": "Cremer_Pople_1975_TableI_sucrose", **value,
                     "source": PAPER, "published_Q_A": .556, "published_theta_deg": 5.2,
                     "published_q2_A": .050, "published_q3_A": .554, "published_phi_deg": 183.7,
                     "validation_tolerances": {"amplitudes_A": .001, "theta_deg": .1, "phi_deg": 1},
                     "table_precision_note": "Recomputing from the printed Cartesian coordinates gives theta5.13193 and phi183.10898 rather than the separately printed5.2 and183.7; residuals are retained, not hidden by exact-equality claims."})
    planar = cremer_pople6(analytic_ring(0, 0))
    assert planar["Q_A"] < 1e-12 and planar["theta_deg"] is None
    return {"status": "passed", "controls": controls, "planar_limit_validated": True,
            "cyclic_shifts_preserve_amplitude_and_chair_pole_distance": True,
            "ring_reversal_preserves_chair_pole_distance": True, "proper_rigid_transform_invariance": True}


