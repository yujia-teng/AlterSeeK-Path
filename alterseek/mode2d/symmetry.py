"""Projected symmetry operations for two-dimensional mode."""

import numpy as np

from ..symmetry import reciprocal_operation_cartesian


def project_point_operations(lattice_2d, rotations, *, add_inversion=True):
    """Project submitted-cell point rotations into the physical k plane."""
    operations = []
    for rotation in rotations:
        cartesian = reciprocal_operation_cartesian(
            rotation, lattice_2d.reciprocal_3d
        )
        leakage = lattice_2d.plane_normal @ cartesian @ lattice_2d.plane_basis
        if np.linalg.norm(leakage) > 2e-6:
            continue
        projected = lattice_2d.plane_basis.T @ cartesian @ lattice_2d.plane_basis
        candidates = (projected, -projected) if add_inversion else (projected,)
        for candidate in candidates:
            if not any(
                np.allclose(candidate, old, atol=2e-6) for old in operations
            ):
                operations.append(candidate)
    if not operations:
        raise RuntimeError("No submitted-cell point operation preserves the 2D plane.")
    return operations
