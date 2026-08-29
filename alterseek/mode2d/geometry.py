"""Two-dimensional lattice and Brillouin-zone geometry.

The submitted calculation cell is the source of truth: this module never
replaces an orthogonal submitted supercell by a smaller hidden primitive cell.
"""

from dataclasses import dataclass
import itertools

import numpy as np

@dataclass(frozen=True)
class Lattice2D:
    """Intrinsic geometry of the submitted two-periodic calculation cell."""

    vacuum_axis: int
    in_plane_axes: tuple
    plane_basis: np.ndarray
    plane_normal: np.ndarray
    direct_2d: np.ndarray
    reciprocal_2d: np.ndarray
    reciprocal_3d: np.ndarray
    lattice_class: str
    canonical_transform: np.ndarray
    canonical_direct_2d: np.ndarray
    centered_branch: str | None


def _relative_close(first, second, tol):
    scale = max(abs(float(first)), abs(float(second)), 1.0)
    return abs(float(first) - float(second)) <= tol * scale


def _candidate_basis_transforms(limit=3):
    """Small unimodular 2D transforms, simplest candidates first."""
    candidates = []
    for values in itertools.product(range(-limit, limit + 1), repeat=4):
        matrix = np.array(values, dtype=int).reshape(2, 2)
        if abs(round(np.linalg.det(matrix))) != 1:
            continue
        score = (
            int(np.sum(np.abs(matrix))),
            int(np.sum(np.abs(matrix - np.eye(2, dtype=int)))),
            tuple(values),
        )
        candidates.append((score, matrix))
    return [matrix for _score, matrix in sorted(candidates, key=lambda item: item[0])]


_BASIS_TRANSFORMS = _candidate_basis_transforms()


def _basis_metric_kind(basis, tol):
    first, second = np.asarray(basis, dtype=float)
    lengths = np.linalg.norm(basis, axis=1)
    if np.any(lengths <= 1e-12):
        return None
    cosine = float(np.dot(first, second) / np.prod(lengths))
    equal = _relative_close(lengths[0], lengths[1], tol)
    orthogonal = abs(cosine) <= tol
    hexagonal = equal and abs(abs(cosine) - 0.5) <= tol
    if orthogonal and equal:
        return "square"
    if hexagonal:
        return "hexagonal"
    if orthogonal:
        return "rectangular"
    if equal:
        return "centered_rectangular"
    return None


def _classify_2d_lattice(direct_2d, tol):
    """Classify a 2D translation lattice without changing the submitted cell."""
    matches = {}
    for transform in _BASIS_TRANSFORMS:
        candidate = transform @ direct_2d
        kind = _basis_metric_kind(candidate, tol)
        if kind is not None and kind not in matches:
            matches[kind] = (transform, candidate)

    # Metric specializations must win over their lower-symmetry descriptions.
    for kind in ("square", "hexagonal", "rectangular", "centered_rectangular"):
        if kind in matches:
            transform, candidate = matches[kind]
            branch = None
            if kind == "centered_rectangular":
                cosine = float(
                    np.dot(candidate[0], candidate[1])
                    / np.prod(np.linalg.norm(candidate, axis=1))
                )
                branch = "acute" if cosine > 0.0 else "obtuse"
            return kind, transform, candidate, branch

    return (
        "oblique",
        np.eye(2, dtype=int),
        np.asarray(direct_2d, dtype=float).copy(),
        None,
    )


def analyze_lattice(
    lattice,
    vacuum_axis,
    *,
    orthogonality_tol=2e-3,
    metric_tol=2e-3,
):
    """Build the intrinsic 2D geometry of the submitted calculation cell.

    ``vacuum_axis`` may be any submitted lattice-vector index, but that vector
    must be perpendicular to both periodic vectors and aligned with a
    Cartesian axis.
    """
    lattice = np.asarray(lattice, dtype=float)
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise ValueError("submitted lattice must be a finite 3x3 matrix")
    if vacuum_axis not in (0, 1, 2):
        raise ValueError("vacuum axis must be 0, 1, or 2")

    lengths = np.linalg.norm(lattice, axis=1)
    if np.any(lengths <= 1e-12):
        raise ValueError("submitted lattice contains a zero-length vector")
    in_plane = tuple(index for index in range(3) if index != vacuum_axis)
    vacuum = lattice[vacuum_axis]
    vacuum_unit = vacuum / lengths[vacuum_axis]
    cosines = [
        abs(float(np.dot(vacuum_unit, lattice[index] / lengths[index])))
        for index in in_plane
    ]
    if any(value > orthogonality_tol for value in cosines):
        axis_name = "abc"[vacuum_axis]
        raise RuntimeError(
            f"2D vacuum axis '{axis_name}' must be perpendicular to both "
            "periodic lattice vectors; "
            f"|cos|={','.join(f'{value:.6g}' for value in cosines)}."
        )

    first = lattice[in_plane[0]]
    second = lattice[in_plane[1]]
    normal = np.cross(first, second)
    normal_norm = np.linalg.norm(normal)
    if normal_norm <= 1e-12:
        raise ValueError("submitted periodic lattice vectors are collinear")
    normal /= normal_norm
    if np.dot(normal, vacuum_unit) < 0.0:
        normal = -normal

    # Screen kx/ky are the submitted structure's own fixed Cartesian axes for
    # the two in-plane directions -- never rotated to align with any lattice
    # vector, which is why the vacuum axis must itself be Cartesian-aligned.
    if abs(abs(vacuum_unit[0]) - 1.0) < 1e-6:
        fixed_axes = (1, 2)
    elif abs(abs(vacuum_unit[1]) - 1.0) < 1e-6:
        fixed_axes = (0, 2)
    elif abs(abs(vacuum_unit[2]) - 1.0) < 1e-6:
        fixed_axes = (0, 1)
    else:
        raise RuntimeError(
            "2D vacuum axis must be aligned with a Cartesian axis (x, y, or "
            "z) to define a fixed screen kx/ky frame; got vacuum direction "
            f"{vacuum_unit.tolist()}."
        )
    e1 = np.zeros(3)
    e1[fixed_axes[0]] = 1.0
    e2 = np.zeros(3)
    e2[fixed_axes[1]] = 1.0
    if np.dot(normal, np.cross(e1, e2)) < 0.0:
        e2 = -e2
    plane_basis = np.column_stack((e1, e2))
    direct_2d = np.array(
        [[np.dot(lattice[index], e1), np.dot(lattice[index], e2)]
         for index in in_plane],
        dtype=float,
    )
    reciprocal_2d = 2.0 * np.pi * np.linalg.inv(direct_2d).T
    reciprocal_3d = 2.0 * np.pi * np.linalg.inv(lattice).T
    lattice_class, transform, canonical, branch = _classify_2d_lattice(
        direct_2d, metric_tol
    )
    return Lattice2D(
        vacuum_axis=vacuum_axis,
        in_plane_axes=in_plane,
        plane_basis=plane_basis,
        plane_normal=normal,
        direct_2d=direct_2d,
        reciprocal_2d=reciprocal_2d,
        reciprocal_3d=reciprocal_3d,
        lattice_class=lattice_class,
        canonical_transform=transform,
        canonical_direct_2d=canonical,
        centered_branch=branch,
    )


def _clip_polygon(poly, normal, bound, tol=1e-11):
    """Clip a convex polygon to ``normal . x <= bound``."""
    poly = np.asarray(poly, dtype=float)
    if len(poly) == 0:
        return poly
    result = []
    for index, current in enumerate(poly):
        previous = poly[index - 1]
        current_value = float(np.dot(normal, current) - bound)
        previous_value = float(np.dot(normal, previous) - bound)
        current_inside = current_value <= tol
        previous_inside = previous_value <= tol
        if current_inside != previous_inside:
            denominator = previous_value - current_value
            if abs(denominator) > 1e-15:
                fraction = previous_value / denominator
                result.append(previous + fraction * (current - previous))
        if current_inside:
            result.append(current)
    return _clean_polygon(result)


def _clean_polygon(points, tol=1e-10):
    """Remove consecutive and closing duplicates introduced by clipping."""
    cleaned = []
    for point in np.asarray(points, dtype=float):
        if not cleaned or not np.allclose(point, cleaned[-1], atol=tol, rtol=0.0):
            cleaned.append(point)
    if len(cleaned) > 1 and np.allclose(
        cleaned[0], cleaned[-1], atol=tol, rtol=0.0
    ):
        cleaned.pop()
    return np.array(cleaned, dtype=float)


def polygon_area_centroid_2d(poly):
    """Return the positive area and Cartesian centroid of a 2D polygon."""
    poly = np.asarray(poly, dtype=float)
    if len(poly) < 3:
        raise ValueError("2D polygon must contain at least three vertices")
    x = poly[:, 0]
    y = poly[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    twice_area = float(np.sum(cross))
    if abs(twice_area) <= 1e-14:
        raise ValueError("2D polygon has zero area")
    centroid = np.array([
        np.sum((x + np.roll(x, -1)) * cross),
        np.sum((y + np.roll(y, -1)) * cross),
    ]) / (3.0 * twice_area)
    return abs(twice_area) / 2.0, centroid


def build_bz(reciprocal_2d, grid_radius=3):
    """Construct the Wigner-Seitz BZ of a submitted 2D reciprocal lattice."""
    reciprocal_2d = np.asarray(reciprocal_2d, dtype=float)
    vectors = []
    for first in range(-grid_radius, grid_radius + 1):
        for second in range(-grid_radius, grid_radius + 1):
            if first == 0 and second == 0:
                continue
            vector = first * reciprocal_2d[0] + second * reciprocal_2d[1]
            if np.linalg.norm(vector) > 1e-12:
                vectors.append(vector)
    span = 2.0 * max(np.linalg.norm(vector) for vector in vectors)
    polygon = np.array(
        [[-span, -span], [span, -span], [span, span], [-span, span]],
        dtype=float,
    )
    for vector in vectors:
        polygon = _clip_polygon(
            polygon, vector, 0.5 * float(np.dot(vector, vector))
        )
        if len(polygon) < 3:
            raise RuntimeError("2D Wigner-Seitz construction became degenerate")
    area, _centroid = polygon_area_centroid_2d(polygon)
    direct_area = abs(float(np.linalg.det(2.0 * np.pi * np.linalg.inv(reciprocal_2d).T)))
    expected = (2.0 * np.pi) ** 2 / direct_area
    if not np.isclose(area, expected, rtol=2e-7, atol=2e-9):
        raise RuntimeError(
            "2D BZ area does not match the submitted translation lattice "
            f"({area:.12g} != {expected:.12g})"
        )
    return polygon


def to_input_fractional(points, lattice_2d):
    """Convert internal 2D Cartesian points to submitted reciprocal fractions."""
    points = np.atleast_2d(np.asarray(points, dtype=float))
    cartesian = points @ lattice_2d.plane_basis.T
    fractional = cartesian @ np.linalg.inv(lattice_2d.reciprocal_3d)
    fractional[np.abs(fractional) < 1e-11] = 0.0
    return fractional
