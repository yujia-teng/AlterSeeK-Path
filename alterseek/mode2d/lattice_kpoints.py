"""Special k-points, paths, and conventional IBZs for 2D mode."""

from dataclasses import replace
import itertools

import numpy as np

from .geometry import (
    _BASIS_TRANSFORMS,
    _clip_polygon,
    build_bz,
    polygon_area_centroid_2d,
    reduced_basis_transform_2d,
)


def _to_submitted_fractional(lattice_2d, first, second):
    canonical = np.array([first, second], dtype=float)
    submitted = canonical @ np.linalg.inv(lattice_2d.canonical_transform).T
    point = np.zeros(3)
    point[lattice_2d.in_plane_axes[0]] = submitted[0]
    point[lattice_2d.in_plane_axes[1]] = submitted[1]
    return point.tolist()


# The default relative metric tolerance of ``analyze_lattice``.
_EQUAL_DISTANCE_TOL = 2e-3


def _fold_oblique_special_point(lattice_2d, first, second, search_limit=3):
    """Return the copy of a Bilbao oblique point that lies nearest Gamma.

    The labels hold only modulo reciprocal lattice vectors, so the literal
    representative need not sit inside the submitted BZ.
    """
    canonical = np.array([first, second], dtype=float)
    base = canonical @ np.linalg.inv(lattice_2d.canonical_transform).T
    reciprocal = lattice_2d.reciprocal_2d
    # The same shifts along the shortest basis reach the BZ copy when the
    # submitted basis is too skewed for the plain range.
    reduction = reduced_basis_transform_2d(reciprocal)
    steps = range(-search_limit, search_limit + 1)
    shifts = list(itertools.product(steps, steps))
    shifts += [
        tuple(int(x) for x in np.array(pair) @ reduction)
        for pair in itertools.product(steps, steps)
    ]
    candidates = []
    for shift_first, shift_second in dict.fromkeys(shifts):
        shift = np.array([shift_first, shift_second], dtype=float)
        fractional = base + shift
        cartesian = fractional @ reciprocal
        # Nearest Gamma, then least shifted.  The remaining parts only
        # fix which of two equally shifted copies wins; shortening them
        # changes the answer on 13 of the 29 2D cases, so they stay.
        score = (
            abs(shift_first) + abs(shift_second),
            abs(shift_first),
            abs(shift_second),
            shift_first,
            shift_second,
        )
        candidates.append((float(np.dot(cartesian, cartesian)), score, fractional))
    # Copies as near Gamma as the nearest within the 2D metric tolerance count
    # as equally near: on a lattice rectangular within that tolerance the four
    # BZ corners are all copies of one point, and rounding must not pick one.
    nearest = min(distance for distance, _score, _fractional in candidates)
    fractional = min(
        (item for item in candidates
         if item[0] <= nearest * (1.0 + _EQUAL_DISTANCE_TOL)),
        key=lambda item: item[1],
    )[2]
    point = np.zeros(3)
    point[lattice_2d.in_plane_axes[0]] = fractional[0]
    point[lattice_2d.in_plane_axes[1]] = fractional[1]
    return point.tolist()


def _cell_rotations(lattice_2d, operations, tol=2e-6):
    """Integer matrix each k-space operation applies to cell coordinates."""
    direct = np.asarray(lattice_2d.direct_2d, dtype=float)
    inverse = np.linalg.inv(direct)
    rotations = []
    for operation in operations:
        # A k-space action C = D^-1 R^-T D moves cell coordinates by R.
        exact = (direct @ np.linalg.inv(operation) @ inverse).T
        rotation = np.rint(exact).astype(int)
        if not np.allclose(exact, rotation, atol=tol, rtol=0.0):
            raise RuntimeError(
                "A projected 2D operation does not map the submitted lattice "
                "onto itself."
            )
        rotations.append(rotation)
    return rotations


def _is_common_axis(vector, rotations):
    """Whether every rotation fixes or reverses an integer cell vector."""
    return all(
        np.array_equal(rotation @ vector, vector)
        or np.array_equal(rotation @ vector, -vector)
        for rotation in rotations
    )


def _positive_screen_direction(vector, tol=1e-10):
    """Choose a deterministic sign for an unoriented in-plane axis."""
    result = np.asarray(vector, dtype=float).copy()
    for component in result:
        if abs(float(component)) <= tol:
            continue
        return result if component > 0.0 else -result
    raise RuntimeError("Cannot orient a zero-length conventional axis.")


def _canonical_hexagonal_centered_basis(lattice_2d, primitive_basis):
    """Orient a hexagonal cell's centred basis by its conventional axes.

    Keying on the axes rather than the submitted primitive pair stops a
    vector swap in the input from rotating the labelled IBZ.
    """
    first, second = np.asarray(primitive_basis, dtype=float)
    axes = [first + second, first - second]
    axes.sort(key=lambda axis: float(np.linalg.norm(axis)))
    conventional_a = _positive_screen_direction(axes[0])
    conventional_b = _positive_screen_direction(axes[1])
    if np.linalg.det(np.vstack((conventional_a, conventional_b))) < 0.0:
        conventional_b = -conventional_b
    canonical = np.array([
        0.5 * (conventional_a - conventional_b),
        0.5 * (conventional_a + conventional_b),
    ])
    transform_float = canonical @ np.linalg.inv(lattice_2d.direct_2d)
    transform = np.rint(transform_float).astype(int)
    if (
        not np.allclose(transform_float, transform, atol=2e-6, rtol=0.0)
        or abs(round(np.linalg.det(transform))) != 1
    ):
        raise RuntimeError(
            "Could not express the conventional-axis-centered basis as a "
            "unimodular change of the submitted 2D lattice."
        )
    return transform, canonical


def _rectangular_basis_candidates(lattice_2d, rotations):
    """The small transforms, or, when none puts its vectors on mirror lines, the same ones applied to the shortest basis."""
    def on_mirrors(transform):
        first, second = transform
        return (
            (_is_common_axis(first, rotations) and _is_common_axis(second, rotations))
            or (_is_common_axis(first + second, rotations)
                and _is_common_axis(first - second, rotations))
        )

    if any(on_mirrors(transform) for transform in _BASIS_TRANSFORMS):
        return _BASIS_TRANSFORMS
    reduction = reduced_basis_transform_2d(lattice_2d.direct_2d)
    return [transform @ reduction for transform in _BASIS_TRANSFORMS]


def _oriented_rectangular_path_lattice(lattice_2d, operations):
    """Return the rectangular or centred basis the two mirrors fix.

    Primitive when both basis vectors lie on mirror lines, otherwise the pair
    whose sum and difference do.
    """
    rotations = _cell_rotations(lattice_2d, operations)
    candidates = []
    for transform in _rectangular_basis_candidates(lattice_2d, rotations):
        first, second = transform
        forms = []
        if (
            _is_common_axis(first, rotations)
            and _is_common_axis(second, rotations)
        ):
            forms.append(("rectangular", 0))
        if (
            _is_common_axis(first + second, rotations)
            and _is_common_axis(first - second, rotations)
        ):
            forms.append(("centered_rectangular", 1))
        for path_class, family_rank in forms:
            # Primitive before centred, then the rewrite closest to the
            # submitted cell.  The distance-from-identity part is load-bearing:
            # dropping it changes the chosen basis on the 2mm cases.
            score = (
                family_rank,
                int(np.sum(np.abs(transform))),
                int(np.sum(np.abs(transform - np.eye(2, dtype=int)))),
                tuple(int(value) for value in transform.ravel()),
            )
            candidates.append((score, path_class, transform))

    if not candidates:
        raise RuntimeError(
            "Could not orient a 2D rectangular path basis from the detected "
            "four-operation in-plane point group."
        )
    _score, path_class, transform = min(candidates, key=lambda item: item[0])
    basis = transform @ np.asarray(lattice_2d.direct_2d, dtype=float)
    if (
        lattice_2d.lattice_class == "hexagonal"
        and path_class == "centered_rectangular"
    ):
        transform, basis = _canonical_hexagonal_centered_basis(
            lattice_2d, basis
        )
    branch = None
    if path_class == "centered_rectangular":
        cosine = float(np.dot(basis[0], basis[1]) / np.prod(
            np.linalg.norm(basis, axis=1)
        ))
        branch = "acute" if cosine > 0.0 else "obtuse"
    return replace(
        lattice_2d,
        lattice_class=path_class,
        canonical_transform=np.asarray(transform, dtype=int),
        canonical_direct_2d=np.asarray(basis, dtype=float),
        centered_branch=branch,
    )


def _select_path_lattice(lattice_2d, projected_operations):
    """Return the temporary path setting; the physical BZ is untouched."""
    if projected_operations is None:
        return lattice_2d
    operations = [
        np.asarray(operation, dtype=float)
        for operation in projected_operations
    ]
    order = len(operations)
    has_reflection = any(
        np.linalg.det(operation) < -0.5 for operation in operations
    )
    if order == 2:
        # The BZ stays; only the labels drop to the p2 convention.
        return replace(
            lattice_2d,
            lattice_class="oblique",
            canonical_transform=np.eye(2, dtype=int),
            canonical_direct_2d=np.asarray(lattice_2d.direct_2d, dtype=float),
            centered_branch=None,
        )
    if order == 4 and has_reflection:
        return _oriented_rectangular_path_lattice(lattice_2d, operations)
    return lattice_2d


def _path_metric_branch(path_lattice):
    if path_lattice.lattice_class != "centered_rectangular":
        return None
    cosine = float(
        np.dot(*path_lattice.canonical_direct_2d)
        / np.prod(np.linalg.norm(path_lattice.canonical_direct_2d, axis=1))
    )
    if abs(cosine) <= 2e-5:
        return "a_eq_b"
    return "a_lt_b" if cosine < 0.0 else "a_gt_b"


def build_path(
    lattice_2d,
    point_group,
    projected_group_order=None,
    projected_operations=None,
):
    """Return special points and the ordinary path in the submitted basis."""
    gamma = "Γ"

    path_lattice = _select_path_lattice(lattice_2d, projected_operations)

    def point(first, second):
        return _to_submitted_fractional(path_lattice, first, second)

    lattice_class = path_lattice.lattice_class
    butterfly_path = None
    butterfly_extra = None
    extra_vertices = []
    square_fourfold = False
    if lattice_class == "rectangular":
        points = {
            gamma: point(0.0, 0.0),
            "X": point(0.5, 0.0),
            "S": point(0.5, 0.5),
            "Y": point(0.0, 0.5),
        }
        path = [(gamma, "X"), ("X", "S"), ("S", "Y"), ("Y", gamma)]
    elif lattice_class == "square":
        points = {
            gamma: point(0.0, 0.0),
            "X": point(0.0, 0.5),
            "M": point(0.5, 0.5),
        }
        path = [(gamma, "X"), ("X", "M"), ("M", gamma)]
        square_fourfold = (
            projected_operations is not None
            and projected_group_order is not None
            and int(projected_group_order) == 4
        ) or str(point_group).replace(" ", "") == "4"
        if square_fourfold:
            points["X_A"] = point(0.5, 0.0)
            butterfly_path = list(path)
            butterfly_extra = [gamma, "X_A"]
    elif lattice_class == "hexagonal":
        points = {
            gamma: point(0.0, 0.0),
            "M": point(0.5, 0.0),
            "K": point(1.0 / 3.0, 1.0 / 3.0),
        }
        path = [(gamma, "M"), ("M", "K"), ("K", gamma)]
        if (projected_group_order is not None
                and int(projected_group_order) == 6):
            points["M_A"] = _hexagonal_ma_point(
                path_lattice, points["M"], points["K"]
            )
            butterfly_path = list(path)
            butterfly_extra = [gamma, "M_A"]
            # The copied vertex joins no ordinary segment, so it reaches the
            # written path as a general-k excursion, as the 3D -3 corners do.
            extra_vertices = ["M_A"]
    elif lattice_class == "centered_rectangular":
        first, second = path_lattice.canonical_direct_2d
        conventional_lengths = sorted((
            np.linalg.norm(first + second),
            np.linalg.norm(first - second),
        ))
        x_value = 0.25 * (
            1.0 + (conventional_lengths[0] / conventional_lengths[1]) ** 2
        )
        if path_lattice.centered_branch == "obtuse":
            # Centered rectangular, a < b.
            points = {
                gamma: point(0.0, 0.0),
                "Y": point(-0.5, 0.5),
                "C_0": point(-x_value, 1.0 - x_value),
                "SIGMA_0": point(x_value, x_value),
                "S": point(0.0, 0.5),
            }
            path = [
                (gamma, "Y"), ("Y", "C_0"),
                ("SIGMA_0", gamma), (gamma, "S"),
            ]
        else:
            # Centered rectangular, a > b.
            points = {
                gamma: point(0.0, 0.0),
                "Y": point(0.5, 0.5),
                "F_0": point(x_value, 1.0 - x_value),
                "DELTA_0": point(-x_value, x_value),
                "S": point(0.0, 0.5),
            }
            path = [
                (gamma, "Y"), ("Y", "F_0"),
                ("DELTA_0", gamma), (gamma, "S"),
            ]
    elif lattice_class == "oblique":
        # Bilbao p2 labels, shared by point groups 1 and 2.
        points = {
            gamma: point(0.0, 0.0),
            "Y": _fold_oblique_special_point(path_lattice, 0.0, 0.5),
            "B": _fold_oblique_special_point(path_lattice, 0.5, 0.0),
            "A": _fold_oblique_special_point(path_lattice, 0.5, -0.5),
        }
        path = [(gamma, "B"), ("Y", gamma), (gamma, "A")]
    else:
        raise RuntimeError(f"Unsupported 2D lattice class: {lattice_class}.")

    return {
        "points": points,
        "path": path,
        "extra_vertices": extra_vertices,
        "butterfly_path": butterfly_path,
        "butterfly_extra_vertices": butterfly_extra,
        "path_lattice": path_lattice,
        "path_lattice_class": lattice_class,
        "path_metric_branch": _path_metric_branch(path_lattice),
        "square_fourfold": lattice_class == "square" and square_fourfold,
    }


# Number of in-plane operations each curated IBZ polygon was drawn for.
_CURATED_IBZ_ORDERS = {
    "hexagonal": 12,
    "square": 8,
    "rectangular": 4,
    "centered_rectangular": 4,
}


def _in_plane_fractional(lattice_2d, point):
    """Submitted in-plane fractional pair of a path point."""
    coords = np.asarray(point, dtype=float)
    return np.array([coords[axis] for axis in lattice_2d.in_plane_axes])


def _half_bz_normal(lattice_2d, path_data):
    """Fractional normal of the line whose positive side is the group-2 IBZ.

    Every line through Gamma cuts an admissible half, so the choice is a
    convention: the 2D form of the rule 3D uses for Laue group -1.
    """
    path_lattice = path_data.get("path_lattice", lattice_2d)
    if path_lattice.lattice_class == "oblique":
        # Bilbao's p2 representation domain is cut through Gamma along the
        # Y direction and keeps the side containing B.  Y is taken where it
        # sits in the BZ, which is b2/2 only when b2 is a shortest vector.
        y = _in_plane_fractional(lattice_2d, path_data["points"]["Y"])
        b = _in_plane_fractional(lattice_2d, path_data["points"]["B"])
        normal = np.array([y[1], -y[0]])
        return normal if float(normal @ b) > 0.0 else -normal

    points = {}
    for label, coords in path_data["points"].items():
        value = _in_plane_fractional(lattice_2d, coords)
        if np.allclose(value, 0.0):
            continue  # Gamma is the cut's own origin, never a side test
        points[label] = value
    for axis in range(2):
        axial = np.zeros(2)
        axial[axis] = 0.5
        axial_label = next(
            (label for label, value in points.items()
             if np.allclose(value, axial)),
            None,
        )
        if axial_label is None:
            continue
        normal = np.zeros(2)
        normal[1 - axis] = 1.0
        if all(
            float(np.dot(normal, value)) > 1e-12
            for label, value in points.items() if label != axial_label
        ):
            return normal
    return np.array([1.0, 0.0])


def _half_bz_corner_labels(lattice_2d, path_data, polygon):
    """Name every half-BZ corner after the special point it copies.

    A corner sitting on a special point keeps that name; the rest take their
    base name, matched by distance from Gamma, plus the next free suffix.
    """
    reciprocal = lattice_2d.reciprocal_2d
    special = []
    for label, coords in path_data["points"].items():
        plane = _in_plane_fractional(lattice_2d, coords) @ reciprocal
        radius = float(np.linalg.norm(plane))
        if radius < 1e-9:
            continue  # Gamma sits on the cut, never a corner
        special.append((radius, label, plane))

    labels = [""] * len(polygon)
    used = set()
    for index, vertex in enumerate(polygon):
        exact = next(
            (label for _, label, plane in special
             if np.allclose(plane, vertex, atol=1e-6)),
            None,
        )
        if exact is not None:
            labels[index] = exact
            used.add(exact)
    for index, vertex in enumerate(polygon):
        if labels[index]:
            continue
        radius = float(np.linalg.norm(vertex))
        base = next(
            (label for known, label, _ in special
             if abs(known - radius) < 1e-6),
            None,
        )
        if base is None:
            continue
        name, suffix = base, 0
        while name in used:
            suffix += 1
            name = f"{base}_{chr(ord('A') + suffix - 1)}"
        used.add(name)
        labels[index] = name
    return labels


def _oblique_half_bz_labels(lattice_2d, path_data, polygon):
    """Label every oblique half-BZ corner without inventing special points.

    The cut endpoints are Bilbao's Y orbit.  Q and its suffixes are
    AlterSeeK-Path geometric labels, not further Bilbao k-vector types.
    When the lattice is rectangular but only C2 survives, Y can sit on a BZ
    corner; the cut then runs corner to corner, the half BZ is a triangle,
    and its third corner is one more copy of Y, named Y_B.
    """
    submitted_fractional = [
        _in_plane_fractional(
            lattice_2d, _submitted_from_plane(lattice_2d, vertex)
        )
        for vertex in polygon
    ]
    path_lattice = path_data.get("path_lattice", lattice_2d)
    transform_t = path_lattice.canonical_transform.T
    fractional = [point @ transform_t for point in submitted_fractional]
    y_submitted = _in_plane_fractional(
        lattice_2d, path_data["points"]["Y"]
    )
    y_point = y_submitted @ transform_t
    y_indices = [
        index for index, point in enumerate(fractional)
        if np.allclose(
            point - y_point - np.rint(point - y_point),
            0.0,
            atol=1e-7,
            rtol=0.0,
        )
    ]
    if len(polygon) == 3 and len(y_indices) == 3:
        ends = next(
            (first, second)
            for first in range(3) for second in range(first + 1, 3)
            if np.allclose(polygon[first], -polygon[second], atol=1e-9)
        )
        labels = [""] * 3
        y_index = min(
            ends,
            key=lambda index: float(np.linalg.norm(fractional[index] - y_point)),
        )
        labels[y_index] = "Y"
        labels[next(index for index in ends if index != y_index)] = "Y_A"
        labels[labels.index("")] = "Y_B"
        return labels
    if len(y_indices) != 2 or len(polygon) < 4:
        raise RuntimeError(
            "The p2 half-BZ must have two Y-orbit cut endpoints and at "
            "least two generic Wigner-Seitz corners."
        )

    labels = [""] * len(polygon)
    y_index = min(
        y_indices,
        key=lambda index: float(np.linalg.norm(fractional[index] - y_point)),
    )
    labels[y_index] = "Y"
    labels[next(index for index in y_indices if index != y_index)] = "Y_A"

    generic = [index for index, label in enumerate(labels) if not label]
    q_index = max(
        generic,
        key=lambda index: (fractional[index][0], fractional[index][1]),
    )
    labels[q_index] = "Q"
    remaining = sorted(
        (index for index in generic if index != q_index),
        key=lambda index: fractional[index][1],
        reverse=True,
    )
    for suffix, index in enumerate(remaining, start=1):
        labels[index] = f"Q_{chr(ord('A') + suffix - 1)}"
    return labels


def _half_bz_ibz(lattice_2d, path_data):
    """Half the BZ, for an in-plane group of order 2."""
    normal_frac = _half_bz_normal(lattice_2d, path_data)
    # n . f is the same test as (R^-1 n) . p for a plane-frame point p = f @ R.
    reciprocal = lattice_2d.reciprocal_2d
    path_lattice = path_data.get("path_lattice", lattice_2d)
    if path_lattice.lattice_class == "oblique":
        reciprocal = (
            np.linalg.inv(path_lattice.canonical_transform).T @ reciprocal
        )
        # A and -A are the same point; keep the copy on B's side of the cut.
        a_point = _in_plane_fractional(lattice_2d, path_data["points"]["A"])
        if float(normal_frac @ a_point) < -1e-12:
            path_data["points"]["A"] = [-float(x) for x in path_data["points"]["A"]]
    normal_plane = np.linalg.inv(reciprocal) @ normal_frac
    bz = build_bz(lattice_2d.reciprocal_2d)
    polygon = _clip_polygon(bz, -normal_plane, 0.0)
    area, centroid = polygon_area_centroid_2d(polygon)
    labels = (
        _oblique_half_bz_labels(lattice_2d, path_data, polygon)
        if path_lattice.lattice_class == "oblique"
        else _half_bz_corner_labels(lattice_2d, path_data, polygon)
    )
    if path_lattice.lattice_class == "oblique":
        first, second = ("Q", "Q_A") if "Q" in labels else ("Y_B", "Y")
        start = labels.index(first)
        polygon = np.roll(polygon, -start, axis=0)
        labels = labels[start:] + labels[:start]
        if labels[1] != second:
            polygon = np.concatenate([polygon[:1], polygon[:0:-1]])
            labels = labels[:1] + labels[:0:-1]
    # Register the copied corners so they are drawn and labeled like any other
    # special point, and reach the written path as general-k excursions.
    for label, vertex in zip(labels, polygon):
        if not label or label in path_data["points"]:
            continue
        path_data["points"][label] = _submitted_from_plane(lattice_2d, vertex)
        path_data["extra_vertices"].append(label)
    return polygon, centroid, area, labels


def _hexagonal_ma_point(lattice_2d, m_point, k_point):
    """M's partner across Gamma-K, the extra corner of the enlarged domain.

    Without the mirrors of 6mm the edge midpoint on each side of K is no
    longer equivalent, as X_A is on the square mesh at point group 4.
    """
    reciprocal = lattice_2d.reciprocal_2d
    m_plane = _in_plane_fractional(lattice_2d, m_point) @ reciprocal
    k_plane = _in_plane_fractional(lattice_2d, k_point) @ reciprocal
    direction = k_plane / np.linalg.norm(k_plane)
    reflected = 2.0 * float(np.dot(m_plane, direction)) * direction - m_plane
    return _submitted_from_plane(lattice_2d, reflected)


def _submitted_from_plane(lattice_2d, plane_point):
    """Submitted-basis fractional 3-vector of an in-plane Cartesian k point."""
    fractional = plane_point @ np.linalg.inv(lattice_2d.reciprocal_2d)
    point = np.zeros(3)
    point[lattice_2d.in_plane_axes[0]] = fractional[0]
    point[lattice_2d.in_plane_axes[1]] = fractional[1]
    return point.tolist()


def build_ibz(lattice_2d, path_data, point_group, projected_group_order):
    """Build the IBZ polygon for the projected in-plane point group.

    A curated polygon is used when the projected group is the one the lattice
    class was drawn for; otherwise the magnetic order has lowered the symmetry
    below the lattice shape and the domain is derived instead.
    """
    gamma = "Γ"
    path_lattice = path_data.get("path_lattice", lattice_2d)
    lattice_class = path_lattice.lattice_class
    order = int(projected_group_order)
    if lattice_class == "oblique":
        if order != 2:
            raise RuntimeError(
                "An oblique 2D path setting must have an "
                f"inversion-extended point-group order of 2, not {order}."
            )
        return _half_bz_ibz(lattice_2d, path_data)
    curated_order = _CURATED_IBZ_ORDERS.get(lattice_class)
    square_fourfold = (
        lattice_class == "square"
        and path_data.get(
            "square_fourfold", str(point_group).replace(" ", "") == "4"
        )
    )
    if square_fourfold:
        curated_order = 4
    if lattice_class == "hexagonal" and order == 6:
        curated_order = 6
    if curated_order is not None and order != curated_order:
        if order == 2:
            return _half_bz_ibz(lattice_2d, path_data)
        raise RuntimeError(
            f"The magnetic order lowers the in-plane point group on a "
            f"{lattice_class} lattice to order {order}, for which no IBZ is "
            "defined yet."
        )
    if lattice_class == "hexagonal":
        labels = [gamma, "M", "K", "M_A"] if order == 6 else [gamma, "M", "K"]
    elif lattice_class == "square":
        labels = (
            [gamma, "X_A", "M", "X"]
            if square_fourfold
            else [gamma, "X", "M"]
        )
    elif lattice_class == "rectangular":
        labels = [gamma, "X", "S", "Y"]
    elif lattice_class == "centered_rectangular":
        labels = (
            ["SIGMA_0", gamma, "Y", "C_0"]
            if path_lattice.centered_branch == "obtuse"
            else ["DELTA_0", gamma, "Y", "F_0"]
        )
    else:
        raise RuntimeError(
            f"No conventional altermagnetic IBZ exists for {lattice_class}."
        )

    fractional = np.array([path_data["points"][label] for label in labels])
    cartesian = fractional @ lattice_2d.reciprocal_3d
    polygon = cartesian @ lattice_2d.plane_basis
    area, centroid = polygon_area_centroid_2d(polygon)
    bz = build_bz(lattice_2d.reciprocal_2d)
    bz_area, _ = polygon_area_centroid_2d(bz)
    expected_ratio = 1.0 / int(projected_group_order)
    ratio = area / bz_area
    if not np.isclose(ratio, expected_ratio, rtol=1e-5, atol=2e-8):
        raise RuntimeError(
            "Conventional 2D IBZ does not tile the BZ under the projected "
            f"point group ({ratio:.12g} != {expected_ratio:.12g})."
        )
    return polygon, centroid, area, labels
