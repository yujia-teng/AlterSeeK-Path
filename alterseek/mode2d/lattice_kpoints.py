"""Special k-points, paths, and conventional IBZs for 2D mode."""

import numpy as np

from .geometry import build_bz, polygon_area_centroid_2d


def _to_submitted_fractional(lattice_2d, first, second):
    canonical = np.array([first, second], dtype=float)
    submitted = canonical @ np.linalg.inv(lattice_2d.canonical_transform).T
    point = np.zeros(3)
    point[lattice_2d.in_plane_axes[0]] = submitted[0]
    point[lattice_2d.in_plane_axes[1]] = submitted[1]
    return point.tolist()


def build_path(lattice_2d, point_group):
    """Return special points and the ordinary path in the submitted basis."""
    gamma = "Γ"

    def point(first, second):
        return _to_submitted_fractional(lattice_2d, first, second)

    lattice_class = lattice_2d.lattice_class
    butterfly_path = None
    butterfly_extra = None
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
        if str(point_group).replace(" ", "") == "4":
            points["X_A"] = point(0.5, 0.0)
            butterfly_path = [(gamma, "X"), ("X", "M")]
            butterfly_extra = [gamma, "X_A"]
    elif lattice_class == "hexagonal":
        points = {
            gamma: point(0.0, 0.0),
            "M": point(0.5, 0.0),
            "K": point(1.0 / 3.0, 1.0 / 3.0),
        }
        path = [(gamma, "M"), ("M", "K"), ("K", gamma)]
    elif lattice_class == "centered_rectangular":
        first, second = lattice_2d.canonical_direct_2d
        conventional_lengths = sorted((
            np.linalg.norm(first + second),
            np.linalg.norm(first - second),
        ))
        x_value = 0.25 * (
            1.0 + (conventional_lengths[0] / conventional_lengths[1]) ** 2
        )
        if lattice_2d.centered_branch == "obtuse":
            # oC1 (a < b) restricted to k_z = 0.
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
            # oC2 (a > b) restricted to k_z = 0.
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
    else:
        raise RuntimeError(
            "The submitted calculation cell has an oblique 2D lattice. "
            "No altermagnetic 2D path exists for this class."
        )

    return {
        "points": points,
        "path": path,
        "extra_vertices": [],
        "butterfly_path": butterfly_path,
        "butterfly_extra_vertices": butterfly_extra,
    }


def build_ibz(lattice_2d, path_data, point_group, projected_group_order):
    """Build the conventional path-compatible IBZ polygon."""
    gamma = "Γ"
    lattice_class = lattice_2d.lattice_class
    if lattice_class == "hexagonal":
        labels = [gamma, "M", "K"]
    elif lattice_class == "square":
        labels = (
            [gamma, "X_A", "M", "X"]
            if str(point_group).replace(" ", "") == "4"
            else [gamma, "X", "M"]
        )
    elif lattice_class == "rectangular":
        labels = [gamma, "X", "S", "Y"]
    elif lattice_class == "centered_rectangular":
        labels = (
            ["SIGMA_0", gamma, "Y", "C_0"]
            if lattice_2d.centered_branch == "obtuse"
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
