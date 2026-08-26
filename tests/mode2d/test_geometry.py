import numpy as np
import pytest
from scipy.spatial import ConvexHull

from alterseek.compute_centroid_hybrid import run as compute_centroid
from alterseek.mode2d.geometry import (
    analyze_lattice,
    build_bz,
    polygon_area_centroid_2d,
    to_input_fractional,
)
from alterseek.mode2d.lattice_kpoints import build_ibz, build_path
from alterseek.plotting_common import GAMMA_LABEL


def _rotation_z(angle_degrees):
    angle = np.radians(angle_degrees)
    return np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])


@pytest.mark.parametrize(
    ("lattice", "expected"),
    [
        (np.diag([4.0, 6.0, 20.0]), "rectangular"),
        (np.diag([4.0, 4.0, 20.0]), "square"),
        (
            np.array([
                [4.0, 0.0, 0.0],
                [2.0, 2.0 * np.sqrt(3.0), 0.0],
                [0.0, 0.0, 20.0],
            ]),
            "hexagonal",
        ),
        (
            np.array([
                [4.0, 0.0, 0.0],
                [4.0 * np.cos(np.radians(73.0)),
                 4.0 * np.sin(np.radians(73.0)), 0.0],
                [0.0, 0.0, 20.0],
            ]),
            "centered_rectangular",
        ),
        (
            np.array([
                [4.0, 0.0, 0.0],
                [1.3, 5.2, 0.0],
                [0.0, 0.0, 20.0],
            ]),
            "oblique",
        ),
    ],
)
def test_classifies_submitted_2d_translation_lattice(lattice, expected):
    lattice_2d = analyze_lattice(lattice, 2)
    assert lattice_2d.lattice_class == expected


def test_arbitrary_cartesian_orientation_and_vacuum_index_are_internal_only():
    canonical = np.diag([4.0, 6.0, 20.0])
    rotated = canonical @ _rotation_z(37.0).T
    reordered = rotated[[2, 0, 1]]
    lattice_2d = analyze_lattice(reordered, 0)

    assert lattice_2d.lattice_class == "rectangular"
    assert lattice_2d.vacuum_axis == 0
    assert lattice_2d.in_plane_axes == (1, 2)
    assert np.allclose(
        lattice_2d.direct_2d @ lattice_2d.reciprocal_2d.T,
        2.0 * np.pi * np.eye(2),
    )


def test_nonperpendicular_vacuum_is_rejected():
    lattice = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 6.0, 0.0],
        [0.5, 0.0, 20.0],
    ])
    with pytest.raises(RuntimeError, match="must be perpendicular"):
        analyze_lattice(lattice, 2)


def test_submitted_orthogonal_supercell_keeps_rectangular_bz():
    lattice_2d = analyze_lattice(np.diag([4.0, 6.0, 20.0]), 2)
    polygon = build_bz(lattice_2d.reciprocal_2d)
    area, _ = polygon_area_centroid_2d(polygon)

    assert lattice_2d.lattice_class == "rectangular"
    assert len(polygon) == 4
    assert np.isclose(area, (2.0 * np.pi) ** 2 / 24.0)


def test_submitted_centered_rectangular_primitive_has_hexagonal_bz():
    lattice = np.array([
        [5.1103299815503842, 0.0, 0.0],
        [-2.5057448085893013, 4.4538427761384352, 0.0],
        [0.0, 0.0, 22.4751439124],
    ])
    lattice_2d = analyze_lattice(lattice, 2)
    polygon = build_bz(lattice_2d.reciprocal_2d)

    assert lattice_2d.lattice_class == "centered_rectangular"
    assert lattice_2d.centered_branch == "obtuse"
    assert len(polygon) == 6


def test_centered_rectangular_obtuse_path_points_stay_inside_bz():
    lattice = np.array([
        [5.1103299815503842, 0.0, 0.0],
        [-2.5057448085893013, 4.4538427761384352, 0.0],
        [0.0, 0.0, 22.4751439124],
    ])
    lattice_2d = analyze_lattice(lattice, 2)
    bz = build_bz(lattice_2d.reciprocal_2d)
    hull = ConvexHull(bz)
    path_data = build_path(lattice_2d, "m")

    assert path_data["path"] == [
        (GAMMA_LABEL, "Y"),
        ("Y", "C_0"),
        ("SIGMA_0", GAMMA_LABEL),
        (GAMMA_LABEL, "S"),
    ]
    for point in path_data["points"].values():
        cartesian = np.asarray(point) @ lattice_2d.reciprocal_3d
        screen = cartesian @ lattice_2d.plane_basis
        assert np.max(hull.equations[:, :2] @ screen + hull.equations[:, 2]) < 2e-8


def test_centered_rectangular_acute_path_points_stay_inside_bz():
    # Case 9's cell with the longer conventional side halved, which flips
    # a < b into a > b and so selects oC2 instead of oC1.
    lattice = np.array([
        [4.4450371550885732, 0.0, 1.5715129063682827],
        [0.0, 25.4549000000000021, 0.0],
        [-0.0336569615565514, 0.0, 4.7145387191048478],
    ])
    lattice_2d = analyze_lattice(lattice, 1)
    bz = build_bz(lattice_2d.reciprocal_2d)
    hull = ConvexHull(bz)
    path_data = build_path(lattice_2d, "mm2")

    assert lattice_2d.lattice_class == "centered_rectangular"
    assert lattice_2d.centered_branch == "acute"
    assert path_data["path"] == [
        (GAMMA_LABEL, "Y"),
        ("Y", "F_0"),
        ("DELTA_0", GAMMA_LABEL),
        (GAMMA_LABEL, "S"),
    ]
    for point in path_data["points"].values():
        cartesian = np.asarray(point) @ lattice_2d.reciprocal_3d
        screen = cartesian @ lattice_2d.plane_basis
        assert np.max(hull.equations[:, :2] @ screen + hull.equations[:, 2]) < 2e-8


def test_centered_rectangular_acute_ibz_is_a_quarter_of_the_bz():
    lattice = np.array([
        [4.4450371550885732, 0.0, 1.5715129063682827],
        [0.0, 25.4549000000000021, 0.0],
        [-0.0336569615565514, 0.0, 4.7145387191048478],
    ])
    lattice_2d = analyze_lattice(lattice, 1)
    path_data = build_path(lattice_2d, "mm2")
    polygon, _centroid, area, labels = build_ibz(lattice_2d, path_data, "mm2", 4)
    bz_area, _ = polygon_area_centroid_2d(build_bz(lattice_2d.reciprocal_2d))

    assert labels == ["DELTA_0", GAMMA_LABEL, "Y", "F_0"]
    assert area / bz_area == pytest.approx(0.25, abs=1e-12)


@pytest.mark.parametrize(
    ("lattice", "point_group", "group_order", "expected"),
    [
        (
            np.array([
                [4.0, 0.0, 0.0],
                [-2.0, 2.0 * np.sqrt(3.0), 0.0],
                [0.0, 0.0, 20.0],
            ]),
            "6mm", 12, [5.0 / 18.0, 1.0 / 9.0, 0.0],
        ),
        (np.diag([4.0, 4.0, 20.0]), "4mm", 8,
         [1.0 / 6.0, 1.0 / 3.0, 0.0]),
        (np.diag([4.0, 4.0, 20.0]), "4", 4, [0.25, 0.25, 0.0]),
        (np.diag([4.0, 6.0, 20.0]), "2mm", 4, [0.25, 0.25, 0.0]),
    ],
)
def test_2d_ibz_reproduces_established_easy_case_centroids(
    lattice, point_group, group_order, expected
):
    lattice_2d = analyze_lattice(lattice, 2)
    path_data = build_path(lattice_2d, point_group)
    _polygon, centroid, _area, _labels = build_ibz(
        lattice_2d, path_data, point_group, group_order
    )
    fractional = to_input_fractional(centroid, lattice_2d)[0]
    assert np.allclose(fractional, expected)


def test_conventional_centered_ibz_reproduces_submitted_basis_centroid():
    lattice = np.array([
        [5.1103299815503842, 0.0, 0.0],
        [-2.5057448085893013, 4.4538427761384352, 0.0],
        [0.0, 0.0, 22.4751439124],
    ])
    lattice_2d = analyze_lattice(lattice, 2)
    path_data = build_path(lattice_2d, "m")
    _polygon, centroid, _area, labels = build_ibz(
        lattice_2d, path_data, "m", 4
    )
    fractional = to_input_fractional(centroid, lattice_2d)[0]

    assert labels == ["SIGMA_0", GAMMA_LABEL, "Y", "C_0"]
    assert np.allclose(fractional, [-0.09162815, 0.35137431, 0.0])


def test_screen_to_fractional_restores_submitted_vacuum_axis():
    lattice = np.diag([20.0, 4.0, 6.0])
    lattice_2d = analyze_lattice(lattice, 0)
    point = 0.25 * lattice_2d.reciprocal_2d[0] + 0.5 * lattice_2d.reciprocal_2d[1]
    fractional = to_input_fractional(point, lattice_2d)[0]

    assert np.allclose(fractional, [0.0, 0.25, 0.5])


def test_rg3_end_to_end_uses_submitted_rectangular_cell(tmp_path):
    lattice = np.diag([4.0, 6.0, 20.0])
    positions = np.array([
        [0.17, 0.23, 0.11],
        [0.17, 0.77, 0.89],
        [0.83, 0.77, 0.89],
        [0.83, 0.23, 0.11],
    ])
    result = compute_centroid(
        "synthetic-rg3.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(lattice, positions, [1, 1, 1, 1]),
        analysis_marker_type=2,
    )

    assert result["path_source_2d"] is True
    assert result["sc_type"] == "rectangular"
    assert result["layer_group_symbol"] == "p2/m11"
    assert result["layer_group_number"] == 14
    assert result["layer_point_group"] == "2/m"
    assert result["band_kpath"] == [
        (GAMMA_LABEL, "X"),
        ("X", "S"),
        ("S", "Y"),
        ("Y", GAMMA_LABEL),
    ]
    assert np.allclose(result["b_matrix"], 2.0 * np.pi * np.linalg.inv(lattice).T)
    assert np.allclose(result["centroid_frac"], [0.25, 0.25, 0.0])
    expected_ibz_area = (2.0 * np.pi) ** 2 / (4.0 * 4.0 * 6.0)
    assert np.isclose(result["ibz_volume"], expected_ibz_area)
