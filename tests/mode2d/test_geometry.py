import numpy as np
import pytest
from scipy.spatial import ConvexHull

from alterseek.compute_centroid_3d import run as compute_centroid
from alterseek.mode2d.geometry import (
    analyze_lattice,
    build_bz,
    polygon_area_centroid_2d,
    to_input_fractional,
)
from alterseek.mode2d.lattice_kpoints import build_ibz, build_path
from alterseek.mode2d.symmetry import project_point_operations
from alterseek.plotting_common import GAMMA_LABEL
from alterseek.symmetry import no_altermagnetism_reason


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


def _oblique_lattice():
    return np.array([
        [4.0, 0.0, 0.0],
        [1.3, 5.2, 0.0],
        [0.0, 0.0, 20.0],
    ])


@pytest.mark.parametrize("point_group", ["1", "2"])
def test_oblique_groups_use_bilbao_labels_and_three_gamma_rays(point_group):
    lattice_2d = analyze_lattice(_oblique_lattice(), 2)
    path_data = build_path(lattice_2d, point_group, projected_group_order=2)
    bz = build_bz(lattice_2d.reciprocal_2d)
    hull = ConvexHull(bz)

    assert set(path_data["points"]) == {GAMMA_LABEL, "Y", "B", "A"}
    assert path_data["path"] == [
        (GAMMA_LABEL, "B"), ("Y", GAMMA_LABEL), (GAMMA_LABEL, "A")
    ]
    assert path_data["extra_vertices"] == []
    for point in path_data["points"].values():
        cartesian = np.asarray(point) @ lattice_2d.reciprocal_3d
        screen = cartesian @ lattice_2d.plane_basis
        assert np.max(
            hull.equations[:, :2] @ screen + hull.equations[:, 2]
        ) < 2e-8


@pytest.mark.parametrize("point_group", ["1", "2"])
def test_oblique_groups_label_every_half_wigner_seitz_vertex(point_group):
    lattice_2d = analyze_lattice(_oblique_lattice(), 2)
    path_data = build_path(lattice_2d, point_group, projected_group_order=2)
    polygon, centroid, area, labels = build_ibz(
        lattice_2d, path_data, point_group, 2
    )
    bz_area, _ = polygon_area_centroid_2d(
        build_bz(lattice_2d.reciprocal_2d)
    )
    centroid_frac = to_input_fractional(centroid, lattice_2d)[0]

    assert len(polygon) >= 4
    assert labels == ["Q", "Q_A", "Y", "Y_A", "Q_B"]
    assert path_data["extra_vertices"] == ["Q", "Q_A", "Y_A", "Q_B"]
    assert set(path_data["points"]) == {
        GAMMA_LABEL, "Y", "B", "A", "Q", "Q_A", "Y_A", "Q_B"
    }
    assert area / bz_area == pytest.approx(0.5, abs=1e-12)
    assert centroid_frac[0] > 0.0


@pytest.mark.parametrize(
    ("case", "positions", "types", "point_group", "layer_symbol"),
    [
        (
            "p1",
            [[0.11, 0.17, 0.21], [0.31, 0.27, 0.39]],
            [1, 2],
            "1",
            "p1",
        ),
        (
            "p2",
            [
                [0.11, 0.17, 0.21], [0.89, 0.83, 0.21],
                [0.23, 0.07, 0.37], [0.77, 0.93, 0.37],
            ],
            [1, 1, 2, 2],
            "2",
            "p112",
        ),
    ],
)
def test_oblique_groups_end_to_end_use_bilbao_path(
    tmp_path, case, positions, types, point_group, layer_symbol
):
    result = compute_centroid(
        f"synthetic-{case}.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(_oblique_lattice(), positions, types),
        analysis_has_markers=True,
        symprec=1e-5,
    )

    assert result["point_group"] == point_group
    assert result["layer_group_symbol"] == layer_symbol
    assert len(result["projected_point_operations_2d"]) == 2
    assert set(result["band_kpoints_frac"]) == {
        GAMMA_LABEL, "Y", "B", "A", "Q", "Q_A", "Y_A", "Q_B"
    }
    assert result["band_kpath"] == [
        (GAMMA_LABEL, "B"), ("Y", GAMMA_LABEL), (GAMMA_LABEL, "A")
    ]
    assert result["extra_general_vertices"] == ["Q", "Q_A", "Y_A", "Q_B"]
    assert result["ibz_polygon_labels"] == [
        "Q", "Q_A", "Y", "Y_A", "Q_B"
    ]
    assert result["color_copied_ibz_sectors"] is False
    if point_group == "1":
        assert result["no_altermagnetism"] == {
            "laue_group": "-1",
            "reason": "No altermagnetism",
        }


@pytest.mark.parametrize(
    "direct",
    [
        [[3.2, 0.0], [1.1, 4.7]],
        [[3.2, 0.0], [-1.1, 4.7]],
        [[4.0, 0.0], [0.02, 5.0]],
        [[4.0, 0.0], [1.2, 3.868]],
        [[2.5, 0.0], [1.15, 8.0]],
    ],
)
def test_oblique_half_bz_is_stable_across_metrics(direct):
    direct = np.asarray(direct, dtype=float)
    lattice = np.array([
        [direct[0, 0], direct[0, 1], 0.0],
        [direct[1, 0], direct[1, 1], 0.0],
        [0.0, 0.0, 20.0],
    ])
    lattice_2d = analyze_lattice(lattice, 2)
    path_data = build_path(lattice_2d, "2", projected_group_order=2)
    polygon, _centroid, area, labels = build_ibz(
        lattice_2d, path_data, "2", 2
    )
    bz_area, _ = polygon_area_centroid_2d(
        build_bz(lattice_2d.reciprocal_2d)
    )

    assert lattice_2d.lattice_class == "oblique"
    assert len(polygon) == 5
    assert labels == ["Q", "Q_A", "Y", "Y_A", "Q_B"]
    assert area / bz_area == pytest.approx(0.5, abs=1e-12)


def test_oblique_bilbao_axes_are_the_submitted_reciprocal_axes():
    lattice_2d = analyze_lattice(_oblique_lattice(), 2)
    path_data = build_path(lattice_2d, "1", projected_group_order=2)

    assert np.allclose(path_data["points"]["B"], [0.5, 0.0, 0.0])
    assert np.allclose(path_data["points"]["Y"], [0.0, 0.5, 0.0])


@pytest.mark.parametrize("vacuum_axis", [0, 1, 2])
def test_oblique_convention_is_independent_of_vacuum_axis(vacuum_axis):
    direct = _oblique_lattice()[:2, :2]
    axes = [index for index in range(3) if index != vacuum_axis]
    lattice = np.zeros((3, 3))
    lattice[vacuum_axis, vacuum_axis] = 20.0
    lattice[axes[0], axes] = direct[0]
    lattice[axes[1], axes] = direct[1]
    lattice_2d = analyze_lattice(lattice, vacuum_axis)
    path_data = build_path(lattice_2d, "1", projected_group_order=2)
    _polygon, _centroid, _area, labels = build_ibz(
        lattice_2d, path_data, "1", 2
    )

    assert labels == ["Q", "Q_A", "Y", "Y_A", "Q_B"]
    assert all(
        abs(np.asarray(point)[vacuum_axis]) < 1e-12
        for point in path_data["points"].values()
    )


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
    # Case 9's cell with the longer conventional side halved, which flips the
    # centered-rectangular metric branch from a < b to a > b.
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


def _hexagonal_lattice():
    return np.array([
        [4.0, 0.0, 0.0],
        [-2.0, 2.0 * np.sqrt(3.0), 0.0],
        [0.0, 0.0, 20.0],
    ])


def _reflection_about(vector):
    direction = np.asarray(vector, dtype=float)
    direction /= np.linalg.norm(direction)
    return 2.0 * np.outer(direction, direction) - np.eye(2)


def _hexagonal_2mm_operations(lattice_2d):
    first, second = lattice_2d.direct_2d
    return [
        np.eye(2),
        -np.eye(2),
        _reflection_about(first + second),
        _reflection_about(first - second),
    ]


def _assert_tiles_bz(polygon, bz, operations):
    """Sample-test that symmetry images cover the BZ once off boundaries."""
    images = [np.asarray(polygon) @ np.asarray(operation).T for operation in operations]
    image_hulls = [ConvexHull(image) for image in images]
    bz_hull = ConvexHull(bz)
    rng = np.random.default_rng(20260902)
    lower = np.min(bz, axis=0)
    upper = np.max(bz, axis=0)
    checked = 0
    while checked < 1000:
        point = rng.uniform(lower, upper)
        if np.max(bz_hull.equations[:, :2] @ point + bz_hull.equations[:, 2]) >= -1e-7:
            continue
        memberships = [
            np.max(hull.equations[:, :2] @ point + hull.equations[:, 2]) < -1e-7
            for hull in image_hulls
        ]
        assert sum(memberships) == 1
        checked += 1


def test_hexagonal_2mm_uses_centered_path_in_physical_hexagonal_bz():
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    operations = _hexagonal_2mm_operations(lattice_2d)
    path_data = build_path(
        lattice_2d,
        "2mm",
        projected_group_order=4,
        projected_operations=operations,
    )
    polygon, centroid, area, labels = build_ibz(
        lattice_2d, path_data, "2mm", 4
    )
    bz = build_bz(lattice_2d.reciprocal_2d)
    bz_area, _ = polygon_area_centroid_2d(bz)

    assert lattice_2d.lattice_class == "hexagonal"
    assert path_data["path_lattice_class"] == "centered_rectangular"
    assert path_data["path_metric_branch"] == "a_lt_b"
    assert np.array_equal(
        path_data["path_lattice"].canonical_transform, np.eye(2, dtype=int)
    )
    assert path_data["path"] == [
        (GAMMA_LABEL, "Y"),
        ("Y", "C_0"),
        ("SIGMA_0", GAMMA_LABEL),
        (GAMMA_LABEL, "S"),
    ]
    assert labels == ["SIGMA_0", GAMMA_LABEL, "Y", "C_0"]
    assert area / bz_area == pytest.approx(0.25, abs=1e-12)
    assert np.allclose(
        to_input_fractional(centroid, lattice_2d)[0],
        [-5.0 / 54.0, 19.0 / 54.0, 0.0],
    )
    _assert_tiles_bz(polygon, bz, operations)


def test_hexagonal_2mm_other_orientation_keeps_centered_a_lt_b_convention():
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    operations = [
        np.eye(2), -np.eye(2),
        np.diag([1.0, -1.0]), np.diag([-1.0, 1.0]),
    ]
    path_data = build_path(
        lattice_2d,
        "2mm",
        projected_group_order=4,
        projected_operations=operations,
    )
    polygon, _centroid, area, labels = build_ibz(
        lattice_2d, path_data, "2mm", 4
    )
    bz = build_bz(lattice_2d.reciprocal_2d)
    bz_area, _ = polygon_area_centroid_2d(bz)

    assert path_data["path_lattice_class"] == "centered_rectangular"
    assert path_data["path_metric_branch"] == "a_lt_b"
    assert labels == ["SIGMA_0", GAMMA_LABEL, "Y", "C_0"]
    assert area / bz_area == pytest.approx(0.25, abs=1e-12)
    _assert_tiles_bz(polygon, bz, operations)


def test_hexagonal_centered_labels_follow_conventional_a_under_basis_rewrites():
    direct = np.array([
        [2.0, -2.0 * np.sqrt(3.0)],
        [2.0, 2.0 * np.sqrt(3.0)],
    ])
    operations = [
        np.eye(2), -np.eye(2),
        np.diag([1.0, -1.0]), np.diag([-1.0, 1.0]),
    ]
    rewrites = [
        np.eye(2, dtype=int),
        np.array([[0, 1], [1, 0]], dtype=int),
        np.array([[1, 1], [0, 1]], dtype=int),
        np.array([[-1, 0], [0, 1]], dtype=int),
        np.array([[1, 0], [1, 1]], dtype=int),
    ]
    reference_points = None
    reference_polygon = None
    for rewrite in rewrites:
        lattice = np.zeros((3, 3))
        lattice[:2, :2] = rewrite @ direct
        lattice[2, 2] = 20.0
        lattice_2d = analyze_lattice(lattice, 2)
        path_data = build_path(
            lattice_2d,
            "2mm",
            projected_group_order=4,
            projected_operations=operations,
        )
        polygon, _centroid, _area, labels = build_ibz(
            lattice_2d, path_data, "2mm", 4
        )
        canonical = path_data["path_lattice"].canonical_direct_2d
        conventional_a = canonical[0] + canonical[1]
        conventional_b = canonical[1] - canonical[0]
        points = {
            label: np.asarray(coords) @ lattice_2d.reciprocal_3d
            @ lattice_2d.plane_basis
            for label, coords in path_data["points"].items()
        }

        assert path_data["path_metric_branch"] == "a_lt_b"
        assert np.allclose(conventional_a[1], 0.0, atol=1e-10)
        assert conventional_a[0] > 0.0
        assert np.allclose(conventional_b[0], 0.0, atol=1e-10)
        assert conventional_b[1] > 0.0
        assert np.linalg.norm(conventional_a) < np.linalg.norm(conventional_b)
        assert labels == ["SIGMA_0", GAMMA_LABEL, "Y", "C_0"]
        assert np.allclose(points["Y"], [0.0, np.pi / (2.0 * np.sqrt(3.0))])
        if reference_points is None:
            reference_points = points
            reference_polygon = polygon
        else:
            assert points.keys() == reference_points.keys()
            for label in points:
                assert np.allclose(points[label], reference_points[label])
            assert np.allclose(polygon, reference_polygon)


@pytest.mark.parametrize(
    ("lattice", "physical_class"),
    [
        (np.diag([4.0, 6.0, 20.0]), "rectangular"),
        (np.diag([4.0, 4.0, 20.0]), "square"),
        (_hexagonal_lattice(), "hexagonal"),
        (
            np.array([
                [5.1103299815503842, 0.0, 0.0],
                [-2.5057448085893013, 4.4538427761384352, 0.0],
                [0.0, 0.0, 22.4751439124],
            ]),
            "centered_rectangular",
        ),
    ],
)
def test_order_two_uses_oblique_path_without_changing_physical_bz(
    lattice, physical_class
):
    lattice_2d = analyze_lattice(lattice, 2)
    operations = [np.eye(2), -np.eye(2)]
    path_data = build_path(
        lattice_2d,
        "2",
        projected_group_order=2,
        projected_operations=operations,
    )
    polygon, _centroid, area, labels = build_ibz(
        lattice_2d, path_data, "2", 2
    )
    bz = build_bz(lattice_2d.reciprocal_2d)
    bz_area, _ = polygon_area_centroid_2d(bz)

    assert lattice_2d.lattice_class == physical_class
    assert path_data["path_lattice_class"] == "oblique"
    assert path_data["path"] == [
        (GAMMA_LABEL, "B"),
        ("Y", GAMMA_LABEL),
        (GAMMA_LABEL, "A"),
    ]
    assert labels[0] == "Q"
    assert "Y" in labels and "Y_A" in labels
    assert area / bz_area == pytest.approx(0.5, abs=1e-12)
    _assert_tiles_bz(polygon, bz, operations)


@pytest.mark.parametrize(
    ("operations", "path_class"),
    [
        (
            [
                np.eye(2), -np.eye(2),
                np.diag([1.0, -1.0]), np.diag([-1.0, 1.0]),
            ],
            "rectangular",
        ),
        (
            [
                np.eye(2), -np.eye(2),
                np.array([[0.0, 1.0], [1.0, 0.0]]),
                np.array([[0.0, -1.0], [-1.0, 0.0]]),
            ],
            "centered_rectangular",
        ),
    ],
)
def test_square_2mm_orientation_selects_primitive_or_centered_path(
    operations, path_class
):
    lattice_2d = analyze_lattice(np.diag([4.0, 4.0, 20.0]), 2)
    path_data = build_path(
        lattice_2d,
        "2mm",
        projected_group_order=4,
        projected_operations=operations,
    )
    polygon, _centroid, area, _labels = build_ibz(
        lattice_2d, path_data, "2mm", 4
    )
    bz = build_bz(lattice_2d.reciprocal_2d)
    bz_area, _ = polygon_area_centroid_2d(bz)

    assert lattice_2d.lattice_class == "square"
    assert path_data["path_lattice_class"] == path_class
    assert area / bz_area == pytest.approx(0.25, abs=1e-12)
    _assert_tiles_bz(polygon, bz, operations)


def test_hexagonal_cmm2_end_to_end_keeps_hexagonal_bz_and_uses_oc_labels(
    tmp_path,
):
    lattice = _hexagonal_lattice()
    fractional_operations = [
        np.eye(2, dtype=int),
        -np.eye(2, dtype=int),
        np.array([[0, 1], [1, 0]], dtype=int),
        np.array([[0, -1], [-1, 0]], dtype=int),
    ]
    positions = []
    types = []
    for seed, orbit_type, height in (
        (np.array([0.113, 0.227]), 1, 0.21),
        (np.array([0.071, 0.319]), 2, 0.37),
    ):
        for operation in fractional_operations:
            point = np.mod(operation @ seed, 1.0)
            positions.append([point[0], point[1], height])
            types.append(orbit_type)

    result = compute_centroid(
        "synthetic-cmm2.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(lattice, positions, types),
        analysis_has_markers=True,
        symprec=1e-5,
    )

    assert result["point_group"] == "mm2"
    assert result["layer_group_symbol"] == "cmm2"
    assert len(result["projected_point_operations_2d"]) == 4
    assert result["lattice_class_2d"] == "hexagonal"
    assert result["path_lattice_class_2d"] == "centered_rectangular"
    assert result["path_metric_branch_2d"] == "a_lt_b"
    assert result["ibz_polygon_labels"] == [
        "SIGMA_0", GAMMA_LABEL, "Y", "C_0"
    ]
    assert result["band_kpath"] == [
        (GAMMA_LABEL, "Y"),
        ("Y", "C_0"),
        ("SIGMA_0", GAMMA_LABEL),
        (GAMMA_LABEL, "S"),
    ]
    assert np.allclose(
        result["centroid_frac"], [-5.0 / 54.0, 19.0 / 54.0, 0.0]
    )


def test_hexagonal_p1_end_to_end_uses_oblique_fallback_convention(tmp_path):
    lattice = _hexagonal_lattice()
    positions = [
        [0.113, 0.227, 0.21],
        [0.071, 0.319, 0.37],
        [0.143, 0.091, 0.44],
    ]
    result = compute_centroid(
        "synthetic-p1-on-hexagonal-lattice.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(lattice, positions, [1, 2, 3]),
        analysis_has_markers=True,
        symprec=1e-5,
    )

    assert result["point_group"] == "1"
    assert result["layer_group_symbol"] == "p1"
    assert len(result["projected_point_operations_2d"]) == 2
    assert result["lattice_class_2d"] == "hexagonal"
    assert result["path_lattice_class_2d"] == "oblique"
    assert result["band_kpath"] == [
        (GAMMA_LABEL, "B"),
        ("Y", GAMMA_LABEL),
        (GAMMA_LABEL, "A"),
    ]
    assert result["no_altermagnetism"] == {
        "laue_group": "-1",
        "reason": "No altermagnetism",
    }


def test_hexagonal_m_end_to_end_uses_same_centered_domain_as_2mm(tmp_path):
    lattice = _hexagonal_lattice()
    mirror = np.array([[0, 1], [1, 0]], dtype=int)
    positions = []
    types = []
    for seed, orbit_type, height in (
        (np.array([0.113, 0.227]), 1, 0.21),
        (np.array([0.071, 0.319]), 2, 0.37),
        (np.array([0.143, 0.091]), 3, 0.44),
    ):
        for operation in (np.eye(2, dtype=int), mirror):
            point = np.mod(operation @ seed, 1.0)
            positions.append([point[0], point[1], height])
            types.append(orbit_type)

    result = compute_centroid(
        "synthetic-cm-on-hexagonal-lattice.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(lattice, positions, types),
        analysis_has_markers=True,
        symprec=1e-5,
    )

    assert result["point_group"] == "m"
    assert result["layer_group_symbol"] == "cm11"
    assert len(result["projected_point_operations_2d"]) == 4
    assert result["lattice_class_2d"] == "hexagonal"
    assert result["path_lattice_class_2d"] == "centered_rectangular"
    assert result["path_metric_branch_2d"] == "a_lt_b"
    assert result["ibz_polygon_labels"] == [
        "SIGMA_0", GAMMA_LABEL, "Y", "C_0"
    ]


def test_threefold_inversion_extension_has_six_in_plane_k_operations():
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    threefold = np.array([
        [0, -1, 0],
        [1, -1, 0],
        [0, 0, 1],
    ])
    rotations = [
        np.eye(3, dtype=int),
        threefold,
        threefold @ threefold,
    ]

    operations = project_point_operations(
        lattice_2d, rotations, add_inversion=True
    )

    assert len(operations) == 6
    assert no_altermagnetism_reason("3") == {
        "laue_group": "-3",
        "reason": "No altermagnetism",
    }


@pytest.mark.parametrize(
    "mirror",
    [
        np.array([[1, -1, 0], [0, -1, 0], [0, 0, 1]], dtype=int),
        np.array([[1, 0, 0], [1, -1, 0], [0, 0, 1]], dtype=int),
    ],
)
def test_both_oriented_threefold_mirror_groups_extend_to_same_hexagonal_k_group(
    mirror,
):
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    threefold = np.array([
        [0, -1, 0],
        [1, -1, 0],
        [0, 0, 1],
    ])
    rotations = []
    rotation = np.eye(3, dtype=int)
    for _ in range(3):
        rotations.extend([rotation, rotation @ mirror])
        rotation = rotation @ threefold

    operations = project_point_operations(
        lattice_2d, rotations, add_inversion=True
    )
    path_data = build_path(
        lattice_2d,
        "3m",
        projected_group_order=len(operations),
        projected_operations=operations,
    )
    polygon, _centroid, area, labels = build_ibz(
        lattice_2d, path_data, "3m", len(operations)
    )
    bz_area, _ = polygon_area_centroid_2d(
        build_bz(lattice_2d.reciprocal_2d)
    )

    assert len(operations) == 12
    assert labels == [GAMMA_LABEL, "M", "K"]
    assert area / bz_area == pytest.approx(1.0 / 12.0, abs=1e-12)


@pytest.mark.parametrize("point_group", ["3", "6"])
def test_threefold_and_sixfold_use_the_same_enlarged_hexagonal_ibz(point_group):
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    path_data = build_path(
        lattice_2d, point_group, projected_group_order=6
    )
    polygon, centroid, area, labels = build_ibz(
        lattice_2d, path_data, point_group, 6
    )
    bz_area, _ = polygon_area_centroid_2d(
        build_bz(lattice_2d.reciprocal_2d)
    )
    centroid_frac = to_input_fractional(centroid, lattice_2d)[0]

    # Gamma-K stays in the reference path.  The centroid sits on it, so the
    # general-k filter in kpoints.py drops it there, not here.
    assert path_data["path"] == [
        (GAMMA_LABEL, "M"), ("M", "K"), ("K", GAMMA_LABEL)
    ]
    assert path_data["butterfly_path"] == path_data["path"]
    assert path_data["butterfly_extra_vertices"] == [
        GAMMA_LABEL, "M_A"
    ]
    assert path_data["extra_vertices"] == ["M_A"]
    assert np.allclose(path_data["points"]["M_A"], [0.0, 0.5, 0.0])
    assert labels == [GAMMA_LABEL, "M", "K", "M_A"]
    assert len(polygon) == 4
    assert area / bz_area == pytest.approx(1.0 / 6.0, abs=1e-12)
    assert np.allclose(centroid_frac, [7.0 / 36.0, 7.0 / 36.0, 0.0])


def test_threefold_mirror_group_keeps_the_full_hexagonal_wedge():
    lattice_2d = analyze_lattice(_hexagonal_lattice(), 2)
    path_data = build_path(lattice_2d, "3m", projected_group_order=12)
    polygon, _centroid, area, labels = build_ibz(
        lattice_2d, path_data, "3m", 12
    )
    bz_area, _ = polygon_area_centroid_2d(
        build_bz(lattice_2d.reciprocal_2d)
    )

    assert path_data["extra_vertices"] == []
    assert labels == [GAMMA_LABEL, "M", "K"]
    assert len(polygon) == 3
    assert area / bz_area == pytest.approx(1.0 / 12.0, abs=1e-12)


def test_p3_end_to_end_uses_enlarged_ibz_before_non_altermagnetic_fallback(
    tmp_path,
):
    lattice = _hexagonal_lattice()
    threefold = np.array([
        [0, -1, 0],
        [1, -1, 0],
        [0, 0, 1],
    ])
    positions = []
    types = []
    for seed, orbit_type in (
        (np.array([0.11, 0.17, 0.21]), 1),
        (np.array([0.23, 0.07, 0.37]), 2),
    ):
        point = seed.copy()
        for _ in range(3):
            positions.append(np.mod(point, 1.0))
            types.append(orbit_type)
            point = point @ threefold.T

    result = compute_centroid(
        "synthetic-p3.vasp",
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        mode_2d=True,
        input_vacuum_axis=2,
        analysis_cell=(lattice, positions, types),
        analysis_has_markers=True,
    )

    assert result["point_group"] == "3"
    assert result["layer_group_symbol"] == "p3"
    assert len(result["projected_point_operations_2d"]) == 6
    assert result["ibz_polygon_labels"] == [
        GAMMA_LABEL, "M", "K", "M_A"
    ]
    assert result["extra_general_vertices"] == ["M_A"]
    assert np.allclose(
        result["centroid_frac"], [7.0 / 36.0, 7.0 / 36.0, 0.0]
    )
    assert result["no_altermagnetism"] == {
        "laue_group": "-3",
        "reason": "No altermagnetism",
    }


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
        analysis_has_markers=True,
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
