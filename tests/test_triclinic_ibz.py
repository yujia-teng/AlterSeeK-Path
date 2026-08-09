"""Triclinic IBZ: the axis-containing half-BZ convention.

Laue group -1 makes the nonmagnetic IBZ half the Brillouin zone, and every
plane through Gamma cuts an admissible half, so the plane is a convention.
These tests pin the convention actually chosen and the numbers it produces.
"""

from fractions import Fraction

import numpy as np
import pytest

from alterseek.compute_centroid_hybrid import _compute_ibz_centroid, _format_plane
from alterseek.geometry import (
    calculate_volume_centroid,
    triclinic_half_bz_cell,
    triclinic_halfspace_normal,
)
from alterseek.lattice_kpoints import LATTICE_DATA

# Real triclinic MnO2, standardized primitive reciprocal basis (aP2).
MNO2_B_MATRIX = np.array([
    [-0.00000000, -0.00000000, 1.31420282],
    [-2.15197256, 1.24268350, -0.40111844],
    [-0.00000000, -2.48464112, -0.45654232],
])
MNO2_CENTROID_FRAC = np.array([0.147777, 0.222625, 0.222704])

# Real all-acute MoS3, the aP3 counterpart of the all-obtuse MnO2 case
# (`1-triclinic/MoS3.poscar`).  Unlike MnO2 it is space group P1 rather than
# P-1, so it also covers the triclinic group that has no inversion of its own
# and reaches the same half-BZ domain only through Laue -1.
AP3_B_MATRIX = np.array([
    [-1.62151525, -0.13499029, -0.10486344],
    [-0.00000000, -0.00000000, -0.58331336],
    [-0.00000000, -0.94455210, -0.15834004],
])
AP3_CENTROID_FRAC = np.array([0.053920, -0.078893, 0.230130])


def _selected_points(lattice_key):
    data = LATTICE_DATA[lattice_key]
    points = {}
    for start, end in data["kpath"]:
        for label in (start, end):
            coords = np.array(
                [float(Fraction(v)) for v in data["points_def"][label]]
            )
            if not np.allclose(coords, 0.0):
                points[label] = coords
    return points


@pytest.mark.parametrize(
    "lattice_key, expected_normal, expected_plane, boundary_label",
    [
        ("aP2", (0, 1, 1), "k2+k3=0", "X"),
        ("aP3", (1, 0, 2), "k1+2*k3=0", "Y"),
    ],
)
def test_axis_containing_normal_is_the_recorded_candidate(
    lattice_key, expected_normal, expected_plane, boundary_label
):
    normal = triclinic_halfspace_normal(lattice_key)

    assert normal == expected_normal
    assert _format_plane(normal) == expected_plane

    points = _selected_points(lattice_key)
    values = {
        label: float(np.dot(normal, point)) for label, point in points.items()
    }
    # Exactly the chosen axis's own path point lies on the cut plane; every
    # other selected path point stays strictly inside the positive half.
    assert values[boundary_label] == pytest.approx(0.0, abs=1e-12)
    assert all(
        value > 1e-12
        for label, value in values.items()
        if label != boundary_label
    )


def test_ap3_admits_no_b1_containing_plane():
    """Y and V_2 need opposite signs once the first component is zero."""
    points = _selected_points("aP3")
    # A b1-containing plane always leaves X on the boundary, so the binding
    # requirement is that every *other* selected point be strictly positive.
    others = {label: p for label, p in points.items() if label != "X"}

    for second in range(-4, 5):
        for third in range(-4, 5):
            normal = (0, second, third)
            if not any(normal):
                continue
            assert not all(
                float(np.dot(normal, point)) > 1e-12
                for point in others.values()
            )

    # The specific recorded obstruction: Y and V_2 can never agree in sign.
    for second in (-3, -1, 1, 2, 4):
        y_value = float(np.dot((0, second, 0), points["Y"]))
        v2_value = float(np.dot((0, second, 0), points["V_2"]))
        assert y_value * v2_value < 0.0


def test_selected_normal_is_the_unique_simplest_admissible_one():
    """The admissible family is not equivalent, so the tie-break matters."""
    points = _selected_points("aP2")
    normal = triclinic_halfspace_normal("aP2")
    weight = sum(abs(component) for component in normal)

    competing = [
        candidate
        for second in range(0, 5)
        for third in range(0, 5)
        for candidate in [(0, second, third)]
        if any(candidate)
        and sum(candidate) <= weight
        and all(
            float(np.dot(candidate, point)) > 1e-12
            for label, point in points.items()
            if label != "X"
        )
    ]

    assert competing == [normal]


def test_half_bz_cell_is_exactly_half_with_the_recorded_centroid():
    normal = triclinic_halfspace_normal("aP2")
    _, hull = triclinic_half_bz_cell(MNO2_B_MATRIX, normal)
    centroid_cart, volume = calculate_volume_centroid(hull)
    centroid_frac = centroid_cart @ np.linalg.inv(MNO2_B_MATRIX)

    assert volume == pytest.approx(
        abs(np.linalg.det(MNO2_B_MATRIX)) / 2.0, rel=1e-9
    )
    assert centroid_frac == pytest.approx(MNO2_CENTROID_FRAC, abs=1e-6)


def test_ap3_half_bz_cell_is_exactly_half_and_holds_its_path_points():
    """The all-acute branch, whose plane contains b2 rather than b1.

    MoS3 is P1, so this is also the real-material check that a triclinic cell
    without its own inversion still reaches the Laue -1 half-BZ domain.
    """
    normal = triclinic_halfspace_normal("aP3")
    vertices, hull = triclinic_half_bz_cell(AP3_B_MATRIX, normal)
    centroid_cart, volume = calculate_volume_centroid(hull)
    centroid_frac = centroid_cart @ np.linalg.inv(AP3_B_MATRIX)
    cart_normal = np.linalg.inv(AP3_B_MATRIX) @ np.array(normal, dtype=float)

    assert volume == pytest.approx(
        abs(np.linalg.det(AP3_B_MATRIX)) / 2.0, rel=1e-9
    )
    assert centroid_frac == pytest.approx(AP3_CENTROID_FRAC, abs=1e-6)
    assert np.all(vertices @ cart_normal > -1e-9)
    # The centroid is a genuine interior point of the chosen half.
    assert float(np.dot(normal, centroid_frac)) > 1e-6


def test_half_bz_cell_keeps_every_selected_path_point_inside():
    normal = triclinic_halfspace_normal("aP2")
    vertices, _ = triclinic_half_bz_cell(MNO2_B_MATRIX, normal)
    cart_normal = np.linalg.inv(MNO2_B_MATRIX) @ np.array(normal, dtype=float)

    assert np.all(vertices @ cart_normal > -1e-9)
    for point in _selected_points("aP2").values():
        assert float(np.dot(normal, point)) >= -1e-12


def test_degenerate_normal_is_rejected():
    with pytest.raises(ValueError):
        triclinic_half_bz_cell(MNO2_B_MATRIX, (0, 0, 0))


def _triclinic_centroid_inputs():
    frac = {
        label: np.array([float(Fraction(v)) for v in coords])
        for label, coords in LATTICE_DATA["aP2"]["points_def"].items()
    }
    cart = {label: point @ MNO2_B_MATRIX for label, point in frac.items()}
    return dict(
        mode_2d=False,
        vacuum_axis=2,
        sc_type="aP2",
        b_matrix=MNO2_B_MATRIX,
        unique_ops=[],
        kpoints_frac_centroid=frac,
        kpoints_cart_centroid=cart,
        verbose=False,
    )


@pytest.mark.parametrize("sg", [1, 2])
def test_workflow_triclinic_branch_uses_the_half_bz_for_p1_and_p_minus_1(sg):
    """Laue -1 covers both triclinic groups, so both get the same domain."""
    result = _compute_ibz_centroid(sg=sg, **_triclinic_centroid_inputs())

    assert result["hull"] is not None
    assert result["ibz_volume"] == pytest.approx(
        abs(np.linalg.det(MNO2_B_MATRIX)) / 2.0, rel=1e-9
    )
    assert result["centroid_frac"] == pytest.approx(
        MNO2_CENTROID_FRAC, abs=1e-6
    )
    # The clipped hull's vertices are Voronoi intersections, so they carry no
    # HPKOT labels and must not be reported as if they did.
    assert result["labels_list"] == []
    assert result["hull_matches_labels"] is False

    # The superseded behaviour was the arithmetic mean of the curated points
    # with no hull and zero volume.
    old_mean = np.mean(
        np.array(list(_triclinic_centroid_inputs()["kpoints_cart_centroid"]
                      .values())),
        axis=0,
    ) @ np.linalg.inv(MNO2_B_MATRIX)
    assert not np.allclose(result["centroid_frac"], old_mean, atol=1e-4)


def test_workflow_triclinic_branch_falls_back_when_the_clip_fails(
    monkeypatch, capsys
):
    """A geometry failure must warn, not cost the run its centroid."""
    def _boom(*_args, **_kwargs):
        raise RuntimeError("synthetic clip failure")

    monkeypatch.setattr(
        "alterseek.compute_centroid_hybrid.triclinic_half_bz_cell", _boom
    )
    result = _compute_ibz_centroid(sg=2, **_triclinic_centroid_inputs())

    assert result["hull"] is None
    assert result["ibz_volume"] == 0.0
    assert result["centroid_frac"] is not None
    assert "synthetic clip failure" in capsys.readouterr().out


def test_format_plane_renders_signs_and_coefficients():
    assert _format_plane((0, 1, 1)) == "k2+k3=0"
    assert _format_plane((1, 0, 2)) == "k1+2*k3=0"
    assert _format_plane((-1, 0, 3)) == "-k1+3*k3=0"
