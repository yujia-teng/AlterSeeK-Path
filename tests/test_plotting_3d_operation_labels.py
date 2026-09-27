"""Check the symbol and drawn turn of three-dimensional spin-flip operations."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from alterseek.plotting_3d import _draw_op_visual
from alterseek.symmetry import describe_spinflip_op


def _rotation_z(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


_BZ_LOOPS = [
    np.array([[-1, -1, z], [1, -1, z], [1, 1, z],
              [-1, 1, z], [-1, -1, z]], dtype=float)
    for z in (-1, 1)
]


@pytest.mark.parametrize("order,barred_order", [(3, 6), (4, 4), (6, 3)])
@pytest.mark.parametrize("turn", [1, -1], ids=["plus", "minus"])
@pytest.mark.parametrize("improper", [False, True], ids=["C", "S"])
def test_3d_rotation_label_and_arc_match_matrix(order, barred_order, turn, improper):
    angle = turn * 2 * np.pi / order
    operation = _rotation_z(angle)
    if improper:
        operation = operation @ np.diag([1.0, 1.0, -1.0])

    # The menu uses Schoenflies signs; the barred figure symbol uses the
    # opposite sign for precisely the same rotoinversion matrix.
    menu_symbol = "S" if improper else "C"
    menu_sign = "+" if turn > 0 else "-"
    assert describe_spinflip_op(operation, np.eye(3)) == (
        f"{menu_symbol}{order}{menu_sign} [0 0 1]"
    )
    figure_digit = rf"\bar{{{barred_order}}}" if improper else str(order)
    figure_sign = ("-" if turn > 0 else "+") if improper else menu_sign
    if improper:
        # Schoenflies Sn = horizontal reflection after Cn; the barred
        # International symbol is inversion after the opposite rotation.
        np.testing.assert_allclose(
            operation, -_rotation_z(-turn * 2 * np.pi / barred_order),
            atol=1e-14,
        )

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    try:
        _draw_op_visual(ax, operation, _BZ_LOOPS, 1.0, np.eye(3))
        labels = [artist.get_text() for artist in ax.texts if artist.get_text()]
        assert labels == [
            rf"$\mathbf{{{figure_digit}^{{{figure_sign}}}_{{001}}}}$"
        ]

        # The curved arrow shows the rotation before inversion in the
        # displayed barred symbol, rather than the opposite Schoenflies turn.
        arc = np.column_stack([ax.lines[-1].get_xdata(),
                               ax.lines[-1].get_ydata()])
        first, second = arc[:2]
        expected_arc_turn = -turn if improper else turn
        assert np.sign(np.linalg.det(np.array([first, second]))) == expected_arc_turn
    finally:
        plt.close(fig)


def test_3d_rotation_label_uses_reciprocal_axis_indices():
    # z = b1 + b3 in this reciprocal basis.
    b_matrix = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                         [-1.0, 0.0, 1.0]])
    operation = _rotation_z(2 * np.pi / 3) @ np.diag([1.0, 1.0, -1.0])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    try:
        _draw_op_visual(ax, operation, _BZ_LOOPS, 1.0, b_matrix)
        labels = [artist.get_text() for artist in ax.texts if artist.get_text()]
        assert labels == [r"$\mathbf{\bar{6}^{-}_{101}}$"]
    finally:
        plt.close(fig)


def test_gdauge_eighth_option_maps_k_to_the_shown_blue_image():
    # The real-space fractional matrix selected by option 8 in GdAuGe.
    operation_frac = np.array([[0, -1, 0], [1, -1, 0], [0, 0, -1]])
    k_frac = np.array([5 / 18, 1 / 9, 1 / 4])
    np.testing.assert_allclose(
        np.linalg.inv(operation_frac).T @ k_frac,
        [-7 / 18, 5 / 18, -1 / 4], atol=1e-14,
    )

    # A 60-degree reciprocal basis makes the same matrix Cartesian.
    b_matrix = np.array([[1.0, 0.0, 0.0],
                         [0.5, np.sqrt(3) / 2, 0.0],
                         [0.0, 0.0, 1.0]])
    operation = (b_matrix.T @ np.linalg.inv(operation_frac).T
                 @ np.linalg.inv(b_matrix.T))
    expected = _rotation_z(2 * np.pi / 3) @ np.diag([1.0, 1.0, -1.0])
    np.testing.assert_allclose(operation, expected, atol=1e-14)
    np.testing.assert_allclose(operation, -_rotation_z(-np.pi / 3), atol=1e-14)
    assert describe_spinflip_op(operation, b_matrix) == "S3+ [0 0 1]"

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    try:
        _draw_op_visual(ax, operation, _BZ_LOOPS, 1.0, b_matrix)
        assert [artist.get_text() for artist in ax.texts if artist.get_text()] == [
            r"$\mathbf{\bar{6}^{-}_{001}}$"
        ]
        arc = np.column_stack([ax.lines[-1].get_xdata(),
                               ax.lines[-1].get_ydata()])
        assert np.sign(np.linalg.det(arc[:2])) == -1
    finally:
        plt.close(fig)
