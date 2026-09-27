"""Check the International labels shared by the two-dimensional figures."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from alterseek.mode2d.plotting import _draw_op_visual_2d


@pytest.mark.parametrize("order,barred_order", [(3, 6), (4, 4), (6, 3)])
@pytest.mark.parametrize("turn", [1, -1], ids=["plus", "minus"])
@pytest.mark.parametrize("improper", [False, True], ids=["C", "S"])
def test_2d_rotation_label_and_arrow_match_international_symbol(
    order, barred_order, turn, improper
):
    angle = turn * 2 * np.pi / order
    c, s = np.cos(angle), np.sin(angle)
    operation = np.array([[c, -s, 0.0], [s, c, 0.0],
                          [0.0, 0.0, -1.0 if improper else 1.0]])
    basis = tuple(np.eye(3))
    bz_poly = np.array([[-1.0, -1.0], [1.0, -1.0],
                        [1.0, 1.0], [-1.0, 1.0]])
    digit = rf"\bar{{{barred_order}}}" if improper else str(order)
    symbol_turn = -turn if improper else turn
    sign = "+" if symbol_turn > 0 else "-"

    fig, ax = plt.subplots()
    try:
        _draw_op_visual_2d(ax, operation, np.eye(3), basis, bz_poly)
        labels = [artist.get_text() for artist in ax.texts if artist.get_text()]
        assert labels == [rf"$\mathbf{{{digit}^{{{sign}}}_{{001}}}}$"]

        arc = np.column_stack([ax.lines[0].get_xdata(),
                               ax.lines[0].get_ydata()])
        assert np.sign(np.linalg.det(arc[:2])) == symbol_turn
    finally:
        plt.close(fig)
