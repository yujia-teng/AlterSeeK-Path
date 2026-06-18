import matplotlib.pyplot as plt
import numpy as np

from plot_alterband import _draw_panel


def test_spin_up_band_is_drawn_above_spin_down_band():
    fig, ax = plt.subplots()
    try:
        _draw_panel(
            ax,
            labels=["GAMMA", "M"],
            positions=[0.0, 1.0],
            tick_labels=["Gamma", "M"],
            kpath=np.array([0.0, 0.5, 1.0]),
            bands_up=np.array([[0.0], [0.1], [0.0]]),
            bands_dw=np.array([[0.0], [0.1], [0.0]]),
            elim=(-1.0, 1.0),
            xlim=(0.0, 1.0),
            gap_half=0.01,
            font_size=12,
            rotate_xtick_labels=False,
            xtick_rotation=45.0,
        )

        spin_down = next(line for line in ax.lines if line.get_zorder() == 2)
        spin_up = next(line for line in ax.lines if line.get_zorder() == 3)
        assert spin_down.get_color() == "red"
        assert spin_down.get_linewidth() == 0.7
        assert spin_up.get_color() == "black"
        assert spin_up.get_linewidth() == 0.7
        assert spin_up.get_zorder() > spin_down.get_zorder()
    finally:
        plt.close(fig)
