import numpy as np
import matplotlib.pyplot as plt

from alterseek.lattice_kpoints import get_hull_kpath, get_hull_kpoints
from alterseek.plotting_3d import plot_ibz
from alterseek.plotting_common import (
    _math_label,
    grouped_point_labels,
    prime_point_label,
)


def test_ervo3_metric_keeps_only_path_labels_at_coincident_vertices():
    points = get_hull_kpoints(
        "mP1", 7.5716, 5.2644, 5.5856, 90.0,
        spacegroup_number=14,
    )
    path = get_hull_kpath("mP1", spacegroup_number=14)
    path_labels = {label for segment in path for label in segment}

    displayed = {
        label: coords
        for label, coords in grouped_point_labels(points, path_labels)
    }

    assert "A" in displayed
    assert "E" in displayed
    assert not {"H", "H_2", "M", "M_2"} & displayed.keys()
    assert all("/" not in label for label in displayed)


def test_coincident_labels_used_by_path_are_slash_combined():
    points = {
        "A": np.array([0.0, 0.0, 0.0]),
        "H": np.array([0.0, 0.0, 0.0]),
        "X": np.array([0.5, 0.0, 0.0]),
    }

    displayed = grouped_point_labels(points, {"A", "H", "X"})

    assert [label for label, _coords in displayed] == ["A/H", "X"]


def test_3d_ibz_draws_one_marker_and_one_path_label_per_coincident_group():
    points = {
        "A": np.array([0.0, 0.0, 0.0]),
        "H": np.array([0.0, 0.0, 0.0]),
        "H_2": np.array([0.0, 0.0, 0.0]),
        "E": np.array([1.0, 0.0, 0.0]),
        "M": np.array([1.0, 0.0, 0.0]),
        "M_2": np.array([1.0, 0.0, 0.0]),
    }
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    try:
        plot_ibz(
            ax, points, [("A", "E")], {}, None,
            np.array([0.5, 0.0, 0.0]),
        )
        assert [text.get_text() for text in ax.texts] == [
            _math_label("A"), _math_label("E")
        ]
        assert len(ax.collections) == 2
    finally:
        plt.close(fig)


def test_combined_label_primes_each_name_independently():
    assert prime_point_label("A/H_2") == "A'/H_2'"
    assert _math_label("A'/H_2'") == (
        r"$\mathbf{A^{\prime}/H_{2}^{\prime}}$"
    )
