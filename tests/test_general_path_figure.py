import numpy as np
from scipy.spatial import ConvexHull

from alterseek.geometry import get_bz_loops
from alterseek.kpoints import KPathBuilder
from alterseek.lattice_kpoints import get_hull_kpath, get_hull_kpoints
from alterseek.mode2d.plotting import plot_2d_general_path_figure
from alterseek.plotting_3d import plot_general_path_figure


def test_red_only_general_path_figure_renders_written_sequence(tmp_path):
    points = get_hull_kpoints("cP1", spacegroup_number=200)
    path = get_hull_kpath("cP1", spacegroup_number=200)
    b_matrix = np.eye(3)
    hull_points = np.array(list(points.values()), dtype=float)
    hull = ConvexHull(hull_points)

    builder = KPathBuilder()
    builder.kpoints_data = [
        [*points[label], label]
        for segment in path
        for label in segment
    ]
    sequence = builder.build_ordinary_path_with_general_k(
        [0.3125, 0.3125, 0.125]
    )
    bz_loops = get_bz_loops(b_matrix)
    output = tmp_path / "POSCAR_generalpath_cP1.png"

    plot_general_path_figure(
        b_matrix=b_matrix,
        bz_loops=bz_loops,
        bz_center=np.mean(np.vstack(bz_loops), axis=0),
        kpoints_data=builder.kpoints_data,
        ibz_kpoints_frac=points,
        hull_pts=hull_points,
        hull_simplices=hull.simplices,
        centroid_frac=[0.3125, 0.3125, 0.125],
        output_path=str(output),
        path_sequence=sequence,
        show_plot=False,
    )

    assert output.is_file()
    assert output.stat().st_size > 10_000


def test_non_altermagnetic_figure_helper_picks_the_matching_plate(
    tmp_path, monkeypatch
):
    from alterseek import kpoints as kpoints_module

    calls = []
    sentinel = object()

    def fake_plot(**kwargs):
        calls.append(kwargs)
        return sentinel

    monkeypatch.setattr(kpoints_module, "OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(kpoints_module, "plot_general_path_figure", fake_plot)
    centroid = {
        "bz_loops": [np.zeros((2, 3))],
        "bz_center": np.zeros(3),
        "b_matrix": np.eye(3),
        "ibz_kpoints_frac": {"Γ": [0.0, 0.0, 0.0]},
        "hull_pts": np.zeros((4, 3)),
        "hull_simplices": [[0, 1, 2]],
        "sc_type": "cP1",
        "elev": 14,
        "azim": 20,
    }
    builder = KPathBuilder()
    builder.kpoints_data = []
    figures = []
    sequence = [[0.0, 0.0, 0.0, "Γ"], [0.1, 0.2, 0.3, "k"]]

    builder.mode_2d = False
    builder._generate_non_altermagnetic_figure(
        centroid, "POSCAR", [0.1, 0.2, 0.3], sequence, figures
    )
    assert figures == [sentinel]
    assert calls[0]["output_path"].endswith("POSCAR_generalpath_cP1.png")

    two_d_calls = []
    monkeypatch.setattr(
        kpoints_module,
        "plot_2d_general_path_figure",
        lambda *args, **kwargs: two_d_calls.append((args, kwargs)),
    )
    builder.mode_2d = True
    builder._generate_non_altermagnetic_figure(
        centroid, "POSCAR", [0.1, 0.2, 0.3], sequence, figures
    )
    # 2D goes to the slab plate, never to the 3D one.
    assert len(calls) == 1
    assert len(two_d_calls) == 1
    assert two_d_calls[0][0][2] == "POSCAR"


def test_2d_general_path_figure_renders_a_slab_plate(tmp_path):
    b_matrix = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.25],
    ])
    points = {
        "GAMMA": np.array([0.0, 0.0, 0.0]),
        "X": np.array([0.5, 0.0, 0.0]),
        "M": np.array([0.5, 0.5, 0.0]),
    }
    centroid = {
        "b_matrix": b_matrix,
        "vacuum_axis": 2,
        "sc_type": "tP1",
        "band_kpath": [("GAMMA", "X"), ("X", "M"), ("M", "GAMMA")],
        "band_kpoints_frac": points,
        "ibz_polygon_frac": list(points.values()),
        "ibz_polygon_labels": list(points),
    }
    builder = KPathBuilder()
    builder.kpoints_data = [
        [*points[label], label]
        for segment in centroid["band_kpath"]
        for label in segment
    ]
    sequence = builder.build_ordinary_path_with_general_k([0.25, 0.15, 0.0])

    saved = plot_2d_general_path_figure(
        centroid,
        np.array([0.25, 0.15, 0.0]),
        "POSCAR",
        output_dir=str(tmp_path),
        path_sequence=sequence,
    )

    assert len(saved) == 1
    output = tmp_path / "POSCAR_2d_generalpath_tP1.png"
    assert output.is_file()
    assert output.stat().st_size > 10_000
