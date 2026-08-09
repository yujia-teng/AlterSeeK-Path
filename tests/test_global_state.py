import json
import subprocess
import sys
import warnings
from pathlib import Path

import matplotlib as mpl


REPO_ROOT = Path(__file__).resolve().parents[1]
STYLE_KEYS = (
    "font.family",
    "font.serif",
    "mathtext.fontset",
    "xtick.direction",
    "ytick.direction",
)
CALLER_STYLE = {
    "font.family": "sans-serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavusans",
    "xtick.direction": "out",
    "ytick.direction": "out",
}


def _style_snapshot():
    return {key: mpl.rcParams[key] for key in STYLE_KEYS}


def test_imports_leave_caller_warning_filters_and_plot_style_unchanged():
    code = r'''
import json
import warnings

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot
import numpy
import scipy
import seekpath
import spglib
import sympy
from ase.io import read
from findspingroup import find_spin_group_basic
from pymatgen.core import Structure

style = {
    "font.family": "sans-serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavusans",
    "xtick.direction": "out",
    "ytick.direction": "out",
}
mpl.rcParams.update(style)
warnings.filterwarnings("error", message="caller-owned sentinel")
before_filters = list(warnings.filters)
before_style = {key: mpl.rcParams[key] for key in style}

import alterseek.compute_centroid_hybrid
import alterseek.find_sf_operations
import plot_alterband
import plot_alterband_qe

after_style = {key: mpl.rcParams[key] for key in style}
print(json.dumps({
    "filters_unchanged": warnings.filters == before_filters,
    "style_unchanged": after_style == before_style,
}))
'''
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    state = json.loads(result.stdout.strip().splitlines()[-1])
    assert state == {
        "filters_unchanged": True,
        "style_unchanged": True,
    }


def test_core_figure_stage_uses_scoped_style_on_repeated_runs(
    tmp_path, monkeypatch
):
    from alterseek import compute_centroid_hybrid as centroid_module

    observed_styles = []

    def fake_setup(*args, **kwargs):
        observed_styles.append(mpl.rcParams["mathtext.fontset"])
        return object(), object()

    monkeypatch.setattr(centroid_module, "get_bz_loops", lambda matrix: [
        centroid_module.np.array([
            [-0.5, -0.5, 0.0],
            [0.5, -0.5, 0.0],
            [0.5, 0.5, 0.0],
            [-0.5, 0.5, 0.0],
            [-0.5, -0.5, 0.0],
        ])
    ])
    monkeypatch.setattr(centroid_module, "setup_3d_ax", fake_setup)
    monkeypatch.setattr(centroid_module, "plot_ibz", lambda *args, **kwargs: None)
    monkeypatch.setattr(centroid_module.plt, "tight_layout", lambda: None)
    monkeypatch.setattr(centroid_module.plt, "close", lambda figure: None)
    monkeypatch.setattr(centroid_module, "_save_figure", lambda *args, **kwargs: [])
    monkeypatch.setattr(centroid_module, "_print_saved_paths", lambda *args, **kwargs: None)

    analysis = {
        "sg": 136,
        "sc_type": "tP1",
        "sc_display": "tP1",
        "b_matrix": centroid_module.np.eye(3),
        "kpath_plot": [],
        "display_labels_plot": {},
        "kpoints_cart_plot": {},
    }
    centroid = {
        "hull": None,
        "centroid_cart": centroid_module.np.zeros(3),
        "points_arr": centroid_module.np.empty((0, 3)),
        "labels_list": [],
    }

    with mpl.rc_context(CALLER_STYLE):
        expected = _style_snapshot()
        for _ in range(2):
            centroid_module._generate_figure1(
                analysis,
                centroid,
                output_dir=tmp_path,
                fig_basename="scope-test",
                show_plot=False,
                defer_show=False,
                mode_2d=False,
                view_elev=None,
                view_azim=None,
                save_pdf=False,
                verbose=False,
            )
            assert _style_snapshot() == expected

    assert observed_styles == ["stix", "stix"]


def test_kspace_analysis_restores_warning_filters_on_repeated_runs():
    from alterseek import compute_centroid_hybrid as centroid_module

    structure = REPO_ROOT / "tests" / "references" / "case12_POSCAR"
    before = list(warnings.filters)
    for _ in range(2):
        centroid_module._analyze_kspace(
            structure,
            seekpath_type_numbers=None,
            mode_2d=False,
            input_vacuum_axis=2,
            symprec=None,
            verbose=False,
        )
        assert warnings.filters == before


def test_seekpath_lattice_tag_suppresses_only_the_third_party_deprecation():
    """The second seekpath call site was left unscoped by the Point 9 pass.

    Removing the process-global filters surfaced SeeK-path 2.1's use of
    spglib's deprecated dict interface, which this project cannot fix
    upstream, so it is suppressed locally like the other call site.
    """
    from alterseek.find_sf_operations import _seekpath_lattice_tag

    lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]
    positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

    before = list(warnings.filters)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tag = _seekpath_lattice_tag(lattice, positions, [1, 1], 1e-3)

    assert tag == "cI1"
    assert not [
        entry for entry in caught
        if "dict interface is deprecated" in str(entry.message)
    ]
    # Suppression must be scoped: the caller's own filters survive the call.
    assert warnings.filters == before


def test_seekpath_lattice_tag_still_reports_unrelated_warnings():
    """Scoping must not turn into a blanket mute of the whole call."""
    from alterseek.find_sf_operations import _seekpath_lattice_tag

    lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warnings.warn("caller-owned sentinel", UserWarning)
        _seekpath_lattice_tag(lattice, [[0.0, 0.0, 0.0]], [1], 1e-3)

    assert [
        entry for entry in caught
        if "caller-owned sentinel" in str(entry.message)
    ]


def test_band_plotters_restore_caller_style_on_repeated_runs(tmp_path):
    from plot_alterband import plot_alterband
    from plot_alterband_qe import plot_alterband_qe

    klabels = tmp_path / "KLABELS"
    klabels.write_text("A 0.0\nB 1.0\n", encoding="utf-8")
    vasp_bands = "k band\n0.0 0.0\n0.5 0.1\n1.0 0.0\n"
    vasp_up = tmp_path / "vasp_up.dat"
    vasp_down = tmp_path / "vasp_down.dat"
    vasp_up.write_text(vasp_bands, encoding="utf-8")
    vasp_down.write_text(vasp_bands, encoding="utf-8")

    qe_kpoints = tmp_path / "KPOINTS_alter_qe"
    qe_kpoints.write_text(
        "K_POINTS crystal_b\n"
        "2\n"
        "0.0 0.0 0.0 2 ! A\n"
        "0.5 0.0 0.0 1 ! B\n",
        encoding="utf-8",
    )
    qe_bands = "0.0 0.0\n0.5 0.1\n1.0 0.0\n"
    qe_up = tmp_path / "qe_up.gnu"
    qe_down = tmp_path / "qe_down.gnu"
    qe_up.write_text(qe_bands, encoding="utf-8")
    qe_down.write_text(qe_bands, encoding="utf-8")

    with mpl.rc_context(CALLER_STYLE):
        expected = _style_snapshot()
        for index in range(2):
            plot_alterband(
                klabels=klabels,
                band_up=vasp_up,
                band_down=vasp_down,
                output=tmp_path / f"vasp-{index}.png",
            )
            assert _style_snapshot() == expected
            plot_alterband_qe(
                band_up=qe_up,
                band_down=qe_down,
                kpoints_qe=qe_kpoints,
                output=tmp_path / f"qe-{index}.png",
            )
            assert _style_snapshot() == expected
