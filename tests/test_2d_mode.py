"""Tests for 2D / slab mode (compute_centroid_hybrid + the op classification).

Covers the pieces that have no 3D analogue:
  * vacuum-axis detection in the standardized frame (incl. axis permutation),
  * the 2D area centroid and in-plane IBZ polygon,
  * the in-plane spin-flip operation filter (Filter 1 + Filter 2) that drives
    the 2D altermagnetism YES/NO verdict,
  * the input-slab sanity check,
  * an end-to-end compute_centroid(mode_2d=True) run on a synthetic tetragonal
    slab, asserting the path/centroid are restricted to the vacuum k=0 plane.
"""

import numpy as np
import pytest

from alterseek import compute_centroid_hybrid as cc
from alterseek import geometry
from alterseek import symmetry
from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
from alterseek.plotting_common import _figure_output_paths


def _diag(vals):
    return np.diag(np.asarray(vals, dtype=float))


# Reference operations with the vacuum on axis 2 (z); in-plane block = rows/cols (0,1).
C2Z = _diag([-1, -1, 1])                       # in-plane -I  -> trivial
MZ = _diag([1, 1, -1])                          # in-plane +I  -> trivial
INV = _diag([-1, -1, -1])                       # in-plane -I  -> trivial
C4Z = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])  # in-plane rotation -> valid
MX = _diag([-1, 1, 1])                          # in-plane diag(-1, 1) -> valid
C2X = _diag([1, -1, -1])                        # in-plane diag(1, -1) -> valid


# ---------------------------------------------------------------------------
# Filter 2: trivial 2D spin-flip (in-plane block +-I)  -- the NO-verdict logic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("op", [C2Z, MZ, INV])
def test_trivial_flip_ops_axis_z(op):
    assert symmetry.is_trivial_2d_spin_flip(op, 2)
    assert not symmetry.is_valid_2d_spin_flip(op, 2)


@pytest.mark.parametrize("op", [C4Z, MX, C2X])
def test_valid_flip_ops_axis_z(op):
    assert not symmetry.is_trivial_2d_spin_flip(op, 2)
    assert symmetry.is_valid_2d_spin_flip(op, 2)


def test_classification_is_axis_aware_permutation():
    # Same physics with the vacuum on axis 0 (x); in-plane block = rows/cols (1,2).
    c4x = np.array([[1., 0, 0], [0, 0, -1], [0, 1, 0]])   # rotation about x -> valid
    c2x = _diag([1, -1, -1])                               # in-plane (y,z) = -I -> trivial
    mx = _diag([-1, 1, 1])                                 # in-plane (y,z) = +I -> trivial
    assert symmetry.is_valid_2d_spin_flip(c4x, 0)
    assert symmetry.is_trivial_2d_spin_flip(c2x, 0)
    assert symmetry.is_trivial_2d_spin_flip(mx, 0)


# ---------------------------------------------------------------------------
# Filter 1: operation must keep the plane
# ---------------------------------------------------------------------------

def test_keeps_plane_guard():
    assert symmetry.keeps_2d_plane(C4Z, 2)
    out_of_plane = np.array([[0., 0, 1], [0, 1, 0], [1, 0, 0]])  # swaps x and z
    assert not symmetry.keeps_2d_plane(out_of_plane, 2)
    assert not symmetry.is_valid_2d_spin_flip(out_of_plane, 2)         # fails Filter 1


def test_cartesian_plane_filter_survives_magnetic_axis_permutation():
    """Operation basis (c,a,b) must not reuse submitted vacuum index c=2."""
    submitted_lattice = _diag([3.0, 3.0, 20.0])
    magnetic_lattice = np.array([
        [0.0, 0.0, 20.0],  # a_magnetic = c_input: physical vacuum
        [3.0, 0.0, 0.0],   # b_magnetic = a_input
        [0.0, 3.0, 0.0],   # c_magnetic = b_input
    ])
    b_magnetic = 2 * np.pi * np.linalg.inv(magnetic_lattice).T
    plane_normal = symmetry.slab_plane_normal_cartesian(
        submitted_lattice, vacuum_axis=2
    )

    # Leaves a_magnetic (the physical vacuum) fixed and swaps the two physical
    # in-plane axes. The old index-2 check rejects it because it mistakes
    # c_magnetic for the vacuum.
    valid_physical_flip = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ])
    assert not symmetry.is_valid_2d_spin_flip(valid_physical_flip, 2)
    assert symmetry.is_valid_2d_spin_flip_cartesian(
        valid_physical_flip, b_magnetic, plane_normal
    )

    # Swaps physical vacuum and in-plane axes. The old index-2 check accepts it
    # as an in-plane swap, while the Cartesian test correctly rejects it.
    mixes_physical_vacuum = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    assert symmetry.is_valid_2d_spin_flip(mixes_physical_vacuum, 2)
    assert not symmetry.is_valid_2d_spin_flip_cartesian(
        mixes_physical_vacuum, b_magnetic, plane_normal
    )

    modifier = KPointsModifier(mode_2d=True, input_vacuum_axis=2)
    centroid_result = {"b_matrix_input": b_magnetic}
    modifier._configure_2d_plane(
        centroid_result,
        submitted_lattice=submitted_lattice,
    )
    assert modifier._is_valid_2d_operation(
        valid_physical_flip, centroid_result
    )
    assert not modifier._is_valid_2d_operation(
        mixes_physical_vacuum, centroid_result
    )


def test_2d_output_projection_uses_physical_plane_not_input_axis_index():
    magnetic_lattice = np.array([
        [0.0, 0.0, 20.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
    ])
    b_magnetic = 2 * np.pi * np.linalg.inv(magnetic_lattice).T
    modifier = KPointsModifier(mode_2d=True, input_vacuum_axis=2)
    modifier.kpoints_basis_matrix = b_magnetic
    modifier.output_basis_matrix = b_magnetic
    modifier.plane_normal_cartesian = np.array([0.0, 0.0, 1.0])

    output = modifier._kpoint_for_output_basis(
        [1e-10, 0.2, 0.3, "k"]
    )
    assert output[:3] == pytest.approx([0.0, 0.2, 0.3], abs=1e-12)


def test_2d_output_rejects_a_genuinely_out_of_plane_point():
    lattice = _diag([3.0, 3.0, 20.0])
    reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
    modifier = KPointsModifier(mode_2d=True, input_vacuum_axis=2)
    modifier.kpoints_basis_matrix = reciprocal
    modifier.output_basis_matrix = reciprocal
    modifier.plane_normal_cartesian = np.array([0.0, 0.0, 1.0])

    with pytest.raises(RuntimeError, match="outside the physical slab plane"):
        modifier._kpoint_for_output_basis([0.2, 0.3, 0.1, "bad"])


# ---------------------------------------------------------------------------
# Vacuum-axis detection
# ---------------------------------------------------------------------------

def test_detect_vacuum_axis_shortest_reciprocal():
    # Reciprocal rows: axis 1 is much shorter -> the vacuum direction.
    b = _diag([2.0, 0.3, 2.0])
    axis, info = geometry.detect_vacuum_axis_2d(b)
    assert axis == 1
    assert info["separated"] and info["orthogonal"]


def test_detect_vacuum_axis_ambiguous_when_cubic():
    axis, info = geometry.detect_vacuum_axis_2d(_diag([1.0, 1.0, 1.0]))
    assert not info["separated"]


# ---------------------------------------------------------------------------
# 2D area centroid + in-plane IBZ polygon
# ---------------------------------------------------------------------------

def test_area_centroid_unit_square():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    cfrac, ccart, area = geometry.area_centroid_2d(pts, 2, np.eye(3))
    assert np.allclose(cfrac, [0.5, 0.5, 0.0])
    assert cfrac[2] == 0.0
    assert abs(area - 1.0) < 1e-9


def test_area_centroid_triangle():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    cfrac, _ccart, area = geometry.area_centroid_2d(pts, 2, np.eye(3))
    assert np.allclose(cfrac, [1 / 3, 1 / 3, 0.0])
    assert abs(area - 0.5) < 1e-9


def test_ordered_polygon_drops_interior_point():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0.5, 0.5, 0]], dtype=float)
    poly = geometry.ordered_2d_polygon_frac(pts, 2)
    assert len(poly) == 4                                  # interior point excluded
    assert all(abs(p[2]) < 1e-12 for p in poly)


# ---------------------------------------------------------------------------
# Input-slab sanity check
# ---------------------------------------------------------------------------

def test_input_slab_clean():
    assert geometry.check_input_slab(_diag([3.0, 3.0, 20.0]), 2) == []


def test_input_slab_tilted_axis_warns():
    tilted = np.array([[3., 0, 0], [0, 3, 0], [1.5, 0, 20.]])  # c not orthogonal to a
    assert geometry.check_input_slab(tilted, 2)                       # non-empty


def test_input_slab_no_vacuum_warns():
    assert geometry.check_input_slab(_diag([3.0, 3.0, 2.0]), 2)       # c not the longest


# ---------------------------------------------------------------------------
# End-to-end: compute_centroid(mode_2d=True) on a synthetic tetragonal slab
# ---------------------------------------------------------------------------

def _write_tetragonal_slab(path):
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.vasp import Poscar
    lattice = Lattice.from_parameters(3.0, 3.0, 20.0, 90, 90, 90)
    structure = Structure(lattice, ["Po"], [[0.0, 0.0, 0.0]])
    Poscar(structure).write_file(str(path))


def test_compute_centroid_2d_tetragonal_slab(tmp_path):
    poscar = tmp_path / "POSCAR"
    _write_tetragonal_slab(poscar)
    result = cc.run(str(poscar), output_dir=str(tmp_path), show_plot=False,
                    verbose=False, mode_2d=True)

    assert result["vacuum_axis"] == 2
    assert result["mode_2d"] is True
    # Centroid is in the physical plane.
    assert abs(result["centroid_frac"][2]) < 1e-12
    # Every band-path label lies in the vacuum k=0 plane.
    for segment in result["band_kpath"]:
        for label in segment:
            frac = result["band_kpoints_frac"].get(label)
            if frac is not None:
                assert abs(frac[2]) < 1e-6, f"{label} is out of plane"
    # The IBZ polygon is a real 2D polygon at z=0.
    poly = result["ibz_polygon_frac"]
    assert poly is not None and len(poly) >= 3
    assert all(abs(vertex[2]) < 1e-12 for vertex in poly)


def test_compute_centroid_2d_differs_from_3d(tmp_path):
    poscar = tmp_path / "POSCAR"
    _write_tetragonal_slab(poscar)
    r2 = cc.run(str(poscar), output_dir=str(tmp_path), show_plot=False,
                verbose=False, mode_2d=True)
    r3 = cc.run(str(poscar), output_dir=str(tmp_path), show_plot=False,
                verbose=False, mode_2d=False)
    # 3D keeps the out-of-plane centroid component; 2D zeroes it.
    assert abs(r3["centroid_frac"][2]) > 1e-6
    assert abs(r2["centroid_frac"][2]) < 1e-12


@pytest.mark.parametrize("stale_spin_log", [None, "old structure sentinel\n"])
def test_interactive_2d_output_stays_in_physical_plane(
    tmp_path, monkeypatch, stale_spin_log
):
    poscar = tmp_path / "POSCAR"
    _write_tetragonal_slab(poscar)
    spin_log = tmp_path / OUTPUT_DIR / "spin_operations.txt"
    if stale_spin_log is not None:
        spin_log.parent.mkdir()
        spin_log.write_text(stale_spin_log, encoding="utf-8")
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = ""\n'
        'path = ""\n'
        'output_code = "vasp"\n',
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    assert KPointsModifier(mode_2d=True, input_vacuum_axis=2).interactive_modify()
    lines = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8").splitlines()
    coordinate_rows = [
        line.split()
        for line in lines[4:]
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert coordinate_rows
    assert all(float(row[2]) == 0.0 for row in coordinate_rows)
    if stale_spin_log is None:
        assert not spin_log.exists()
    else:
        assert spin_log.read_text(encoding="utf-8") == stale_spin_log


def test_save_pdf_adds_pdf_without_leaking_to_next_run(monkeypatch):
    monkeypatch.delenv("ALTERSEEK_BZ_FORMATS", raising=False)
    monkeypatch.delenv("ALTERSEEK_BZ_EXTRA_FORMATS", raising=False)

    assert _figure_output_paths(
        "first_2d_ibz_tP1.png",
        extra_formats=("pdf",),
    ) == [
        "first_2d_ibz_tP1.png",
        "first_2d_ibz_tP1.pdf",
    ]
    assert _figure_output_paths("second_2d_ibz_tP1.png") == [
        "second_2d_ibz_tP1.png",
    ]


def test_environment_extra_figure_formats_remain_supported(monkeypatch):
    monkeypatch.delenv("ALTERSEEK_BZ_FORMATS", raising=False)
    monkeypatch.setenv("ALTERSEEK_BZ_EXTRA_FORMATS", "pdf")
    assert _figure_output_paths("slab_2d_ibz_tP1.png") == [
        "slab_2d_ibz_tP1.png",
        "slab_2d_ibz_tP1.pdf",
    ]
