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

import compute_centroid_hybrid as cc
import geometry
import symmetry


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
