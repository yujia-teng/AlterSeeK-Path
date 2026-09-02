"""Regression tests for the curated HPKOT tables in lattice_kpoints.py.

Golden values below were captured on 2026-06-12 from the state that passed
the full 54-case manual validation (MSG and SeeK-path path audits 54/54).
If a test fails after a code change, the change altered validated behavior:
either fix the code or consciously update the golden value here.
"""

import numpy as np
import pytest
from scipy.spatial import ConvexHull

from alterseek import lattice_kpoints as lk
from alterseek.compute_centroid_3d import _prepare_enlarged_ibz_path
from alterseek.geometry import calculate_volume_centroid

GAMMA = "Γ"


# ---------------------------------------------------------------------------
# Structural invariants: every HPKOT type must stay self-consistent.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lattice_type", lk.HPKOT_LATTICE_TYPES)
def test_kpath_labels_exist_in_point_table(lattice_type):
    points = lk.get_kpoints(lattice_type)
    path = lk.get_kpath(lattice_type)
    assert path, f"{lattice_type}: empty k-path"
    assert all(len(segment) == 2 for segment in path)
    missing = [label for segment in path for label in segment
               if label not in points]
    assert not missing, f"{lattice_type}: path labels missing from points: {missing}"


@pytest.mark.parametrize("lattice_type", lk.HPKOT_LATTICE_TYPES)
def test_gamma_present_and_at_origin(lattice_type):
    points = lk.get_kpoints(lattice_type)
    assert GAMMA in points
    assert list(points[GAMMA]) == [0.0, 0.0, 0.0]


@pytest.mark.parametrize("lattice_type", lk.HPKOT_LATTICE_TYPES)
def test_hull_kpoints_nonempty(lattice_type):
    hull = lk.get_hull_kpoints(lattice_type)
    assert hull
    assert GAMMA in hull


# ---------------------------------------------------------------------------
# Convention goldens (non-negotiable project rules).
# ---------------------------------------------------------------------------

def test_cubic_caption_segment_only_for_listed_spacegroups():
    # Hinuma et al. Table 69 caption: M-X_1 only for SG {195,198,200,201,205}.
    assert lk.get_kpath("cP1", 195)[-1] == ("M", "X_1")
    assert ("M", "X_1") not in lk.get_kpath("cP1", 221)


def test_hp1_trigonal_keeps_seekpath_K_H2_segment():
    # hP1 trigonal middle-row convention (2026-06-11 decision).
    assert lk.get_kpath("hP1", 149)[-1] == ("K", "H_2")
    assert ("K", "H_2") not in lk.get_kpath("hP1", 191)


def test_hp1_trigonal_hull_doubles_in_plane_with_KA_HA():
    # SG 149 family: the second half is the in-plane image across Gamma-M, so
    # the copies are K_A/H_A and H_2 at kz = -1/2 is outside the hull.
    keys_149 = sorted(lk.get_hull_kpoints("hP1", spacegroup_number=149))
    assert keys_149 == sorted(["A", "H", "H_A", "K", "K_A", "L", "M", GAMMA])
    keys_191 = sorted(lk.get_hull_kpoints("hP1", spacegroup_number=191))
    assert keys_191 == sorted(["A", "H", "K", "L", "M", GAMMA])


def test_hp2_trigonal_and_hexagonal_share_MA_LA():
    # hP2 uses the same in-plane copied sector across Gamma-K for trigonal
    # -3m (e.g. SG 165) and hexagonal 6/m (e.g. SG 170).
    expected = sorted(["A", "H", "K", "L", "L_A", "M", "M_A", GAMMA])
    assert sorted(lk.get_hull_kpoints("hP2", spacegroup_number=170)) == expected
    assert sorted(lk.get_hull_kpoints("hP2", spacegroup_number=165)) == expected


def test_tp1_4m_hull_uses_RA_XA():
    # tP1 4/m doubled-IBZ adds R_A and X_A (path addition A-R_A-X_A-M).
    keys = sorted(lk.get_hull_kpoints("tP1", spacegroup_number=80))
    assert keys == sorted(["A", "M", "R", "R_A", "X", "X_A", "Z", GAMMA])


def test_ti1_4m_hull_uses_XA_PA():
    # tI1 4/m doubled-IBZ keeps every copy a corner of the doubled wedge:
    # X_A/P_A (path addition N-P_A-X_A-M), not the M_A/N_A/Z_0A half.
    keys = sorted(lk.get_hull_kpoints("tI1", 1.0, c=2.0, spacegroup_number=87))
    assert keys == sorted(
        ["M", "N", "P", "P_A", "X", "X_A", "Z", "Z_0", GAMMA]
    )


def test_ti2_literal_G_is_not_gamma():
    # In tI2 the literal label "G" is a distinct high-symmetry point.
    points = lk.get_kpoints("tI2")
    assert "G" in points
    assert GAMMA in points
    assert list(points["G"]) != list(points[GAMMA])


@pytest.mark.parametrize(
    "lattice_type,spacegroup,kwargs,factor,centroid,extra_segments,stubs",
    [
        ("cP1", 200, {}, 2, [0.3125, 0.3125, 0.125],
         [("M", "X_A")], []),
        ("cF1", 202, {}, 2, [0.431640625, 0.23046875, 0.431640625],
         [("X", "W_A")], ["K_A"]),
        ("cI1", 204, {}, 2, [0.25, -0.0625, 0.25],
         [], ["N_A"]),
        ("hP1", 147, {"c": 1.2}, 4, [7 / 24, 0.0, 0.25],
         [("H_A", "K_A")], ["M_A", "L_A", "M_B", "L_B"]),
        ("hR1", 148, {"a": 5.141007422395145,
                       "c": 14.67537993611477}, 2,
         [0.32317619, 0.01285480, 0.32317619],
         [], ["L_A", "S_0A", "S_2A", "H_0A", "H_2A", "M_4A"]),
        ("hR2", 148, {"a": 10.847751571259092,
                       "c": 9.47659577265194}, 2,
         [0.375, -0.04093461, -0.04093461],
         [], ["T_A", "P_0A", "F_A"]),
    ],
)
def test_fixed_minus3_mbar3_enlargements(
    lattice_type, spacegroup, kwargs, factor, centroid, extra_segments, stubs,
):
    base = lk.get_hull_kpoints(lattice_type, **kwargs)
    enlarged = lk.get_hull_kpoints(
        lattice_type, spacegroup_number=spacegroup, **kwargs)
    base_hull = ConvexHull(np.array(list(base.values())))
    enlarged_hull = ConvexHull(np.array(list(enlarged.values())))
    enlarged_centroid, _ = calculate_volume_centroid(enlarged_hull)

    assert enlarged_hull.volume / base_hull.volume == pytest.approx(factor)
    assert enlarged_centroid == pytest.approx(centroid, abs=1e-6)

    ordinary_path = lk.get_hull_kpath(lattice_type)
    enlarged_path = lk.get_hull_kpath(lattice_type, spacegroup)
    assert enlarged_path[:len(ordinary_path)] == ordinary_path
    assert enlarged_path[len(ordinary_path):] == extra_segments

    band_path, band_points, general_vertices = _prepare_enlarged_ibz_path(
        lattice_type, spacegroup, ordinary_path, enlarged_path,
        base, enlarged,
    )
    assert band_path == enlarged_path
    assert general_vertices == stubs
    assert set(lk.get_project_hull_extra_labels(
        lattice_type, spacegroup)).issubset(band_points)
