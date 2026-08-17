"""Regression tests for the curated HPKOT tables in lattice_kpoints.py.

Golden values below were captured on 2026-06-12 from the state that passed
the full 54-case manual validation (MSG and SeeK-path path audits 54/54).
If a test fails after a code change, the change altered validated behavior:
either fix the code or consciously update the golden value here.
"""

import pytest

from alterseek import lattice_kpoints as lk

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


def test_aliases_resolve_to_hpkot_keys():
    assert lk.canonical_lattice_type("HEX") == "hP2"
    assert lk.canonical_lattice_type("BCT2") == "tI2"
    assert lk.canonical_lattice_type("ORCI") == "oI1"
    # Canonical keys resolve to themselves.
    for key in lk.HPKOT_LATTICE_TYPES:
        assert lk.canonical_lattice_type(key) == key


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
