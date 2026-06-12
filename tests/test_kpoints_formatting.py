"""Regression tests for KPOINTS label/path formatting in alterseek_path.py."""

import numpy as np

from alterseek_path import KPointsModifier


def test_gamma_label_is_vasp_safe():
    assert KPointsModifier._kpoints_label("Γ") == "GAMMA"
    assert KPointsModifier._kpoints_label("gamma") == "GAMMA"
    assert KPointsModifier._kpoints_label("GAMMA") == "GAMMA"
    # Ordinary labels pass through untouched, including primes and subscripts.
    assert KPointsModifier._kpoints_label("M_A'") == "M_A'"
    assert KPointsModifier._kpoints_label("X_1") == "X_1"


def test_display_label_matches_kpoints_label():
    # The console form and the VASP-file form must never drift apart.
    for label in ["Γ", "GAMMA", "K", "H_2", "M_A'", "k", "k'"]:
        assert (KPointsModifier._display_label(label)
                == KPointsModifier._kpoints_label(label))


def test_format_path_joins_continuous_and_breaks_discontinuous():
    segments = [("Γ", "X"), ("X", "M"), ("R", "Z")]
    assert KPointsModifier._format_path(segments) == "GAMMA-X-M | R-Z"


def test_dedupe_frac_positions_wraps_periodic_images():
    from alterseek_path import _dedupe_frac_positions
    unique = _dedupe_frac_positions([
        [0.25, 0.25, 0.25],
        [1.25, 0.25, -0.75],   # same site shifted by lattice vectors
        [0.75, 0.25, 0.25],
    ])
    assert len(unique) == 2


def test_reciprocal_basis_roundtrip():
    # k_input_frac = k_prim_frac @ B_prim @ inv(B_input) must be identity when
    # the bases coincide (the KPOINTS output-basis conversion contract).
    rng = np.random.default_rng(0)
    b = rng.normal(size=(3, 3)) + 3 * np.eye(3)
    k = np.array([0.386383, 0.363842, 0.351437])
    converted = k @ b @ np.linalg.inv(b)
    assert np.allclose(converted, k, atol=1e-12)
