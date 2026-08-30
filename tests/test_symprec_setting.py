"""The symmetry tolerance is a physics setting, not a file-format detail.

The tolerance used to be applied only to .mcif input, so the *same structure
with the same moments* reported different symmetry as .mcif and as POSCAR:
MnSe2 as .mcif correctly reported its cubic Pa-3 parent, while as POSCAR it
read as orthorhombic.

The configured tolerance now applies to every format. Only the *extra*
validated recovery stays mcif-specific, since only an mcif declares a parent
space group to validate a loosened tolerance against.
"""
import numpy as np
import pytest

from alterseek import find_sf_operations as F
from alterseek import compute_centroid_3d as C
from alterseek.kpoints import _validate_input_config


def test_default_tolerance_is_1e_3_in_both_modules():
    """Deposited coordinates are routinely rounded to 5 decimals; spglib's own
    1e-5 default hides symmetry that is really present."""
    assert F._DEFAULT_SYMPREC == 1e-3
    assert C._DEFAULT_SYMPREC == 1e-3


def _cubic_cell_with_rounded_coordinates():
    """A cell whose symmetry is only visible above ~1e-4."""
    lattice = np.eye(3) * 4.0
    # Two sites that are equivalent under a centering translation, but with
    # the second nudged by 5e-4 fractional (2e-3 A) -- invisible at 1e-5.
    positions = np.array([[0.0, 0.0, 0.0], [0.5005, 0.5, 0.5]])
    numbers = np.array([26, 26])
    return lattice, positions, numbers


def test_tolerance_is_not_gated_on_file_extension():
    """Same coordinates, same tolerance -> same answer for either format."""
    lattice, positions, numbers = _cubic_cell_with_rounded_coordinates()
    # "x.mcif" declares no parent, so no validated recovery can apply and both
    # calls must fall through to the configured tolerance.
    as_mcif = F._non_magnetic_symmetry("x.mcif", lattice, positions, numbers,
                                       True, symprec=1e-2)
    as_poscar = F._non_magnetic_symmetry("POSCAR", lattice, positions, numbers,
                                         False, symprec=1e-2)
    assert as_mcif == as_poscar


def test_configured_tolerance_actually_reaches_spglib():
    """A tight tolerance and a loose one must disagree on this cell, or the
    test above would pass vacuously."""
    lattice, positions, numbers = _cubic_cell_with_rounded_coordinates()
    tight = F._non_magnetic_symmetry("POSCAR", lattice, positions, numbers,
                                     False, symprec=1e-6)
    loose = F._non_magnetic_symmetry("POSCAR", lattice, positions, numbers,
                                     False, symprec=1e-2)
    assert tight['spacegroup_number'] != loose['spacegroup_number']


def test_validated_mcif_recovery_still_wins_over_the_configured_tolerance():
    """When a declared parent exists it is checked against ground truth, so it
    is preferred over a bare tolerance."""
    lattice, positions, numbers = _cubic_cell_with_rounded_coordinates()
    called = {}

    def fake_hint(filename):
        called['hint'] = filename
        return None  # no declaration -> must fall back to the configured value

    original = F._declared_mcif_parent_hint
    try:
        F._declared_mcif_parent_hint = fake_hint
        chosen = F._select_mcif_symprec_for_non_magnetic_label(
            "x.mcif", lattice, positions, numbers, fallback=0.25)
    finally:
        F._declared_mcif_parent_hint = original
    assert called['hint'] == "x.mcif"
    assert chosen == 0.25


@pytest.mark.parametrize("value", [1e-5, 1e-3, 0.01, 1, 2.5])
def test_symprec_accepts_positive_numbers(value):
    assert _validate_input_config({"symprec": value}) == {"symprec": value}


@pytest.mark.parametrize("value", [0, -1, -1e-9])
def test_symprec_rejects_non_positive(value):
    with pytest.raises(ValueError, match="symprec must be positive"):
        _validate_input_config({"symprec": value})


@pytest.mark.parametrize("value", [np.inf, -np.inf, np.nan])
def test_symprec_rejects_non_finite_numbers(value):
    with pytest.raises(ValueError, match="symprec must be finite"):
        _validate_input_config({"symprec": value})


@pytest.mark.parametrize("value", ["1e-3", True, False, None, [1e-3]])
def test_symprec_rejects_non_numbers(value):
    with pytest.raises(ValueError, match="symprec must be a number"):
        _validate_input_config({"symprec": value})
