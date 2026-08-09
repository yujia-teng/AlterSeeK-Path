"""MCIF inputs must match AlterSeeK-Path's collinear spin model."""

import numpy as np
import pytest

from alterseek import find_sf_operations as spin_ops
from alterseek.io import _load_magnetic_input_data, _write_magnetic_mcif
from alterseek.mcif import _validate_collinear_moments


def _write_test_mcif(path, moments):
    _write_magnetic_mcif(
        str(path),
        "MCIF collinearity test",
        np.eye(3) * 4.0,
        ["Fe", "Fe", "O"],
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.5, 0.5, 0.5]),
            np.array([0.25, 0.25, 0.25]),
        ],
        np.asarray(moments, dtype=float),
    )


def test_collinearity_validator_accepts_parallel_antiparallel_and_zero_moments():
    _validate_collinear_moments([
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 3.0],
        [-2.0, -4.0, -6.0],
    ])


def test_magnetic_input_loader_accepts_rounded_collinear_mcif(tmp_path):
    mcif = tmp_path / "rounded-collinear.mcif"
    expected_moments = np.array([
        [0.34, 0.42, 0.0],
        [1.32, 1.63, 0.0],
        [0.0, 0.0, 0.0],
    ])
    _write_test_mcif(mcif, expected_moments)

    _, _, _, moments, setting = _load_magnetic_input_data(str(mcif), "", None)

    assert setting == "cartesian"
    assert np.allclose(moments, expected_moments)


def test_collinearity_validator_uses_absolute_moment_tolerance():
    _validate_collinear_moments([
        [1.0, 0.0, 0.0],
        [1.0, 0.02, 0.0],
    ])

    with pytest.raises(
        ValueError,
        match="Only collinear magnetic structures are supported",
    ):
        _validate_collinear_moments([
            [1.0, 0.0, 0.0],
            [1.0, 0.021, 0.0],
        ])


def test_magnetic_input_loader_rejects_noncollinear_mcif(tmp_path):
    mcif = tmp_path / "noncollinear.mcif"
    _write_test_mcif(
        mcif,
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
    )

    with pytest.raises(
        ValueError,
        match="Only collinear magnetic structures are supported",
    ):
        _load_magnetic_input_data(str(mcif), "", None)


def test_spin_symmetry_route_rejects_noncollinear_mcif_before_findspingroup(
    tmp_path, monkeypatch
):
    mcif = tmp_path / "noncollinear.mcif"
    _write_test_mcif(
        mcif,
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
    )
    monkeypatch.setattr(
        spin_ops,
        "_non_magnetic_symmetry",
        lambda *args, **kwargs: {
            "non_mag_label": "Pm-3m (221)",
            "point_group": "m-3m",
            "laue_group": "m-3m",
            "spacegroup_number": 221,
            "sites": 3,
            "lattice": "cP",
        },
    )

    def unexpected_findspingroup_call(*args, **kwargs):
        pytest.fail("FindSpinGroup must not run for a noncollinear MCIF")

    monkeypatch.setattr(
        spin_ops,
        "find_spin_group_basic",
        unexpected_findspingroup_call,
    )

    with pytest.raises(
        spin_ops.SpinSymmetryError,
        match="Only collinear magnetic structures are supported",
    ):
        spin_ops.run(str(mcif), "", verbose=False, output_dir=str(tmp_path))

    assert not (tmp_path / "spin_operations.txt").exists()
