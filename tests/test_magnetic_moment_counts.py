"""Manual moment counts must match sites without hiding excess values."""

import numpy as np
import pytest

from alterseek import find_sf_operations as spin_ops
from alterseek.io import _load_magnetic_input_data


def _write_two_atom_poscar(path):
    path.write_text(
        "two-site count validation\n"
        "1.0\n"
        "3.0 0.0 0.0\n"
        "0.0 3.0 0.0\n"
        "0.0 0.0 3.0\n"
        "Fe\n"
        "2\n"
        "Direct\n"
        "0.0 0.0 0.0\n"
        "0.5 0.5 0.5\n",
        encoding="utf-8",
    )


def test_missing_manual_moments_still_fill_with_zero():
    values = spin_ops.fit_magmoms_to_structure(
        spin_ops.parse_magmoms("1 -1"),
        4,
    )
    assert values == [1.0, -1.0, 0.0, 0.0]


@pytest.mark.parametrize(
    ("moments", "provided", "atoms", "excess"),
    [
        ("1 -1 0 0 7", 5, 4, 1),
        ("2*1 3*-1", 5, 4, 1),
        ("6*0", 6, 4, 2),
    ],
)
def test_excess_manual_moments_are_rejected_after_expansion(
    moments, provided, atoms, excess
):
    values = spin_ops.parse_magmoms(moments)
    noun = "value" if excess == 1 else "values"
    message = (
        rf"{provided} magnetic moments were provided for a {atoms}-atom "
        rf"structure; remove {excess} excess {noun}"
    )
    with pytest.raises(ValueError, match=message):
        spin_ops.fit_magmoms_to_structure(values, atoms)


def test_magnetic_primitive_loader_uses_the_same_count_validation(tmp_path):
    poscar = tmp_path / "POSCAR"
    _write_two_atom_poscar(poscar)

    with pytest.raises(
        ValueError,
        match="3 magnetic moments were provided for a 2-atom structure",
    ):
        _load_magnetic_input_data(str(poscar), "2*1 1*-1", "0 0 1")

    _, _, _, moments, _ = _load_magnetic_input_data(
        str(poscar), "1", "0 0 1"
    )
    assert np.allclose(moments, [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])


def test_spin_symmetry_route_reports_excess_count(tmp_path, monkeypatch, capsys):
    poscar = tmp_path / "POSCAR"
    _write_two_atom_poscar(poscar)
    monkeypatch.setattr(
        spin_ops,
        "_non_magnetic_symmetry",
        lambda *args, **kwargs: {
            "non_mag_label": "Im-3m (229)",
            "point_group": "m-3m",
            "laue_group": "m-3m",
            "spacegroup_number": 229,
            "sites": 1,
            "lattice": "cI",
        },
    )

    result = spin_ops.run(
        str(poscar),
        "1 -1 7",
        verbose=False,
        spin_axis_cart="0 0 1",
        output_dir=str(tmp_path / "output"),
    )

    assert result is False
    assert "3 magnetic moments were provided for a 2-atom structure" in (
        capsys.readouterr().out
    )
    assert not (tmp_path / "output" / "spin_operations.txt").exists()
