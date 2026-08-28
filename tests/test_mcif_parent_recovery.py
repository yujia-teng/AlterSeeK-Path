import numpy as np

from alterseek.compute_centroid_3d import _select_mcif_parent_symprec
from alterseek.mcif import _parent_hint_from_cif_block


def test_parent_hint_reads_declared_child_transform():
    hint = _parent_hint_from_cif_block({
        "_parent_space_group.child_transform_Pp_abc": "2a+b,-a+b,c;0,0,0",
        "_parent_space_group.IT_number": "194",
        "_parent_space_group.name_H-M_alt": "P 6_3/m m c",
    })

    assert hint == {
        "index": 3,
        "spacegroup_number": 194,
        "spacegroup_symbol": "P 6_3/m m c",
        "transform": "2a+b,-a+b,c;0,0,0",
    }


def test_parent_symprec_uses_smallest_value_matching_declared_parent(monkeypatch):
    lattice = np.diag([2.0, 1.0, 1.0])
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.5001, 0.0, 0.0],
    ])
    numbers = [1, 1]
    hint = {
        "index": 2,
        "spacegroup_number": 221,
        "spacegroup_symbol": "P m -3 m",
        "transform": "2a,b,c;0,0,0",
    }
    monkeypatch.setattr(
        "alterseek.compute_centroid_3d._declared_mcif_parent_hint",
        lambda _filename: hint,
    )

    symprec, recovered = _select_mcif_parent_symprec(
        "rounded-supercell.mcif",
        lattice,
        positions,
        numbers,
    )

    assert symprec == 1e-3
    assert recovered["primitive_sites"] == 1
    assert recovered["input_sites"] == 2
    assert recovered["detected_spacegroup_number"] == 221
