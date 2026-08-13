"""The submitted translation lattice defines the analyzed Brillouin zone."""

import itertools

import numpy as np
import spglib
from ase import Atoms
from ase.io import write

from alterseek.compute_centroid_hybrid import run as compute_centroid
from alterseek.ssg_setting import prepare_submitted_cell_analysis
from alterseek.ssg_setting import build_submitted_analysis_cell


def _signed_permutation_operations():
    operations = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            rotation = np.zeros((3, 3), dtype=int)
            for row, column in enumerate(permutation):
                rotation[row, column] = signs[row]
            operations.append(rotation)
    return operations


def _orthorhombic_operations():
    return [
        np.diag(signs).astype(int)
        for signs in itertools.product((-1, 1), repeat=3)
    ]


def _change_operation_basis(rotation, source_lattice, target_lattice):
    cartesian = source_lattice.T @ rotation @ np.linalg.inv(source_lattice.T)
    transformed = (
        np.linalg.inv(target_lattice.T) @ cartesian @ target_lattice.T
    )
    return np.rint(transformed).astype(int)


def test_primitive_input_is_not_reduced():
    lattice = np.diag([4.0, 5.0, 6.0])
    result = build_submitted_analysis_cell(
        lattice,
        real_type_numbers=[26, 26],
        point_operations=_orthorhombic_operations(),
    )

    assert np.allclose(result["cell"][0], lattice)
    assert result["volume_original_wrt_prim"] == 1.0
    assert result["intended_point_operation_count"] == 8
    assert result["detected_point_operation_count"] == 8


def test_true_222_supercell_keeps_its_native_bz():
    lattice = np.diag([8.0, 8.0, 8.0])
    result = build_submitted_analysis_cell(
        lattice,
        real_type_numbers=[26] * 8,
        point_operations=_signed_permutation_operations(),
    )

    assert np.allclose(result["cell"][0], lattice)
    assert result["volume_original_wrt_prim"] == 1.0
    assert result["seekpath_bravais"] == "cP2"
    assert {"GAMMA", "X", "M", "R"}.issubset(result["point_coords"])
    expected_bz_volume = (2 * np.pi) ** 3 / abs(np.linalg.det(lattice))
    assert np.isclose(result["submitted_bz_volume"], expected_bz_volume)


def test_determinant_one_basis_change_keeps_submitted_lattice():
    orthorhombic = np.diag([4.0, 5.0, 6.0])
    transform = np.array([
        [0, -1, 0],
        [1, 1, 0],
        [0, 0, 1],
    ])
    submitted = transform @ orthorhombic
    submitted_operations = [
        _change_operation_basis(rotation, orthorhombic, submitted)
        for rotation in _orthorhombic_operations()
    ]

    result = build_submitted_analysis_cell(
        submitted,
        real_type_numbers=[64, 79, 32],
        point_operations=submitted_operations,
    )

    assert round(np.linalg.det(transform)) == 1
    assert np.allclose(result["cell"][0], submitted)
    assert result["volume_original_wrt_prim"] == 1.0
    assert result["seekpath_bravais"].startswith("o")


def test_artificial_type_is_positive_and_distinct_from_real_types():
    result = build_submitted_analysis_cell(
        np.array([
            [4.0, 0.0, 0.0],
            [0.3, 5.0, 0.0],
            [0.2, 0.4, 6.0],
        ]),
        real_type_numbers=[1, 2, 3, 118],
        point_operations=[np.eye(3, dtype=int)],
    )

    marker_types = set(result["cell"][2])
    assert marker_types == {result["marker_type"]}
    assert result["marker_type"] > 0
    assert result["marker_type"] not in {1, 2, 3, 118}


def test_anisotropic_supercell_keeps_its_full_submitted_volume():
    lattice = np.diag([8.0, 15.0, 6.0])
    result = build_submitted_analysis_cell(
        lattice,
        real_type_numbers=[26] * 6,
        point_operations=_orthorhombic_operations(),
    )

    assert result["seekpath_bravais"] == "oP1"
    assert result["volume_original_wrt_prim"] == 1.0
    assert np.isclose(
        result["submitted_bz_volume"],
        (2 * np.pi) ** 3 / np.linalg.det(lattice),
    )


def test_q_nonzero_enlarged_cell_is_not_returned_to_a_parent_period():
    lattice = np.diag([12.0, 4.0, 4.0])
    result = build_submitted_analysis_cell(
        lattice,
        real_type_numbers=[25, 25, 25],
        point_operations=_orthorhombic_operations(),
    )

    assert np.allclose(result["cell"][0], lattice)
    assert result["volume_original_wrt_prim"] == 1.0
    assert np.isclose(
        result["submitted_bz_volume"],
        (2 * np.pi) ** 3 / (12.0 * 4.0 * 4.0),
    )


def test_rhombohedral_primitive_and_hexagonal_conventional_inputs_keep_distinct_bzs():
    a_hex = np.array([
        [5.0, 0.0, 0.0],
        [-2.5, 2.5 * np.sqrt(3.0), 0.0],
        [0.0, 0.0, 13.0],
    ])
    to_rhombohedral = np.array([
        [2 / 3, 1 / 3, 1 / 3],
        [-1 / 3, 1 / 3, 1 / 3],
        [-1 / 3, -2 / 3, 1 / 3],
    ])
    a_rhombohedral = to_rhombohedral @ a_hex
    hex_positions = np.array([
        [0.0, 0.0, 0.0],
        [2 / 3, 1 / 3, 1 / 3],
        [1 / 3, 2 / 3, 2 / 3],
    ])

    def structural_rotations(lattice, positions):
        symmetry = spglib.get_symmetry(
            (lattice, positions, [26] * len(positions)),
            symprec=1e-5,
        )
        return symmetry["rotations"]

    primitive = build_submitted_analysis_cell(
        a_rhombohedral,
        real_type_numbers=[26],
        point_operations=structural_rotations(
            a_rhombohedral, np.zeros((1, 3))
        ),
        symprec=1e-5,
    )
    conventional = build_submitted_analysis_cell(
        a_hex,
        real_type_numbers=[26, 26, 26],
        point_operations=structural_rotations(a_hex, hex_positions),
        symprec=1e-5,
    )

    assert np.isclose(np.linalg.det(a_hex) / np.linalg.det(a_rhombohedral), 3)
    assert np.isclose(primitive["volume_original_wrt_prim"], 1.0)
    assert np.isclose(conventional["volume_original_wrt_prim"], 1.0)
    assert np.isclose(
        primitive["submitted_bz_volume"],
        3 * conventional["submitted_bz_volume"],
    )


def test_no_moments_supercell_uses_the_same_submitted_cell_helper(tmp_path):
    positions = list(itertools.product((0.0, 0.5), repeat=3))
    structure = Atoms(
        symbols=["Fe"] * 8,
        scaled_positions=positions,
        cell=np.eye(3) * 8.0,
        pbc=True,
    )
    path = tmp_path / "POSCAR"
    write(path, structure, format="vasp", direct=True, vasp5=True)

    preparation = prepare_submitted_cell_analysis(str(path))
    result = compute_centroid(
        str(path),
        output_dir=str(tmp_path),
        show_plot=False,
        verbose=False,
        analysis_cell=preparation["analysis_cell"],
        analysis_marker_type=preparation["analysis_marker_type"],
    )

    assert result["sc_type"] == "cP2"
    assert np.allclose(result["b_matrix_input"], np.eye(3) * np.pi / 4)
    assert not (tmp_path / "POSCAR_seekpath_standard.vasp").exists()
    assert (tmp_path / "POSCAR_seekpath_basis_mapping.txt").exists()
    assert not (tmp_path / "POSCAR_magnetic_primitive.mcif").exists()
