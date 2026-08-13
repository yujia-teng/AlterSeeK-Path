"""Build validated submitted-cell marker inputs for SeeK-path."""

import os
import warnings

import numpy as np
import seekpath
import spglib
from findspingroup import find_spin_group_acc_primitive_from_data

from .io import (
    _dedupe_frac_positions,
    _group_poscar_sites,
    _load_magnetic_input_data,
    _min_periodic_cart_distance,
    _write_magnetic_mcif,
)


_MARKER_SEEDS = [
    np.array([0.11000000, 0.12000000, 0.15000001]),
    np.array([0.13000000, 0.17000000, 0.23000001]),
    np.array([0.07100000, 0.19300000, 0.31700001]),
    np.array([0.21100000, 0.13700000, 0.29300001]),
]

_MARKER_SEED_SETS = (
    (_MARKER_SEEDS[0], _MARKER_SEEDS[1]),
    (_MARKER_SEEDS[0], _MARKER_SEEDS[2]),
    (_MARKER_SEEDS[1], _MARKER_SEEDS[3]),
    (_MARKER_SEEDS[0], _MARKER_SEEDS[1], _MARKER_SEEDS[2]),
)


def _wrapped_translation(translation, tol=1e-7):
    wrapped = np.mod(np.asarray(translation, dtype=float), 1.0)
    wrapped[np.isclose(wrapped, 0.0, atol=tol, rtol=0.0)] = 0.0
    wrapped[np.isclose(wrapped, 1.0, atol=tol, rtol=0.0)] = 0.0
    return wrapped


def _space_operation_key(rotation, translation, tol=1e-7):
    decimals = max(0, int(np.ceil(-np.log10(tol))))
    return (
        tuple(np.rint(rotation).astype(int).ravel()),
        tuple(np.round(_wrapped_translation(translation, tol), decimals)),
    )


def _validated_space_operations(space_operations, tol=1e-7):
    """Return distinct integral Seitz operations in the submitted basis."""
    operations = []
    keys = set()
    for operation in space_operations:
        if isinstance(operation, dict):
            value = operation.get("real_rotation")
            translation = operation.get("translation", np.zeros(3))
        else:
            value = operation
            translation = np.zeros(3)
        rotation = np.asarray(value, dtype=float)
        translation = np.asarray(translation, dtype=float)
        if rotation.shape != (3, 3) or not np.all(np.isfinite(rotation)):
            continue
        if translation.shape != (3,) or not np.all(np.isfinite(translation)):
            continue
        rounded = np.rint(rotation).astype(int)
        if not np.allclose(rotation, rounded, atol=tol):
            continue
        if not np.isclose(abs(np.linalg.det(rounded)), 1.0, atol=tol):
            continue
        wrapped = _wrapped_translation(translation, tol)
        key = _space_operation_key(rounded, wrapped, tol)
        if key not in keys:
            keys.add(key)
            operations.append({
                "real_rotation": rounded,
                "translation": wrapped,
            })
    if not operations:
        raise RuntimeError(
            "No spatial symmetry operation is compatible with the submitted "
            "cell lattice."
        )
    return operations


def _point_operation_keys(rotations):
    return {
        tuple(np.asarray(rotation, dtype=int).ravel())
        for rotation in rotations
    }


def _space_operation_keys(rotations, translations, tol=1e-7):
    return {
        _space_operation_key(rotation, translation, tol)
        for rotation, translation in zip(rotations, translations)
    }


def _one_full_operation_per_rotation(space_operations):
    """Remove supercell translation copies without discarding Seitz shifts."""
    selected = []
    rotation_keys = set()
    for operation in _validated_space_operations(space_operations):
        rotation_key = tuple(operation["real_rotation"].ravel())
        if rotation_key not in rotation_keys:
            rotation_keys.add(rotation_key)
            selected.append(operation)
    return selected


def _magnetic_primitive_nssg_operations(result):
    views = (
        result.get("operation_views", {})
        .get("magnetic_primitive_cartesian", {})
        .get("views", {})
    )
    nssg = views.get("nssg", {})
    operations = nssg.get("ops", []) if isinstance(nssg, dict) else []
    if not operations:
        raise RuntimeError(
            "FindSpinGroup did not return magnetic-primitive nssg operations."
        )
    return operations


def _magnetic_operations_in_submitted_basis(result):
    """Transform FindSpinGroup's full primitive Seitz operations to input."""
    transform = result.get("T_input_to_acc_primitive")
    if not transform or len(transform) != 2:
        raise RuntimeError(
            "FindSpinGroup did not return the input-to-magnetic-primitive "
            "setting transformation."
        )
    input_to_primitive = np.asarray(transform[0], dtype=float)
    origin_shift = np.asarray(transform[1], dtype=float)
    if input_to_primitive.shape != (3, 3) or origin_shift.shape != (3,):
        raise RuntimeError(
            "FindSpinGroup returned an invalid magnetic-primitive setting "
            "transformation."
        )
    primitive_to_input = np.linalg.inv(input_to_primitive)
    operations = []
    for operation in _magnetic_primitive_nssg_operations(result):
        rotation_primitive = np.asarray(
            operation["real_rotation"], dtype=float
        )
        translation_primitive = np.asarray(
            operation.get("translation", np.zeros(3)), dtype=float
        )
        rotation_input = (
            primitive_to_input
            @ rotation_primitive
            @ input_to_primitive
        )
        translation_input = primitive_to_input @ (
            rotation_primitive @ origin_shift
            + translation_primitive
            - origin_shift
        )
        operations.append({
            "real_rotation": rotation_input,
            "translation": translation_input,
        })
    return _validated_space_operations(operations)


def _seekpath_lattice_tag(lattice, symprec):
    """Return the HPKOT tag of a primitive lattice from its metric alone."""
    metric_cell = (
        np.asarray(lattice, dtype=float).tolist(),
        [[0.11000000, 0.12000000, 0.15000001]],
        [1],
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=r".*dict interface is deprecated.*"
        )
        warnings.filterwarnings(
            "ignore",
            category=DeprecationWarning,
            module=r"seekpath\.hpkot(\..*)?",
        )
        result = seekpath.get_path(
            metric_cell,
            with_time_reversal=True,
            symprec=symprec,
        )
    return result["bravais_lattice_extended"]


def build_submitted_analysis_cell(
    lattice,
    real_positions,
    real_type_numbers,
    space_operations,
    *,
    symprec=1e-3,
):
    """Build and validate an in-memory marker cell in the submitted lattice.

    This is the old marker-helper construction in memory: every generic seed
    is transformed by the complete compatible Seitz operations ``R r + t``.
    The artificial type replaces He chemistry; the symmetry construction is
    otherwise unchanged.
    """
    lattice = np.asarray(lattice, dtype=float)
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise ValueError("Submitted lattice must be a finite 3x3 matrix.")
    if abs(np.linalg.det(lattice)) < 1e-12:
        raise ValueError("Submitted lattice is singular.")

    real_positions = np.mod(np.asarray(real_positions, dtype=float), 1.0)
    real_type_numbers = [int(value) for value in real_type_numbers]
    if real_positions.shape != (len(real_type_numbers), 3):
        raise ValueError(
            "Real positions and type numbers must describe the same sites."
        )
    if not np.all(np.isfinite(real_positions)):
        raise ValueError("Real positions must be finite fractional coordinates.")
    real_types = set(real_type_numbers)
    marker_type = max(real_types, default=0) + 1
    while marker_type in real_types or marker_type <= 0:
        marker_type += 1

    operations = _validated_space_operations(space_operations)
    rotations = [operation["real_rotation"] for operation in operations]
    translations = [operation["translation"] for operation in operations]
    intended_keys = _point_operation_keys(rotations)
    intended_space_keys = _space_operation_keys(rotations, translations)
    failures = []
    for seeds in _MARKER_SEED_SETS:
        markers = _dedupe_frac_positions([
            seed @ rotation.T + translation
            for seed in seeds
            for rotation, translation in zip(rotations, translations)
        ])
        helper_positions = [*real_positions.tolist(), *markers]
        cell = (
            lattice.tolist(),
            [np.asarray(position, dtype=float).tolist()
             for position in helper_positions],
            [*real_type_numbers, *([marker_type] * len(markers))],
        )
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        seed_label = [seed.tolist() for seed in seeds]
        if dataset is None:
            failures.append(f"seeds {seed_label}: spglib found no symmetry")
            continue
        detected_keys = _point_operation_keys(dataset.rotations)
        detected_space_keys = _space_operation_keys(
            dataset.rotations, dataset.translations
        )
        if detected_space_keys != intended_space_keys:
            failures.append(
                f"seeds {seed_label}: intended "
                f"{len(intended_space_keys)} full spatial operations but "
                f"detected {len(detected_space_keys)}"
            )
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message=r".*dict interface is deprecated.*"
            )
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                module=r"seekpath\.hpkot(\..*)?",
            )
            sp_result = seekpath.get_path(
                cell,
                with_time_reversal=True,
                symprec=symprec,
            )
        volume_ratio = float(sp_result["volume_original_wrt_prim"])
        if not np.isclose(volume_ratio, 1.0, atol=1e-7, rtol=0.0):
            failures.append(
                f"seeds {seed_label}: volume_original_wrt_prim="
                f"{volume_ratio:.12g}"
            )
            continue
        reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
        return {
            "cell": cell,
            "marker_type": marker_type,
            "marker_seeds": [seed.copy() for seed in seeds],
            "marker_count": len(markers),
            "marker_min_distance": _min_periodic_cart_distance(
                helper_positions, lattice
            ),
            "point_operations": rotations,
            "space_operations": operations,
            "intended_point_operation_count": len(intended_keys),
            "detected_point_operation_count": len(detected_keys),
            "intended_space_operation_count": len(intended_space_keys),
            "detected_space_operation_count": len(detected_space_keys),
            "volume_original_wrt_prim": volume_ratio,
            "submitted_bz_volume": abs(float(np.linalg.det(reciprocal))),
            "seekpath_bravais": sp_result["bravais_lattice_extended"],
            "point_coords": dict(sp_result["point_coords"]),
            "analysis_spacegroup_number": int(dataset.number),
            "analysis_spacegroup_symbol": str(dataset.international),
            "analysis_point_group": str(dataset.pointgroup),
        }

    raise RuntimeError(
        "Could not validate a submitted-cell marker helper without primitive "
        f"reduction. {'; '.join(failures)}"
    )


def prepare_submitted_cell_analysis(
    structure_file,
    *,
    moments_str="",
    spin_axis_cart=None,
    output_dir=".",
    symprec=1e-3,
    write_magnetic_diagnostic=False,
):
    """Prepare the validated in-memory analysis cell in the submitted lattice."""
    from ase.data import atomic_numbers

    lattice, positions, elements, moments, spin_setting = (
        _load_magnetic_input_data(structure_file, moments_str, spin_axis_cart)
    )
    real_types = [atomic_numbers[str(element)] for element in elements]
    fsg_result = None
    has_magnetic_moments = bool(np.any(
        np.linalg.norm(np.asarray(moments, dtype=float), axis=1) > 1e-10
    ))
    if has_magnetic_moments:
        fsg_result = find_spin_group_acc_primitive_from_data(
            structure_file,
            lattice,
            positions,
            elements,
            [1.0] * len(elements),
            moments,
            input_spin_setting=spin_setting,
        )
        space_operations = _one_full_operation_per_rotation(
            _magnetic_operations_in_submitted_basis(fsg_result)
        )
    else:
        symmetry = spglib.get_symmetry(
            (lattice, positions, real_types),
            symprec=symprec,
        )
        if symmetry is None:
            raise RuntimeError(
                "Could not determine submitted-cell structural operations."
            )
        space_operations = _one_full_operation_per_rotation([
            {
                "real_rotation": rotation,
                "translation": translation,
            }
            for rotation, translation in zip(
                symmetry["rotations"], symmetry["translations"]
            )
        ])

    helper = build_submitted_analysis_cell(
        lattice,
        positions,
        real_types,
        space_operations,
        symprec=symprec,
    )
    basename = os.path.splitext(os.path.basename(structure_file))[0]
    result = {
        "analysis_cell": helper["cell"],
        "analysis_marker_type": helper["marker_type"],
        "submitted_lattice": np.array(lattice, dtype=float),
        "submitted_sites": len(elements),
        "operation_basis_label": (
            f"submitted structure '{os.path.basename(structure_file)}'"
        ),
        "analysis_symmetry": {
            "number": helper["analysis_spacegroup_number"],
            "symbol": helper["analysis_spacegroup_symbol"],
            "point_group": helper["analysis_point_group"],
            "seekpath_bravais": helper["seekpath_bravais"],
        },
        "summary": {
            "marker_seeds": [
                seed.tolist() for seed in helper["marker_seeds"]
            ],
            "marker_count": helper["marker_count"],
            "marker_type": helper["marker_type"],
            "marker_min_distance": helper["marker_min_distance"],
            "intended_point_operations": helper[
                "intended_point_operation_count"
            ],
            "detected_point_operations": helper[
                "detected_point_operation_count"
            ],
            "intended_space_operations": helper[
                "intended_space_operation_count"
            ],
            "detected_space_operations": helper[
                "detected_space_operation_count"
            ],
            "volume_original_wrt_prim": helper[
                "volume_original_wrt_prim"
            ],
        },
    }

    if write_magnetic_diagnostic:
        if fsg_result is None:
            raise RuntimeError(
                "A magnetic primitive diagnostic was requested for a "
                "structure without nonzero magnetic moments."
            )
        magnetic_cell = fsg_result["acc_primitive_cell_detail"]
        magnetic_lattice = np.asarray(
            magnetic_cell["lattice"], dtype=float
        )
        magnetic_positions = [
            np.asarray(position, dtype=float)
            for position in magnetic_cell["positions"]
        ]
        magnetic_elements = [
            str(value) for value in magnetic_cell["elements"]
        ]
        magnetic_moments = np.asarray(
            magnetic_cell.get(
                "moments", np.zeros((len(magnetic_elements), 3))
            ),
            dtype=float,
        )
        _, _, grouped_positions, ordered_indices = _group_poscar_sites(
            magnetic_elements,
            magnetic_positions,
        )
        ordered_elements = [
            magnetic_elements[index] for index in ordered_indices
        ]
        ordered_moments = magnetic_moments[ordered_indices]
        os.makedirs(output_dir, exist_ok=True)
        mcif_path = os.path.join(
            output_dir, f"{basename}_magnetic_primitive.mcif"
        )
        _write_magnetic_mcif(
            mcif_path,
            f"{basename}_magnetic_primitive",
            magnetic_lattice,
            ordered_elements,
            grouped_positions,
            ordered_moments,
        )
        result.update({
            "mcif_path": mcif_path,
            "magnetic_primitive_lattice": magnetic_lattice,
            "magnetic_primitive_sites": len(magnetic_elements),
            "magnetic_primitive_lattice_tag": _seekpath_lattice_tag(
                magnetic_lattice, symprec
            ),
            "magnetic_summary": {
                "index": fsg_result.get("index"),
                "acc_symbol": fsg_result.get("acc_symbol"),
                "setting": fsg_result.get("acc_primitive_cell_setting"),
            },
        })
    return result
