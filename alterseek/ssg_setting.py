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


def _compatible_point_operations(lattice, space_operations, tol=1e-7):
    """Return the closed point group compatible with the submitted lattice.

    The submitted fractional basis represents a lattice automorphism only when
    its rotation matrix is integral and unimodular.  ``_validated_space_operations``
    performs that compatibility filter.  Translation columns are deliberately
    discarded here: they remain part of the physical symmetry data, but the
    conventional/supercell BZ proxy is the symmorphic group ``(R, 0)`` built on
    the submitted edge-vector translation lattice.
    """
    lattice = np.asarray(lattice, dtype=float)
    metric = lattice @ lattice.T
    candidates = _validated_space_operations(space_operations, tol=tol)
    compatible = []
    metric_atol = tol * max(1.0, float(np.max(np.abs(metric))))
    for operation in candidates:
        rotation = np.asarray(operation["real_rotation"], dtype=int)
        if np.allclose(
            rotation.T @ metric @ rotation,
            metric,
            atol=metric_atol,
            rtol=0.0,
        ):
            compatible.append(operation)
    if not compatible:
        raise RuntimeError(
            "No point operation preserves the submitted cell metric."
        )
    rotations = {}
    for operation in compatible:
        rotation = np.asarray(operation["real_rotation"], dtype=int)
        rotations.setdefault(tuple(rotation.ravel()), rotation)

    identity_key = tuple(np.eye(3, dtype=int).ravel())
    if identity_key not in rotations:
        raise RuntimeError(
            "Submitted-cell compatible point operations do not contain the "
            "identity."
        )

    keys = set(rotations)
    for left in rotations.values():
        for right in rotations.values():
            product_key = tuple((left @ right).ravel())
            if product_key not in keys:
                raise RuntimeError(
                    "Submitted-cell compatible rotations do not form a "
                    "closed point group."
                )
    return list(rotations.values()), compatible, candidates


def _database_operations_in_input_basis(
    hall_number,
    input_to_standard,
    origin_shift,
    *,
    tol=1e-6,
):
    """Transform one complete standard Hall operation set to the input."""
    input_to_standard = np.asarray(input_to_standard, dtype=float)
    origin_shift = np.asarray(origin_shift, dtype=float)
    if (
        input_to_standard.shape != (3, 3)
        or origin_shift.shape != (3,)
        or not np.all(np.isfinite(input_to_standard))
        or not np.all(np.isfinite(origin_shift))
        or abs(np.linalg.det(input_to_standard)) < 1e-12
    ):
        raise RuntimeError("Invalid input-to-standard setting transformation.")

    database = spglib.get_symmetry_from_database(int(hall_number))
    if database is None:
        raise RuntimeError(
            f"spglib has no operation set for Hall number {hall_number}."
        )
    standard_to_input = np.linalg.inv(input_to_standard)
    transformed = []
    for rotation_standard, translation_standard in zip(
        database["rotations"], database["translations"]
    ):
        rotation_standard = np.asarray(rotation_standard, dtype=float)
        translation_standard = np.asarray(translation_standard, dtype=float)
        rotation_input = (
            standard_to_input
            @ rotation_standard
            @ input_to_standard
        )
        if not np.allclose(
            rotation_input, np.rint(rotation_input), atol=tol, rtol=0.0
        ):
            raise RuntimeError(
                f"Hall {hall_number} produces a nonintegral input rotation."
            )
        translation_input = standard_to_input @ (
            rotation_standard @ origin_shift
            + translation_standard
            - origin_shift
        )
        transformed.append({
            "real_rotation": rotation_input,
            "translation": translation_input,
        })
    return _validated_space_operations(transformed, tol=tol)


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
    """Transform FindSpinGroup's primitive Seitz representatives to input."""
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


def _g0_spacegroup_number(result):
    details = result.get("identify_index_details")
    number = details.get("G0_id") if isinstance(details, dict) else None
    try:
        number = int(number)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            "FindSpinGroup did not return a valid G0 space-group number."
        ) from exc
    if not 1 <= number <= 230:
        raise RuntimeError(
            f"FindSpinGroup returned invalid G0 space-group number {number}."
        )
    return number


def _input_to_g0_standard_transform(result):
    """Compose input -> magnetic primitive -> G0 standard coordinates."""
    transforms = []
    for key in ("T_input_to_acc_primitive", "T_acc_primitive_to_G0std"):
        transform = result.get(key)
        if not transform or len(transform) != 2:
            raise RuntimeError(
                f"FindSpinGroup did not return the required {key} transform."
            )
        matrix = np.asarray(transform[0], dtype=float)
        shift = np.asarray(transform[1], dtype=float)
        if (
            matrix.shape != (3, 3)
            or shift.shape != (3,)
            or not np.all(np.isfinite(matrix))
            or not np.all(np.isfinite(shift))
            or abs(np.linalg.det(matrix)) < 1e-12
        ):
            raise RuntimeError(
                f"FindSpinGroup returned an invalid {key} transform."
            )
        transforms.append((matrix, shift))

    (input_to_primitive, input_shift), (
        primitive_to_standard,
        primitive_shift,
    ) = transforms
    return (
        primitive_to_standard @ input_to_primitive,
        primitive_to_standard @ input_shift + primitive_shift,
    )


def _complete_magnetic_operations_in_submitted_basis(result, tol=1e-6):
    """Return the complete G0 Seitz set in the submitted setting.

    FindSpinGroup's magnetic-primitive view contains one operation per
    primitive-cell representative.  Transforming that list alone cannot
    recover a centered conventional target.  Use those representatives only
    to identify the matching standard Hall setting, then transform the full
    crystallographic database operation set into the submitted basis.
    """
    g0_number = _g0_spacegroup_number(result)
    representatives = _magnetic_operations_in_submitted_basis(result)
    representative_keys = _space_operation_keys(
        [operation["real_rotation"] for operation in representatives],
        [operation["translation"] for operation in representatives],
        tol=tol,
    )
    input_to_standard, origin_shift = _input_to_g0_standard_transform(result)
    distinct_candidates = {}
    failures = []
    for hall_number in range(1, 531):
        spacegroup_type = spglib.get_spacegroup_type(hall_number)
        if (
            spacegroup_type is None
            or int(spacegroup_type.number) != g0_number
        ):
            continue
        try:
            operations = _database_operations_in_input_basis(
                hall_number,
                input_to_standard,
                origin_shift,
                tol=tol,
            )
        except RuntimeError as exc:
            failures.append(str(exc))
            continue
        operation_keys = _space_operation_keys(
            [operation["real_rotation"] for operation in operations],
            [operation["translation"] for operation in operations],
            tol=tol,
        )
        if not representative_keys.issubset(operation_keys):
            failures.append(
                f"Hall {hall_number}: does not contain the transformed "
                "magnetic-primitive representatives"
            )
            continue
        distinct_candidates.setdefault(
            frozenset(operation_keys), (hall_number, operations)
        )

    if not distinct_candidates:
        raise RuntimeError(
            "Could not reconstruct the complete G0 operation set in the "
            f"submitted basis for space group {g0_number}. "
            f"{' ; '.join(failures)}"
        )
    if len(distinct_candidates) != 1:
        hall_numbers = sorted(
            hall_number
            for hall_number, _operations in distinct_candidates.values()
        )
        raise RuntimeError(
            "Multiple inequivalent standard Hall settings match the "
            f"submitted G0 operation representatives: {hall_numbers}."
        )
    return next(iter(distinct_candidates.values()))[1]


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


def _physical_symmetry_from_operations(
    lattice,
    operations,
    expected_number,
    symprec,
):
    """Identify the unchanged physical Seitz set separately from the BZ proxy."""
    spacegroup_type = spglib.get_spacegroup_type_from_symmetry(
        [operation["real_rotation"] for operation in operations],
        [operation["translation"] for operation in operations],
        lattice=np.asarray(lattice, dtype=float),
        symprec=symprec,
    )
    if spacegroup_type is None:
        raise RuntimeError(
            "Could not identify the complete physical space-group operation "
            "set in the submitted basis."
        )
    if int(spacegroup_type.number) != int(expected_number):
        raise RuntimeError(
            "Physical operation-set identification changed the expected "
            f"space-group number {int(expected_number)} to "
            f"{int(spacegroup_type.number)}."
        )
    return {
        "number": int(spacegroup_type.number),
        "symbol": str(spacegroup_type.international_short),
        "point_group": str(spacegroup_type.pointgroup_international),
        "hall_number": int(spacegroup_type.hall_number),
    }


def _standard_physical_symmetry(spacegroup_number):
    """Return a conventional database label when a supercell hides the setting."""
    fallback = None
    for hall_number in range(1, 531):
        spacegroup_type = spglib.get_spacegroup_type(hall_number)
        if (
            spacegroup_type is None
            or int(spacegroup_type.number) != int(spacegroup_number)
        ):
            continue
        candidate = {
            "number": int(spacegroup_type.number),
            "symbol": str(spacegroup_type.international_short),
            "point_group": str(spacegroup_type.pointgroup_international),
            "hall_number": int(spacegroup_type.hall_number),
        }
        if fallback is None:
            fallback = candidate
        if str(spacegroup_type.choice) == "":
            return candidate
    if fallback is None:
        raise RuntimeError(
            f"No spglib setting exists for space group {spacegroup_number}."
        )
    return fallback


def build_submitted_analysis_cell(
    lattice,
    real_positions,
    real_type_numbers,
    space_operations,
    *,
    symprec=1e-3,
):
    """Build the conventional/supercell-BZ proxy in the submitted lattice.

    The proxy is marker-only.  Real sites generally obey nonsymmorphic
    ``(R, t)`` operations and therefore cannot be mixed into an artificial
    symmorphic ``(R, 0)`` orbit.  The physical sites and complete Seitz set are
    retained by the caller for magnetic and diagnostic uses.
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

    rotations, compatible_space_operations, source_space_operations = (
        _compatible_point_operations(
            lattice, space_operations
        )
    )
    translations = [np.zeros(3) for _rotation in rotations]
    operations = [
        {"real_rotation": rotation, "translation": np.zeros(3)}
        for rotation in rotations
    ]
    intended_keys = _point_operation_keys(rotations)
    intended_space_keys = _space_operation_keys(rotations, translations)
    failures = []
    for seeds in _MARKER_SEED_SETS:
        helper_positions = []
        helper_types = []
        marker_types = []
        used_types = set(real_types)
        next_marker_type = marker_type
        for seed in seeds:
            orbit = _dedupe_frac_positions([
                seed @ rotation.T for rotation in rotations
            ])
            while next_marker_type in used_types or next_marker_type <= 0:
                next_marker_type += 1
            orbit_type = next_marker_type
            used_types.add(orbit_type)
            next_marker_type += 1
            marker_types.append(orbit_type)
            helper_positions.extend(orbit)
            helper_types.extend([orbit_type] * len(orbit))
        cell = (
            lattice.tolist(),
            [np.asarray(position, dtype=float).tolist()
             for position in helper_positions],
            helper_types,
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
        if not np.isfinite(volume_ratio) or volume_ratio <= 0.0:
            failures.append(
                f"seeds {seed_label}: invalid volume_original_wrt_prim="
                f"{volume_ratio!r}"
            )
            continue
        if not np.isclose(volume_ratio, 1.0, atol=1e-6, rtol=0.0):
            failures.append(
                f"seeds {seed_label}: helper reduced the submitted "
                f"translation lattice, volume_original_wrt_prim="
                f"{volume_ratio!r}"
            )
            continue
        reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
        return {
            "cell": cell,
            "marker_type": marker_type,
            "marker_types": marker_types,
            "marker_seeds": [seed.copy() for seed in seeds],
            "marker_count": len(helper_positions),
            "marker_min_distance": _min_periodic_cart_distance(
                helper_positions, lattice
            ),
            "point_operations": rotations,
            "space_operations": operations,
            "source_space_operations": source_space_operations,
            "source_space_operation_count": len(source_space_operations),
            "compatible_space_operation_count": len(
                compatible_space_operations
            ),
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
        "Could not validate the submitted-cell marker helper. "
        f"{'; '.join(failures)}"
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
        expected_spacegroup_number = _g0_spacegroup_number(fsg_result)
        # Magnetic-primitive representatives contain every G0 point part.
        # Transform them directly and keep only rotations compatible with the
        # submitted lattice; arbitrary supercells need not admit all parent
        # rotations as integral lattice automorphisms.
        space_operations = _magnetic_operations_in_submitted_basis(fsg_result)
        try:
            physical_operations = (
                _complete_magnetic_operations_in_submitted_basis(fsg_result)
            )
            physical_symmetry = _physical_symmetry_from_operations(
                lattice,
                physical_operations,
                expected_spacegroup_number,
                symprec,
            )
            physical_operation_set_verified = True
        except RuntimeError:
            # The complete conventional Hall set need not embed integrally in
            # an anisotropic or non-diagonal supercell.  This does not affect
            # the BZ helper, which needs only compatible point representatives.
            physical_operations = space_operations
            physical_symmetry = _standard_physical_symmetry(
                expected_spacegroup_number
            )
            physical_operation_set_verified = False
    else:
        dataset = spglib.get_symmetry_dataset(
            (lattice, positions, real_types),
            symprec=symprec,
        )
        if dataset is None:
            raise RuntimeError(
                "Could not determine submitted-cell structural operations."
            )
        expected_spacegroup_number = int(dataset.number)
        symmetry = spglib.get_symmetry(
            (lattice, positions, real_types),
            symprec=symprec,
        )
        if symmetry is None:
            raise RuntimeError(
                "Could not determine submitted-cell structural operations."
            )
        space_operations = _validated_space_operations([
            {
                "real_rotation": rotation,
                "translation": translation,
            }
            for rotation, translation in zip(
                symmetry["rotations"], symmetry["translations"]
            )
        ])
        physical_operations = space_operations
        physical_operation_set_verified = True
        physical_symmetry = {
            "number": int(dataset.number),
            "symbol": str(dataset.international),
            "point_group": str(dataset.pointgroup),
            "hall_number": int(dataset.hall_number),
        }

    helper = build_submitted_analysis_cell(
        lattice,
        positions,
        real_types,
        space_operations,
        symprec=symprec,
    )
    bz_helper_symmetry = {
        "number": helper["analysis_spacegroup_number"],
        "symbol": helper["analysis_spacegroup_symbol"],
        "point_group": helper["analysis_point_group"],
        "seekpath_bravais": helper["seekpath_bravais"],
    }
    basename = os.path.splitext(os.path.basename(structure_file))[0]
    result = {
        "analysis_cell": helper["cell"],
        "analysis_marker_type": helper["marker_type"],
        "submitted_lattice": np.array(lattice, dtype=float),
        "submitted_sites": len(elements),
        "operation_basis_label": (
            f"submitted structure '{os.path.basename(structure_file)}'"
        ),
        "physical_symmetry": physical_symmetry,
        "bz_helper_symmetry": bz_helper_symmetry,
        # Compatibility alias for internal callers during this redesign.
        "analysis_symmetry": bz_helper_symmetry,
        "summary": {
            "marker_seeds": [
                seed.tolist() for seed in helper["marker_seeds"]
            ],
            "marker_count": helper["marker_count"],
            "marker_type": helper["marker_type"],
            "marker_types": helper["marker_types"],
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
            "source_space_operations": helper[
                "source_space_operation_count"
            ],
            "compatible_space_operations": helper[
                "compatible_space_operation_count"
            ],
            "volume_original_wrt_prim": helper[
                "volume_original_wrt_prim"
            ],
        },
    }
    result["summary"]["physical_space_operations"] = len(
        physical_operations
    )
    result["summary"]["physical_operation_set_verified"] = (
        physical_operation_set_verified
    )

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
