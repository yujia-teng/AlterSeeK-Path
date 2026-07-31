"""SSG-setting (He-marker) workflow: build spin-space-group setting inputs.

Extracted from alterseek_path.py (restructuring phase 4).
"""
import os
import shutil
import numpy as np

try:
    from findspingroup import (
        find_spin_group_acc_primitive_from_data,
    )
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = True
except ImportError:  # pragma: no cover
    find_spin_group_acc_primitive_from_data = None
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = False

from .io import (
    _atomic_write_text, _group_poscar_sites, _write_poscar, _write_without_species,
    _dedupe_frac_positions,
    _min_periodic_cart_distance, _write_magnetic_mcif,
    _load_magnetic_input_data,
)


_MARKER_SEEDS = [
    np.array([0.11000000, 0.12000000, 0.15000001]),
    np.array([0.13000000, 0.17000000, 0.23000001]),
    np.array([0.07100000, 0.19300000, 0.31700001]),
    np.array([0.21100000, 0.13700000, 0.29300001]),
]


def _marker_orbit(seed, operations):
    positions = []
    for op in operations:
        rotation = np.array(op["real_rotation"], dtype=float)
        translation = np.array(op.get("translation", [0.0, 0.0, 0.0]), dtype=float)
        positions.append(seed @ rotation.T + translation)
    return _dedupe_frac_positions(positions)


def _build_marker_helper(lattice, symbols, counts, positions, operations):
    best = None
    for seed in _MARKER_SEEDS:
        markers = _marker_orbit(seed, operations)
        helper_positions = list(positions) + markers
        candidate = {
            "seed": seed,
            "symbols": list(symbols) + ["He"],
            "counts": list(counts) + [len(markers)],
            "positions": helper_positions,
            "markers": markers,
            "min_distance": _min_periodic_cart_distance(helper_positions, lattice),
        }
        if best is None or candidate["min_distance"] > best["min_distance"]:
            best = candidate
    if best is None:
        raise RuntimeError("Could not build SSG-setting marker helper.")
    return best


def _spin_axis_from_moments(moments):
    for moment in np.asarray(moments, dtype=float):
        norm = np.linalg.norm(moment)
        if norm > 1e-10:
            return moment / norm
    raise ValueError("No nonzero magnetic moment found.")


def _operation_class_indices(operations, spin_axis, flip):
    indices = []
    for i, op in enumerate(operations):
        spin_rotation = np.array(op["spin_rotation"], dtype=float)
        mapped_axis = spin_rotation @ spin_axis
        is_flip = np.allclose(mapped_axis, -spin_axis, atol=1e-7)
        is_preserve = np.allclose(mapped_axis, spin_axis, atol=1e-7)
        if not is_flip and not is_preserve:
            continue
        if is_flip == flip:
            indices.append(i)
    return indices


def _collect_point_ops_from_payload(operations, indices, include_inversion=True):
    point_ops = []
    source_indices = []
    for idx in indices:
        rotation = np.array(operations[idx]["real_rotation"], dtype=float)
        variants = [rotation]
        if include_inversion:
            variants.append(-rotation)
        for candidate in variants:
            rounded = np.rint(candidate).astype(int)
            compare = rounded if np.allclose(candidate, rounded, atol=1e-7) else candidate
            if not any(np.allclose(compare, existing, atol=1e-7) for existing in point_ops):
                point_ops.append(compare)
                source_indices.append(int(operations[idx].get("index", idx + 1)))
    return point_ops, source_indices


def _write_operation_file(filename, rotations, source_indices, label, basis_label):
    with open(filename, "w", encoding="utf-8", newline="\n") as f:
        f.write(
            f"# Basis: {basis_label} real-space fractional basis "
            "(a1, a2, a3).\n"
        )
        f.write(
            "# k mapping: k' = R^(-T) k (mod G) in the corresponding "
            "reciprocal basis (b1, b2, b3).\n"
        )
        f.write(f"# Found {len(rotations)} inversion-extended spin-{label} point operations\n")
        f.write(f"# Original Indices: {source_indices}\n")
        for i, rotation in enumerate(rotations):
            f.write(f"Operation_{i + 1}\n")
            for row in rotation:
                f.write(" ".join(f"{value:.10g}" for value in row) + "\n")
            f.write("\n")
    return len(rotations)


def _write_magnetic_setting_operation_files(
    operations,
    spin_axis,
    basis_label,
    output_dir=".",
):
    flip_indices = _operation_class_indices(operations, spin_axis, flip=True)
    preserve_indices = _operation_class_indices(operations, spin_axis, flip=False)
    flip_ops, flip_sources = _collect_point_ops_from_payload(operations, flip_indices)
    preserve_ops, preserve_sources = _collect_point_ops_from_payload(operations, preserve_indices)
    flip_count = _write_operation_file(
        os.path.join(output_dir, "spin_flip_operations.txt"),
        flip_ops,
        flip_sources,
        "flipping",
        basis_label,
    )
    preserve_count = _write_operation_file(
        os.path.join(output_dir, "spin_preserve_operations.txt"),
        preserve_ops,
        preserve_sources,
        "preserving",
        basis_label,
    )
    return flip_count, preserve_count


def _magnetic_primitive_ssg_operations(result):
    views = (
        result.get("operation_views", {})
        .get("magnetic_primitive_cartesian", {})
        .get("views", {})
    )
    nssg_view = views.get("nssg")
    if nssg_view and nssg_view.get("ops"):
        return nssg_view["ops"]
    all_view = views.get("all")
    if all_view and all_view.get("ops"):
        return all_view["ops"]
    operations = result.get("acc_primitive_ssg_ops_cartesian")
    if operations is None:
        operations = result.get("acc_primitive_ssg_operation_matrices", [])
    return operations


def prepare_magnetic_setting_files(structure_file, moments_str="", spin_axis_cart=None, output_dir="."):
    """Write real magnetic primitive POSCAR/MCIF files from FindSpinGroup."""
    if not FIND_SG_MAGNETIC_SETTING_AVAILABLE:
        raise RuntimeError("FindSpinGroup accurate primitive API is not available.")

    # Use the same pymatgen-backed MCIF loader as the rest of AlterSeeK-Path
    # before entering FindSpinGroup's public from-data ACC-primitive route.
    # Pymatgen restores near-special fractional coordinates (for example,
    # 0.33333 -> 1/3) that a direct MCIF parse can otherwise leave just outside
    # SeeK-path/spglib's default symmetry tolerance.
    lattice_in, positions_in, elements_in, moments_in, spin_setting = _load_magnetic_input_data(
        structure_file, moments_str, spin_axis_cart
    )
    result = find_spin_group_acc_primitive_from_data(
        structure_file,
        lattice_in,
        positions_in,
        elements_in,
        [1.0] * len(elements_in),
        moments_in,
        input_spin_setting=spin_setting,
    )
    cell = result["acc_primitive_cell_detail"]
    lattice = np.array(cell["lattice"], dtype=float)
    positions = [np.array(pos, dtype=float) for pos in cell["positions"]]
    elements = [str(el) for el in cell["elements"]]
    moments = np.array(cell.get("moments", np.zeros((len(elements), 3))), dtype=float)

    basename = os.path.splitext(os.path.basename(structure_file))[0]
    temp_dir = os.path.join(output_dir, ".alterseek_magnetic_tmp")
    os.makedirs(temp_dir, exist_ok=True)
    real_path = os.path.join(temp_dir, f"{basename}_magnetic_primitive.vasp")
    mcif_path = os.path.join(output_dir, f"{basename}_magnetic_primitive.mcif")
    magmom_path = os.path.join(
        temp_dir, f"{basename}_magnetic_primitive_MAGMOM.txt")
    helper_path = os.path.join(
        temp_dir, f"{basename}_magnetic_marker_input.vasp")

    symbols, counts, grouped_positions, ordered_indices = _group_poscar_sites(
        elements, positions)
    _write_poscar(
        real_path,
        f"{basename} magnetic primitive from FindSpinGroup; "
        f"moments: {basename}_magnetic_primitive_MAGMOM.txt",
        lattice,
        symbols,
        counts,
        grouped_positions,
    )

    ordered_moments = moments[ordered_indices]
    ordered_elements = [elements[idx] for idx in ordered_indices]
    _write_magnetic_mcif(
        mcif_path,
        f"{basename}_magnetic_primitive",
        lattice,
        ordered_elements,
        grouped_positions,
        ordered_moments,
    )
    axis = None
    for moment in ordered_moments:
        norm = np.linalg.norm(moment)
        if norm > 1e-10:
            axis = moment / norm
            break
    with open(magmom_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Magnetic primitive POSCAR atom order matches:\n")
        f.write(f"# {real_path}\n")
        f.write(f"# Species order: {' '.join(symbols)}\n")
        f.write(f"# Counts: {' '.join(str(count) for count in counts)}\n")
        if axis is not None:
            scalars = [float(np.dot(moment, axis)) for moment in ordered_moments]
            f.write("# Collinear spin axis (Cartesian): " + " ".join(
                f"{value:.8g}" for value in axis) + "\n")
            f.write("# Collinear scalar MAGMOM along first nonzero moment axis:\n")
            f.write("MAGMOM = " + " ".join(
                f"{value:.8g}" for value in scalars) + "\n")

    operations = _magnetic_primitive_ssg_operations(result)
    marker_helper = _build_marker_helper(
        lattice,
        symbols,
        counts,
        grouped_positions,
        operations,
    )
    _write_poscar(
        helper_path,
        f"{basename} SSG-setting marker helper from FindSpinGroup",
        lattice,
        marker_helper["symbols"],
        marker_helper["counts"],
        marker_helper["positions"],
    )
    spin_axis = _spin_axis_from_moments(moments)
    magnetic_basis = (
        f"magnetic primitive cell '{os.path.basename(mcif_path)}'"
    )
    flip_count, preserve_count = _write_magnetic_setting_operation_files(
        operations,
        spin_axis,
        magnetic_basis,
        output_dir=output_dir,
    )

    return {
        "real_poscar_path": real_path,
        "helper_path": helper_path,
        "mcif_path": mcif_path,
        "magmom_path": magmom_path,
        "temp_dir": temp_dir,
        "basename": basename,
        "seekpath_type_numbers": None,
        "magnetic_cell_sites": len(elements),
        # Both lattices are kept so finalize_magnetic_setting_outputs() can
        # decide whether the magnetic cell is genuinely a different cell or
        # merely the submitted one with its axes relabelled.
        "magnetic_primitive_lattice": lattice,
        "submitted_lattice": np.array(lattice_in, dtype=float),
        "magnetic_symbols": list(symbols),
        "magnetic_counts": list(counts),
        "spin_flip_operations": flip_count,
        "spin_preserve_operations": preserve_count,
        "operation_basis_label": magnetic_basis,
        "summary": {
            "index": result.get("index"),
            "acc_symbol": result.get("acc_symbol"),
            "setting": result.get("acc_primitive_cell_setting"),
            "marker_seed": marker_helper["seed"].tolist(),
            "marker_count": len(marker_helper["markers"]),
            "marker_min_distance": marker_helper["min_distance"],
        },
    }


def _is_axis_relabelling(magnetic_lattice, submitted_lattice, tol=1e-6):
    """True when the magnetic cell is the submitted cell with its axes merely
    reordered and/or reversed, rather than a genuinely different cell.

    FindSpinGroup returns the magnetic primitive cell in its own convention,
    which routinely permutes and flips the axes even when the cell is the one
    that was submitted. MnSe2 (MAGNDATA 1.0.47) is the clean example: the
    submitted 36-site cell comes back as the same three vectors in the order
    (c, a, b), so the transformation is a bare permutation. Re-expressing the
    path in it would rewrite every coordinate for no physical reason.

    GdAuGe AFM5 is the contrasting case: its transformation combines two
    submitted vectors, which is what turns the 120-degree hexagonal cell into
    the 90-degree orthorhombic one. That is a real change of cell and the
    calculation has to follow it.
    """
    submitted = np.asarray(submitted_lattice, dtype=float)
    magnetic = np.asarray(magnetic_lattice, dtype=float)
    try:
        transform = magnetic @ np.linalg.inv(submitted)
    except np.linalg.LinAlgError:
        return False
    if not np.allclose(transform, np.rint(transform), atol=tol):
        return False
    rounded = np.rint(transform)
    if not np.all(np.isin(rounded, (-1.0, 0.0, 1.0))):
        return False
    return (np.all(np.count_nonzero(rounded, axis=0) == 1)
            and np.all(np.count_nonzero(rounded, axis=1) == 1))


def _is_proper_magnetic_supercell(
    magnetic_lattice,
    submitted_lattice,
    tol=1e-6,
):
    """True when the submitted calculation cell is an integer supercell of
    the magnetic primitive cell.

    A deliberate calculation supercell remains a valid output basis.  The
    magnetic primitive cell is still the right cell for the G0/IBZ analysis,
    but the resulting k-points must be transformed back to the reciprocal
    basis of the cell VASP will actually use.  GdAuGe SUPERCELL_221 is the
    regression case: its submitted cell has index four over the six-site
    magnetic primitive cell.

    Index one is deliberately excluded.  A nontrivial unimodular basis change
    is not a supercell and can be physically important for making the magnetic
    symmetry explicit (GdAuGe SUPERCELL_211, hexagonal-looking input to the
    SSG-adapted orthorhombic setting).
    """
    submitted = np.asarray(submitted_lattice, dtype=float)
    magnetic = np.asarray(magnetic_lattice, dtype=float)
    try:
        transform = submitted @ np.linalg.inv(magnetic)
    except np.linalg.LinAlgError:
        return False
    rounded = np.rint(transform)
    if not np.allclose(transform, rounded, atol=tol):
        return False
    determinant = abs(float(np.linalg.det(rounded)))
    index = int(round(determinant))
    return index > 1 and np.isclose(determinant, index, atol=tol)


def finalize_magnetic_setting_outputs(
    mag_setting,
    centroid_result,
    output_dir=".",
    verbose_output=False,
    calculation_cell_dir=".",
):
    helper_source = centroid_result.get("standardized_structure_path")
    if not helper_source or not os.path.exists(helper_source):
        return {}

    basename = mag_setting["basename"]
    real_final = os.path.join(
        output_dir, f"{basename}_seekpath_standard.vasp"
    )

    _write_without_species(
        helper_source,
        real_final,
        {"He"},
        f"{basename} SeeK-path standard cell from SSG G0 analysis",
    )
    helper_final = None
    if verbose_output:
        helper_final = os.path.join(
            output_dir, f"{basename}_seekpath_marker_helper.vasp"
        )
        if os.path.abspath(helper_source) != os.path.abspath(helper_final):
            shutil.copyfile(helper_source, helper_final)

    # The basis mapping is the only record of how the path's reciprocal basis
    # relates to the submitted cell, and this route is where that is least
    # obvious, since the cell being standardized is one the user never supplied.
    # Keep it under the same name the ordinary route uses.
    mapping_source = centroid_result.get("standard_mapping_path")
    mapping_final = None
    if mapping_source and os.path.exists(mapping_source):
        mapping_final = os.path.join(output_dir, f"{basename}_seekpath_basis_mapping.txt")
        if os.path.abspath(mapping_source) != os.path.abspath(mapping_final):
            shutil.copyfile(mapping_source, mapping_final)

    # Output basis. Reaching here means the analysis had to run in the magnetic
    # cell, but that on its own does not mean the calculation cell changed:
    # FindSpinGroup routinely returns the submitted cell with its axes reordered
    # or reversed (MnSe2), and a deliberately submitted calculation supercell
    # remains valid (GdAuGe SUPERCELL_221). Only a same-volume, nontrivial basis
    # change -- GdAuGe AFM5's hexagonal-looking 211 cell to the SSG-adapted
    # orthorhombic setting -- moves the calculation off the user's own cell.
    #
    # Either way the basis is never taken from the standardized cell written
    # above: that is spglib's standardized *conventional* cell, so for a centred
    # lattice it is 2x or 4x the primitive, and taking the basis from it used to
    # hand the user a cell several times larger than the magnetism required.
    #
    # Runs before the temp-dir cleanup below, because the magnetic primitive
    # POSCAR promoted here still lives in it.
    magnetic_lattice = np.asarray(
        mag_setting["magnetic_primitive_lattice"], dtype=float)
    submitted_lattice = mag_setting.get("submitted_lattice")
    submitted_cell_is_usable = (
        submitted_lattice is not None
        and (
            _is_axis_relabelling(magnetic_lattice, submitted_lattice)
            or _is_proper_magnetic_supercell(
                magnetic_lattice, submitted_lattice)
        )
    )
    cell_changed = not submitted_cell_is_usable
    if cell_changed:
        b_matrix_output = 2 * np.pi * np.linalg.inv(magnetic_lattice).T
        output_lattice = magnetic_lattice
    else:
        output_lattice = np.asarray(submitted_lattice, dtype=float)
        b_matrix_output = 2 * np.pi * np.linalg.inv(output_lattice).T
    # Written next to KPOINTS_alter rather than into the output folder: it is
    # not a record of the run but an input the user has to calculate with, and
    # it only appears when the cell genuinely changed.
    calculation_cell_path = None
    calculation_magmom_path = None
    source = mag_setting.get("real_poscar_path") if cell_changed else None
    magmom_source = mag_setting.get("magmom_path") if cell_changed else None
    if cell_changed:
        if not source or not os.path.exists(source):
            raise RuntimeError(
                "Magnetic calculation POSCAR is missing; refusing an incomplete "
                "calculation-cell handoff."
            )
        if not magmom_source or not os.path.exists(magmom_source):
            raise RuntimeError(
                "Matching magnetic moments are missing; refusing to hand over a "
                "reordered calculation POSCAR without its MAGMOM data."
            )
        os.makedirs(calculation_cell_dir, exist_ok=True)
        calculation_cell_path = os.path.normpath(
            os.path.join(
                calculation_cell_dir,
                f"{basename}_magnetic_primitive.vasp",
            )
        )
        calculation_magmom_path = os.path.normpath(
            os.path.join(
                calculation_cell_dir,
                f"{basename}_magnetic_primitive_MAGMOM.txt",
            )
        )

        with open(source, "r", encoding="utf-8") as f:
            _atomic_write_text(calculation_cell_path, f.read())

        with open(magmom_source, "r", encoding="utf-8") as f:
            magmom_lines = f.read().splitlines()
        if len(magmom_lines) >= 2 and magmom_lines[0] == \
                "# Magnetic primitive POSCAR atom order matches:":
            magmom_lines[1] = f"# {os.path.basename(calculation_cell_path)}"
        else:
            magmom_lines = [
                "# Magnetic primitive POSCAR atom order matches:",
                f"# {os.path.basename(calculation_cell_path)}",
                *magmom_lines,
            ]
        _atomic_write_text(
            calculation_magmom_path,
            "\n".join(magmom_lines) + "\n",
        )

    if mapping_final:
        if calculation_cell_path:
            output_label = os.path.basename(calculation_cell_path)
        else:
            output_label = "submitted structure (calculation cell unchanged)"
        with open(mapping_final, "r", encoding="utf-8") as f:
            mapping_text = f.read().rstrip()
        output_lines = [
            "",
            f"# kpoints_output_lattice is the direct lattice of {output_label}.",
            "# KPOINTS_alter fractional coordinates are written in its "
            "reciprocal basis.",
            "kpoints_output_lattice:",
        ]
        output_lines.extend(
            "  " + " ".join(f"{float(value): .10f}" for value in row)
            for row in output_lattice
        )
        _atomic_write_text(
            mapping_final,
            mapping_text + "\n" + "\n".join(output_lines) + "\n",
        )

    # Remove or relocate low-level seekpath artifacts for the hidden marker helper.
    # The clean final standardized structure above is the user-facing record.
    temp_dir = mag_setting.get("temp_dir")
    for path in (
        helper_source,
        mapping_source,
    ):
        if not path or not os.path.exists(path):
            continue
        try:
            if verbose_output and temp_dir and os.path.isdir(temp_dir):
                shutil.move(path, os.path.join(temp_dir, os.path.basename(path)))
            else:
                os.remove(path)
        except OSError:
            pass

    if not verbose_output and temp_dir and os.path.isdir(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except OSError:
            pass

    result = {
        "standard_real_path": real_final,
        "b_matrix_output": b_matrix_output,
        "cell_changed": cell_changed,
    }
    if calculation_cell_path:
        result["calculation_cell_path"] = calculation_cell_path
    if calculation_magmom_path:
        result["calculation_magmom_path"] = calculation_magmom_path
    if mag_setting.get("magnetic_symbols"):
        result["calculation_species_order"] = list(
            mag_setting["magnetic_symbols"])
    if mapping_final:
        result["standard_mapping_path"] = mapping_final
    if helper_final:
        result["standard_with_helper_path"] = helper_final
    if verbose_output:
        result["intermediate_dir"] = mag_setting.get("temp_dir")
    return result
