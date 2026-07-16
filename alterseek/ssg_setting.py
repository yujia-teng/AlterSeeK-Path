"""SSG-setting (He-marker) workflow: build spin-space-group setting inputs.

Extracted from alterseek_path.py (restructuring phase 4).
"""
import os
import shutil
import numpy as np

try:
    from findspingroup import (
        find_spin_group_acc_primitive,
        find_spin_group_acc_primitive_from_data,
    )
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = True
except ImportError:  # pragma: no cover
    find_spin_group_acc_primitive = None
    find_spin_group_acc_primitive_from_data = None
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = False

from .io import (
    _group_poscar_sites, _write_poscar, _write_without_species,
    _reciprocal_from_poscar, _dedupe_frac_positions,
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


def _write_operation_file(filename, rotations, source_indices, label):
    with open(filename, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# Found {len(rotations)} inversion-extended spin-{label} point operations\n")
        f.write(f"# Original Indices: {source_indices}\n")
        for i, rotation in enumerate(rotations):
            f.write(f"Operation_{i + 1}\n")
            for row in rotation:
                f.write(" ".join(f"{value:.10g}" for value in row) + "\n")
            f.write("\n")
    return len(rotations)


def _write_magnetic_setting_operation_files(operations, spin_axis, output_dir="."):
    flip_indices = _operation_class_indices(operations, spin_axis, flip=True)
    preserve_indices = _operation_class_indices(operations, spin_axis, flip=False)
    flip_ops, flip_sources = _collect_point_ops_from_payload(operations, flip_indices)
    preserve_ops, preserve_sources = _collect_point_ops_from_payload(operations, preserve_indices)
    flip_count = _write_operation_file(
        os.path.join(output_dir, "spin_flip_operations.txt"),
        flip_ops,
        flip_sources,
        "flipping",
    )
    preserve_count = _write_operation_file(
        os.path.join(output_dir, "spin_preserve_operations.txt"),
        preserve_ops,
        preserve_sources,
        "preserving",
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

    is_mcif = structure_file.lower().endswith(".mcif")
    if is_mcif:
        result = find_spin_group_acc_primitive(structure_file)
    else:
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
    temp_dir = os.path.join(output_dir, ".alterseek_ssgstd_tmp")
    os.makedirs(temp_dir, exist_ok=True)
    real_path = os.path.join(temp_dir, f"{basename}_ssgprim.vasp")
    mcif_path = os.path.join(output_dir, f"{basename}_ssgprim.mcif")
    magmom_path = os.path.join(temp_dir, f"{basename}_ssgprim_MAGMOM.txt")
    helper_path = os.path.join(temp_dir, f"{basename}_ssgstd.vasp")

    symbols, counts, grouped_positions, ordered_indices = _group_poscar_sites(
        elements, positions)
    _write_poscar(
        real_path,
        f"{basename} magnetic primitive from FindSpinGroup",
        lattice,
        symbols,
        counts,
        grouped_positions,
    )

    ordered_moments = moments[ordered_indices]
    ordered_elements = [elements[idx] for idx in ordered_indices]
    _write_magnetic_mcif(
        mcif_path,
        f"{basename}_ssgprim",
        lattice,
        ordered_elements,
        grouped_positions,
        ordered_moments,
    )
    with open(magmom_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Magnetic primitive POSCAR atom order matches:\n")
        f.write(f"# {real_path}\n")
        f.write("# Vector moments from FindSpinGroup acc_primitive_cell_detail:\n")
        for moment in ordered_moments:
            f.write(f"{moment[0]: .10f} {moment[1]: .10f} {moment[2]: .10f}\n")
        axis = None
        for moment in ordered_moments:
            norm = np.linalg.norm(moment)
            if norm > 1e-10:
                axis = moment / norm
                break
        if axis is not None:
            scalars = [float(np.dot(moment, axis)) for moment in ordered_moments]
            f.write("# Collinear scalar MAGMOM along first nonzero moment axis:\n")
            f.write("MAGMOM = " + " ".join(f"{value:.8g}" for value in scalars) + "\n")

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
    flip_count, preserve_count = _write_magnetic_setting_operation_files(
        operations, spin_axis, output_dir=output_dir
    )

    return {
        "real_poscar_path": real_path,
        "helper_path": helper_path,
        "mcif_path": mcif_path,
        "magmom_path": magmom_path,
        "temp_dir": temp_dir,
        "basename": basename,
        "seekpath_type_numbers": None,
        "spin_flip_operations": flip_count,
        "spin_preserve_operations": preserve_count,
        "summary": {
            "index": result.get("index"),
            "acc_symbol": result.get("acc_symbol"),
            "setting": result.get("acc_primitive_cell_setting"),
            "marker_seed": marker_helper["seed"].tolist(),
            "marker_count": len(marker_helper["markers"]),
            "marker_min_distance": marker_helper["min_distance"],
        },
    }


def finalize_magnetic_setting_outputs(
    mag_setting,
    centroid_result,
    output_dir=".",
    verbose_output=False,
):
    helper_source = centroid_result.get("standardized_structure_path")
    if not helper_source or not os.path.exists(helper_source):
        return {}

    basename = mag_setting["basename"]
    real_final = os.path.join(
        output_dir, f"{basename}_ssgstd.vasp"
    )

    _write_without_species(
        helper_source,
        real_final,
        {"He"},
        f"{basename} magnetic setting standardized",
    )
    helper_final = None
    if verbose_output:
        helper_final = os.path.join(
            output_dir, f"{basename}_ssgstd_helper.vasp"
        )
        if os.path.abspath(helper_source) != os.path.abspath(helper_final):
            shutil.copyfile(helper_source, helper_final)

    # Remove or relocate low-level seekpath artifacts for the hidden marker helper.
    # The clean final standardized structure above is the user-facing record.
    temp_dir = mag_setting.get("temp_dir")
    for path in (
        helper_source,
        centroid_result.get("standard_mapping_path"),
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
        "b_matrix_output": _reciprocal_from_poscar(real_final),
    }
    if helper_final:
        result["standard_with_helper_path"] = helper_final
    if verbose_output:
        result["intermediate_dir"] = mag_setting.get("temp_dir")
    return result
