"""Complete submitted-cell analysis for AlterSeeK-Path 2D mode."""

import warnings

import numpy as np
import spglib
from pymatgen.core import Structure

from ..geometry import get_symmetry_operations
from ..symmetry import laue_group_from_point_group, no_altermagnetism_reason
from .geometry import analyze_lattice, build_bz, to_input_fractional
from .lattice_kpoints import build_ibz, build_path
from .symmetry import project_point_operations


DEFAULT_SYMPREC = 1e-3


def run(
    filename,
    output_dir=None,
    show_plot=True,
    defer_show=False,
    verbose=True,
    seekpath_type_numbers=None,
    mode_2d=True,
    input_vacuum_axis=2,
    view_elev=None,
    view_azim=None,
    symprec=None,
    figure_basename=None,
    save_pdf=False,
    analysis_cell=None,
    analysis_marker_type=None,
):
    """Analyze one two-periodic submitted calculation cell."""
    del (
        output_dir, show_plot, defer_show, view_elev, view_azim,
        figure_basename, save_pdf, analysis_marker_type,
    )
    if not mode_2d:
        raise ValueError("mode2d.compute_centroid.run requires mode_2d=True")
    if analysis_cell is None:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="We strongly encourage explicit.*encoding"
            )
            structure = Structure.from_file(filename)
        lattice = np.asarray(structure.lattice.matrix, dtype=float)
        positions = np.asarray(structure.frac_coords, dtype=float)
        numbers = (
            [site.Z for site in structure.species]
            if seekpath_type_numbers is None
            else [int(value) for value in seekpath_type_numbers]
        )
    else:
        lattice = np.asarray(analysis_cell[0], dtype=float)
        positions = np.asarray(analysis_cell[1], dtype=float)
        numbers = [int(value) for value in analysis_cell[2]]
        if seekpath_type_numbers is not None:
            raise ValueError(
                "seekpath_type_numbers cannot be combined with analysis_cell"
            )
    if len(numbers) != len(positions):
        raise ValueError("Analysis-cell type numbers must match its sites")

    tolerance = DEFAULT_SYMPREC if symprec is None else float(symprec)
    dataset = spglib.get_symmetry_dataset(
        (lattice, positions, numbers), symprec=tolerance
    )
    if dataset is None:
        raise RuntimeError("Could not determine submitted-cell symmetry")
    layer_dataset = spglib.get_symmetry_layerdataset(
        (lattice, positions, numbers),
        aperiodic_dir=input_vacuum_axis,
        symprec=tolerance,
    )
    if layer_dataset is None:
        raise RuntimeError("Could not determine the submitted-cell layer group")
    shortest = max(float(np.min(np.linalg.norm(lattice, axis=1))), 1e-12)
    lattice_2d = analyze_lattice(
        lattice,
        input_vacuum_axis,
        orthogonality_tol=max(2e-3, tolerance / shortest),
        metric_tol=max(2e-3, 2.0 * tolerance / shortest),
    )
    path_data = build_path(lattice_2d, dataset.pointgroup)
    point_operations = project_point_operations(
        lattice_2d, dataset.rotations, add_inversion=True
    )
    bz_polygon = build_bz(lattice_2d.reciprocal_2d)
    ibz_polygon, centroid_2d, ibz_area, ibz_labels = build_ibz(
        lattice_2d,
        path_data,
        dataset.pointgroup,
        len(point_operations),
    )
    centroid_frac = to_input_fractional(centroid_2d, lattice_2d)[0]
    ibz_polygon_frac = to_input_fractional(ibz_polygon, lattice_2d).tolist()
    centroid_cart = centroid_frac @ lattice_2d.reciprocal_3d
    path_points = dict(path_data["points"])
    path = list(path_data["path"])
    path_points_cart = {
        label: np.asarray(coords) @ lattice_2d.reciprocal_3d
        for label, coords in path_points.items()
    }
    sym_ops_cart, unique_ops = get_symmetry_operations(
        lattice_2d.reciprocal_3d, dataset
    )
    lattice_name = lattice_2d.lattice_class.replace("_", "-")
    laue_group = laue_group_from_point_group(dataset.pointgroup)
    no_altermag = no_altermagnetism_reason(dataset.pointgroup, dataset.number)
    if verbose:
        branch = (
            f", {lattice_2d.centered_branch}"
            if lattice_2d.centered_branch else ""
        )
        print(
            f"[2D mode] submitted lattice: {lattice_name}{branch}; "
            f"projected point operations: {len(point_operations)}"
        )
        print(
            "[2D mode] centroid (submitted reciprocal basis): "
            f"[{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, "
            f"{centroid_frac[2]:.6f}]"
        )

    return {
        "sc_type": lattice_name,
        "lattice_key": lattice_name,
        "seekpath_bravais": None,
        "spacegroup": int(dataset.number),
        "sg_symbol": str(dataset.international),
        "point_group": str(dataset.pointgroup),
        "layer_group_number": int(layer_dataset.number),
        "layer_group_symbol": str(layer_dataset.international),
        "layer_point_group": str(layer_dataset.pointgroup),
        "laue_group": laue_group,
        "no_altermagnetism": no_altermag,
        "kpoints_frac": path_points,
        "centroid_cart": centroid_cart,
        "centroid_frac": centroid_frac,
        "ibz_volume": float(ibz_area),
        "n_symmetry_ops": len(unique_ops),
        "sp_path": path,
        "sp_point_coords": path_points,
        "b_matrix": lattice_2d.reciprocal_3d,
        "bz_loops": None,
        "bz_center": None,
        "bz_span": None,
        "elev": 90.0,
        "azim": 0.0,
        "ibz_kpoints_frac": path_points,
        "path_kpoints_frac": path_points,
        "ibz_kpath": path,
        "band_kpoints_frac": path_points,
        "band_kpath": path,
        "extra_general_vertices": list(path_data["extra_vertices"]),
        "butterfly_kpath": path_data["butterfly_path"],
        "butterfly_extra_vertices": path_data["butterfly_extra_vertices"],
        "hull_pts": ibz_polygon @ lattice_2d.plane_basis.T,
        "hull_simplices": None,
        "hull_labels": ibz_labels,
        "sym_ops_cart": sym_ops_cart,
        "unique_ops": unique_ops,
        "b_matrix_conv": np.linalg.inv(lattice).T,
        "b_matrix_input": lattice_2d.reciprocal_3d,
        "mode_2d": True,
        "vacuum_axis": input_vacuum_axis,
        "ibz_polygon_frac": ibz_polygon_frac,
        "ibz_polygon_labels": ibz_labels,
        "seekpath_rotation_matrix": np.eye(3),
        "standardized_structure_path": None,
        "standard_mapping_path": None,
        "symprec": tolerance,
        "mcif_parent_recovery": None,
        "path_source_2d": True,
        "lattice_class_2d": lattice_2d.lattice_class,
        "metric_branch_2d": lattice_2d.centered_branch,
        "bz_polygon_2d": bz_polygon,
        "projected_point_operations_2d": point_operations,
        "display_figures": [],
        "path_points_cart_2d": path_points_cart,
    }
