#!/usr/bin/env python3
"""
IBZ Centroid Calculator (Hybrid: seekpath + HPKOT)
================================================================
Uses seekpath for:  lattice type detection, cell standardization
Uses our own data:  curated IRBZ k-point vertices (HPKOT kP convention)

This ensures the IRBZ shape is consistent for all space groups within
the same extended Bravais lattice type while preserving paper-defined
optional path points such as H_2.

Supports all HPKOT extended Bravais lattice variations.

Usage:
    python compute_centroid_hybrid.py <structure_file>
    python compute_centroid_hybrid.py <structure_file> <output_dir>

Requires:
    pip install seekpath pymatgen spglib numpy scipy matplotlib sympy
"""

import sys
import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import seekpath
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
import spglib

from .atomic_write import _atomic_open_text
from .mcif import _MCIF_PARENT_SYMPREC_CANDIDATES, _declared_mcif_parent_hint

# See find_sf_operations._DEFAULT_SYMPREC for why 1e-3 is used rather than spglib's 1e-5 default.
# Override it per run with `symprec` in alterseek_input.toml.
_DEFAULT_SYMPREC = 1e-3


def _mC2d_centered_rectangular_path(submitted_lattice, vacuum_axis):
    """Return the canonical 2D centred-rectangular path for an mC slab."""
    in_plane = [index for index in range(3) if index != vacuum_axis]
    first, second = np.asarray(submitted_lattice, dtype=float)[in_plane]
    rectangular_axes = (first + second, first - second)
    lengths = sorted(np.linalg.norm(vector) for vector in rectangular_axes)
    if lengths[0] < 1e-10:
        raise RuntimeError("Could not determine the mC 2D in-plane metric")
    x_2d = 0.25 * (1.0 + (lengths[0] / lengths[1]) ** 2)
    def point(first_value, second_value):
        result = [0.0, 0.0, 0.0]
        result[in_plane[0]] = first_value
        result[in_plane[1]] = second_value
        return result
    points = {
        GAMMA_LABEL: point(0.0, 0.0), 'Y': point(0.5, 0.5),
        'S': point(0.0, 0.5), 'C': point(x_2d, 1.0 - x_2d),
        'SIGMA': point(-x_2d, x_2d),
    }
    return x_2d, points


def _select_mcif_parent_symprec(filename, cell, positions, numbers, fallback=None):
    """Use the smallest conservative tolerance that recovers a declared parent.

    Preferred whenever the file declares a parent, since it is validated
    against the structure's own ground truth. `fallback` (the configured
    symprec) applies when there is nothing to validate against.
    """
    if fallback is None:
        fallback = _DEFAULT_SYMPREC
    hint = _declared_mcif_parent_hint(filename)
    if hint is None or len(positions) % hint["index"] != 0:
        return fallback, None

    expected_sites = len(positions) // hint["index"]
    structure_cell = (cell, positions, numbers)
    for symprec in _MCIF_PARENT_SYMPREC_CANDIDATES:
        dataset = spglib.get_symmetry_dataset(structure_cell, symprec=symprec)
        primitive = spglib.find_primitive(structure_cell, symprec=symprec)
        if dataset is None or primitive is None or len(primitive[1]) != expected_sites:
            continue
        parent_number = hint.get("spacegroup_number")
        if parent_number is not None and int(dataset.number) != parent_number:
            continue
        recovered = dict(hint)
        recovered.update({
            "symprec": symprec,
            "input_sites": len(positions),
            "primitive_sites": expected_sites,
            "detected_spacegroup_number": int(dataset.number),
            "detected_spacegroup_symbol": str(dataset.international),
        })
        return symprec, recovered
    return fallback, None


def _write_seekpath_standard_poscar(lattice, positions, types, output_path, source_name):
    """Write the standardized input-cell setting used by SeeK-path."""
    lattice = np.array(lattice, dtype=float)
    positions = np.array(positions, dtype=float)
    types = list(types)

    species_order = []
    grouped = {}
    for pos, atomic_number in zip(positions, types):
        symbol = Element.from_Z(int(atomic_number)).symbol
        if symbol not in grouped:
            species_order.append(symbol)
            grouped[symbol] = []
        grouped[symbol].append(pos)

    lines = [
        f"SeeK-path standardized cell from {source_name}",
        "1.0",
    ]
    for row in lattice:
        lines.append("   " + " ".join(f"{x:22.16f}" for x in row))
    lines.append(" ".join(species_order))
    lines.append(" ".join(str(len(grouped[symbol])) for symbol in species_order))
    lines.append("Direct")
    for symbol in species_order:
        for pos in grouped[symbol]:
            frac = np.mod(pos, 1.0)
            lines.append("   " + " ".join(f"{x:22.16f}" for x in frac) + f" {symbol}")

    with _atomic_open_text(output_path) as f:
        f.write("\n".join(lines) + "\n")


def _write_seekpath_basis_mapping(
    input_lattice,
    primitive_lattice,
    conventional_lattice,
    rotation_matrix,
    output_path,
    source_name,
):
    """Record the submitted analysis lattice to SeeK-path basis chain."""
    def _fmt_matrix(mat):
        return "\n".join(
            "  " + " ".join(f"{float(x): .10f}" for x in row)
            for row in np.array(mat, dtype=float)
        )

    input_lattice = np.array(input_lattice, dtype=float)
    primitive_lattice = np.array(primitive_lattice, dtype=float)
    conventional_lattice = np.array(conventional_lattice, dtype=float)
    rotation_matrix = np.array(rotation_matrix, dtype=float)
    lines = [
        "# SeeK-path standardization mapping",
        f"# analysis_input_lattice is the lattice of {source_name}, "
        "the cell given to SeeK-path.",
        "# seekpath_standard_primitive_lattice is the internal HPKOT path basis.",
        "# seekpath_standard_conventional_lattice is SeeK-path's conventional "
        "diagnostic setting.",
        "# rotation_matrix is reported by SeeK-path for the input-to-standard orientation.",
        "",
        "analysis_input_lattice:",
        _fmt_matrix(input_lattice),
        "",
        "seekpath_standard_primitive_lattice:",
        _fmt_matrix(primitive_lattice),
        "",
        "seekpath_standard_conventional_lattice:",
        _fmt_matrix(conventional_lattice),
        "",
        "seekpath_rotation_matrix:",
        _fmt_matrix(rotation_matrix),
        "",
        "# kpoints_output_lattice is the submitted direct lattice. Final VASP "
        "and QE fractional coordinates use its reciprocal basis.",
        "kpoints_output_lattice:",
        _fmt_matrix(input_lattice),
        "",
    ]
    with _atomic_open_text(output_path) as f:
        f.write("\n".join(lines))


from .lattice_kpoints import (
    get_kpoints, get_hull_kpoints, get_hull_kpath,
    get_display_labels, get_params, _normalize_label,
)
from .symmetry import (
    laue_group_from_point_group,
    no_altermagnetism_reason,
    seekpath_to_hpkot_type,
)
from .geometry import (
    get_symmetry_operations,
    calculate_volume_centroid,
    detect_vacuum_axis_2d,
    trace_vacuum_axis_2d,
    area_centroid_2d,
    ordered_2d_polygon_frac,
    check_input_slab,
    compute_symbolic_centroid,
    get_bz_loops,
    build_symmetry_ibz_cell,
    triclinic_halfspace_normal,
    triclinic_half_bz_cell,
)
from .plotting_common import (
    GAMMA_LABEL,
    _save_figure,
    _print_saved_paths,
    alterseek_plot_style,
)
from .plotting_3d import setup_3d_ax, plot_ibz


def _analyze_kspace(
    filename,
    *,
    analysis_cell=None,
    seekpath_type_numbers,
    mode_2d,
    input_vacuum_axis,
    symprec,
    verbose,
):
    """Perform the required structure, symmetry, IBZ, and path analysis."""
    if analysis_cell is None:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="We strongly encourage explicit.*encoding",
            )
            struct = Structure.from_file(filename)
        a_matrix = np.asarray(struct.lattice.matrix, dtype=float)
        cell = a_matrix.tolist()
        positions = struct.frac_coords.tolist()
        if seekpath_type_numbers is None:
            numbers = [site.Z for site in struct.species]
        else:
            numbers = [int(number) for number in seekpath_type_numbers]
    else:
        a_matrix = np.asarray(analysis_cell[0], dtype=float)
        cell = a_matrix.tolist()
        positions = np.asarray(analysis_cell[1], dtype=float).tolist()
        numbers = [int(number) for number in analysis_cell[2]]
        if seekpath_type_numbers is not None:
            raise ValueError(
                "seekpath_type_numbers cannot be combined with analysis_cell."
            )

    if len(numbers) != len(positions):
        raise ValueError(
            "Analysis-cell type numbers must match the number of sites."
        )
    requested_symprec = _DEFAULT_SYMPREC if symprec is None else float(symprec)
    if analysis_cell is None:
        symprec, mcif_parent_recovery = _select_mcif_parent_symprec(
            filename, cell, positions, numbers, fallback=requested_symprec
        )
    else:
        symprec = requested_symprec
        mcif_parent_recovery = None

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(r".*dict interface is deprecated.*"
                     r"Use attribute interface instead.*"),
        )
        warnings.filterwarnings(
            "ignore",
            message=r".*dict interface is deprecated.*",
        )
        warnings.filterwarnings(
            "ignore",
            category=DeprecationWarning,
            module=r"seekpath\.hpkot(\..*)?",
        )
        sp_result = seekpath.get_path(
            (cell, positions, numbers),
            with_time_reversal=True,
            symprec=symprec,
        )

    input_dataset = spglib.get_symmetry_dataset(
        (cell, positions, numbers),
        symprec=symprec,
    )
    spg_cell = (
        np.array(sp_result['primitive_lattice']),
        np.array(sp_result['primitive_positions']),
        sp_result['primitive_types'],
    )
    dataset = spglib.get_symmetry_dataset(spg_cell, symprec=symprec)
    b_matrix = np.array(sp_result['reciprocal_primitive_lattice'])
    b_matrix_input = 2 * np.pi * np.linalg.inv(np.array(a_matrix)).T
    conventional_lattice = np.array(
        sp_result.get('conv_lattice', sp_result['primitive_lattice'])
    )
    b_matrix_conv = np.linalg.inv(conventional_lattice).T
    b1, b2, b3 = b_matrix

    vacuum_axis = None
    if mode_2d:
        # The declared input axis governs the slice; detection only cross-checks it.
        detected_axis, vacuum_info = detect_vacuum_axis_2d(b_matrix)
        vacuum_axis, trace_info = trace_vacuum_axis_2d(
            a_matrix, input_vacuum_axis, b_matrix,
            sp_result['rotation_matrix'],
        )
        if not trace_info['traced']:
            vacuum_axis = detected_axis
            print("[2D mode][Warning] input vacuum axis "
                  f"'{'abc'[input_vacuum_axis]}' matches no single "
                  "standardized axis; using the longest axis "
                  f"'{'abc'[detected_axis]}' instead.")
        elif vacuum_axis != detected_axis:
            print("[2D mode][Warning] input vacuum axis "
                  f"'{'abc'[input_vacuum_axis]}' maps to standardized "
                  f"'{'abc'[vacuum_axis]}', but the longest standardized axis "
                  f"is '{'abc'[detected_axis]}'; verify the structure is a "
                  "proper slab.")
        if verbose:
            print(f"\n[2D mode] vacuum axis (standardized frame): "
                  f"{vacuum_axis} ('{'abc'[vacuum_axis]}'); reciprocal norms "
                  f"{[round(value, 4) for value in vacuum_info['reciprocal_norms']]}")
        if not (vacuum_info['separated'] and vacuum_info['orthogonal']):
            print("[2D mode][Warning] vacuum-axis detection is ambiguous "
                  f"(separated={vacuum_info['separated']}, "
                  f"orthogonal={vacuum_info['orthogonal']}); verify the "
                  "structure is a proper slab.")
        for warning_text in check_input_slab(a_matrix, input_vacuum_axis):
            print(f"[2D mode][Warning] {warning_text}")

    sg = dataset.number
    laue_group = laue_group_from_point_group(dataset.pointgroup)
    no_altermag = no_altermagnetism_reason(dataset.pointgroup, sg)
    if verbose:
        print(f"\nSpace Group: {sg} ({dataset.international})")
        print(f"Point Group: {dataset.pointgroup}")
        print(f"Laue Group: "
              f"{laue_group if laue_group is not None else 'Unknown'}")
        if no_altermag:
            print(f"[Note] {no_altermag['reason']} for Laue group "
                  f"{no_altermag['laue_group']}")
        print(f"Seekpath Bravais: {sp_result['bravais_lattice_extended']}")

    sc_type, conv_params = seekpath_to_hpkot_type(sp_result)
    sc_display = sc_type
    centroid_type = sc_type
    kpath = [
        (_normalize_label(start), _normalize_label(end))
        for start, end in sp_result['path']
    ]
    seekpath_point_coords = {
        _normalize_label(label): list(coords)
        for label, coords in sp_result['point_coords'].items()
    }
    kpoints_frac = get_kpoints(
        sc_type,
        conv_params['a'],
        conv_params.get('b'),
        conv_params.get('c'),
        conv_params.get('alpha'),
    )
    path_kpoints_frac = {
        label: seekpath_point_coords[label]
        for label in {label for segment in kpath for label in segment}
        if label in seekpath_point_coords
    }
    # In the k_z = 0 plane the hP1/hP2 wedge already tiles exactly, so the 3D
    # doubling would double-cover it. Tetragonal 75-88 still needs its doubling.
    doubling_sg = None if (mode_2d and sc_type in ('hP1', 'hP2')) else sg
    kpoints_frac_centroid = get_hull_kpoints(
        sc_type,
        conv_params['a'],
        conv_params.get('b'),
        conv_params.get('c'),
        conv_params.get('alpha'),
        spacegroup_number=doubling_sg,
    )
    hull_kpath = get_hull_kpath(sc_type, spacegroup_number=doubling_sg)
    display_labels = get_display_labels(sc_type)
    params = get_params(
        sc_type,
        conv_params['a'],
        conv_params.get('b'),
        conv_params.get('c'),
        conv_params.get('alpha'),
    )
    if params and verbose:
        print(f"Parameters: "
              f"{', '.join(f'{key}={value:.6f}' for key, value in params.items())}")

    # A slab with an mC HPKOT setting has the same physical centred-
    # rectangular BZ as oC/oA, but its 3D mC path does not contain the two
    # metric-dependent 2D vertices.  Define the 2D-only table in the native
    # mC primitive reciprocal basis; never change the 3D HPKOT table itself.
    if mode_2d and sc_type in {'mC1', 'mC2', 'mC3'}:
        x_2d, mC_2d_points = _mC2d_centered_rectangular_path(
            a_matrix, input_vacuum_axis
        )
        kpath = [
            (GAMMA_LABEL, 'Y'), ('Y', 'C'),
            ('SIGMA', GAMMA_LABEL), (GAMMA_LABEL, 'S'),
        ]
        path_kpoints_frac = dict(mC_2d_points)
        kpoints_frac_centroid = dict(mC_2d_points)
        display_labels.update({'C': 'C', 'SIGMA': r'$\\Sigma$'})
        if verbose:
            print(f"[2D mode] mC centred-rectangular parameter: X_2D={x_2d:.6f}")

    path_kpoints_cart = {
        key: value[0] * b1 + value[1] * b2 + value[2] * b3
        for key, value in path_kpoints_frac.items()
    }
    kpoints_cart_centroid = {
        key: value[0] * b1 + value[1] * b2 + value[2] * b3
        for key, value in kpoints_frac_centroid.items()
    }
    if sc_type == 'mP1':
        for label in ('Y', 'C'):
            kpoints_frac_centroid.pop(label, None)
            kpoints_cart_centroid.pop(label, None)

    sym_ops_cart, unique_ops = get_symmetry_operations(b_matrix, dataset)
    if verbose:
        print(f"\nSymmetry operations: {len(sym_ops_cart)}")
        print(f"With time-reversal: {len(unique_ops)}")

    doubled_ibz_case = (
        (75 <= sg <= 88 and sc_type in {'tP1', 'tI1', 'tI2'})
        or (149 <= sg <= 176 and sc_type == 'hP2')
        or (sc_type == 'hP1' and sg in {149, 151, 153, 157, 159, 162, 163})
    )
    band_kpath = list(hull_kpath if doubled_ibz_case else kpath)
    band_kpoints_frac = dict(path_kpoints_frac)
    extra_general_vertices = []
    if sc_type == 'hP1' and sg in {149, 151, 153, 157, 159, 162, 163}:
        extra_general_vertices = ["K_A", "H_A"]
        for label in extra_general_vertices:
            if label in kpoints_frac_centroid:
                band_kpoints_frac[label] = kpoints_frac_centroid[label]
    elif sc_type == 'hP2' and 149 <= sg <= 176:
        extra_general_vertices = ["L_A", "M_A"]
        for label in ("L_A", "M_A"):
            if label in kpoints_frac_centroid:
                band_kpoints_frac[label] = kpoints_frac_centroid[label]
    elif 75 <= sg <= 88 and sc_type in {'tP1', 'tI1', 'tI2'}:
        tetragonal_extra = {
            'tP1': ["X_A", "R_A"],
            'tI1': ["X_A", "P_A"],
            'tI2': ["R_A", "S_0A", "S_A", "N_A"],
        }
        extra_general_vertices = [
            label for label in tetragonal_extra[sc_type]
            if label in kpoints_frac_centroid
        ]
        for label in extra_general_vertices:
            band_kpoints_frac[label] = kpoints_frac_centroid[label]

    path_labels = {label for segment in band_kpath for label in segment}
    extra_general_vertices = [
        label for label in extra_general_vertices if label not in path_labels
    ]
    if mode_2d:
        def _in_plane_label(label):
            value = (
                band_kpoints_frac.get(label)
                or kpoints_frac_centroid.get(label)
                or path_kpoints_frac.get(label)
            )
            return value is not None and abs(value[vacuum_axis]) < 1e-4

        band_kpath = [
            segment for segment in band_kpath
            if all(_in_plane_label(label) for label in segment)
        ]
        extra_general_vertices = [
            label for label in extra_general_vertices
            if _in_plane_label(label)
        ]
        path_labels = {
            label for segment in band_kpath for label in segment
        }
        if doubled_ibz_case:
            def _is_doubled_label(label):
                base = label[1:] if label.startswith('_') else label
                return (
                    '_A' in base
                    or '_B' in base
                    or base.endswith('_0A')
                )

            connected_to_gamma = {
                segment[1] if segment[0] == GAMMA_LABEL else segment[0]
                for segment in band_kpath if GAMMA_LABEL in segment
            }
            doubled_labels = [
                label for label in path_labels if _is_doubled_label(label)
            ]
            for label in doubled_labels:
                if label not in connected_to_gamma:
                    band_kpath.append((GAMMA_LABEL, label))
            path_labels = {
                label for segment in band_kpath for label in segment
            }
        if verbose:
            print(f"[2D mode] in-plane band path: {len(band_kpath)} "
                  f"segments, labels {sorted(path_labels)}")

    kpath_plot = band_kpath
    kpoints_cart_plot = dict(kpoints_cart_centroid)
    for label in {label for segment in kpath_plot for label in segment}:
        if label in path_kpoints_cart:
            kpoints_cart_plot[label] = path_kpoints_cart[label]
    kpoints_frac_for_output = dict(path_kpoints_frac)
    for label in {label for segment in kpath_plot for label in segment}:
        if label in kpoints_frac_centroid:
            kpoints_frac_for_output[label] = kpoints_frac_centroid[label]

    return {
        'a_matrix': a_matrix,
        'input_dataset': input_dataset,
        'sp_result': sp_result,
        'dataset': dataset,
        'b_matrix': b_matrix,
        'b_matrix_input': b_matrix_input,
        'b_matrix_conv': b_matrix_conv,
        'vacuum_axis': vacuum_axis,
        'sg': sg,
        'laue_group': laue_group,
        'no_altermagnetism': no_altermag,
        'sc_type': sc_type,
        'sc_display': sc_display,
        'centroid_type': centroid_type,
        'conv_params': conv_params,
        'kpoints_frac': kpoints_frac,
        'kpoints_frac_centroid': kpoints_frac_centroid,
        'kpoints_cart_centroid': kpoints_cart_centroid,
        'sym_ops_cart': sym_ops_cart,
        'unique_ops': unique_ops,
        'band_kpath': band_kpath,
        'band_kpoints_frac': band_kpoints_frac,
        'extra_general_vertices': extra_general_vertices,
        'kpath_plot': kpath_plot,
        'display_labels_plot': display_labels,
        'kpoints_cart_plot': kpoints_cart_plot,
        'kpoints_frac_for_output': kpoints_frac_for_output,
        'symprec': symprec,
        'mcif_parent_recovery': mcif_parent_recovery,
    }


def _write_optional_diagnostics(
    analysis,
    centroid,
    *,
    filename,
    output_dir,
    basename,
    verbose,
    spin_log_current_run,
    analysis_marker_type,
):
    """Write optional standardization and symbolic-centroid diagnostics."""
    standardized_structure_output = os.path.join(
        output_dir, f"{basename}_seekpath_standard.vasp"
    )
    standard_mapping_output = os.path.join(
        output_dir, f"{basename}_seekpath_basis_mapping.txt"
    )
    input_dataset = analysis['input_dataset']
    sp_result = analysis['sp_result']

    standardized_structure_path = None
    if analysis_marker_type is None:
        try:
            os.makedirs(output_dir, exist_ok=True)
            _write_seekpath_standard_poscar(
                np.array(input_dataset.std_lattice),
                np.array(input_dataset.std_positions),
                list(input_dataset.std_types),
                standardized_structure_output,
                os.path.basename(filename),
            )
            standardized_structure_path = standardized_structure_output
            if verbose:
                print(f"Saved standardized structure: {standardized_structure_path}")
        except Exception as exc:
            print(f"[Warning] Could not write SeeK-path standardized structure: {exc}")

    standard_mapping_path = None
    try:
        os.makedirs(output_dir, exist_ok=True)
        _write_seekpath_basis_mapping(
            analysis['a_matrix'],
            np.array(sp_result['primitive_lattice']),
            np.array(input_dataset.std_lattice),
            np.array(sp_result['rotation_matrix']),
            standard_mapping_output,
            os.path.basename(filename),
        )
        standard_mapping_path = standard_mapping_output
        if verbose:
            print(f"Saved SeeK-path basis mapping: {standard_mapping_path}")
    except Exception as exc:
        print(f"[Warning] Could not write SeeK-path basis mapping: {exc}")

    if spin_log_current_run and centroid['hull_matches_labels']:
        try:
            sym_centroid, _ = compute_symbolic_centroid(
                analysis['kpoints_frac_centroid'],
                centroid['hull'],
                centroid['labels_list'],
                analysis['centroid_type'],
                analysis['conv_params'],
            )
            if sym_centroid is not None:
                sym_lines = "\n".join(
                    f"  {axis_name} = {sym_centroid[index]}"
                    for index, axis_name in enumerate(['k1', 'k2', 'k3'])
                )
                try:
                    with open(
                        os.path.join(output_dir, "spin_operations.txt"),
                        "a",
                        encoding="utf-8",
                        newline="\n",
                    ) as stream:
                        stream.write(
                            "\nSymbolic IBZ centroid (fractional):\n"
                            f"{sym_lines}\n"
                        )
                except OSError as exc:
                    print("[Warning] Could not record the symbolic IBZ "
                          f"centroid: {exc}")
        except Exception as exc:
            print(f"[Warning] Symbolic IBZ centroid unavailable: {exc}")

    return {
        'standardized_structure_path': standardized_structure_path,
        'standard_mapping_path': standard_mapping_path,
    }


def _format_plane(normal_frac):
    """Render an integer fractional normal as a readable plane equation."""
    terms = []
    for index, coefficient in enumerate(normal_frac):
        if not coefficient:
            continue
        axis = f"k{index + 1}"
        magnitude = abs(coefficient)
        term = axis if magnitude == 1 else f"{magnitude}*{axis}"
        sign = "-" if coefficient < 0 else ("+" if terms else "")
        terms.append(f"{sign}{term}")
    return f"{''.join(terms)}=0"


def _compute_ibz_centroid(
    *,
    mode_2d,
    vacuum_axis,
    sg,
    sc_type,
    b_matrix,
    unique_ops,
    kpoints_frac_centroid,
    kpoints_cart_centroid,
    verbose,
):
    """Compute the required numerical centroid from analyzed IBZ geometry."""
    labels_list = list(kpoints_cart_centroid.keys())
    points_arr = np.array([kpoints_cart_centroid[k] for k in labels_list])
    ibz_polygon_frac = None
    ibz_polygon_labels = None
    hull_matches_labels = False

    if mode_2d:
        # In 2D slab mode, the physical IBZ is the k[vacuum_axis] = 0 cross-section.
        # Restrict the curated hull points to that plane and take the 2D area centroid because the 3D volume centroid is meaningless and ConvexHull crashes on coplanar input.
        in_plane_labels = [
            lab for lab in labels_list
            if abs(kpoints_frac_centroid[lab][vacuum_axis]) < 1e-4
        ]
        # An HPKOT <base>_2 point that is a distinct 3D hull vertex can become a redundant mirror-image duplicate after slicing to this 2D plane.
        # This reduction is 2D-only.
        base_labels = {lab for lab in in_plane_labels if not lab.endswith('_2')}
        in_plane_labels = [
            lab for lab in in_plane_labels
            if not (lab.endswith('_2') and lab[:-2] in base_labels)
        ]
        frac_pts = np.array([
            kpoints_frac_centroid[lab] for lab in in_plane_labels
        ])
        hull = None
        centroid_frac, centroid_cart, ibz_volume = area_centroid_2d(
            frac_pts, vacuum_axis, b_matrix
        )
        ibz_polygon_frac, ibz_polygon_labels = ordered_2d_polygon_frac(
            frac_pts, vacuum_axis, labels=in_plane_labels
        )
        if verbose:
            print(f"[2D mode] in-plane hull labels: {len(in_plane_labels)} / "
                  f"{len(labels_list)} ({', '.join(in_plane_labels)})")
            print(f"[2D mode] area centroid (frac): "
                  f"[{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, "
                  f"{centroid_frac[2]:.6f}]  area={ibz_volume:.6f}")
    elif sg in (1, 2):
        # The curated triclinic points are not a closed IBZ hull, so use the Wigner-Seitz BZ instead.
        # Laue group -1 applies to both P1 and P-1, making the nonmagnetic IBZ half the BZ; the selected axis-containing plane fixes which half.
        hull = None
        normal_frac = triclinic_halfspace_normal(sc_type)
        if normal_frac is not None:
            try:
                points_arr, hull = triclinic_half_bz_cell(b_matrix, normal_frac)
                labels_list = []
                centroid_cart, ibz_volume = calculate_volume_centroid(hull)
            except Exception as exc:
                hull = None
                print("[Warning] Could not build the triclinic half-BZ "
                      f"({exc}); using the mean of the curated k-points.")
        if hull is None:
            centroid_cart = np.mean(points_arr, axis=0)
            ibz_volume = 0.0
        centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)
        if verbose:
            if hull is not None:
                print(f"\n[Note] Triclinic: IBZ = half BZ cut by "
                      f"{_format_plane(normal_frac)} (Laue -1)")
                print(f"IBZ centroid: [{centroid_frac[0]:.6f}, "
                      f"{centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]  "
                      f"volume={ibz_volume:.6f}")
            else:
                print("\n[Note] Triclinic: IBZ shading skipped")
                print(f"Centroid (mean of k-points): [{centroid_frac[0]:.6f}, "
                      f"{centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]")
    else:
        hull = ConvexHull(points_arr)
        centroid_cart, ibz_volume = calculate_volume_centroid(hull)
        centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)
        hull_matches_labels = True

        # The mC2/mC3 HPKOT boundary-label hull is slightly larger than the true C2/m fundamental domain.
        if sc_type in {'mC2', 'mC3'}:
            mono_pts, mono_simplices = build_symmetry_ibz_cell(
                b_matrix, unique_ops, centroid_cart
            )
            if mono_pts is not None and mono_simplices is not None:
                points_arr = np.array(mono_pts, dtype=float)
                hull = ConvexHull(points_arr)
                centroid_cart, ibz_volume = calculate_volume_centroid(hull)
                centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)
                hull_matches_labels = False
                if verbose:
                    print("[Note] Using symmetry/Voronoi IBZ cell for "
                          "monoclinic hull.")

    return {
        'labels_list': labels_list,
        'points_arr': points_arr,
        'hull': hull,
        'centroid_cart': centroid_cart,
        'centroid_frac': centroid_frac,
        'ibz_volume': ibz_volume,
        'ibz_polygon_frac': ibz_polygon_frac,
        'ibz_polygon_labels': ibz_polygon_labels,
        'hull_matches_labels': hull_matches_labels,
    }


@alterseek_plot_style
def _generate_figure1(
    analysis,
    centroid,
    *,
    output_dir,
    fig_basename,
    show_plot,
    defer_show,
    mode_2d,
    view_elev,
    view_azim,
    save_pdf,
    verbose,
):
    """Generate optional BZ geometry and the ordinary 3D Figure 1."""
    sg = analysis['sg']
    sc_type = analysis['sc_type']
    sc_display = analysis['sc_display']
    b_matrix = analysis['b_matrix']
    kpath_plot = analysis['kpath_plot']
    display_labels_plot = analysis['display_labels_plot']
    kpoints_cart_plot = analysis['kpoints_cart_plot']
    hull = centroid['hull']
    centroid_cart = centroid['centroid_cart']
    points_arr = centroid['points_arr']
    labels_list = centroid['labels_list']

    fig1_title = (
        f"BZ: {fig_basename} ({sc_display})"
        if sg in (1, 2) and hull is None
        else f"IBZ + BZ: {fig_basename} ({sc_display})"
    )
    display_figures = []
    default_elev, default_azim = (
        (view_elev, view_azim)
        if view_elev is not None and view_azim is not None
        else (14, 20)
    )
    elev1, azim1 = default_elev, default_azim
    fig1_path = os.path.join(
        output_dir, f'{fig_basename}_ibz_{sc_display}.png'
    )
    bz_loops = None
    bz_center = None
    bz_span = None
    fig1 = None
    fig1s = None
    try:
        # The BZ wireframe is presentation data used only by Figures 1-4.
        bz_loops = get_bz_loops(b_matrix)
        all_bz_pts = np.vstack(bz_loops)
        bz_center = np.mean(all_bz_pts, axis=0)
        bz_span = np.max(all_bz_pts) - np.min(all_bz_pts)

        if show_plot and not mode_2d:
            fig1, ax1 = setup_3d_ax(
                fig1_title,
                bz_loops,
                b_matrix,
                bz_center,
                bz_span,
                elev=default_elev,
                azim=default_azim,
            )
            plot_ibz(
                ax1,
                kpoints_cart_plot,
                kpath_plot,
                display_labels_plot,
                hull,
                centroid_cart,
                hull_pts=points_arr,
                lattice_type=sc_type,
                hull_labels=labels_list,
            )
            plt.tight_layout()
            if defer_show:
                def _save_fig1_after_show(fig=fig1, ax=ax1):
                    save_figure = None
                    try:
                        save_figure, save_ax = setup_3d_ax(
                            fig1_title,
                            bz_loops,
                            b_matrix,
                            bz_center,
                            bz_span,
                            elev=ax.elev,
                            azim=ax.azim,
                            dashed_back=True,
                        )
                        plot_ibz(
                            save_ax,
                            kpoints_cart_plot,
                            kpath_plot,
                            display_labels_plot,
                            hull,
                            centroid_cart,
                            hull_pts=points_arr,
                            lattice_type=sc_type,
                            hull_labels=labels_list,
                        )
                        plt.tight_layout()
                        saved_paths = _save_figure(
                            save_figure,
                            fig1_path,
                            extra_formats=("pdf",) if save_pdf else (),
                            dpi=300,
                            bbox_inches='tight',
                        )
                        _print_saved_paths(saved_paths)
                    finally:
                        if save_figure is not None:
                            plt.close(save_figure)
                        plt.close(fig)
                fig1._alterseek_save_after_show = _save_fig1_after_show
                display_figures.append(fig1)
            else:
                plt.show()
                elev1, azim1 = ax1.elev, ax1.azim

        if not (show_plot and defer_show) and not mode_2d:
            fig1s, ax1s = setup_3d_ax(
                fig1_title,
                bz_loops,
                b_matrix,
                bz_center,
                bz_span,
                elev=elev1,
                azim=azim1,
                dashed_back=True,
            )
            plot_ibz(
                ax1s,
                kpoints_cart_plot,
                kpath_plot,
                display_labels_plot,
                hull,
                centroid_cart,
                hull_pts=points_arr,
                lattice_type=sc_type,
                hull_labels=labels_list,
            )
            plt.tight_layout()
            saved_paths = _save_figure(
                fig1s,
                fig1_path,
                extra_formats=("pdf",) if save_pdf else (),
                dpi=300,
                bbox_inches='tight',
            )
            _print_saved_paths(saved_paths, verbose=verbose)
            plt.close(fig1s)
    except Exception as exc:
        display_figures.clear()
        for figure in (fig1s, fig1):
            if figure is not None:
                plt.close(figure)
        elev1, azim1 = default_elev, default_azim
        print(f"[Warning] Could not generate Figure 1: {exc}")

    return {
        'bz_loops': bz_loops,
        'bz_center': bz_center,
        'bz_span': bz_span,
        'elev': elev1,
        'azim': azim1,
        'display_figures': display_figures,
    }


def run(
    filename,
    output_dir=None,
    show_plot=True,
    defer_show=False,
    verbose=True,
    seekpath_type_numbers=None,
    mode_2d=False,
    input_vacuum_axis=2,
    view_elev=None,
    view_azim=None,
    symprec=None,
    figure_basename=None,
    save_pdf=False,
    spin_log_current_run=False,
    analysis_cell=None,
    analysis_marker_type=None,
):
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(filename))
    basename = os.path.splitext(os.path.basename(filename))[0]
    # Keep all figures named after the submitted structure.
    fig_basename = figure_basename or basename

    if verbose:
        print("=" * 60)
        print(f"Processing: {filename}")
        print("=" * 60)

    analysis_result = _analyze_kspace(
        filename,
        analysis_cell=analysis_cell,
        seekpath_type_numbers=seekpath_type_numbers,
        mode_2d=mode_2d,
        input_vacuum_axis=input_vacuum_axis,
        symprec=symprec,
        verbose=verbose,
    )

    centroid_result = _compute_ibz_centroid(
        mode_2d=mode_2d,
        vacuum_axis=analysis_result['vacuum_axis'],
        sg=analysis_result['sg'],
        sc_type=analysis_result['sc_type'],
        b_matrix=analysis_result['b_matrix'],
        unique_ops=analysis_result['unique_ops'],
        kpoints_frac_centroid=analysis_result['kpoints_frac_centroid'],
        kpoints_cart_centroid=analysis_result['kpoints_cart_centroid'],
        verbose=verbose,
    )

    diagnostic_result = _write_optional_diagnostics(
        analysis_result,
        centroid_result,
        filename=filename,
        output_dir=output_dir,
        basename=basename,
        verbose=verbose,
        spin_log_current_run=spin_log_current_run,
        analysis_marker_type=analysis_marker_type,
    )

    figure_result = _generate_figure1(
        analysis_result,
        centroid_result,
        output_dir=output_dir,
        fig_basename=fig_basename,
        show_plot=show_plot,
        defer_show=defer_show,
        mode_2d=mode_2d,
        view_elev=view_elev,
        view_azim=view_azim,
        save_pdf=save_pdf,
        verbose=verbose,
    )
    bz_loops = figure_result['bz_loops']
    bz_center = figure_result['bz_center']
    bz_span = figure_result['bz_span']
    elev1 = figure_result['elev']
    azim1 = figure_result['azim']
    display_figures = figure_result['display_figures']

    return {
        'sc_type': analysis_result['sc_display'],
        'lattice_key': analysis_result['sc_type'],
        'seekpath_bravais': analysis_result[
            'sp_result'
        ]['bravais_lattice_extended'],
        'spacegroup': analysis_result['sg'],
        'sg_symbol': analysis_result['dataset'].international,
        'point_group': analysis_result['dataset'].pointgroup,
        'laue_group': analysis_result['laue_group'],
        'no_altermagnetism': analysis_result['no_altermagnetism'],
        'kpoints_frac': analysis_result['kpoints_frac'],
        'centroid_cart': centroid_result['centroid_cart'],
        'centroid_frac': centroid_result['centroid_frac'],
        'ibz_volume': centroid_result['ibz_volume'],
        'n_symmetry_ops': len(analysis_result['unique_ops']),
        'sp_path': analysis_result['sp_result']['path'],
        'sp_point_coords': analysis_result['sp_result']['point_coords'],
        'b_matrix': analysis_result['b_matrix'],
        'bz_loops': bz_loops,
        'bz_center': bz_center,
        'bz_span': bz_span,
        'elev': elev1,
        'azim': azim1,
        'ibz_kpoints_frac': (
            analysis_result['kpoints_frac_centroid']
            if analysis_result['sg'] not in (1, 2)
            else analysis_result['kpoints_frac']
        ),
        'path_kpoints_frac': analysis_result['kpoints_frac_for_output'],
        'ibz_kpath': analysis_result['kpath_plot'],
        'band_kpoints_frac': analysis_result['band_kpoints_frac'],
        'band_kpath': analysis_result['band_kpath'],
        'extra_general_vertices': analysis_result['extra_general_vertices'],
        # The clipped triclinic half-BZ has no curated hull labels.
        'hull_pts': (
            centroid_result['points_arr']
            if (
                analysis_result['sg'] not in (1, 2)
                or centroid_result['hull'] is not None
            )
            else None
        ),
        'hull_simplices': (
            centroid_result['hull'].simplices.tolist()
            if centroid_result['hull'] is not None
            else None
        ),
        'hull_labels': (
            centroid_result['labels_list']
            if analysis_result['sg'] not in (1, 2)
            else None
        ),
        'sym_ops_cart': analysis_result['sym_ops_cart'],
        'unique_ops': analysis_result['unique_ops'],
        'b_matrix_conv': analysis_result['b_matrix_conv'],
        'b_matrix_input': analysis_result['b_matrix_input'],
        'mode_2d': mode_2d,
        'vacuum_axis': analysis_result['vacuum_axis'],
        'ibz_polygon_frac': centroid_result['ibz_polygon_frac'],
        'ibz_polygon_labels': centroid_result['ibz_polygon_labels'],
        'seekpath_rotation_matrix': np.array(
            analysis_result['sp_result']['rotation_matrix']
        ),
        'standardized_structure_path': diagnostic_result[
            'standardized_structure_path'
        ],
        'standard_mapping_path': diagnostic_result['standard_mapping_path'],
        'symprec': analysis_result['symprec'],
        'mcif_parent_recovery': analysis_result['mcif_parent_recovery'],
        'display_figures': display_figures,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python compute_centroid_hybrid.py <structure_file> [output_dir]")
        sys.exit(1)

    structure_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None
    results = run(structure_file, out_dir)

