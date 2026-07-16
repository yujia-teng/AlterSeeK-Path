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
warnings.filterwarnings("ignore", message="We strongly encourage explicit.*encoding")
warnings.filterwarnings("ignore", message="dict interface is deprecated")
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module=r"seekpath\.hpkot(\..*)-",
)


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import seekpath
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
import spglib

plt.rcParams["mathtext.fontset"] = "stix"


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

    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def _write_seekpath_basis_mapping(input_lattice, standard_lattice, rotation_matrix, output_path):
    """Record the input-cell to SeeK-path-standard axis mapping for users."""
    def _fmt_matrix(mat):
        return "\n".join(
            "  " + " ".join(f"{float(x): .10f}" for x in row)
            for row in np.array(mat, dtype=float)
        )

    input_lattice = np.array(input_lattice, dtype=float)
    standard_lattice = np.array(standard_lattice, dtype=float)
    rotation_matrix = np.array(rotation_matrix, dtype=float)
    lines = [
        "# SeeK-path standardization mapping",
        "# input_lattice is the lattice from the submitted structure file.",
        "# seekpath_standard_lattice is the standardized cell written to *_seekpath_standard.vasp.",
        "# rotation_matrix is reported by SeeK-path for the input-to-standard orientation.",
        "",
        "input_lattice:",
        _fmt_matrix(input_lattice),
        "",
        "seekpath_standard_lattice:",
        _fmt_matrix(standard_lattice),
        "",
        "seekpath_rotation_matrix:",
        _fmt_matrix(rotation_matrix),
        "",
    ]
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
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
    area_centroid_2d,
    ordered_2d_polygon_frac,
    check_input_slab,
    compute_symbolic_centroid,
    get_bz_loops,
    build_symmetry_ibz_cell,
)
from .plotting_common import (
    GAMMA_LABEL,
    _save_figure,
    _print_saved_paths,
)
from .plotting_3d import setup_3d_ax, plot_ibz


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
):
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(filename))
    basename = os.path.splitext(os.path.basename(filename))[0]

    if verbose:
        print("=" * 60)
        print(f"Processing: {filename}")
        print("=" * 60)

    struct = Structure.from_file(filename)
    a_matrix = struct.lattice.matrix
    cell = a_matrix.tolist()
    positions = struct.frac_coords.tolist()
    if seekpath_type_numbers is None:
        numbers = [s.Z for s in struct.species]
    else:
        numbers = [int(n) for n in seekpath_type_numbers]
        if len(numbers) != len(positions):
            raise ValueError(
                "seekpath_type_numbers length must match the number of structure sites."
            )

    # ---- seekpath: lattice detection & standardization ----
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*dict interface is deprecated.*Use attribute interface instead.*",
        )
        warnings.filterwarnings(
            "ignore",
            message=r".*dict interface is deprecated.*",
        )
        sp_result = seekpath.get_path((cell, positions, numbers), with_time_reversal=True)

    input_dataset = spglib.get_symmetry_dataset((cell, positions, numbers))
    os.makedirs(output_dir, exist_ok=True)
    standardized_structure_path = os.path.join(output_dir, f"{basename}_seekpath_standard.vasp")
    standard_mapping_path = os.path.join(output_dir, f"{basename}_seekpath_basis_mapping.txt")
    try:
        _write_seekpath_standard_poscar(
            np.array(input_dataset.std_lattice),
            np.array(input_dataset.std_positions),
            list(input_dataset.std_types),
            standardized_structure_path,
            os.path.basename(filename),
        )
        _write_seekpath_basis_mapping(
            a_matrix,
            np.array(input_dataset.std_lattice),
            np.array(sp_result["rotation_matrix"]),
            standard_mapping_path,
        )
        if verbose:
            print(f"Saved standardized structure: {standardized_structure_path}")
            print(f"Saved SeeK-path basis mapping: {standard_mapping_path}")
    except Exception as exc:
        standardized_structure_path = None
        standard_mapping_path = None
        if verbose:
            print(f"[Warning] Could not write SeeK-path standardization files: {exc}")

    spg_cell = (
        np.array(sp_result['primitive_lattice']),
        np.array(sp_result['primitive_positions']),
        sp_result['primitive_types'],
    )
    dataset = spglib.get_symmetry_dataset(spg_cell)
    b_matrix = np.array(sp_result['reciprocal_primitive_lattice'])
    # Reciprocal lattice of the user-provided input cell. Step 0 rotations
    # are written in this basis by find_sf_operations.py.
    b_matrix_input = 2 * np.pi * np.linalg.inv(np.array(a_matrix)).T

    # Conventional-cell reciprocal lattice (no 2闁?needed --cancels in formula).
    # Used to correctly convert seed flip ops that were written using the
    # conventional cell (input-cell POSCAR and spin_flip_operations.txt).
    _conv_lat = np.array(sp_result.get('conv_lattice', sp_result['primitive_lattice']))
    b_matrix_conv = np.linalg.inv(_conv_lat).T
    b1, b2, b3 = b_matrix

    # ---- 2D / slab mode setup ----
    # Detect the vacuum axis in the *standardized* frame (seekpath may permute
    # axes, e.g. monoclinic -> unique axis b), and sanity-check the input cell.
    vacuum_axis = None
    if mode_2d:
        vacuum_axis, vac_info = detect_vacuum_axis_2d(b_matrix)
        if verbose:
            print(f"\n[2D mode] vacuum axis (standardized frame): "
                  f"{vacuum_axis} ('{'abc'[vacuum_axis]}'); reciprocal norms "
                  f"{[round(x, 4) for x in vac_info['reciprocal_norms']]}")
        if not (vac_info['separated'] and vac_info['orthogonal']):
            print("[2D mode][Warning] vacuum-axis detection is ambiguous "
                  f"(separated={vac_info['separated']}, orthogonal={vac_info['orthogonal']}); "
                  "verify the structure is a proper slab.")
        for w in check_input_slab(a_matrix, input_vacuum_axis):
            print(f"[2D mode][Warning] {w}")

    sg = dataset.number
    laue_group = laue_group_from_point_group(dataset.pointgroup)
    no_altermag = no_altermagnetism_reason(dataset.pointgroup, sg)
    if verbose:
        print(f"\nSpace Group: {sg} ({dataset.international})")
        print(f"Point Group: {dataset.pointgroup}")
        print(f"Laue Group: {laue_group if laue_group is not None else 'Unknown'}")
        if no_altermag:
            print(f"[Note] {no_altermag['reason']} for Laue group {no_altermag['laue_group']}")
        print(f"Seekpath Bravais: {sp_result['bravais_lattice_extended']}")

    # ---- Map to SeeK-path/HPKOT extended Bravais key ----
    sc_type, conv_params = seekpath_to_hpkot_type(sp_result)
    sc_display = sc_type
    centroid_type = sc_type

    # ---- Get SeeK-path band path plus curated HPKOT/project hull points ----
    # SeeK-path is the source of the ordinary band path.  The local HPKOT table
    # is still used for project-curated IBZ hull vertices and copied-sector
    # points that SeeK-path does not define.
    kpath = [
        (_normalize_label(start), _normalize_label(end))
        for start, end in sp_result['path']
    ]
    seekpath_point_coords = {
        _normalize_label(label): list(coords)
        for label, coords in sp_result['point_coords'].items()
    }
    kpoints_frac = get_kpoints(sc_type,
                               conv_params['a'], conv_params.get('b'),
                               conv_params.get('c'), conv_params.get('alpha'))
    path_kpoints_frac = {
        label: seekpath_point_coords[label]
        for label in {label for segment in kpath for label in segment}
        if label in seekpath_point_coords
    }
    kpoints_frac_centroid = get_hull_kpoints(sc_type,
                                             conv_params['a'], conv_params.get('b'),
                                             conv_params.get('c'), conv_params.get('alpha'),
                                             spacegroup_number=sg)
    hull_kpath = get_hull_kpath(sc_type, spacegroup_number=sg)
    display_labels = get_display_labels(sc_type)
    params = get_params(sc_type,
                        conv_params['a'], conv_params.get('b'),
                        conv_params.get('c'), conv_params.get('alpha'))
    if params and verbose:
        print(f"Parameters: {', '.join(f'{k}={v:.6f}' for k, v in params.items())}")

    path_kpoints_cart = {
        k: v[0]*b1 + v[1]*b2 + v[2]*b3
        for k, v in path_kpoints_frac.items()
    }
    kpoints_cart_centroid = {
        k: v[0]*b1 + v[1]*b2 + v[2]*b3
        for k, v in kpoints_frac_centroid.items()
    }

    if sc_type == 'mP1':
        # HPKOT mP1 includes Y and C as path labels, but they are not vertices
        # of the selected simple-monoclinic IBZ hull.
        for label in ('Y', 'C'):
            kpoints_frac_centroid.pop(label, None)
            kpoints_cart_centroid.pop(label, None)
    elif sc_type == 'mC1':
        # Keep project-only hidden closure vertices for the mC1 hull. They are
        # present in kpoints_frac_centroid but absent from display labels/path.
        pass

    # ---- Symmetry operations ----
    sym_ops_cart, unique_ops = get_symmetry_operations(b_matrix, dataset)
    if verbose:
        print(f"\nSymmetry operations: {len(sym_ops_cart)}")
        print(f"With time-reversal: {len(unique_ops)}")

    # ---- Convex Hull & Centroid ----
    # Use the project-curated HPKOT hull point set.  This is distinct from the
    # public HPKOT band path: lower-symmetry classes can add copied boundary
    # vertices (e.g. M_A/L_A for hexagonal 6/m) to preserve the standard doubled
    # IBZ wedge used by the project.
    labels_list = list(kpoints_cart_centroid.keys())
    points_arr = np.array([kpoints_cart_centroid[k] for k in labels_list])

    ibz_polygon_frac = None
    if mode_2d:
        # 2D slab: the physical IBZ is the k[vacuum_axis]=0 cross-section.
        # Restrict the curated hull points to that plane and take the 2D area
        # centroid (the 3D volume centroid/ConvexHull are meaningless here and
        # ConvexHull crashes on coplanar input).
        in_plane_labels = [
            lab for lab in labels_list
            if abs(kpoints_frac_centroid[lab][vacuum_axis]) < 1e-4
        ]
        # A "<base>_2" HPKOT point that is genuinely a distinct 3D hull
        # vertex (e.g. mP1's B_2/D_2, needed for the real 3D IBZ volume) can
        # still be a redundant mirror-image duplicate once sliced to this 2D
        # plane, if both "<base>" and "<base>_2" happen to survive the
        # in-plane filter. Keeping both would double the in-plane hull's
        # extent along whichever coordinate differs between them, pushing
        # Gamma off being a hull corner and skewing the 2D area centroid.
        # This is a 2D-only reduction -- the curated 3D hull table itself is
        # untouched, since 3D genuinely needs both points.
        base_labels = {lab for lab in in_plane_labels if not lab.endswith('_2')}
        in_plane_labels = [
            lab for lab in in_plane_labels
            if not (lab.endswith('_2') and lab[:-2] in base_labels)
        ]
        frac_pts = np.array([kpoints_frac_centroid[lab] for lab in in_plane_labels])
        hull = None
        centroid_frac, centroid_cart, ibz_vol = area_centroid_2d(
            frac_pts, vacuum_axis, b_matrix)
        ibz_polygon_frac = ordered_2d_polygon_frac(frac_pts, vacuum_axis)
        if verbose:
            print(f"[2D mode] in-plane hull labels: {len(in_plane_labels)} / "
                  f"{len(labels_list)} ({', '.join(in_plane_labels)})")
            print(f"[2D mode] area centroid (frac): "
                  f"[{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]"
                  f"  area={ibz_vol:.6f}")
    elif sg in (1, 2):
        # Triclinic: IBZ boundary is hard to define on Wigner-Seitz BZ.
        # Skip hull/centroid --not needed (no altermagnetic splitting).
        hull = None
        centroid_cart = np.mean(points_arr, axis=0)
        centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)
        ibz_vol = 0.0
        if verbose:
            print(f"\n[Note] Triclinic: IBZ shading skipped (IBZ = {'full BZ' if sg == 1 else 'half BZ'})")
            print(f"Centroid (mean of k-points): [{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]")
    else:
        hull = ConvexHull(points_arr)
        centroid_cart, ibz_vol = calculate_volume_centroid(hull)
        centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)

        # The mC2/mC3 HPKOT tables include distinctive boundary labels whose
        # convex hull is slightly larger than the true C2/m fundamental domain.
        # For these branches only, keep the labels/path from HPKOT but shade and
        # use the centroid of the symmetry/Voronoi IBZ cell so four images tile
        # the BZ without overlap.
        hull_matches_labels = True
        if sc_type in {'mC2', 'mC3'}:
            mono_pts, mono_simplices = build_symmetry_ibz_cell(
                b_matrix, unique_ops, centroid_cart)
            if mono_pts is not None and mono_simplices is not None:
                points_arr = np.array(mono_pts, dtype=float)
                hull = ConvexHull(points_arr)
                centroid_cart, ibz_vol = calculate_volume_centroid(hull)
                centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)
                hull_matches_labels = False
                if verbose:
                    print("[Note] Using symmetry/Voronoi IBZ cell for monoclinic hull.")

        # ---- Symbolic Centroid (saved to file, not printed) ----
        if hull_matches_labels:
            try:
                sym_centroid, param_syms = compute_symbolic_centroid(
                    kpoints_frac_centroid, hull, labels_list, centroid_type, conv_params)
                if sym_centroid is not None:
                    sym_lines = "\n".join(
                        f"  {ax_name} = {sym_centroid[i]}"
                        for i, ax_name in enumerate(['k1', 'k2', 'k3'])
                    )
                    try:
                        with open("spin_operations.txt", "a", encoding="utf-8", newline="\n") as f:
                            f.write(f"\nSymbolic IBZ centroid (fractional):\n{sym_lines}\n")
                    except Exception:
                        pass
            except Exception:
                pass

    # ---- Plotting ----
    bz_loops = get_bz_loops(b_matrix)
    all_bz_pts = np.vstack(bz_loops)
    bz_center = np.mean(all_bz_pts, axis=0)
    bz_span = np.max(all_bz_pts) - np.min(all_bz_pts)

    # Selected band path. For doubled-IBZ cases, copied vertices may either be
    # included in project-defined path segments or appended as isolated
    # general-point anchors when no nonredundant high-symmetry edge is needed.
    #
    # The doubled-IBZ construction (SG 75-88 tP1/tI1/tI2, SG 168-176 hP2)
    # doubles the wedge because the *actual* point group (4, -4, 4/m, 6, -6,
    # 6/m, 3, -3) has no vertical mirrors -- that deficiency is an in-plane
    # fact, not a k_z/bulk-only one, so it still applies in 2D. The curated
    # hull (kpoints_frac_centroid, from get_hull_kpoints below) always
    # includes the doubled _A points regardless of mode_2d, and the 2D-mode
    # centroid/ibz_polygon_frac block further down already uses that doubled
    # hull. So band_kpath must also use the doubled hull_kpath in 2D mode --
    # otherwise the *displayed* wedge is only half of what the (correctly
    # doubled) centroid was computed from, and the plotted k point ends up
    # looking off-center/outside a too-small triangle.
    doubled_ibz_case = (75 <= sg <= 88 and sc_type in {'tP1', 'tI1', 'tI2'}) or (
        168 <= sg <= 176 and sc_type == 'hP2'
    )
    if doubled_ibz_case:
        band_kpath = list(hull_kpath)
    else:
        band_kpath = list(kpath)
    band_kpoints_frac = dict(path_kpoints_frac)
    extra_general_vertices = []
    if sc_type == 'hP1' and sg in {149, 151, 153, 157, 159, 162, 163}:
        extra_general_vertices = ["H_2", "A_2", "L_A"]
        for label in extra_general_vertices:
            if label in kpoints_frac_centroid:
                band_kpoints_frac[label] = kpoints_frac_centroid[label]
    elif sc_type == 'hP2' and 149 <= sg <= 167:
        extra_general_vertices = ["H_2", "A_2", "L_A"]
        for label in extra_general_vertices:
            if label in kpoints_frac_centroid:
                band_kpoints_frac[label] = kpoints_frac_centroid[label]
    elif sc_type == 'hP2' and 168 <= sg <= 176:
        extra_general_vertices = ["M_A", "L_A"]
        for label in ("L_A", "M_A"):
            if label in kpoints_frac_centroid:
                band_kpoints_frac[label] = kpoints_frac_centroid[label]
    elif 75 <= sg <= 88 and sc_type in {'tP1', 'tI1', 'tI2'}:
        tetragonal_extra = {
            'tP1': ["X_A", "R_A"],
            'tI1': ["M_A", "N_A", "Z_0A"],
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
        label for label in extra_general_vertices
        if label not in path_labels
    ]

    if mode_2d:
        # Keep only path segments and anchors that live in the physical plane
        # (both endpoints have fractional k[vacuum_axis] ~ 0). Out-of-plane
        # high-symmetry points (e.g. Z, A, R) are the same in-plane points
        # raised along the dead vacuum reciprocal direction and are dropped.
        def _in_plane_label(lbl):
            v = (band_kpoints_frac.get(lbl)
                 or kpoints_frac_centroid.get(lbl)
                 or path_kpoints_frac.get(lbl))
            return v is not None and abs(v[vacuum_axis]) < 1e-4
        band_kpath = [
            seg for seg in band_kpath if all(_in_plane_label(l) for l in seg)
        ]
        extra_general_vertices = [
            l for l in extra_general_vertices if _in_plane_label(l)
        ]
        path_labels = {label for segment in band_kpath for label in segment}

        if doubled_ibz_case:
            # The 3D doubled-IBZ chain (e.g. tP1's A-R_A-X_A-M) only reaches
            # its copied _A vertex via out-of-plane anchors (A, R_A, ...),
            # which the plane filter above just dropped. That leaves the
            # surviving in-plane tail (e.g. X_A-M) as a disconnected stub
            # with no edge back to Gamma, even though X_A is a genuine
            # vertex of the physical in-plane doubled wedge (Gamma-X-M-X_A).
            # Close the loop explicitly: any doubled _A/_B-style vertex that
            # is in-plane but not directly connected to Gamma gets a new
            # Gamma-<vertex> segment.
            def _is_doubled_label(lbl):
                base = lbl[1:] if lbl.startswith('_') else lbl
                return '_A' in base or '_B' in base or base.endswith('_0A')

            connected_to_gamma = {
                seg[1] if seg[0] == GAMMA_LABEL else seg[0]
                for seg in band_kpath if GAMMA_LABEL in seg
            }
            doubled_labels = [l for l in path_labels if _is_doubled_label(l)]
            for label in doubled_labels:
                if label not in connected_to_gamma:
                    band_kpath.append((GAMMA_LABEL, label))
            path_labels = {label for segment in band_kpath for label in segment}

        if verbose:
            print(f"[2D mode] in-plane band path: {len(band_kpath)} segments, "
                  f"labels {sorted(path_labels)}")

    # For plotting: draw the selected band path on top of the project IBZ hull.
    # Add path-only optional points (e.g. H_2) only when the selected path
    # actually uses them; do not let them affect the hull or centroid.
    kpath_plot = band_kpath
    display_labels_plot = display_labels
    kpoints_cart_plot = dict(kpoints_cart_centroid)
    for label in {lbl for segment in kpath_plot for lbl in segment}:
        if label in path_kpoints_cart:
            kpoints_cart_plot[label] = path_kpoints_cart[label]

    kpoints_frac_for_output = dict(path_kpoints_frac)
    for label in {lbl for segment in kpath_plot for lbl in segment}:
        if label in kpoints_frac_centroid:
            kpoints_frac_for_output[label] = kpoints_frac_centroid[label]

    # ---- Plotting ----
    fig1_title = (f"BZ: {basename} ({sc_display})" if sg in (1, 2)
                  else f"IBZ + BZ: {basename} ({sc_display})")

    display_figures = []
    default_elev, default_azim = (
        (view_elev, view_azim) if view_elev is not None and view_azim is not None
        else (14, 20)
    )
    elev1, azim1 = default_elev, default_azim
    fig1_path = os.path.join(output_dir, f'{basename}_ibz_{sc_display}.png')
    if show_plot and not mode_2d:
        # Interactive mode: create the figure now. alterseek_path can defer
        # the actual plt.show() call until all prompts and file writes finish.
        # Opens at (default_elev, default_azim) so a fixed camera angle (e.g.
        # matched across cases) is used unless the user rotates it manually.
        fig1, ax1 = setup_3d_ax(fig1_title,
                                bz_loops, b_matrix, bz_center, bz_span,
                                elev=default_elev, azim=default_azim)
        plot_ibz(ax1, kpoints_cart_plot, kpath_plot, display_labels_plot,
                 hull, centroid_cart, hull_pts=points_arr, lattice_type=sc_type,
                 hull_labels=labels_list)
        plt.tight_layout()
        if defer_show:
            def _save_fig1_after_show(fig=fig1, ax=ax1):
                print(f"[View] elev={ax.elev:.2f}, azim={ax.azim:.2f}")
                fig1s, ax1s = setup_3d_ax(fig1_title,
                                          bz_loops, b_matrix, bz_center, bz_span,
                                          elev=ax.elev, azim=ax.azim,
                                          dashed_back=True)
                plot_ibz(ax1s, kpoints_cart_plot, kpath_plot,
                         display_labels_plot, hull, centroid_cart,
                         hull_pts=points_arr, lattice_type=sc_type,
                         hull_labels=labels_list)
                plt.tight_layout()
                saved_paths = _save_figure(fig1s, fig1_path, dpi=300, bbox_inches='tight')
                plt.close(fig1s)
                plt.close(fig)
                _print_saved_paths(saved_paths)
            fig1._alterseek_save_after_show = _save_fig1_after_show
            display_figures.append(fig1)
            elev1, azim1 = default_elev, default_azim
        else:
            plt.show()
            elev1, azim1 = ax1.elev, ax1.azim
    else:
        # Automated mode (called from alterseek_path): use default angles, no window
        elev1, azim1 = default_elev, default_azim

    # Render with dashed back-edges and save unless deferred post-show saving is active
    # (skipped in 2D mode: Figures stay 3D-only and are not produced for slabs yet)
    if not (show_plot and defer_show) and not mode_2d:
        fig1s, ax1s = setup_3d_ax(fig1_title,
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  elev=elev1, azim=azim1, dashed_back=True)
        plot_ibz(ax1s, kpoints_cart_plot, kpath_plot, display_labels_plot,
                 hull, centroid_cart, hull_pts=points_arr, lattice_type=sc_type,
                 hull_labels=labels_list)
        plt.tight_layout()
        saved_paths = _save_figure(fig1s, fig1_path, dpi=300, bbox_inches='tight')
        _print_saved_paths(saved_paths, verbose=verbose)
        plt.close(fig1s)

    return {
        'sc_type': sc_display,
        'lattice_key': sc_type,
        'seekpath_bravais': sp_result['bravais_lattice_extended'],
        'spacegroup': sg,
        'sg_symbol': dataset.international,
        'point_group': dataset.pointgroup,
        'laue_group': laue_group,
        'no_altermagnetism': no_altermag,
        'kpoints_frac': kpoints_frac,
        'centroid_cart': centroid_cart,
        'centroid_frac': centroid_frac,
        'ibz_volume': ibz_vol,
        'n_symmetry_ops': len(unique_ops),
        'sp_path': sp_result['path'],
        'sp_point_coords': sp_result['point_coords'],
        'b_matrix': b_matrix,
        'bz_loops': bz_loops,
        'bz_center': bz_center,
        'bz_span': bz_span,
        'elev': elev1,
        'azim': azim1,
        'ibz_kpoints_frac': kpoints_frac_centroid if sg not in (1, 2) else kpoints_frac,
        'path_kpoints_frac': kpoints_frac_for_output,
        'ibz_kpath': kpath_plot,
        'band_kpoints_frac': band_kpoints_frac,
        'band_kpath': band_kpath,
        'extra_general_vertices': extra_general_vertices,
        'hull_pts': points_arr if sg not in (1, 2) else None,
        'hull_simplices': hull.simplices.tolist() if (sg not in (1, 2) and hull is not None) else None,
        'hull_labels': labels_list if sg not in (1, 2) else None,
        'sym_ops_cart': sym_ops_cart,
        'unique_ops': unique_ops,
        'b_matrix_conv': b_matrix_conv,
        'b_matrix_input': b_matrix_input,
        'mode_2d': mode_2d,
        'vacuum_axis': vacuum_axis,
        'ibz_polygon_frac': ibz_polygon_frac,
        'seekpath_rotation_matrix': np.array(sp_result['rotation_matrix']),
        'standardized_structure_path': standardized_structure_path,
        'standard_mapping_path': standard_mapping_path,
        'display_figures': display_figures,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python compute_centroid_hybrid.py <structure_file> [output_dir]")
        sys.exit(1)

    structure_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None
    results = run(structure_file, out_dir)

