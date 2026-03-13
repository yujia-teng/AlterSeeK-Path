#!/usr/bin/env python3
"""
IBZ Centroid Calculator (Hybrid: seekpath + Setyawan-Curtarolo)
================================================================
Uses seekpath for:  lattice type detection, cell standardization
Uses our own data:  IRBZ k-point vertices (Setyawan-Curtarolo convention)

This ensures the IRBZ shape is consistent for all space groups within
the same Bravais lattice type — no extra points like H_2 in hexagonal.

Supports ALL 25 Bravais lattice variations (including triclinic).

Usage:
    python compute_centroid_hybrid.py <structure_file>
    python compute_centroid_hybrid.py <structure_file> <output_dir>

Requires:
    pip install seekpath pymatgen spglib numpy scipy matplotlib sympy
"""

import sys
import os
import warnings
import threading
import atexit
warnings.filterwarnings("ignore", message="We strongly encourage explicit.*encoding")
warnings.filterwarnings("ignore", message="dict interface is deprecated")
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    module=r"seekpath\.hpkot(\..*)?",
)


def _suppress_stderr_lines(tokens):
    """Filter selected stderr lines (including native C-level writes)."""
    try:
        orig_fd = os.dup(2)
        r_fd, w_fd = os.pipe()
        os.dup2(w_fd, 2)
    except Exception:
        return

    running = {"on": True}

    def _pump():
        buf = ""
        while running["on"]:
            try:
                chunk = os.read(r_fd, 4096)
                if not chunk:
                    break
                buf += chunk.decode("utf-8", errors="replace")
                while "\n" in buf:
                    line, buf = buf.split("\n", 1)
                    if not any(tok in line for tok in tokens):
                        os.write(orig_fd, (line + "\n").encode("utf-8", errors="replace"))
            except Exception:
                break

    def _stop():
        running["on"] = False
        try:
            os.dup2(orig_fd, 2)
        except Exception:
            pass
        for fd in (w_fd, r_fd, orig_fd):
            try:
                os.close(fd)
            except Exception:
                pass

    threading.Thread(target=_pump, daemon=True).start()
    atexit.register(_stop)


_suppress_stderr_lines(("libpng warning: iCCP: known incorrect sRGB profile",))


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sympy as sp
import seekpath
from pymatgen.core import Structure
import spglib

from lattice_kpoints import (
    LATTICE_DATA, get_kpoints, get_kpath, get_display_labels, get_params,
)
# ============================================================================
# Seekpath Bravais → Setyawan-Curtarolo BZ variation mapping
# ============================================================================
def seekpath_to_sc_type(sp_result):
    """
    Map seekpath's Bravais label to Setyawan-Curtarolo BZ variation.

    Returns
    -------
    sc_type : str (e.g. 'CUB', 'FCC', 'BCT1', 'TRI1a', ...)
    conv_params : dict with 'a', 'b', 'c', 'alpha' for parametric types
    """
    bravais = sp_result['bravais_lattice']  # e.g. 'cF', 'hP', 'oI', 'aP'

    # Get conventional cell parameters
    conv_lattice = np.array(sp_result.get('conv_lattice',
                            sp_result.get('primitive_lattice')))
    va, vb, vc = conv_lattice[0], conv_lattice[1], conv_lattice[2]
    a = np.linalg.norm(va)
    b = np.linalg.norm(vb)
    c = np.linalg.norm(vc)
    alpha = np.degrees(np.arccos(np.clip(np.dot(vb, vc)/(b*c), -1, 1)))
    beta = np.degrees(np.arccos(np.clip(np.dot(va, vc)/(a*c), -1, 1)))

    conv_params = {'a': a, 'b': b, 'c': c, 'alpha': alpha}

    # --- Simple mappings ---
    simple = {
        'cP': 'CUB', 'cF': 'FCC', 'cI': 'BCC',
        'tP': 'TET', 'oP': 'ORC', 'hP': 'HEX',
    }
    if bravais in simple:
        return simple[bravais], conv_params

    # --- Triclinic ---
    if bravais == 'aP':
        recip_latt = np.array(sp_result['reciprocal_primitive_lattice'])
        rb1, rb2, rb3 = recip_latt
        k_alpha = np.degrees(np.arccos(np.clip(
            np.dot(rb2, rb3) / (np.linalg.norm(rb2)*np.linalg.norm(rb3)), -1, 1)))
        k_beta = np.degrees(np.arccos(np.clip(
            np.dot(rb1, rb3) / (np.linalg.norm(rb1)*np.linalg.norm(rb3)), -1, 1)))
        k_gamma = np.degrees(np.arccos(np.clip(
            np.dot(rb1, rb2) / (np.linalg.norm(rb1)*np.linalg.norm(rb2)), -1, 1)))

        tol = 0.1  # degree tolerance
        all_leq_90 = (k_alpha <= 90+tol) and (k_beta <= 90+tol) and (k_gamma <= 90+tol)
        all_geq_90 = (k_alpha >= 90-tol) and (k_beta >= 90-tol) and (k_gamma >= 90-tol)
        any_eq_90 = (abs(k_alpha-90) < tol) or (abs(k_beta-90) < tol) or (abs(k_gamma-90) < tol)

        if all_leq_90:
            return ('TRI1b' if any_eq_90 else 'TRI1a'), conv_params
        elif all_geq_90:
            return ('TRI2b' if any_eq_90 else 'TRI2a'), conv_params
        else:
            print(f"  WARNING: Triclinic reciprocal angles mixed "
                  f"({k_alpha:.1f}°, {k_beta:.1f}°, {k_gamma:.1f}°)")
            return 'TRI2a', conv_params

    # --- Rhombohedral ---
    if bravais == 'hR':
        prim_latt = np.array(sp_result['primitive_lattice'])
        pa = np.linalg.norm(prim_latt[0])
        rhomb_alpha = np.degrees(np.arccos(np.clip(
            np.dot(prim_latt[0], prim_latt[1]) / (pa**2), -1, 1)))
        conv_params['alpha'] = rhomb_alpha
        conv_params['a'] = pa
        return ('RHL1' if rhomb_alpha < 90 else 'RHL2'), conv_params

    # --- Body-centered tetragonal ---
    if bravais == 'tI':
        return ('BCT1' if c < a else 'BCT2'), conv_params

    # --- Face-centered orthorhombic ---
    if bravais == 'oF':
        a_s, b_s, c_s = sorted([a, b, c])
        conv_params['a'], conv_params['b'], conv_params['c'] = a_s, b_s, c_s
        inv_a2, inv_b2, inv_c2 = 1/a_s**2, 1/b_s**2, 1/c_s**2
        if abs(inv_a2 - (inv_b2 + inv_c2)) < 1e-8:
            return 'ORCF3', conv_params
        elif inv_a2 > inv_b2 + inv_c2:
            return 'ORCF1', conv_params
        else:
            return 'ORCF2', conv_params

    # --- Body-centered orthorhombic ---
    if bravais == 'oI':
        a_s, b_s, c_s = sorted([a, b, c])
        conv_params['a'], conv_params['b'], conv_params['c'] = a_s, b_s, c_s
        return 'ORCI', conv_params

    # --- Base-centered orthorhombic (seekpath variants: oS, oC, oA, oB) ---
    if bravais in ('oS', 'oC', 'oA', 'oB'):
        if a > b: a, b = b, a
        conv_params['a'], conv_params['b'] = a, b
        return 'ORCC', conv_params

    # --- Simple monoclinic ---
    if bravais == 'mP':
        # MCL params in lattice_kpoints use monoclinic beta (angle a-c).
        conv_params['alpha'] = beta
        return 'MCL', conv_params

    # --- C-centered monoclinic ---
    # HPKOT/seekpath: mC1/mC2/mC3 (3 cases)
    # Setyawan-Curtarolo: MCLC1-5 (5 cases, adding boundary cases MCLC2 and MCLC4)
    if bravais in ('mS', 'mC'):
        alpha_rad = np.radians(beta)
        conv_params['alpha'] = beta
        tol_mc = 1e-6
        b_sin = a * np.sin(alpha_rad)
        cond = -a * np.cos(alpha_rad) / c + a**2 * np.sin(alpha_rad)**2 / b**2
        if abs(b - b_sin) < tol_mc:
            return 'MCLC2_SC', conv_params  # SC: MCLC2 (k_gamma = 90, boundary)
        if b < b_sin:
            return 'MCLC1', conv_params
        if abs(cond - 1.0) < tol_mc:
            return 'MCLC4_SC', conv_params  # SC: MCLC4 (cond = 1, boundary)
        if cond < 1.0:
            return 'MCLC2', conv_params
        return 'MCLC3', conv_params

    print(f"WARNING: Unknown Bravais type '{bravais}', defaulting to CUB")
    return 'CUB', conv_params


def _seekpath_label_to_internal(label):
    return 'Γ' if label == 'GAMMA' else label


def _display_label_from_internal(label):
    if label == 'Γ':
        return r'$\Gamma$'
    if '_' in label:
        base, sub = label.split('_', 1)
        return rf'${base}_{sub}$'
    return label


# ============================================================================
# Symmetry Operations
# ============================================================================
def get_symmetry_operations(b_matrix, dataset):
    """Convert real-space rotations to k-space, add time-reversal."""
    b_mat_T = b_matrix.T
    b_mat_T_inv = np.linalg.inv(b_mat_T)

    sym_ops_cart = [b_mat_T @ np.linalg.inv(R).T @ b_mat_T_inv
                    for R in dataset.rotations]

    all_ops = [op for R in sym_ops_cart for op in (R, -R)]
    unique_ops = []
    for op in all_ops:
        if not any(np.allclose(op, ex, atol=1e-6) for ex in unique_ops):
            unique_ops.append(op)

    return sym_ops_cart, unique_ops


# ============================================================================
# Centroid Calculation
# ============================================================================
def calculate_volume_centroid(hull):
    """Compute volume centroid via signed tetrahedra decomposition."""
    ref = np.mean(hull.points[hull.vertices], axis=0)
    total_vol = 0.0
    w_cent = np.zeros(3)
    for simplex in hull.simplices:
        a, b, c = hull.points[simplex[0]], hull.points[simplex[1]], hull.points[simplex[2]]
        vol = np.abs(np.dot(a - ref, np.cross(b - ref, c - ref))) / 6.0
        total_vol += vol
        w_cent += vol * (ref + a + b + c) / 4.0
    return w_cent / total_vol, total_vol


def compute_symbolic_centroid(kpoints_frac, hull, labels_list, lattice_type, conv_params):
    """Compute symbolic centroid (exact fractions or parametric)."""
    data = LATTICE_DATA[lattice_type]

    if 'kpoints' in data:
        kp_sym = {k: [sp.nsimplify(c, rational=True) for c in v]
                  for k, v in data['kpoints'].items()}
        param_symbols = {}
    elif 'params_func' in data:
        actual = data['params_func'](
            conv_params['a'], conv_params.get('b', conv_params['a']),
            conv_params.get('c', conv_params['a']),
            conv_params.get('alpha', 90.0))
        param_symbols = {p: sp.Symbol(p, real=True, positive=True) for p in actual}
        kp_from_func = data['kpoints_func'](param_symbols)
        kp_sym = {k: [sp.nsimplify(c, rational=True) if isinstance(c, (int, float)) else c
                       for c in v] for k, v in kp_from_func.items()}
    else:
        return None, {}

    sym_points = [sp.Matrix(kp_sym[k]) for k in labels_list]
    sym_ref = sum([sym_points[i] for i in hull.vertices], sp.zeros(3, 1)) / len(hull.vertices)

    sym_total_vol = sp.Integer(0)
    sym_weighted_centroid = sp.zeros(3, 1)

    if 'params_func' in data:
        num_params = data['params_func'](
            conv_params['a'], conv_params.get('b', conv_params['a']),
            conv_params.get('c', conv_params['a']),
            conv_params.get('alpha', 90.0))
        subs_list = [(param_symbols[k], num_params[k]) for k in param_symbols]
    else:
        subs_list = []

    for simplex in hull.simplices:
        a_s, b_s, c_s = sym_points[simplex[0]], sym_points[simplex[1]], sym_points[simplex[2]]
        det_val = sp.Matrix([(a_s-sym_ref).T, (b_s-sym_ref).T, (c_s-sym_ref).T]).det()
        num_det = float(det_val.subs(subs_list)) if subs_list else float(det_val)
        sign = 1 if num_det >= 0 else -1
        vol = sign * det_val / 6
        sym_total_vol += vol
        sym_weighted_centroid += vol * (sym_ref + a_s + b_s + c_s) / 4

    raw_centroid = sp.Matrix(sym_weighted_centroid / sym_total_vol)
    sym_centroid = simplify_symbolic_centroid(raw_centroid, lattice_type, param_symbols)
    return sym_centroid, param_symbols
def _relation_candidates(lattice_type, param_symbols):
    """
    Return substitution candidates used to eliminate dependent symbols.
    Extend this map for other parametric lattice types as needed.
    """
    eta = param_symbols.get('eta')
    nu = param_symbols.get('nu')

    candidates = []
    if lattice_type in ('RHL1', 'RHL2') and eta is not None and nu is not None:
        candidates.append({nu: sp.Rational(3, 4) - eta / 2})
        candidates.append({eta: sp.Rational(3, 2) - 2 * nu})
    return candidates


def _expr_complexity(expr):
    """Lower is simpler."""
    return (sp.count_ops(expr), len(str(expr)))


def simplify_symbolic_centroid(expr_vec, lattice_type, param_symbols):
    """
    Simplify centroid expressions and optionally apply known parameter relations.
    Chooses the least complex equivalent form.
    """
    base = sp.Matrix([sp.simplify(sp.together(e)) for e in expr_vec])
    best = base
    best_score = sum(_expr_complexity(e)[0] for e in base), sum(_expr_complexity(e)[1] for e in base)

    for sub_map in _relation_candidates(lattice_type, param_symbols):
        cand = sp.Matrix([sp.simplify(sp.together(e.subs(sub_map))) for e in expr_vec])
        score = sum(_expr_complexity(e)[0] for e in cand), sum(_expr_complexity(e)[1] for e in cand)
        if score < best_score:
            best, best_score = cand, score

    return best


# ============================================================================
# BZ Boundary & Plotting
# ============================================================================
def get_bz_loops(b_matrix):
    grid = np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1])).T.reshape(-1,3)
    points = grid @ b_matrix
    vor = Voronoi(points)
    origin_idx = 13
    loops = []
    for i, pair in enumerate(vor.ridge_points):
        if origin_idx not in pair: continue
        idx = vor.ridge_vertices[i]
        if -1 in idx: continue
        pts = vor.vertices[idx]
        center = np.mean(pts, axis=0)
        neighbor = points[pair[0] if pair[1] == origin_idx else pair[1]]
        normal = neighbor - points[origin_idx]
        normal /= np.linalg.norm(normal)
        ref = np.array([0.,0.,1.]) if np.abs(normal[2]) < 0.9 else np.array([0.,1.,0.])
        u = np.cross(normal, ref); u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        angles = np.arctan2((pts-center)@v, (pts-center)@u)
        loop = pts[np.argsort(angles)]
        loops.append(np.vstack([loop, loop[0]]))
    return loops


def find_bz_exit(vec, b_matrix):
    grid = np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1])).T.reshape(-1,3)
    G_vectors = grid @ b_matrix
    t_min = np.inf
    for G in G_vectors:
        dot = np.dot(vec, G)
        if dot > 1e-10:
            t = np.dot(G, G) / (2 * dot)
            if t < t_min: t_min = t
    return t_min


def _get_view_direction(ax):
    """Get the unit vector pointing from the scene toward the camera."""
    elev = np.radians(ax.elev)
    azim = np.radians(ax.azim)
    # Camera direction in data coordinates
    view = np.array([
        np.cos(elev) * np.cos(azim),
        np.cos(elev) * np.sin(azim),
        np.sin(elev),
    ])
    return view


def _classify_bz_edges(bz_loops, view_dir):
    """
    Classify BZ edges as front or back based on adjacent face normals.

    Each loop is one BZ face.  An edge shared by two faces is 'back' if
    BOTH adjacent face normals point away from the viewer.
    An edge belonging to only one face is 'back' if that face normal
    points away from the viewer.
    """
    # Build edge → face-normal map
    from collections import defaultdict
    edge_normals = defaultdict(list)

    face_normals = []
    for loop in bz_loops:
        pts = loop[:-1]  # remove closing duplicate
        center = np.mean(pts, axis=0)
        # Face normal (cross product of two edge vectors)
        if len(pts) >= 3:
            n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
            # Orient outward (away from origin)
            if np.dot(n, center) < 0:
                n = -n
            n = n / (np.linalg.norm(n) + 1e-15)
        else:
            n = np.array([0., 0., 0.])
        face_normals.append(n)

        for i in range(len(pts)):
            p1 = tuple(np.round(pts[i], 8))
            p2 = tuple(np.round(pts[(i+1) % len(pts)], 8))
            edge_key = (min(p1, p2), max(p1, p2))
            edge_normals[edge_key].append(n)

    front_edges = []
    back_edges = []

    for edge_key, normals in edge_normals.items():
        # Edge is front if ANY adjacent face is front-facing
        is_front = any(np.dot(n, view_dir) > 1e-6 for n in normals)
        seg = np.array([list(edge_key[0]), list(edge_key[1])])
        if is_front:
            front_edges.append(seg)
        else:
            back_edges.append(seg)

    return front_edges, back_edges


def draw_bz_edges(ax, bz_loops, dashed_back=False):
    """Draw BZ edges. If dashed_back, use view-dependent solid/dashed."""
    if dashed_back:
        view_dir = _get_view_direction(ax)
        front, back = _classify_bz_edges(bz_loops, view_dir)
        for seg in back:
            ax.plot(seg[:,0], seg[:,1], seg[:,2],
                    c='black', ls='--', lw=1.0, alpha=0.4)
        for seg in front:
            ax.plot(seg[:,0], seg[:,1], seg[:,2],
                    c='black', ls='-', lw=1.5, alpha=0.7)
    else:
        for loop in bz_loops:
            ax.plot(loop[:,0], loop[:,1], loop[:,2],
                    c='black', ls='-', lw=1.5, alpha=0.6)


def setup_3d_ax(title, bz_loops, b_matrix, bz_center, bz_span,
                elev=25, azim=-55, dashed_back=False):
    b1, b2, b3 = b_matrix[0], b_matrix[1], b_matrix[2]
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=elev, azim=azim)
    draw_bz_edges(ax, bz_loops, dashed_back=dashed_back)
    vec_labels = [r'$\mathbf{b}_1$', r'$\mathbf{b}_2$', r'$\mathbf{b}_3$']
    for i, vec in enumerate([b1, b2, b3]):
        t_exit = find_bz_exit(vec, b_matrix)
        t_exit = min(t_exit, 1.0)
        exit_pt = vec * t_exit
        # Dotted line inside BZ
        ax.plot([0,exit_pt[0]], [0,exit_pt[1]], [0,exit_pt[2]],
                color='black', ls=':', lw=1.5, alpha=0.6, zorder=100)
        # Short arrow segment outside BZ (30% of outside length)
        outside = vec - exit_pt
        arrow_frac = 0.35
        arrow_end = exit_pt + outside * arrow_frac
        ax.quiver(exit_pt[0], exit_pt[1], exit_pt[2],
                  outside[0]*arrow_frac, outside[1]*arrow_frac, outside[2]*arrow_frac,
                  color='black', arrow_length_ratio=0.4, lw=2.0, zorder=100)
        ax.text(arrow_end[0] + outside[0]*0.08, arrow_end[1] + outside[1]*0.08,
                arrow_end[2] + outside[2]*0.08, vec_labels[i],
                color='black', fontsize=20, fontweight='bold', zorder=101)
    # Use per-axis ranges so the BZ isn't squashed along short axes
    all_pts = np.vstack([np.array(loop) for loop in bz_loops])
    ranges = np.ptp(all_pts, axis=0)  # [dx, dy, dz]
    pad = 0.25  # fractional padding around BZ
    for i, (set_lim, r) in enumerate(zip(
            [ax.set_xlim, ax.set_ylim, ax.set_zlim], ranges)):
        half = r / 2 * (1 + pad)
        set_lim(bz_center[i] - half, bz_center[i] + half)
    ax.set_box_aspect(ranges / ranges.max())
    ax.set_axis_off()
    ax.set_title(title, fontsize=14)
    return fig, ax


def plot_ibz(ax, kpoints_cart, kpath, display_labels, hull, centroid_cart):
    points_list = list(kpoints_cart.values())
    ibz_faces = [[points_list[s] for s in simplex] for simplex in hull.simplices]
    ax.add_collection3d(Poly3DCollection(
        ibz_faces, facecolor='steelblue', alpha=0.15, edgecolor='none'))
    for k1, k2 in kpath:
        if k1 in kpoints_cart and k2 in kpoints_cart:
            p1, p2 = kpoints_cart[k1], kpoints_cart[k2]
            ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]],
                    c='red', ls='-', lw=2.5, alpha=0.9)
    ibz_center = np.mean(points_list, axis=0)
    ibz_span = np.max(np.ptp(np.array(points_list), axis=0))
    label_offset = ibz_span * 0.1  # scale offset to IBZ size
    for label, coords in kpoints_cart.items():
        if label.startswith('_'):
            continue  # hidden vertex used only for IRBZ hull
        ax.scatter(coords[0], coords[1], coords[2],
                   c='red', s=80, zorder=110, edgecolors='darkred', linewidths=0.5)
        direction = coords - ibz_center
        norm_dir = np.linalg.norm(direction)
        offset = direction / norm_dir * label_offset if norm_dir > 1e-8 else np.array([0, 0, label_offset])
        ax.text(coords[0]+offset[0], coords[1]+offset[1], coords[2]+offset[2],
                display_labels.get(label, label),
                fontsize=20, color='black',
                zorder=111, ha='center', va='center')
    ax.scatter(*centroid_cart, c='gold', marker='*', s=400,
               edgecolors='k', zorder=112, label="Vol. Centroid")
    ax.legend(loc='upper right')


def plot_mapped_bz(ax, points_arr, hull, centroid_cart, unique_ops):
    colormap = plt.colormaps["nipy_spectral"]
    num_ops = len(unique_ops)
    for i, R in enumerate(unique_ops):
        mapped_pts = (R @ points_arr.T).T
        ax.plot_trisurf(mapped_pts[:,0], mapped_pts[:,1], mapped_pts[:,2],
                        triangles=hull.simplices, color=colormap(i/num_ops),
                        edgecolor='none', alpha=0.2, shade=False)
        mc = R @ centroid_cart
        ax.scatter(mc[0], mc[1], mc[2], c='gold', marker='*', s=250,
                   edgecolors='k', zorder=200,
                   label="Mapped Centroids" if i == 0 else None, depthshade=False)
        ax.text(mc[0], mc[1], mc[2], f"  {i+1}", fontsize=10, fontweight='bold', zorder=201)
    avg_pt = np.mean(points_arr, axis=0)
    ax.scatter(*avg_pt, c='cyan', marker='D', s=100, edgecolors='k', zorder=200, label='Avg Point')
    ax.legend(loc='upper right')


# ============================================================================
# Main Pipeline
# ============================================================================
def run(filename, output_dir=None, show_plot=True):
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(filename))
    basename = "BZ"

    print("=" * 60)
    print(f"Processing: {filename}")
    print("=" * 60)

    struct = Structure.from_file(filename)
    a_matrix = struct.lattice.matrix
    cell = a_matrix.tolist()
    positions = struct.frac_coords.tolist()
    numbers = [s.Z for s in struct.species]

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

    spg_cell = (
        np.array(sp_result['primitive_lattice']),
        np.array(sp_result['primitive_positions']),
        sp_result['primitive_types'],
    )
    dataset = spglib.get_symmetry_dataset(spg_cell)
    b_matrix = np.array(sp_result['reciprocal_primitive_lattice'])
    b1, b2, b3 = b_matrix

    sg = dataset.number
    print(f"\nSpace Group: {sg} ({dataset.international})")
    print(f"Point Group: {dataset.pointgroup}")
    print(f"Seekpath Bravais: {sp_result['bravais_lattice_extended']}")

    # ---- Map to Setyawan-Curtarolo BZ type ----
    sc_type, conv_params = seekpath_to_sc_type(sp_result)

    # Override for lower-symmetry point groups with doubled IBZ
    # 6/m, -6, 6 (SG 168-176): 60° sector instead of 30°
    # 4/m, -4, 4 (SG 75-88, P-centered): quadrant instead of octant
    if sc_type == 'HEX' and 168 <= sg <= 176:
        sc_type = 'HEX2'
    elif sc_type == 'CUB' and 195 <= sg <= 206:
        sc_type = 'CUB2'
    elif sc_type == 'FCC' and 195 <= sg <= 206:
        sc_type = 'FCC2'
    elif sc_type == 'BCC' and 195 <= sg <= 206:
        sc_type = 'BCC2'
    elif sc_type == 'TET' and 75 <= sg <= 88:
        sc_type = 'TET2'

    # Map internal HPKOT labels to Setyawan-Curtarolo labels for display
    _hpkot_to_sc = {
        'MCLC1': 'MCLC1', 'MCLC2_SC': 'MCLC2',
        'MCLC2': 'MCLC3', 'MCLC4_SC': 'MCLC4',
        'MCLC3': 'MCLC5',
    }
    sc_display = _hpkot_to_sc.get(sc_type, sc_type)
    print(f"Setyawan-Curtarolo type: {sc_display}")

    # ---- Get IRBZ k-points ----
    # Monoclinic uses the HPKOT tables in lattice_kpoints.py.
    kpoints_frac = get_kpoints(sc_type,
                               conv_params['a'], conv_params.get('b'),
                               conv_params.get('c'), conv_params.get('alpha'))
    kpath = get_kpath(sc_type)
    display_labels = get_display_labels(sc_type)
    params = get_params(sc_type,
                        conv_params['a'], conv_params.get('b'),
                        conv_params.get('c'), conv_params.get('alpha'))
    if params:
        print(f"Parameters: {', '.join(f'{k}={v:.6f}' for k, v in params.items())}")

    print(f"\nHigh-symmetry k-points ({len(kpoints_frac)}):")
    for label, coords in kpoints_frac.items():
        print(f"  {label:8s}: [{coords[0]:8.4f}, {coords[1]:8.4f}, {coords[2]:8.4f}]")

    kpoints_cart = {k: v[0]*b1 + v[1]*b2 + v[2]*b3 for k, v in kpoints_frac.items()}

    # ---- Symmetry operations ----
    sym_ops_cart, unique_ops = get_symmetry_operations(b_matrix, dataset)
    print(f"\nSymmetry operations: {len(sym_ops_cart)}")
    print(f"With time-reversal: {len(unique_ops)}")

    # ---- Convex Hull & Centroid ----
    labels_list = list(kpoints_cart.keys())
    points_arr = np.array([kpoints_cart[k] for k in labels_list])
    hull = ConvexHull(points_arr)
    centroid_cart, ibz_vol = calculate_volume_centroid(hull)
    centroid_frac = centroid_cart @ np.linalg.inv(b_matrix)

    print(f"\n{'='*50}")
    print("NUMERICAL VOLUME CENTROID")
    print(f"{'='*50}")
    print(f"Cartesian:  [{centroid_cart[0]:.6f}, {centroid_cart[1]:.6f}, {centroid_cart[2]:.6f}]")
    print(f"Fractional: [{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]")
    print(f"IBZ Volume: {ibz_vol:.6e}")

    # ---- Symbolic Centroid ----
    print(f"\n{'='*50}")
    print("SYMBOLIC VOLUME CENTROID")
    print(f"{'='*50}")
    try:
        sym_centroid, param_syms = compute_symbolic_centroid(
            kpoints_frac, hull, labels_list, sc_type, conv_params)
        if sym_centroid is not None:
            for i, ax_name in enumerate(['k1', 'k2', 'k3']):
                print(f"  {ax_name} = {sym_centroid[i]}")
            if param_syms:
                subs = [(param_syms[k], params[k]) for k in param_syms if k in params]
                verify = [float(sym_centroid[i].subs(subs)) for i in range(3)]
            else:
                verify = [float(sym_centroid[i]) for i in range(3)]
            print(f"\nVerification:")
            print(f"  Symbolic:   [{verify[0]:.6f}, {verify[1]:.6f}, {verify[2]:.6f}]")
            print(f"  Numerical:  [{centroid_frac[0]:.6f}, {centroid_frac[1]:.6f}, {centroid_frac[2]:.6f}]")
    except Exception as e:
        print(f"  Symbolic computation failed: {e}")
    print(f"{'='*50}")

    # ---- Plotting ----
    bz_loops = get_bz_loops(b_matrix)
    all_bz_pts = np.vstack(bz_loops)
    bz_center = np.mean(all_bz_pts, axis=0)
    bz_span = np.max(all_bz_pts) - np.min(all_bz_pts)

    if show_plot:
        # ---- Interactive plots (all-solid BZ edges, rotate freely) ----
        fig1, ax1 = setup_3d_ax(f"IBZ + BZ: {basename} ({sc_display})",
                                bz_loops, b_matrix, bz_center, bz_span)
        plot_ibz(ax1, kpoints_cart, kpath, display_labels, hull, centroid_cart)
        plt.tight_layout()

        fig2, ax2 = setup_3d_ax(f"Mapped BZ: {basename} — {len(unique_ops)} ops",
                                bz_loops, b_matrix, bz_center, bz_span)
        plot_mapped_bz(ax2, points_arr, hull, centroid_cart, unique_ops)
        plt.tight_layout()

        plt.show()

        # ---- Capture view angles & re-render with dashed back-edges ----
        elev1, azim1 = ax1.elev, ax1.azim
        elev2, azim2 = ax2.elev, ax2.azim
        print(f"\nCaptured view angles:")
        print(f"  Fig1 (IBZ):    elev={elev1:.1f}, azim={azim1:.1f}")
        print(f"  Fig2 (Mapped): elev={elev2:.1f}, azim={azim2:.1f}")

        # Re-render Fig1 with dashed back-edges
        fig1s, ax1s = setup_3d_ax(f"IBZ + BZ: {basename} ({sc_display})",
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  elev=elev1, azim=azim1, dashed_back=True)
        plot_ibz(ax1s, kpoints_cart, kpath, display_labels, hull, centroid_cart)
        plt.tight_layout()
        fig1_path = os.path.join(output_dir, f'{basename}_ibz_{sc_type}.png')
        plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {fig1_path}")
        plt.close(fig1s)

        # Re-render Fig2 with dashed back-edges
        fig2s, ax2s = setup_3d_ax(f"Mapped BZ: {basename} — {len(unique_ops)} ops",
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  elev=elev2, azim=azim2, dashed_back=True)
        plot_mapped_bz(ax2s, points_arr, hull, centroid_cart, unique_ops)
        plt.tight_layout()
        fig2_path = os.path.join(output_dir, f'{basename}_mapped_{sc_type}.png')
        plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {fig2_path}")
        plt.close(fig2s)
    else:
        # ---- Save plots directly without opening a window ----
        plt.switch_backend('Agg')
        fig1s, ax1s = setup_3d_ax(f"IBZ + BZ: {basename} ({sc_display})",
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  dashed_back=True)
        plot_ibz(ax1s, kpoints_cart, kpath, display_labels, hull, centroid_cart)
        plt.tight_layout()
        fig1_path = os.path.join(output_dir, f'{basename}_ibz_{sc_type}.png')
        plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {fig1_path}")
        plt.close(fig1s)

        fig2s, ax2s = setup_3d_ax(f"Mapped BZ: {basename} — {len(unique_ops)} ops",
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  dashed_back=True)
        plot_mapped_bz(ax2s, points_arr, hull, centroid_cart, unique_ops)
        plt.tight_layout()
        fig2_path = os.path.join(output_dir, f'{basename}_mapped_{sc_type}.png')
        plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {fig2_path}")
        plt.close(fig2s)

    return {
        'sc_type': sc_type,
        'seekpath_bravais': sp_result['bravais_lattice_extended'],
        'spacegroup': sg,
        'sg_symbol': dataset.international,
        'point_group': dataset.pointgroup,
        'kpoints_frac': kpoints_frac,
        'centroid_cart': centroid_cart,
        'centroid_frac': centroid_frac,
        'ibz_volume': ibz_vol,
        'n_symmetry_ops': len(unique_ops),
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python compute_centroid_hybrid.py <structure_file> [output_dir]")
        sys.exit(1)

    structure_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None
    results = run(structure_file, out_dir)
