#!/usr/bin/env python3
"""
IBZ Centroid Calculator (Hybrid: seekpath + Setyawan-Curtarolo)
================================================================
Uses seekpath for:  lattice type detection, cell standardization
Uses our own data:  IRBZ k-point vertices (Setyawan-Curtarolo convention)

This ensures the IRBZ shape is consistent for all space groups within
the same Bravais lattice type --no extra points like H_2 in hexagonal.

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
    module=r"seekpath\.hpkot(\..*)-",
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

NO_ALTERMAGNETISM_LAUE_GROUPS = {'-1', '-3', 'm-3'}


def laue_group_from_point_group(point_group):
    """Return the centrosymmetric Laue group associated with a point group."""
    pg = str(point_group).strip().replace(' ', '')
    pg = pg.replace('−', '-').replace('bar', '-')
    mapping = {
        '1': '-1', '-1': '-1',
        '2': '2/m', 'm': '2/m', '2/m': '2/m',
        '222': 'mmm', 'mm2': 'mmm', 'mmm': 'mmm',
        '4': '4/m', '-4': '4/m', '4/m': '4/m',
        '422': '4/mmm', '4mm': '4/mmm', '-42m': '4/mmm',
        '-4m2': '4/mmm', '4/mmm': '4/mmm',
        '3': '-3', '-3': '-3',
        '32': '-3m', '3m': '-3m', '-3m': '-3m',
        '6': '6/m', '-6': '6/m', '6/m': '6/m',
        '622': '6/mmm', '6mm': '6/mmm', '-6m2': '6/mmm',
        '-62m': '6/mmm', '6/mmm': '6/mmm',
        '23': 'm-3', 'm-3': 'm-3',
        '432': 'm-3m', '-43m': 'm-3m', 'm-3m': 'm-3m',
    }
    return mapping.get(pg)


def no_altermagnetism_reason(point_group=None, spacegroup=None):
    """Return a reason dict when the Laue group cannot support altermagnetism."""
    laue_group = laue_group_from_point_group(point_group) if point_group else None

    if laue_group is None and spacegroup is not None:
        sg = int(spacegroup)
        if 1 <= sg <= 2:
            laue_group = '-1'
        elif 143 <= sg <= 148:
            laue_group = '-3'
        elif 195 <= sg <= 206:
            laue_group = 'm-3'

    if laue_group in NO_ALTERMAGNETISM_LAUE_GROUPS:
        return {'laue_group': laue_group, 'reason': 'No altermagnetism'}
    return None

# ============================================================================
# Seekpath Bravais --Setyawan-Curtarolo BZ variation mapping
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
                  f"({k_alpha:.1f}闁? {k_beta:.1f}闁? {k_gamma:.1f}闁?")
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
    # ORCC1: a < b (Table 82, 闁?= (1 + a闁?b闁?/4)
    # ORCC2: a > b (Table 83, 闁?= (1 + b闁?a闁?/4)
    if bravais in ('oS', 'oC', 'oA', 'oB'):
        conv_params['a'], conv_params['b'] = a, b
        return ('ORCC1' if a < b else 'ORCC2'), conv_params

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
    if label == 'GAMMA':
        return '螕'
    # seekpath/HPKOT often writes subscripted labels as B_1, P_1, etc.;
    # lattice_kpoints.py stores the same S-C vertices as B1, P1, etc.
    return label.replace('_', '')


def _display_label_from_internal(label):
    if label == '螕':
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
    # Build edge --face-normal map
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
    ax.set_title(title, fontsize=20)
    return fig, ax


def plot_ibz(ax, kpoints_cart, kpath, display_labels, hull, centroid_cart):
    points_list = list(kpoints_cart.values())
    # Draw IBZ faces (skip for triclinic where hull is None)
    if hull is not None:
        ibz_faces = [[points_list[s] for s in simplex] for simplex in hull.simplices]
        ax.add_collection3d(Poly3DCollection(
            ibz_faces, facecolor='salmon', alpha=0.20, edgecolor='none'))
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
    if hull is not None:
        ax.scatter(*centroid_cart, c='gold', marker='*', s=400,
                   edgecolors='k', zorder=112, label="Vol. Centroid")
        ax.legend(loc='upper right')


def plot_mapped_bz(ax, points_arr, hull, centroid_cart, unique_ops):
    colormap = plt.colormaps["nipy_spectral"]
    num_ops = len(unique_ops)
    for i, R in enumerate(unique_ops):
        mapped_pts = (R @ points_arr.T).T
        if hull is not None:
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
# IBZ frame edge helper
# ============================================================================
def _get_ibz_frame_edges(hull_pts, hull_simplices):
    """Return only the non-coplanar edges of the IBZ hull as (pt1, pt2) pairs.

    Filters out internal triangulation diagonals within flat faces by checking
    whether adjacent face normals are nearly parallel (|cos 闁煎啿鍓?--0.99).
    """
    from collections import defaultdict
    hull_pts = np.array(hull_pts)
    edge_faces = defaultdict(list)
    face_normals = []

    for i, tri in enumerate(hull_simplices):
        a, b, c = hull_pts[tri[0]], hull_pts[tri[1]], hull_pts[tri[2]]
        n = np.cross(b - a, c - a)
        nn = np.linalg.norm(n)
        face_normals.append(n / nn if nn > 1e-10 else np.zeros(3))
        for e in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])]:
            edge_faces[tuple(sorted(e))].append(i)

    edges = []
    for (i, j), faces in edge_faces.items():
        if len(faces) < 2:
            edges.append((hull_pts[i], hull_pts[j]))
        else:
            cos_a = abs(np.dot(face_normals[faces[0]], face_normals[faces[1]]))
            if cos_a < 0.99:
                edges.append((hull_pts[i], hull_pts[j]))
    return edges


# ============================================================================
# Spin-flip Connection Figure  (replaces the old "mapped BZ" Fig 2)
# ============================================================================
def plot_spin_flip_figure(b_matrix, bz_loops, bz_center, bz_span,
                          kpoints_data, ibz_kpoints_frac,
                          hull_pts, hull_simplices,
                          centroid_frac, R,
                          output_path, elev=25, azim=-55, show_plot=True,
                          block=True, path_sequence=None, R_cart=None,
                          defer_show=False):
    """
    Generate the spin-flip connection figure (replaces Fig 2 / mapped-BZ figure).

    Shows the 3D BZ with:
      - Salmon shading   : spin-up IBZ (original)
      - Blue shading     : spin-down IBZ (R-mapped)
      - Red solid lines  : spin-up segments of the generated path (non-primed / k)
      - Blue solid lines : spin-down segments of the generated path (primed / k')
      - Gold star / label: k  (IBZ centroid) with dashed lines to original high-sym points
      - Blue circle      : k' (spin-flip partner) with dashed lines to mapped high-sym points
      - Pink plane       : mirror plane (only when R is a pure mirror)

    path_sequence : list returned by KPointsModifier.insert_general_kpoints(), optional.
        Each entry is [kx, ky, kz, label] (fractional coords) or None (segment break).
        When provided, only the generated-path segments are drawn, coloured by spin side.
        When None, falls back to drawing all raw KPOINTS file segments in red.
    """
    b1, b2, b3 = b_matrix[0], b_matrix[1], b_matrix[2]
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)   # maps Cartesian k -> fractional
    R_inv_T = np.linalg.inv(R).T
    if R_cart is None:
        R_cart = b_T @ R_inv_T @ b_T_inv
    else:
        R_cart = np.array(R_cart, dtype=float)

    # k and k' in Cartesian
    k_frac = np.array(centroid_frac[:3])
    k_cart = k_frac[0] * b1 + k_frac[1] * b2 + k_frac[2] * b3
    kp_cart = R_cart @ k_cart
    kp_frac = b_T_inv @ kp_cart

    # Original IBZ high-sym points (from lattice_kpoints, not KPOINTS file)
    ibz_orig = {}
    for lbl, frac in ibz_kpoints_frac.items():
        if not lbl.startswith('_'):
            f = np.array(frac)
            ibz_orig[lbl] = f[0] * b1 + f[1] * b2 + f[2] * b3

    # Spin-down IBZ high-sym points: apply R^{-T} to each original point.
    # ibz_mapped      : for labels only --skip points that coincide with an
    #                   original high-sym point (e.g. 闁剧粯娲熼崺?- to avoid stacked labels.
    # ibz_mapped_lines: for dashed-line drawing --include ALL mapped points so
    #                   that the k'闂?-line (and any other self-mapped point) is
    #                   still drawn even when it carries no new label.
    ibz_mapped       = {}
    ibz_mapped_lines = {}
    for lbl, frac in ibz_kpoints_frac.items():
        if lbl.startswith('_'):
            continue
        f = np.array(frac)
        pt0 = f[0] * b1 + f[1] * b2 + f[2] * b3
        pt = R_cart @ pt0
        ibz_mapped_lines[lbl + "'"] = pt
        coincident = any(
            np.linalg.norm(pt - orig_pt) < 1e-6
            for orig_pt in ibz_orig.values()
        )
        if not coincident:
            ibz_mapped[lbl + "'"] = pt

    # Spin-down IBZ hull vertices (apply same fractional transform to each vertex)
    mapped_hull_pts = None
    if hull_pts is not None:
        mapped_hull_pts = (R_cart @ np.array(hull_pts).T).T

    # Filter threshold: skip connections shorter than 5% of BZ radius
    bz_radius = np.max(np.linalg.norm(np.vstack(bz_loops), axis=1))
    threshold = 0.05 * bz_radius

    def _disp(lbl):
        prime = lbl.endswith("'")
        base = lbl.rstrip("'")
        suffix = r"$'$" if prime else ''
        if base.upper() == 'GAMMA':
            return r"$\Gamma$" + (r"$'$" if prime else '')
        if '_' in base:
            b, s = base.split('_', 1)
            return rf"${b}_{{{s}}}${suffix}"
        return f"${base}${suffix}" if prime else base

    def _draw(ax):
        # Spin-up IBZ: salmon shading (no triangulation diagonals)
        if hull_pts is not None and hull_simplices is not None:
            ax.plot_trisurf(hull_pts[:, 0], hull_pts[:, 1], hull_pts[:, 2],
                            triangles=hull_simplices, color='salmon',
                            edgecolor='none', alpha=0.20, shade=False)
            for p1, p2 in _get_ibz_frame_edges(hull_pts, hull_simplices):
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                        c='darkred', lw=1.8, alpha=0.85, zorder=10)

        # Spin-down IBZ: blue shading (no triangulation diagonals)
        if mapped_hull_pts is not None and hull_simplices is not None:
            ax.plot_trisurf(mapped_hull_pts[:, 0], mapped_hull_pts[:, 1],
                            mapped_hull_pts[:, 2],
                            triangles=hull_simplices, color='cornflowerblue',
                            edgecolor='none', alpha=0.20, shade=False)
            for p1, p2 in _get_ibz_frame_edges(mapped_hull_pts, hull_simplices):
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                        c='navy', lw=1.8, alpha=0.85, zorder=10)

        # Generated path segments (red = spin-up side, navy = spin-down side).
        # Use path_sequence for labels/order only. The coordinates are looked up from
        # ibz_orig/ibz_mapped_lines so Figure 2 uses the same high-symmetry convention
        # as Figure 1, even when seekpath and S-C labels share names but differ by basis.
        def _path_point(label):
            if label in ('k', "k'"):
                return None
            is_prime = label.endswith("'")
            base = _seekpath_label_to_internal(label.rstrip("'"))
            if is_prime:
                return ibz_mapped_lines.get(base + "'")
            return ibz_orig.get(base)

        if path_sequence is not None:
            for i in range(len(path_sequence) - 1):
                A = path_sequence[i]
                B = path_sequence[i + 1]
                if A is None or B is None:
                    continue
                la, lb = A[3], B[3]
                # Skip only k<->k' spin-flip jumps and k endpoint helper segments.
                # Mixed primed/unprimed high-symmetry transitions are real line-mode
                # KPOINTS segments and must be drawn.
                if la in ('k', "k'") or lb in ('k', "k'"):
                    continue
                a_prime = la.endswith("'")
                b_prime = lb.endswith("'")
                pa_c = _path_point(la)
                pb_c = _path_point(lb)
                if pa_c is None or pb_c is None:
                    continue
                col = 'navy' if (a_prime and b_prime) else ('red' if not (a_prime or b_prime) else 'purple')
                ax.plot([pa_c[0], pb_c[0]], [pa_c[1], pb_c[1]], [pa_c[2], pb_c[2]],
                        c=col, lw=4.0, alpha=0.9, zorder=50)
        else:
            # Fallback: draw all raw KPOINTS file segments in red
            for i in range(0, len(kpoints_data) - 1, 2):
                p1c = (kpoints_data[i][0] * b1 + kpoints_data[i][1] * b2
                       + kpoints_data[i][2] * b3)
                p2c = (kpoints_data[i+1][0] * b1 + kpoints_data[i+1][1] * b2
                       + kpoints_data[i+1][2] * b3)
                ax.plot([p1c[0], p2c[0]], [p1c[1], p2c[1]], [p1c[2], p2c[2]],
                        c='red', lw=2.5, alpha=0.9, zorder=50)

        # Label helpers
        def _label_pts(pts_dict, color, edgecolor):
            arr = np.array(list(pts_dict.values()))
            center = np.mean(arr, axis=0) if len(arr) else np.zeros(3)
            span = max(np.max(np.ptp(arr, axis=0)), 1e-8) if len(arr) else 1.0
            off_sc = span * 0.18
            for lbl, hpt in pts_dict.items():
                ax.scatter(*hpt, c=color, s=60, zorder=110,
                           edgecolors=edgecolor, linewidths=0.5)
                direction = hpt - center
                nd = np.linalg.norm(direction)
                off = direction / nd * off_sc if nd > 1e-8 else np.array([0, 0, off_sc])
                ax.text(*(hpt + off), _disp(lbl), fontsize=20, color=edgecolor,
                        zorder=111, ha='center', va='center')

        _label_pts(ibz_orig,   color='salmon',          edgecolor='darkred')
        _label_pts(ibz_mapped, color='cornflowerblue',  edgecolor='navy')

        # k --original high-sym points (dashed blue)
        for hpt in ibz_orig.values():
            if np.linalg.norm(hpt - k_cart) > threshold:
                ax.plot([hpt[0], k_cart[0]], [hpt[1], k_cart[1]], [hpt[2], k_cart[2]],
                        c='deepskyblue', lw=2.0, ls='--', alpha=0.75, zorder=40)

        # k' --spin-down high-sym points (dashed blue)
        # Use ibz_mapped_lines (includes self-mapped 闁? so k'闂?-is always drawn.
        for hpt in ibz_mapped_lines.values():
            if np.linalg.norm(hpt - kp_cart) > threshold:
                ax.plot([hpt[0], kp_cart[0]], [hpt[1], kp_cart[1]], [hpt[2], kp_cart[2]],
                        c='deepskyblue', lw=2.0, ls='--', alpha=0.75, zorder=40)

        # k and k' markers
        ax.scatter(*k_cart, c='gold', s=300, marker='*',
                   edgecolors='k', linewidths=0.8, zorder=120, label=r'$k$')
        ax.scatter(*kp_cart, c='cornflowerblue', s=150, marker='o',
                   edgecolors='k', linewidths=0.8, zorder=120, label=r"$k'$")
        ax.legend(loc='upper right', fontsize=16)

    fig, ax = setup_3d_ax("Spin-flip path connections",
                          bz_loops, b_matrix, bz_center, bz_span,
                          elev=elev, azim=azim, dashed_back=False)
    _draw(ax)
    plt.tight_layout()

    display_fig = fig if show_plot and defer_show else None
    if display_fig is not None:
        def _save_after_show(fig=fig, ax=ax):
            fig_save, ax_save = setup_3d_ax("Spin-flip path connections",
                                            bz_loops, b_matrix, bz_center, bz_span,
                                            elev=ax.elev, azim=ax.azim,
                                            dashed_back=True)
            _draw(ax_save)
            plt.tight_layout()
            fig_save.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig_save)
            plt.close(fig)
            print(f"Saved: {output_path}")
        display_fig._alterseek_save_after_show = _save_after_show
        return display_fig

    if show_plot and not defer_show:
        plt.show(block=block)
        if block:
            # Blocking mode: window already closed, capture adjusted view
            elev, azim = ax.elev, ax.azim
            plt.close(fig)
        # Non-blocking mode: window stays open alongside the next figure;
        # save uses the original elev/azim (view adjustments not captured).
    elif not show_plot:
        plt.close(fig)

    # Re-render with dashed back-edges for saved PNG
    fig_save, ax_save = setup_3d_ax("Spin-flip path connections",
                                    bz_loops, b_matrix, bz_center, bz_span,
                                    elev=elev, azim=azim, dashed_back=True)
    _draw(ax_save)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig_save)
    print(f"Saved: {output_path}")
    return display_fig


# ============================================================================
# Spin-up / Spin-down Full-BZ Coloring Figure  (replaces old rainbow Fig 2)
# ============================================================================
def plot_spin_bz_figure(b_matrix, bz_loops, bz_center, bz_span,
                        unique_ops, centroid_cart,
                        hull_pts, hull_simplices,
                        R, output_path,
                        flip_ops_frac=None,
                        elev=25, azim=-55, show_plot=True,
                        defer_show=False):
    """
    Show the full BZ colored by spin channel (replaces old rainbow Fig 2).

    Each orbit of the IBZ is colored by proximity: an op g is "spin-down" if
    g @ centroid_cart lands closer to the spin-down centroid (R^{-T} @ centroid)
    than to the spin-up centroid. This works correctly for all lattice types,
    including high-symmetry groups (e.g. Oh) where many ops rotate the spin
    axis to a third direction rather than cleanly flipping it.
      - RED   (salmon)         : spin-up  regions
      - BLUE  (cornflowerblue) : spin-down regions
    """
    if hull_pts is None or hull_simplices is None or not len(unique_ops):
        print("[Note] Skipping spin-BZ figure (no hull or symmetry ops available).")
        return

    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)

    hull_pts = np.array(hull_pts)
    centroid_cart = np.array(centroid_cart)
    hull_simplices_arr = np.array(hull_simplices)

    # Classify spin channel for each op in unique_ops.
    #
    # Primary method (when flip_ops_frac is available):
    #   Convert each g_cart --g_frac (primitive) and check if it matches any
    #   entry in flip_ops_frac.  This is exact and R-independent: the coloring
    #   is the same regardless of which spin-flip R the user picked.
    #
    # Fallback (proximity, R-dependent --kept only as safety net):
    #   An op is "spin-down" if g @ centroid is closer to kp = R^{-T}@centroid
    #   than to centroid.  This can give different results for different R choices
    #   because kp depends on R.
    if flip_ops_frac is not None and len(flip_ops_frac):
        flip_set = [np.array(f, dtype=float) for f in flip_ops_frac]
        spin_down_mask = np.zeros(len(unique_ops), dtype=bool)
        for i, g_cart in enumerate(unique_ops):
            # g_cart = b_T @ inv(g_frac).T @ b_T_inv  -- g_frac = inv((b_T_inv@g_cart@b_T).T)
            M = b_T_inv @ g_cart @ b_T          # = inv(g_frac).T
            g_frac = np.linalg.inv(M.T)
            spin_down_mask[i] = any(np.allclose(g_frac, f, atol=1e-6) for f in flip_set)
    else:
        # Fallback: proximity against a single kp reference (R-dependent)
        R_inv_T = np.linalg.inv(np.array(R)).T
        R_cart  = b_T @ R_inv_T @ b_T_inv
        kp_cart = R_cart @ centroid_cart

        def _proximity_mask(c_pt, kp_pt):
            return np.array([
                np.linalg.norm(g @ c_pt - kp_pt) < np.linalg.norm(g @ c_pt - c_pt)
                for g in unique_ops
            ])

        spin_down_mask = _proximity_mask(centroid_cart, kp_cart)
        n_expected = len(unique_ops) // 2
        if spin_down_mask.sum() != n_expected:
            eps_scale = np.linalg.norm(centroid_cart) * 3e-4
            for trial in range(30):
                rng = np.random.default_rng(trial)
                c_pert  = centroid_cart + rng.standard_normal(3) * eps_scale
                kp_pert = R_cart @ c_pert
                mask_try = _proximity_mask(c_pert, kp_pert)
                if mask_try.sum() == n_expected:
                    spin_down_mask = mask_try
                    break

    def _draw(ax):
        for i, g in enumerate(unique_ops):
            mapped_pts = (g @ hull_pts.T).T
            if spin_down_mask[i]:
                color, alpha = 'cornflowerblue', 0.2
            else:
                color, alpha = 'salmon', 0.2
            ax.plot_trisurf(mapped_pts[:, 0], mapped_pts[:, 1], mapped_pts[:, 2],
                            triangles=hull_simplices_arr, color=color,
                            edgecolor='none', alpha=alpha, shade=False)
        # Gold star at the original IBZ centroid
        ax.scatter(*centroid_cart, c='gold', s=350, marker='*',
                   edgecolors='k', linewidths=0.8, zorder=200,
                   label=r'$k$ (IBZ centroid)', depthshade=False)
        ax.legend(loc='upper right', fontsize=10)

    fig, ax = setup_3d_ax("Spin-up (red) / Spin-down (blue) BZ",
                          bz_loops, b_matrix, bz_center, bz_span,
                          elev=elev, azim=azim, dashed_back=False)
    _draw(ax)
    plt.tight_layout()

    display_fig = fig if show_plot and defer_show else None
    if display_fig is not None:
        def _save_after_show(fig=fig, ax=ax):
            fig_save, ax_save = setup_3d_ax("Spin-up (red) / Spin-down (blue) BZ",
                                            bz_loops, b_matrix, bz_center, bz_span,
                                            elev=ax.elev, azim=ax.azim,
                                            dashed_back=True)
            _draw(ax_save)
            plt.tight_layout()
            fig_save.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig_save)
            plt.close(fig)
            print(f"Saved: {output_path}")
        display_fig._alterseek_save_after_show = _save_after_show
        return display_fig

    if show_plot and not defer_show:
        plt.show()
        elev, azim = ax.elev, ax.azim
        plt.close(fig)
    elif not show_plot:
        plt.close(fig)

    # Re-render with dashed back-edges for saved PNG
    fig, ax = setup_3d_ax("Spin-up (red) / Spin-down (blue) BZ",
                          bz_loops, b_matrix, bz_center, bz_span,
                          elev=elev, azim=azim, dashed_back=True)
    _draw(ax)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")
    return display_fig



def _classify_spin_down_ops(b_matrix, unique_ops, centroid_cart, R, flip_ops_frac=None):
    """Return a boolean mask selecting the spin-down symmetry images."""
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)
    centroid_cart = np.array(centroid_cart)

    if flip_ops_frac is not None and len(flip_ops_frac):
        flip_set = [np.array(f, dtype=float) for f in flip_ops_frac]
        spin_down_mask = np.zeros(len(unique_ops), dtype=bool)
        for i, g_cart in enumerate(unique_ops):
            M = b_T_inv @ g_cart @ b_T
            g_frac = np.linalg.inv(M.T)
            spin_down_mask[i] = any(np.allclose(g_frac, f, atol=1e-6) for f in flip_set)
        return spin_down_mask

    R_inv_T = np.linalg.inv(np.array(R)).T
    R_cart  = b_T @ R_inv_T @ b_T_inv
    kp_cart = R_cart @ centroid_cart

    def _proximity_mask(c_pt, kp_pt):
        return np.array([
            np.linalg.norm(g @ c_pt - kp_pt) < np.linalg.norm(g @ c_pt - c_pt)
            for g in unique_ops
        ])

    spin_down_mask = _proximity_mask(centroid_cart, kp_cart)
    n_expected = len(unique_ops) // 2
    if spin_down_mask.sum() != n_expected:
        eps_scale = np.linalg.norm(centroid_cart) * 3e-4
        for trial in range(30):
            rng = np.random.default_rng(trial)
            c_pert  = centroid_cart + rng.standard_normal(3) * eps_scale
            kp_pert = R_cart @ c_pert
            mask_try = _proximity_mask(c_pert, kp_pert)
            if mask_try.sum() == n_expected:
                spin_down_mask = mask_try
                break
    return spin_down_mask


def _points_on_kz_plane(points, simplices, z0=0.0, tol=1e-8):
    """Return the 2D convex section of a triangular hull with the kz=z0 plane."""
    section = []
    points = np.array(points, dtype=float)

    for tri in np.array(simplices, dtype=int):
        verts = points[tri]
        for a, b in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            da = a[2] - z0
            db = b[2] - z0
            if abs(da) <= tol:
                section.append(a[:2])
            if abs(db) <= tol:
                section.append(b[:2])
            if da * db < -tol * tol:
                t = da / (da - db)
                section.append((a + t * (b - a))[:2])

    if len(section) < 3:
        return None

    pts = np.unique(np.round(np.array(section), 10), axis=0)
    if len(pts) < 3:
        return None

    try:
        hull = ConvexHull(pts)
        return pts[hull.vertices]
    except Exception:
        center = pts.mean(axis=0)
        angles = np.arctan2(pts[:, 1] - center[1], pts[:, 0] - center[0])
        ordered = pts[np.argsort(angles)]
        return ordered if len(ordered) >= 3 else None


def _bz_kz_plane_outline(bz_loops, z0=0.0, tol=1e-8):
    """Return the top-view outline where the BZ boundary cuts kz=z0."""
    section = []
    for loop in bz_loops:
        pts = np.array(loop, dtype=float)
        for a, b in zip(pts[:-1], pts[1:]):
            da = a[2] - z0
            db = b[2] - z0
            if abs(da) <= tol:
                section.append(a[:2])
            if abs(db) <= tol:
                section.append(b[:2])
            if da * db < -tol * tol:
                t = da / (da - db)
                section.append((a + t * (b - a))[:2])

    if len(section) < 3:
        return None

    pts = np.unique(np.round(np.array(section), 10), axis=0)
    if len(pts) < 3:
        return None

    try:
        hull = ConvexHull(pts)
        return pts[hull.vertices]
    except Exception:
        center = pts.mean(axis=0)
        angles = np.arctan2(pts[:, 1] - center[1], pts[:, 0] - center[0])
        return pts[np.argsort(angles)]


def plot_spin_bz_top_view_figure(b_matrix, bz_loops,
                                 unique_ops, centroid_cart,
                                 hull_pts, hull_simplices,
                                 R, output_path,
                                 flip_ops_frac=None,
                                 show_plot=True, defer_show=False,
                                 z0=0.0):
    """Draw a darker top-view kz=0 slice of the Figure 3 spin-BZ coloring."""
    if hull_pts is None or hull_simplices is None or not len(unique_ops):
        print("[Note] Skipping spin-BZ top-view figure (no hull or symmetry ops available).")
        return

    hull_pts = np.array(hull_pts, dtype=float)
    hull_simplices_arr = np.array(hull_simplices, dtype=int)
    centroid_cart = np.array(centroid_cart, dtype=float)
    spin_down_mask = _classify_spin_down_ops(
        b_matrix, unique_ops, centroid_cart, R, flip_ops_frac)

    z_span = np.ptp(np.vstack(bz_loops)[:, 2])
    z_eps = max(float(z_span) * 1e-6, 1e-8)
    side = np.sign(centroid_cart[2] - z0)
    if abs(side) < 1e-12:
        side = 1.0
    section_z = z0 + side * z_eps

    fig, ax = plt.subplots(figsize=(9, 9))
    up_labeled = False
    down_labeled = False

    for i, g in enumerate(unique_ops):
        mapped_pts = (g @ hull_pts.T).T
        poly = _points_on_kz_plane(mapped_pts, hull_simplices_arr, z0=section_z)
        if poly is None:
            continue

        is_down = spin_down_mask[i]
        if is_down:
            color = '#1f4e9e'
            label = 'spin-down' if not down_labeled else None
            down_labeled = True
        else:
            color = '#b22222'
            label = 'spin-up' if not up_labeled else None
            up_labeled = True

        closed = np.vstack([poly, poly[0]])
        ax.fill(poly[:, 0], poly[:, 1], facecolor=color, alpha=0.68,
                edgecolor='none', label=label)
        ax.plot(closed[:, 0], closed[:, 1], color=color, lw=0.9, alpha=0.95)

    outline = _bz_kz_plane_outline(bz_loops, z0=z0)
    if outline is not None:
        closed = np.vstack([outline, outline[0]])
        ax.plot(closed[:, 0], closed[:, 1], color='black', lw=2.0, label='BZ boundary')

    ax.set_aspect('equal', adjustable='box')
    ax.set_title(r'Spin-up / Spin-down BZ top view ($k_z = 0$)', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    ax.legend(loc='upper right', fontsize=10)
    fig.tight_layout()

    display_fig = fig if show_plot and defer_show else None
    if display_fig is not None:
        def _save_after_show(fig=fig):
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_path}")
        display_fig._alterseek_save_after_show = _save_after_show
        return display_fig

    if show_plot and not defer_show:
        plt.show()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close(fig)
    return display_fig
# ============================================================================
# Main Pipeline
# ============================================================================
def run(filename, output_dir=None, show_plot=True, defer_show=False, verbose=True):
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
    # Reciprocal lattice of the user-provided input cell. spinspg rotations
    # are written in this basis by find_sf_operations.py.
    b_matrix_input = 2 * np.pi * np.linalg.inv(np.array(a_matrix)).T

    # Conventional-cell reciprocal lattice (no 2闁?needed --cancels in formula).
    # Used to correctly convert seed flip ops that were written using the
    # conventional cell (ASE-read POSCAR --spinspg --flip_spin_operations.txt).
    _conv_lat = np.array(sp_result.get('conv_lattice', sp_result['primitive_lattice']))
    b_matrix_conv = np.linalg.inv(_conv_lat).T
    b1, b2, b3 = b_matrix

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

    # ---- Map to Setyawan-Curtarolo BZ type ----
    sc_type, conv_params = seekpath_to_sc_type(sp_result)

    # Override for lower-symmetry point groups with enlarged IBZ
    # Hexagonal 6/m, -6, 6 (SG 168-176): 60闁?sector (2闁?HEX)
    # Trigonal 32, 3m, -3m on P lattice (SG 149-167): Laue=-3m (12 ops) --2闁?HEX = HEX2
    # Trigonal 3, -3 on P lattice (SG 143-148): Laue=-3 (6 ops) --4闁?HEX = HEX4
    # 4/m, -4, 4 (SG 75-88, P-centered): quadrant instead of octant
    # 4/m, -4, 4 (SG 75-88, I-centered): doubled BCT IBZ (闁挎粎绂哾(110) absent in C4h)
    if sc_type in ('RHL1', 'RHL2') and 143 <= sg <= 148:
        sc_type = sc_type + '_2'
    elif sc_type == 'HEX' and 168 <= sg <= 176:
        sc_type = 'HEX2'
    elif sc_type == 'HEX' and 149 <= sg <= 167:
        sc_type = 'HEX2'
    elif sc_type == 'HEX' and 143 <= sg <= 148:
        sc_type = 'HEX4'
    elif sc_type == 'CUB' and 195 <= sg <= 206:
        sc_type = 'CUB2'
    elif sc_type == 'FCC' and 195 <= sg <= 206:
        sc_type = 'FCC2'
    elif sc_type == 'BCC' and 195 <= sg <= 206:
        sc_type = 'BCC2'
    elif sc_type == 'TET' and 75 <= sg <= 88:
        sc_type = 'TET2'
    elif sc_type == 'BCT1' and 75 <= sg <= 88:
        sc_type = 'BCT1_2'
    elif sc_type == 'BCT2' and 75 <= sg <= 88:
        sc_type = 'BCT2_2'

    # Map internal HPKOT labels to Setyawan-Curtarolo labels for display
    _hpkot_to_sc = {
        'MCLC1': 'MCLC1', 'MCLC2_SC': 'MCLC2',
        'MCLC2': 'MCLC3', 'MCLC4_SC': 'MCLC4',
        'MCLC3': 'MCLC5',
    }
    sc_display = _hpkot_to_sc.get(sc_type, sc_type)

    # ---- Get IRBZ k-points ----
    # Use the real (possibly enlarged) IBZ for centroid computation.
    # For primitive trigonal -3m/3m/32 the enlarged HEX2 IBZ must preserve
    # the standard hexagonal irreducible triangle Gamma-M-K and add the
    # neighboring sector; do not rotate it to the Gamma-K-K2 sector.
    centroid_type = sc_type

    if sg in (1, 2):
        # Triclinic: use seekpath k-points directly
        # (lattice_kpoints.py TRI types have convention issues)
        if verbose:
            print("[Note] Triclinic --using seekpath k-points for visualization")
        sp_coords = sp_result['point_coords']
        kpoints_frac = {k: list(v) for k, v in sp_coords.items()}
        kpath = list(sp_result['path'])
        display_labels = {k: ('Γ' if k == 'GAMMA' else k) for k in sp_coords}
        sc_display = sp_result['bravais_lattice_extended']
        kpoints_cart = {k: v[0]*b1 + v[1]*b2 + v[2]*b3
                        for k, v in sp_coords.items()}
        kpoints_frac_centroid = kpoints_frac
        kpoints_cart_centroid = dict(kpoints_cart)
    else:
        # Full IBZ k-points for plotting (sc_type, may include extra vertices)
        kpoints_frac = get_kpoints(sc_type,
                                   conv_params['a'], conv_params.get('b'),
                                   conv_params.get('c'), conv_params.get('alpha'))
        kpath = get_kpath(sc_type)
        display_labels = get_display_labels(sc_type)
        params = get_params(sc_type,
                            conv_params['a'], conv_params.get('b'),
                            conv_params.get('c'), conv_params.get('alpha'))
        if params and verbose:
            print(f"Parameters: {', '.join(f'{k}={v:.6f}' for k, v in params.items())}")

        kpoints_cart = {k: v[0]*b1 + v[1]*b2 + v[2]*b3 for k, v in kpoints_frac.items()}

        # IBZ k-points for centroid
        kpoints_frac_centroid = get_kpoints(centroid_type,
                                   conv_params['a'], conv_params.get('b'),
                                   conv_params.get('c'), conv_params.get('alpha'))
        kpoints_cart_centroid = {k: v[0]*b1 + v[1]*b2 + v[2]*b3
                                 for k, v in kpoints_frac_centroid.items()}

    # ---- Symmetry operations ----
    sym_ops_cart, unique_ops = get_symmetry_operations(b_matrix, dataset)
    if verbose:
        print(f"\nSymmetry operations: {len(sym_ops_cart)}")
        print(f"With time-reversal: {len(unique_ops)}")

    # ---- Convex Hull & Centroid (using base IBZ) ----
    labels_list = list(kpoints_cart_centroid.keys())
    points_arr = np.array([kpoints_cart_centroid[k] for k in labels_list])

    if sg in (1, 2):
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

        # ---- Symbolic Centroid (saved to file, not printed) ----
        try:
            sym_centroid, param_syms = compute_symbolic_centroid(
                kpoints_frac_centroid, hull, labels_list, centroid_type, conv_params)
            if sym_centroid is not None:
                sym_lines = "\n".join(
                    f"  {ax_name} = {sym_centroid[i]}"
                    for i, ax_name in enumerate(['k1', 'k2', 'k3'])
                )
                try:
                    with open("spin_operations.txt", "a") as f:
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

    # For plotting: use base type k-points, kpath and display labels
    # (hull was computed from centroid type, so faces must use same point set)
    if sg in (1, 2):
        kpath_plot = kpath
        display_labels_plot = display_labels
    else:
        kpath_plot = get_kpath(centroid_type)
        display_labels_plot = get_display_labels(centroid_type)

    # ---- Plotting ----
    fig1_title = (f"BZ: {basename} ({sc_display})" if sg in (1, 2)
                  else f"IBZ + BZ: {basename} ({sc_display})")

    display_figures = []
    fig1_path = os.path.join(output_dir, f'{basename}_ibz_{sc_type}.png')
    if show_plot:
        # Interactive mode: create the figure now. alterseek_path can defer
        # the actual plt.show() call until all prompts and file writes finish.
        fig1, ax1 = setup_3d_ax(fig1_title,
                                bz_loops, b_matrix, bz_center, bz_span)
        plot_ibz(ax1, kpoints_cart_centroid, kpath_plot, display_labels_plot, hull, centroid_cart)
        plt.tight_layout()
        if defer_show:
            def _save_fig1_after_show(fig=fig1, ax=ax1):
                fig1s, ax1s = setup_3d_ax(fig1_title,
                                          bz_loops, b_matrix, bz_center, bz_span,
                                          elev=ax.elev, azim=ax.azim, dashed_back=True)
                plot_ibz(ax1s, kpoints_cart_centroid, kpath_plot, display_labels_plot, hull, centroid_cart)
                plt.tight_layout()
                fig1s.savefig(fig1_path, dpi=300, bbox_inches='tight')
                plt.close(fig1s)
                plt.close(fig)
                print(f"Saved: {fig1_path}")
            fig1._alterseek_save_after_show = _save_fig1_after_show
            display_figures.append(fig1)
            elev1, azim1 = 25, -55
        else:
            plt.show()
            elev1, azim1 = ax1.elev, ax1.azim
    else:
        # Automated mode (called from alterseek_path): use default angles, no window
        elev1, azim1 = 25, -55

    # Render with dashed back-edges and save unless deferred post-show saving is active
    if not (show_plot and defer_show):
        fig1s, ax1s = setup_3d_ax(fig1_title,
                                  bz_loops, b_matrix, bz_center, bz_span,
                                  elev=elev1, azim=azim1, dashed_back=True)
        plot_ibz(ax1s, kpoints_cart_centroid, kpath_plot, display_labels_plot, hull, centroid_cart)
        plt.tight_layout()
        fig1s.savefig(fig1_path, dpi=300, bbox_inches='tight')
        if verbose:
            print(f"Saved: {fig1_path}")
        plt.close(fig1s)

    return {
        'sc_type': sc_type,
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
        'ibz_kpath': kpath_plot,
        'hull_pts': points_arr if sg not in (1, 2) else None,
        'hull_simplices': hull.simplices.tolist() if (sg not in (1, 2) and hull is not None) else None,
        'sym_ops_cart': sym_ops_cart,
        'unique_ops': unique_ops,
        'b_matrix_conv': b_matrix_conv,
        'b_matrix_input': b_matrix_input,
        'display_figures': display_figures,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python compute_centroid_hybrid.py <structure_file> [output_dir]")
        sys.exit(1)

    structure_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None
    results = run(structure_file, out_dir)

