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
from scipy.spatial import Voronoi, ConvexHull, HalfspaceIntersection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
import sympy as sp
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

    with open(output_path, "w") as f:
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
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


from lattice_kpoints import (
    LATTICE_DATA, get_kpoints, get_hull_kpoints,
    get_hull_kpath, get_display_labels, get_params,
    canonical_lattice_type, _normalize_label,
)

NO_ALTERMAGNETISM_LAUE_GROUPS = {'-1', '-3', 'm-3'}
GAMMA_LABEL = "\u0393"
BZ_SPECIAL_COLORS = {
    "orange": "#e68613",
    "purple": "#6b5596",
}
IBZ_FACE_COLORS = {
    "up_main": "salmon",
    "up_extra": "#e98f8f",
    "down_main": "cornflowerblue",
    "down_extra": "#91b2e8",
}
BZ_PATH_STYLE_OVERRIDES = {
    "cP1": {("M", "X_1"): {"color": "red", "ls": "--"}},
    "cF1": {("X", "W_2"): {"color": "red", "ls": "--"}},
    "tI1": {(GAMMA_LABEL, "N"): "orange"},
    "tI2": {(GAMMA_LABEL, "N"): "orange"},
    "oF1": {(GAMMA_LABEL, "L"): "orange"},
    "oF2": {(GAMMA_LABEL, "L"): "orange"},
    "oF3": {(GAMMA_LABEL, "L"): "orange"},
    "oI1": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oI2": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oI3": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oA1": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oA2": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oC1": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oC2": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "hR1": {
        (GAMMA_LABEL, "L"): "orange",
        (GAMMA_LABEL, "F"): "orange",
    },
    "hR2": {(GAMMA_LABEL, "L"): "orange"},
    "mP1": {
        (GAMMA_LABEL, "B"): "orange",
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "Y_2"): "orange",
        ("Z", "D"): "purple",
        ("Z", "E"): "purple",
        ("Z", "C_2"): "purple",
    },
    "mC1": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "M_2"): "orange",
        (GAMMA_LABEL, "Y_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
    },
    "mC2": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
        ("M", "Y"): "purple",
    },
    "mC3": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "M_2"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
    },
}


def _is_doubled_ibz_extra_label(label):
    label = str(label)
    suffix = label.rsplit("_", 1)[-1] if "_" in label else ""
    return suffix.endswith("A") or label.endswith("_2")


def _doubled_ibz_extra_flags(hull_labels):
    labels = [str(label) for label in hull_labels]
    # Project doubled sectors always include at least one copied _A label.
    # Ordinary HPKOT labels ending in _2 must not trigger sector coloring.
    if not any(
        "_" in label and label.rsplit("_", 1)[-1].endswith("A")
        for label in labels
    ):
        return [False] * len(labels)
    return [_is_doubled_ibz_extra_label(label) for label in labels]


def _split_hull_faces_by_extra_labels(hull_pts, hull_simplices, hull_labels=None):
    points = np.array(hull_pts)
    main_faces = []
    extra_faces = []
    for simplex in hull_simplices:
        face = [points[s] for s in simplex]
        if hull_labels is not None and any(
            _is_doubled_ibz_extra_label(hull_labels[s])
            for s in simplex
            if s < len(hull_labels)
        ):
            extra_faces.append(face)
        else:
            main_faces.append(face)
    return main_faces, extra_faces


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
# SeeK-path/HPKOT extended Bravais key and conventional parameters
# ============================================================================
def seekpath_to_hpkot_type(sp_result):
    """Return the HPKOT extended key and parameters for lattice_kpoints.py."""
    lattice_key = sp_result.get('bravais_lattice_extended',
                                sp_result.get('bravais_lattice', 'cP'))

    conv_lattice = np.array(sp_result.get('conv_lattice',
                            sp_result.get('primitive_lattice')), dtype=float)
    va, vb, vc = conv_lattice[0], conv_lattice[1], conv_lattice[2]
    a = np.linalg.norm(va)
    b = np.linalg.norm(vb)
    c = np.linalg.norm(vc)
    alpha = np.degrees(np.arccos(np.clip(np.dot(vb, vc)/(b*c), -1, 1)))
    beta = np.degrees(np.arccos(np.clip(np.dot(va, vc)/(a*c), -1, 1)))
    gamma = np.degrees(np.arccos(np.clip(np.dot(va, vb)/(a*b), -1, 1)))

    conv_params = {
        'a': a, 'b': b, 'c': c,
        'alpha': alpha, 'beta': beta, 'gamma': gamma,
    }
    if lattice_key.startswith('m'):
        # HPKOT monoclinic table expressions use beta between a and c.
        conv_params['alpha'] = beta
    return lattice_key, conv_params


# Backward-compatible import name.
def seekpath_to_sc_type(sp_result):
    return seekpath_to_hpkot_type(sp_result)


def _seekpath_label_to_internal(label):
    if label == 'GAMMA':
        return '\u0393'
    return label


def _display_label_from_internal(label):
    return _math_label(label)


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


# ---------------------------------------------------------------------------
# 2D / slab helpers
#
# For a 2D material computed as a 3D slab (large vacuum along one axis) the
# whole symmetry analysis is reused unchanged; only the geometry is restricted
# to the physical reciprocal plane.  Two facts make this robust:
#   * The reciprocal vector dual to the real vacuum axis is always along the
#     layer normal, so the in-plane condition is *fractional* k=0 along that
#     axis (not Cartesian kz=0).
#   * seekpath may permute axes when it standardizes (monoclinic -> unique
#     axis b), so the vacuum axis is detected in the standardized frame rather
#     than taken from the user's --vacuum-axis flag.
# ---------------------------------------------------------------------------

def detect_vacuum_axis_2d(b_matrix, sep_ratio=0.8, ortho_tol=0.05):
    """Identify the slab vacuum axis in the standardized reciprocal frame.

    The vacuum (out-of-plane) real axis is the longest, so its reciprocal
    vector has the smallest norm.  The choice is cross-checked for
    orthogonality to the in-plane reciprocal vectors (the vacuum reciprocal
    vector lies along the layer normal).  For a proper slab the smallest
    reciprocal norm is well separated and both checks pass; the flags let the
    caller warn on an ambiguous cell.

    Returns ``(axis_index, info)`` where ``info`` carries the diagnostics.
    """
    b = np.array(b_matrix, dtype=float)
    norms = np.linalg.norm(b, axis=1)
    axis = int(np.argmin(norms))
    others = [i for i in range(3) if i != axis]
    bn = b / norms[:, None]
    dots = [abs(float(bn[axis] @ bn[j])) for j in others]
    sorted_norms = np.sort(norms)
    separated = bool(sorted_norms[0] < sep_ratio * sorted_norms[1])
    orthogonal = bool(all(d < ortho_tol for d in dots))
    info = {
        "reciprocal_norms": [float(x) for x in norms],
        "dots_to_in_plane": [float(d) for d in dots],
        "separated": separated,
        "orthogonal": orthogonal,
    }
    return axis, info


def area_centroid_2d(frac_points, vacuum_axis, b_matrix):
    """Area centroid of the in-plane IBZ polygon.

    ``frac_points`` are fractional coordinates already restricted to the
    physical plane (``k[vacuum_axis] == 0``).  The centroid is computed in the
    in-plane fractional coordinates and mapped back through ``b_matrix``; this
    is exact because the area centroid is affine-equivariant, so a monoclinic
    in-plane shear does not require any special handling.

    Returns ``(centroid_frac, centroid_cart, area_frac)`` with the vacuum
    fractional component set to exactly 0.
    """
    in_plane = [i for i in range(3) if i != vacuum_axis]
    pts2 = np.array(frac_points, dtype=float)[:, in_plane]
    uniq = []
    for p in pts2:
        if not any(np.allclose(p, q, atol=1e-8) for q in uniq):
            uniq.append(p)
    pts2 = np.array(uniq, dtype=float)
    area = 0.0
    if len(pts2) < 3:
        c2 = pts2.mean(axis=0) if len(pts2) else np.zeros(2)
    else:
        hull = ConvexHull(pts2)
        poly = pts2[hull.vertices]
        x, y = poly[:, 0], poly[:, 1]
        cross = x * np.roll(y, -1) - np.roll(x, -1) * y
        a2 = float(np.sum(cross))
        if abs(a2) < 1e-12:
            c2 = poly.mean(axis=0)
        else:
            cx = float(np.sum((x + np.roll(x, -1)) * cross) / (3.0 * a2))
            cy = float(np.sum((y + np.roll(y, -1)) * cross) / (3.0 * a2))
            c2 = np.array([cx, cy], dtype=float)
            area = abs(a2) / 2.0
    centroid_frac = np.zeros(3)
    centroid_frac[in_plane[0]] = c2[0]
    centroid_frac[in_plane[1]] = c2[1]
    centroid_frac[vacuum_axis] = 0.0
    centroid_cart = centroid_frac @ np.array(b_matrix, dtype=float)
    return centroid_frac, centroid_cart, area


def ordered_2d_polygon_frac(frac_points, vacuum_axis):
    """Return the in-plane IBZ polygon vertices (3D fractional, vacuum=0),
    ordered around the convex hull.  Used by the 2D spin-pattern figure."""
    in_plane = [i for i in range(3) if i != vacuum_axis]
    pts = np.array(frac_points, dtype=float)
    uniq = []
    for p in pts:
        if not any(np.allclose(p, q, atol=1e-8) for q in uniq):
            uniq.append(p)
    uniq = np.array(uniq, dtype=float)
    if len(uniq) < 3:
        return uniq.tolist()
    hull = ConvexHull(uniq[:, in_plane])
    return uniq[hull.vertices].tolist()


def keeps_2d_plane(R, vacuum_axis, tol=1e-6):
    """True if operation R maps in-plane k back into the plane (Filter 1).

    For a true layer group every operation passes; this guards against spurious
    3D bulk operations (C4x, C3[111], ...) that would send an in-plane k out of
    the vacuum k=0 plane.
    """
    R = np.array(R, dtype=float)
    in_plane = [i for i in range(3) if i != vacuum_axis]
    return all(abs(R[vacuum_axis, i]) < tol for i in in_plane)


def is_trivial_2d_spin_flip(R, vacuum_axis, tol=1e-6):
    """True if R's in-plane 2x2 block is +-I (Filter 2 -- a trivial 2D flip).

    +I  (e.g. mz)            : k_par -> k_par, forces E_up = E_down everywhere.
    -I  (e.g. C2z, inversion): k_par -> -k_par, PT-like for collinear spins.
    Either way there is no in-plane altermagnetic splitting.
    """
    R = np.array(R, dtype=float)
    in_plane = [i for i in range(3) if i != vacuum_axis]
    block = R[np.ix_(in_plane, in_plane)]
    return (np.allclose(block, np.eye(2), atol=tol)
            or np.allclose(block, -np.eye(2), atol=tol))


def is_valid_2d_spin_flip(R, vacuum_axis, tol=1e-6):
    """True if R is a genuine 2D spin-flip: keeps the plane and is non-trivial."""
    return (keeps_2d_plane(R, vacuum_axis, tol)
            and not is_trivial_2d_spin_flip(R, vacuum_axis, tol))


def check_input_slab(a_matrix, declared_axis, ortho_tol=0.02):
    """Sanity-check that the INPUT cell looks like a proper 2D slab.

    A real 2D slab has its vacuum axis orthogonal to the two in-plane axes and
    longer than them.  A tilted/elongated bulk cell (the common mistake of
    just stretching c of a 3D structure) fails these.  Returns a list of
    human-readable warning strings (empty when the cell looks fine).
    """
    A = np.array(a_matrix, dtype=float)
    lengths = np.linalg.norm(A, axis=1)
    others = [i for i in range(3) if i != declared_axis]
    An = A / lengths[:, None]
    dots = [abs(float(An[declared_axis] @ An[j])) for j in others]
    name = "abc"[declared_axis]
    warns = []
    if any(d > ortho_tol for d in dots):
        warns.append(
            f"input vacuum axis '{name}' is not orthogonal to the in-plane axes "
            f"(|cos| = {[round(d, 3) for d in dots]}); this does not look like a "
            "proper 2D slab (a tilted or elongated bulk cell?)."
        )
    if lengths[declared_axis] < max(lengths[j] for j in others):
        warns.append(
            f"input vacuum axis '{name}' (length {lengths[declared_axis]:.2f}) is not "
            "the longest axis; is there real vacuum along it?"
        )
    return warns


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


# ============================================================================
# IBZ frame edge helper
# ============================================================================
def _get_ibz_frame_edges(hull_pts, hull_simplices, hull_labels=None):
    """Return only the non-coplanar edges of the IBZ hull as (pt1, pt2) pairs.

    Filters out internal triangulation diagonals within flat faces by checking
    whether adjacent face normals are nearly parallel (|cos 闁煎啿鍓?--0.99).
    """
    from collections import defaultdict
    hull_pts = np.array(hull_pts)
    hull_labels = list(hull_labels) if hull_labels is not None else None
    forbidden_label_edges = {
        tuple(sorted(edge))
        for edge in (
            ("G", "G_6"), ("G_2", "G_4"),
            ("N", "N_6"), ("N_2", "N_4"),
        )
    }

    def _is_forbidden_label_edge(i, j):
        if hull_labels is None or i >= len(hull_labels) or j >= len(hull_labels):
            return False
        edge = tuple(sorted((hull_labels[i], hull_labels[j])))
        return edge in forbidden_label_edges

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
        if _is_forbidden_label_edge(i, j):
            continue
        if len(faces) < 2:
            edges.append((hull_pts[i], hull_pts[j]))
        else:
            cos_a = max(
                abs(np.dot(face_normals[a], face_normals[b]))
                for a in faces for b in faces if a != b
            )
            if cos_a < 0.97:
                edges.append((hull_pts[i], hull_pts[j]))
    return edges


# ============================================================================
# Spin-flip operation classification and Figure 2 geometric visuals
# ============================================================================

def _perp_unit(v):
    """Return a unit vector perpendicular to v (for any nonzero v)."""
    v = np.asarray(v, dtype=float)
    idx = int(np.argmin(np.abs(v)))
    w = np.zeros(3)
    w[idx] = 1.0
    w = w - (w @ v) * v
    return w / np.linalg.norm(w)


def _axis_bz_exit(axis, bz_loops):
    """
    Return the parameter t where the ray origin + t*axis first exits the BZ.
    Uses the convex-hull half-space equations of the BZ vertices.
    Falls back to bz_radius if the hull cannot be computed.
    """
    all_pts = np.vstack([np.asarray(loop, dtype=float) for loop in bz_loops])
    fallback = float(np.max(np.linalg.norm(all_pts, axis=1)))
    try:
        hull = ConvexHull(all_pts)
        t_vals = []
        for eq in hull.equations:
            n, d = eq[:3], eq[3]
            denom = float(n @ axis)
            if denom > 1e-10:
                t = float(-d / denom)
                if t > 0:
                    t_vals.append(t)
        return min(t_vals) if t_vals else fallback
    except Exception:
        return fallback


def _classify_spinflip_op(R_cart):
    """
    Classify a Cartesian orthogonal matrix as a crystallographic point-group
    operation using det and trace.

    Returns a dict:
      'type'  : 'identity' | 'rotation' | 'mirror' | 'inversion' | 'rotoreflection'
      'axis'  : unit 3-vector or None
                  rotation -> rotation axis (real eigenvec with eigenvalue +1)
                  mirror   -> plane normal  (real eigenvec with eigenvalue -1)
                  Sn       -> rotation axis (real eigenvec with eigenvalue -1)
      'order' : int n for Cn/Sn, or None
    """
    det = int(np.round(np.linalg.det(R_cart)))
    tr  = int(np.round(np.trace(R_cart)))
    eigvals, eigvecs = np.linalg.eig(R_cart)

    def _real_eigvec(target):
        dists = np.where(
            np.abs(eigvals.imag) < 1e-5,
            np.abs(eigvals.real - target),
            np.inf,
        )
        idx = int(np.argmin(dists))
        if dists[idx] > 0.15:
            return None
        vec = eigvecs[:, idx].real
        norm = np.linalg.norm(vec)
        return vec / norm if norm > 1e-10 else None

    if det == 1:
        if tr == 3:
            return {'type': 'identity', 'axis': None, 'order': 1}
        order = {-1: 2, 0: 3, 1: 4, 2: 6}.get(tr)
        return {'type': 'rotation', 'axis': _real_eigvec(+1.0), 'order': order}

    # det == -1 (improper)
    if tr == -3:
        return {'type': 'inversion', 'axis': None, 'order': None}
    if tr == 1:
        return {'type': 'mirror', 'axis': _real_eigvec(-1.0), 'order': None}
    # Rotoreflections Sn: S3 tr=-2, S4 tr=-1, S6 tr=0
    order = {-2: 3, -1: 4, 0: 6}.get(tr)
    return {'type': 'rotoreflection', 'axis': _real_eigvec(-1.0), 'order': order}


def _mirror_plane_bz_polygon(normal, bz_loops):
    """
    Return the polygon (N,3) where the mirror plane n·k=0 cuts the BZ edges,
    sorted angularly around the centroid.  Returns None if fewer than 3 pts.
    """
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    pts = []
    for loop in bz_loops:
        loop_pts = np.asarray(loop, dtype=float)
        for a, b in zip(loop_pts[:-1], loop_pts[1:]):
            da = float(n @ a)
            db = float(n @ b)
            if abs(da) < 1e-8:
                pts.append(a.copy())
            if abs(db) < 1e-8:
                pts.append(b.copy())
            if da * db < -1e-14:
                t = da / (da - db)
                pts.append(a + t * (b - a))

    if len(pts) < 3:
        return None
    pts = np.unique(np.round(np.array(pts, dtype=float), 10), axis=0)
    if len(pts) < 3:
        return None

    centroid = pts.mean(axis=0)
    u = _perp_unit(n)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    uv = np.column_stack([u, v])
    coords_2d = pts @ uv
    c2d = centroid @ uv
    angles = np.arctan2(coords_2d[:, 1] - c2d[1], coords_2d[:, 0] - c2d[0])
    return pts[np.argsort(angles)]


def _reduce_int_vector(vec):
    """Reduce a direction to smallest integer indices, with the sign convention
    that the first non-zero index is positive. Used for both reciprocal (hkl)
    plane normals and direct [uvw] axes (caller supplies the projected vector)."""
    v = np.asarray(vec, dtype=float)
    m = np.max(np.abs(v))
    if m < 1e-9:
        return [0, 0, 0]
    v = v / m
    ints = None
    for denom in range(1, 13):
        scaled = v * denom
        if np.all(np.abs(scaled - np.round(scaled)) < 1e-3):
            ints = np.round(scaled).astype(int)
            break
    if ints is None:
        ints = np.round(v * 12).astype(int)
    g = int(np.gcd.reduce(np.abs(ints)))
    if g > 0:
        ints = ints // g
    for x in ints:
        if x != 0:
            if x < 0:
                ints = -ints
            break
    return ints.tolist()


def _format_miller(letter, idx):
    """Bold mathtext label, e.g. 'm' + [1,-1,0] -> $\\mathbf{m_{1\\bar{1}0}}$."""
    sub = "".join(rf"\bar{{{abs(i)}}}" if i < 0 else f"{i}" for i in idx)
    return rf"$\mathbf{{{letter}_{{{sub}}}}}$"


def _rotation_sense(R_cart, axis):
    """Sign of the rotation angle about `axis` for a proper or improper rotation:
    +1 (counter-clockwise about +axis), -1, or 0 when undefined (order 2, where
    sin(theta)=0 and +/- coincide). Works for Cn and Sn: the antisymmetric part
    of R equals sin(theta)*axis in both cases (the on-axis +-1 is symmetric)."""
    R = np.asarray(R_cart, dtype=float)
    a = np.asarray(axis, dtype=float)
    a = a / (np.linalg.norm(a) or 1.0)
    w = 0.5 * np.array([R[2, 1] - R[1, 2],
                        R[0, 2] - R[2, 0],
                        R[1, 0] - R[0, 1]])
    s = float(w @ a)
    if abs(s) < 1e-6:
        return 0
    return 1 if s > 0 else -1


def describe_spinflip_op(R_cart, b_matrix=None):
    """Plain-text crystallographic name for a Cartesian point operation.

    type/order come from det/trace (basis-invariant). The axis/normal is read as
    the Cartesian eigenvector and expressed as reduced integer components in the
    reciprocal (b1,b2,b3) frame shown in the figures (vec @ inv(b_matrix)):
      rotation Cn   -> axis [hkl] (rotation-axis direction in b1,b2,b3)
      Sn            -> axis [hkl]
      mirror m      -> plane (hkl) (plane normal in b1,b2,b3)
      inversion / identity -> no axis
    Square brackets mark a direction (rotation axis), parentheses a plane
    (mirror). If b_matrix is None the axis is omitted (type/order still correct).
    """
    op = _classify_spinflip_op(np.asarray(R_cart, dtype=float))
    t, order = op['type'], op['order']
    if t == 'identity':
        return 'identity (E)'
    if t == 'inversion':
        return 'inversion (i)'
    if t == 'mirror':
        if b_matrix is None or op['axis'] is None:
            return 'mirror m'
        hkl = _reduce_int_vector(np.asarray(op['axis']) @ np.linalg.inv(b_matrix))
        return f"mirror m ({' '.join(str(i) for i in hkl)})"
    if t in ('rotation', 'rotoreflection'):
        sym = 'C' if t == 'rotation' else 'S'
        if b_matrix is None or op['axis'] is None:
            return f"{sym}{order}"
        axis = np.asarray(op['axis'], dtype=float)
        # Axis components in the reciprocal (b1,b2,b3) frame, same projection the
        # mirror normal uses, so both labels read against the drawn axes.
        hkl_raw = axis @ np.linalg.inv(np.asarray(b_matrix))
        # Match the first-non-zero-positive convention of the reported indices,
        # and measure the rotation sense about that same reported axis direction.
        flip = 1.0
        for val in hkl_raw:
            if abs(val) > 1e-9:
                flip = 1.0 if val > 0 else -1.0
                break
        sense = _rotation_sense(R_cart, flip * axis)
        sgn = '' if (order < 3 or sense == 0) else ('+' if sense > 0 else '-')
        hkl = _reduce_int_vector(hkl_raw)
        return f"{sym}{order}{sgn} [{' '.join(str(i) for i in hkl)}]"
    return t


# ============================================================================
# Spin-flip Connection Figure  (replaces the old "mapped BZ" Fig 2)
# ============================================================================


# ============================================================================
# Spin-up / Spin-down Full-BZ Coloring Figure  (replaces old rainbow Fig 2)
# ============================================================================


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


def _bz_halfspaces(b_matrix, grid_radius=2):
    """Return halfspaces for the Wigner-Seitz BZ in Cartesian k coordinates."""
    halfspaces = []
    for h in range(-grid_radius, grid_radius + 1):
        for k in range(-grid_radius, grid_radius + 1):
            for l in range(-grid_radius, grid_radius + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                G = h * b_matrix[0] + k * b_matrix[1] + l * b_matrix[2]
                norm2 = float(np.dot(G, G))
                if norm2 < 1e-12:
                    continue
                # Points closer to 0 than to G satisfy G.x <= |G|^2/2.
                halfspaces.append(np.r_[G, -0.5 * norm2])
    return np.array(halfspaces, dtype=float)


def _dedupe_points(points, decimals=10):
    if len(points) == 0:
        return np.empty((0, 3))
    return np.unique(np.round(np.array(points, dtype=float), decimals), axis=0)


def _spin_bz_cells(b_matrix, unique_ops, centroid_cart):
    """Build non-overlapping symmetry/Voronoi cells clipped to the BZ.

    The high-symmetry point convex hull can be slightly too large for skew
    monoclinic cells.  For full spin-BZ coloring, construct the actual partition
    generated by the symmetry images of the IBZ centroid instead.
    """
    b_matrix = np.array(b_matrix, dtype=float)
    centers = np.array([g @ centroid_cart for g in unique_ops], dtype=float)
    bz_hs = _bz_halfspaces(b_matrix)
    cells = []

    for i, ci in enumerate(centers):
        hs = [*bz_hs]
        for j, cj in enumerate(centers):
            if i == j or np.linalg.norm(ci - cj) < 1e-10:
                continue
            # Closer to ci than to cj:
            # ||x-ci||^2 <= ||x-cj||^2
            # 2(cj-ci).x + |ci|^2 - |cj|^2 <= 0
            normal = 2.0 * (cj - ci)
            offset = float(np.dot(ci, ci) - np.dot(cj, cj))
            hs.append(np.r_[normal, offset])
        hs = np.array(hs, dtype=float)

        interior = ci.copy()
        vals = hs[:, :3] @ interior + hs[:, 3]
        if np.max(vals) >= -1e-9:
            # Move very slightly toward the BZ center if the centroid image lies
            # on a numerical boundary.  This preserves the intended cell.
            interior = ci * (1.0 - 1e-7)
        vals = hs[:, :3] @ interior + hs[:, 3]
        if np.max(vals) >= -1e-9:
            cells.append((None, None))
            continue

        try:
            hs_int = HalfspaceIntersection(hs, interior)
            verts = _dedupe_points(hs_int.intersections)
            if len(verts) < 4:
                cells.append((None, None))
                continue
            hull = ConvexHull(verts)
            cells.append((verts, hull.simplices))
        except Exception:
            cells.append((None, None))
    return cells


def _fractional_real_op_to_cart_k(b_matrix, operation):
    """Convert a real-space fractional operation to its Cartesian k action."""
    b_t = np.array(b_matrix, dtype=float).T
    operation = np.array(operation, dtype=float)
    return b_t @ np.linalg.inv(operation).T @ np.linalg.inv(b_t)


def _mapped_spin_hulls(b_matrix, hull_pts, hull_simplices, preserve_ops_frac, flip_ops_frac):
    """Map the exact Figure 1 hull by labeled FindSpinGroup operations."""
    labeled_ops = []
    for ops, is_down in ((preserve_ops_frac, False), (flip_ops_frac, True)):
        for operation in ops:
            cart_op = _fractional_real_op_to_cart_k(b_matrix, operation)
            matching = [
                old_down for old_op, old_down in labeled_ops
                if np.allclose(cart_op, old_op, atol=1e-7)
            ]
            if matching:
                if matching[0] != is_down:
                    print("[Warning] A k-space operation is both spin-preserving and spin-flipping.")
                    return None
                continue
            labeled_ops.append((cart_op, is_down))

    return [
        ((operation @ np.asarray(hull_pts, dtype=float).T).T, hull_simplices, is_down)
        for operation, is_down in labeled_ops
    ]


def build_symmetry_ibz_cell(b_matrix, unique_ops, seed_cart):
    """Build the fundamental BZ cell selected by a generic seed k-point."""
    seed_cart = np.array(seed_cart, dtype=float)
    cells = _spin_bz_cells(b_matrix, unique_ops, seed_cart)
    centers = np.array([g @ seed_cart for g in unique_ops], dtype=float)
    if not len(centers):
        return None, None

    order = np.argsort(np.linalg.norm(centers - seed_cart, axis=1))
    for idx in order:
        pts, simplices = cells[int(idx)]
        if pts is not None and simplices is not None:
            return np.array(pts, dtype=float), np.array(simplices, dtype=int)
    return None, None


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


# ============================================================================
# Main Pipeline
# ============================================================================
def run(
    filename,
    output_dir=None,
    show_plot=True,
    defer_show=False,
    verbose=True,
    seekpath_type_numbers=None,
    mode_2d=False,
    input_vacuum_axis=2,
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
    elev1, azim1 = 14, 20
    fig1_path = os.path.join(output_dir, f'{basename}_ibz_{sc_display}.png')
    if show_plot and not mode_2d:
        # Interactive mode: create the figure now. alterseek_path can defer
        # the actual plt.show() call until all prompts and file writes finish.
        fig1, ax1 = setup_3d_ax(fig1_title,
                                bz_loops, b_matrix, bz_center, bz_span)
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
            elev1, azim1 = 14, 20
        else:
            plt.show()
            elev1, azim1 = ax1.elev, ax1.azim
    else:
        # Automated mode (called from alterseek_path): use default angles, no window
        elev1, azim1 = 14, 20

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


# ==========================================================================
# Restructuring phase 1: figure code moved to plotting_common / plotting_3d.
# These re-exports keep `from compute_centroid_hybrid import <fn>` working for
# existing callers and tests. Placed at end so cc is fully defined first;
# removed in a later phase once geometry/symmetry are also extracted.
# ==========================================================================
from plotting_common import (
    _get_bz_path_style,
    _figure_output_paths,
    _save_figure,
    _print_saved_paths,
    _math_label,
)
from plotting_3d import (
    _Arrow3D,
    _draw_ibz_faces_by_sector,
    _get_view_direction,
    _classify_bz_edges,
    draw_bz_edges,
    setup_3d_ax,
    plot_ibz,
    plot_mapped_bz,
    _screen_xy,
    _best_label_anchor,
    _mirror_label_candidates,
    _draw_op_visual,
    _connect_arc_view_follow,
    plot_spin_flip_figure,
    plot_spin_bz_figure,
    draw_kz0_helper_plane,
    draw_projected_reciprocal_axes,
    plot_spin_bz_top_view_figure,
)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python compute_centroid_hybrid.py <structure_file> [output_dir]")
        sys.exit(1)

    structure_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None
    results = run(structure_file, out_dir)

