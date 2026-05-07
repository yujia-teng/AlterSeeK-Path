#!/usr/bin/env python3
"""
2D/slab altermagnetic k-path generator.

This experimental entry point keeps the existing collinear spin-flip workflow,
but restricts the band path and general k point to the physical 2D reciprocal
plane.  The ordinary 3D workflow in alterseek_path.py is intentionally left
untouched until the 2D behavior is validated.
"""

import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from scipy.spatial import ConvexHull
except ImportError:
    ConvexHull = None

from alterseek_path import (
    FIND_SF_AVAILABLE,
    KPointsModifier,
    find_sf_run,
)
from compute_centroid_hybrid import run as compute_centroid
from compute_centroid_hybrid import (
    _points_on_kz_plane,
    _bz_kz_plane_outline,
)


VACUUM_AXIS = 2
VACUUM_AXIS_NAME = "c"


def reciprocal_lattice_rows(lattice) -> np.ndarray:
    """Return reciprocal basis vectors as rows for a real-space row-vector lattice."""
    return 2 * np.pi * np.linalg.inv(np.array(lattice, dtype=float)).T


def convert_real_frac_op_to_primitive_k_basis(
    R, source_lattice, primitive_b_matrix, standardization_rotation=None
):
    """Convert a real-space fractional operation to the seekpath primitive basis."""
    b_input = reciprocal_lattice_rows(source_lattice)
    b_prim = np.array(primitive_b_matrix, dtype=float)
    R_arr = np.array(R, dtype=float)
    R_cart = b_input.T @ np.linalg.inv(R_arr).T @ np.linalg.inv(b_input.T)
    if standardization_rotation is not None:
        q = np.array(standardization_rotation, dtype=float)
        R_cart = q @ R_cart @ np.linalg.inv(q)
    R_prim_inv_T = np.linalg.inv(b_prim.T) @ R_cart @ b_prim.T
    return np.linalg.inv(R_prim_inv_T.T), R_cart


def _plane_value(value: float) -> float:
    """Return distance from the k_i = 0 reciprocal plane, modulo integers."""
    return abs(((float(value) + 0.5) % 1.0) - 0.5)


def _on_plane(coords, axis: int, tol: float) -> bool:
    return _plane_value(coords[axis]) <= tol


def _project_to_plane(coords, axis: int) -> List[float]:
    projected = [float(coords[0]), float(coords[1]), float(coords[2])]
    projected[axis] = 0.0
    return projected


def _is_monoclinic_result(result: dict) -> bool:
    key = str(result.get("lattice_key", result.get("sc_type", "")))
    return key in {"mP1", "mC1", "mC2", "mC3"}


def _apply_2d_monoclinic_seekpath_basis(result: dict) -> dict:
    """Record monoclinic 2D basis context without overriding the slab plane."""
    if not _is_monoclinic_result(result):
        return result
    result["coordinate_basis"] = "monoclinic_seekpath"
    result["coordinate_basis_note"] = (
        "2D monoclinic: seekpath may permute the input cell axes; "
        "the 2D reciprocal plane is chosen from the slab vacuum direction."
    )
    return result


def _frac_point_between_bases(point, from_b_matrix, to_b_matrix):
    """Convert fractional reciprocal coordinates between two reciprocal bases."""
    cart = np.array(point, dtype=float) @ np.array(from_b_matrix, dtype=float)
    return (cart @ np.linalg.inv(np.array(to_b_matrix, dtype=float))).tolist()


def _frac_dict_between_bases(coords: dict, from_b_matrix, to_b_matrix) -> dict:
    return {
        label: _frac_point_between_bases(point, from_b_matrix, to_b_matrix)
        for label, point in coords.items()
    }


def _cart_ops_between_reciprocal_bases(ops, from_b_matrix, to_b_matrix):
    """Convert Cartesian k-ops expressed through one reciprocal basis to another."""
    b_from = np.array(from_b_matrix, dtype=float)
    b_to = np.array(to_b_matrix, dtype=float)
    # Row fractional coordinates obey f_to = f_from @ p.
    # Column k-operation matrices therefore transform as K_to = p.T K_from inv(p.T).
    p = b_from @ np.linalg.inv(b_to)
    p_t = p.T
    p_t_inv = np.linalg.inv(p_t)
    from_t = b_from.T
    from_t_inv = np.linalg.inv(from_t)
    to_t = b_to.T
    to_t_inv = np.linalg.inv(to_t)
    converted = []
    for op in ops:
        k_from_frac = from_t_inv @ np.array(op, dtype=float) @ from_t
        k_to_frac = p_t @ k_from_frac @ p_t_inv
        converted.append(to_t @ k_to_frac @ to_t_inv)
    return converted


def _convert_result_to_input_slab_basis(result: dict) -> dict:
    """Use the input reciprocal basis as the user-facing 2D slab convention."""
    b_seek = np.array(result["b_matrix"], dtype=float)
    b_input = np.array(result.get("b_matrix_input", b_seek), dtype=float)
    result["seekpath_b_matrix"] = b_seek
    result["b_matrix"] = b_input
    result["coordinate_basis"] = "input_slab"
    result["coordinate_basis_note"] = (
        "2D slab output uses the input reciprocal basis; vacuum is input c, "
        "so the selected plane is input kz=0."
    )
    for key in (
        "band_kpoints_frac",
        "path_kpoints_frac",
        "ibz_kpoints_frac",
        "kpoints_frac",
        "sp_point_coords",
    ):
        if isinstance(result.get(key), dict):
            result[key] = _frac_dict_between_bases(result[key], b_seek, b_input)
    if result.get("centroid_frac") is not None:
        result["centroid_frac"] = _frac_point_between_bases(
            result["centroid_frac"], b_seek, b_input)
        result["centroid_cart"] = np.array(result["centroid_frac"]) @ b_input
    for key in ("unique_ops", "sym_ops_cart"):
        if result.get(key):
            result[key] = _cart_ops_between_reciprocal_bases(
                result[key], b_seek, b_input)
    return result


def _infer_vacuum_axis_from_reciprocal_metric(result: dict, default_axis: int) -> int:
    """Choose the slab normal as the shortest reciprocal basis vector."""
    b_matrix = result.get("b_matrix")
    result["vacuum_axis_detection"] = {
        "method": "fallback",
        "default_axis": int(default_axis),
        "axis": int(default_axis),
        "reciprocal_norms": None,
    }
    if b_matrix is None:
        return int(default_axis)
    norms = np.linalg.norm(np.array(b_matrix, dtype=float), axis=1)
    result["vacuum_axis_detection"]["reciprocal_norms"] = [float(x) for x in norms]
    if len(norms) != 3 or not np.all(np.isfinite(norms)):
        return int(default_axis)
    order = np.argsort(norms)
    # A slab vacuum produces a clear small reciprocal vector.  When the metric is
    # ambiguous, keep the caller's requested axis instead of guessing.
    if norms[order[0]] < 0.75 * norms[order[1]]:
        result["vacuum_axis_detection"].update({
            "method": "shortest_reciprocal_vector",
            "axis": int(order[0]),
        })
        return int(order[0])
    return int(default_axis)


def _forbidden_2d_altermagnetic_spin_flip_operation(ops, axis: int = VACUUM_AXIS, tol=1e-8):
    """Return spin-flip C2 or mirror normal to the selected 2D plane."""
    c2 = np.eye(3)
    c2[(axis + 1) % 3, (axis + 1) % 3] = -1.0
    c2[(axis + 2) % 3, (axis + 2) % 3] = -1.0
    mirror = np.eye(3)
    mirror[axis, axis] = -1.0
    axis_label = ("x", "y", "z")[axis]
    forbidden = [
        (f"C2{axis_label}", c2),
        (f"m{axis_label}", mirror),
    ]
    for op in ops:
        op = np.array(op, dtype=float)
        for name, matrix in forbidden:
            if np.allclose(op, matrix, atol=tol):
                return name, op
    return None, None


def _operation_identity_on_plane(op, axis: int, tol=1e-8) -> bool:
    """Return True if a full 3D operation leaves all in-plane k coordinates fixed."""
    op = np.array(op, dtype=float)
    in_plane_axes = [i for i in range(3) if i != axis]
    for i in in_plane_axes:
        row = np.zeros(3)
        row[i] = 1.0
        if not np.allclose(op[i], row, atol=tol):
            return False
    return True


def _same_in_plane_operation(a, b, axis: int, tol=1e-6) -> bool:
    """Compare only the selected 2D reciprocal-plane action of two operations."""
    in_plane_axes = [i for i in range(3) if i != axis]
    a2 = np.array(a, dtype=float)[np.ix_(in_plane_axes, in_plane_axes)]
    b2 = np.array(b, dtype=float)[np.ix_(in_plane_axes, in_plane_axes)]
    return np.allclose(a2, b2, atol=tol)


def _classify_spin_down_ops_2d(b_matrix, unique_ops, R, flip_ops_frac=None,
                               axis: int = VACUUM_AXIS):
    """Classify spin-down domains by matching the in-plane operation only."""
    if not unique_ops:
        return np.zeros(0, dtype=bool)

    b_matrix = np.array(b_matrix, dtype=float)
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)
    if flip_ops_frac is None or not len(flip_ops_frac):
        flip_ops_frac = [R]
    flip_set = [np.array(f, dtype=float) for f in flip_ops_frac]

    spin_down_mask = np.zeros(len(unique_ops), dtype=bool)
    for i, g_cart in enumerate(unique_ops):
        M = b_T_inv @ np.array(g_cart, dtype=float) @ b_T
        g_frac = np.linalg.inv(M.T)
        spin_down_mask[i] = any(
            _same_in_plane_operation(g_frac, f, axis)
            for f in flip_set
        )
    return spin_down_mask


def _dedupe_in_plane_ops_2d(b_matrix, unique_ops, spin_down_mask,
                            axis: int = VACUUM_AXIS, decimals: int = 8):
    """Keep one 3D symmetry representative for each distinct in-plane action."""
    b_matrix = np.array(b_matrix, dtype=float)
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)
    in_plane_axes = [i for i in range(3) if i != axis]
    kept_ops = []
    kept_mask = []
    seen = {}
    conflicts = 0
    for i, g_cart in enumerate(unique_ops):
        M = b_T_inv @ np.array(g_cart, dtype=float) @ b_T
        g_frac = np.linalg.inv(M.T)
        block = g_frac[np.ix_(in_plane_axes, in_plane_axes)]
        key = tuple(np.round(block, decimals=decimals).ravel())
        is_down = bool(spin_down_mask[i])
        if key in seen:
            if seen[key] != is_down:
                conflicts += 1
            continue
        seen[key] = is_down
        kept_ops.append(g_cart)
        kept_mask.append(is_down)
    if conflicts:
        print(
            "[Warning] Some 3D operations collapse to the same 2D operation "
            "but have different spin labels; those 2D sectors may be spin-degenerate."
        )
    return kept_ops, np.array(kept_mask, dtype=bool)


def _default_2d_axis_from_symmetry(sf_result: dict) -> int:
    """Return the legacy c-vacuum default used before metric-based detection."""
    return VACUUM_AXIS


def _dedupe_segments(path, coords: Dict[str, List[float]], axis: int, tol: float):
    filtered = []
    seen = set()
    for start, end in path:
        if start not in coords or end not in coords:
            continue
        if not (_on_plane(coords[start], axis, tol) and _on_plane(coords[end], axis, tol)):
            continue
        a = _project_to_plane(coords[start], axis)
        b = _project_to_plane(coords[end], axis)
        if np.allclose(a, b, atol=tol):
            continue
        key = (start, end)
        if key in seen:
            continue
        seen.add(key)
        filtered.append((start, end))
    return filtered


def _segments_from_kpoints_data(kpoints_data, axis: int, tol: float):
    """Build 2D path/coordinate dictionaries from a line-mode KPOINTS file."""
    coords: Dict[str, List[float]] = {}
    path = []
    if len(kpoints_data) < 2:
        raise ValueError("manual 2D KPATH must contain at least one segment")
    for i in range(0, len(kpoints_data) - 1, 2):
        start = kpoints_data[i]
        end = kpoints_data[i + 1]
        if start is None or end is None:
            continue
        start_label = str(start[3])
        end_label = str(end[3])
        start_point = _project_to_plane(start[:3], axis)
        end_point = _project_to_plane(end[:3], axis)
        if np.allclose(start_point, end_point, atol=tol):
            print(f"  [Note] Manual segment {start_label}-{end_label} skipped: endpoints coincide")
            continue
        coords.setdefault(start_label, start_point)
        coords.setdefault(end_label, end_point)
        path.append((start_label, end_label))
    if not path:
        raise ValueError("manual 2D KPATH did not contain any non-degenerate segments")
    return path, coords


def _apply_manual_2d_path(centroid_result: dict, modifier: KPointsModifier,
                          axis: int, tol: float):
    """Replace the automatic 2D band path with a user-supplied KPATH file."""
    path, coords = _segments_from_kpoints_data(modifier.kpoints_data, axis, tol)
    centroid_result["band_kpath"] = path
    centroid_result["band_kpoints_frac"] = coords
    centroid_result["extra_general_vertices"] = []
    labels = _dedupe_2d_labels(list(coords.keys()), coords, axis, tol)
    if len(labels) >= 3:
        centroid = _compute_2d_centroid(coords, labels, axis, tol)
        centroid_result["centroid_frac"] = centroid
        centroid_result["centroid_cart"] = np.array(centroid) @ np.array(centroid_result["b_matrix"])
        centroid_result["ibz_polygon_frac"] = _ordered_2d_polygon_frac(coords, labels, axis, tol)
        print("[Note] Recomputed 2D centroid from the manual KPATH labels.")
    else:
        print("[Warning] Manual KPATH has fewer than three unique labels; keeping the automatic 2D centroid.")
    return path, coords


def _axis_name(axis: int) -> str:
    return ("kx", "ky", "kz")[axis]


def _polygon_centroid_2d(points_2d: np.ndarray) -> np.ndarray:
    """Centroid of the convex polygon spanned by 2D points."""
    if len(points_2d) == 0:
        raise ValueError("no 2D points available")
    unique = []
    for point in points_2d:
        if not any(np.allclose(point, old, atol=1e-8) for old in unique):
            unique.append(point)
    pts = np.array(unique, dtype=float)
    if len(pts) < 3 or ConvexHull is None:
        return np.mean(pts, axis=0)

    hull = ConvexHull(pts)
    poly = pts[hull.vertices]
    x = poly[:, 0]
    y = poly[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    area2 = np.sum(cross)
    if abs(area2) < 1e-12:
        return np.mean(poly, axis=0)
    cx = np.sum((x + np.roll(x, -1)) * cross) / (3.0 * area2)
    cy = np.sum((y + np.roll(y, -1)) * cross) / (3.0 * area2)
    return np.array([cx, cy], dtype=float)


def _compute_2d_centroid(
    coords: Dict[str, List[float]],
    labels: List[str],
    axis: int,
    tol: float,
) -> List[float]:
    in_plane_axes = [i for i in range(3) if i != axis]
    points = []
    for label in labels:
        if label in coords and _on_plane(coords[label], axis, tol):
            projected = _project_to_plane(coords[label], axis)
            points.append([projected[in_plane_axes[0]], projected[in_plane_axes[1]]])
    centroid_2d = _polygon_centroid_2d(np.array(points, dtype=float))
    centroid = [0.0, 0.0, 0.0]
    centroid[in_plane_axes[0]] = float(centroid_2d[0])
    centroid[in_plane_axes[1]] = float(centroid_2d[1])
    centroid[axis] = 0.0
    return centroid


def _base_2d_label(label: str) -> str:
    """Return the physical 2D label, dropping project-only doubled-IBZ suffixes."""
    label = str(label)
    if label.endswith("_A"):
        return label[:-2]
    if label.endswith("_0A"):
        return label[:-1]
    return label


def _dedupe_2d_labels(labels: List[str], coords: Dict[str, List[float]],
                      axis: int, tol: float) -> List[str]:
    """Keep one representative for each physical 2D label/coordinate."""
    selected = []
    seen = set()
    for label in labels:
        if label not in coords or not _on_plane(coords[label], axis, tol):
            continue
        point = _project_to_plane(coords[label], axis)
        coord_key = tuple(round(float(x) % 1.0, 10) for x in point)
        key = (_base_2d_label(label), coord_key)
        if key in seen:
            continue
        seen.add(key)
        selected.append(label)
    return selected


def _ordered_2d_polygon_frac(coords: Dict[str, List[float]],
                             labels: List[str],
                             axis: int,
                             tol: float):
    in_plane_axes = [i for i in range(3) if i != axis]
    points = []
    for label in labels:
        if label in coords and _on_plane(coords[label], axis, tol):
            projected = _project_to_plane(coords[label], axis)
            points.append(projected)
    unique = []
    for point in points:
        if not any(np.allclose(point, old, atol=1e-8) for old in unique):
            unique.append(point)
    if len(unique) < 3:
        return unique
    pts = np.array(unique, dtype=float)
    pts2 = pts[:, in_plane_axes]
    if ConvexHull is None:
        center = pts2.mean(axis=0)
        angles = np.arctan2(pts2[:, 1] - center[1], pts2[:, 0] - center[0])
        return pts[np.argsort(angles)].tolist()
    hull = ConvexHull(pts2)
    return pts[hull.vertices].tolist()


def _ordered_2d_polygon_labels(coords: Dict[str, List[float]],
                               labels: List[str],
                               axis: int,
                               tol: float) -> List[str]:
    in_plane_axes = [i for i in range(3) if i != axis]
    entries = []
    for label in labels:
        if label in coords and _on_plane(coords[label], axis, tol):
            projected = _project_to_plane(coords[label], axis)
            entries.append((label, projected))
    unique = []
    for label, point in entries:
        if not any(np.allclose(point, old_point, atol=tol) for _old_label, old_point in unique):
            unique.append((label, point))
    if len(unique) < 3:
        return [label for label, _point in unique]
    labels_unique = [label for label, _point in unique]
    pts = np.array([point for _label, point in unique], dtype=float)
    pts2 = pts[:, in_plane_axes]
    if ConvexHull is None:
        center = pts2.mean(axis=0)
        angles = np.arctan2(pts2[:, 1] - center[1], pts2[:, 0] - center[0])
        return [labels_unique[i] for i in np.argsort(angles)]
    hull = ConvexHull(pts2)
    return [labels_unique[i] for i in hull.vertices]


def _closed_path_from_labels(labels: List[str]):
    if len(labels) < 2:
        return []
    return [
        (labels[i], labels[(i + 1) % len(labels)])
        for i in range(len(labels))
    ]


def _clip_polygon_halfplane(poly: np.ndarray, normal, offset: float, tol=1e-10):
    """Clip a 2D polygon to normal.x >= offset."""
    normal = np.array(normal, dtype=float)
    poly = np.array(poly, dtype=float)
    if len(poly) == 0:
        return poly
    out = []
    for idx, cur in enumerate(poly):
        prev = poly[idx - 1]
        cur_val = float(np.dot(normal, cur) - offset)
        prev_val = float(np.dot(normal, prev) - offset)
        cur_in = cur_val >= -tol
        prev_in = prev_val >= -tol
        denom = float(np.dot(normal, cur - prev))
        if cur_in:
            if not prev_in and abs(denom) > tol:
                t = (offset - np.dot(normal, prev)) / denom
                out.append(prev + t * (cur - prev))
            out.append(cur)
        elif prev_in and abs(denom) > tol:
            t = (offset - np.dot(normal, prev)) / denom
            out.append(prev + t * (cur - prev))
    return np.array(out, dtype=float)


def _cart2_to_frac(point_2d, b_matrix, axis: int, basis):
    in_plane_axes = [i for i in range(3) if i != axis]
    b2 = np.array([
        _to_2d(np.array(b_matrix[in_plane_axes[0]], dtype=float), basis),
        _to_2d(np.array(b_matrix[in_plane_axes[1]], dtype=float), basis),
    ])
    frac2 = np.array(point_2d, dtype=float) @ np.linalg.inv(b2)
    frac = [0.0, 0.0, 0.0]
    frac[in_plane_axes[0]] = float(frac2[0])
    frac[in_plane_axes[1]] = float(frac2[1])
    frac[axis] = 0.0
    return frac


def _in_plane_cart_ops(result: dict, axis: int, basis):
    e1, e2, _ = basis
    ops = []
    for g_cart in result.get("unique_ops") or []:
        g_arr = np.array(g_cart, dtype=float)
        ops.append(np.column_stack([
            _to_2d(g_arr @ e1, basis),
            _to_2d(g_arr @ e2, basis),
        ]))
    return ops


def _apply_geometric_2d_ibz(result: dict, axis: int, tol=1e-8) -> bool:
    """Build a layer-plane IBZ by clipping the actual 2D BZ with mirror lines."""
    if ConvexHull is None or not result.get("unique_ops"):
        return False
    b_matrix = np.array(result["b_matrix"], dtype=float)
    try:
        bz_poly, basis = _bz_polygon_2d(b_matrix, axis)
    except Exception:
        return False

    seed = np.mean(bz_poly, axis=0) + np.array([
        0.173 * max(np.ptp(bz_poly[:, 0]), tol),
        0.219 * max(np.ptp(bz_poly[:, 1]), tol),
    ])
    ibz_poly = np.array(bz_poly, dtype=float)
    used_boundaries = 0
    ident = np.eye(2)
    for op in _in_plane_cart_ops(result, axis, basis):
        if np.allclose(op, ident, atol=tol):
            continue
        fixed = op - ident
        if np.linalg.matrix_rank(fixed, tol=tol) != 1:
            continue
        normal = fixed[np.argmax(np.linalg.norm(fixed, axis=1))]
        norm = np.linalg.norm(normal)
        if norm < tol:
            continue
        normal = normal / norm
        if np.dot(normal, seed) < 0:
            normal = -normal
        clipped = _clip_polygon_halfplane(ibz_poly, normal, 0.0, tol=tol)
        if len(clipped) >= 3 and abs(ConvexHull(clipped).volume) > tol:
            ibz_poly = clipped
            used_boundaries += 1

    if used_boundaries == 0 or len(ibz_poly) < 3:
        return False
    center = np.mean(ibz_poly, axis=0)
    order = np.argsort(np.arctan2(ibz_poly[:, 1] - center[1],
                                  ibz_poly[:, 0] - center[0]))
    ibz_poly = ibz_poly[order]

    origin_idx = int(np.argmin(np.linalg.norm(ibz_poly, axis=1)))
    ibz_poly = np.roll(ibz_poly, -origin_idx, axis=0)
    coords = {}
    labels = []
    for idx, point in enumerate(ibz_poly):
        label = "GAMMA" if idx == 0 and np.linalg.norm(point) < 1e-7 else f"V{idx}"
        labels.append(label)
        coords[label] = _cart2_to_frac(point, b_matrix, axis, basis)
    path = _closed_path_from_labels(labels)
    centroid_2d = _polygon_centroid_2d(ibz_poly)
    centroid = _cart2_to_frac(centroid_2d, b_matrix, axis, basis)

    result["band_kpath"] = path
    result["band_kpoints_frac"] = coords
    result["centroid_frac"] = centroid
    result["centroid_cart"] = np.array(centroid) @ b_matrix
    result["ibz_polygon_frac"] = [coords[label] for label in labels]
    result["extra_general_vertices"] = []
    result["coordinate_basis_note"] = (
        "2D slab output uses the input reciprocal basis. The 2D IBZ was "
        "built geometrically from the input-kz=0 BZ and in-plane symmetry "
        "boundaries, not from projected 3D HPKOT labels."
    )
    return True


def _build_2d_path_result(struct_file: str, axis: int, tol: float):
    result = compute_centroid(
        struct_file, output_dir=".", show_plot=False, defer_show=False, verbose=False
    )
    result = _apply_2d_monoclinic_seekpath_basis(result)
    result = _convert_result_to_input_slab_basis(result)
    axis = VACUUM_AXIS
    result["vacuum_axis_detection"] = {
        "method": "input_slab_c_axis",
        "axis": int(axis),
        "reciprocal_norms": [
            float(x) for x in np.linalg.norm(np.array(result["b_matrix"], dtype=float), axis=1)
        ],
    }
    result["dimension"] = 2
    result["vacuum_axis"] = axis
    result["vacuum_axis_name"] = ("a", "b", "c")[axis]
    basename = os.path.basename(struct_file)
    for key in ("mC1", "mC2", "mC3", "mP"):
        if key in basename:
            result["sc_type"] = key
            break
    if _is_monoclinic_result(result) and _apply_geometric_2d_ibz(result, axis, tol=tol):
        return result
    coords = result.get(
        "band_kpoints_frac",
        result.get("path_kpoints_frac", result.get("ibz_kpoints_frac", {})),
    )
    section_coords_source = result.get(
        "ibz_kpoints_frac",
        result.get("kpoints_frac", coords),
    )
    path = result.get("band_kpath", result.get("ibz_kpath", result.get("sp_path", [])))
    plane_path = _dedupe_segments(path, coords, axis, tol)
    labels = sorted({label for segment in plane_path for label in segment})
    section_labels = sorted(
        label for label, point in section_coords_source.items()
        if _on_plane(point, axis, tol) and not str(label).startswith("_")
    )
    section_coords_all = dict(section_coords_source)
    section_coords_all.update(coords)
    polygon_labels = _ordered_2d_polygon_labels(
        section_coords_all, section_labels, axis, tol)
    if len(polygon_labels) >= 3:
        polygon_edges = {
            frozenset(edge) for edge in _closed_path_from_labels(polygon_labels)
        }
        plane_edges = {frozenset(edge) for edge in plane_path}
        if _is_monoclinic_result(result) and not polygon_edges.issubset(plane_edges):
            plane_path = _closed_path_from_labels(polygon_labels)
            labels = polygon_labels
    if len(labels) < 2:
        labels = section_labels
        plane_path = []
    if len(labels) < 2:
        raise RuntimeError(
            f"Only {len(labels)} high-symmetry label(s) remain on {_axis_name(axis)}=0. "
            "This structure may need a manual 2D path."
        )
    # In 2D mode, do not inherit the project-specific doubled 3D IBZ hull
    # (labels such as M_A/L_A/R_A).  Those copied vertices are useful as 3D
    # anchors, but in a true 2D BZ they duplicate equivalent boundaries and can
    # place the general k point on a spin-flip mirror line.
    centroid_labels = _dedupe_2d_labels(labels, section_coords_all, axis, tol)
    centroid_coords = section_coords_all
    centroid = _compute_2d_centroid(centroid_coords, centroid_labels, axis, tol)
    ibz_polygon_frac = _ordered_2d_polygon_frac(
        centroid_coords, centroid_labels, axis, tol)
    result["band_kpath"] = plane_path
    section_coords = {
        label: _project_to_plane(section_coords_all[label], axis)
        for label in labels
        if label in section_coords_all
    }
    result["band_kpoints_frac"] = section_coords
    result["centroid_frac"] = centroid
    result["centroid_cart"] = np.array(centroid) @ np.array(result["b_matrix"])
    result["ibz_polygon_frac"] = ibz_polygon_frac
    result["extra_general_vertices"] = []
    return result


def _prompt_structure_and_spin(BOLD: str, RESET: str):
    print(f"\n{BOLD}>>> Step 0: Spin symmetry{RESET}")
    print("Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): ", end="", flush=True)
    struct_file = input().strip() or "POSCAR"
    step0_wrote_flip_file = None
    actual_flip_ops_for_plot = None
    spin_operation_lattice = None
    flip_file = "spin_flip_operations.txt"

    if not os.path.exists(struct_file):
        print(f"[Note] '{struct_file}' not found. Skipping Step 0, using existing spin_flip_operations.txt")
        return None, step0_wrote_flip_file, actual_flip_ops_for_plot, spin_operation_lattice
    if not FIND_SF_AVAILABLE:
        print("[Note] find_sf_operations.py not found. Skipping Step 0.")
        return struct_file, step0_wrote_flip_file, actual_flip_ops_for_plot, spin_operation_lattice

    if struct_file.lower().endswith(".mcif"):
        print("Detected .mcif file --magnetic moments will be read from file.")
        moments_str = ""
    else:
        print("Magnetic moments (atom order, trailing atoms auto-fill to 0): ", end="", flush=True)
        moments_str = input().strip()

    mtime_before = os.path.getmtime(flip_file) if os.path.exists(flip_file) else None
    sf_result = find_sf_run(struct_file, moments_str, verbose=False)
    mtime_after = os.path.getmtime(flip_file) if os.path.exists(flip_file) else None
    step0_wrote_flip_file = mtime_after is not None and mtime_after != mtime_before

    if isinstance(sf_result, dict):
        step0_wrote_flip_file = sf_result.get("spin_flip_operations", 0) > 0
        spin_operation_lattice = sf_result.get("operation_lattice")
        print(f"\nStructure: {sf_result['structure_file']}, atoms: {sf_result['num_atoms']}")
        print(f"Space Group: {sf_result['space_group']}")
        print(f"Point Group: {sf_result['point_group']}")
        print(f"Laue Group: {sf_result['laue_group']}")
        print(f"Magnetic SG: {sf_result['magnetic_space_group']}")
        print(f"Spin group: {sf_result['spin_group']}")
        print("Actual point operations: "
              f"{sf_result['actual_spin_flip_point_operations']} spin-flip, "
              f"{sf_result['actual_spin_preserve_point_operations']} spin-preserving")
        print("Inversion-extended k operations: "
              f"{sf_result['extended_spin_flip_point_operations']} spin-flip, "
              f"{sf_result['extended_spin_preserve_point_operations']} spin-preserving")
        if sf_result["spin_split_diagnostic"]:
            print(f"{BOLD}Warning! {sf_result['spin_split_diagnostic']}{RESET}")
        actual_flip_ops = sf_result.get("actual_spin_flip_point_matrices", [])
        actual_flip_ops_for_plot = actual_flip_ops
        diagnostic_axis = _default_2d_axis_from_symmetry(sf_result)
        forbidden_name, _ = _forbidden_2d_altermagnetic_spin_flip_operation(
            actual_flip_ops, axis=diagnostic_axis
        )
        if forbidden_name is not None:
            print(
                f"{BOLD}Warning! The actual spin-flip point operations contain "
                f"{forbidden_name}, which a 2D altermagnet should break.{RESET}"
            )
        print(f"Saved: {', '.join(sf_result['saved_files'])}")

    return struct_file, step0_wrote_flip_file, actual_flip_ops_for_plot, spin_operation_lattice


def _choose_spin_flip_operation(
    modifier: KPointsModifier,
    step0_wrote_flip_file: Optional[bool],
    centroid_result: dict,
):
    axis = int(centroid_result.get("vacuum_axis", VACUUM_AXIS))
    if step0_wrote_flip_file is False:
        print("[Note] Step 0 found no spin-flip operations --skipping any existing spin_flip_operations.txt.")
        flip_ops = []
    else:
        flip_ops = modifier.load_flip_operations()

    R = None
    label = None
    if flip_ops:
        print(f"Found {len(flip_ops)} spin-flip operations.")
        print("Default R: Option 1")
        while R is None:
            print("Press [Enter] to use default, type a number, 'list' to show matrices, or 'manual': ",
                  end="", flush=True)
            choice = input().strip().lower()
            if not choice:
                R = flip_ops[0]
                label = "Option 1"
                print("Selected: Option 1")
            elif choice == "list":
                for i, op in enumerate(flip_ops):
                    print(f"\n  Option {i + 1}:")
                    print(modifier._format_matrix(op))
                print()
            elif choice == "manual":
                break
            else:
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(flip_ops):
                        R = flip_ops[idx]
                        label = f"Option {idx + 1}"
                        print(f"Selected: Option {idx + 1}")
                    else:
                        print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")
                except ValueError:
                    print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")

    if R is None:
        print("Enter custom transformation matrix R.")
        print("Enter row by row (3 numbers per row, space-separated):")
        rows = []
        for i in range(3):
            while True:
                try:
                    print(f"Row {i + 1}: ", end="", flush=True)
                    values = input().strip().split()
                    if len(values) == 3:
                        rows.append([float(x) for x in values])
                        break
                    print("Please enter exactly 3 numbers.")
                except ValueError:
                    print("Invalid input.")
        R = np.array(rows)
        label = "manual"

    return np.array(R, dtype=float), label, flip_ops


def _convert_input_R_to_primitive(centroid_result: dict, R, source_lattice=None):
    if centroid_result.get("coordinate_basis") == "input_slab":
        return np.array(R, dtype=float)
    if source_lattice is not None:
        converted, _ = convert_real_frac_op_to_primitive_k_basis(
            R,
            source_lattice,
            centroid_result["b_matrix"],
            centroid_result.get("seekpath_rotation_matrix"),
        )
        return converted
    b_input = np.array(
        centroid_result.get("b_matrix_input", centroid_result["b_matrix_conv"]),
        dtype=float,
    )
    b_prim = np.array(centroid_result["b_matrix"], dtype=float)
    R_arr = np.array(R, dtype=float)
    R_cart = b_input.T @ np.linalg.inv(R_arr).T @ np.linalg.inv(b_input.T)
    R_prim_inv_T = np.linalg.inv(b_prim.T) @ R_cart @ b_prim.T
    return np.linalg.inv(R_prim_inv_T.T)


def _cart_from_frac(frac, b_matrix):
    f = np.array(frac[:3], dtype=float)
    return f[0] * b_matrix[0] + f[1] * b_matrix[1] + f[2] * b_matrix[2]


def _plane_projector(b_matrix, axis: int):
    in_plane_axes = [i for i in range(3) if i != axis]
    b_matrix = np.array(b_matrix, dtype=float)
    g1 = b_matrix[in_plane_axes[0]]
    g2 = b_matrix[in_plane_axes[1]]
    e1_norm = np.linalg.norm(g1)
    if e1_norm < 1e-14:
        raise ValueError("invalid zero-length in-plane reciprocal vector")
    e1 = g1 / e1_norm
    g2_perp = g2 - np.dot(g2, e1) * e1
    e2_norm = np.linalg.norm(g2_perp)
    if e2_norm < 1e-14:
        raise ValueError("collinear in-plane reciprocal vectors")
    e2 = g2_perp / e2_norm
    return (e1, e2, in_plane_axes)


def _use_cartesian_xy_plot(result: dict) -> bool:
    """For 2D figures, keep screen axes aligned with Cartesian kx, ky."""
    return int(result.get("vacuum_axis", VACUUM_AXIS)) == 2


def _to_2d(point, basis):
    """Project Cartesian k-space point to the selected 2D reciprocal plane."""
    point = np.array(point, dtype=float)
    if basis is not None:
        e1, e2, _ = basis
        return np.array([np.dot(point, e1), np.dot(point, e2)], dtype=float)
    return np.array([point[0], point[1]], dtype=float)


def _bz_polygon_2d(b_matrix, axis: int, radius=2, cartesian_xy=False):
    """Return the 2D Wigner-Seitz BZ polygon for the selected reciprocal plane."""
    if ConvexHull is None:
        raise RuntimeError("scipy is required for 2D BZ polygon construction")
    if cartesian_xy:
        if axis != 2:
            raise ValueError("Cartesian x/y 2D plotting requires the kz=0 plane")
        basis = None
        in_plane_axes = [0, 1]
    else:
        basis = _plane_projector(b_matrix, axis)
        _, _, in_plane_axes = basis
    vectors = []
    for i in range(-radius, radius + 1):
        for j in range(-radius, radius + 1):
            if i == 0 and j == 0:
                continue
            g3 = i * np.array(b_matrix[in_plane_axes[0]]) + j * np.array(b_matrix[in_plane_axes[1]])
            vectors.append(_to_2d(g3, basis))

    # Half-plane intersection by clipping a large square against |k| <= |k-G|.
    span = max(np.linalg.norm(v) for v in vectors) * 2.0
    poly = np.array([
        [-span, -span],
        [ span, -span],
        [ span,  span],
        [-span,  span],
    ], dtype=float)

    def clip_polygon(poly_pts, normal, offset):
        if len(poly_pts) == 0:
            return poly_pts
        out = []
        for idx, cur in enumerate(poly_pts):
            prev = poly_pts[idx - 1]
            cur_in = np.dot(normal, cur) <= offset + 1e-10
            prev_in = np.dot(normal, prev) <= offset + 1e-10
            denom = np.dot(normal, cur - prev)
            if cur_in:
                if not prev_in and abs(denom) > 1e-14:
                    t = (offset - np.dot(normal, prev)) / denom
                    out.append(prev + t * (cur - prev))
                out.append(cur)
            elif prev_in and abs(denom) > 1e-14:
                t = (offset - np.dot(normal, prev)) / denom
                out.append(prev + t * (cur - prev))
        return np.array(out, dtype=float)

    for g in sorted(vectors, key=lambda v: np.linalg.norm(v)):
        poly = clip_polygon(poly, g, 0.5 * np.dot(g, g))
        if len(poly) == 0:
            raise RuntimeError("2D BZ polygon clipping failed")

    hull = ConvexHull(poly)
    poly = poly[hull.vertices]
    center = np.mean(poly, axis=0)
    angles = np.arctan2(poly[:, 1] - center[1], poly[:, 0] - center[0])
    return poly[np.argsort(angles)], basis


def _points_on_frac_axis_plane(points, simplices, b_matrix, axis: int,
                               basis, value=0.0, tol=1e-8):
    """Return a 2D convex section of a cartesian hull at fractional k_axis=value."""
    points = np.array(points, dtype=float)
    b_inv = np.linalg.inv(np.array(b_matrix, dtype=float))
    frac = points @ b_inv
    section = []
    for simplex in np.array(simplices, dtype=int):
        tri_cart = points[simplex]
        tri_frac = frac[simplex]
        vals = tri_frac[:, axis] - value
        for i in range(3):
            j = (i + 1) % 3
            vi, vj = vals[i], vals[j]
            if abs(vi) <= tol:
                section.append(tri_cart[i])
            if vi * vj < -tol * tol:
                t = -vi / (vj - vi)
                section.append(tri_cart[i] + t * (tri_cart[j] - tri_cart[i]))
            elif abs(vj) <= tol:
                section.append(tri_cart[j])
    if len(section) < 3:
        return None
    pts2 = np.unique(
        np.round(np.array([_to_2d(point, basis) for point in section]), 10),
        axis=0,
    )
    if len(pts2) < 3:
        return None
    try:
        hull = ConvexHull(pts2)
        return pts2[hull.vertices]
    except Exception:
        center = pts2.mean(axis=0)
        angles = np.arctan2(pts2[:, 1] - center[1], pts2[:, 0] - center[0])
        ordered = pts2[np.argsort(angles)]
        return ordered if len(ordered) >= 3 else None


def _setup_2d_ax(title, bz_poly):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))
    closed = np.vstack([bz_poly, bz_poly[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="0.25", lw=2.2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    ax.set_title(title, fontsize=20)
    span = np.ptp(bz_poly, axis=0)
    pad = 0.22 * max(span.max(), 1e-8)
    ax.set_xlim(np.min(bz_poly[:, 0]) - pad, np.max(bz_poly[:, 0]) + pad)
    ax.set_ylim(np.min(bz_poly[:, 1]) - pad, np.max(bz_poly[:, 1]) + pad)
    return fig, ax


def _label_offset(points):
    pts = np.array(list(points), dtype=float)
    center = np.mean(pts, axis=0)
    span = max(np.max(np.ptp(pts, axis=0)), 1e-8)
    return center, 0.12 * span


def _draw_labeled_points(ax, points, color, edgecolor, labels=True,
                         prime=False, label_color=None):
    if not points:
        return
    center, offset_scale = _label_offset(points.values())
    for label, point in points.items():
        ax.scatter(point[0], point[1], s=85, c=color, edgecolors=edgecolor,
                   linewidths=0.7, zorder=5)
        if not labels:
            continue
        direction = point - center
        norm = np.linalg.norm(direction)
        offset = direction / norm * offset_scale if norm > 1e-8 else np.array([offset_scale, 0.0])
        text = _display_label(label.rstrip("'"), prime=prime)
        ax.text(*(point + offset), text, fontsize=24,
                color=label_color if label_color is not None else edgecolor,
                ha="center", va="center", zorder=6)


def _plot_spin_pattern_top_view_2d(centroid_result, R_for_kpts,
                                   flip_ops_for_plot,
                                   output_path, show_title=False,
                                   close_figure=True):
    """Color all 2D symmetry images of the IBZ by spin-preserving/flipping ops."""
    import matplotlib.pyplot as plt
    axis = int(centroid_result.get("vacuum_axis", VACUUM_AXIS))

    unique_ops = centroid_result.get("unique_ops")
    ibz_polygon_frac = centroid_result.get("ibz_polygon_frac")
    if not ibz_polygon_frac or len(ibz_polygon_frac) < 3 or not unique_ops:
        print("[Note] Skipping 2D spin-pattern figure (no 2D IBZ polygon or symmetry ops available).")
        return None

    b_matrix = np.array(centroid_result["b_matrix"], dtype=float)
    cartesian_xy = _use_cartesian_xy_plot(centroid_result)
    bz_poly, basis = _bz_polygon_2d(b_matrix, axis, cartesian_xy=cartesian_xy)
    centroid_frac = np.array(centroid_result["centroid_frac"], dtype=float)
    centroid_frac[axis] = 0.0
    centroid_cart = _cart_from_frac(centroid_frac, b_matrix)

    spin_down_mask = _classify_spin_down_ops_2d(
        b_matrix, unique_ops, R_for_kpts,
        flip_ops_frac=flip_ops_for_plot if flip_ops_for_plot else None,
        axis=axis,
    )
    unique_ops, spin_down_mask = _dedupe_in_plane_ops_2d(
        b_matrix, unique_ops, spin_down_mask, axis=axis
    )
    if spin_down_mask.sum() in {0, len(spin_down_mask)}:
        print(
            "[Warning] 2D spin-BZ coloring did not split the symmetry images. "
            "Check the 3D spin-flip operation list and operation basis."
            )
    if flip_ops_for_plot:
        identity_like = [
            op for op in flip_ops_for_plot
            if _operation_identity_on_plane(op, axis)
        ]
        if identity_like:
            print(
                "[Warning] A 3D spin-flip operation leaves the selected 2D "
                f"plane {_axis_name(axis)}=0 unchanged. This protects "
                "spin degeneracy at the same in-plane k and is not a "
                "nontrivial 2D altermagnetic splitting pattern."
            )
    else:
        identity_like = []

    fig, ax = plt.subplots(figsize=(9, 9))
    fill_alpha = 0.68
    up_labeled = False
    down_labeled = False
    if identity_like or spin_down_mask.sum() == 0:
        ax.fill(bz_poly[:, 0], bz_poly[:, 1], facecolor="#b22222",
                alpha=fill_alpha, edgecolor="none", label="spin-up")
        up_labeled = True
    elif spin_down_mask.sum() == len(spin_down_mask):
        ax.fill(bz_poly[:, 0], bz_poly[:, 1], facecolor="#1f4e9e",
                alpha=fill_alpha, edgecolor="none", label="spin-down")
        down_labeled = True
    ibz_polygon_cart = np.array([
        _cart_from_frac(point, b_matrix)
        for point in ibz_polygon_frac
    ], dtype=float)
    if not (identity_like or spin_down_mask.sum() in {0, len(spin_down_mask)}):
        for i, g in enumerate(unique_ops):
            cell_pts = (g @ ibz_polygon_cart.T).T
            poly = np.array([_to_2d(point, basis) for point in cell_pts], dtype=float)
            is_down = spin_down_mask[i]
            color = "#1f4e9e" if is_down else "#b22222"
            if is_down:
                label = "spin-down" if not down_labeled else None
                down_labeled = True
            else:
                label = "spin-up" if not up_labeled else None
                up_labeled = True
            closed = np.vstack([poly, poly[0]])
            ax.fill(poly[:, 0], poly[:, 1], facecolor=color, alpha=fill_alpha,
                    edgecolor="none", label=label)
            ax.plot(closed[:, 0], closed[:, 1], color=color, lw=0.9, alpha=0.95)

    closed = np.vstack([bz_poly, bz_poly[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="black", lw=2.0, label="BZ boundary")

    ax.set_aspect("equal", adjustable="box")
    if show_title:
        ax.set_title(r"Spin-up / Spin-down BZ top view ($k_z=0$)", fontsize=18)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    ax.legend(loc="upper right", fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    if close_figure:
        plt.close(fig)
    return output_path


def _display_label(label, prime=False):
    label = str(label)
    greek = {
        "GAMMA": r"\Gamma",
        "\u0393": r"\Gamma",
        "DELTA": r"\Delta",
        "\u0394": r"\Delta",
        "LAMBDA": r"\Lambda",
        "\u039b": r"\Lambda",
        "SIGMA": r"\Sigma",
        "\u03a3": r"\Sigma",
    }
    prime_text = "'" if prime else ""
    if "_" in label:
        base, sub = label.split("_", 1)
        math_base = greek.get(base.strip().upper(), base)
        return rf"${math_base}_{{{sub}}}{prime_text}$"
    math_label = greek.get(label.strip().upper())
    if math_label is not None:
        return rf"${math_label}{prime_text}$"
    text = KPointsModifier._display_label(label)
    if "_" in text:
        base, sub = text.split("_", 1)
        text = rf"${base}_{{{sub}}}$"
    return f"{text}'" if prime else text


def _plot_2d_figures(centroid_result, modifier, general_kpoint,
                     R_for_kpts, flip_ops_for_plot, path_sequence, basename):
    """Save the 2D Figure 1/2/3 stack for the selected reciprocal plane."""
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[Warning] Could not import matplotlib for 2D figures: {exc}")
        return []

    b_matrix = np.array(centroid_result["b_matrix"], dtype=float)
    axis = int(centroid_result.get("vacuum_axis", VACUUM_AXIS))
    bz_poly, basis = _bz_polygon_2d(
        b_matrix, axis, cartesian_xy=_use_cartesian_xy_plot(centroid_result)
    )
    sc_type = centroid_result.get("sc_type", "BZ")
    kpath = centroid_result["band_kpath"]
    coords = centroid_result["band_kpoints_frac"]
    kpoints_cart = {
        label: _to_2d(_cart_from_frac(frac, b_matrix), basis)
        for label, frac in coords.items()
        if not str(label).startswith("_")
    }
    centroid_cart3 = _cart_from_frac(general_kpoint, b_matrix)
    centroid_xy = _to_2d(centroid_cart3, basis)
    R_inv_T = np.linalg.inv(np.array(R_for_kpts, dtype=float)).T
    k_prime = R_inv_T @ np.array(general_kpoint, dtype=float)
    k_prime[axis] = 0.0
    k_prime_xy = _to_2d(_cart_from_frac(k_prime, b_matrix), basis)

    saved = []

    # Figure 1: true 2D BZ + in-plane path + centroid.
    fig1, ax1 = _setup_2d_ax(f"2D BZ: {basename} ({sc_type})", bz_poly)
    for start, end in kpath:
        if start in kpoints_cart and end in kpoints_cart:
            p1, p2 = kpoints_cart[start], kpoints_cart[end]
            ax1.plot([p1[0], p2[0]], [p1[1], p2[1]],
                     c="red", lw=3.0, alpha=0.9, zorder=3)
    _draw_labeled_points(ax1, kpoints_cart, "red", "darkred",
                         label_color="black")
    ax1.scatter(*centroid_xy, c="gold", marker="*", s=420,
                edgecolors="k", zorder=112, label="2D centroid")
    ax1.legend(loc="upper right")
    fig1.tight_layout()
    fig1_path = f"{basename}_2d_ibz_{sc_type}.png"
    fig1.savefig(fig1_path, dpi=300, bbox_inches="tight")
    saved.append(fig1_path)

    mapped_cart = {}
    mapped_cart_lines = {}
    for label, frac in coords.items():
        if str(label).startswith("_"):
            continue
        mapped_frac = R_inv_T @ np.array(frac, dtype=float)
        mapped_frac[axis] = 0.0
        mapped_point = _to_2d(_cart_from_frac(mapped_frac, b_matrix), basis)
        mapped_cart_lines[label + "'"] = mapped_point
        if not any(np.allclose(mapped_point, orig, atol=1e-8)
                   for orig in kpoints_cart.values()):
            mapped_cart[label + "'"] = mapped_point

    # Figure 2: 2D mapped spin-flip figure.
    fig2, ax2 = _setup_2d_ax("2D spin-flip path connections", bz_poly)
    orig_poly = np.array(list(kpoints_cart.values()), dtype=float)
    mapped_poly = np.array(list(mapped_cart_lines.values()), dtype=float)
    if len(orig_poly) >= 3 and ConvexHull is not None:
        hull = ConvexHull(orig_poly)
        hp = orig_poly[hull.vertices]
        ax2.fill(hp[:, 0], hp[:, 1], color="salmon", alpha=0.20, zorder=1)
    if len(mapped_poly) >= 3 and ConvexHull is not None:
        hull = ConvexHull(mapped_poly)
        hp = mapped_poly[hull.vertices]
        ax2.fill(hp[:, 0], hp[:, 1], color="cornflowerblue", alpha=0.20, zorder=1)
    for start, end in kpath:
        if start in kpoints_cart and end in kpoints_cart:
            p1, p2 = kpoints_cart[start], kpoints_cart[end]
            ax2.plot([p1[0], p2[0]], [p1[1], p2[1]],
                     c="red", lw=4.0, alpha=0.9, zorder=50)
    for start, end in kpath:
        sp = start + "'"
        ep = end + "'"
        if sp in mapped_cart_lines and ep in mapped_cart_lines:
            p1, p2 = mapped_cart_lines[sp], mapped_cart_lines[ep]
            ax2.plot([p1[0], p2[0]], [p1[1], p2[1]],
                     c="navy", lw=4.0, alpha=0.9, zorder=50)
    _draw_labeled_points(ax2, kpoints_cart, "salmon", "darkred", prime=False)
    _draw_labeled_points(ax2, mapped_cart, "cornflowerblue", "navy", prime=True)
    threshold = 0.05 * max(np.max(np.linalg.norm(bz_poly, axis=1)), 1e-8)
    for point in kpoints_cart.values():
        if np.linalg.norm(point - centroid_xy) > threshold:
            ax2.plot([point[0], centroid_xy[0]], [point[1], centroid_xy[1]],
                     c="deepskyblue",
                     lw=2.0, ls="--", alpha=0.75, zorder=40)
    for point in mapped_cart_lines.values():
        if np.linalg.norm(point - k_prime_xy) > threshold:
            ax2.plot([point[0], k_prime_xy[0]], [point[1], k_prime_xy[1]],
                     c="deepskyblue",
                     lw=2.0, ls="--", alpha=0.75, zorder=40)

    ax2.scatter(*centroid_xy, c="gold", s=300, marker="*",
                edgecolors="k", linewidths=0.8, zorder=120, label=r"$k$")
    ax2.scatter(*k_prime_xy, c="cornflowerblue", s=150, marker="o",
                edgecolors="k", linewidths=0.8, zorder=120, label=r"$k'$")
    ax2.legend(loc="upper right", fontsize=18)
    fig2.tight_layout()
    fig2_path = f"{basename}_2d_spinflip_{sc_type}.png"
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
    saved.append(fig2_path)

    fig3_path = f"{basename}_2d_spinbz_{sc_type}.png"
    fig3_saved = _plot_spin_pattern_top_view_2d(
        centroid_result, R_for_kpts, flip_ops_for_plot, fig3_path,
        show_title=False, close_figure=False
    )
    if fig3_saved is not None:
        saved.append(fig3_saved)

    print(f"Saved 2D figures: {', '.join(saved)}")
    plt.show()
    plt.close("all")
    return saved


def interactive_2d():
    BOLD = "\033[1m"
    RESET = "\033[0m"
    print("=== Altermagnetic 2D K-Path Generator ===")

    (
        struct_file,
        step0_wrote_flip_file,
        actual_flip_ops_for_plot,
        spin_operation_lattice,
    ) = _prompt_structure_and_spin(BOLD, RESET)
    if struct_file is None:
        print("A structure file is required for the automatic 2D path.")
        return

    print(f"\n{BOLD}>>> Step 1: 2D high-symmetry k-path{RESET}")
    axis = VACUUM_AXIS
    tol = 1e-6
    centroid_result = _build_2d_path_result(struct_file, axis, tol)
    axis = int(centroid_result.get("vacuum_axis", axis))
    modifier = KPointsModifier()
    modifier.header_lines = [
        f"K-Path generated by AlterSeeK-Path 2D ({_axis_name(axis)}=0, HPKOT {centroid_result.get('sc_type', 'unknown')})",
        "20",
        "Line-Mode",
        "Reciprocal",
    ]

    coords = centroid_result["band_kpoints_frac"]
    path = centroid_result["band_kpath"]
    print(f"Source: {centroid_result.get('sc_type', 'unknown')} | Plane: {_axis_name(axis)} = 0")
    if path:
        print(f"2D path: {modifier._format_path(path)}")
    else:
        print("[Warning] No automatic in-plane path was found.")
    print("Press [Enter] to use this 2D path, or type a filename to load your own: ",
          end="", flush=True)
    path_choice = input().strip()
    if path_choice:
        if not modifier.read_kpoints_file(path_choice):
            return
        try:
            path, coords = _apply_manual_2d_path(centroid_result, modifier, axis, tol)
        except Exception as exc:
            print(f"[Error] Could not use manual 2D KPATH: {exc}")
            return
        modifier.header_lines = [
            f"Manual 2D KPATH ({_axis_name(axis)}=0, HPKOT {centroid_result.get('sc_type', 'unknown')})",
            "20",
            "Line-Mode",
            "Reciprocal",
        ]
        print(f"Using manual 2D path ({len(path)} segments, {len(modifier.kpoints_data)} k-points).")
    else:
        if not path:
            print("Automatic 2D path is empty. Please rerun and provide a manual KPATH file.")
            return
        modifier.kpoints_data = []
        for start, end in path:
            for label in (start, end):
                point = coords[label]
                modifier.kpoints_data.append([point[0], point[1], point[2], label])

    modifier.extra_general_points = []
    for label in centroid_result.get("extra_general_vertices", []):
        if label in coords:
            point = coords[label]
            modifier.extra_general_points.append([point[0], point[1], point[2], label])

    print(f"\n{BOLD}>>> Step 2: General k-point{RESET}")
    general_kpoint = centroid_result["centroid_frac"]
    print(f"2D IBZ centroid: [{general_kpoint[0]:.6f}, {general_kpoint[1]:.6f}, {general_kpoint[2]:.6f}]")
    try:
        with open("spin_operations.txt", "a") as f:
            f.write(
                "\nGeneral k-point (2D IBZ centroid, fractional): "
                f"[{general_kpoint[0]:.6f}, {general_kpoint[1]:.6f}, {general_kpoint[2]:.6f}]\n"
            )
    except Exception:
        pass

    print(f"\n{BOLD}>>> Step 3: Spin-flip operation{RESET}")
    R, selected_label, _flip_ops = _choose_spin_flip_operation(
        modifier, step0_wrote_flip_file, centroid_result
    )

    print(f"\n{BOLD}>>> Step 4: Build 2D altermagnetic path{RESET}")
    R_for_kpts = _convert_input_R_to_primitive(
        centroid_result, R, spin_operation_lattice)
    plot_flip_source = _flip_ops or actual_flip_ops_for_plot
    flip_ops_for_plot = [
        _convert_input_R_to_primitive(centroid_result, op, spin_operation_lattice)
        for op in plot_flip_source
    ] if plot_flip_source else None
    print("[Basis] Converted R from input-cell basis to primitive basis.")
    print("Primitive-basis R used for KPOINTS:")
    print(modifier._format_matrix(R_for_kpts))

    k_prime = modifier.transform_kpoint(general_kpoint, R_for_kpts)
    if abs(k_prime[axis]) > 1e-6:
        print(
            f"[Warning] k' has nonzero vacuum-axis coordinate ({k_prime[axis]:.6f}). "
            "Projecting k' back to the 2D plane for KPOINTS output."
        )
    k_prime[axis] = 0.0
    new_kpoints = modifier.insert_general_kpoints(
        general_kpoint, R_for_kpts, modifier.extra_general_points
    )
    for point in new_kpoints:
        if point is not None:
            point[axis] = 0.0

    print(f"\n{BOLD}>>> Step 5: Save{RESET}")
    print("Enter output filename (default: KPOINTS_2D): ", end="", flush=True)
    output_file = input().strip() or "KPOINTS_2D"
    modifier.write_kpoints_file(new_kpoints, output_file, R, selected_label)

    basename = os.path.splitext(os.path.basename(struct_file))[0] if struct_file else "POSCAR"
    _plot_2d_figures(
        centroid_result,
        modifier,
        general_kpoint,
        R_for_kpts,
        flip_ops_for_plot,
        new_kpoints,
        basename,
    )
    print("\nDone.")


def main():
    if len(sys.argv) > 1 and sys.argv[1] in {"-h", "--help"}:
        print("usage: python alterseek_path_2d.py")
        print("Run interactively to generate a 2D/slab altermagnetic KPOINTS path.")
        return
    interactive_2d()


if __name__ == "__main__":
    main()
