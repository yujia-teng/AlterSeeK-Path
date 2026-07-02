"""Top-down 2D Brillouin-zone figures for AlterSeeK-Path slab mode.

For a 2D material computed as a slab the physical reciprocal space is the
``k[vacuum_axis] = 0`` plane.  These figures render that plane straight-on as a
flat top-down view, instead of the tilted 3D BZ plate used for bulk crystals:

  * Figure 1 (``*_2d_ibz_*``)      -- 2D BZ outline + in-plane band path + k.
  * Figure 2 (``*_2d_spinflip_*``) -- spin-up IBZ (red) and its spin-flip image
                                      (blue) with the general point k and its
                                      spin-flip partner k'.

The vacuum axis is taken from ``centroid_result['vacuum_axis']`` (detected in
the standardized frame by ``compute_centroid_hybrid``).  For the common
``kz = 0`` slab (vacuum on axis c) the screen axes are the Cartesian kx, ky;
when seekpath permutes the axes (e.g. monoclinic -> unique axis b) the in-plane
reciprocal vectors define an orthonormal screen basis instead.

Ported from the original dev-branch ``alterseek_path_2d.py`` 2D figure stack.
"""

import os
import numpy as np

try:
    from scipy.spatial import ConvexHull
except Exception:  # pragma: no cover - scipy is a hard dependency in practice
    ConvexHull = None

try:
    from .symmetry import (
        _classify_spinflip_op, _reduce_int_vector, _format_miller,
        _rotation_sense,
    )
except Exception:  # pragma: no cover - op-visual is best-effort only
    _classify_spinflip_op = None


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

def _cart_from_frac(frac, b_matrix):
    f = np.array(frac[:3], dtype=float)
    return f[0] * b_matrix[0] + f[1] * b_matrix[1] + f[2] * b_matrix[2]


def _plane_projector(b_matrix, axis):
    """Orthonormal 2D screen basis spanning the in-plane reciprocal vectors."""
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


def _to_2d(point, basis):
    """Project a Cartesian k-space point onto the selected 2D plane."""
    point = np.array(point, dtype=float)
    if basis is not None:
        e1, e2, _ = basis
        return np.array([np.dot(point, e1), np.dot(point, e2)], dtype=float)
    return np.array([point[0], point[1]], dtype=float)


def _bz_polygon_2d(b_matrix, axis, radius=2, cartesian_xy=False):
    """2D Wigner-Seitz BZ polygon for the selected reciprocal plane."""
    if ConvexHull is None:
        raise RuntimeError("scipy is required for 2D BZ polygon construction")
    if cartesian_xy:
        if axis != 2:
            raise ValueError("Cartesian kx/ky 2D plotting requires the kz=0 plane")
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
            g = (i * np.array(b_matrix[in_plane_axes[0]])
                 + j * np.array(b_matrix[in_plane_axes[1]]))
            vectors.append(_to_2d(g, basis))

    span = max(np.linalg.norm(v) for v in vectors) * 2.0
    poly = np.array([[-span, -span], [span, -span],
                     [span, span], [-span, span]], dtype=float)

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


# ---------------------------------------------------------------------------
# drawing helpers
# ---------------------------------------------------------------------------

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


_GREEK = {
    "GAMMA": r"\Gamma", "Γ": r"\Gamma",
    "DELTA": r"\Delta", "Δ": r"\Delta",
    "LAMBDA": r"\Lambda", "Λ": r"\Lambda",
    "SIGMA": r"\Sigma", "Σ": r"\Sigma",
}


def _fig_label(label, prime=False):
    label = str(label).rstrip("'")
    pt = "'" if prime else ""
    if "_" in label:
        base, sub = label.split("_", 1)
        math_base = _GREEK.get(base.strip().upper(), base)
        return rf"${math_base}_{{{sub}}}{pt}$"
    math_label = _GREEK.get(label.strip().upper())
    if math_label is not None:
        return rf"${math_label}{pt}$"
    return f"{label}{pt}"


def _label_offset(points):
    pts = np.array(list(points), dtype=float)
    center = np.mean(pts, axis=0)
    span = max(np.max(np.ptp(pts, axis=0)), 1e-8)
    return center, 0.12 * span


def _draw_labeled_points(ax, points, color, edgecolor, labels=True,
                         prime=False, label_color=None, center=None,
                         offset_scale=None):
    if not points:
        return
    auto_center, auto_scale = _label_offset(points.values())
    # A degenerate (single-point) set has no meaningful own centroid/span;
    # fall back to the caller-supplied shared reference in that case only,
    # so a well-formed multi-point set keeps its own (already-correct) offset.
    if len(points) <= 1:
        if center is None:
            center = auto_center
        if offset_scale is None:
            offset_scale = auto_scale
    else:
        center = auto_center
        offset_scale = auto_scale
    for label, point in points.items():
        ax.scatter(point[0], point[1], s=85, c=color, edgecolors=edgecolor,
                   linewidths=0.7, zorder=5)
        if not labels:
            continue
        direction = point - center
        norm = np.linalg.norm(direction)
        offset = (direction / norm * offset_scale if norm > 1e-8
                  else np.array([offset_scale, 0.0]))
        ax.text(*(point + offset), _fig_label(label, prime=prime), fontsize=24,
                color=label_color if label_color is not None else edgecolor,
                ha="center", va="center", zorder=6)


# ---------------------------------------------------------------------------
# in-plane operation classification (Figure 3)
# ---------------------------------------------------------------------------

def _axis_name(axis):
    return ("kx", "ky", "kz")[axis]


def _operation_identity_on_plane(op, axis, tol=1e-8):
    """True if a full 3D operation leaves every in-plane k coordinate fixed."""
    op = np.array(op, dtype=float)
    for i in [a for a in range(3) if a != axis]:
        row = np.zeros(3)
        row[i] = 1.0
        if not np.allclose(op[i], row, atol=tol):
            return False
    return True


def _same_in_plane_operation(a, b, axis, tol=1e-6):
    """Compare only the in-plane 2x2 action of two fractional operations."""
    ip = [i for i in range(3) if i != axis]
    a2 = np.array(a, dtype=float)[np.ix_(ip, ip)]
    b2 = np.array(b, dtype=float)[np.ix_(ip, ip)]
    return np.allclose(a2, b2, atol=tol)


def _classify_spin_down_ops_2d(b_matrix, unique_ops, R, flip_ops_frac=None,
                               axis=2):
    """Mark each structural op spin-down if its in-plane action matches a flip op."""
    if not unique_ops:
        return np.zeros(0, dtype=bool)
    b_matrix = np.array(b_matrix, dtype=float)
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)
    flip_set = [np.array(f, dtype=float)
                for f in (flip_ops_frac if flip_ops_frac else [R])]
    mask = np.zeros(len(unique_ops), dtype=bool)
    for i, g_cart in enumerate(unique_ops):
        M = b_T_inv @ np.array(g_cart, dtype=float) @ b_T
        g_frac = np.linalg.inv(M.T)
        mask[i] = any(_same_in_plane_operation(g_frac, f, axis) for f in flip_set)
    return mask


def _dedupe_in_plane_ops_2d(b_matrix, unique_ops, spin_down_mask, axis=2,
                            decimals=8):
    """Keep one representative 3D op per distinct in-plane action."""
    b_matrix = np.array(b_matrix, dtype=float)
    b_T = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)
    ip = [i for i in range(3) if i != axis]
    kept_ops, kept_mask, seen, conflicts = [], [], {}, 0
    for i, g_cart in enumerate(unique_ops):
        M = b_T_inv @ np.array(g_cart, dtype=float) @ b_T
        g_frac = np.linalg.inv(M.T)
        block = g_frac[np.ix_(ip, ip)]
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
        print("[Warning] Some 3D operations collapse to the same 2D operation "
              "but carry different spin labels; those 2D sectors may be "
              "spin-degenerate.")
    return kept_ops, np.array(kept_mask, dtype=bool)


def _polygon_ray_exit(poly, unit_vec, origin=None):
    """Distance from ``origin`` to where a ray along ``unit_vec`` exits the
    closed polygon ``poly``, or None if it never crosses an edge."""
    origin = np.array([0.0, 0.0]) if origin is None else np.asarray(origin, dtype=float)
    closed = np.vstack([poly, poly[0]])
    t_hits = []
    for p1, p2 in zip(closed[:-1], closed[1:]):
        edge = p2 - p1
        mat = np.column_stack([unit_vec, -edge])
        det = np.linalg.det(mat)
        if abs(det) < 1e-12:
            continue
        t, u = np.linalg.solve(mat, p1 - origin)
        if t > 1e-10 and -1e-10 <= u <= 1.0 + 1e-10:
            t_hits.append(float(t))
    return min(t_hits) if t_hits else None


def _draw_reciprocal_axes_2d(ax, b_matrix, axis, basis, bz_poly):
    """Draw the two in-plane reciprocal axis arrows (b_i, i != vacuum axis).

    Same convention as ``draw_projected_reciprocal_axes`` in
    ``compute_centroid_hybrid.py`` (Figure 4): a dotted lead-in inside the BZ,
    a solid arrow outside, and a bold label at the tip. The vacuum-axis
    reciprocal vector is always perpendicular to this plane and carries no
    in-plane information, so it is skipped entirely (unlike the 3D top view,
    where all three b_i are shown for orientation against the full 3D BZ).
    """
    span = max(float(np.max(np.ptp(bz_poly, axis=0))), 1e-8)
    target = 0.60 * span
    origin = np.array([0.0, 0.0])
    color = "#202020"
    all_labels = [r"$\mathbf{b}_1$", r"$\mathbf{b}_2$", r"$\mathbf{b}_3$"]

    in_plane = [(vec, all_labels[i]) for i, vec in enumerate(np.array(b_matrix, dtype=float))
               if i != axis]
    for vec, label in in_plane:
        projected = _to_2d(vec, basis)
        length = float(np.linalg.norm(projected))
        if length < 1e-10:
            continue

        unit = projected / length
        exit_t = _polygon_ray_exit(bz_poly, unit)
        exit_len = exit_t if exit_t is not None else target * 0.65
        exit_pt = origin + unit * min(exit_len, target * 0.86)
        end = origin + unit * max(target, exit_len * 1.16)
        ax.plot([origin[0], exit_pt[0]], [origin[1], exit_pt[1]],
                color=color, ls=":", lw=1.5, alpha=0.6, zorder=219,
                clip_on=False)
        ann = ax.annotate(
            "",
            xy=end, xytext=exit_pt,
            arrowprops=dict(arrowstyle="->", color=color, lw=2.2,
                            mutation_scale=28, shrinkA=0, shrinkB=0),
            zorder=220,
            annotation_clip=False,
        )
        ann.set_clip_on(False)
        if ann.arrow_patch is not None:
            ann.arrow_patch.set_clip_on(False)
        offset = projected / length * (0.04 * span)
        ax.text(end[0] + offset[0], end[1] + offset[1], label,
                fontsize=24, fontweight="bold", ha="center", va="center",
                color=color, zorder=221, clip_on=False)


def _line_label_anchor(p_pos, p_neg, avoid_pts):
    """Pick whichever end of a through-Gamma line is farthest from the
    already-drawn high-symmetry points/labels, so the op label doesn't land
    on top of them. Same idea as the 3D ``_best_label_anchor``, simplified
    since there is no camera/projection in the 2D figures."""
    candidates = [p_pos, p_neg]
    if avoid_pts is None or len(avoid_pts) == 0:
        return candidates[0]
    avoid = np.asarray(avoid_pts, dtype=float)
    best, best_score = candidates[0], -np.inf
    for c in candidates:
        d = float(np.min(np.linalg.norm(avoid - c, axis=1)))
        if d > best_score:
            best_score, best = d, c
    return best


def _draw_op_visual_2d(ax, R_frac, b_matrix, basis, bz_poly, avoid_pts=None):
    """Draw the geometric visual for the selected 2D spin-flip operation.

    Reuses the same op classification/labeling as the 3D Figure 2
    (``_classify_spinflip_op``/``describe_spinflip_op`` in
    ``compute_centroid_hybrid.py``) and keeps the same 3-index (b1,b2,b3)
    reciprocal labeling convention, since the underlying operation is still
    the full 3x3 matrix -- only the drawing is 2D-specific:

      Cn/Sn rotation (det=+1, or det=-1 rotoreflection): the physical
        rotation axis is always exactly perpendicular to this plane (the
        vacuum direction), the same for every case, so there is nothing to
        project -- just a small curved arrow glyph centered on Gamma (the
        plot origin) with a Cn(+/-) label, instead of the 3D camera-aware
        ring + long axis arrow.
      Mirror (det=-1, order None): the mirror's fixed set is a LINE through
        Gamma (not a plane), drawn straight across the BZ, labeled by the
        plane normal's (b1,b2,b3) indices exactly as the 3D m_{hkl} label.

    Trivial/degenerate ops (identity, inversion, or an axis with no in-plane
    component) draw nothing, matching the 3D convention.
    """
    if _classify_spinflip_op is None:
        return
    b_matrix = np.array(b_matrix, dtype=float)
    R_frac = np.array(R_frac, dtype=float)
    b_T = b_matrix.T
    R_cart = b_T @ np.linalg.inv(R_frac).T @ np.linalg.inv(b_T)

    op = _classify_spinflip_op(R_cart)
    op_type = op['type']
    if op_type in ('identity', 'inversion') or op['axis'] is None:
        return

    span = max(float(np.max(np.ptp(bz_poly, axis=0))), 1e-8)
    origin = np.array([0.0, 0.0])
    COLOR = os.environ.get('ALTERSEEK_OP_AXIS_COLOR', '#00c853')

    if op_type in ('rotation', 'rotoreflection') and op['order'] and op['order'] >= 3:
        # Canonical "out of the page" direction for this screen basis, so the
        # drawn sense (+/- = CCW/CW) matches what's actually on screen.
        e1, e2, _ = basis if basis is not None else (
            np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), None)
        normal_dir = np.cross(e1, e2)
        axis = np.array(op['axis'], dtype=float)
        if np.dot(axis, normal_dir) < 0:
            axis = -axis
        sense = _rotation_sense(R_cart, axis)
        order = op['order']
        # International/Bilbao numeral convention (same as the 3D figure):
        # a proper rotation Cn is just its order n; an improper Sn is
        # relabeled by its rotoinversion order n-bar (S3<->6bar, S4<->4bar,
        # S6<->3bar). The axis subscript here is always the trivial vacuum
        # direction (this glyph only fires for the perpendicular-axis case),
        # kept for consistency with the 3D label format.
        if op_type == 'rotation':
            display_order = order
            improper = False
        else:
            display_order = {3: 6, 4: 4, 6: 3}.get(order, order)
            improper = True
        sgn = '' if sense == 0 else ('+' if sense > 0 else '-')
        idx = _reduce_int_vector(axis @ np.linalg.inv(b_matrix))
        axis_sub = "".join(rf"\bar{{{abs(i)}}}" if i < 0 else f"{i}" for i in idx)
        digit = rf"\bar{{{display_order}}}" if improper else f"{display_order}"
        label = (rf"$\mathbf{{{digit}^{{{sgn}}}_{{{axis_sub}}}}}$" if sgn
                 else rf"$\mathbf{{{digit}_{{{axis_sub}}}}}$")

        r_arc = 0.14 * span
        theta = np.radians(np.linspace(30.0, 330.0, 100))  # 300 deg, gap at bottom
        arc = origin + r_arc * np.column_stack([np.cos(theta), np.sin(theta)])
        # Solid shaft stops at the arrow base; only the final segment is the
        # arrowhead (same convention as the 3D rotation-axis arrow), so the
        # plain curved line never overlaps/co-terminates with the separately
        # drawn straight arrowhead triangle at the tip.
        head_cut = 4
        shaft = arc[:-head_cut]
        ax.plot(shaft[:, 0], shaft[:, 1], color=COLOR, lw=2.4, alpha=0.95, zorder=210)
        ann = ax.annotate(
            "", xy=arc[-1], xytext=arc[-head_cut],
            arrowprops=dict(arrowstyle="-|>", color=COLOR, lw=2.4,
                            mutation_scale=22, shrinkA=0, shrinkB=0),
            zorder=211, annotation_clip=False,
        )
        ann.set_clip_on(False)
        label_pt = origin + r_arc * 1.35 * np.array([0.0, -1.0])
        ax.text(label_pt[0], label_pt[1], label, fontsize=22, fontweight='bold',
                color=COLOR, ha='center', va='center', zorder=212, clip_on=False)

    elif op_type == 'rotation' and op['order'] == 2:
        # A C2 whose 3D axis lies IN this plane (common for 6/mmm-family
        # in-plane two-fold axes) restricts to a genuine in-plane mirror when
        # drawn flat: the axis direction itself is the fixed line (+1
        # eigenvalue), while the in-plane-perpendicular direction flips sign.
        # A C2 axis exactly along the vacuum direction (the excluded trivial
        # -I case) never reaches here since it's filtered out upstream.
        axis_full = np.array(op['axis'], dtype=float)
        axis_2d = _to_2d(axis_full, basis)
        axis_2d_len = float(np.linalg.norm(axis_2d))
        if axis_2d_len < 1e-8:
            return
        line_dir = axis_2d / axis_2d_len

        exit_pos = _polygon_ray_exit(bz_poly, line_dir, origin) or (0.5 * span)
        exit_neg = _polygon_ray_exit(bz_poly, -line_dir, origin) or (0.5 * span)
        p_pos = origin + line_dir * exit_pos
        p_neg = origin - line_dir * exit_neg
        ax.plot([p_neg[0], p_pos[0]], [p_neg[1], p_pos[1]],
                color=COLOR, lw=2.4, alpha=0.95, zorder=210, clip_on=False)

        idx = _reduce_int_vector(axis_full @ np.linalg.inv(b_matrix))
        label = _format_miller('2', idx)
        offset = line_dir * (0.05 * span)
        anchor = _line_label_anchor(p_pos, p_neg, avoid_pts)
        sign = 1.0 if np.array_equal(anchor, p_pos) else -1.0
        label_pt = anchor + sign * offset
        ax.text(label_pt[0], label_pt[1], label, fontsize=22, fontweight='bold',
                color=COLOR, ha='center', va='center', zorder=212, clip_on=False)

    elif op_type == 'mirror':
        normal3 = np.array(op['axis'], dtype=float)
        n2d = _to_2d(normal3, basis)
        n2d_len = float(np.linalg.norm(n2d))
        if n2d_len < 1e-8:
            return
        n2d_unit = n2d / n2d_len
        line_dir = np.array([-n2d_unit[1], n2d_unit[0]])

        exit_pos = _polygon_ray_exit(bz_poly, line_dir, origin) or (0.5 * span)
        exit_neg = _polygon_ray_exit(bz_poly, -line_dir, origin) or (0.5 * span)
        p_pos = origin + line_dir * exit_pos
        p_neg = origin - line_dir * exit_neg
        ax.plot([p_neg[0], p_pos[0]], [p_neg[1], p_pos[1]],
                color=COLOR, lw=2.4, alpha=0.95, zorder=210, clip_on=False)

        idx = _reduce_int_vector(normal3 @ np.linalg.inv(b_matrix))
        label = _format_miller('m', idx)
        offset = line_dir * (0.05 * span)
        anchor = _line_label_anchor(p_pos, p_neg, avoid_pts)
        sign = 1.0 if np.array_equal(anchor, p_pos) else -1.0
        label_pt = anchor + sign * offset
        ax.text(label_pt[0], label_pt[1], label, fontsize=22, fontweight='bold',
                color=COLOR, ha='center', va='center', zorder=212, clip_on=False)


def _plot_spin_pattern_top_view_2d(centroid_result, R_for_kpts,
                                   flip_ops_for_plot, output_path,
                                   show_title=False, show_legend=False):
    """Figure 3: color every 2D symmetry image of the IBZ red/blue by spin."""
    import matplotlib.pyplot as plt
    axis = int(centroid_result.get("vacuum_axis", 2))
    unique_ops = centroid_result.get("unique_ops")
    ibz_polygon_frac = centroid_result.get("ibz_polygon_frac")
    if not ibz_polygon_frac or len(ibz_polygon_frac) < 3 or not unique_ops:
        print("[Note] Skipping 2D spin-pattern figure (no 2D IBZ polygon or "
              "symmetry ops available).")
        return None

    b_matrix = np.array(centroid_result["b_matrix"], dtype=float)
    bz_poly, basis = _bz_polygon_2d(b_matrix, axis, cartesian_xy=(axis == 2))

    spin_down_mask = _classify_spin_down_ops_2d(
        b_matrix, unique_ops, R_for_kpts,
        flip_ops_frac=flip_ops_for_plot if flip_ops_for_plot else None,
        axis=axis)
    unique_ops, spin_down_mask = _dedupe_in_plane_ops_2d(
        b_matrix, unique_ops, spin_down_mask, axis=axis)

    identity_like = [op for op in (flip_ops_for_plot or [])
                     if _operation_identity_on_plane(op, axis)]
    if identity_like:
        print(f"[Warning] A spin-flip operation leaves the plane "
              f"{_axis_name(axis)}=0 unchanged; this protects spin degeneracy "
              "at the same in-plane k (not a nontrivial 2D splitting).")

    ibz_cart = np.array([_cart_from_frac(p, b_matrix) for p in ibz_polygon_frac],
                        dtype=float)

    fig, ax = plt.subplots(figsize=(9, 9))
    fill_alpha = 0.68
    up_labeled = down_labeled = False
    if identity_like or spin_down_mask.sum() == 0:
        ax.fill(bz_poly[:, 0], bz_poly[:, 1], facecolor="#b22222",
                alpha=fill_alpha, edgecolor="none", label="spin-up")
        up_labeled = True
    elif spin_down_mask.sum() == len(spin_down_mask):
        ax.fill(bz_poly[:, 0], bz_poly[:, 1], facecolor="#1f4e9e",
                alpha=fill_alpha, edgecolor="none", label="spin-down")
        down_labeled = True
    else:
        for i, g in enumerate(unique_ops):
            cell_pts = (np.array(g, dtype=float) @ ibz_cart.T).T
            poly = np.array([_to_2d(p, basis) for p in cell_pts], dtype=float)
            is_down = bool(spin_down_mask[i])
            color = "#1f4e9e" if is_down else "#b22222"
            if is_down:
                label = None if down_labeled else "spin-down"
                down_labeled = True
            else:
                label = None if up_labeled else "spin-up"
                up_labeled = True
            closed = np.vstack([poly, poly[0]])
            ax.fill(poly[:, 0], poly[:, 1], facecolor=color, alpha=fill_alpha,
                    edgecolor="none", label=label)
            ax.plot(closed[:, 0], closed[:, 1], color=color, lw=0.9, alpha=0.95)

    closed = np.vstack([bz_poly, bz_poly[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="black", lw=2.0, label="BZ boundary")
    _draw_reciprocal_axes_2d(ax, b_matrix, axis, basis, bz_poly)
    ax.set_aspect("equal", adjustable="box")
    if show_title:
        ax.set_title(r"Spin-up / Spin-down BZ ($k_{\rm vac}=0$)", fontsize=18)
    ax.set_axis_off()
    if show_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=12,
                 borderaxespad=0, frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------

def plot_2d_figures(centroid_result, general_kpoint, R_for_kpts, basename,
                    output_dir=".", flip_ops_for_plot=None):
    """Save the 2D Figure 1 (IBZ), 2 (spin-flip), and 3 (spin pattern) views.

    ``R_for_kpts`` is the selected spin-flip operation in the standardized
    primitive basis (the same matrix used to build the KPOINTS path), so the
    figure's k' and primed points coincide with the written path.
    ``flip_ops_for_plot`` (standardized primitive basis) drives the Figure 3
    spin-up/down coloring.  Returns the list of saved file paths.
    """
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"[Warning] Could not import matplotlib for 2D figures: {exc}")
        return []

    b_matrix = np.array(centroid_result["b_matrix"], dtype=float)
    axis = int(centroid_result.get("vacuum_axis", 2))
    sc_type = centroid_result.get("sc_type", "BZ")
    cartesian_xy = (axis == 2)
    bz_poly, basis = _bz_polygon_2d(b_matrix, axis, cartesian_xy=cartesian_xy)

    kpath = centroid_result["band_kpath"]
    coords = centroid_result["band_kpoints_frac"]
    # Only physical in-plane labels appear in the top-down view.
    kpoints_cart = {
        label: _to_2d(_cart_from_frac(frac, b_matrix), basis)
        for label, frac in coords.items()
        if not str(label).startswith("_") and abs(np.array(frac)[axis]) < 1e-4
    }

    centroid_xy = _to_2d(_cart_from_frac(general_kpoint, b_matrix), basis)
    R_inv_T = np.linalg.inv(np.array(R_for_kpts, dtype=float)).T
    k_prime = R_inv_T @ np.array(general_kpoint, dtype=float)
    k_prime[axis] = 0.0
    k_prime_xy = _to_2d(_cart_from_frac(k_prime, b_matrix), basis)

    saved = []

    # ----- Figure 1: 2D BZ + in-plane path + centroid -----
    fig1, ax1 = _setup_2d_ax(f"2D BZ: {basename} ({sc_type})", bz_poly)
    for start, end in kpath:
        if start in kpoints_cart and end in kpoints_cart:
            p1, p2 = kpoints_cart[start], kpoints_cart[end]
            ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], c="red", lw=3.0,
                     alpha=0.9, zorder=3)
    _draw_labeled_points(ax1, kpoints_cart, "red", "darkred", label_color="black")
    ax1.scatter(*centroid_xy, c="gold", marker="*", s=420, edgecolors="k",
                zorder=112, label=r"$k$")
    # Legend placed outside the axes (bbox_to_anchor) rather than in a
    # corner, so it never overlaps plotted points/labels and the BZ content
    # itself doesn't need extra padding to make room for it.
    ax1.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=18,
              borderaxespad=0, frameon=True)
    fig1.tight_layout()
    fig1_path = os.path.join(output_dir, f"{basename}_2d_ibz_{sc_type}.png")
    fig1.savefig(fig1_path, dpi=300, bbox_inches="tight")
    plt.close(fig1)
    saved.append(fig1_path)

    # ----- Figure 2: spin-up IBZ + spin-flip image -----
    mapped_cart_lines = {}
    mapped_cart = {}
    for label, frac in coords.items():
        if str(label).startswith("_") or abs(np.array(frac)[axis]) >= 1e-4:
            continue
        mapped_frac = R_inv_T @ np.array(frac, dtype=float)
        mapped_frac[axis] = 0.0
        mapped_point = _to_2d(_cart_from_frac(mapped_frac, b_matrix), basis)
        mapped_cart_lines[label + "'"] = mapped_point
        if not any(np.allclose(mapped_point, orig, atol=1e-8)
                   for orig in kpoints_cart.values()):
            mapped_cart[label + "'"] = mapped_point

    fig2, ax2 = _setup_2d_ax("2D spin-flip path connections", bz_poly)
    avoid_pts = (list(kpoints_cart.values()) + list(mapped_cart_lines.values())
                 + [centroid_xy, k_prime_xy])
    _draw_op_visual_2d(ax2, R_for_kpts, b_matrix, basis, bz_poly, avoid_pts=avoid_pts)
    orig_poly = np.array(list(kpoints_cart.values()), dtype=float)
    mapped_poly = np.array(list(mapped_cart_lines.values()), dtype=float)
    if len(orig_poly) >= 3 and ConvexHull is not None:
        hp = orig_poly[ConvexHull(orig_poly).vertices]
        ax2.fill(hp[:, 0], hp[:, 1], color="salmon", alpha=0.20, zorder=1)
    if len(mapped_poly) >= 3 and ConvexHull is not None:
        hp = mapped_poly[ConvexHull(mapped_poly).vertices]
        ax2.fill(hp[:, 0], hp[:, 1], color="cornflowerblue", alpha=0.20, zorder=1)
    for start, end in kpath:
        if start in kpoints_cart and end in kpoints_cart:
            p1, p2 = kpoints_cart[start], kpoints_cart[end]
            ax2.plot([p1[0], p2[0]], [p1[1], p2[1]], c="red", lw=4.0,
                     alpha=0.9, zorder=50)
    for start, end in kpath:
        sp, ep = start + "'", end + "'"
        if sp in mapped_cart_lines and ep in mapped_cart_lines:
            p1, p2 = mapped_cart_lines[sp], mapped_cart_lines[ep]
            ax2.plot([p1[0], p2[0]], [p1[1], p2[1]], c="navy", lw=4.0,
                     alpha=0.9, zorder=50)
    shared_points = list(kpoints_cart.values()) + list(mapped_cart_lines.values())
    shared_center, shared_scale = _label_offset(shared_points)
    _draw_labeled_points(ax2, kpoints_cart, "salmon", "darkred", prime=False,
                         center=shared_center, offset_scale=shared_scale)
    _draw_labeled_points(ax2, mapped_cart, "cornflowerblue", "navy", prime=True,
                         center=shared_center, offset_scale=shared_scale)
    threshold = 0.05 * max(np.max(np.linalg.norm(bz_poly, axis=1)), 1e-8)
    for point in kpoints_cart.values():
        if np.linalg.norm(point - centroid_xy) > threshold:
            ax2.plot([point[0], centroid_xy[0]], [point[1], centroid_xy[1]],
                     c="deepskyblue", lw=2.0, ls="--", alpha=0.75, zorder=40)
    for point in mapped_cart_lines.values():
        if np.linalg.norm(point - k_prime_xy) > threshold:
            ax2.plot([point[0], k_prime_xy[0]], [point[1], k_prime_xy[1]],
                     c="deepskyblue", lw=2.0, ls="--", alpha=0.75, zorder=40)
    ax2.scatter(*centroid_xy, c="gold", s=300, marker="*", edgecolors="k",
                linewidths=0.8, zorder=120, label=r"$k$")
    ax2.scatter(*k_prime_xy, c="cornflowerblue", s=150, marker="o",
                edgecolors="k", linewidths=0.8, zorder=120, label=r"$k'$")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=18,
              borderaxespad=0, frameon=True)
    fig2.tight_layout()
    fig2_path = os.path.join(output_dir, f"{basename}_2d_spinflip_{sc_type}.png")
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)
    saved.append(fig2_path)

    # ----- Figure 3: spin-up / spin-down BZ pattern -----
    fig3_path = os.path.join(output_dir, f"{basename}_2d_spinbz_{sc_type}.png")
    fig3_saved = _plot_spin_pattern_top_view_2d(
        centroid_result, R_for_kpts, flip_ops_for_plot, fig3_path)
    if fig3_saved is not None:
        saved.append(fig3_saved)

    print(f"Saved 2D figures: {', '.join(saved)}")
    return saved
