"""Top-down 2D Brillouin-zone figures for AlterSeeK-Path slab mode.

For a 2D material computed as a slab the physical reciprocal space is the
``k[vacuum_axis] = 0`` plane.  These figures render that plane straight-on as a
flat top-down view, instead of the tilted 3D BZ plate used for bulk crystals:

  * Figure 1 (``*_2d_ibz_*``)      -- 2D BZ outline + in-plane band path + k.
  * Figure 2 (``*_2d_spinflip_*``) -- spin-up IBZ (red) and its spin-flip image
                                      (blue) with the general point k and its
                                      spin-flip partner k'.
  * Figure 3 (``*_2d_spinbz_*``)   -- every symmetry image of the IBZ tiling
                                      the BZ, colored red/blue by spin.

The vacuum axis is taken from ``centroid_result['vacuum_axis']``; the screen
basis for the remaining two in-plane directions is chosen by
``_plane_projector``.
"""

import os
import numpy as np

from ..plotting_common import (
    IBZ_FACE_COLORS, _figure_output_paths, _math_label, _print_saved_paths,
    _save_figure, combine_point_labels, generated_plain_path_segments,
    grouped_point_labels, label_aliases, prime_point_label,
)

from scipy.spatial import ConvexHull

from ..symmetry import (
    _classify_spinflip_op, _doubled_ibz_extra_flags, _reduce_int_vector,
    _format_miller, _rotation_sense,
)


# --- Geometry helpers ---

def _cart_from_frac(frac, b_matrix):
    f = np.array(frac[:3], dtype=float)
    return f[0] * b_matrix[0] + f[1] * b_matrix[1] + f[2] * b_matrix[2]


def _plane_projector(b_matrix, axis, lattice_class=None):
    """2D screen basis (kx, ky) for the in-plane reciprocal vectors.

    Screen x/y are the structure's own Cartesian axes when the vacuum
    direction is axis-aligned, and follow the first in-plane reciprocal
    vector otherwise.  A centred rectangular lattice instead uses its
    conventional axes A = a1 + a2 and B = a1 - a2, so the figure does not
    depend on how the input file happened to be oriented.
    """
    in_plane_axes = [i for i in range(3) if i != axis]
    b_matrix = np.array(b_matrix, dtype=float)
    g1 = b_matrix[in_plane_axes[0]]
    g2 = b_matrix[in_plane_axes[1]]
    normal = np.cross(g1, g2)
    normal_norm = np.linalg.norm(normal)
    if normal_norm < 1e-14:
        raise ValueError("collinear in-plane reciprocal vectors")
    normal /= normal_norm

    if lattice_class == "centered-rectangular":
        # The conventional axis A = a1 + a2 is parallel to g1 + g2, so screen x
        # is the conventional a direction in both metric branches.
        conv_a = g1 + g2
        e1 = conv_a / np.linalg.norm(conv_a)
        # Derive screen y from the plane normal so the frame stays right-handed
        # and the drawn sense of every rotation operation is preserved.
        e2 = np.cross(normal, e1)
        e2 /= np.linalg.norm(e2)
        return (e1, e2, in_plane_axes)

    axis_aligned = next(
        (i for i in range(3) if abs(abs(normal[i]) - 1.0) < 1e-6), None
    )
    if axis_aligned is not None:
        e1 = np.zeros(3)
        e2 = np.zeros(3)
        fixed_axes = [i for i in range(3) if i != axis_aligned]
        e1[fixed_axes[0]] = 1.0
        e2[fixed_axes[1]] = 1.0
        if np.dot(normal, np.cross(e1, e2)) < 0.0:
            e2 = -e2
        return (e1, e2, in_plane_axes)

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


def _bz_polygon_2d(b_matrix, axis, radius=2, cartesian_xy=False,
                   lattice_class=None):
    """2D Wigner-Seitz BZ polygon for the selected reciprocal plane."""
    if cartesian_xy:
        if axis != 2:
            raise ValueError("Cartesian kx/ky 2D plotting requires the kz=0 plane")
        basis = None
        in_plane_axes = [0, 1]
    else:
        basis = _plane_projector(b_matrix, axis, lattice_class=lattice_class)
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


# --- Drawing helpers ---

# Zoom: the drawing half-width is 1.22x the Gamma-to-furthest-vertex distance,
# leaving a 22% margin around the BZ.  Measuring from Gamma rather than the
# polygon's bounding box keeps the zoom unchanged when the screen frame rotates.
_AX_LIMIT_FACTOR = 1.22


def _ax_limit(bz_poly):
    """Half-width of the drawing area, shared by the axes and the b_i arrows."""
    radius = max(float(np.max(np.linalg.norm(bz_poly, axis=1))), 1e-8)
    return _AX_LIMIT_FACTOR * radius


def _setup_2d_ax(title, bz_poly):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))
    closed = np.vstack([bz_poly, bz_poly[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="0.25", lw=2.2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    # A b_i arrow pointing straight up nearly reaches the top of the axes, so
    # the title needs padding to clear it.
    ax.set_title(title, fontsize=20, pad=18)
    limit = _ax_limit(bz_poly)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    return fig, ax


# How far a label sits from its marker, in points.  The x gap is larger because
# a text box is much wider than tall, so equal gaps would look uneven.
_LABEL_OFFSET_POINTS_X = 14.0
_LABEL_OFFSET_POINTS_Y = 9.0
# Extra gap for a label whose point also carries a reciprocal-axis arrow.
_AXIS_LABEL_GAP_SCALE = 1.7


def _bz_outward_normal(point, bz_poly):
    """Outward normal of the BZ edge(s) a point sits on.

    Boundary labels read best pushed straight out of the zone; the radial
    direction from the wedge centre aims diagonally and drops the text on the
    boundary line.  At a corner the two edge normals are averaged.  Returns
    None for a point that is not on the boundary.
    """
    if bz_poly is None or len(bz_poly) < 3:
        return None
    poly = np.asarray(bz_poly, dtype=float)
    interior = poly.mean(axis=0)
    tol = 1e-3 * max(np.max(np.ptp(poly, axis=0)), 1e-12)
    normals = []
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        edge = b - a
        length = np.linalg.norm(edge)
        if length < 1e-12:
            continue
        t = np.clip(np.dot(point - a, edge) / length ** 2, 0.0, 1.0)
        if np.linalg.norm(point - (a + t * edge)) > tol:
            continue
        normal = np.array([edge[1], -edge[0]], dtype=float) / length
        if np.dot(normal, 0.5 * (a + b) - interior) < 0:
            normal = -normal
        if not any(np.allclose(normal, kept, atol=1e-6) for kept in normals):
            normals.append(normal)
    if not normals:
        return None
    bisector = np.sum(normals, axis=0)
    length = np.linalg.norm(bisector)
    if length < 1e-12:
        return normals[0]
    return bisector / length


def _label_offset(points):
    pts = np.array(list(points), dtype=float)
    center = np.mean(pts, axis=0)
    span = max(np.max(np.ptp(pts, axis=0)), 1e-8)
    return center, 0.12 * span


def _steer_off_line(direction, line_dir, point, other_points, offset_scale):
    """Turn a label off a line drawn through it, e.g. the mirror/axis visual.

    A label pushed along the operation line lands on top of it.  Both
    perpendiculars clear the line equally, so pick the one that ends up
    farther from the other labelled points.
    """
    line_dir = np.asarray(line_dir, dtype=float)[:2]
    norm = np.linalg.norm(line_dir)
    if norm < 1e-12:
        return direction
    line_dir = line_dir / norm
    if abs(float(np.dot(direction, line_dir))) < np.cos(np.radians(35.0)):
        return direction
    candidates = [np.array([-line_dir[1], line_dir[0]]),
                  np.array([line_dir[1], -line_dir[0]])]
    others = [np.asarray(other, dtype=float)[:2] for other in other_points
              if not np.allclose(np.asarray(other, dtype=float)[:2], point[:2])]
    if not others:
        return candidates[0]

    def clearance(candidate):
        probe = point[:2] + candidate * max(offset_scale, 1e-12)
        return min(float(np.linalg.norm(probe - other)) for other in others)

    return max(candidates, key=clearance)


def _draw_mixed_spin_label(ax, name, point, offset, alignment,
                           unprimed_color, primed_color, zorder):
    """Draw a slash label whose primed aliases use the spin-down color."""
    aliases = label_aliases(name)
    primed = [alias.endswith("'") for alias in aliases]
    if len(aliases) < 2 or not any(primed) or all(primed):
        return None

    from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea

    children = []
    alias_colors = []
    for index, alias in enumerate(aliases):
        if index:
            children.append(TextArea(
                "/", textprops={"fontsize": 24, "color": unprimed_color},
            ))
        color = primed_color if alias.endswith("'") else unprimed_color
        children.append(TextArea(
            _math_label(alias),
            textprops={"fontsize": 24, "color": color},
        ))
        alias_colors.append((alias, color))

    packed = HPacker(children=children, align="center", pad=0, sep=0)
    annotation = AnnotationBbox(
        packed,
        (point[0], point[1]),
        xybox=offset,
        xycoords="data",
        boxcoords="offset points",
        box_alignment=alignment,
        frameon=False,
        annotation_clip=False,
        zorder=zorder,
    )
    annotation.set_clip_on(False)
    annotation._alterseek_label = name
    annotation._alterseek_alias_colors = tuple(alias_colors)
    ax.add_artist(annotation)
    return annotation


def _draw_labeled_points(ax, points, color, edgecolor, labels=True,
                          prime=False, label_color=None, center=None,
                          offset_scale=None, path_labels=(), bz_poly=None,
                          avoid_dir=None, avoid_dirs=(), avoid_points=None,
                          prime_label_color=None, zorder=5):
    if not points:
        return
    auto_center, auto_scale = _label_offset(points.values())
    # A single-point set has no meaningful centroid or span, so only that case
    # falls back to the caller's shared reference; a normal multi-point set
    # keeps its own offset.
    if len(points) <= 1:
        if center is None:
            center = auto_center
        if offset_scale is None:
            offset_scale = auto_scale
    else:
        center = auto_center
        offset_scale = auto_scale
    layout_points = points.values() if avoid_points is None else avoid_points
    for label, point in grouped_point_labels(points, path_labels):
        ax.scatter(point[0], point[1], s=85, c=color, edgecolors=edgecolor,
                   linewidths=0.7, zorder=zorder)
        if not labels:
            continue
        direction = point - center
        norm = np.linalg.norm(direction)
        direction = (direction / norm if norm > 1e-8
                     else np.array([1.0, 0.0]))
        # A label that has to share its point with a b_i arrow needs a wider
        # berth than the standard gap, since the arrow shaft runs through the
        # marker itself.
        gap_scale = 1.0
        boundary_normal = _bz_outward_normal(point, bz_poly)
        if boundary_normal is not None:
            direction = boundary_normal
            # If a reciprocal axis exits the BZ through this point, retain the
            # outward component and add a tangent component away from nearby
            # labels.  The axis remains visible without putting the label back
            # onto the BZ boundary.
            for line_dir in avoid_dirs:
                line_dir = np.asarray(line_dir, dtype=float)[:2]
                line_norm = np.linalg.norm(line_dir)
                if line_norm < 1e-12:
                    continue
                line_dir = line_dir / line_norm
                scale = max(float(np.linalg.norm(point)), 1.0)
                cross = line_dir[0] * point[1] - line_dir[1] * point[0]
                if abs(float(cross)) > 1e-8 * scale:
                    continue
                perpendiculars = [np.array([-line_dir[1], line_dir[0]]),
                                  np.array([line_dir[1], -line_dir[0]])]
                if abs(line_dir[0]) >= 2.0 * abs(line_dir[1]):
                    tangent = min(perpendiculars, key=lambda side: side[1])
                elif abs(line_dir[1]) >= 2.0 * abs(line_dir[0]):
                    tangent = min(perpendiculars, key=lambda side: side[0])
                else:
                    tangent = _steer_off_line(
                        line_dir, line_dir, point, layout_points, offset_scale,
                    )
                # Weighted towards the tangent rather than the outward
                # normal: the outward direction here is exactly where the b_i
                # arrow and its own bold label sit, so a mostly-radial label
                # lands on the arrow shaft.
                direction = direction + 1.4 * tangent
                direction = direction / np.linalg.norm(direction)
                gap_scale = _AXIS_LABEL_GAP_SCALE
                break
        elif avoid_dir is not None:
            line_dir = np.asarray(avoid_dir, dtype=float)[:2]
            line_norm = np.linalg.norm(line_dir)
            crossing_dir = next((axis_dir for axis_dir in avoid_dirs
                                 if line_norm >= 1e-12
                                 and np.linalg.norm(axis_dir) >= 1e-12
                                 and abs(float(np.dot(line_dir / line_norm,
                                     axis_dir / np.linalg.norm(axis_dir)))) < 0.2), None)
            if np.linalg.norm(point[:2]) < 1e-8 and crossing_dir is not None:
                side_a = _steer_off_line(line_dir, line_dir, point,
                                         layout_points, offset_scale)
                side_b = _steer_off_line(crossing_dir, crossing_dir, point,
                                         layout_points, offset_scale)
                candidates = [a * side_a + b * side_b
                              for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
                candidates = [candidate / np.linalg.norm(candidate)
                              for candidate in candidates]
                all_dirs = [line_dir, *avoid_dirs]
                direction = max(candidates, key=lambda candidate: (
                    min(1.0 - abs(float(np.dot(candidate, axis_dir / np.linalg.norm(axis_dir))))
                        for axis_dir in all_dirs if np.linalg.norm(axis_dir) >= 1e-12),
                    min(np.linalg.norm(point + candidate * offset_scale - other)
                        for other in layout_points if not np.allclose(other, point)),
                ))
            else:
                direction = _steer_off_line(direction, avoid_dir, point,
                                            points.values(), offset_scale)
        # Share the 3D figures' label typography instead of a second formatter.
        name = str(label)
        if prime and not name.endswith("'"):
            name += "'"
        # Offset in typographic units and align the near edge of the text box,
        # so the gap does not shrink as the label gets wider or the wedge smaller.
        dx, dy = float(direction[0]), float(direction[1])
        text_color = label_color if label_color is not None else edgecolor
        offset = (_LABEL_OFFSET_POINTS_X * gap_scale * dx,
                  _LABEL_OFFSET_POINTS_Y * gap_scale * dy)
        alignment = (
            0.0 if dx > 0.3 else 1.0 if dx < -0.3 else 0.5,
            0.0 if dy > 0.3 else 1.0 if dy < -0.3 else 0.5,
        )
        if prime_label_color is not None and _draw_mixed_spin_label(
            ax,
            name,
            point,
            offset,
            alignment,
            text_color,
            prime_label_color,
            zorder + 1,
        ) is not None:
            continue
        ax.annotate(
            _math_label(name), xy=(point[0], point[1]),
            textcoords="offset points",
            xytext=offset,
            fontsize=24,
            color=text_color,
            ha="left" if dx > 0.3 else "right" if dx < -0.3 else "center",
            va="bottom" if dy > 0.3 else "top" if dy < -0.3 else "center",
            zorder=zorder + 1, annotation_clip=False)


# --- In-plane operation classification ---

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


# --- Reciprocal-axis arrows and figure output ---

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


def _draw_reciprocal_axes_2d(ax, b_matrix, axis, basis, bz_poly,
                             zorder=219):
    """Draw the two in-plane reciprocal axis arrows (b_i, i != vacuum axis).

    Same convention as ``draw_projected_reciprocal_axes`` in ``plotting_3d.py``.
    The vacuum-axis b_i is perpendicular to the plane and carries no in-plane
    information, so it is skipped.
    """
    span = max(float(np.max(np.ptp(bz_poly, axis=0))), 1e-8)
    target = 0.60 * span
    # A b_i pointing straight up would run its tip and label into the title, so
    # cap each arrow against the axes edge it actually approaches rather than by
    # one shared length -- a near-horizontal b_i keeps its full reach.
    label_gap = 0.04 * span
    tip_budget = max(_ax_limit(bz_poly) * 0.94 - label_gap, 0.1 * span)
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
        max_tip = tip_budget / max(float(np.max(np.abs(unit))), 1e-8)
        tip_len = min(max(target, exit_len * 1.16), max_tip)
        exit_pt = origin + unit * min(exit_len, tip_len * 0.86)
        end = origin + unit * tip_len
        ax.plot([origin[0], exit_pt[0]], [origin[1], exit_pt[1]],
                color=color, ls=":", lw=1.5, alpha=0.6, zorder=zorder,
                clip_on=False)
        ann = ax.annotate(
            "",
            xy=end, xytext=exit_pt,
            arrowprops=dict(arrowstyle="->", color=color, lw=2.2,
                            mutation_scale=28, shrinkA=0, shrinkB=0),
            zorder=zorder + 1,
            annotation_clip=False,
        )
        ann.set_clip_on(False)
        if ann.arrow_patch is not None:
            ann.arrow_patch.set_clip_on(False)
        offset = unit * label_gap
        ax.text(end[0] + offset[0], end[1] + offset[1], label,
                fontsize=24, fontweight="bold", ha="center", va="center",
                color=color, zorder=zorder + 2, clip_on=False)


def _finish_2d_figure(fig, output_path, save_pdf, deferred_figures=None):
    """Save now for direct calls or defer display/save to the workflow."""
    import matplotlib.pyplot as plt
    extra_formats = ("pdf",) if save_pdf else ()
    if deferred_figures is None:
        saved_paths = _save_figure(
            fig,
            output_path,
            extra_formats=extra_formats,
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)
        return saved_paths

    expected_paths = _figure_output_paths(
        output_path, extra_formats=extra_formats
    )

    def _save_after_show(fig=fig):
        try:
            saved_paths = _save_figure(
                fig,
                output_path,
                extra_formats=extra_formats,
                dpi=300,
                bbox_inches="tight",
            )
            _print_saved_paths(saved_paths)
        finally:
            plt.close(fig)

    fig._alterseek_save_after_show = _save_after_show
    deferred_figures.append(fig)
    return expected_paths


# --- Operation visuals ---

def _line_operation_label_position(p_pos, p_neg, avoid_pts):
    """Return the line endpoint farthest from the avoid points."""
    positions = [p_pos, p_neg]
    if avoid_pts is None or len(avoid_pts) == 0:
        return positions[0]
    avoid = np.asarray(avoid_pts, dtype=float)
    best, best_score = positions[0], -np.inf
    for position in positions:
        d = float(np.min(np.linalg.norm(avoid - position, axis=1)))
        if d > best_score:
            best_score, best = d, position
    return best


def _operation_line_endpoints(bz_poly, line_dir, span, origin=None):
    """Return a through-Gamma line that overhangs the BZ on both sides."""
    origin = (np.array([0.0, 0.0]) if origin is None
              else np.asarray(origin, dtype=float))
    line_dir = np.asarray(line_dir, dtype=float)
    extension = 0.04 * span
    exit_pos = _polygon_ray_exit(bz_poly, line_dir, origin) or (0.5 * span)
    exit_neg = _polygon_ray_exit(bz_poly, -line_dir, origin) or (0.5 * span)
    p_pos = origin + line_dir * (exit_pos + extension)
    p_neg = origin - line_dir * (exit_neg + extension)
    return p_pos, p_neg


def _draw_op_visual_2d(ax, R_frac, b_matrix, basis, bz_poly, avoid_pts=None):
    """Draw the geometric visual for the selected 2D spin-flip operation.

    An order >= 3 rotation or rotoreflection becomes a curved arrow at Gamma
    labeled with its International symbol and sense.  An in-plane C2 or a mirror
    becomes a line across the BZ, labeled by the (b1,b2,b3) indices of its axis
    or plane normal.  Identity, inversion, and axes with no in-plane component
    draw nothing.  Classification is shared with the 3D figures
    (``_classify_spinflip_op`` in ``symmetry.py``).

    Returns the unit direction of the line drawn through Gamma, so point labels
    can be kept off it, or None when the visual is an arc or nothing at all.
    """
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
        # Use the screen basis's canonical out-of-page direction so the
        # displayed +/- sense matches the visible counterclockwise or clockwise
        # rotation.
        e1, e2, _ = basis if basis is not None else (
            np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), None)
        normal_dir = np.cross(e1, e2)
        axis = np.array(op['axis'], dtype=float)
        if np.dot(axis, normal_dir) < 0:
            axis = -axis
        sense = _rotation_sense(R_cart, axis)
        order = op['order']
        # Use International/Bilbao labels: a proper Cn is shown as n, while S3,
        # S4, and S6 become 6-bar, 4-bar, and 3-bar.
        # Retain the trivial vacuum-axis subscript because this glyph is used
        # only for the perpendicular-axis case and should match the 3D label
        # format.
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
        # Match the arc direction to the displayed rotation sense so n+ is
        # counterclockwise and n- is clockwise.
        dir_sign = -1 if sense < 0 else 1
        theta = np.radians(np.linspace(180.0 - dir_sign * 150.0,
                                        180.0 + dir_sign * 150.0, 100))  # 300 deg, gap at bottom
        arc = origin + r_arc * np.column_stack([np.cos(theta), np.sin(theta)])
        # Stop the solid shaft at the arrow base and draw only the final segment
        # as the arrowhead so the two strokes do not overlap.
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
        # A C2 whose 3D axis lies in the plane restricts to a genuine mirror in
        # the flat view: the axis is fixed while the in-plane perpendicular
        # direction changes sign.
        # A C2 along the vacuum direction would reduce to the excluded trivial
        # -I action and is filtered upstream.
        axis_full = np.array(op['axis'], dtype=float)
        axis_2d = _to_2d(axis_full, basis)
        axis_2d_len = float(np.linalg.norm(axis_2d))
        if axis_2d_len < 1e-8:
            return
        line_dir = axis_2d / axis_2d_len

        p_pos, p_neg = _operation_line_endpoints(
            bz_poly, line_dir, span, origin,
        )
        ax.plot([p_neg[0], p_pos[0]], [p_neg[1], p_pos[1]],
                color=COLOR, lw=2.4, alpha=0.95, zorder=45, clip_on=False)

        idx = _reduce_int_vector(axis_full @ np.linalg.inv(b_matrix))
        label = _format_miller('2', idx)
        label_position = _line_operation_label_position(p_pos, p_neg, avoid_pts)
        sign = 1.0 if np.array_equal(label_position, p_pos) else -1.0
        # The label position sits on the BZ boundary, so offset the text box's near edge
        # outward rather than centring it there, which straddles the boundary.
        dx, dy = float(sign * line_dir[0]), float(sign * line_dir[1])
        ax.annotate(
            label, xy=(label_position[0], label_position[1]), textcoords="offset points",
            xytext=(_LABEL_OFFSET_POINTS_X * dx,
                    _LABEL_OFFSET_POINTS_Y * dy),
            fontsize=22, fontweight='bold', color=COLOR,
            ha="left" if dx > 0.3 else "right" if dx < -0.3 else "center",
            va="bottom" if dy > 0.3 else "top" if dy < -0.3 else "center",
            zorder=212, annotation_clip=False)
        return line_dir

    elif op_type == 'mirror':
        normal3 = np.array(op['axis'], dtype=float)
        n2d = _to_2d(normal3, basis)
        n2d_len = float(np.linalg.norm(n2d))
        if n2d_len < 1e-8:
            return
        n2d_unit = n2d / n2d_len
        line_dir = np.array([-n2d_unit[1], n2d_unit[0]])

        p_pos, p_neg = _operation_line_endpoints(
            bz_poly, line_dir, span, origin,
        )
        ax.plot([p_neg[0], p_pos[0]], [p_neg[1], p_pos[1]],
                color=COLOR, lw=2.4, alpha=0.95, zorder=45, clip_on=False)

        idx = _reduce_int_vector(normal3 @ np.linalg.inv(b_matrix))
        label = _format_miller('m', idx)
        label_position = _line_operation_label_position(p_pos, p_neg, avoid_pts)
        sign = 1.0 if np.array_equal(label_position, p_pos) else -1.0
        # The label position sits on the BZ boundary, so offset the text box's near edge
        # outward rather than centring it there, which straddles the boundary.
        dx, dy = float(sign * line_dir[0]), float(sign * line_dir[1])
        ax.annotate(
            label, xy=(label_position[0], label_position[1]), textcoords="offset points",
            xytext=(_LABEL_OFFSET_POINTS_X * dx,
                    _LABEL_OFFSET_POINTS_Y * dy),
            fontsize=22, fontweight='bold', color=COLOR,
            ha="left" if dx > 0.3 else "right" if dx < -0.3 else "center",
            va="bottom" if dy > 0.3 else "top" if dy < -0.3 else "center",
            zorder=212, annotation_clip=False)
        return line_dir


# --- Spin-flip labels and the spin-pattern figure ---

def _order_mixed_label_left_to_right(label, point, spin_up_center,
                                     spin_down_center, bz_span):
    """Match mixed-label order to a clear left/right spin-sector layout."""
    aliases = label_aliases(label)
    primed = [alias.endswith("'") for alias in aliases]
    if len(aliases) < 2 or not any(primed) or all(primed):
        return label
    if spin_up_center is None or spin_down_center is None:
        return label

    point = np.asarray(point, dtype=float)
    up_delta = np.asarray(spin_up_center, dtype=float) - point
    down_delta = np.asarray(spin_down_center, dtype=float) - point
    min_horizontal = 0.05 * max(float(bz_span), 1e-8)
    if (abs(up_delta[0]) < min_horizontal
            or abs(down_delta[0]) < min_horizontal
            or up_delta[0] * down_delta[0] >= 0.0):
        return label

    center_separation = down_delta - up_delta
    if abs(center_separation[0]) <= 1.25 * abs(center_separation[1]):
        return label

    ordered = sorted(
        enumerate(aliases),
        key=lambda item: (
            down_delta[0] if item[1].endswith("'") else up_delta[0],
            item[0],
        ),
    )
    return combine_point_labels(*(alias for _index, alias in ordered))


def _spinflip_display_points_2d(original_points, mapped_points, path_labels,
                                spin_up_center=None, spin_down_center=None,
                                bz_span=1.0):
    """Return 2D spin-up/down point maps with coincident labels combined.

    Coincident original and mapped vertices are drawn once with the spin-up
    marker style, while every name used by the generated path is retained in
    one slash-combined label, such as ``Y/Y'``. The 3D figure suppresses
    coincident mapped labels.
    """
    combined_points = dict(original_points)
    for label, point in mapped_points.items():
        if label not in combined_points:
            combined_points[label] = point

    grouped = grouped_point_labels(combined_points, path_labels)
    original_coords = list(original_points.values())
    spin_up, spin_down = {}, {}
    for label, point in grouped:
        label = _order_mixed_label_left_to_right(
            label,
            point,
            spin_up_center,
            spin_down_center,
            bz_span,
        )
        target = spin_up if any(
            np.allclose(point, original, atol=1e-8, rtol=0.0)
            for original in original_coords
        ) else spin_down
        target[label] = point
    return spin_up, spin_down


def _plot_spin_pattern_top_view_2d(centroid_result, R_for_kpts,
                                   flip_ops_for_plot, output_path,
                                   show_title=False, show_legend=False,
                                   save_pdf=False, deferred_figures=None):
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
    # A standardized fractional c axis is not necessarily Cartesian z.  For
    # example, oA2 can place its physical vacuum reciprocal vector along kx.
    # Always construct the screen basis from the physical reciprocal plane.
    bz_poly, basis = _bz_polygon_2d(
        b_matrix, axis, lattice_class=centroid_result.get("sc_type"))

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

    # A project-doubled IBZ (the square lattice's X_A at point group 4)
    # carries an extra sector beyond the conventional wedge; one solid fill
    # would hide it.  Fill the whole polygon in a lighter color and overlay
    # the genuine (non-"_A") sub-polygon at full strength, mirroring
    # plotting_3d.plot_spin_bz_top_view_figure.  ibz_polygon_frac's vertices
    # are already in order around the boundary, so a subset can be filled
    # directly, with no convex-hull step.
    ibz_polygon_labels = centroid_result.get("ibz_polygon_labels")
    extra_flags = None
    if ibz_polygon_labels and len(ibz_polygon_labels) == len(ibz_polygon_frac):
        flags = _doubled_ibz_extra_flags(ibz_polygon_labels)
        if any(flags):
            extra_flags = flags
    original_idx = (
        [j for j, is_extra in enumerate(extra_flags) if not is_extra]
        if extra_flags is not None else None
    )

    # Same canvas and same limits as Figures 1 and 2 (_setup_2d_ax), so the BZ
    # is drawn at an identical scale across all three views of a case.  Left to
    # autoscale, this figure hugged the BZ and rendered it much larger than the
    # other two.
    fig, ax = plt.subplots(figsize=(8, 8))
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
            extra_color = (IBZ_FACE_COLORS["down_extra"] if is_down
                           else IBZ_FACE_COLORS["up_extra"])
            if is_down:
                label = None if down_labeled else "spin-down"
                down_labeled = True
            else:
                label = None if up_labeled else "spin-up"
                up_labeled = True
            closed = np.vstack([poly, poly[0]])
            if original_idx is not None and len(original_idx) >= 3:
                ax.fill(poly[:, 0], poly[:, 1], facecolor=extra_color,
                        alpha=0.46, edgecolor="none", label=label)
                sub_poly = poly[original_idx]
                ax.fill(sub_poly[:, 0], sub_poly[:, 1], facecolor=color,
                        alpha=fill_alpha, edgecolor="none")
            else:
                ax.fill(poly[:, 0], poly[:, 1], facecolor=color,
                        alpha=fill_alpha, edgecolor="none", label=label)
            ax.plot(closed[:, 0], closed[:, 1], color=color, lw=0.9, alpha=0.95)

    closed = np.vstack([bz_poly, bz_poly[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="black", lw=2.0, label="BZ boundary")
    _draw_reciprocal_axes_2d(ax, b_matrix, axis, basis, bz_poly)
    ax.set_aspect("equal", adjustable="box")
    limit = _ax_limit(bz_poly)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    if show_title:
        ax.set_title(r"Spin-up / Spin-down BZ ($k_{\rm vac}=0$)", fontsize=18)
    ax.set_axis_off()
    if show_legend:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=12,
                 borderaxespad=0, frameon=True)
    fig.tight_layout()
    return _finish_2d_figure(
        fig, output_path, save_pdf, deferred_figures=deferred_figures
    )


# --- Main plotting entry point ---

def _physical_lattice_title(centroid_result):
    """Return a physical 2D lattice name without an HPKOT setting tag."""
    lattice_class = centroid_result.get("lattice_class_2d")
    if lattice_class:
        return str(lattice_class).replace("_", " ")
    sc_type = str(centroid_result.get("sc_type", "BZ"))
    if sc_type.startswith(("oA", "oC", "mC")):
        return "centered rectangular"
    return sc_type.replace("_", " ").replace("-", " ")


def plot_2d_figures(centroid_result, general_kpoint, R_for_kpts, basename,
                    output_dir=".", flip_ops_for_plot=None, save_pdf=False,
                    deferred_figures=None, path_sequence=None):
    """Save the 2D Figure 1 (IBZ), 2 (spin-flip), and 3 (spin pattern) views.

    ``R_for_kpts`` is the selected spin-flip operation in the standardized
    primitive basis (the same matrix used to build the KPOINTS path), so the
    figure's k' and primed points coincide with the written path.
    ``flip_ops_for_plot`` (standardized primitive basis) drives the Figure 3
    spin-up/down coloring.  Returns the list of saved file paths.
    """
    import matplotlib.pyplot as plt

    b_matrix = np.array(centroid_result["b_matrix"], dtype=float)
    axis = int(centroid_result.get("vacuum_axis", 2))
    sc_type = centroid_result.get("sc_type", "BZ")
    lattice_title = _physical_lattice_title(centroid_result)
    # Do not infer a Cartesian kx/ky section from a fractional-axis index.
    # The standardized basis may permute Cartesian axes.
    bz_poly, basis = _bz_polygon_2d(b_matrix, axis, lattice_class=sc_type)
    reciprocal_axis_dirs = [
        _to_2d(vec, basis)
        for i, vec in enumerate(b_matrix)
        if i != axis
    ]

    kpath = centroid_result["band_kpath"]
    path_labels = {label for segment in kpath for label in segment}
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
    fig1, ax1 = _setup_2d_ax(
        f"2D BZ: {basename} ({lattice_title})", bz_poly
    )
    # Shade the IBZ itself: without it the wedge is invisible and the centroid
    # star looks misplaced relative to the plain high-symmetry path.  A
    # project-doubled IBZ (the square lattice's X_A at point group 4) carries an
    # extra sector beyond the conventional wedge, so fill it in a lighter color
    # and the genuine sub-polygon at full strength, as the spin-BZ figure does.
    ibz_polygon_frac = centroid_result.get("ibz_polygon_frac")
    if ibz_polygon_frac is not None and len(ibz_polygon_frac) >= 3:
        ibz_xy = np.array(
            [_to_2d(_cart_from_frac(point, b_matrix), basis)
             for point in ibz_polygon_frac],
            dtype=float,
        )
        ibz_polygon_labels_fig1 = centroid_result.get("ibz_polygon_labels")
        extra_flags_fig1 = None
        if ibz_polygon_labels_fig1 and \
                len(ibz_polygon_labels_fig1) == len(ibz_polygon_frac):
            flags = _doubled_ibz_extra_flags(ibz_polygon_labels_fig1)
            if any(flags):
                extra_flags_fig1 = flags
        if extra_flags_fig1 is not None:
            original_idx_fig1 = [j for j, is_extra in enumerate(extra_flags_fig1)
                                  if not is_extra]
            ax1.fill(ibz_xy[:, 0], ibz_xy[:, 1],
                     color=IBZ_FACE_COLORS["up_extra"], alpha=0.20, zorder=1)
            if len(original_idx_fig1) >= 3:
                sub_xy = ibz_xy[original_idx_fig1]
                ax1.fill(sub_xy[:, 0], sub_xy[:, 1],
                         color=IBZ_FACE_COLORS["up_main"], alpha=0.35, zorder=1)
        else:
            ax1.fill(ibz_xy[:, 0], ibz_xy[:, 1], color="salmon", alpha=0.20,
                     zorder=1)
    for start, end in kpath:
        if start in kpoints_cart and end in kpoints_cart:
            p1, p2 = kpoints_cart[start], kpoints_cart[end]
            ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], c="red", lw=3.0,
                     alpha=0.9, zorder=3)
    _draw_reciprocal_axes_2d(
        ax1, b_matrix, axis, basis, bz_poly, zorder=4,
    )
    _draw_labeled_points(ax1, kpoints_cart, "red", "darkred", label_color="black",
                         path_labels=path_labels, bz_poly=bz_poly,
                         avoid_dirs=reciprocal_axis_dirs)
    ax1.scatter(*centroid_xy, c="gold", marker="*", s=420, edgecolors="k",
                zorder=112, label=r"$k$")
    # Place the legend outside the axes so it cannot overlap plotted points or
    # labels and the BZ itself needs no extra padding.
    ax1.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=18,
              borderaxespad=0, frameon=True)
    fig1.tight_layout()
    fig1_path = os.path.join(output_dir, f"{basename}_2d_ibz_{sc_type}.png")
    saved.extend(_finish_2d_figure(
        fig1, fig1_path, save_pdf, deferred_figures=deferred_figures
    ))

    # ----- Figure 2: spin-up IBZ + spin-flip image -----
    mapped_cart_lines = {}
    for label, frac in coords.items():
        if str(label).startswith("_") or abs(np.array(frac)[axis]) >= 1e-4:
            continue
        mapped_frac = R_inv_T @ np.array(frac, dtype=float)
        mapped_frac[axis] = 0.0
        mapped_point = _to_2d(_cart_from_frac(mapped_frac, b_matrix), basis)
        mapped_cart_lines[prime_point_label(label)] = mapped_point

    if path_sequence is not None:
        figure2_path_labels = {
            alias
            for point in path_sequence
            if point is not None
            for alias in label_aliases(point[3])
        }
    else:
        figure2_path_labels = set(path_labels)
        figure2_path_labels.update(
            prime_point_label(label) for label in path_labels
        )
    up_display_points, down_display_points = _spinflip_display_points_2d(
        kpoints_cart,
        mapped_cart_lines,
        figure2_path_labels,
        spin_up_center=centroid_xy,
        spin_down_center=k_prime_xy,
        bz_span=max(float(np.max(np.ptp(bz_poly, axis=0))), 1e-8),
    )

    fig2, ax2 = _setup_2d_ax("2D spin-flip path connections", bz_poly)
    avoid_pts = (list(kpoints_cart.values()) + list(mapped_cart_lines.values())
                 + [centroid_xy, k_prime_xy])
    op_line_dir = _draw_op_visual_2d(ax2, R_for_kpts, b_matrix, basis, bz_poly,
                                     avoid_pts=avoid_pts)
    orig_labels = list(kpoints_cart.keys())
    mapped_labels = list(mapped_cart_lines.keys())
    orig_poly = np.array(list(kpoints_cart.values()), dtype=float)
    mapped_poly = np.array(list(mapped_cart_lines.values()), dtype=float)

    def _fill_doubled_aware(ax, poly, labels, main_color, extra_color):
        # Same lighter-color split as Figure 1.  The spin-flip image carries the
        # doubling on the same vertices but with primed names, which the "_A"
        # check does not recognize, so remove the prime before checking.
        if len(poly) < 3:
            return
        hull_idx = ConvexHull(poly).vertices
        hp = poly[hull_idx]
        hp_labels = [str(labels[i]).rstrip("'") for i in hull_idx]
        flags = _doubled_ibz_extra_flags(hp_labels)
        if any(flags):
            original_idx = [j for j, is_extra in enumerate(flags) if not is_extra]
            ax.fill(hp[:, 0], hp[:, 1], color=extra_color, alpha=0.20, zorder=1)
            if len(original_idx) >= 3:
                sub = hp[original_idx]
                ax.fill(sub[:, 0], sub[:, 1], color=main_color, alpha=0.35, zorder=1)
        else:
            ax.fill(hp[:, 0], hp[:, 1], color=main_color, alpha=0.20, zorder=1)

    _fill_doubled_aware(ax2, orig_poly, orig_labels,
                        IBZ_FACE_COLORS["up_main"], IBZ_FACE_COLORS["up_extra"])
    _fill_doubled_aware(ax2, mapped_poly, mapped_labels,
                        IBZ_FACE_COLORS["down_main"], IBZ_FACE_COLORS["down_extra"])
    if path_sequence is not None:
        for first, second, spin_side in generated_plain_path_segments(
                path_sequence):
            first_xy = _to_2d(
                _cart_from_frac(first, b_matrix), basis
            )
            second_xy = _to_2d(
                _cart_from_frac(second, b_matrix), basis
            )
            color = "navy" if spin_side == "down" else "red"
            ax2.plot(
                [first_xy[0], second_xy[0]],
                [first_xy[1], second_xy[1]],
                c=color, lw=4.0, alpha=0.9, zorder=50,
            )
    else:
        for start, end in kpath:
            if start in kpoints_cart and end in kpoints_cart:
                p1, p2 = kpoints_cart[start], kpoints_cart[end]
                ax2.plot([p1[0], p2[0]], [p1[1], p2[1]], c="red", lw=4.0,
                         alpha=0.9, zorder=50)
    shared_points = list(kpoints_cart.values()) + list(mapped_cart_lines.values())
    _draw_reciprocal_axes_2d(
        ax2, b_matrix, axis, basis, bz_poly, zorder=51,
    )
    shared_center, shared_scale = _label_offset(shared_points)
    _draw_labeled_points(ax2, up_display_points, "salmon", "darkred",
                         center=shared_center, offset_scale=shared_scale,
                         path_labels=figure2_path_labels, bz_poly=bz_poly,
                         avoid_dir=op_line_dir,
                         avoid_dirs=reciprocal_axis_dirs,
                         avoid_points=shared_points,
                         prime_label_color="navy", zorder=213)
    _draw_labeled_points(ax2, down_display_points, "cornflowerblue", "navy",
                         bz_poly=bz_poly, avoid_dir=op_line_dir,
                         center=shared_center, offset_scale=shared_scale,
                         path_labels=figure2_path_labels,
                         avoid_dirs=reciprocal_axis_dirs,
                         avoid_points=shared_points, zorder=213)
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
    ax2.scatter(*k_prime_xy, c="cornflowerblue", s=180, marker="D",
                edgecolors="k", linewidths=0.8, zorder=120, label=r"$k'$")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=18,
              borderaxespad=0, frameon=True)
    fig2.tight_layout()
    fig2_path = os.path.join(output_dir, f"{basename}_2d_spinflip_{sc_type}.png")
    saved.extend(_finish_2d_figure(
        fig2, fig2_path, save_pdf, deferred_figures=deferred_figures
    ))

    # ----- Figure 3: spin-up / spin-down BZ pattern -----
    fig3_path = os.path.join(output_dir, f"{basename}_2d_spinbz_{sc_type}.png")
    fig3_saved = _plot_spin_pattern_top_view_2d(
        centroid_result,
        R_for_kpts,
        flip_ops_for_plot,
        fig3_path,
        save_pdf=save_pdf,
        deferred_figures=deferred_figures,
    )
    if fig3_saved is not None:
        saved.extend(fig3_saved)

    return saved
