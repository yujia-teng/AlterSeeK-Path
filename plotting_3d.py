"""3D Brillouin-zone / IBZ / spin-flip figures (matplotlib).

Extracted from compute_centroid_hybrid.py (restructuring phase 1).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

from plotting_common import _get_bz_path_style, _save_figure, _print_saved_paths, _math_label
from symmetry import (
    _axis_bz_exit,
    _classify_spin_down_ops,
    _classify_spinflip_op,
    _doubled_ibz_extra_flags,
    _format_miller,
    _is_doubled_ibz_extra_label,
    _mirror_plane_bz_polygon,
    _perp_unit,
    _reduce_int_vector,
    _rotation_sense,
    _seekpath_label_to_internal,
)
from geometry import (
    _bz_kz_plane_outline,
    _get_ibz_frame_edges,
    _mapped_spin_hulls,
    _points_on_kz_plane,
    find_bz_exit,
)

IBZ_FACE_COLORS = {
    "up_main": "salmon",
    "up_extra": "#e98f8f",
    "down_main": "cornflowerblue",
    "down_extra": "#91b2e8",
}


class _Arrow3D(FancyArrowPatch):
    """FancyArrowPatch projected into 3D — paper-quality arrow from any view angle."""
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def _draw_ibz_faces_by_sector(ax, hull_pts, hull_simplices, hull_labels,
                              main_color, extra_color,
                              main_alpha=0.20, extra_alpha=0.14):
    points = np.array(hull_pts)
    labels = list(hull_labels) if hull_labels is not None else None
    extra_flags = (
        _doubled_ibz_extra_flags(labels) if labels is not None
        else [False] * len(points)
    )
    if len(extra_flags) != len(points) or not any(extra_flags):
        faces = [[points[s] for s in simplex] for simplex in hull_simplices]
        ax.add_collection3d(Poly3DCollection(
            faces, facecolor=main_color, alpha=main_alpha, edgecolor='none'))
        return

    full_faces = [[points[s] for s in simplex] for simplex in hull_simplices]
    ax.add_collection3d(Poly3DCollection(
        full_faces, facecolor=extra_color, alpha=extra_alpha, edgecolor='none'))

    original_points = np.array([
        point for point, is_extra in zip(points, extra_flags)
        if not is_extra
    ])
    if len(original_points) < 4:
        return
    try:
        original_hull = ConvexHull(original_points)
    except Exception:
        return

    original_faces = [
        [original_points[s] for s in simplex]
        for simplex in original_hull.simplices
    ]
    ax.add_collection3d(Poly3DCollection(
        original_faces, facecolor=main_color,
        alpha=main_alpha, edgecolor='none'))


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
                elev=14, azim=20, dashed_back=False):
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
        # Arrow outside BZ: FancyArrowPatch projected into 3D
        outside = vec - exit_pt
        arrow_frac = 0.45
        tip = exit_pt + outside * arrow_frac
        arrow = _Arrow3D(
            [exit_pt[0], tip[0]], [exit_pt[1], tip[1]], [exit_pt[2], tip[2]],
            arrowstyle='->', mutation_scale=28,
            color='black', lw=2.5, shrinkA=0, shrinkB=0, zorder=100,
        )
        ax.add_artist(arrow)
        ax.text(tip[0] + outside[0]*0.08, tip[1] + outside[1]*0.08,
                tip[2] + outside[2]*0.08, vec_labels[i],
                color='black', fontsize=24, fontweight='bold', zorder=101)
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


def plot_ibz(ax, kpoints_cart, kpath, display_labels, hull, centroid_cart,
             hull_pts=None, lattice_type=None, hull_labels=None):
    points_list = list(kpoints_cart.values())
    # Draw IBZ faces (skip for triclinic where hull is None)
    if hull is not None:
        face_points = np.array(hull_pts) if hull_pts is not None else np.array(points_list)
        _draw_ibz_faces_by_sector(
            ax, face_points, hull.simplices, hull_labels,
            IBZ_FACE_COLORS["up_main"], IBZ_FACE_COLORS["up_extra"])
        for p1, p2 in _get_ibz_frame_edges(face_points, hull.simplices, hull_labels):
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                    c='darkred', lw=1.8, alpha=0.85, zorder=10)
    for k1, k2 in kpath:
        if k1 in kpoints_cart and k2 in kpoints_cart:
            p1, p2 = kpoints_cart[k1], kpoints_cart[k2]
            style = _get_bz_path_style(lattice_type, k1, k2)
            ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]],
                    c=style["color"], ls=style["ls"], lw=style["lw"],
                    alpha=style["alpha"], zorder=50)
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
                _math_label(label),
                fontsize=22, color='black',
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


def _screen_xy(ax, pt3):
    """Project a 3D point to the axes' 2D projected (screen) coordinates."""
    x, y, _ = proj3d.proj_transform(pt3[0], pt3[1], pt3[2], ax.get_proj())
    return np.array([x, y])


def _best_label_anchor(ax, candidates, avoid_pts):
    """Pick the candidate 3D anchor whose projected screen position is farthest
    from every avoid point. Generic placement so a label (e.g. a mirror m_{hkl})
    clears the drawn high-symmetry points/labels in any case, not one structure.
    """
    candidates = np.asarray(candidates, dtype=float)
    if avoid_pts is None or len(avoid_pts) == 0:
        return candidates[0]
    avoid_screen = np.array([_screen_xy(ax, p) for p in avoid_pts])
    best, best_score = candidates[0], -np.inf
    for c in candidates:
        d = float(np.min(np.linalg.norm(avoid_screen - _screen_xy(ax, c), axis=1)))
        if d > best_score:
            best_score, best = d, c
    return best


def _mirror_label_candidates(poly, bz_radius):
    """Anchor candidates around a mirror-plane polygon: each rim vertex and edge
    midpoint pushed outward (in-plane) and lifted slightly, so the label sits
    just off the plane rim. The caller screen-tests these to avoid overlaps."""
    poly = np.asarray(poly, dtype=float)
    c0 = poly.mean(axis=0)
    margin = bz_radius * 0.18
    lift = np.array([0.0, 0.0, 1.0]) * bz_radius * 0.05
    n = len(poly)
    cands = []
    for j in range(n):
        for p in (poly[j], 0.5 * (poly[j] + poly[(j + 1) % n])):
            o = p - c0
            no = np.linalg.norm(o)
            o = o / no if no > 1e-9 else np.array([0.0, 0.0, 1.0])
            cands.append(p + o * margin + lift)
    return np.array(cands)


def _draw_op_visual(ax, R_cart, bz_loops, bz_radius, b_matrix,
                    avoid_pts=None, elev=None, azim=None):
    """
    Draw the geometric visual for the chosen spin-flip operation on ax.

      Mirror  (det=-1, tr=1)  : grey shaded plane through Γ + m_{hkl} label
      Cn rotation (det=+1)    : bold colored axis arrow + curved rotation arrow
      Inversion / Sn / identity : nothing drawn

    The mirror is labelled by its plane normal expressed as reduced integer
    (hkl) indices in the figure's reciprocal (b1,b2,b3) basis.

    elev/azim are the Matplotlib 3D view angles; they orient the curved
    rotation arrow so its back gap stays away from the camera.
    """
    op = _classify_spinflip_op(R_cart)
    op_type = op['type']

    if op_type == 'mirror' and op['axis'] is not None:
        poly = _mirror_plane_bz_polygon(op['axis'], bz_loops)
        if poly is not None and len(poly) >= 3:
            poly = np.asarray(poly, dtype=float)
            ax.add_collection3d(Poly3DCollection(
                [poly],
                alpha=0.30,
                facecolor='#606060',
                edgecolor='#303030',
                linewidth=0.8,
                zorder=5,
            ))

            # m_{hkl} label: mirror-plane normal in the reciprocal (b1,b2,b3)
            # basis, reduced to integers. Grey to match the plane, bold,
            # screen-projected so it never stacks/rotates in saved views.
            MIRROR_COLOR = os.environ.get('ALTERSEEK_OP_MIRROR_COLOR', '#606060')
            n = np.asarray(op['axis'], dtype=float)
            idx = _reduce_int_vector(n @ np.linalg.inv(b_matrix))
            m_label = _format_miller('m', idx)
            # Place the label at the rim anchor whose screen position is farthest
            # from the drawn high-symmetry points/labels (avoid_pts), so it never
            # lands on a BZ-boundary point such as K. General for all cases.
            candidates = _mirror_label_candidates(poly, bz_radius)
            label_pt = _best_label_anchor(ax, candidates, avoid_pts)
            mlx, mly, _ = proj3d.proj_transform(*label_pt, ax.get_proj())
            m_artist = ax.text2D(
                mlx, mly, m_label, transform=ax.transData,
                fontsize=17, fontweight='bold', color=MIRROR_COLOR, zorder=120,
                ha='center', va='center',
            )

            def _update_mirror_label():
                lx, ly, _ = proj3d.proj_transform(*label_pt, ax.get_proj())
                m_artist.set_position((lx, ly))

            ax._alterseek_arc_updater = _update_mirror_label

    elif op_type in ('rotation', 'rotoreflection') and op['axis'] is not None:
        axis = np.array(op['axis'], dtype=float)
        # Orient tip toward upper hemisphere for consistent appearance
        if axis[2] < -1e-6 or (abs(axis[2]) < 1e-6 and axis[1] < -1e-6):
            axis = -axis

        # Stand-out color so the rotation axis is easy to spot against the
        # red/navy IBZ and black BZ frame. Override with ALTERSEEK_OP_AXIS_COLOR.
        AXIS_COLOR = os.environ.get('ALTERSEEK_OP_AXIS_COLOR', '#00c853')

        # Axis exits BZ at bz_exit along +axis and bz_exit_neg along -axis.
        # Extend beyond the BZ so the arrow is always visibly longer than the BZ
        # (e.g. for flat hexagonal where c* << a*).
        bz_exit     = _axis_bz_exit( axis, bz_loops)
        bz_exit_neg = _axis_bz_exit(-axis, bz_loops)
        extension   = max(bz_radius * 0.35, bz_exit * 0.18)
        tip  =  axis * (bz_exit     + extension)
        base = -axis * (bz_exit_neg + extension)

        # Lower half (through the cell): thin dashed, so the axis is seen to pass
        # through Γ without overpowering the figure.
        ax.plot(
            [base[0], 0.0], [base[1], 0.0], [base[2], 0.0],
            color=AXIS_COLOR, lw=1.6, ls='--', alpha=0.55, zorder=6,
        )
        # Solid arrowhead at the tip (like the c_n axis arrow in the reference).
        head_len = bz_radius * 0.18
        shaft_pt = tip - axis * head_len
        # Upper shaft stops at the arrow base; the final segment is the arrow.
        ax.plot(
            [0.0, shaft_pt[0]], [0.0, shaft_pt[1]], [0.0, shaft_pt[2]],
            color=AXIS_COLOR, lw=2.8, ls='-', alpha=0.95, zorder=6,
        )
        axis_arrow = _Arrow3D(
            [shaft_pt[0], tip[0]], [shaft_pt[1], tip[1]], [shaft_pt[2], tip[2]],
            arrowstyle='-|>', mutation_scale=26,
            color=AXIS_COLOR, lw=2.8, shrinkA=0, shrinkB=0, zorder=8,
        )
        ax.add_artist(axis_arrow)

        # Curved rotation arrow: a 3D ring in the plane perpendicular to the axis
        # (wraps around it). Drawn once as a plain 3D object, so it rotates WITH the
        # cell. The 60° gap is placed at the BACK (away from the camera for the
        # starting view) and the arrowhead sits on the FRONT of the ring.
        order = op['order'] or 2
        # International/Bilbao numeral convention instead of Schoenflies: a
        # proper rotation Cn is just its order n; an improper rotoreflection
        # Sn is relabeled by its rotoinversion-axis order n-bar, using the
        # standard correspondence S3<->6bar, S4<->4bar, S6<->3bar (S1/S2 are
        # the separately-handled mirror/inversion cases and never reach here).
        if op_type == 'rotation':
            display_order = order
            improper = False
        else:
            display_order = {3: 6, 4: 4, 6: 3}.get(order, order)
            improper = True
        # Distinguish the rotation sense (n^+ vs n^-) by a +/- superscript,
        # measured about the drawn (upper-hemisphere) axis. None for order 2.
        sense = _rotation_sense(R_cart, axis)
        sgn = '' if (order < 3 or sense == 0) else ('+' if sense > 0 else '-')
        idx = _reduce_int_vector(np.asarray(axis) @ np.linalg.inv(np.asarray(b_matrix)))
        axis_sub = "".join(rf"\bar{{{abs(i)}}}" if i < 0 else f"{i}" for i in idx)
        digit = rf"\bar{{{display_order}}}" if improper else f"{display_order}"
        order_str = (rf"$\mathbf{{{digit}^{{{sgn}}}_{{{axis_sub}}}}}$" if sgn
                     else rf"$\mathbf{{{digit}_{{{axis_sub}}}}}$")
        u = _perp_unit(axis)
        v = np.cross(axis, u)
        v = v / np.linalg.norm(v)

        arc_center = tip - axis * (bz_radius * 0.05)
        r_arc = bz_radius * 0.20
        span  = np.radians(300.0)

        # Angle of the ring point nearest the camera (the front) for this view.
        elev = ax.elev if elev is None else elev
        azim = ax.azim if azim is None else azim
        el0, az0 = np.radians(elev), np.radians(azim)
        view0 = np.array([np.cos(el0) * np.cos(az0),
                          np.cos(el0) * np.sin(az0), np.sin(el0)])
        cu, cv = float(view0 @ u), float(view0 @ v)
        phi_front = np.arctan2(cv, cu) if (abs(cu) + abs(cv)) > 1e-9 else 0.0

        # Centre the 300° arc on the front, so the 60° gap lands at the back.
        theta = np.linspace(phi_front - span / 2, phi_front + span / 2, 121)
        arc_pts = (arc_center[None, :]
                   + r_arc * (np.cos(theta)[:, None] * u + np.sin(theta)[:, None] * v))
        arc_head_start = -10
        arc_line = ax.plot(arc_pts[:arc_head_start, 0],
                           arc_pts[:arc_head_start, 1],
                           arc_pts[:arc_head_start, 2],
                           color=AXIS_COLOR, lw=3.0, alpha=0.95, zorder=9)[0]

        # Draw only the final arc segment as an arrow. The plain arc stops at
        # the arrow base, so the stroke does not run through the arrowhead.
        arc_arrow = _Arrow3D(
            [arc_pts[arc_head_start, 0], arc_pts[-1, 0]],
            [arc_pts[arc_head_start, 1], arc_pts[-1, 1]],
            [arc_pts[arc_head_start, 2], arc_pts[-1, 2]],
            arrowstyle='->', mutation_scale=22,
            color=AXIS_COLOR, lw=3.0, shrinkA=0, shrinkB=0, zorder=9,
        )
        ax.add_artist(arc_arrow)

        # Bold label, kept general for any axis orientation:
        #  - vertical axes (C6/S3/C2∥c): push the label straight up past the tip.
        #  - in-plane / tilted axes: lift modestly in +z so it clears the b1/b2
        #    labels without floating far above the figure.
        if abs(axis[2]) > 0.9:
            # Tie the gap to the axis-direction BZ extent (bz_exit), not the wide
            # in-plane bz_radius, so on flat cells (hexagonal c* << a*) the label
            # hugs the arrowhead instead of floating far above. Small bz_radius
            # floor keeps it clear of the arrowhead when bz_exit is tiny.
            label_pt = tip + axis * max(bz_exit * 0.30, bz_radius * 0.05)
        else:
            label_pt = tip + axis * bz_radius * 0.10 + np.array([0.0, 0.0, 1.0]) * bz_radius * 0.16
        label_x, label_y, _ = proj3d.proj_transform(*label_pt, ax.get_proj())
        label_artist = ax.text2D(
            label_x, label_y, order_str, transform=ax.transData,
            fontsize=17, fontweight='bold', color=AXIS_COLOR, zorder=120,
            ha='center', va='center',
        )

        def _update_camera_facing_arc():
            el_now, az_now = np.radians(ax.elev), np.radians(ax.azim)
            view_now = np.array([
                np.cos(el_now) * np.cos(az_now),
                np.cos(el_now) * np.sin(az_now),
                np.sin(el_now),
            ])
            cu_now, cv_now = float(view_now @ u), float(view_now @ v)
            phi_now = (
                np.arctan2(cv_now, cu_now)
                if (abs(cu_now) + abs(cv_now)) > 1e-9 else 0.0
            )
            theta_now = np.linspace(phi_now - span / 2, phi_now + span / 2, 121)
            arc_pts_now = (
                arc_center[None, :]
                + r_arc * (
                    np.cos(theta_now)[:, None] * u
                    + np.sin(theta_now)[:, None] * v
                )
            )
            arc_line.set_data_3d(
                arc_pts_now[:arc_head_start, 0],
                arc_pts_now[:arc_head_start, 1],
                arc_pts_now[:arc_head_start, 2],
            )
            arc_arrow._verts3d = (
                [arc_pts_now[arc_head_start, 0], arc_pts_now[-1, 0]],
                [arc_pts_now[arc_head_start, 1], arc_pts_now[-1, 1]],
                [arc_pts_now[arc_head_start, 2], arc_pts_now[-1, 2]],
            )
            label_x_now, label_y_now, _ = proj3d.proj_transform(
                *label_pt, ax.get_proj()
            )
            label_artist.set_position((label_x_now, label_y_now))

        ax._alterseek_arc_updater = _update_camera_facing_arc


def _connect_arc_view_follow(fig, ax):
    """
    Make the rotation arc follow the camera in an interactive window: whenever
    the user drags to a new view, re-orient the arc so its arrowhead stays on
    the near (visible) side instead of swinging behind the BZ. No-op if the axis
    has no rotation arc (mirror/identity/inversion ops).
    """
    updater = getattr(ax, '_alterseek_arc_updater', None)
    if updater is None:
        return
    state = {'view': (ax.elev, ax.azim)}

    def _on_move(event):
        cur = (ax.elev, ax.azim)
        if cur != state['view']:
            state['view'] = cur
            updater()
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('motion_notify_event', _on_move)


def plot_spin_flip_figure(b_matrix, bz_loops, bz_center, bz_span,
                          kpoints_data, ibz_kpoints_frac,
                          hull_pts, hull_simplices,
                          centroid_frac, R,
                          output_path, elev=14, azim=20, show_plot=True,
                          block=True, path_sequence=None, R_cart=None,
                          defer_show=False, unique_ops=None):
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

    # Original IBZ high-sym points (from lattice_kpoints, not KPOINTS file)
    ibz_orig = {}
    for lbl, frac in ibz_kpoints_frac.items():
        if not lbl.startswith('_'):
            f = np.array(frac)
            ibz_orig[lbl] = f[0] * b1 + f[1] * b2 + f[2] * b3
    if not ibz_orig and kpoints_data:
        # Defensive fallback for older callers: rebuild high-symmetry labels
        # from the generated KPOINTS rows.  Keep k/k' as guide points only.
        for row in kpoints_data:
            if len(row) < 4:
                continue
            lbl = str(row[3])
            if lbl in ('k', "k'") or lbl.endswith("'"):
                continue
            if lbl not in ibz_orig:
                ibz_orig[lbl] = row[0] * b1 + row[1] * b2 + row[2] * b3

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
    hull_labels = (
        list(ibz_kpoints_frac.keys())
        if hull_pts is not None and len(ibz_kpoints_frac) == len(hull_pts)
        else None
    )

    # Filter threshold: skip connections shorter than 5% of BZ radius
    bz_radius = np.max(np.linalg.norm(np.vstack(bz_loops), axis=1))
    threshold = 0.05 * bz_radius

    def _draw(ax):
        # Mirror plane or rotation axis visual (drawn first so IBZ sits on top).
        # Pass the labeled high-symmetry points + reciprocal-axis tips so the
        # mirror m_{hkl} label can be placed clear of them (screen-space).
        avoid = list(ibz_orig.values()) + list(ibz_mapped.values())
        avoid += [b1 * 0.5, b2 * 0.5, b3 * 0.5]
        avoid_pts = np.array(avoid) if avoid else None
        _draw_op_visual(ax, R_cart, bz_loops, bz_radius, b_matrix,
                        avoid_pts=avoid_pts, elev=ax.elev, azim=ax.azim)

        # Spin-up IBZ: use the same curated HPKOT/project hull as Figure 1.
        up_pts, up_simplices = hull_pts, hull_simplices
        if up_pts is not None and up_simplices is not None:
            up_pts = np.array(up_pts)
            _draw_ibz_faces_by_sector(
                ax, up_pts, up_simplices, hull_labels,
                IBZ_FACE_COLORS["up_main"], IBZ_FACE_COLORS["up_extra"])
            for p1, p2 in _get_ibz_frame_edges(up_pts, up_simplices, hull_labels):
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                        c='darkred', lw=1.8, alpha=0.85, zorder=10)

        # Spin-down IBZ: the same Figure-1 hull mapped by the chosen spin-flip op.
        down_pts, down_simplices = mapped_hull_pts, hull_simplices
        if down_pts is not None and down_simplices is not None:
            down_pts = np.array(down_pts)
            _draw_ibz_faces_by_sector(
                ax, down_pts, down_simplices, hull_labels,
                IBZ_FACE_COLORS["down_main"], IBZ_FACE_COLORS["down_extra"])
            for p1, p2 in _get_ibz_frame_edges(down_pts, down_simplices, hull_labels):
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                        c='navy', lw=1.8, alpha=0.85, zorder=10)

        # Generated path segments (red = spin-up side, navy = spin-down side).
        # Use path_sequence for labels/order only. The coordinates are looked up from
        # ibz_orig/ibz_mapped_lines so Figure 2 uses the same high-symmetry convention
        # as Figure 1, even when labels share names across different bases.
        def _path_point(label):
            if label in ('k', "k'"):
                return None
            is_prime = label.endswith("'")
            raw_base = label.rstrip("'")
            base = _seekpath_label_to_internal(raw_base)
            if is_prime:
                pt = ibz_mapped_lines.get(raw_base + "'")
                return pt if pt is not None else ibz_mapped_lines.get(base + "'")
            pt = ibz_orig.get(raw_base)
            return pt if pt is not None else ibz_orig.get(base)

        def _is_gamma_label(label):
            base = str(label).rstrip("'")
            return base in ('Γ', 'GAMMA', 'ŚŁ')

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
                gamma_path = _is_gamma_label(la) or _is_gamma_label(lb)
                col = 'navy' if (
                    (a_prime and b_prime) or
                    (gamma_path and (a_prime or b_prime))
                ) else 'red'
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

        # Label helpers — use combined span/center so offset scale matches Fig 1
        _all_pts = np.array(list(ibz_orig.values()) + list(ibz_mapped.values()))
        _lbl_center = np.mean(_all_pts, axis=0) if len(_all_pts) else np.zeros(3)
        # Base offset scale on original IBZ only — mapped sector can be far away
        # (e.g. C2 flips IBZ to opposite side), which would bloat _lbl_span and
        # push all labels far from their points.
        _orig_pts = np.array(list(ibz_orig.values())) if ibz_orig else _all_pts
        _lbl_span = max(np.max(np.ptp(_orig_pts, axis=0)), 1e-8) if len(_orig_pts) else 1.0
        _off_sc = _lbl_span * 0.10

        def _label_pts(pts_dict, color, edgecolor):
            for lbl, hpt in pts_dict.items():
                ax.scatter(*hpt, c=color, s=60, zorder=110,
                           edgecolors=edgecolor, linewidths=0.5)
                direction = hpt - _lbl_center
                nd = np.linalg.norm(direction)
                off = direction / nd * _off_sc if nd > 1e-8 else np.array([0, 0, _off_sc])
                ax.text(*(hpt + off), _math_label(lbl), fontsize=22, color=edgecolor,
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
        ax.legend(loc='upper right', fontsize=18)

    fig, ax = setup_3d_ax("Spin-flip path connections",
                          bz_loops, b_matrix, bz_center, bz_span,
                          elev=elev, azim=azim, dashed_back=False)
    _draw(ax)
    # In an interactive window, keep the rotation arc facing the camera on drag.
    if show_plot:
        _connect_arc_view_follow(fig, ax)
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
            saved_paths = _save_figure(fig_save, output_path, dpi=300, bbox_inches='tight')
            plt.close(fig_save)
            plt.close(fig)
            _print_saved_paths(saved_paths)
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
    saved_paths = _save_figure(fig_save, output_path, dpi=300, bbox_inches='tight')
    plt.close(fig_save)
    _print_saved_paths(saved_paths)
    return display_fig


def plot_spin_bz_figure(b_matrix, bz_loops, bz_center, bz_span,
                        unique_ops, centroid_cart,
                        hull_pts, hull_simplices,
                        R, output_path,
                        flip_ops_frac=None,
                        preserve_ops_frac=None,
                        elev=14, azim=20, show_plot=True,
                        defer_show=False, z0=0.0,
                        show_helper_plane=True,
                        hull_labels=None):
    """
    Show Figure 1's IBZ hull mapped by spin-preserving/spin-flipping operations.

    The outer BZ remains the structural BZ. Each colored copy is exactly the
    Figure 1 hull; reduced spin symmetry can therefore leave part of the
    structural BZ uncolored.
      - RED   (salmon)         : spin-up  regions
      - BLUE  (cornflowerblue) : spin-down regions
    """
    if hull_pts is None or hull_simplices is None or not len(unique_ops):
        print("[Note] Skipping spin-BZ figure (no hull or symmetry ops available).")
        return

    hull_pts = np.array(hull_pts)
    centroid_cart = np.array(centroid_cart)
    hull_simplices_arr = np.array(hull_simplices)

    mapped_spin_hulls = None
    if (
        preserve_ops_frac is not None and len(preserve_ops_frac)
        and flip_ops_frac is not None and len(flip_ops_frac)
    ):
        mapped_spin_hulls = _mapped_spin_hulls(
            b_matrix, hull_pts, hull_simplices_arr,
            preserve_ops_frac, flip_ops_frac
        )
        if mapped_spin_hulls is None:
            print("[Warning] Skipping spin-BZ figure: inconsistent spin-operation classes.")
            return
    else:
        # Compatibility fallback for workflows that supply no preserve file.
        spin_down_mask = _classify_spin_down_ops(
            b_matrix, unique_ops, centroid_cart, R, flip_ops_frac
        )

    def _draw(ax):
        if show_helper_plane:
            draw_kz0_helper_plane(ax, bz_loops, z0=z0)

        if mapped_spin_hulls is not None:
            cells_to_draw = mapped_spin_hulls
        else:
            cells_to_draw = [
                ((g @ hull_pts.T).T, hull_simplices_arr, bool(spin_down_mask[i]))
                for i, g in enumerate(unique_ops)
            ]
        for cell_pts, cell_simplices, is_down in cells_to_draw:
            alpha = 0.2
            main_key = "down_main" if is_down else "up_main"
            extra_key = "down_extra" if is_down else "up_extra"
            _draw_ibz_faces_by_sector(
                ax, cell_pts, cell_simplices, hull_labels,
                IBZ_FACE_COLORS[main_key], IBZ_FACE_COLORS[extra_key],
                main_alpha=alpha, extra_alpha=0.10)
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
            saved_paths = _save_figure(fig_save, output_path, dpi=300, bbox_inches='tight')
            plt.close(fig_save)
            plt.close(fig)
            _print_saved_paths(saved_paths)
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
    saved_paths = _save_figure(fig, output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    _print_saved_paths(saved_paths)
    return display_fig


def draw_kz0_helper_plane(ax, bz_loops, z0=0.0, pad=0.08):
    """Draw the kz=0 BZ-section outline without changing any BZ geometry."""
    outline = _bz_kz_plane_outline(bz_loops, z0=z0)
    if outline is not None:
        outline3d = np.column_stack([
            outline[:, 0],
            outline[:, 1],
            np.full(len(outline), z0),
        ])
        closed_outline = np.vstack([outline3d, outline3d[0]])
        ax.plot(closed_outline[:, 0], closed_outline[:, 1], closed_outline[:, 2],
                color='#3f5268', lw=3.0, ls='--', alpha=0.95, zorder=90)


def draw_projected_reciprocal_axes(ax, b_matrix, bz_loops, z0=0.0):
    """Project b1, b2, b3 onto kz=0 and draw the projected in-plane arrows."""
    all_pts = np.vstack([np.array(loop, dtype=float) for loop in bz_loops])
    span_xy = np.ptp(all_pts[:, :2], axis=0)
    span = max(float(np.max(span_xy)), 1e-8)
    target = 0.60 * span
    origin = np.array([0.0, 0.0])
    colors = ['#202020', '#202020', '#202020']
    labels = [r'$\mathbf{b}_1$', r'$\mathbf{b}_2$', r'$\mathbf{b}_3$']
    outline = _bz_kz_plane_outline(bz_loops, z0=z0)

    def _ray_outline_exit(unit_vec):
        if outline is None or len(outline) < 3:
            return None
        t_hits = []
        closed = np.vstack([outline, outline[0]])
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

    for vec, label, color in zip(np.array(b_matrix, dtype=float), labels, colors):
        projected = np.array([vec[0], vec[1]], dtype=float)
        length = float(np.linalg.norm(projected))
        if length < 1e-10:
            ax.scatter(origin[0], origin[1], s=36, c=color, zorder=220)
            ax.text(origin[0] + 0.025 * span, origin[1] + 0.025 * span,
                    label, fontsize=24, fontweight='bold',
                    ha='left', va='bottom', color=color, zorder=221,
                    clip_on=False)
            continue

        unit = projected / length
        exit_t = _ray_outline_exit(unit)
        exit_len = exit_t if exit_t is not None else target * 0.65
        exit_pt = origin + unit * min(exit_len, target * 0.86)
        end = origin + unit * max(target, exit_len * 1.16)
        ax.plot([origin[0], exit_pt[0]], [origin[1], exit_pt[1]],
                color=color, ls=':', lw=1.5, alpha=0.6, zorder=219,
                clip_on=False)
        ann = ax.annotate(
            '',
            xy=end, xytext=exit_pt,
            arrowprops=dict(arrowstyle='->', color=color, lw=2.2,
                            mutation_scale=28, shrinkA=0, shrinkB=0),
            zorder=220,
            annotation_clip=False,
        )
        ann.set_clip_on(False)
        if ann.arrow_patch is not None:
            ann.arrow_patch.set_clip_on(False)
        offset = projected / length * (0.04 * span)
        ax.text(end[0] + offset[0], end[1] + offset[1], label,
                fontsize=24, fontweight='bold', ha='center', va='center',
                color=color, zorder=221, clip_on=False)


def plot_spin_bz_top_view_figure(b_matrix, bz_loops,
                                 unique_ops, centroid_cart,
                                 hull_pts, hull_simplices,
                                 R, output_path,
                                 flip_ops_frac=None,
                                 preserve_ops_frac=None,
                                 show_plot=True, defer_show=False,
                                 z0=0.0, show_title=False,
                                 show_projected_axes=True,
                                 show_legend=False,
                                 hull_labels=None):
    """Draw a darker top-view kz=0 slice of the Figure 3 spin-BZ coloring."""
    if hull_pts is None or hull_simplices is None or not len(unique_ops):
        print("[Note] Skipping spin-BZ top-view figure (no hull or symmetry ops available).")
        return

    hull_pts = np.array(hull_pts, dtype=float)
    hull_simplices_arr = np.array(hull_simplices, dtype=int)
    centroid_cart = np.array(centroid_cart, dtype=float)
    mapped_spin_hulls = None
    if (
        preserve_ops_frac is not None and len(preserve_ops_frac)
        and flip_ops_frac is not None and len(flip_ops_frac)
    ):
        mapped_spin_hulls = _mapped_spin_hulls(
            b_matrix, hull_pts, hull_simplices_arr,
            preserve_ops_frac, flip_ops_frac
        )
        if mapped_spin_hulls is None:
            print("[Warning] Skipping spin-BZ top-view figure: inconsistent spin-operation classes.")
            return
    else:
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

    if mapped_spin_hulls is not None:
        cells_to_draw = mapped_spin_hulls
    else:
        cells_to_draw = [
            ((g @ hull_pts.T).T, hull_simplices_arr, bool(spin_down_mask[i]))
            for i, g in enumerate(unique_ops)
        ]
    for cell_pts, cell_simplices, is_down in cells_to_draw:
        poly = _points_on_kz_plane(cell_pts, cell_simplices, z0=section_z)
        if poly is None:
            continue

        has_extra = (
            hull_labels is not None
            and any(_is_doubled_ibz_extra_label(lbl) for lbl in hull_labels)
        )
        if is_down:
            color = '#1f4e9e'
            extra_color = IBZ_FACE_COLORS["down_extra"]
            label = 'spin-down' if not down_labeled else None
            down_labeled = True
        else:
            color = '#b22222'
            extra_color = IBZ_FACE_COLORS["up_extra"]
            label = 'spin-up' if not up_labeled else None
            up_labeled = True

        closed = np.vstack([poly, poly[0]])
        if has_extra:
            ax.fill(poly[:, 0], poly[:, 1], facecolor=extra_color, alpha=0.46,
                    edgecolor='none', label=label)
            original_pts = np.array([
                point for point, lbl in zip(cell_pts, hull_labels)
                if not _is_doubled_ibz_extra_label(lbl)
            ])
            if len(original_pts) >= 4:
                try:
                    original_hull = ConvexHull(original_pts)
                    original_poly = _points_on_kz_plane(
                        original_pts, original_hull.simplices, z0=section_z)
                except Exception:
                    original_poly = None
                if original_poly is not None:
                    ax.fill(original_poly[:, 0], original_poly[:, 1],
                            facecolor=color, alpha=0.68, edgecolor='none')
        else:
            ax.fill(poly[:, 0], poly[:, 1], facecolor=color, alpha=0.68,
                    edgecolor='none', label=label)
        ax.plot(closed[:, 0], closed[:, 1], color=color, lw=0.9, alpha=0.95)

    outline = _bz_kz_plane_outline(bz_loops, z0=z0)
    if outline is not None:
        closed = np.vstack([outline, outline[0]])
        ax.plot(closed[:, 0], closed[:, 1], color='black', lw=2.0, label='BZ boundary')

    if show_projected_axes:
        draw_projected_reciprocal_axes(ax, b_matrix, bz_loops, z0=z0)

    ax.set_aspect('equal', adjustable='box')
    if show_title:
        ax.set_title(r'Spin-up / Spin-down BZ top view ($k_z = 0$)', fontsize=18)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()
    if show_legend:
        ax.legend(loc='upper right', fontsize=12)
    fig.tight_layout()

    display_fig = fig if show_plot and defer_show else None
    if display_fig is not None:
        def _save_after_show(fig=fig):
            saved_paths = _save_figure(fig, output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            _print_saved_paths(saved_paths)
        display_fig._alterseek_save_after_show = _save_after_show
        return display_fig

    if show_plot and not defer_show:
        plt.show()
    saved_paths = _save_figure(fig, output_path, dpi=300, bbox_inches='tight')
    _print_saved_paths(saved_paths)
    plt.close(fig)
    return display_fig
