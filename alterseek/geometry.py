"""Brillouin-zone / IBZ hull / centroid geometry.

Extracted from compute_centroid_hybrid.py (restructuring phase 3).
"""
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection, Voronoi
import sympy as sp
from .symmetry import _is_doubled_ibz_extra_label
from .lattice_kpoints import LATTICE_DATA


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
