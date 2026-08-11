"""Point-operation classification, 2D spin-flip filters, Laue/type mapping.

Extracted from compute_centroid_hybrid.py (restructuring phase 2). Leaf module:
depends only on numpy + scipy.
"""
import numpy as np
from scipy.spatial import ConvexHull


NO_ALTERMAGNETISM_LAUE_GROUPS = {'-1', '-3', 'm-3'}


def _is_doubled_ibz_extra_label(label):
    label = str(label)
    suffix = label.rsplit("_", 1)[-1] if "_" in label else ""
    return suffix.endswith("A") or label.endswith("_2")


def _doubled_ibz_extra_flags(hull_labels):
    labels = [str(label) for label in hull_labels]
    # Project-doubled sectors always include a copied _A label, while ordinary HPKOT _2 labels must not trigger sector coloring.
    if not any(
        "_" in label and label.rsplit("_", 1)[-1].endswith("A")
        for label in labels
    ):
        return [False] * len(labels)
    return [_is_doubled_ibz_extra_label(label) for label in labels]


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


_LAUE_GROUP_BY_SPACEGROUP_RANGE = (
    (2, '-1'), (15, '2/m'), (74, 'mmm'), (88, '4/m'), (142, '4/mmm'),
    (148, '-3'), (167, '-3m'), (176, '6/m'), (194, '6/mmm'),
    (206, 'm-3'), (230, 'm-3m'),
)

# Verified to reproduce spglib's pointgroup_international for all 230 groups.
_POINT_GROUP_BY_SPACEGROUP_RANGE = (
    (1, '1'), (2, '-1'), (5, '2'), (9, 'm'), (15, '2/m'), (24, '222'),
    (46, 'mm2'), (74, 'mmm'), (80, '4'), (82, '-4'), (88, '4/m'), (98, '422'),
    (110, '4mm'), (122, '-42m'), (142, '4/mmm'), (146, '3'), (148, '-3'),
    (155, '32'), (161, '3m'), (167, '-3m'), (173, '6'), (174, '-6'),
    (176, '6/m'), (182, '622'), (186, '6mm'), (190, '-6m2'), (194, '6/mmm'),
    (199, '23'), (206, 'm-3'), (214, '432'), (220, '-43m'), (230, 'm-3m'),
)


def _by_spacegroup_range(spacegroup_number, table):
    try:
        number = int(spacegroup_number)
    except (TypeError, ValueError):
        return None
    if not 1 <= number <= 230:
        return None
    for upper, value in table:
        if number <= upper:
            return value
    return None


def point_group_from_spacegroup_number(spacegroup_number):
    """Return the crystallographic point group of a space-group number.

    Companion to laue_group_from_spacegroup_number, used for symmetry known as
    a group rather than as coordinates (notably the SSG's G0). Printing the
    point group alongside the Laue group makes the reduction visible -- mm2
    reduces to mmm -- instead of leaving the reader to supply it.
    """
    return _by_spacegroup_range(spacegroup_number, _POINT_GROUP_BY_SPACEGROUP_RANGE)


def laue_group_from_spacegroup_number(spacegroup_number):
    """Return the Laue group of an international space-group number.

    Used for symmetry already known as a space group rather than as a set of
    coordinates -- notably the spatial part G0 of a spin space group, which
    FindSpinGroup reports directly. Reading the Laue group off the group itself
    avoids re-detecting symmetry from atomic positions, where the answer
    depends on a tolerance.
    """
    return _by_spacegroup_range(spacegroup_number, _LAUE_GROUP_BY_SPACEGROUP_RANGE)


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


def _seekpath_label_to_internal(label):
    if label == 'GAMMA':
        return '\u0393'
    return label


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


def slab_plane_normal_cartesian(lattice, vacuum_axis):
    """Return the physical slab-plane normal from a direct lattice.

    ``vacuum_axis`` identifies the submitted lattice vector containing the
    vacuum. The physical plane is defined by the other two direct vectors; its
    Cartesian normal remains meaningful after any cell-basis permutation or
    unimodular change.
    """
    lattice = np.asarray(lattice, dtype=float)
    if lattice.shape != (3, 3):
        raise ValueError("lattice must be a 3x3 matrix")
    if vacuum_axis not in (0, 1, 2):
        raise ValueError("vacuum axis must be 0, 1, or 2")
    in_plane = [index for index in range(3) if index != vacuum_axis]
    normal = np.cross(lattice[in_plane[0]], lattice[in_plane[1]])
    norm = np.linalg.norm(normal)
    if not np.isfinite(norm) or norm <= 1e-12:
        raise ValueError("submitted in-plane lattice vectors are collinear")
    return normal / norm


def reciprocal_operation_cartesian(R, b_matrix):
    """Convert real-space fractional R to its Cartesian reciprocal action."""
    R = np.asarray(R, dtype=float)
    b_matrix = np.asarray(b_matrix, dtype=float)
    if R.shape != (3, 3) or b_matrix.shape != (3, 3):
        raise ValueError("operation and reciprocal basis must both be 3x3")
    b_transpose = b_matrix.T
    return b_transpose @ np.linalg.inv(R).T @ np.linalg.inv(b_transpose)


def _cartesian_plane_basis(plane_normal):
    normal = np.asarray(plane_normal, dtype=float)
    if normal.shape != (3,) or not np.all(np.isfinite(normal)):
        raise ValueError("plane normal must contain three finite components")
    norm = np.linalg.norm(normal)
    if norm <= 1e-12:
        raise ValueError("plane normal cannot be zero")
    normal = normal / norm
    first = _perp_unit(normal)
    second = np.cross(normal, first)
    second /= np.linalg.norm(second)
    return normal, np.column_stack((first, second))


def keeps_2d_plane_cartesian(R, b_matrix, plane_normal, tol=1e-6):
    """True when R preserves the physical Cartesian slab plane."""
    operation = reciprocal_operation_cartesian(R, b_matrix)
    normal, plane_basis = _cartesian_plane_basis(plane_normal)
    mapped_plane = operation @ plane_basis
    return np.all(np.abs(normal @ mapped_plane) < tol)


def is_trivial_2d_spin_flip_cartesian(
    R,
    b_matrix,
    plane_normal,
    tol=1e-6,
):
    """True when R restricts to +I or -I on the physical slab plane."""
    operation = reciprocal_operation_cartesian(R, b_matrix)
    _, plane_basis = _cartesian_plane_basis(plane_normal)
    block = plane_basis.T @ operation @ plane_basis
    return (np.allclose(block, np.eye(2), atol=tol)
            or np.allclose(block, -np.eye(2), atol=tol))


def is_valid_2d_spin_flip_cartesian(R, b_matrix, plane_normal, tol=1e-6):
    """True for a plane-preserving, nontrivial physical 2D spin flip."""
    return (
        keeps_2d_plane_cartesian(R, b_matrix, plane_normal, tol)
        and not is_trivial_2d_spin_flip_cartesian(
            R, b_matrix, plane_normal, tol
        )
    )


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
    Return a flat rectangle (4,3) representing the mirror plane n·k=0, centered
    on Gamma (the plane always passes through the origin) and sized to the
    axis-aligned bounding box (in-plane) of where the plane actually cuts the
    BZ edges. This reads as a plain textbook mirror-plane rectangle rather than
    tracing the BZ's own cross-section outline at that cut (which for e.g. a
    horizontal mirror in a hexagonal BZ would itself be a hexagon), while
    staying tightly sized to the true local cut extent instead of a uniform
    whole-BZ radius. Returns None if fewer than 3 intersection points.
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

    u = _perp_unit(n)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    uv = np.column_stack([u, v])
    coords_2d = pts @ uv
    umin, vmin = coords_2d.min(axis=0)
    umax, vmax = coords_2d.max(axis=0)
    corners_2d = [(umax, vmax), (umin, vmax), (umin, vmin), (umax, vmin)]
    return np.array([c[0] * u + c[1] * v for c in corners_2d])


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
        # Express axis components in the reciprocal b1/b2/b3 frame used by the mirror normal so both labels are read against the drawn axes.
        hkl_raw = axis @ np.linalg.inv(np.asarray(b_matrix))
        # Match the reported indices' first-nonzero-positive convention and measure rotation sense about that same reported axis direction.
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
