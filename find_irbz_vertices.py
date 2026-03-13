#!/usr/bin/env python3
"""Find the complete IRBZ vertices for MCLC1 (C2/m) and identify missing ones."""
import numpy as np
from scipy.spatial import Voronoi, ConvexHull
import seekpath
from pymatgen.core import Structure
import spglib
import warnings
warnings.filterwarnings('ignore')

struct = Structure.from_file('../monoclinic/C2-m.vasp')
cell = struct.lattice.matrix.tolist()
positions = struct.frac_coords.tolist()
numbers = [s.Z for s in struct.species]
sp_result = seekpath.get_path((cell, positions, numbers), with_time_reversal=True)
prim_latt = np.array(sp_result['primitive_lattice'])
prim_pos = np.array(sp_result['primitive_positions'])
prim_types = sp_result['primitive_types']
b_matrix = np.array(sp_result['reciprocal_primitive_lattice'])
b1, b2, b3 = b_matrix
b_inv = np.linalg.inv(b_matrix)

spg_cell = (prim_latt, prim_pos, prim_types)
dataset = spglib.get_symmetry_dataset(spg_cell)
print(f"Space group: {dataset.number} ({dataset.international})")

# C2/m point group in fractional reciprocal coords:
# R0: E   (k1,k2,k3) -> (k1,k2,k3)
# R1: i   (k1,k2,k3) -> (-k1,-k2,-k3)
# R2: C2  (k1,k2,k3) -> (k2,k1,-k3)
# R3: m   (k1,k2,k3) -> (-k2,-k1,k3)
#
# Group order = 4, so IRBZ = BZ/4
#
# Two independent constraints:
# 1. From C2 (R2): k1 <= k2  (when k1=k2, tiebreaker: k3 >= 0)
# 2. From mirror (R3): k1+k2 >= 0

# BZ construction
grid = np.array(np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1])).T.reshape(-1,3)
points_cart = grid @ b_matrix
vor = Voronoi(points_cart)
origin_idx = 13

bz_faces = []
for i, pair in enumerate(vor.ridge_points):
    if origin_idx not in pair:
        continue
    idx = vor.ridge_vertices[i]
    if -1 in idx:
        continue
    pts = vor.vertices[idx]
    center = np.mean(pts, axis=0)
    neighbor = points_cart[pair[0] if pair[1] == origin_idx else pair[1]]
    normal = neighbor - points_cart[origin_idx]
    normal /= np.linalg.norm(normal)
    ref = np.array([0.,0.,1.]) if abs(normal[2]) < 0.9 else np.array([0.,1.,0.])
    u = np.cross(normal, ref); u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    angles = np.arctan2((pts-center)@v, (pts-center)@u)
    ordered = pts[np.argsort(angles)]
    bz_faces.append(ordered)

bz_edges = []
for face in bz_faces:
    n = len(face)
    for j in range(n):
        bz_edges.append((face[j], face[(j+1) % n]))


def inside_bz(pt):
    G_all = grid @ b_matrix
    for G in G_all:
        if np.linalg.norm(G) < 1e-10:
            continue
        if np.dot(pt, G) > np.dot(G, G) / 2 + 1e-8:
            return False
    return True


# Two symmetry plane normals in Cartesian
n_C2 = b_inv[:, 1] - b_inv[:, 0]       # positive when k2 > k1
n_mirror = b_inv[:, 0] + b_inv[:, 1]    # positive when k1+k2 > 0
sym_normals = [n_C2, n_mirror]


def in_fundamental_domain(pt_cart):
    pf = pt_cart @ b_inv
    k1, k2, k3 = pf
    if k1 > k2 + 1e-8:
        return False
    if k1 + k2 < -1e-8:
        return False
    return True


def line_plane_intersect(p1, p2, normal):
    d1 = np.dot(p1, normal)
    d2 = np.dot(p2, normal)
    if abs(d2 - d1) < 1e-12:
        return None
    t = -d1 / (d2 - d1)
    if -1e-10 <= t <= 1 + 1e-10:
        return p1 + t * (p2 - p1)
    return None


candidates = [np.zeros(3)]  # Gamma

# 1. BZ vertices in FD
bz_verts = set()
for face in bz_faces:
    for v in face:
        bz_verts.add(tuple(np.round(v, 10)))
bz_verts = np.array(list(bz_verts))
for v in bz_verts:
    if in_fundamental_domain(v):
        candidates.append(v.copy())

# 2. BZ edge / symmetry plane intersections
for p1, p2 in bz_edges:
    for n in sym_normals:
        pt = line_plane_intersect(p1, p2, n)
        if pt is not None and inside_bz(pt) and in_fundamental_domain(pt):
            candidates.append(pt)

# 3. Line = intersection of 2 sym planes, find where it crosses BZ faces
line_dir = np.cross(n_C2, n_mirror)
if np.linalg.norm(line_dir) > 1e-10:
    line_dir_n = line_dir / np.linalg.norm(line_dir)
    G_all = grid @ b_matrix
    for G in G_all:
        if np.linalg.norm(G) < 1e-10:
            continue
        denom = np.dot(line_dir_n, G)
        if abs(denom) < 1e-10:
            continue
        t = np.dot(G, G) / (2 * denom)
        for sign in [1, -1]:
            pt = line_dir_n * t * sign
            if inside_bz(pt) and in_fundamental_domain(pt):
                candidates.append(pt)

# 4. Triple intersections: 1 sym plane + 2 BZ planes
G_all = grid @ b_matrix
G_nz = [G for G in G_all if np.linalg.norm(G) > 1e-10]
for n_sym in sym_normals:
    for i_g, G1 in enumerate(G_nz):
        for j_g, G2 in enumerate(G_nz):
            if j_g <= i_g:
                continue
            A = np.array([n_sym, G1, G2])
            if abs(np.linalg.det(A)) < 1e-10:
                continue
            b_rhs = np.array([0, np.dot(G1, G1) / 2, np.dot(G2, G2) / 2])
            pt = np.linalg.solve(A, b_rhs)
            if inside_bz(pt) and in_fundamental_domain(pt):
                candidates.append(pt)

# Remove duplicates
unique = [candidates[0]]
for c in candidates[1:]:
    if not any(np.allclose(c, u, atol=1e-6) for u in unique):
        unique.append(c)

candidates_arr = np.array(unique)
candidates_frac = candidates_arr @ b_inv

# Sort by k3 then k2
order = np.lexsort((candidates_frac[:, 1], candidates_frac[:, 2]))
candidates_arr = candidates_arr[order]
candidates_frac = candidates_frac[order]

# Get HPKOT points
from lattice_kpoints import get_kpoints, _mclc1_params

conv_lattice = np.array(sp_result.get('conv_lattice', sp_result.get('primitive_lattice')))
va, vb, vc = conv_lattice
a = np.linalg.norm(va)
b_param = np.linalg.norm(vb)
c_param = np.linalg.norm(vc)
beta = np.degrees(np.arccos(np.clip(np.dot(va, vc) / (a * c_param), -1, 1)))

params = _mclc1_params(a, b_param, c_param, beta)
zeta, eta, psi, phi = params['zeta'], params['eta'], params['psi'], params['phi']
print(f"\nConventional: a={a:.6f}, b={b_param:.6f}, c={c_param:.6f}, beta={beta:.4f}")
print(f"Params: zeta={zeta:.6f}, eta={eta:.6f}, psi={psi:.6f}, phi={phi:.6f}")

kpoints_frac_table = get_kpoints('MCLC1', a, b_param, c_param, beta)

print(f"\n{'='*70}")
print(f"IRBZ vertices: {len(unique)}")
print(f"{'='*70}")
for i, (cart, frac) in enumerate(zip(candidates_arr, candidates_frac)):
    match = ""
    for label, coords in kpoints_frac_table.items():
        kc = coords[0] * b1 + coords[1] * b2 + coords[2] * b3
        if np.allclose(kc, cart, atol=1e-4):
            match = f" = {label}"
            break
    on = []
    if abs(frac[0] - frac[1]) < 1e-5:
        on.append("k1=k2")
    if abs(frac[0] + frac[1]) < 1e-5:
        on.append("k1=-k2")
    if abs(frac[2]) < 1e-5 or abs(frac[2] - 1) < 1e-5:
        on.append("k3=0")
    bnd = f" [{', '.join(on)}]" if on else ""
    print(f"  V{i:2d}: [{frac[0]:9.5f}, {frac[1]:9.5f}, {frac[2]:9.5f}]{match}{bnd}")

# Volumes
hull_irbz = ConvexHull(candidates_arr)
all_bz_verts = np.vstack(bz_faces)
bz_hull = ConvexHull(all_bz_verts)
print(f"\nIRBZ volume: {hull_irbz.volume:.8f}")
print(f"BZ volume:   {bz_hull.volume:.8f}")
print(f"BZ/4:        {bz_hull.volume/4:.8f}")
print(f"IRBZ/(BZ/4): {hull_irbz.volume/(bz_hull.volume/4):.6f}")

# HPKOT points status
print(f"\n{'='*70}")
print("HPKOT points status:")
print(f"{'='*70}")
for label, coords in kpoints_frac_table.items():
    kc = coords[0] * b1 + coords[1] * b2 + coords[2] * b3
    fd = in_fundamental_domain(kc)
    bz = inside_bz(kc)
    is_vertex = any(np.allclose(kc, c, atol=1e-4) for c in candidates_arr)
    if not bz:
        status = "OUTSIDE_BZ!"
    elif is_vertex:
        status = "VERTEX"
    elif fd:
        status = "interior/edge"
    else:
        status = "outside_FD"
    print(f"  {label:5s}: [{coords[0]:8.4f}, {coords[1]:8.4f}, {coords[2]:8.4f}]  {status}")

# Identify unlabeled vertices
print(f"\n{'='*70}")
print("UNLABELED IRBZ vertices (need formulas):")
print(f"{'='*70}")
for j, (cand, cf) in enumerate(zip(candidates_arr, candidates_frac)):
    matched = False
    for label, coords in kpoints_frac_table.items():
        kc = coords[0] * b1 + coords[1] * b2 + coords[2] * b3
        if np.allclose(kc, cand, atol=1e-4):
            matched = True
            break
    if not matched:
        # Try to express in terms of parameters
        print(f"  V{j}: [{cf[0]:9.5f}, {cf[1]:9.5f}, {cf[2]:9.5f}]")
        # Check relations
        print(f"       k1+k2 = {cf[0]+cf[1]:.5f}")
        print(f"       k1-k2 = {cf[0]-cf[1]:.5f}")
        print(f"       k3    = {cf[2]:.5f}")
        # Check against known params
        for pname, pval in [('zeta', zeta), ('eta', eta), ('psi', psi), ('phi', phi),
                            ('1-zeta', 1-zeta), ('1-eta', 1-eta), ('1-psi', 1-psi), ('1-phi', 1-phi),
                            ('-1+2*psi', -1+2*psi), ('-1+2*phi', -1+2*phi),
                            ('-1+zeta+psi', -1+zeta+psi), ('-1-zeta+psi', -1-zeta+psi),
                            ('zeta+psi', zeta+psi), ('psi-zeta', psi-zeta)]:
            for coord_name, coord_val in [('k1', cf[0]), ('k2', cf[1]), ('k3', cf[2]),
                                           ('k1+k2', cf[0]+cf[1]), ('k1-k2', cf[0]-cf[1])]:
                if abs(coord_val - pval) < 1e-4:
                    print(f"       {coord_name} ~ {pname} = {pval:.5f}")
                if abs(coord_val + pval) < 1e-4:
                    print(f"       {coord_name} ~ -{pname} = {-pval:.5f}")
                if abs(coord_val - pval/2) < 1e-4:
                    print(f"       {coord_name} ~ {pname}/2 = {pval/2:.5f}")
