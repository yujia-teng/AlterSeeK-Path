"""
Batch Figure-2 generator — tests every spin-flip operation without user interaction.

Usage (run from the structure folder that contains POSCAR):

    python "G:\\...\\worktrees\\op-type-visualization\\tools\\batch_test_fig2.py" [moments]

  moments  (optional) space-separated VASP MAGMOM string, default "1 -1"

Output: fig2_op01.png … fig2_opNN.png  in the current working directory,
        one image per spin-flip operation (including inversion-extended partners).
"""

import sys
import os
import argparse

import matplotlib
matplotlib.use('Agg')          # no GUI window needed
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Locate the worktree so relative imports work regardless of cwd
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)   # worktree root
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import find_sf_operations
import compute_centroid_hybrid
from compute_centroid_hybrid import (
    plot_spin_flip_figure,
    _classify_spinflip_op,
)

# ---------------------------------------------------------------------------
# Helpers (mirrors the logic in alterseek_path.py)
# ---------------------------------------------------------------------------

def _load_ops(filename):
    """Parse spin_flip_operations.txt — left side only for '|' format."""
    matrices, unique, cur = [], [], []
    if not os.path.exists(filename):
        return matrices
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('Operation'):
                continue
            if '|' in line:
                line = line.split('|', 1)[0].strip()
            parts = line.split()
            if len(parts) == 3:
                cur.append([float(x) for x in parts])
            if len(cur) == 3:
                mat = np.array(cur)
                if not any(np.allclose(mat, ex, atol=1e-8) for ex in unique):
                    unique.append(mat)
                    matrices.append(mat)
                cur = []
    return matrices


def _inversion_extended(ops):
    expanded = list(ops)
    for op in ops:
        neg = -np.array(op, dtype=float)
        if not any(np.allclose(neg, ex, atol=1e-8) for ex in expanded):
            expanded.append(neg)
    return expanded


def _to_prim(R_input, b_input, b_prim):
    """Convert R from input-cell fractional to primitive-cell fractional."""
    R_cart = b_input.T @ np.linalg.inv(R_input).T @ np.linalg.inv(b_input.T)
    R_prim_inv_T = np.linalg.inv(b_prim.T) @ R_cart @ b_prim.T
    return np.linalg.inv(R_prim_inv_T.T)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Batch Figure-2 test for all flip ops.')
    parser.add_argument('moments', nargs='?', default='1 -1',
                        help='Magnetic moments string (default: "1 -1")')
    parser.add_argument('--poscar', default='POSCAR',
                        help='Structure file (default: POSCAR)')
    parser.add_argument('--spin-axis', default=None,
                        help='Cartesian spin axis, e.g. "0 0 1" (default: 0 0 1)')
    args = parser.parse_args()

    poscar = args.poscar
    moments = args.moments
    spin_axis = args.spin_axis

    if not os.path.exists(poscar):
        print(f"[Error] {poscar} not found in {os.getcwd()}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 0: find spin-flip operations
    # ------------------------------------------------------------------
    print('=== Step 0: finding spin-flip operations ===')
    sf_result = find_sf_operations.run(poscar, moments,
                                       verbose=True, spin_axis_cart=spin_axis)
    if not sf_result:
        print('[Error] find_sf_operations failed.')
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 1: compute BZ / centroid data (no plot)
    # ------------------------------------------------------------------
    print('\n=== Step 1: computing BZ / centroid ===')
    centroid_result = compute_centroid_hybrid.run(poscar, '.', show_plot=False, verbose=False)
    if centroid_result is None:
        print('[Error] compute_centroid_hybrid failed.')
        sys.exit(1)

    b_matrix  = np.array(centroid_result['b_matrix'])
    b_input   = np.array(centroid_result.get('b_matrix_input',
                         centroid_result.get('b_matrix_conv', b_matrix)))
    rot_mat   = np.array(centroid_result.get('seekpath_rotation_matrix', np.eye(3)))
    b_prim    = b_matrix @ rot_mat

    # ------------------------------------------------------------------
    # Step 2: load and convert operations
    # ------------------------------------------------------------------
    raw_ops = _load_ops('spin_flip_operations.txt')
    if not raw_ops:
        print('[Error] spin_flip_operations.txt is empty or missing.')
        sys.exit(1)

    flip_ops = _inversion_extended(raw_ops)
    print(f'\n{len(flip_ops)} spin-flip operations (including inversion-extended).')

    # Convert each to primitive-cell fractional basis
    prim_ops = [_to_prim(op, b_input, b_prim) for op in flip_ops]

    # ------------------------------------------------------------------
    # Step 3: generate Figure 2 for every operation
    # ------------------------------------------------------------------
    print('\n=== Step 2: generating figures ===')

    # R_cart computation mirrors plot_spin_flip_figure internals
    b_T     = b_matrix.T
    b_T_inv = np.linalg.inv(b_T)

    for i, R_prim in enumerate(prim_ops, start=1):
        R_cart = b_T @ np.linalg.inv(R_prim).T @ b_T_inv
        op_info = _classify_spinflip_op(R_cart)
        op_label = op_info['type']
        if op_info['order']:
            op_label += str(op_info['order'])

        out_path = f'fig2_op{i:02d}_{op_label}.png'
        print(f'  Op {i:2d}: {op_label:18s} -> {out_path}')

        plot_spin_flip_figure(
            b_matrix        = b_matrix,
            bz_loops        = centroid_result['bz_loops'],
            bz_center       = centroid_result['bz_center'],
            bz_span         = centroid_result['bz_span'],
            kpoints_data    = [],          # no path segments needed for visual test
            ibz_kpoints_frac= centroid_result.get('ibz_kpoints_frac', {}),
            hull_pts        = centroid_result.get('hull_pts'),
            hull_simplices  = centroid_result.get('hull_simplices'),
            centroid_frac   = centroid_result['centroid_frac'],
            R               = R_prim,
            output_path     = out_path,
            elev            = centroid_result.get('elev', 14),
            azim            = centroid_result.get('azim', 20),
            show_plot       = False,
            path_sequence   = None,
            unique_ops      = centroid_result.get('unique_ops'),
        )

    print(f'\nDone — {len(prim_ops)} figures saved in {os.getcwd()}')


if __name__ == '__main__':
    main()
