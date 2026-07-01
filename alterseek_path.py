# Script for generating general k-path for band structure calculation
# Author: Yujia Teng, Rutgers University
# 06/14/2025 - Creation of the script
# 01/12/2026 - Automatically reads the spin-flip operation based on outputs of find_sf_operations.py
# 03/2026    - Integrated find_sf_operations and compute_centroid_hybrid into one workflow
import numpy as np
import os
import shutil
import sys
from typing import List, Optional

# Try to import spin-flip operations finder
try:
    from find_sf_operations import run as find_sf_run
    from find_sf_operations import parse_cartesian_spin_axis, parse_magmoms
    FIND_SF_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] find_sf_operations unavailable ({_exc}); "
          "Step 0 spin-operation detection is disabled.")
    FIND_SF_AVAILABLE = False
    parse_cartesian_spin_axis = None
    parse_magmoms = None

try:
    from findspingroup import (
        find_spin_group_acc_primitive,
        find_spin_group_acc_primitive_from_data,
    )
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] findspingroup unavailable ({_exc}); --ssg-setting is disabled. "
          "Install with: pip install \"findspingroup>=0.15.6,<0.16\"")
    find_spin_group_acc_primitive = None
    find_spin_group_acc_primitive_from_data = None
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = False

# Try to import IBZ centroid calculator
try:
    from compute_centroid_hybrid import (run as compute_centroid,
                                         no_altermagnetism_reason,
                                         is_valid_2d_spin_flip,
                                         is_trivial_2d_spin_flip,
                                         plot_spin_flip_figure,
                                         plot_spin_bz_figure,
                                         plot_spin_bz_top_view_figure,
                                         describe_spinflip_op)
    from plot_2d_figures import plot_2d_figures
    import matplotlib.pyplot as plt
    CENTROID_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] compute_centroid_hybrid/matplotlib unavailable ({_exc}); "
          "centroid and figure generation are disabled.")
    CENTROID_AVAILABLE = False
    no_altermagnetism_reason = None
    plot_spin_flip_figure = None
    plot_spin_bz_figure = None
    plot_spin_bz_top_view_figure = None
    describe_spinflip_op = None
    plot_2d_figures = None
    plt = None


# ==========================================================================
# Restructuring phase 4: structure/MCIF I/O + band-plot config moved to io_vasp.py.
# Re-exported so existing `from alterseek_path import <fn>` keeps working.
# ==========================================================================
from io_vasp import (
    _group_poscar_sites,
    _write_poscar,
    _read_grouped_poscar,
    _write_poscar_from_sites,
    _write_without_species,
    _reciprocal_from_poscar,
    _dedupe_frac_positions,
    _min_periodic_cart_distance,
    _lattice_lengths_angles,
    _write_magnetic_mcif,
    _load_magnetic_input_data,
    write_bandplot_lattice_config,
)


# ==========================================================================
# Restructuring phase 4: SSG-setting (He-marker) workflow moved to ssg_setting.py.
# Re-exported so existing `from alterseek_path import <fn>` keeps working.
# ==========================================================================
from ssg_setting import (
    _marker_orbit,
    _build_marker_helper,
    _spin_axis_from_moments,
    _operation_class_indices,
    _collect_point_ops_from_payload,
    _write_operation_file,
    _write_magnetic_setting_operation_files,
    _magnetic_primitive_ssg_operations,
    prepare_magnetic_setting_files,
    finalize_magnetic_setting_outputs,
)


# ==========================================================================
# Restructuring phase 4: KPointsModifier moved to kpoints.py.
# Re-exported so existing `from alterseek_path import KPointsModifier` keeps working.
# ==========================================================================
from kpoints import KPointsModifier


def main():
    import argparse

    argv = sys.argv[1:]
    # Forward the bandplot subcommand untouched so plot_alterband's own
    # parser handles its options (including --help).
    if argv and argv[0].lower() in {"bandplot", "plot-band", "plot"}:
        from plot_alterband import main as plot_alterband_main
        plot_alterband_main(argv[1:])
        return

    parser = argparse.ArgumentParser(
        prog="alterseek-path",
        description=(
            "Generate an altermagnetic KPOINTS path interactively. "
            "Subcommand: 'bandplot' plots spin-resolved bands from KLABELS "
            "and reformatted band data (run: alterseek-path bandplot --help)."
        ),
    )
    parser.add_argument(
        "--ssg-setting",
        action="store_true",
        help="Experimental: generate Figure 1/KPOINTS from FindSpinGroup SSG setting.",
    )
    parser.add_argument(
        "--output",
        choices=["verbose"],
        help="verbose: keep SSG-setting intermediate/helper structures for debugging.",
    )
    parser.add_argument(
        "--2d",
        dest="mode_2d",
        action="store_true",
        help="2D/slab mode: restrict the k-path and IBZ centroid to the physical "
             "in-plane (vacuum k=0) reciprocal plane.",
    )
    parser.add_argument(
        "--vacuum-axis",
        choices=["a", "b", "c"],
        default="c",
        help="Input-cell vacuum axis for the 2D-slab sanity check (default: c). "
             "The slicing axis is auto-detected in the standardized frame "
             "regardless of this flag.",
    )
    args = parser.parse_args(argv)

    modifier = KPointsModifier(
        magnetic_setting=args.ssg_setting,
        output_verbose=args.output == "verbose",
        mode_2d=args.mode_2d,
        input_vacuum_axis={"a": 0, "b": 1, "c": 2}[args.vacuum_axis],
    )
    modifier.interactive_modify()


if __name__ == "__main__":
    main()

