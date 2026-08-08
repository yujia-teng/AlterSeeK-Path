# Script for generating general k-path for band structure calculation
# Author: Yujia Teng, Rutgers University
# 06/14/2025 - Creation of the script
# 01/12/2026 - Automatically reads the spin-flip operation based on outputs of find_sf_operations.py
# 03/2026    - Integrated find_sf_operations and compute_centroid_hybrid into one workflow
# 07/2026    - Restructured into focused modules; this file is now just the CLI entry point.
#              The workflow lives in kpoints.KPointsModifier (see io / ssg_setting /
#              compute_centroid_hybrid / plotting_* / symmetry / geometry).
import sys

from alterseek.kpoints import KPointsModifier


def main():
    import argparse

    argv = sys.argv[1:]
    # Forward the bandplot subcommand untouched so plot_alterband's own
    # parser handles its options (including --help).
    if argv and argv[0].lower() in {"bandplot", "plot-band", "plot"}:
        from plot_alterband import main as plot_alterband_main
        plot_alterband_main(argv[1:])
        return 0

    parser = argparse.ArgumentParser(
        prog="alterseek-path",
        description=(
            "Generate an altermagnetic KPOINTS path interactively. "
            "Subcommand: 'bandplot' plots spin-resolved bands from KLABELS "
            "and reformatted band data (run: alterseek-path bandplot --help)."
        ),
    )
    parser.add_argument(
        "--parent-setting",
        action="store_true",
        help="Build Figure 1/KPOINTS in the nonmagnetic parent cell instead of "
             "the magnetic primitive cell. Useful for comparing several magnetic "
             "orders against one fixed reference path, but that path does not "
             "respect the symmetry of the magnetic state.",
    )
    parser.add_argument(
        "--output",
        choices=["verbose"],
        help="verbose: keep intermediate/helper structures for debugging.",
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
        magnetic_setting=not args.parent_setting,
        output_verbose=args.output == "verbose",
        mode_2d=args.mode_2d,
        input_vacuum_axis={"a": 0, "b": 1, "c": 2}[args.vacuum_axis],
    )
    try:
        success = modifier.interactive_modify()
    except Exception as exc:
        # This is the command-line workflow's final safety boundary. Required
        # calculations raise to here; expected input failures and optional
        # output warnings remain the responsibility of their own subsystems.
        # Never continue after an unexpected failure or reinterpret it as a
        # request for different scientific input.
        print(f"[Error] AlterSeeK-Path failed: {exc}", file=sys.stderr)
        return 1
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
