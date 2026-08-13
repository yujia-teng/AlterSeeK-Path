"""Command-line entry point for AlterSeeK-Path."""
import sys

from alterseek.kpoints import KPointsModifier


def main():
    import argparse

    argv = sys.argv[1:]
    # Forward the bandplot subcommand untouched so plot_alterband parses its own options, including --help.
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
        output_verbose=args.output == "verbose",
        mode_2d=args.mode_2d,
        input_vacuum_axis={"a": 0, "b": 1, "c": 2}[args.vacuum_axis],
    )
    try:
        success = modifier.interactive_modify()
    except Exception as exc:
        # This is the command-line workflow's final failure boundary.
        # Required calculations raise here, while expected input failures and optional-output warnings remain with their own subsystems.
        # Never continue after an unexpected failure or reinterpret it as a request for different scientific input.
        print(f"[Error] AlterSeeK-Path failed: {exc}", file=sys.stderr)
        return 1
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
