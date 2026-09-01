"""Command-line entry point for AlterSeeK-Path."""
import os
import sys

from alterseek.kpoints import KPathBuilder, OUTPUT_DIR
from alterseek.mode2d.kpoints import KPathBuilder2D
from alterseek.run_log import RUN_LOG_FILENAME, run_log


def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="alterseek-path",
        description="Generate a general-k path.",
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
        default=None,
        help="Input-cell vacuum axis for the 2D-slab sanity check (default: c, "
             "or vacuum_axis in alterseek_input.toml). The slicing axis is "
             "auto-detected in the standardized frame regardless of this flag.",
    )
    args = parser.parse_args(sys.argv[1:])

    vacuum_axis = (
        None if args.vacuum_axis is None
        else {"a": 0, "b": 1, "c": 2}[args.vacuum_axis]
    )
    if args.mode_2d:
        builder = KPathBuilder2D(input_vacuum_axis=vacuum_axis)
    else:
        builder = KPathBuilder()
    with run_log(os.path.join(OUTPUT_DIR, RUN_LOG_FILENAME)) as log:
        try:
            success = builder.interactive_build()
        except Exception as exc:
            # This is the command-line workflow's final failure boundary.
            # Required calculations raise here, while expected input failures and optional-output warnings remain with their own subsystems.
            # Never continue after an unexpected failure or reinterpret it as a request for different scientific input.
            print(f"[Error] AlterSeeK-Path failed: {exc}", file=sys.stderr)
            return 1
        if success:
            log.mark_success()
        return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
