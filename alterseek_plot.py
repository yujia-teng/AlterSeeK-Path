"""Command-line dispatcher for the VASP, QE, and ABINIT band plotters."""
from pathlib import Path
import sys


_CONFIG_FILENAMES = {
    "vasp": "alterseek_plot_vasp.toml",
    "qe": "alterseek_plot_qe.toml",
    "abinit": "alterseek_plot_abinit.toml",
}


def _print_help(file=None):
    stream = sys.stdout if file is None else file
    print("usage: alterseek-plot {vasp,qe,abinit} [options]", file=stream)
    print(file=stream)
    print("Plot spin-resolved bands for VASP, Quantum ESPRESSO, or ABINIT.", file=stream)
    print("If exactly one AlterSeeK-Path plot config is present, the code", file=stream)
    print("argument may be omitted.", file=stream)
    print(file=stream)
    print("codes:", file=stream)
    print("  vasp      VASPKIT reformatted band files", file=stream)
    print("  qe        Quantum ESPRESSO bands.x .gnu files", file=stream)
    print("  abinit    ABINIT plain-text _EIG output", file=stream)


def _detected_codes(directory=Path(".")):
    """Return codes identified by AlterSeeK-Path plot config filenames."""
    return [
        code
        for code, filename in _CONFIG_FILENAMES.items()
        if (directory / filename).is_file()
    ]


def _run_plotter(code, argv):
    if code == "vasp":
        from plotting.plot_alterband import main as plot_main
    elif code == "qe":
        from plotting.plot_alterband_qe import main as plot_main
    else:
        from plotting.plot_alterband_abinit import main as plot_main
    return plot_main(argv, prog=f"alterseek-plot {code}")


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] in {"-h", "--help"}:
        _print_help()
        return 0

    if args and args[0].lower() in _CONFIG_FILENAMES:
        code = args.pop(0).lower()
        return _run_plotter(code, args)

    if args and not args[0].startswith("-"):
        print(f"[Error] Unknown code: {args[0]}", file=sys.stderr)
        _print_help(file=sys.stderr)
        return 2

    detected = _detected_codes()
    if len(detected) == 1:
        return _run_plotter(detected[0], args)
    if len(detected) > 1:
        names = ", ".join(detected)
        print(
            f"[Error] Multiple plot configurations found ({names}); "
            "specify a code.",
            file=sys.stderr,
        )
        return 2

    if args:
        print(
            "[Error] No plot configuration found; specify vasp, qe, or abinit.",
            file=sys.stderr,
        )
        return 2
    _print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
