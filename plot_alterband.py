#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from VASPKIT reformatted data."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10 fallback
    tomllib = None

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = [
    "Times New Roman",
    "STIX Two Text",
    "Nimbus Roman",
    "Liberation Serif",
    "DejaVu Serif",
]
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator


DEFAULT_ELIM = (-2.0, 2.0)
DEFAULT_GAP_FRAC = 0.004
DEFAULT_FIG_SIZE = (10.0, 5.0)
GREY_COLOR = "0.65"
BAND_LW = 0.7
BAND_UP_COLOR = "black"
BAND_DOWN_COLOR = "red"
VLINE_LW = 0.8
VLINE_COLOR = "black"
FERMI_LW = 1.2
FERMI_COLOR = "0"
FONT_SIZE = 14

GREEK_LABELS = {
    "GAMMA": r"$\Gamma$",
    "DELTA": r"$\Delta$",
    "LAMBDA": r"$\Lambda$",
    "SIGMA": r"$\Sigma$",
    "\u0393": r"$\Gamma$",
    "\u0394": r"$\Delta$",
    "\u039b": r"$\Lambda$",
    "\u03a3": r"$\Sigma$",
    "\u8795": r"$\Gamma$",
}

def _read_klabels(path: Path) -> tuple[list[str], list[float]]:
    labels: list[str] = []
    positions: list[float] = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                positions.append(float(parts[1]))
            except ValueError:
                continue
            labels.append(parts[0])

    if not labels:
        raise ValueError(f"No KLABELS entries found in {path}")
    return labels, positions


def _parse_simple_toml_value(value: str) -> Any:
    value = value.strip()
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Unsupported TOML value: {value}") from exc


def _read_plot_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    if tomllib is not None:
        with path.open("rb") as f:
            data = tomllib.load(f)
        if not isinstance(data, dict):
            raise ValueError(f"Config file must contain key-value settings: {path}")
        return data

    config: dict[str, Any] = {}
    with path.open() as f:
        for raw_line in f:
            line = raw_line.split("#", 1)[0].strip()
            if not line or line.startswith("["):
                continue
            if "=" not in line:
                raise ValueError(f"Invalid config line in {path}: {raw_line.rstrip()}")
            key, value = line.split("=", 1)
            config[key.strip()] = _parse_simple_toml_value(value)
    return config


def _format_tick_label(label: str) -> str:
    """Return labels with mathtext only where real Greek/subscripts are needed."""
    if "|" in label:
        return "|".join(_format_tick_label(part) for part in label.split("|"))

    prime_count = 0
    while label.endswith("'"):
        prime_count += 1
        label = label[:-1]

    if "_" in label:
        base, subscript = label.split("_", 1)
    else:
        base, subscript = label, None

    greek = GREEK_LABELS.get(base.upper(), GREEK_LABELS.get(base))
    if subscript is not None:
        sub_body = subscript if subscript.isdigit() else rf"\mathrm{{{subscript}}}"
        if greek is not None:
            body = greek[1:-1]
        else:
            body = rf"\mathrm{{{base}}}"
        return rf"${body}_{{{sub_body}}}$" + "'" * prime_count

    if greek is not None:
        return greek + "'" * prime_count

    if prime_count:
        return base + "'" * prime_count
    return label


def plot_alterband(
    *,
    klabels: str | Path = "KLABELS",
    band_up: str | Path = "REFORMATTED_BAND_UP.dat",
    band_down: str | Path = "REFORMATTED_BAND_DW.dat",
    output: str | Path = "alterband.png",
    elim: tuple[float, float] = DEFAULT_ELIM,
    fig_size: tuple[float, float] = DEFAULT_FIG_SIZE,
    gap_frac: float = DEFAULT_GAP_FRAC,
    rotate_xtick_labels: bool = False,
    xtick_rotation: float = 45.0,
) -> Path:
    """Create the spin-resolved band plot and return the output path."""
    klabels_path = Path(klabels)
    band_up_path = Path(band_up)
    band_down_path = Path(band_down)
    output_path = Path(output)

    labels, positions = _read_klabels(klabels_path)
    x_total = positions[-1] - positions[0]
    gap_half = x_total * gap_frac

    tick_lab = [_format_tick_label(label) for label in labels]
    alt_gap = {"k|k'", "k'|k"}
    non_gap_pos = [p for label, p in zip(labels, positions) if label not in alt_gap]
    gap_pos = [p for label, p in zip(labels, positions) if label in alt_gap]

    up = np.loadtxt(band_up_path, skiprows=1)
    dw = np.loadtxt(band_down_path, skiprows=1)
    kpath = up[:, 0]
    bands_up = up[:, 1:]
    bands_dw = dw[:, 1:]

    in_window = (
        ((bands_up.max(axis=0) >= elim[0]) & (bands_up.min(axis=0) <= elim[1]))
        | ((bands_dw.max(axis=0) >= elim[0]) & (bands_dw.min(axis=0) <= elim[1]))
    )
    bands_up = bands_up[:, in_window]
    bands_dw = bands_dw[:, in_window]

    fig, ax = plt.subplots(figsize=fig_size)

    for i in range(len(positions) - 1):
        if "k" not in labels[i] and "k" not in labels[i + 1]:
            ax.axvspan(
                positions[i], positions[i + 1], color=GREY_COLOR, lw=0, zorder=0
            )

    for ib in range(bands_up.shape[1]):
        ax.plot(kpath, bands_dw[:, ib], color=BAND_DOWN_COLOR, lw=BAND_LW, zorder=2)
    for ib in range(bands_up.shape[1]):
        ax.plot(kpath, bands_up[:, ib], color=BAND_UP_COLOR, lw=BAND_LW, zorder=3)

    for pos in gap_pos:
        ax.axvspan(pos - gap_half, pos + gap_half, color="white", zorder=4, lw=0)

    for pos in non_gap_pos:
        ax.axvline(x=pos, color=VLINE_COLOR, lw=VLINE_LW, zorder=5)

    for pos in gap_pos:
        ax.axvline(x=pos - gap_half, color=VLINE_COLOR, lw=VLINE_LW, zorder=5)
        ax.axvline(x=pos + gap_half, color=VLINE_COLOR, lw=VLINE_LW, zorder=5)

    ax.axhline(y=0, color=FERMI_COLOR, lw=FERMI_LW, ls="--", zorder=1)
    ax.set_xticks(positions)
    xtick_label_kwargs: dict[str, Any] = {"fontsize": FONT_SIZE}
    if rotate_xtick_labels:
        xtick_label_kwargs.update(
            {
                "rotation": xtick_rotation,
                "ha": "right",
                "rotation_mode": "anchor",
            }
        )
    ax.set_xticklabels(tick_lab, **xtick_label_kwargs)
    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(elim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.set_ylabel(r"E - E$_\mathrm{F}$ (eV)", fontsize=FONT_SIZE + 1)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot spin-resolved AlterSeeK band output from VASPKIT files."
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Optional TOML config file. Defaults to alterband.toml if present.",
    )
    parser.add_argument("--klabels", default=None, help="KLABELS file path.")
    parser.add_argument(
        "--up", default=None, help="Spin-up band data file."
    )
    parser.add_argument(
        "--down", default=None, help="Spin-down band data file."
    )
    parser.add_argument("-o", "--output", default=None, help="Output PNG.")
    parser.add_argument(
        "--emin", type=float, default=None, help="Minimum plotted energy."
    )
    parser.add_argument(
        "--emax", type=float, default=None, help="Maximum plotted energy."
    )
    parser.add_argument(
        "--fig-width", type=float, default=None, help="Figure width in inches."
    )
    parser.add_argument(
        "--fig-height", type=float, default=None, help="Figure height in inches."
    )
    parser.add_argument(
        "--gap-frac",
        type=float,
        default=None,
        help="Half-width of k|k' gap as a fraction of the total k-path.",
    )
    parser.add_argument(
        "--rotate-xtick-labels",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Rotate all x-axis tick labels. Can also be set in TOML.",
    )
    parser.add_argument(
        "--xtick-rotation",
        type=float,
        default=None,
        help="X-axis tick label rotation angle in degrees when rotation is enabled.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    config_path = Path(args.config) if args.config else Path("alterband.toml")
    config = _read_plot_config(config_path) if args.config or config_path.exists() else {}

    def option(name: str, default: Any) -> Any:
        arg_value = getattr(args, name)
        if arg_value is not None:
            return arg_value
        return config.get(name, default)

    emin = float(option("emin", DEFAULT_ELIM[0]))
    emax = float(option("emax", DEFAULT_ELIM[1]))
    fig_width = float(option("fig_width", DEFAULT_FIG_SIZE[0]))
    fig_height = float(option("fig_height", DEFAULT_FIG_SIZE[1]))
    rotate_xtick_labels = bool(option("rotate_xtick_labels", False))
    xtick_rotation = float(option("xtick_rotation", 45.0))
    output = plot_alterband(
        klabels=option("klabels", "KLABELS"),
        band_up=option("up", config.get("band_up", "REFORMATTED_BAND_UP.dat")),
        band_down=option("down", config.get("band_down", "REFORMATTED_BAND_DW.dat")),
        output=option("output", "alterband.png"),
        elim=(emin, emax),
        fig_size=(fig_width, fig_height),
        gap_frac=float(option("gap_frac", DEFAULT_GAP_FRAC)),
        rotate_xtick_labels=rotate_xtick_labels,
        xtick_rotation=xtick_rotation,
    )
    print(f"Band plot written to: {output}")


if __name__ == "__main__":
    main()
