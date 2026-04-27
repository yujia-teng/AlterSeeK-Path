#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from VASPKIT reformatted data."""

from __future__ import annotations

import argparse
from pathlib import Path

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


DEFAULT_ELIM = (-2.0, 2.0)
DEFAULT_GAP_FRAC = 0.004
GREY_COLOR = "0.65"
BAND_LW = 0.7
BAND_UP_COLOR = "black"
BAND_DOWN_COLOR = "red"
VLINE_LW = 0.8
VLINE_COLOR = "black"
FERMI_LW = 1.2
FERMI_COLOR = "0"
FONT_SIZE = 14


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


def plot_alterband(
    *,
    klabels: str | Path = "KLABELS",
    band_up: str | Path = "REFORMATTED_BAND_UP.dat",
    band_down: str | Path = "REFORMATTED_BAND_DW.dat",
    output: str | Path = "alterband.png",
    elim: tuple[float, float] = DEFAULT_ELIM,
    gap_frac: float = DEFAULT_GAP_FRAC,
) -> Path:
    """Create the spin-resolved band plot and return the output path."""
    klabels_path = Path(klabels)
    band_up_path = Path(band_up)
    band_down_path = Path(band_down)
    output_path = Path(output)

    labels, positions = _read_klabels(klabels_path)
    x_total = positions[-1] - positions[0]
    gap_half = x_total * gap_frac

    tick_lab = [r"$\Gamma$" if label == "GAMMA" else label for label in labels]
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

    fig, ax = plt.subplots(figsize=(10, 5))

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
    ax.set_xticklabels(tick_lab, fontsize=FONT_SIZE)
    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(elim)
    ax.set_yticks(range(int(elim[0]), int(elim[1]) + 1))
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
    parser.add_argument("--klabels", default="KLABELS", help="KLABELS file path.")
    parser.add_argument(
        "--up", default="REFORMATTED_BAND_UP.dat", help="Spin-up band data file."
    )
    parser.add_argument(
        "--down", default="REFORMATTED_BAND_DW.dat", help="Spin-down band data file."
    )
    parser.add_argument("-o", "--output", default="alterband.png", help="Output PNG.")
    parser.add_argument(
        "--emin", type=float, default=DEFAULT_ELIM[0], help="Minimum plotted energy."
    )
    parser.add_argument(
        "--emax", type=float, default=DEFAULT_ELIM[1], help="Maximum plotted energy."
    )
    parser.add_argument(
        "--gap-frac",
        type=float,
        default=DEFAULT_GAP_FRAC,
        help="Half-width of k|k' gap as a fraction of the total k-path.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    output = plot_alterband(
        klabels=args.klabels,
        band_up=args.up,
        band_down=args.down,
        output=args.output,
        elim=(args.emin, args.emax),
        gap_frac=args.gap_frac,
    )
    print(f"Band plot written to: {output}")


if __name__ == "__main__":
    main()
