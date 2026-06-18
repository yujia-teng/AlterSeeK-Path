#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from Quantum ESPRESSO bands.x .gnu files."""

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
mpl.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator


DEFAULT_ELIM = (-2.0, 2.0)
DEFAULT_FIG_SIZE = (12.0, 5.0)
DEFAULT_GAP_FRAC = 0.003
GREY_COLOR = "0.65"
BAND_LW = 0.7
BAND_UP_COLOR = "black"
BAND_DOWN_COLOR = "red"
VLINE_COLOR = "black"
VLINE_LW = 0.8
FERMI_LW = 0.5
FERMI_COLOR = "0.5"
FONT_SIZE = 14
GAP_LABELS = {"k|k'", "k'|k"}
HELPER_LABELS = {"k", "k'", *GAP_LABELS}

GREEK_LABELS = {
    "GAMMA": r"$\Gamma$",
    "Γ": r"$\Gamma$",
    "螕": r"$\Gamma$",
    "DELTA": r"$\Delta$",
    "Δ": r"$\Delta$",
    "LAMBDA": r"$\Lambda$",
    "Λ": r"$\Lambda$",
    "SIGMA": r"$\Sigma$",
    "Σ": r"$\Sigma$",
}


def _parse_simple_toml_value(value: str) -> Any:
    value = value.strip()
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if value.lower() in {"none", "null"}:
        return None
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


def _read_gnu_bands(path: Path, fermi_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """Parse QE bands.x .gnu file (band-major, blank-line separated).

    Returns (kpath, bands) where bands.shape = (nk, nbands), energies shifted by fermi_ev.
    """
    blocks: list[np.ndarray] = []
    current: list[list[float]] = []
    with path.open() as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                if current:
                    blocks.append(np.array(current, dtype=float))
                    current = []
            else:
                k, e = map(float, stripped.split())
                current.append([k, e])
    if current:
        blocks.append(np.array(current, dtype=float))
    if not blocks:
        raise ValueError(f"No band data found in {path}")
    kpath = blocks[0][:, 0]
    bands = np.column_stack([b[:, 1] - fermi_ev for b in blocks])
    return kpath, bands


def _parse_kpoints_qe(path: Path) -> list[tuple[str, int, int]]:
    """Parse KPOINTS_modified_qe written by write_kpoints_file_qe.

    Returns list of (label, ninterp, k_index) for each waypoint.
    k_index is the 0-based index into the .gnu k-point array.
    """
    with path.open(encoding="utf-8-sig") as f:
        all_lines = f.readlines()

    start = -1
    for i, line in enumerate(all_lines):
        if line.strip().upper().startswith("K_POINTS"):
            start = i
            break
    if start == -1:
        raise ValueError(f"K_POINTS card not found in {path}")

    n_waypoints = int(all_lines[start + 1].strip().split()[0])

    waypoints_raw: list[tuple[str, int]] = []
    for line in all_lines[start + 2:]:
        if len(waypoints_raw) >= n_waypoints:
            break
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        # format: x y z  ni  ! LABEL
        parts = stripped.split("!", 1)
        ni = int(parts[0].split()[3])
        label = parts[1].strip() if len(parts) > 1 else ""
        waypoints_raw.append((label, ni))

    # cumulative k-indices: waypoint i starts at sum of preceding ninterps
    result: list[tuple[str, int, int]] = []
    cum = 0
    for label, ni in waypoints_raw:
        result.append((label, ni, cum))
        cum += ni
    return result


def _format_tick_label(label: str) -> str:
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
        body = greek[1:-1] if greek is not None else rf"\mathrm{{{base}}}"
        return rf"${body}_{{{sub_body}}}$" + "'" * prime_count
    if greek is not None:
        return greek + "'" * prime_count
    if prime_count:
        return base + "'" * prime_count
    return label


def _build_tick_data(
    waypoints: list[tuple[str, int, int]], kpath: np.ndarray
) -> tuple[list[str], list[float]]:
    """Convert waypoints to tick labels and x-positions.

    Adjacent k / k' pairs are merged into a single k|k' tick at their midpoint.
    """
    labels: list[str] = []
    positions: list[float] = []
    i = 0
    while i < len(waypoints):
        label, _ni, k_idx = waypoints[i]
        if (
            label == "k"
            and i + 1 < len(waypoints)
            and waypoints[i + 1][0] == "k'"
        ):
            x_k = kpath[waypoints[i][2]]
            x_kp = kpath[waypoints[i + 1][2]]
            labels.append("k|k'")
            positions.append((x_k + x_kp) / 2.0)
            i += 2
            continue
        if k_idx >= len(kpath):
            raise ValueError(
                f"Waypoint '{label}' has k-index {k_idx} but .gnu file has only "
                f"{len(kpath)} k-points — KPOINTS_modified_qe and .gnu file mismatch"
            )
        labels.append(label)
        positions.append(float(kpath[k_idx]))
        i += 1
    return labels, positions


def plot_alterband_qe(
    *,
    band_up: str | Path = "band_up.gnu",
    band_down: str | Path = "band_down.gnu",
    kpoints_qe: str | Path = "KPOINTS_modified_qe",
    output: str | Path = "alterband_qe.png",
    fermi_ev: float = 0.0,
    elim: tuple[float, float] = DEFAULT_ELIM,
    fig_size: tuple[float, float] = DEFAULT_FIG_SIZE,
) -> Path:
    """Create the QE spin-resolved band plot and return the output path."""
    band_up_path = Path(band_up)
    band_down_path = Path(band_down)
    kpoints_path = Path(kpoints_qe)
    output_path = Path(output)

    kpath_up, bands_up = _read_gnu_bands(band_up_path, fermi_ev)
    kpath_dw, bands_dw = _read_gnu_bands(band_down_path, fermi_ev)
    if len(kpath_up) != len(kpath_dw) or not np.allclose(kpath_up, kpath_dw, atol=1e-5):
        raise ValueError(
            "Spin-up and spin-down .gnu files have different k-path grids"
        )
    kpath = kpath_up

    waypoints = _parse_kpoints_qe(kpoints_path)
    labels, positions = _build_tick_data(waypoints, kpath)

    in_window = (
        ((bands_up.max(axis=0) >= elim[0]) & (bands_up.min(axis=0) <= elim[1]))
        | ((bands_dw.max(axis=0) >= elim[0]) & (bands_dw.min(axis=0) <= elim[1]))
    )
    bands_up = bands_up[:, in_window]
    bands_dw = bands_dw[:, in_window]

    x_total = float(kpath[-1] - kpath[0])
    gap_half = x_total * DEFAULT_GAP_FRAC

    boundary_pos = [p for lbl, p in zip(labels, positions) if lbl not in HELPER_LABELS]
    gap_pos = [p for lbl, p in zip(labels, positions) if lbl in GAP_LABELS]
    tick_labels = [_format_tick_label(lbl) for lbl in labels]

    fig, ax = plt.subplots(figsize=fig_size, constrained_layout=True)

    # grey shading for HPKOT path intervals (skip helper-label boundaries)
    for i in range(len(positions) - 1):
        if labels[i] in HELPER_LABELS or labels[i + 1] in HELPER_LABELS:
            continue
        ax.axvspan(positions[i], positions[i + 1], color=GREY_COLOR, lw=0, zorder=0)

    # spin-down (red) below spin-up (black)
    for ib in range(bands_dw.shape[1]):
        ax.plot(kpath, bands_dw[:, ib], color=BAND_DOWN_COLOR, lw=BAND_LW, zorder=2)
    for ib in range(bands_up.shape[1]):
        ax.plot(kpath, bands_up[:, ib], color=BAND_UP_COLOR, lw=BAND_LW, zorder=3)

    # white gap strip for k|k' discontinuities
    for pos in gap_pos:
        ax.axvspan(pos - gap_half, pos + gap_half, color="white", zorder=4, lw=0)

    # vertical lines at high-sym points
    for pos in boundary_pos:
        ax.axvline(x=pos, color=VLINE_COLOR, lw=VLINE_LW, zorder=5)

    # double vertical lines at gap edges
    for pos in gap_pos:
        ax.axvline(x=pos - gap_half, color="black", lw=0.8, zorder=5)
        ax.axvline(x=pos + gap_half, color="black", lw=0.8, zorder=5)

    ax.axhline(y=0, color=FERMI_COLOR, lw=FERMI_LW, ls="--", zorder=1)

    visible = [(p, t) for p, t in zip(positions, tick_labels) if kpath[0] <= p <= kpath[-1]]
    ax.set_xticks([p for p, _ in visible])
    ax.set_xticklabels([t for _, t in visible], fontsize=FONT_SIZE)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)
    ax.set_xlim(kpath[0], kpath[-1])
    ax.set_ylim(elim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.set_ylabel(r"E - E$_\mathrm{F}$ (eV)", fontsize=FONT_SIZE + 1)

    fig.savefig(output_path, dpi=800, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Plot spin-resolved AlterSeeK band output from QE bands.x files."
    )
    parser.add_argument(
        "--config",
        default=None,
        help="TOML config file. Defaults to alterband_qe.toml if present.",
    )
    parser.add_argument("-o", "--output", default=None, help="Override output file path.")
    args = parser.parse_args(argv)

    config_path = Path(args.config) if args.config else Path("alterband_qe.toml")
    config = _read_plot_config(config_path) if args.config or config_path.exists() else {}

    emin = float(config.get("emin", DEFAULT_ELIM[0]))
    emax = float(config.get("emax", DEFAULT_ELIM[1]))
    fig_width = float(config.get("fig_width", DEFAULT_FIG_SIZE[0]))
    fig_height = float(config.get("fig_height", DEFAULT_FIG_SIZE[1]))
    out_path = args.output or config.get("output", "alterband_qe.png")

    output = plot_alterband_qe(
        band_up=config.get("band_up", "band_up.gnu"),
        band_down=config.get("band_down", "band_down.gnu"),
        kpoints_qe=config.get("kpoints_qe", "KPOINTS_modified_qe"),
        output=out_path,
        fermi_ev=float(config.get("fermi_ev", 0.0)),
        elim=(emin, emax),
        fig_size=(fig_width, fig_height),
    )
    print(f"Band plot written to: {output}")


if __name__ == "__main__":
    main()
