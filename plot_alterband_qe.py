#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from Quantum ESPRESSO bands.x .gnu files."""

from __future__ import annotations

import argparse
from functools import wraps
import math
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10 fallback
    tomllib = None

import matplotlib as mpl

mpl.use("Agg")

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator


DEFAULT_ELIM = (-2.0, 2.0)
DEFAULT_FIG_SIZE = (12.0, 5.0)
DEFAULT_GAP_WIDTH_INCHES = 0.05
DEFAULT_PANEL_GAP = 0.08
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

_PLOT_STYLE = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "xtick.direction": "in",
    "ytick.direction": "in",
}


def _plot_style(func):
    @wraps(func)
    def wrapped(*args, **kwargs):
        with mpl.rc_context(_PLOT_STYLE):
            return func(*args, **kwargs)
    return wrapped

_PLOT_CONFIG_KEYS = {
    "band_up", "band_down", "kpoints_qe", "output", "fermi_ev",
    "emin", "emax", "fig_width", "fig_height",
    "gap_width_inches", "split_panels", "rotate_xtick_labels",
    "xtick_rotation", "save_pdf",
}
_STRING_CONFIG_KEYS = {
    "band_up", "band_down", "kpoints_qe", "output",
}
_NUMBER_CONFIG_KEYS = {
    "fermi_ev", "emin", "emax", "fig_width", "fig_height",
    "gap_width_inches", "xtick_rotation",
}

GREEK_LABELS = {
    "GAMMA": r"$\Gamma$",
    "Γ": r"$\Gamma$",
    "\u8795": r"$\Gamma$",  # Legacy mis-encoded Gamma.
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
    with path.open(encoding="utf-8-sig") as f:
        for raw_line in f:
            line = raw_line.split("#", 1)[0].strip()
            if not line or line.startswith("["):
                continue
            if "=" not in line:
                raise ValueError(f"Invalid config line in {path}: {raw_line.rstrip()}")
            key, value = line.split("=", 1)
            config[key.strip()] = _parse_simple_toml_value(value)
    return config


# Keep this function's body identical to the copy in plot_alterband.py;
# only the module-level key sets differ between the two plotters.
def _validate_plot_config(config: dict[str, Any], path: Path) -> dict[str, Any]:
    unknown = sorted(set(config) - _PLOT_CONFIG_KEYS)
    if unknown:
        raise ValueError(f"Unknown setting(s) in {path}: {', '.join(unknown)}")
    for key in _STRING_CONFIG_KEYS:
        if key in config and not isinstance(config[key], str):
            raise ValueError(f"{key} in {path} must be a string")
    for key in _NUMBER_CONFIG_KEYS:
        if key in config and (isinstance(config[key], bool) or
                              not isinstance(config[key], (int, float))):
            raise ValueError(f"{key} in {path} must be a number")
    if "split_panels" in config and (
            isinstance(config["split_panels"], bool) or
            not isinstance(config["split_panels"], int) or
            config["split_panels"] not in {1, 2, 3}):
        raise ValueError(f"split_panels in {path} must be 1, 2, or 3")
    if ("rotate_xtick_labels" in config and
            not isinstance(config["rotate_xtick_labels"], bool)):
        raise ValueError(f"rotate_xtick_labels in {path} must be true or false")
    if "save_pdf" in config and not isinstance(config["save_pdf"], bool):
        raise ValueError(f"save_pdf in {path} must be true or false")
    return config


def _read_gnu_bands(path: Path, fermi_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """Parse QE bands.x .gnu file (band-major, blank-line separated).

    Returns (kpath, bands) where bands.shape = (nk, nbands), energies shifted by fermi_ev.
    """
    blocks: list[np.ndarray] = []
    current: list[list[float]] = []
    with path.open(encoding="utf-8-sig") as f:
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
    """Parse KPOINTS_alter_qe written by write_kpoints_file_qe.

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
    if start + 1 >= len(all_lines):
        raise ValueError(f"{path} is truncated right after the K_POINTS card")
    if not all_lines[start + 1].split():
        raise ValueError(f"{path}: missing waypoint-count line after the K_POINTS card")

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
        fields = parts[0].split()
        if len(fields) < 4:
            raise ValueError(
                f"Malformed waypoint line in {path}: '{stripped}' "
                "(expected 'x y z ni  ! LABEL')"
            )
        ni = int(fields[3])
        label = parts[1].strip() if len(parts) > 1 else ""
        waypoints_raw.append((label, ni))

    if len(waypoints_raw) != n_waypoints:
        raise ValueError(
            f"{path} declares {n_waypoints} waypoints but contains only "
            f"{len(waypoints_raw)}"
        )

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


_COINCIDENT_POSITION_TOL = 1e-4


def _build_tick_data(
    waypoints: list[tuple[str, int, int]], kpath: np.ndarray
) -> tuple[list[str], list[float]]:
    """Convert waypoints to tick labels and x-positions.

    Adjacent k / k' pairs (in either order) are merged into a single k|k' or
    k'|k gap tick at their midpoint. Any other pair of adjacent waypoints
    that land at the same k-path position (the segment revisits a point
    that is physically identical up to a reciprocal lattice vector, e.g.
    a path landing back on a zone-boundary point under a different label)
    is merged into a single "A|B" tick -- otherwise matplotlib silently
    drops one of the two same-position labels (mirrors
    _fix_klabels_missing_merge in plot_alterband.py, which handles the
    same coincident-point case for VASP/VASPKIT-derived KLABELS).
    """
    for label, _ni, k_idx in waypoints:
        if k_idx >= len(kpath):
            raise ValueError(
                f"Waypoint '{label}' has k-index {k_idx} but .gnu file has only "
                f"{len(kpath)} k-points — KPOINTS_alter_qe and .gnu file mismatch"
            )

    labels: list[str] = []
    positions: list[float] = []
    i = 0
    while i < len(waypoints):
        label, _ni, k_idx = waypoints[i]
        if (
            label in ("k", "k'")
            and i + 1 < len(waypoints)
            and waypoints[i + 1][0] in ("k", "k'")
            and waypoints[i + 1][0] != label
        ):
            x_first = kpath[waypoints[i][2]]
            x_second = kpath[waypoints[i + 1][2]]
            labels.append(f"{label}|{waypoints[i + 1][0]}")
            positions.append((x_first + x_second) / 2.0)
            i += 2
            continue
        x_pos = float(kpath[k_idx])
        next_label, _next_ni, next_k_idx = (
            waypoints[i + 1] if i + 1 < len(waypoints) else (None, None, None)
        )
        if (
            next_label is not None
            and label not in ("k", "k'")
            and next_label not in ("k", "k'")
            and next_k_idx < len(kpath)
            and abs(float(kpath[next_k_idx]) - x_pos) < _COINCIDENT_POSITION_TOL
        ):
            labels.append(label if next_label == label else f"{label}|{next_label}")
            positions.append(x_pos)
            i += 2
            continue
        labels.append(label)
        positions.append(x_pos)
        i += 1
    return labels, positions


def _is_valid_split_label(label: str) -> bool:
    return label not in HELPER_LABELS


def _split_indices(labels: list[str], split_panels: int) -> list[int]:
    if split_panels in (None, 1):
        return []
    if split_panels not in {2, 3}:
        raise ValueError("split_panels must be 1, 2, or 3")

    candidates = [
        i for i in range(1, len(labels) - 1)
        if _is_valid_split_label(labels[i])
    ]
    if not candidates:
        return []

    targets = [math.ceil(len(labels) / split_panels * i) for i in range(1, split_panels)]
    selected: list[int] = []
    for target in targets:
        options = [idx for idx in candidates if idx not in selected]
        if not options:
            break

        def score(idx: int) -> tuple[float, float, int]:
            trial = sorted(selected + [idx])
            bounds = [0, *trial, len(labels) - 1]
            widths = [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]
            return (abs(idx - target), max(widths) - min(widths), idx)

        selected.append(min(options, key=score))

    return sorted(selected)


def _panel_ranges(labels: list[str], positions: list[float], split_panels: int):
    splits = _split_indices(labels, split_panels)
    bounds = [0, *splits, len(labels) - 1]
    return [(positions[bounds[i]], positions[bounds[i + 1]]) for i in range(len(bounds) - 1)]


def _draw_panel(
    ax,
    *,
    labels: list[str],
    positions: list[float],
    tick_labels: list[str],
    kpath: np.ndarray,
    bands_up: np.ndarray,
    bands_dw: np.ndarray,
    elim: tuple[float, float],
    xlim: tuple[float, float],
    gap_half: float,
    font_size: int,
    rotate_xtick_labels: bool,
    xtick_rotation: float,
) -> None:
    boundary_pos = [p for label, p in zip(labels, positions) if label not in HELPER_LABELS]
    gap_pos = [p for label, p in zip(labels, positions) if label in GAP_LABELS]

    for i in range(len(positions) - 1):
        left, right = positions[i], positions[i + 1]
        if right < xlim[0] or left > xlim[1]:
            continue
        if labels[i] in HELPER_LABELS or labels[i + 1] in HELPER_LABELS:
            continue
        ax.axvspan(left, right, color=GREY_COLOR, lw=0, zorder=0)

    for ib in range(bands_dw.shape[1]):
        ax.plot(kpath, bands_dw[:, ib], color=BAND_DOWN_COLOR, lw=BAND_LW, zorder=2)
    for ib in range(bands_up.shape[1]):
        ax.plot(kpath, bands_up[:, ib], color=BAND_UP_COLOR, lw=BAND_LW, zorder=3)

    for pos in gap_pos:
        ax.axvspan(pos - gap_half, pos + gap_half, color="white", zorder=4, lw=0)

    for pos in boundary_pos:
        if xlim[0] <= pos <= xlim[1]:
            ax.axvline(x=pos, color=VLINE_COLOR, lw=VLINE_LW, zorder=5)

    for pos in gap_pos:
        if xlim[0] <= pos <= xlim[1]:
            ax.axvline(x=pos - gap_half, color="black", lw=0.8, zorder=5)
            ax.axvline(x=pos + gap_half, color="black", lw=0.8, zorder=5)

    visible_ticks = [
        (pos, lab) for pos, lab in zip(positions, tick_labels)
        if xlim[0] <= pos <= xlim[1]
    ]
    ax.axhline(y=0, color=FERMI_COLOR, lw=FERMI_LW, ls="--", zorder=1)
    ax.set_xticks([pos for pos, _ in visible_ticks])
    xtick_label_kwargs: dict[str, Any] = {"fontsize": font_size}
    if rotate_xtick_labels:
        xtick_label_kwargs.update(
            {"rotation": xtick_rotation, "ha": "right", "rotation_mode": "anchor"}
        )
    ax.set_xticklabels([lab for _, lab in visible_ticks], **xtick_label_kwargs)
    ax.set_xlim(xlim)
    ax.set_ylim(elim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=font_size)


@_plot_style
def plot_alterband_qe(
    *,
    band_up: str | Path = "band_up.gnu",
    band_down: str | Path = "band_down.gnu",
    kpoints_qe: str | Path = "KPOINTS_alter_qe",
    output: str | Path = "alterband_qe.png",
    fermi_ev: float = 0.0,
    elim: tuple[float, float] = DEFAULT_ELIM,
    fig_size: tuple[float, float] = DEFAULT_FIG_SIZE,
    gap_width_inches: float = DEFAULT_GAP_WIDTH_INCHES,
    split_panels: int = 1,
    rotate_xtick_labels: bool = False,
    xtick_rotation: float = 45.0,
    save_pdf: bool = False,
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
    x_total = positions[-1] - positions[0]
    if x_total <= 0:
        raise ValueError("KPOINTS_alter_qe positions must increase from first to last entry")
    if gap_width_inches < 0:
        raise ValueError("gap_width_inches must be non-negative")

    tick_labels = [_format_tick_label(lbl) for lbl in labels]

    in_window = (
        ((bands_up.max(axis=0) >= elim[0]) & (bands_up.min(axis=0) <= elim[1]))
        | ((bands_dw.max(axis=0) >= elim[0]) & (bands_dw.min(axis=0) <= elim[1]))
    )
    bands_up = bands_up[:, in_window]
    bands_dw = bands_dw[:, in_window]

    ranges = _panel_ranges(labels, positions,
                           1 if split_panels is None else int(split_panels))
    n_panels = len(ranges)
    total_size = fig_size
    if n_panels > 1:
        total_size = (fig_size[0], fig_size[1] * n_panels)
    font_size = FONT_SIZE + (2 if n_panels > 1 else 0)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=total_size,
        sharey=True,
        squeeze=False,
        constrained_layout=True,
        gridspec_kw={"hspace": DEFAULT_PANEL_GAP},
    )
    flat_axes = list(axes[:, 0])

    # Every panel is rendered at the same physical figure width regardless of
    # how much k-path range it covers, so a gap sized relative to the GLOBAL
    # x_total looks wrong on any panel that only shows a slice of it (a
    # narrower panel is effectively "zoomed in", so the same absolute gap
    # covers more of its visible width). Size each panel's gap relative to
    # its own local range instead, so the physical/printed gap width is the
    # same across all panels and panel counts.
    axis_width_inches = total_size[0]

    def _gap_half_for(xlim: tuple[float, float]) -> float:
        panel_span = xlim[1] - xlim[0]
        return 0.5 * gap_width_inches * panel_span / axis_width_inches

    for ax, xlim in zip(flat_axes, ranges):
        _draw_panel(
            ax,
            labels=labels,
            positions=positions,
            tick_labels=tick_labels,
            kpath=kpath,
            bands_up=bands_up,
            bands_dw=bands_dw,
            elim=elim,
            xlim=xlim,
            gap_half=_gap_half_for(xlim),
            font_size=font_size,
            rotate_xtick_labels=rotate_xtick_labels,
            xtick_rotation=xtick_rotation,
        )

    fig.supylabel(r"E - E$_\mathrm{F}$ (eV)", fontsize=font_size + 1)
    fig.savefig(output_path, dpi=800, bbox_inches="tight")
    if save_pdf and output_path.suffix.lower() != ".pdf":
        fig.savefig(output_path.with_suffix(".pdf"), dpi=800, bbox_inches="tight")
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
    try:
        config = (
            _validate_plot_config(_read_plot_config(config_path), config_path)
            if args.config or config_path.exists() else {}
        )
        emin = float(config.get("emin", DEFAULT_ELIM[0]))
        emax = float(config.get("emax", DEFAULT_ELIM[1]))
        fig_width = float(config.get("fig_width", DEFAULT_FIG_SIZE[0]))
        fig_height = float(config.get("fig_height", DEFAULT_FIG_SIZE[1]))
        out_path = args.output or config.get("output", "alterband_qe.png")
        output = plot_alterband_qe(
            band_up=config.get("band_up", "band_up.gnu"),
            band_down=config.get("band_down", "band_down.gnu"),
            kpoints_qe=config.get("kpoints_qe", "KPOINTS_alter_qe"),
            output=out_path,
            fermi_ev=float(config.get("fermi_ev", 0.0)),
            elim=(emin, emax),
            fig_size=(fig_width, fig_height),
            gap_width_inches=float(config.get("gap_width_inches", DEFAULT_GAP_WIDTH_INCHES)),
            split_panels=int(config.get("split_panels", 1)),
            rotate_xtick_labels=bool(config.get("rotate_xtick_labels", False)),
            xtick_rotation=float(config.get("xtick_rotation", 45.0)),
            save_pdf=bool(config.get("save_pdf", False)),
        )
    except (ValueError, OSError) as exc:
        raise SystemExit(f"[Error] {exc}") from exc
    print(f"Band plot written to: {output}")


if __name__ == "__main__":
    main()
