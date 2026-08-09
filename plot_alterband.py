#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from VASPKIT reformatted data."""

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
DEFAULT_GAP_WIDTH_INCHES = 0.05
DEFAULT_FIG_SIZE = (12.0, 5.0)
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
    "klabels", "up", "down", "band_up", "band_down", "output",
    "emin", "emax", "fig_width", "fig_height",
    "gap_width_inches", "lattice_type", "split_panels",
    "rotate_xtick_labels", "xtick_rotation", "save_pdf",
}
_STRING_CONFIG_KEYS = {
    "klabels", "up", "down", "band_up", "band_down", "output", "lattice_type",
}
_NUMBER_CONFIG_KEYS = {
    "emin", "emax", "fig_width", "fig_height",
    "gap_width_inches", "xtick_rotation",
}

GREEK_LABELS = {
    "GAMMA": r"$\Gamma$",
    "\u0393": r"$\Gamma$",
    "\u8795": r"$\Gamma$",
    "DELTA": r"$\Delta$",
    "\u0394": r"$\Delta$",
    "LAMBDA": r"$\Lambda$",
    "\u039b": r"$\Lambda$",
    "SIGMA": r"$\Sigma$",
    "\u03a3": r"$\Sigma$",
}


def _read_klabels(path: Path) -> tuple[list[str], list[float]]:
    labels: list[str] = []
    positions: list[float] = []
    with path.open(encoding="utf-8-sig") as f:
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


VASPKIT_TRUNCATED_LABEL_FIXES = {
    # VASPKIT can drop the trailing zero in combined labels such as
    # LAMBDA_0'|G_.  Repair only the affected HPKOT labels.
    "oI2": {"G_": "G_2"},
    "oI3": {"G_": "G_0"},
    "oF2": {"Q_": "Q_0"},
}


def _fix_klabels_missing_merge(labels: list[str], kpoints_path: Path | None) -> list[str]:
    """Fix vaspkit 1.5.0 bug: boundary label written as start-only instead of 'end|start'.

    Triggered when adjacent segment endpoints share kx and kz but differ in ky.
    Detects the condition from KPOINTS coordinates, so works for any label names.
    """
    if kpoints_path is None or not kpoints_path.exists():
        return labels
    segments: list[list] = []
    buf: list = []
    with kpoints_path.open(encoding="utf-8-sig") as f:
        for i, line in enumerate(f):
            if i < 4:
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    buf.append((parts[3], (float(parts[0]), float(parts[1]), float(parts[2]))))
                    if len(buf) == 2:
                        segments.append(buf)
                        buf = []
                except ValueError:
                    pass
    if not segments:
        return labels
    # Build expected KLABELS sequence from segment structure.
    expected: list[str] = [segments[0][0][0]]
    for i in range(len(segments) - 1):
        (el, ek), (sl, sk) = segments[i][1], segments[i + 1][0]
        dx, dy, dz = abs(ek[0] - sk[0]), abs(ek[1] - sk[1]), abs(ek[2] - sk[2])
        all_same = dx < 1e-4 and dy < 1e-4 and dz < 1e-4
        expected.append(el if (all_same and el == sl) else f"{el}|{sl}")
    expected.append(segments[-1][1][0])
    if len(expected) != len(labels):
        return labels  # length mismatch means our parsing diverges; don't corrupt
    # Replace only where actual is the start-label tail of a buggy expected merge.
    return [exp if (act != exp and exp.endswith(f"|{act}")) else act
            for act, exp in zip(labels, expected)]


def _fix_vaspkit_truncated_label(label: str, lattice_type: str | None) -> str:
    if lattice_type is None:
        return label
    fixes = VASPKIT_TRUNCATED_LABEL_FIXES.get(_canonical_lattice_type(lattice_type), {})
    if not fixes:
        return label

    fixed_parts = []
    for part in str(label).split("|"):
        prime_count = len(part) - len(part.rstrip("'"))
        base = part[:-prime_count] if prime_count else part
        fixed_parts.append(fixes.get(base, base) + "'" * prime_count)
    return "|".join(fixed_parts)


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


# Keep this function's body identical to the copy in plot_alterband_qe.py;
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


def _format_tick_label(label: str) -> str:
    """Return labels with mathtext only where Greek/subscripts are needed."""
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


def _canonical_lattice_type(lattice_type: str) -> str:
    lattice_type = str(lattice_type).strip()
    if lattice_type.upper().startswith("ORCF"):
        return lattice_type
    try:
        from alterseek.lattice_kpoints import canonical_lattice_type

        return canonical_lattice_type(lattice_type)
    except Exception:
        return lattice_type


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
def plot_alterband(
    *,
    klabels: str | Path = "KLABELS",
    band_up: str | Path = "REFORMATTED_BAND_UP.dat",
    band_down: str | Path = "REFORMATTED_BAND_DW.dat",
    output: str | Path = "alterband.png",
    elim: tuple[float, float] = DEFAULT_ELIM,
    fig_size: tuple[float, float] = DEFAULT_FIG_SIZE,
    gap_width_inches: float = DEFAULT_GAP_WIDTH_INCHES,
    lattice_type: str | None = None,
    split_panels: int = 1,
    rotate_xtick_labels: bool = False,
    xtick_rotation: float = 45.0,
    save_pdf: bool = False,
) -> Path:
    """Create the spin-resolved band plot and return the output path."""
    klabels_path = Path(klabels)
    band_up_path = Path(band_up)
    band_down_path = Path(band_down)
    output_path = Path(output)

    labels, positions = _read_klabels(klabels_path)
    original_labels = labels
    labels = _fix_klabels_missing_merge(labels, klabels_path.parent / "KPOINTS")
    labels = [_fix_vaspkit_truncated_label(label, lattice_type) for label in labels]
    corrections = [
        f"{old} -> {new}"
        for old, new in zip(original_labels, labels)
        if old != new
    ]
    if corrections:
        print(
            "[Note] Corrected VASPKIT label(s) for plotting only; KLABELS was "
            f"not modified: {', '.join(corrections)}"
        )
    x_total = positions[-1] - positions[0]
    if x_total <= 0:
        raise ValueError("KLABELS positions must increase from first to last entry")
    if gap_width_inches < 0:
        raise ValueError("gap_width_inches must be non-negative")

    tick_labels = [_format_tick_label(label) for label in labels]

    up = np.loadtxt(band_up_path, skiprows=1, ndmin=2)
    dw = np.loadtxt(band_down_path, skiprows=1, ndmin=2)
    kpath = up[:, 0]
    bands_up = up[:, 1:]
    bands_dw = dw[:, 1:]

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
    parser.add_argument("--up", default=None, help="Spin-up band data file.")
    parser.add_argument("--down", default=None, help="Spin-down band data file.")
    parser.add_argument("-o", "--output", default=None, help="Output file.")
    parser.add_argument("--emin", type=float, default=None, help="Minimum plotted energy.")
    parser.add_argument("--emax", type=float, default=None, help="Maximum plotted energy.")
    parser.add_argument("--fig-width", type=float, default=None, help="Figure width in inches.")
    parser.add_argument("--fig-height", type=float, default=None, help="Figure height in inches.")
    parser.add_argument("--lattice-type", default=None, help="HPKOT lattice type such as tI1 or mC2.")
    parser.add_argument("--split-panels", type=int, default=None, help="Use 1, 2, or 3 stacked panels.")
    parser.add_argument(
        "--gap-width-inches",
        type=float,
        default=None,
        help="Full visual width of each k|k' gap in inches.",
    )
    parser.add_argument(
        "--rotate-xtick-labels",
        action="store_true",
        default=None,
        help="Rotate x-axis tick labels.",
    )
    parser.add_argument("--xtick-rotation", type=float, default=None, help="X tick rotation angle.")
    parser.add_argument(
        "--save-pdf",
        action="store_true",
        default=None,
        help="Also save a PDF copy alongside the primary output.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    config_path = Path(args.config) if args.config else Path("alterband.toml")
    def option(name: str, default: Any) -> Any:
        arg_value = getattr(args, name)
        if arg_value is not None:
            return arg_value
        return config.get(name, default)

    try:
        config = (
            _validate_plot_config(_read_plot_config(config_path), config_path)
            if args.config or config_path.exists() else {}
        )
        emin = float(option("emin", DEFAULT_ELIM[0]))
        emax = float(option("emax", DEFAULT_ELIM[1]))
        fig_width = float(option("fig_width", DEFAULT_FIG_SIZE[0]))
        fig_height = float(option("fig_height", DEFAULT_FIG_SIZE[1]))
        output = plot_alterband(
            klabels=option("klabels", "KLABELS"),
            band_up=option("up", config.get("band_up", "REFORMATTED_BAND_UP.dat")),
            band_down=option("down", config.get("band_down", "REFORMATTED_BAND_DW.dat")),
            output=option("output", "alterband.png"),
            elim=(emin, emax),
            fig_size=(fig_width, fig_height),
            gap_width_inches=float(option("gap_width_inches", DEFAULT_GAP_WIDTH_INCHES)),
            lattice_type=option("lattice_type", None),
            split_panels=int(option("split_panels", 1)),
            rotate_xtick_labels=bool(option("rotate_xtick_labels", False)),
            xtick_rotation=float(option("xtick_rotation", 45.0)),
            save_pdf=bool(option("save_pdf", False)),
        )
    except (ValueError, OSError) as exc:
        raise SystemExit(f"[Error] {exc}") from exc
    print(f"Band plot written to: {output}")


if __name__ == "__main__":
    main()
