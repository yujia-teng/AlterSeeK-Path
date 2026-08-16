#!/usr/bin/env python3
"""Plot spin-resolved AlterSeeK band output from an ABINIT _EIG file."""

from __future__ import annotations

import argparse
from functools import wraps
import math
from pathlib import Path
import re
import tomllib
from typing import Any

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
HARTREE_TO_EV = 27.211386245988

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
    "eig", "kpoints_abinit", "poscar", "abo", "output", "fermi_ev",
    "emin", "emax", "fig_width", "fig_height",
    "gap_width_inches", "split_panels", "rotate_xtick_labels",
    "xtick_rotation", "save_pdf",
}
_STRING_CONFIG_KEYS = {
    "eig", "kpoints_abinit", "poscar", "abo", "output",
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


def _read_plot_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("rb") as f:
        data = tomllib.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Config file must contain key-value settings: {path}")
    return data


# Keep this function identical to plot_alterband.py/plot_alterband_qe.py except for the module-level key sets.
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


def _read_poscar_lattice(path: Path) -> np.ndarray:
    """Read the 3x3 real-space lattice (Angstrom) from a VASP POSCAR."""
    with path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()
    if len(lines) < 5:
        raise ValueError(f"{path} is too short to contain a lattice")
    scale = float(lines[1].split()[0])
    lattice = np.array(
        [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)],
        dtype=float,
    )
    return scale * lattice


def _reciprocal_lattice(real_lattice: np.ndarray) -> np.ndarray:
    """Rows are reciprocal lattice vectors b1,b2,b3 (1/Angstrom), the same
    2*pi convention ABINIT's own abinit_eignc_to_bandstructure.py post-
    processing script uses to convert fractional k-points to Cartesian."""
    return 2.0 * np.pi * np.linalg.inv(real_lattice).T


_FERMI_HARTREE_RE = re.compile(
    r"Fermi \(or HOMO\) energy \(hartree\)\s*=\s*([-\d.Ee+]+)"
)


def _read_fermi_ev(path: Path) -> float:
    """Read the Fermi (or HOMO) energy from an ABINIT .abo output file and
    convert it to eV. When a .abo covers multiple datasets, ABINIT only
    computes and prints this line once, during the self-consistent (SCF)
    dataset -- a fixed-density NSCF band-path run has no physically
    meaningful occupation to redetermine it from, so it just carries the
    SCF value forward. The _EIG file itself carries no reference energy at
    all, unlike VASPKIT's reformatted output or QE's bands.x .gnu files."""
    with path.open(encoding="utf-8-sig") as f:
        text = f.read()
    match = _FERMI_HARTREE_RE.search(text)
    if not match:
        raise ValueError(f"{path}: no 'Fermi (or HOMO) energy' line found")
    return float(match.group(1)) * HARTREE_TO_EV


_SPIN_HEADER_RE = re.compile(
    r"Eigenvalues \(hartree\) for nkpt=\s*\d+\s*k points,\s*SPIN\s+(UP|DOWN):"
)
_KPT_HEADER_RE = re.compile(
    r"kpt#\s*\d+,\s*nband=\s*(\d+),.*kpt=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)"
)


def _read_eig_bands(
    path: Path, fermi_ev: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse an ABINIT plain-text _EIG file (both spin channels required).

    Returns (kpoints_frac, bands_up, bands_dw): kpoints_frac has shape
    (nkpt, 3) in reduced coordinates; each bands_* array has shape
    (nkpt, nband), energies converted to eV and shifted by fermi_ev.
    """
    with path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()

    per_spin_kpts: dict[str, list] = {"UP": [], "DOWN": []}
    per_spin_bands: dict[str, list] = {"UP": [], "DOWN": []}
    current_spin: str | None = None
    current_energies: list[float] = []
    nband: int | None = None

    def flush():
        if not current_energies:
            return
        if nband is not None and len(current_energies) != nband:
            raise ValueError(
                f"{path}: expected {nband} eigenvalues per k-point, got "
                f"{len(current_energies)}"
            )
        per_spin_bands[current_spin].append(list(current_energies))

    for line in lines:
        stripped = line.strip()
        spin_header = _SPIN_HEADER_RE.match(stripped)
        if spin_header:
            if current_spin is not None:
                flush()
            current_spin = spin_header.group(1)
            current_energies = []
            continue
        kpt_header = _KPT_HEADER_RE.match(stripped)
        if kpt_header:
            if current_spin is None:
                raise ValueError(f"{path}: eigenvalue block found before a SPIN header")
            flush()
            current_energies = []
            nband = int(kpt_header.group(1))
            per_spin_kpts[current_spin].append(
                [float(kpt_header.group(i)) for i in (2, 3, 4)]
            )
            continue
        if stripped and current_spin is not None and nband is not None:
            current_energies.extend(float(x) for x in stripped.split())
    flush()

    if not per_spin_bands["UP"] or not per_spin_bands["DOWN"]:
        raise ValueError(f"{path} must contain both SPIN UP and SPIN DOWN blocks")
    if len(per_spin_kpts["UP"]) != len(per_spin_kpts["DOWN"]):
        raise ValueError(f"{path}: SPIN UP and SPIN DOWN have different k-point counts")
    if not np.allclose(per_spin_kpts["UP"], per_spin_kpts["DOWN"], atol=1e-6):
        raise ValueError(f"{path}: SPIN UP and SPIN DOWN k-point grids differ")

    kpoints_frac = np.array(per_spin_kpts["UP"], dtype=float)
    bands_up = np.array(per_spin_bands["UP"], dtype=float) * HARTREE_TO_EV - fermi_ev
    bands_dw = np.array(per_spin_bands["DOWN"], dtype=float) * HARTREE_TO_EV - fermi_ev
    return kpoints_frac, bands_up, bands_dw


def _parse_kpoints_abinit(path: Path) -> list[tuple[str, int, int]]:
    """Parse KPOINTS_alter_abinit written by write_kpoints_file_abinit.

    Returns list of (label, ndivk, k_index) for each waypoint, mirroring
    _parse_kpoints_qe's (label, ninterp, k_index) shape. k_index is the
    0-based index into the _EIG file's k-point array.
    """
    with path.open(encoding="utf-8-sig") as f:
        lines = f.readlines()

    kptopt_value: int | None = None
    ndivk_values: list[int] = []
    kptbounds_start: int | None = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("kptopt"):
            kptopt_value = int(stripped.split()[1])
        elif stripped.startswith("ndivk"):
            ndivk_values = [int(x) for x in stripped.split()[1:]]
        elif stripped.startswith("kptbounds"):
            kptbounds_start = i + 1
            break
    if kptopt_value is None or not ndivk_values or kptbounds_start is None:
        raise ValueError(f"{path}: missing kptopt/ndivk/kptbounds")
    if kptopt_value >= 0:
        raise ValueError(f"{path}: kptopt must be negative (band-structure path mode)")

    n_segments = -kptopt_value
    if len(ndivk_values) != n_segments:
        raise ValueError(
            f"{path}: kptopt implies {n_segments} segments but ndivk has "
            f"{len(ndivk_values)} values"
        )

    waypoints_raw: list[str] = []
    for line in lines[kptbounds_start:]:
        if len(waypoints_raw) >= n_segments + 1:
            break
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split("#", 1)
        fields = parts[0].split()
        if len(fields) < 3:
            break  # end of the kptbounds block
        waypoints_raw.append(parts[1].strip() if len(parts) > 1 else "")

    if len(waypoints_raw) != n_segments + 1:
        raise ValueError(
            f"{path} declares {n_segments} segments ({n_segments + 1} waypoints) "
            f"but contains only {len(waypoints_raw)}"
        )

    result: list[tuple[str, int, int]] = []
    cum = 0
    for idx, label in enumerate(waypoints_raw):
        ni = ndivk_values[idx] if idx < len(ndivk_values) else 0
        result.append((label, ni, cum))
        cum += ni
    return result


def _kpath_distance(kpoints_frac: np.ndarray, reciprocal_lattice: np.ndarray) -> np.ndarray:
    """Cumulative Cartesian path distance (1/Angstrom) through the k-points,
    in file order -- including k<->k' "break" jumps, exactly as VASPKIT and
    QE's bands.x already do for the VASP/QE plotters (the jump still gets a
    real physical width; the k|k' gap marker is drawn on top of it)."""
    cart = kpoints_frac @ reciprocal_lattice
    deltas = np.linalg.norm(np.diff(cart, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(deltas)])


def _format_tick_label(label: str) -> str:
    if "|" in label:
        return "|".join(_format_tick_label(part) for part in label.split("|"))
    if "/" in label:
        return "/".join(_format_tick_label(part) for part in label.split("/"))
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
    that land at the same k-path position is merged into a single "A|B"
    tick -- otherwise matplotlib silently drops one of the two same-position
    labels (mirrors plot_alterband.py/plot_alterband_qe.py).
    """
    for label, _ni, k_idx in waypoints:
        if k_idx >= len(kpath):
            raise ValueError(
                f"Waypoint '{label}' has k-index {k_idx} but the _EIG file has only "
                f"{len(kpath)} k-points — KPOINTS_alter_abinit and _EIG file mismatch"
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
def plot_alterband_abinit(
    *,
    eig: str | Path = "EIG",
    kpoints_abinit: str | Path = "KPOINTS_alter_abinit",
    poscar: str | Path = "POSCAR",
    abo: str | Path | None = "abo",
    output: str | Path = "alterband_abinit.png",
    fermi_ev: float | None = None,
    elim: tuple[float, float] = DEFAULT_ELIM,
    fig_size: tuple[float, float] = DEFAULT_FIG_SIZE,
    gap_width_inches: float = DEFAULT_GAP_WIDTH_INCHES,
    split_panels: int = 1,
    rotate_xtick_labels: bool = False,
    xtick_rotation: float = 45.0,
    save_pdf: bool = False,
) -> Path:
    """Create the ABINIT spin-resolved band plot and return the output path.

    fermi_ev, if not given directly, is read from the ABINIT .abo output
    file named by abo and converted from hartree -- the user only needs to
    have that file present, not compute or type the shift by hand."""
    eig_path = Path(eig)
    kpoints_path = Path(kpoints_abinit)
    poscar_path = Path(poscar)
    output_path = Path(output)

    if fermi_ev is None:
        if abo is None:
            raise ValueError(
                "fermi_ev was not given and no .abo file was supplied to read it from"
            )
        abo_path = Path(abo)
        if not abo_path.exists():
            raise ValueError(
                f"{abo_path} not found -- supply fermi_ev directly, or point "
                "'abo' at the ABINIT .abo output file to read it from"
            )
        fermi_ev = _read_fermi_ev(abo_path)

    kpoints_frac, bands_up, bands_dw = _read_eig_bands(eig_path, fermi_ev)
    reciprocal_lattice = _reciprocal_lattice(_read_poscar_lattice(poscar_path))
    kpath = _kpath_distance(kpoints_frac, reciprocal_lattice)

    waypoints = _parse_kpoints_abinit(kpoints_path)
    labels, positions = _build_tick_data(waypoints, kpath)
    x_total = positions[-1] - positions[0]
    if x_total <= 0:
        raise ValueError("KPOINTS_alter_abinit positions must increase from first to last entry")
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

    # Every panel has the same physical width even when it covers a different k-path range, so a gap scaled to the global x_total looks too wide in a panel showing only a narrow slice.
    # Scale each gap to its panel's local range so its physical printed width stays constant across panel counts.
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
        description="Plot spin-resolved AlterSeeK band output from an ABINIT _EIG file."
    )
    parser.add_argument(
        "--config",
        default=None,
        help="TOML config file. Defaults to alterband_abinit.toml if present.",
    )
    parser.add_argument("-o", "--output", default=None, help="Override output file path.")
    args = parser.parse_args(argv)

    config_path = Path(args.config) if args.config else Path("alterband_abinit.toml")
    try:
        config = (
            _validate_plot_config(_read_plot_config(config_path), config_path)
            if args.config or config_path.exists() else {}
        )
        emin = float(config.get("emin", DEFAULT_ELIM[0]))
        emax = float(config.get("emax", DEFAULT_ELIM[1]))
        fig_width = float(config.get("fig_width", DEFAULT_FIG_SIZE[0]))
        fig_height = float(config.get("fig_height", DEFAULT_FIG_SIZE[1]))
        out_path = args.output or config.get("output", "alterband_abinit.png")
        output = plot_alterband_abinit(
            eig=config.get("eig", "EIG"),
            kpoints_abinit=config.get("kpoints_abinit", "KPOINTS_alter_abinit"),
            poscar=config.get("poscar", "POSCAR"),
            abo=config.get("abo", "abo"),
            output=out_path,
            fermi_ev=(float(config["fermi_ev"]) if "fermi_ev" in config else None),
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
