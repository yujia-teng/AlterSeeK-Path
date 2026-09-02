"""Shared figure I/O, label, and style helpers."""
from functools import wraps
import os

import matplotlib as mpl
import numpy as np

from .atomic_write import _atomic_open_text
from .lattice_kpoints import canonical_lattice_type


ALTERSEEK_PLOT_STYLE = {
    "mathtext.fontset": "stix",
}


def alterseek_plot_style(func):
    """Run one AlterSeeK figure stage without changing caller rcParams."""
    @wraps(func)
    def wrapped(*args, **kwargs):
        with mpl.rc_context(ALTERSEEK_PLOT_STYLE):
            return func(*args, **kwargs)
    return wrapped

GAMMA_LABEL = "\u0393"
POINT_COINCIDENCE_ATOL = 1e-6
BZ_SPECIAL_COLORS = {
    "orange": "#e68613",
    "purple": "#6b5596",
}
IBZ_FACE_COLORS = {
    "up_main": "salmon",
    "up_extra": "#e98f8f",
    "down_main": "cornflowerblue",
    "down_extra": "#91b2e8",
}
IBZ_UP_EXTRA_SECTOR_COLORS = ("#8b0000", "#d6604d", "#f4a582")
BZ_PATH_STYLE_OVERRIDES = {
    "cP1": {("M", "X_1"): {"color": "red", "ls": "--"}},
    "cF1": {("X", "W_2"): {"color": "red", "ls": "--"}},
    "tI1": {(GAMMA_LABEL, "N"): "orange"},
    "tI2": {(GAMMA_LABEL, "N"): "orange"},
    "oF1": {(GAMMA_LABEL, "L"): "orange"},
    "oF2": {(GAMMA_LABEL, "L"): "orange"},
    "oF3": {(GAMMA_LABEL, "L"): "orange"},
    "oI1": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oI2": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oI3": {
        (GAMMA_LABEL, "T"): "orange",
        (GAMMA_LABEL, "R"): "orange",
        (GAMMA_LABEL, "S"): "orange",
    },
    "oA1": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oA2": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oC1": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "oC2": {
        (GAMMA_LABEL, "S"): "orange",
        ("Z", "R"): "purple",
    },
    "hR1": {
        (GAMMA_LABEL, "L"): "orange",
        (GAMMA_LABEL, "F"): "orange",
    },
    "hR2": {(GAMMA_LABEL, "L"): "orange"},
    "mP1": {
        (GAMMA_LABEL, "B"): "orange",
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "Y_2"): "orange",
        ("Z", "D"): "purple",
        ("Z", "E"): "purple",
        ("Z", "C_2"): "purple",
    },
    "mC1": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "M_2"): "orange",
        (GAMMA_LABEL, "Y_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
    },
    "mC2": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
        ("M", "Y"): "purple",
    },
    "mC3": {
        (GAMMA_LABEL, "A"): "orange",
        (GAMMA_LABEL, "M_2"): "orange",
        (GAMMA_LABEL, "L_2"): "orange",
        (GAMMA_LABEL, "V_2"): "orange",
    },
}


def _get_bz_path_style(lattice_type, k1, k2):
    """Return the display style for a recommended HPKOT path segment."""
    style = {"color": "red", "ls": "-", "lw": 4.0, "alpha": 0.9}
    if not lattice_type:
        return style

    try:
        lattice_key = canonical_lattice_type(lattice_type)
    except Exception:
        lattice_key = lattice_type

    if lattice_key in {"aP2", "aP3"}:
        style["color"] = BZ_SPECIAL_COLORS["orange"]
        style["alpha"] = 0.95
        return style

    overrides = BZ_PATH_STYLE_OVERRIDES.get(lattice_key, {})
    override = overrides.get((k1, k2), overrides.get((k2, k1)))
    if override:
        if isinstance(override, dict):
            style.update(override)
            style["color"] = BZ_SPECIAL_COLORS.get(style["color"], style["color"])
        else:
            style["color"] = BZ_SPECIAL_COLORS.get(override, override)
        style["alpha"] = 0.95
    return style


def _figure_output_paths(output_path, extra_formats=None):
    """Return the requested figure output paths (``bz.png`` ->
    ``["bz.png", "bz.pdf"]``). Set ALTERSEEK_BZ_FORMATS (e.g. "png,pdf") to
    override the format list entirely, or ALTERSEEK_BZ_EXTRA_FORMATS (e.g.
    "pdf") to add formats on top of the default. ``extra_formats`` adds
    formats for one call without modifying process-wide environment state."""
    root, ext = os.path.splitext(output_path)
    default_fmt = ext[1:] if ext else 'png'
    raw_formats = os.environ.get('ALTERSEEK_BZ_FORMATS', default_fmt)
    raw_formats += ',' + os.environ.get('ALTERSEEK_BZ_EXTRA_FORMATS', '')
    if isinstance(extra_formats, str):
        extra_formats = [extra_formats]
    raw_formats += ',' + ','.join(extra_formats or [])
    formats = []
    for item in raw_formats.replace(';', ',').split(','):
        fmt = item.strip().lower().lstrip('.')
        if fmt and fmt not in formats:
            formats.append(fmt)
    if not formats:
        formats = [default_fmt]
    return [f"{root}.{fmt}" for fmt in formats]


def _save_figure(fig, output_path, extra_formats=None, **kwargs):
    saved_paths = _figure_output_paths(output_path, extra_formats=extra_formats)
    for path in saved_paths:
        fig.savefig(path, **kwargs)
    return saved_paths


def _print_saved_paths(saved_paths, verbose=True, view=None):
    if not verbose:
        return
    suffix = ""
    if view is not None:
        suffix = f"  (view_elev = {view[0]:.2f}, view_azim = {view[1]:.2f})"
    for path in saved_paths:
        print(f"Saved: {path}{suffix}")


CAMERA_ANGLE_FILENAME = "figure_camera_angle.txt"


def write_camera_angle_file(entries):
    """Record the closing camera angle of each 3D figure beside the figures."""
    if not entries:
        return None
    output_path = os.path.join(
        os.path.dirname(entries[0][0]), CAMERA_ANGLE_FILENAME
    )
    lines = [
        "# Camera angle captured when each interactive 3D window was closed.",
        "# Paste one pair into alterseek_input.toml to reopen at that angle.",
        "",
    ]
    for figure_path, elev, azim in entries:
        lines.append(os.path.basename(figure_path))
        lines.append(f"view_elev = {elev:.2f}")
        lines.append(f"view_azim = {azim:.2f}")
        lines.append("")
    with _atomic_open_text(output_path) as handle:
        handle.write("\n".join(lines))
    return output_path


_GREEK_MATH = {
    "DELTA": r"\Delta",
    "LAMBDA": r"\Lambda",
    "SIGMA": r"\Sigma",
}


def label_aliases(label):
    """Return the individual names in a slash-combined point label."""
    return tuple(part for part in str(label).split("/") if part)


def combine_point_labels(*labels):
    """Combine point names with slashes while preserving their first occurrence."""
    combined = []
    for label in labels:
        for alias in label_aliases(label):
            if alias not in combined:
                combined.append(alias)
    return "/".join(combined)


def prime_point_label(label):
    """Prime every non-Gamma name in a possibly slash-combined label."""
    primed = []
    for alias in label_aliases(label):
        if alias == GAMMA_LABEL or alias.strip().upper() == "GAMMA":
            primed.append(alias)
        else:
            primed.append(f"{alias}'")
    return combine_point_labels(*primed)


def point_label_is_primed(label):
    """Return whether every name in a combined label is primed."""
    aliases = label_aliases(label)
    return bool(aliases) and all(alias.endswith("'") for alias in aliases)


def unprime_point_label(label):
    """Remove primes independently from every name in a combined label."""
    return combine_point_labels(*(alias.rstrip("'") for alias in label_aliases(label)))


def generated_plain_path_segments(path_sequence):
    """Yield generated high-symmetry segments with their spin side."""
    if path_sequence is None:
        return
    for first, second in zip(path_sequence, path_sequence[1:]):
        if first is None or second is None:
            continue
        first_label, second_label = first[3], second[3]
        if first_label in ('k', "k'") or second_label in ('k', "k'"):
            continue
        first_prime = point_label_is_primed(first_label)
        second_prime = point_label_is_primed(second_label)
        first_base = unprime_point_label(first_label)
        second_base = unprime_point_label(second_label)
        gamma_path = any(
            alias == GAMMA_LABEL or alias.strip().upper() == 'GAMMA'
            for label in (first_base, second_base)
            for alias in label_aliases(label)
        )
        spin_side = 'down' if (
            (first_prime and second_prime)
            or (gamma_path and (first_prime or second_prime))
        ) else 'up'
        yield first, second, spin_side


def grouped_point_labels(points, path_labels, atol=POINT_COINCIDENCE_ATOL):
    """Merge high-symmetry points that land on the same spot into one label.

    A, B and C at the same coordinates are drawn as "A/B/C". Names missing
    from the plotted path are dropped, so it shows "A/B" when C is not in
    the path. If none of the names is in the path, the first is kept anyway
    so the vertex still gets a label. Names starting with "_" are skipped.
    """
    preferred = {
        alias
        for label in path_labels
        for alias in label_aliases(label)
    }
    groups = []
    for label, coords in points.items():
        if str(label).startswith("_"):
            continue
        coords = np.asarray(coords, dtype=float)
        for group in groups:
            if np.allclose(coords, group["coords"], atol=atol, rtol=0.0):
                group["labels"].append(str(label))
                break
        else:
            groups.append({"coords": coords, "labels": [str(label)]})

    result = []
    for group in groups:
        selected = [label for label in group["labels"] if label in preferred]
        if not selected:
            selected = group["labels"][:1]
        result.append((combine_point_labels(*selected), group["coords"]))
    return result


def _math_symbol(label):
    """Return one high-symmetry name without the outer mathtext delimiters."""
    prime_count = len(label) - len(label.rstrip("'"))
    base = label[:-prime_count] if prime_count else label
    if base == GAMMA_LABEL or base.upper() == "GAMMA":
        symbol = r"\Gamma"
    elif '_' in base:
        head, sub = base.split('_', 1)
        symbol = rf"{_GREEK_MATH.get(head.upper(), head)}_{{{sub}}}"
    else:
        symbol = _GREEK_MATH.get(base.upper(), base)
    if prime_count:
        prime_marks = " ".join([r"\prime"] * prime_count)
        symbol = rf"{symbol}^{{{prime_marks}}}"
    return symbol


def _math_label(label):
    """Return a bold mathtext label for high-symmetry point names."""
    symbol = "/".join(_math_symbol(alias) for alias in label_aliases(label))
    return rf"$\mathbf{{{symbol}}}$"
