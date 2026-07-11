"""Shared figure I/O + label/style helpers.

Extracted from compute_centroid_hybrid.py (restructuring phase 1).
"""
import os
from .lattice_kpoints import canonical_lattice_type

GAMMA_LABEL = "\u0393"
BZ_SPECIAL_COLORS = {
    "orange": "#e68613",
    "purple": "#6b5596",
}
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


def _figure_output_paths(output_path):
    """Return the requested figure output paths, preserving the requested
    extension as the default. Set ALTERSEEK_BZ_FORMATS (e.g. "png,pdf") to
    override the format list entirely, or ALTERSEEK_BZ_EXTRA_FORMATS (e.g.
    "pdf") to add formats on top of the default -- the latter is what the
    toml `save_pdf` option sets."""
    root, ext = os.path.splitext(output_path)
    default_fmt = ext[1:] if ext else 'png'
    raw_formats = os.environ.get('ALTERSEEK_BZ_FORMATS', default_fmt)
    raw_formats += ',' + os.environ.get('ALTERSEEK_BZ_EXTRA_FORMATS', '')
    formats = []
    for item in raw_formats.replace(';', ',').split(','):
        fmt = item.strip().lower().lstrip('.')
        if fmt and fmt not in formats:
            formats.append(fmt)
    if not formats:
        formats = [default_fmt]
    return [f"{root}.{fmt}" for fmt in formats]


def _save_figure(fig, output_path, **kwargs):
    saved_paths = _figure_output_paths(output_path)
    for path in saved_paths:
        fig.savefig(path, **kwargs)
    return saved_paths


def _print_saved_paths(saved_paths, verbose=True):
    if not verbose:
        return
    for path in saved_paths:
        print(f"Saved: {path}")


def _math_label(label):
    """Return a bold mathtext label for high-symmetry point names."""
    label = str(label)
    prime = label.endswith("'")
    base = label.rstrip("'")
    if base == '\u0393' or base.upper() == "GAMMA":
        symbol = r"\Gamma"
    elif '_' in base:
        head, sub = base.split('_', 1)
        greek = {
            "DELTA": r"\Delta",
            "LAMBDA": r"\Lambda",
            "SIGMA": r"\Sigma",
        }
        symbol = rf"{greek.get(head.upper(), head)}_{{{sub}}}"
    else:
        greek = {
            "DELTA": r"\Delta",
            "LAMBDA": r"\Lambda",
            "SIGMA": r"\Sigma",
        }
        symbol = greek.get(base.upper(), base)
    if prime:
        # Attach the prime as a real mathtext superscript (inside the same
        # braces as any subscript) so it stacks tightly next to the base
        # symbol, matching e.g. 6_{001}^{+}, instead of floating off as a
        # trailing plain-text character.
        symbol = rf"{symbol}^{{\prime}}"
    return rf"$\mathbf{{{symbol}}}$"
