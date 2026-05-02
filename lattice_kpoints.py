"""
Curated high-symmetry k-point data in the SeeK-path/HPKOT convention.

Primary source:

  Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka,
  Computational Materials Science 128, 140-184 (2017), Tables 69-92.

Coordinates returned by this module are kP coordinates: fractional
coefficients in the crystallographic primitive reciprocal basis used by
SeeK-path.  Gamma is stored internally as the Greek label ``Γ``.  Other
Greek and subscripted labels are preserved semantically, e.g. ``H_2`` is
displayed as ``$H_2$``.

The public HPKOT band path and the project-curated IBZ hull are related
but distinct objects.  This module stores all HPKOT points and can add
project-only hidden closure vertices (labels beginning with ``_``) for
centroid/hull construction where the visible band-path labels do not close
the selected irreducible domain.
"""

from __future__ import annotations

import math
from copy import deepcopy

import numpy as np

from seekpath.hpkot.tools import (
    eval_expr,
    eval_expr_simple,
    extend_kparam,
    get_path_data,
)


# SeeK-path/HPKOT extended Bravais symbols with bundled table data.
HPKOT_LATTICE_TYPES = [
    "aP2", "aP3",
    "cF1", "cF2", "cI1", "cP1", "cP2",
    "hP1", "hP2", "hR1", "hR2",
    "mC1", "mC2", "mC3", "mP1",
    "oA1", "oA2", "oC1", "oC2", "oF1", "oF2", "oF3",
    "oI1", "oI2", "oI3", "oP1",
    "tI1", "tI2", "tP1",
]


GREEK_LABELS = {
    "GAMMA": "\u0393",
    "DELTA": "\u0394",
    "LAMBDA": "\u039b",
    "SIGMA": "\u03a3",
}

GREEK_INTERNAL_LABELS = set(GREEK_LABELS.values())


EXTRA_PATH_RULES = {
    # Hinuma et al. Fig. 2 / Table 69 caption.
    "cP1": [
        {"spacegroups": {195, 198, 200, 201, 205}, "segments": [("M", "X_1")]},
    ],
    "cP2": [
        {"spacegroups": {195, 198, 200, 201, 205}, "segments": [("M", "X_1")]},
    ],
    # Hinuma et al. Fig. 3 / Table 70 caption.
    "cF1": [
        {"spacegroups": {196, 202, 203}, "segments": [("X", "W_2")]},
    ],
    "cF2": [
        {"spacegroups": {196, 202, 203}, "segments": [("X", "W_2")]},
    ],
    # Hinuma et al. Fig. 17 / Table 85 caption.
    "hP1": [
        {
            "spacegroups": {143, 144, 145, 147, 149, 151, 153, 157, 159, 162, 163},
            "segments": [("K", "H_2")],
        },
    ],
    "hP2": [
        {
            "spacegroups": {143, 144, 145, 147, 149, 151, 153, 157, 159, 162, 163},
            "segments": [("K", "H_2")],
        },
    ],
}


# HPKOT caption-only/additional path labels are stored in the point tables so
# optional segments can be generated, but they are not vertices of the selected
# IBZ hull used for centroid construction.
HULL_EXCLUDED_POINTS = {
    "cP1": {"X_1"},
    "cP2": {"X_1"},
    "cF1": {"W_2"},
    "cF2": {"W_2"},
    "hP1": {"H_2"},
    "hP2": {"H_2"},
    # seekpath includes all symmetry-equivalent copies across the full BZ.
    # Only the IBZ-wedge vertices (those that appear in the kpath) are kept.
    "hR1": {"L_4", "F_2", "S_4", "S_6", "H_4", "H_6",
            "M_8", "M_6"}, 
    "hR2": {"M", "M_2"}
}


PROJECT_HULL_EXTRA_POINTS_BY_SG = {
    # Cubic 23 and m-3 point groups: doubled project IBZ relative to m-3m.
    ("cP1", range(195, 207)): {
        "X_A": ("1/2", "0", "0"),
    },
    ("cF1", range(195, 207)): {
        "U_A": ("1/4", "5/8", "5/8"),
        "W_A": ("1/4", "1/2", "3/4"),
        "X_A": ("0", "1/2", "1/2"),
    },
    ("cI1", range(195, 207)): {
        "H_A": ("-1/2", "1/2", "1/2"),
    },

    # Tetragonal 4, -4, and 4/m point groups: doubled project IBZ.
    ("tP1", range(75, 89)): {
        "X_A": ("1/2", "0", "0"),
        "R_A": ("1/2", "0", "1/2"),
    },
    ("tI1", range(75, 89)): {
        "M_A": ("1/2", "-1/2", "1/2"),
        "N_A": ("1/2", "0", "0"),
        "Z_0A": ("1-H", "-H", "H"),
    },
    ("tI2", range(75, 89)): {
        "N_A": ("1/2", "0", "0"),
        "S_0A": ("H", "-H", "H"),
        "S_A": ("1-H", "H", "-H"),
        "R_A": ("Z", "-Z", "1/2"),
    },

    # Trigonal 3 and -3 point groups on a primitive hexagonal lattice:
    # quadrupled 120-degree project IBZ relative to hP holohedry.
    ("hP1", range(143, 149)): {
        "M_A": ("0", "1/2", "0"),
        "L_A": ("0", "1/2", "1/2"),
        "K_A": ("-1/3", "2/3", "0"),
        "H_A": ("-1/3", "2/3", "1/2"),
        "M_B": ("-1/2", "1/2", "0"),
        "L_B": ("-1/2", "1/2", "1/2"),
    },

    # Trigonal 32, 3m, -3m and hexagonal 6, -6, 6/m point groups:
    # doubled 60-degree project IBZ relative to the highest-symmetry hP wedge.
    ("hP1", range(149, 168)): {
        "M_A": ("0", "1/2", "0"),
        "L_A": ("0", "1/2", "1/2"),
    },
    ("hP2", range(149, 168)): {
        "M_A": ("0", "1/2", "0"),
        "L_A": ("0", "1/2", "1/2"),
    },
    ("hP2", range(168, 177)): {
        "M_A": ("0", "1/2", "0"),
        "L_A": ("0", "1/2", "1/2"),
    },
}


PROJECT_HULL_PATH_BY_SG = {
    ("cP1", range(195, 207)): [
        ("\u0393", "X"), ("X", "M"), ("M", "X_A"), ("X_A", "\u0393"),
        ("\u0393", "R"), ("R", "X"), ("R", "M"), ("R", "X_A"),
    ],
    ("cF1", range(195, 207)): [
        ("\u0393", "X"), ("X", "W"), ("W", "K"), ("K", "\u0393"),
        ("\u0393", "L"), ("L", "U"), ("U", "W"), ("W", "L"),
        ("L", "K"), ("U", "X"),
        ("\u0393", "X_A"), ("X_A", "W_A"), ("W_A", "K"),
        ("L", "U_A"), ("U_A", "W_A"), ("U_A", "X_A"),
    ],
    ("cI1", range(195, 207)): [
        ("\u0393", "H"), ("H", "N"), ("N", "\u0393"), ("\u0393", "P"),
        ("P", "H"), ("P", "N"),
        ("\u0393", "H_A"), ("H_A", "N"), ("P", "H_A"),
    ],
    ("tP1", range(75, 89)): [
        ("\u0393", "X"), ("X", "M"), ("M", "X_A"), ("X_A", "\u0393"),
        ("\u0393", "Z"),
        ("Z", "R"), ("R", "A"), ("A", "R_A"), ("R_A", "Z"),
        ("X", "R"), ("M", "A"), ("X_A", "R_A"),
    ],
    ("tI1", range(75, 89)): [
        ("\u0393", "X"), ("X", "M"), ("M", "\u0393"), ("\u0393", "Z"),
        ("Z", "P"), ("P", "N"), ("N", "Z_0"), ("Z_0", "M"), ("X", "P"),
        ("M_A", "\u0393"), ("P", "N_A"), ("N_A", "Z_0A"), ("Z_0A", "M_A"),
    ],
    ("tI2", range(75, 89)): [
        ("\u0393", "X"), ("X", "R"), ("R", "S_0"), ("S_0", "\u0393"),
        ("\u0393", "M"), ("M", "S"), ("S", "N"), ("N", "P"),
        ("P", "G"), ("G", "M"), ("X", "P"),
        ("X", "R_A"), ("R_A", "S_0A"), ("S_0A", "\u0393"),
        ("M", "S_A"), ("S_A", "N_A"), ("N_A", "P"),
    ],
    ("hP1", range(143, 149)): [
        ("\u0393", "M"), ("M", "K"), ("K", "M_A"), ("M_A", "K_A"),
        ("K_A", "M_B"), ("M_B", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "L_A"), ("L_A", "H_A"),
        ("H_A", "L_B"), ("L_B", "A"),
        ("L", "M"), ("H", "K"), ("L_A", "M_A"), ("H_A", "K_A"),
        ("L_B", "M_B"),
    ],
    ("hP1", range(149, 168)): [
        ("\u0393", "M"), ("M", "K"), ("K", "M_A"), ("M_A", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "L_A"), ("L_A", "A"),
        ("L", "M"), ("H", "K"), ("M_A", "L_A"),
    ],
    ("hP2", range(149, 168)): [
        ("\u0393", "M"), ("M", "K"), ("K", "M_A"), ("M_A", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "L_A"), ("L_A", "A"),
        ("L", "M"), ("H", "K"), ("M_A", "L_A"),
    ],
    ("hP2", range(168, 177)): [
        ("\u0393", "M"), ("M", "K"), ("K", "M_A"), ("M_A", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "L_A"), ("L_A", "A"),
        ("L", "M"), ("H", "K"), ("M_A", "L_A"),
    ],
}


def _normalize_label(label: str) -> str:
    """Convert SeeK-path table labels to the local display-preserving labels."""
    return GREEK_LABELS.get(label, label)


def _denormalize_label(label: str) -> str:
    """Convert local labels back to SeeK-path table labels."""
    if label == "\u0393":
        return "GAMMA"
    for raw, normalized in GREEK_LABELS.items():
        if label == normalized:
            return raw
    return label


def _normalize_path(path):
    return [(_normalize_label(a), _normalize_label(b)) for a, b in path]


def _strip_path_segments(path, segments_to_remove):
    remove = {tuple(seg) for seg in segments_to_remove}
    return [seg for seg in path if tuple(seg) not in remove]


def _format_display_label(label: str) -> str:
    if label.startswith("_"):
        label = label[1:]
    if label in GREEK_INTERNAL_LABELS:
        return rf"${label}$"
    if "_" in label:
        base, sub = label.split("_", 1)
        base = GREEK_LABELS.get(base, base)
        return rf"${base}_{{{sub.replace('_', '')}}}$"
    return label


def _hpkot_table(ext_bravais: str):
    kparam_def, points_def, path = get_path_data(ext_bravais)
    points_def = {
        _normalize_label(label): tuple(exprs)
        for label, exprs in points_def.items()
    }
    path = _normalize_path(path)

    # The SeeK-path bundled cP1/cF1/hP1 path files contain the optional
    # HPKOT caption segments unconditionally.  Keep them as points, but add
    # the segments only when the paper's space-group list applies.
    if ext_bravais in {"cP1", "cP2"}:
        path = _strip_path_segments(path, [("M", "X_1")])
    elif ext_bravais in {"cF1", "cF2"}:
        path = _strip_path_segments(path, [("X", "W_2")])
    elif ext_bravais in {"hP1", "hP2"}:
        path = _strip_path_segments(path, [("K", "H_2")])

    return {
        "source": "HPKOT",
        "kparam_def": kparam_def,
        "points_def": points_def,
        "kpath": path,
        "display_labels": {label: _format_display_label(label) for label in points_def},
        "hull_excluded_points": set(HULL_EXCLUDED_POINTS.get(ext_bravais, set())),
    }


def _build_lattice_data():
    data = {key: _hpkot_table(key) for key in HPKOT_LATTICE_TYPES}

    # Project-curated hidden closure vertices for the mC1 selected IBZ hull.
    # These are not public HPKOT path labels, so they are excluded from the
    # display label map and are never added to k-paths.
    data["mC1"]["hidden_points_def"] = {
        "_F4": ("-1/2+Z+S", "-1/2-Z+S", "1-H"),
        "_Q1": ("-1/2-Z+S", "-1/2+Z+S", "H"),
        "_Q2": ("1/2+Z-S", "3/2-Z-S", "1-H"),
        "_P1": ("-3/2+Z+S", "1/2-Z+S", "1-H"),
    }

    # Backward-compatible aliases for older callers and scripts.  New code
    # should use HPKOT extended symbols directly.
    aliases = {
        "CUB": "cP2", "CUB2": "cP1",
        "FCC": "cF2", "FCC2": "cF1",
        "BCC": "cI1", "BCC2": "cI1",
        "TET": "tP1", "TET2": "tP1",
        "BCT1": "tI1", "BCT1_2": "tI1",
        "BCT2": "tI2", "BCT2_2": "tI2",
        "ORC": "oP1",
        "ORCF1": "oF1", "ORCF2": "oF2", "ORCF3": "oF3",
        "ORCI": "oI1",
        "ORCC1": "oC1", "ORCC2": "oC2", "ORCC": "oC2",
        "HEX": "hP2", "HEX2": "hP1", "HEX4": "hP1",
        "RHL1": "hR1", "RHL1_2": "hR1",
        "RHL2": "hR2", "RHL2_2": "hR2",
        "MCL": "mP1",
        "MCLC1": "mC1", "MCLC2_SC": "mC1",
        "MCLC2": "mC2", "MCLC4_SC": "mC2", "MCLC4": "mC2",
        "MCLC3": "mC3", "MCLC5": "mC3",
        "TRI1a": "aP2", "TRI1b": "aP2",
        "TRI2a": "aP3", "TRI2b": "aP3",
    }
    for alias, target in aliases.items():
        data[alias] = {"alias_for": target}

    return data


LATTICE_DATA = _build_lattice_data()


def canonical_lattice_type(lattice_type: str) -> str:
    """Return the canonical HPKOT key for a lattice type or old alias."""
    data = LATTICE_DATA[lattice_type]
    while "alias_for" in data:
        lattice_type = data["alias_for"]
        data = LATTICE_DATA[lattice_type]
    return lattice_type


def _cos_from_degrees(angle, default=90.0):
    if angle is None:
        angle = default
    return math.cos(math.radians(float(angle)))


def _evaluate_kparams(data, a, b, c, alpha, beta, gamma):
    a = 1.0 if a is None else float(a)
    b = a if b is None else float(b)
    c = a if c is None else float(c)

    cosalpha = _cos_from_degrees(alpha, 90.0)
    cosbeta = _cos_from_degrees(beta if beta is not None else alpha, 90.0)
    cosgamma = _cos_from_degrees(gamma, 90.0)

    kparam = {}
    for name, expr in data.get("kparam_def", []):
        kparam[name] = eval_expr(expr, a, b, c, cosalpha, cosbeta, cosgamma, kparam)
    return extend_kparam(kparam)


def _evaluate_point_exprs(point_exprs, kparam):
    values = []
    safe_globals = {"__builtins__": {}}
    safe_locals = dict(kparam)
    for expr in point_exprs:
        try:
            values.append(eval_expr_simple(expr, kparam))
        except ValueError:
            # Project-only hidden closure vertices can use small arithmetic
            # combinations of HPKOT parameters, e.g. "-1/2+Z+S".
            values.append(float(eval(expr, safe_globals, safe_locals)))
    return values


def get_kpoints(
    lattice_type,
    a=None,
    b=None,
    c=None,
    alpha=None,
    beta=None,
    gamma=None,
    include_hidden=True,
):
    """
    Return HPKOT/project k-points in primitive reciprocal coordinates kP.

    Parameters keep the historical call shape used by this project.  For
    monoclinic cells, ``alpha`` is treated as the monoclinic beta angle when
    ``beta`` is not provided.
    """
    key = canonical_lattice_type(lattice_type)
    data = LATTICE_DATA[key]
    kparam = _evaluate_kparams(data, a, b, c, alpha, beta, gamma)

    points = {
        label: _evaluate_point_exprs(exprs, kparam)
        for label, exprs in data["points_def"].items()
    }
    if include_hidden:
        points.update({
            label: _evaluate_point_exprs(exprs, kparam)
            for label, exprs in data.get("hidden_points_def", {}).items()
        })
    return points


def get_hull_kpoints(
    lattice_type,
    a=None,
    b=None,
    c=None,
    alpha=None,
    beta=None,
    gamma=None,
    include_hidden=True,
    spacegroup_number=None,
):
    """Return k-points used for the project IBZ hull/centroid."""
    key = canonical_lattice_type(lattice_type)
    points = get_kpoints(
        key, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
        include_hidden=include_hidden,
    )
    for label in LATTICE_DATA[key].get("hull_excluded_points", set()):
        points.pop(label, None)
    if spacegroup_number is not None:
        sg = int(spacegroup_number)
        kparam = _evaluate_kparams(
            LATTICE_DATA[key], a, b, c, alpha, beta, gamma)
        for (rule_key, sg_range), extra_points in PROJECT_HULL_EXTRA_POINTS_BY_SG.items():
            if key == rule_key and sg in sg_range:
                points.update({
                    label: _evaluate_point_exprs(exprs, kparam)
                    for label, exprs in extra_points.items()
                })
    return points


def get_hull_kpath(lattice_type, spacegroup_number=None):
    """Return the project hull/display path, distinct from the HPKOT band path."""
    key = canonical_lattice_type(lattice_type)
    if spacegroup_number is not None:
        sg = int(spacegroup_number)
        for (rule_key, sg_range), path in PROJECT_HULL_PATH_BY_SG.items():
            if key == rule_key and sg in sg_range:
                return list(path)
    return list(LATTICE_DATA[key]["kpath"])


def get_path_kpoints(
    lattice_type,
    path_segments,
    a=None,
    b=None,
    c=None,
    alpha=None,
    beta=None,
    gamma=None,
):
    """Return just the point coordinates required by a band path."""
    points = get_kpoints(
        lattice_type, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
        include_hidden=False,
    )
    labels = {label for segment in path_segments for label in segment}
    return {label: points[label] for label in labels if label in points}


def get_kpath(lattice_type, spacegroup_number=None, with_time_reversal=True):
    """
    Return the HPKOT base path plus applicable paper-defined extra segments.

    ``with_time_reversal`` is accepted to keep the API explicit; path doubling
    for no-time-reversal workflows is handled by the caller because it also
    needs operation-specific primed coordinates.
    """
    key = canonical_lattice_type(lattice_type)
    path = list(LATTICE_DATA[key]["kpath"])
    if spacegroup_number is not None:
        sg = int(spacegroup_number)
        for rule in EXTRA_PATH_RULES.get(key, []):
            if sg in rule["spacegroups"]:
                path.extend(rule["segments"])
    return path


def get_display_labels(lattice_type, include_hidden=False):
    """Return label -> Matplotlib/LaTeX display text."""
    key = canonical_lattice_type(lattice_type)
    labels = dict(LATTICE_DATA[key]["display_labels"])
    if include_hidden:
        labels.update({
            label: _format_display_label(label)
            for label in LATTICE_DATA[key].get("hidden_points_def", {})
        })
    for extra_points in PROJECT_HULL_EXTRA_POINTS_BY_SG.values():
        labels.update({
            label: _format_display_label(label)
            for label in extra_points
        })
    return labels


def get_params(
    lattice_type,
    a=None,
    b=None,
    c=None,
    alpha=None,
    beta=None,
    gamma=None,
):
    """Return evaluated HPKOT k-vector parameters for a parametric table."""
    key = canonical_lattice_type(lattice_type)
    data = LATTICE_DATA[key]
    if not data.get("kparam_def"):
        return {}
    raw = _evaluate_kparams(data, a, b, c, alpha, beta, gamma)
    return {name: raw[name] for name, _expr in data["kparam_def"]}


def seekpath_label_to_internal(label: str) -> str:
    """Normalize a SeeK-path label for lookup in this module."""
    return _normalize_label(label)


def internal_label_to_seekpath(label: str) -> str:
    """Convert a local label to the corresponding SeeK-path label."""
    return _denormalize_label(label)


# Historical detection helper retained for scripts that import it directly.
def get_bravais_type(
    spacegroup_number,
    conv_a,
    conv_b,
    conv_c,
    conv_alpha=90.0,
    conv_beta=90.0,
    conv_gamma=90.0,
    centering="P",
):
    """Best-effort HPKOT extended key from space group and cell parameters."""
    sg = int(spacegroup_number)
    centering = str(centering).upper()

    if 195 <= sg <= 230:
        if centering == "F":
            return "cF1" if sg <= 206 else "cF2"
        if centering == "I":
            return "cI1"
        return "cP1" if sg <= 206 else "cP2"
    if 75 <= sg <= 142:
        return "tI1" if centering == "I" and conv_c < conv_a else (
            "tI2" if centering == "I" else "tP1")
    if 143 <= sg <= 194:
        if centering == "R":
            return "hR1" if math.sqrt(3) * conv_a < math.sqrt(2) * conv_c else "hR2"
        return "hP1" if sg <= 163 else "hP2"
    if 16 <= sg <= 74:
        if centering == "F":
            a, b, c = sorted([conv_a, conv_b, conv_c])
            inv_a2, inv_b2, inv_c2 = 1 / a**2, 1 / b**2, 1 / c**2
            if abs(inv_a2 - (inv_b2 + inv_c2)) < 1e-8:
                return "oF3"
            return "oF1" if inv_a2 > inv_b2 + inv_c2 else "oF2"
        if centering == "I":
            largest = max((conv_a, "a"), (conv_b, "b"), (conv_c, "c"))[1]
            return {"a": "oI2", "b": "oI3", "c": "oI1"}[largest]
        if centering in {"C", "A", "B"}:
            return "oC1" if conv_a < conv_b else "oC2"
        return "oP1"
    if 3 <= sg <= 15:
        if centering in {"C", "A", "B"}:
            beta_rad = math.radians(conv_beta)
            if conv_b < conv_a * math.sin(beta_rad):
                return "mC1"
            cond = (
                -conv_a * math.cos(beta_rad) / conv_c
                + conv_a**2 * math.sin(beta_rad) ** 2 / conv_b**2
            )
            return "mC2" if cond < 1.0 else "mC3"
        return "mP1"
    return "aP2"


ALL_LATTICE_TYPES = list(LATTICE_DATA.keys())
CANONICAL_LATTICE_TYPES = HPKOT_LATTICE_TYPES
FIXED_LATTICES = [
    key for key in HPKOT_LATTICE_TYPES
    if not LATTICE_DATA[key].get("kparam_def")
]
PARAMETRIC_LATTICES = [
    key for key in HPKOT_LATTICE_TYPES
    if LATTICE_DATA[key].get("kparam_def")
]


if __name__ == "__main__":
    print("=" * 60)
    print("HPKOT/SeeK-path k-point database")
    print("=" * 60)
    print(f"Canonical lattice types: {len(CANONICAL_LATTICE_TYPES)}")
    for lt in CANONICAL_LATTICE_TYPES:
        kp = get_kpoints(lt, a=5.0, b=6.0, c=7.0, alpha=100.0)
        print(f"{lt:4s}: {len(kp):2d} k-points, {len(get_kpath(lt))} path segments")
