"""Load HPKOT k-points and paths from SeeK-path.

Coordinates are fractional coefficients in SeeK-path's crystallographic
primitive reciprocal basis. This module normalizes labels and adds the
project-specific data needed to construct the selected IBZ hull.
"""

from __future__ import annotations

import math

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


# HPKOT caption-only and additional path labels stay in the point tables so optional segments can be generated, but they are not vertices of the selected IBZ hull used for centroid construction.
HULL_EXCLUDED_POINTS = {
    "cP1": {"X_1"},
    "cP2": {"X_1"},
    "cF1": {"W_2"},
    "cF2": {"W_2"},
    "hP1": {"H_2"},
    "hP2": {"H_2"},
    # SeeK-path includes symmetry-equivalent copies across the full BZ, so keep only the IBZ-wedge vertices that occur in the selected path.
    "hR1": {"L_4", "F_2", "S_4", "S_6", "H_4", "H_6",
            "M_8", "M_6"}, 
    "hR2": {"M", "M_2"},
    # In oC2, HPKOT _2 labels are mirror copies at kz < 0 in the bottom half of the BZ, outside the selected IBZ wedge, and must not enter its hull.
    "oC2": {"B_2", "G_2", "R_2", "T_2", "Z_2"},
}


# Cubic 23/m-3 and trigonal 3/-3 use no project-only copied vertices.
PROJECT_HULL_EXTRA_POINTS_BY_SG = {
    # Tetragonal 4, -4, and 4/m point groups: doubled project IBZ.
    ("tP1", range(75, 89)): {
        "X_A": ("1/2", "0", "0"),
        "R_A": ("1/2", "0", "1/2"),
    },
    ("tI1", range(75, 89)): {
        "X_A": ("-1/2", "1/2", "0"),
        "P_A": ("-1/4", "3/4", "-1/4"),
    },
    ("tI2", range(75, 89)): {
        "N_A": ("1/2", "0", "0"),
        "S_0A": ("H", "-H", "H"),
        "S_A": ("1-H", "H", "-H"),
        "R_A": ("Z", "-Z", "1/2"),
    },

    # The trigonal hP1 second half is the in-plane image across the Gamma-M plane, which is the only in-plane doubling the group allows here: M and L stay put and K, H are copied.
    # H_2 of the SeeK-path k-vector table sits at kz = -1/2, outside this IBZ, and is symmetry-equivalent to the copied H_A.
    ("hP1", frozenset({149, 151, 153, 157, 159, 162, 163})): {
        "K_A": ("2/3", "-1/3", "0"),
        "H_A": ("2/3", "-1/3", "1/2"),
    },

    # Trigonal 32, 3m, -3m and hexagonal 6, -6, 6/m share the same in-plane side copied sector across the Gamma-K plane.
    ("hP2", range(149, 177)): {
        "M_A": ("0", "1/2", "0"),
        "L_A": ("0", "1/2", "1/2"),
    },
}


# Reference paths for doubled IBZs; cubic 23/m-3 and trigonal 3/-3 use the ordinary path.
PROJECT_HULL_PATH_BY_SG = {
    ("tP1", range(75, 89)): [
        ("\u0393", "X"), ("X", "M"), ("M", "\u0393"),
        ("\u0393", "Z"),
        ("Z", "R"), ("R", "A"), ("A", "Z"),
        ("X", "R"), ("M", "A"),
    ],
    ("tI1", range(75, 89)): [
        ("\u0393", "X"), ("X", "M"), ("M", "\u0393"), ("\u0393", "Z"),
        ("Z_0", "M"), ("X", "P"), ("P", "N"), ("N", "\u0393"),
        ("N", "P_A"),
    ],
    ("tI2", range(75, 89)): [
        ("\u0393", "X"), ("X", "P"), ("P", "N"), ("N", "\u0393"),
        ("\u0393", "M"), ("M", "S"),
        ("S_0", "\u0393"), ("X", "R"), ("G", "M"),
        ("P", "N_A"),
    ],
    ("hP1", frozenset({149, 151, 153, 157, 159, 162, 163})): [
        ("\u0393", "M"), ("M", "K"), ("K", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "A"),
        ("L", "M"), ("H", "K"),
        ("H_A", "K_A"),
    ],
    ("hP2", range(149, 177)): [
        ("\u0393", "M"), ("M", "K"), ("K", "\u0393"),
        ("\u0393", "A"),
        ("A", "L"), ("L", "H"), ("H", "A"),
        ("L", "M"), ("H", "K"),
    ],
}


def _normalize_label(label: str) -> str:
    """Convert SeeK-path table labels to the local display-preserving labels."""
    return GREEK_LABELS.get(label, label)


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


def _load_hpkot_table(ext_bravais: str):
    """Load and normalize one SeeK-path HPKOT table."""
    kparam_def, points_def, path = get_path_data(ext_bravais)
    points_def = {
        _normalize_label(label): tuple(exprs)
        for label, exprs in points_def.items()
    }
    path = _normalize_path(path)

    # SeeK-path includes segments specified only for particular space groups in the HPKOT table captions without applying those conditions.
    # Keep their point definitions, but let EXTRA_PATH_RULES decide when to include the segments.
    if ext_bravais in {"cP1", "cP2"}:
        path = _strip_path_segments(path, [("M", "X_1")])
    elif ext_bravais in {"cF1", "cF2"}:
        path = _strip_path_segments(path, [("X", "W_2")])
    elif ext_bravais in {"hP1", "hP2"}:
        path = _strip_path_segments(path, [("K", "H_2")])

    return {
        "kparam_def": kparam_def,
        "points_def": points_def,
        "kpath": path,
        "display_labels": {label: _format_display_label(label) for label in points_def},
        "hull_excluded_points": set(HULL_EXCLUDED_POINTS.get(ext_bravais, set())),
    }


def _load_lattice_data():
    data = {key: _load_hpkot_table(key) for key in HPKOT_LATTICE_TYPES}

    # Project-curated mC1 closure vertices define the selected IBZ hull but are not public HPKOT labels, so exclude them from display labels and k-paths.
    data["mC1"]["hidden_points_def"] = {
        "_F4": ("-1/2+Z+S", "-1/2-Z+S", "1-H"),
        "_Q1": ("-1/2-Z+S", "-1/2+Z+S", "H"),
        "_Q2": ("1/2+Z-S", "3/2-Z-S", "1-H"),
        "_P1": ("-3/2+Z+S", "1/2-Z+S", "1-H"),
    }

    return data


LATTICE_DATA = _load_lattice_data()


def canonical_lattice_type(lattice_type: str) -> str:
    """Return a validated HPKOT extended lattice key."""
    if lattice_type not in LATTICE_DATA:
        raise KeyError(lattice_type)
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
            # Hidden closure vertices may use small arithmetic combinations of HPKOT parameters, such as -1/2+Z+S.
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
    """Return evaluated HPKOT and project-specific k-points in the primitive reciprocal basis.

    For monoclinic lattices, ``alpha`` supplies the beta angle when ``beta`` is omitted.
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
    """Return k-points used for the AlterSeeK-Path IBZ hull and centroid."""
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
    """Return the IBZ hull/display path, distinct from the HPKOT band path."""
    key = canonical_lattice_type(lattice_type)
    if spacegroup_number is not None:
        sg = int(spacegroup_number)
        for (rule_key, sg_range), path in PROJECT_HULL_PATH_BY_SG.items():
            if key == rule_key and sg in sg_range:
                return list(path)
    return list(LATTICE_DATA[key]["kpath"])


def get_kpath(lattice_type, spacegroup_number=None):
    """Return the HPKOT path with any space-group-specific extra segments."""
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
