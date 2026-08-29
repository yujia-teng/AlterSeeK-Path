"""MCIF-specific moment validation and nonmagnetic-parent metadata.

MAGNDATA MCIF coordinates are commonly rounded. The declared parent
transformation and space group allow a looser spglib tolerance when it
reproduces that parent.
"""
import warnings

import numpy as np


# This is a conservative ladder because a looser tolerance is accepted only when it reproduces the structure's own declared parent, preventing a genuinely lower-symmetry structure from being over-symmetrized.
_MCIF_PARENT_SYMPREC_CANDIDATES = (1e-5, 1e-4, 1e-3)
_MCIF_COLLINEAR_MOMENT_TOLERANCE = 0.02


def _validate_collinear_moments(
    moments,
    moment_tolerance=_MCIF_COLLINEAR_MOMENT_TOLERANCE,
):
    """Reject MCIF moments that are not collinear.

    Zero moments are ignored, and parallel and antiparallel moments are treated
    as collinear. The tolerance accommodates rounded moment components.
    """
    moments = np.asarray(moments, dtype=float)
    norms = np.linalg.norm(moments, axis=1)
    nonzero = moments[norms > 1e-10]
    if len(nonzero) < 2:
        return

    axis = nonzero[0] / np.linalg.norm(nonzero[0])
    transverse = np.linalg.norm(
        np.cross(nonzero, axis),
        axis=1,
    )
    if np.any(transverse > moment_tolerance):
        raise ValueError(
            "Only collinear magnetic structures are supported; noncollinear "
            "moment directions were detected in the MCIF."
        )


def _cif_scalar(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _parent_hint_from_cif_block(block):
    """Return the declared nonmagnetic parent-cell index and SG, if present."""
    transform = _cif_scalar(block.get("_parent_space_group.child_transform_Pp_abc"))
    if not transform:
        return None
    try:
        from pymatgen.symmetry.settings import JonesFaithfulTransformation

        parsed = JonesFaithfulTransformation.from_transformation_str(str(transform))
        index = int(round(abs(float(np.linalg.det(np.asarray(parsed.P, dtype=float))))))
    except Exception:
        return None
    if index <= 1:
        return None

    parent_number = _cif_scalar(block.get("_parent_space_group.IT_number"))
    try:
        parent_number = int(parent_number)
    except (TypeError, ValueError):
        parent_number = None
    parent_symbol = _cif_scalar(block.get("_parent_space_group.name_H-M_alt"))
    return {
        "index": index,
        "spacegroup_number": parent_number,
        "spacegroup_symbol": str(parent_symbol) if parent_symbol else None,
        "transform": str(transform),
    }


def _declared_mcif_parent_hint(filename):
    if not str(filename).lower().endswith(".mcif"):
        return None
    try:
        from pymatgen.io.cif import CifParser

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            blocks = CifParser(filename).as_dict()
        if not blocks:
            return None
        return _parent_hint_from_cif_block(next(iter(blocks.values())))
    except Exception:
        return None
