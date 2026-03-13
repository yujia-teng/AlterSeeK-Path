"""
Standardized high-symmetry k-point data for all 14 Bravais lattice types
(and their sub-variations) from:

  Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299-312 (2010)
  Tables 2-21, Appendix A

Coordinates are given as fractions of reciprocal lattice vectors b1, b2, b3.
For parametric lattices, a function returns k-points given conventional cell
parameters (a, b, c, alpha in the paper's convention).

Covers ALL Bravais lattices including triclinic (TRI1a/1b/2a/2b).
"""

import numpy as np

# ============================================================================
# Lattice type detection from space group + lattice parameters
# ============================================================================
def get_bravais_type(spacegroup_number, conv_a, conv_b, conv_c,
                     conv_alpha=90.0, conv_beta=90.0, conv_gamma=90.0,
                     centering='P'):
    """
    Determine the Setyawan-Curtarolo BZ variation from space group info.

    Parameters
    ----------
    spacegroup_number : int
    conv_a, conv_b, conv_c : float
        Conventional cell axial lengths.
    conv_alpha, conv_beta, conv_gamma : float
        Conventional cell interaxial angles in degrees.
    centering : str
        Lattice centering symbol: 'P', 'F', 'I', 'C', 'A', 'R'

    Returns
    -------
    str : BZ variation label (e.g. 'CUB', 'FCC', 'BCT1', 'MCLC3', ...)
    """
    sg = spacegroup_number

    # --- Cubic m-3m (SG 207-230) ---
    if 207 <= sg <= 230:
        if centering == 'F':
            return 'FCC'
        elif centering == 'I':
            return 'BCC'
        else:
            return 'CUB'
    
    # --- Cubic m-3 (SG 195-206) ---
    if 195 <= sg <= 206:
        if centering == 'F':
            return 'FCC2'
        elif centering == 'I':
            return 'BCC2'
        else:
            return 'CUB2'

    # --- Hexagonal 6/mmm (SG 177-194) ---
    if 177 <= sg <= 194:
        return 'HEX'
    
    # --- Hexagonal 6/m (SG 168-176) —doubled IBZ ---
    if 168 <= sg <= 176:
        return 'HEX2'

    # --- Trigonal -3m (SG 149-167) ---
    if 149 <= sg <= 167:
        if centering == 'R':
            # Rhombohedral: need primitive cell angle
            # Convention: use alpha of the rhombohedral primitive cell
            # RHL1: alpha < 90, RHL2: alpha > 90
            # Note: conv_alpha here should be the rhombohedral angle
            if conv_alpha < 90.0:
                return 'RHL1'
            else:
                return 'RHL2'
        else:
            return 'HEX'
    
    # --- Trigonal -3 (SG 143-148) ---
    if 143 <= sg <= 148:
        if centering == 'R':
            # Rhombohedral: need primitive cell angle
            # Convention: use alpha of the rhombohedral primitive cell
            # RHL1: alpha < 90, RHL2: alpha > 90
            # Note: conv_alpha here should be the rhombohedral angle
            if conv_alpha < 90.0:
                return 'RHL1'
            else:
                return 'RHL2'
        else:
            return 'HEX'

    # --- Tetragonal 4/mmm (SG 89-142) ---
    if 89 <= sg <= 142:
        if centering == 'I':
            if conv_c < conv_a:
                return 'BCT1'
            else:
                return 'BCT2'
        else:
            return 'TET'

    # --- Tetragonal 4/m (SG 75-88) —doubled IBZ ---
    if 75 <= sg <= 88:
        if centering == 'I':
            if conv_c < conv_a:
                return 'BCT1'
            else:
                return 'BCT2'
        else:
            return 'TET2'

    # --- Orthorhombic (SG 16-74) ---
    if 16 <= sg <= 74:
        a, b, c = sorted([conv_a, conv_b, conv_c])
        if centering == 'F':
            inv_a2 = 1.0 / a**2
            inv_b2 = 1.0 / b**2
            inv_c2 = 1.0 / c**2
            if abs(inv_a2 - (inv_b2 + inv_c2)) < 1e-8:
                return 'ORCF3'
            elif inv_a2 > inv_b2 + inv_c2:
                return 'ORCF1'
            else:
                return 'ORCF2'
        elif centering == 'I':
            return 'ORCI'
        elif centering in ('C', 'A'):
            return 'ORCC'
        else:
            return 'ORC'

    # --- Monoclinic (SG 3-15) ---
    if 3 <= sg <= 15:
        if centering in ('C', 'A'):
            # C-centered monoclinic (HPKOT/seekpath 3-case convention)
            a, b, c = conv_a, conv_b, conv_c
            alpha = np.radians(conv_alpha)
            if b < a * np.sin(alpha):
                return 'MCLC1'
            cond = -a * np.cos(alpha) / c + a**2 * np.sin(alpha)**2 / b**2
            if cond < 1.0:
                return 'MCLC2'
            return 'MCLC3'
        else:
            return 'MCL'

    # --- Triclinic (SG 1-2) ---
    return 'TRI'


# ============================================================================
# K-point data for each Bravais lattice variation
# ============================================================================
# Each entry is a dict with:
#   'kpoints': dict of label -> [b1, b2, b3] (for fixed-coordinate lattices)
#   'kpoints_func': callable(params) -> dict (for parametric lattices)
#   'params_func': callable(a, b, c, alpha) -> dict of parameter values
#   'kpath': list of (label1, label2) tuples
#   'display_labels': dict of label -> LaTeX string

LATTICE_DATA = {}

# -------------------------------------------------------
# CUB - Simple Cubic (Table 2)
# -------------------------------------------------------
LATTICE_DATA['CUB'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'M': [1/2, 1/2, 0],
        'R': [1/2, 1/2, 1/2],
        'X': [0, 1/2, 0],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'M'), ('M', 'Γ'), ('Γ', 'R'),
        ('R', 'X'), ('M', 'R'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'M': 'M', 'R': 'R', 'X': 'X',
    },
}

# -------------------------------------------------------
# CUB2 - Simple Cubic with doubled IBZ (23, m-3 point groups)
# For SG 195-206 (P-centered): include mirror-related X_A point
# -------------------------------------------------------
LATTICE_DATA['CUB2'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'X': [0, 1/2, 0],
        'XA': [1/2, 0, 0],
        'M': [1/2, 1/2, 0],
        'R': [1/2, 1/2, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'M'), ('M', 'XA'), ('XA', 'Γ'),
        ('Γ', 'R'), ('R', 'X'), ('R', 'M'), ('R', 'XA'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'X': 'X', 'XA': r'$X_A$', 'M': 'M', 'R': 'R',
    },
}

# -------------------------------------------------------
# FCC - Face-Centered Cubic (Table 3)
# -------------------------------------------------------
LATTICE_DATA['FCC'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'K': [3/8, 3/8, 3/4],
        'L': [1/2, 1/2, 1/2],
        'U': [5/8, 1/4, 5/8],
        'W': [1/2, 1/4, 3/4],
        'X': [1/2, 0, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'W'), ('W', 'K'), ('K', 'Γ'),
        ('Γ', 'L'), ('L', 'U'), ('U', 'W'), ('W', 'L'),
        ('L', 'K'), ('U', 'X'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'K': 'K', 'L': 'L',
        'U': 'U', 'W': 'W', 'X': 'X',
    },
}

# -------------------------------------------------------
# FCC2 - Face-Centered Cubic with doubled IBZ (23, m-3 point groups)
# For SG 195-206 (F-centered): include mirror-related X/U/W counterparts
# -------------------------------------------------------
LATTICE_DATA['FCC2'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'K': [3/8, 3/8, 3/4],
        'L': [1/2, 1/2, 1/2],
        'U': [5/8, 1/4, 5/8],
        'UA': [1/4, 5/8, 5/8],
        'W': [1/2, 1/4, 3/4],
        'WA': [1/4, 1/2, 3/4],
        'X': [1/2, 0, 1/2],
        'XA': [0, 1/2, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'W'), ('W', 'K'), ('K', 'Γ'),
        ('Γ', 'L'), ('L', 'U'), ('U', 'W'), ('W', 'L'),
        ('L', 'K'), ('U', 'X'),
        ('Γ', 'XA'), ('XA', 'WA'), ('WA', 'K'),
        ('L', 'UA'), ('UA', 'WA'), ('UA', 'XA'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'K': 'K', 'L': 'L',
        'U': 'U', 'UA': r'$U_A$',
        'W': 'W', 'WA': r'$W_A$',
        'X': 'X', 'XA': r'$X_A$',
    },
}

# -------------------------------------------------------
# BCC - Body-Centered Cubic (Table 4)
# -------------------------------------------------------
LATTICE_DATA['BCC'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'H': [1/2, -1/2, 1/2],
        'N': [0, 0, 1/2],
        'P': [1/4, 1/4, 1/4],
    },
    'kpath': [
        ('Γ', 'H'), ('H', 'N'), ('N', 'Γ'), ('Γ', 'P'),
        ('P', 'H'), ('P', 'N'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'H': 'H', 'N': 'N', 'P': 'P',
    },
}

# -------------------------------------------------------
# BCC2 - Body-Centered Cubic with doubled IBZ (23, m-3 point groups)
# For SG 195-206 (I-centered): include mirror-related H_A point
# -------------------------------------------------------
LATTICE_DATA['BCC2'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'H': [1/2, -1/2, 1/2],
        'HA': [-1/2, 1/2, 1/2],
        'N': [0, 0, 1/2],
        'P': [1/4, 1/4, 1/4],
    },
    'kpath': [
        ('Γ', 'H'), ('H', 'N'), ('N', 'Γ'), ('Γ', 'P'),
        ('P', 'H'), ('P', 'N'),
        ('Γ', 'HA'), ('HA', 'N'), ('P', 'HA'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'H': 'H', 'HA': r'$H_A$', 'N': 'N', 'P': 'P',
    },
}

# -------------------------------------------------------
# TET - Tetragonal (Table 5)
# -------------------------------------------------------
LATTICE_DATA['TET'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'A': [1/2, 1/2, 1/2],
        'M': [1/2, 1/2, 0],
        'R': [0, 1/2, 1/2],
        'X': [0, 1/2, 0],
        'Z': [0, 0, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'M'), ('M', 'Γ'), ('Γ', 'Z'),
        ('Z', 'R'), ('R', 'A'), ('A', 'Z'),
        ('X', 'R'), ('M', 'A'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'M': 'M',
        'R': 'R', 'X': 'X', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# TET2 - Tetragonal with doubled IBZ (4/m, -4, 4 point groups)
# For SG 75-88 (P-centered): the IBZ is a quadrant (double the TET octant)
# 8 vertices forming a rectangular prism
# -------------------------------------------------------
LATTICE_DATA['TET2'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'X': [0, 1/2, 0],
        'M': [1/2, 1/2, 0],
        'XA': [1/2, 0, 0],
        'Z': [0, 0, 1/2],
        'R': [0, 1/2, 1/2],
        'A': [1/2, 1/2, 1/2],
        'RA': [1/2, 0, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'M'), ('M', 'XA'), ('XA', 'Γ'),
        ('Γ', 'Z'),
        ('Z', 'R'), ('R', 'A'), ('A', 'RA'), ('RA', 'Z'),
        ('X', 'R'), ('M', 'A'), ('XA', 'RA'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'M': 'M',
        'R': 'R', 'X': 'X', 'Z': 'Z',
        'XA': r'$X_A$', 'RA': r'$R_A$',
    },
}

# -------------------------------------------------------
# BCT1 - Body-Centered Tetragonal 1, c < a (Table 6)
# -------------------------------------------------------
def _bct1_params(a, b=None, c=None, alpha=None):
    """BCT1: conventional a, c with c < a."""
    eta = (1 + c**2 / a**2) / 4
    return {'eta': eta}


def _bct1_kpoints(p):
    eta = p['eta']
    return {
        'Γ': [0, 0, 0],
        'M': [-1/2, 1/2, 1/2],
        'N': [0, 1/2, 0],
        'P': [1/4, 1/4, 1/4],
        'X': [0, 0, 1/2],
        'Z': [eta, eta, -eta],
        'Z1': [-eta, 1 - eta, eta],
    }

LATTICE_DATA['BCT1'] = {
    'params_func': _bct1_params,
    'kpoints_func': _bct1_kpoints,
    'kpath': [
        ('Γ', 'X'), ('X', 'M'), ('M', 'Γ'), ('Γ', 'Z'),
        ('Z', 'P'), ('P', 'N'), ('N', 'Z1'), ('Z1', 'M'),
        ('X', 'P'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'M': 'M', 'N': 'N', 'P': 'P',
        'X': 'X', 'Z': 'Z', 'Z1': r'$Z_1$',
    },
}

# -------------------------------------------------------
# BCT2 - Body-Centered Tetragonal 2, c > a (Table 7)
# -------------------------------------------------------
def _bct2_params(a, b=None, c=None, alpha=None):
    eta = (1 + a**2 / c**2) / 4
    zeta = a**2 / (2 * c**2)
    return {'eta': eta, 'zeta': zeta}


def _bct2_kpoints(p):
    eta, zeta = p['eta'], p['zeta']
    return {
        'Γ': [0, 0, 0],
        'N': [0, 1/2, 0],
        'P': [1/4, 1/4, 1/4],
        'Σ': [-eta, eta, eta],
        'Σ1': [eta, 1 - eta, -eta],
        'X': [0, 0, 1/2],
        'Y': [-zeta, zeta, 1/2],
        'Y1': [1/2, 1/2, -zeta],
        'Z': [1/2, 1/2, -1/2],
    }

LATTICE_DATA['BCT2'] = {
    'params_func': _bct2_params,
    'kpoints_func': _bct2_kpoints,
    'kpath': [
        ('Γ', 'X'), ('X', 'Y'), ('Y', 'Σ'), ('Σ', 'Γ'),
        ('Γ', 'Z'), ('Z', 'Σ1'), ('Σ1', 'N'), ('N', 'P'),
        ('P', 'Y1'), ('Y1', 'Z'), ('X', 'P'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'N': 'N', 'P': 'P',
        'Σ': r'$\Sigma$', 'Σ1': r'$\Sigma_1$',
        'X': 'X', 'Y': 'Y', 'Y1': r'$Y_1$', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORC - Simple Orthorhombic (Table 8)
# Convention: a < b < c
# -------------------------------------------------------
LATTICE_DATA['ORC'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'R': [1/2, 1/2, 1/2],
        'S': [1/2, 1/2, 0],
        'T': [0, 1/2, 1/2],
        'U': [1/2, 0, 1/2],
        'X': [1/2, 0, 0],
        'Y': [0, 1/2, 0],
        'Z': [0, 0, 1/2],
    },
    'kpath': [
        ('Γ', 'X'), ('X', 'S'), ('S', 'Y'), ('Y', 'Γ'),
        ('Γ', 'Z'), ('Z', 'U'), ('U', 'R'), ('R', 'T'),
        ('T', 'Z'), ('Y', 'T'), ('U', 'X'), ('S', 'R'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'R': 'R', 'S': 'S', 'T': 'T',
        'U': 'U', 'X': 'X', 'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORCF1 - Face-Centered Orthorhombic 1 (Table 9)
# Convention: a < b < c, 1/a^2 > 1/b^2 + 1/c^2
# -------------------------------------------------------
def _orcf1_params(a, b, c, alpha=None):
    zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4
    eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
    return {'zeta': zeta, 'eta': eta}


def _orcf1_kpoints(p):
    zeta, eta = p['zeta'], p['eta']
    return {
        'Γ': [0, 0, 0],
        'A': [1/2, 1/2 + zeta, zeta],
        'A1': [1/2, 1/2 - zeta, 1 - zeta],
        'L': [1/2, 1/2, 1/2],
        'T': [1, 1/2, 1/2],
        'X': [0, eta, eta],
        'X1': [1, 1 - eta, 1 - eta],
        'Y': [1/2, 0, 1/2],
        'Z': [1/2, 1/2, 0],
    }

LATTICE_DATA['ORCF1'] = {
    'params_func': _orcf1_params,
    'kpoints_func': _orcf1_kpoints,
    'kpath': [
        ('Γ', 'Y'), ('Y', 'T'), ('T', 'Z'), ('Z', 'Γ'),
        ('Γ', 'X'), ('X', 'A1'), ('A1', 'Y'),
        ('T', 'X1'), ('X', 'A'), ('A', 'Z'), ('L', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'A1': r'$A_1$',
        'L': 'L', 'T': 'T', 'X': 'X', 'X1': r'$X_1$',
        'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORCF2 - Face-Centered Orthorhombic 2 (Table 10)
# Convention: a < b < c, 1/a^2 < 1/b^2 + 1/c^2
# -------------------------------------------------------
def _orcf2_params(a, b, c, alpha=None):
    eta = (1 + a**2 / b**2 - a**2 / c**2) / 4
    delta = (1 + b**2 / a**2 - b**2 / c**2) / 4
    phi = (1 + c**2 / b**2 - c**2 / a**2) / 4
    return {'eta': eta, 'delta': delta, 'phi': phi}


def _orcf2_kpoints(p):
    eta, delta, phi = p['eta'], p['delta'], p['phi']
    return {
        'Γ': [0, 0, 0],
        'C': [1/2, 1/2 - eta, 1 - eta],
        'C1': [1/2, 1/2 + eta, eta],
        'D': [1/2 - delta, 1/2, 1 - delta],
        'D1': [1/2 + delta, 1/2, delta],
        'L': [1/2, 1/2, 1/2],
        'H': [1 - phi, 1/2 - phi, 1/2],
        'H1': [phi, 1/2 + phi, 1/2],
        'X': [0, 1/2, 1/2],
        'Y': [1/2, 0, 1/2],
        'Z': [1/2, 1/2, 0],
    }

LATTICE_DATA['ORCF2'] = {
    'params_func': _orcf2_params,
    'kpoints_func': _orcf2_kpoints,
    'kpath': [
        ('Γ', 'Y'), ('Y', 'C'), ('C', 'D'), ('D', 'X'),
        ('X', 'Γ'), ('Γ', 'Z'), ('Z', 'D1'), ('D1', 'H'),
        ('H', 'C'), ('C1', 'Z'), ('X', 'H1'), ('H', 'Y'),
        ('L', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'C': 'C', 'C1': r'$C_1$',
        'D': 'D', 'D1': r'$D_1$', 'L': 'L',
        'H': 'H', 'H1': r'$H_1$',
        'X': 'X', 'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORCF3 - Face-Centered Orthorhombic 3 (Table 9)
# Convention: a < b < c, 1/a^2 = 1/b^2 + 1/c^2
# Same k-points as ORCF1, different path
# -------------------------------------------------------
LATTICE_DATA['ORCF3'] = {
    'params_func': _orcf1_params,
    'kpoints_func': _orcf1_kpoints,
    'kpath': [
        ('Γ', 'Y'), ('Y', 'T'), ('T', 'Z'), ('Z', 'Γ'),
        ('Γ', 'X'), ('X', 'A1'), ('A1', 'Y'),
        ('X', 'A'), ('A', 'Z'), ('L', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'A1': r'$A_1$',
        'L': 'L', 'T': 'T', 'X': 'X', 'X1': r'$X_1$',
        'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORCI - Body-Centered Orthorhombic (Table 11)
# Convention: a < b < c
# -------------------------------------------------------
def _orci_params(a, b, c, alpha=None):
    zeta = (1 + a**2 / c**2) / 4
    eta = (1 + b**2 / c**2) / 4
    delta = (b**2 - a**2) / (4 * c**2)
    mu = (a**2 + b**2) / (4 * c**2)
    return {'zeta': zeta, 'eta': eta, 'delta': delta, 'mu': mu}


def _orci_kpoints(p):
    zeta, eta, delta, mu = p['zeta'], p['eta'], p['delta'], p['mu']
    return {
        'Γ': [0, 0, 0],
        'L': [-mu, mu, 1/2 - delta],
        'L1': [mu, -mu, 1/2 + delta],
        'L2': [1/2 - delta, 1/2 + delta, -mu],
        'R': [0, 1/2, 0],
        'S': [1/2, 0, 0],
        'T': [0, 0, 1/2],
        'W': [1/4, 1/4, 1/4],
        'X': [-zeta, zeta, zeta],
        'X1': [zeta, 1 - zeta, -zeta],
        'Y': [eta, -eta, eta],
        'Y1': [1 - eta, eta, -eta],
        'Z': [1/2, 1/2, -1/2],
    }

LATTICE_DATA['ORCI'] = {
    'params_func': _orci_params,
    'kpoints_func': _orci_kpoints,
    'kpath': [
        ('Γ', 'X'), ('X', 'L'), ('L', 'T'), ('T', 'W'),
        ('W', 'R'), ('R', 'X1'), ('X1', 'Z'), ('Z', 'Γ'),
        ('Γ', 'Y'), ('Y', 'S'), ('S', 'W'),
        ('L1', 'Y'), ('Y1', 'Z'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'L': 'L', 'L1': r'$L_1$', 'L2': r'$L_2$',
        'R': 'R', 'S': 'S', 'T': 'T', 'W': 'W',
        'X': 'X', 'X1': r'$X_1$',
        'Y': 'Y', 'Y1': r'$Y_1$', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# ORCC - C-Centered Orthorhombic (Table 12)
# Convention: a < b
# -------------------------------------------------------
def _orcc_params(a, b, c=None, alpha=None):
    zeta = (1 + a**2 / b**2) / 4
    return {'zeta': zeta}


def _orcc_kpoints(p):
    zeta = p['zeta']
    return {
        'Γ': [0, 0, 0],
        'A': [zeta, zeta, 1/2],
        'A1': [-zeta, 1 - zeta, 1/2],
        'R': [0, 1/2, 1/2],
        'S': [0, 1/2, 0],
        'T': [-1/2, 1/2, 1/2],
        'X': [zeta, zeta, 0],
        'X1': [-zeta, 1 - zeta, 0],
        'Y': [-1/2, 1/2, 0],
        'Z': [0, 0, 1/2],
    }

LATTICE_DATA['ORCC'] = {
    'params_func': _orcc_params,
    'kpoints_func': _orcc_kpoints,
    'kpath': [
        ('Γ', 'X'), ('X', 'S'), ('S', 'R'), ('R', 'A'),
        ('A', 'Z'), ('Z', 'Γ'), ('Γ', 'Y'), ('Y', 'X1'),
        ('X1', 'A1'), ('A1', 'T'), ('T', 'Y'), ('Z', 'T'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'A1': r'$A_1$',
        'R': 'R', 'S': 'S', 'T': 'T',
        'X': 'X', 'X1': r'$X_1$', 'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# HEX - Hexagonal (Table 13)
# -------------------------------------------------------
LATTICE_DATA['HEX'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'A': [0, 0, 1/2],
        'H': [1/3, 1/3, 1/2],
        'K': [1/3, 1/3, 0],
        'L': [1/2, 0, 1/2],
        'M': [1/2, 0, 0],
    },
    'kpath': [
        ('Γ', 'M'), ('M', 'K'), ('K', 'Γ'), ('Γ', 'A'),
        ('A', 'L'), ('L', 'H'), ('H', 'A'),
        ('L', 'M'), ('K', 'H'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'H': 'H',
        'K': 'K', 'L': 'L', 'M': 'M',
    },
}

# -------------------------------------------------------
# HEX2 - Hexagonal with doubled IBZ (6/m, -6, 6 point groups)
# For SG 168-176: the IBZ is a 60° sector (double the 30° HEX wedge)
# 8 vertices forming a quadrilateral prism
# -------------------------------------------------------
LATTICE_DATA['HEX2'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'A': [0, 0, 1/2],
        'K': [1/3, 1/3, 0],
        'H': [1/3, 1/3, 1/2],
        'M': [1/2, 0, 0],
        'L': [1/2, 0, 1/2],
        'MA': [0, 1/2, 0],
        'LA': [0, 1/2, 1/2],
    },
    'kpath': [
        ('Γ', 'M'), ('M', 'K'), ('K', 'MA'), ('MA', 'Γ'),
        ('Γ', 'A'),
        ('A', 'L'), ('L', 'H'), ('H', 'LA'), ('LA', 'A'),
        ('L', 'M'), ('K', 'H'), ('MA', 'LA'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'A': 'A', 'H': 'H',
        'K': 'K', 'L': 'L', 'M': 'M',
        'MA': r'$M_A$', 'LA': r'$L_A$',
    },
}

# -------------------------------------------------------
# RHL1 - Rhombohedral 1, alpha < 90° (Table 14)
# -------------------------------------------------------
def _rhl1_params(a, b=None, c=None, alpha=None):
    """alpha is the rhombohedral angle in degrees."""
    cos_a = np.cos(np.radians(alpha))
    eta = (1 + 4 * cos_a) / (2 + 4 * cos_a)
    nu = 3/4 - eta / 2
    return {'eta': eta, 'nu': nu}


def _rhl1_kpoints(p):
    eta, nu = p['eta'], p['nu']
    return {
        'Γ': [0, 0, 0],
        'B': [eta, 1/2, 1 - eta],
        'B1': [1/2, 1 - eta, eta - 1],
        'F': [1/2, 1/2, 0],
        'L': [1/2, 0, 0],
        'L1': [0, 0, -1/2],
        'P': [eta, nu, nu],
        'P1': [1 - nu, 1 - nu, 1 - eta],
        'P2': [nu, nu, eta - 1],
        'Q': [1 - nu, nu, 0],
        'X': [nu, 0, -nu],
        'Z': [1/2, 1/2, 1/2],
    }

LATTICE_DATA['RHL1'] = {
    'params_func': _rhl1_params,
    'kpoints_func': _rhl1_kpoints,
    'kpath': [
        ('Γ', 'L'), ('L', 'B1'),
        ('B', 'Z'), ('Z', 'Γ'), ('Γ', 'X'),
        ('Q', 'F'), ('F', 'P1'), ('P1', 'Z'),
        ('L', 'P'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'B': 'B', 'B1': r'$B_1$',
        'F': 'F', 'L': 'L', 'L1': r'$L_1$',
        'P': 'P', 'P1': r'$P_1$', 'P2': r'$P_2$',
        'Q': 'Q', 'X': 'X', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# RHL2 - Rhombohedral 2, alpha > 90° (Table 15)
# -------------------------------------------------------
def _rhl2_params(a, b=None, c=None, alpha=None):
    """alpha is the rhombohedral angle in degrees."""
    alpha_rad = np.radians(alpha)
    eta = 1 / (2 * np.tan(alpha_rad / 2)**2)
    nu = 3/4 - eta / 2
    return {'eta': eta, 'nu': nu}


def _rhl2_kpoints(p):
    eta, nu = p['eta'], p['nu']
    return {
        'Γ': [0, 0, 0],
        'F': [1/2, -1/2, 0],
        'L': [1/2, 0, 0],
        'P': [1 - nu, -nu, 1 - nu],
        'P1': [nu, nu - 1, nu - 1],
        'Q': [eta, eta, eta],
        'Q1': [1 - eta, -eta, -eta],
        'Z': [1/2, -1/2, 1/2],
    }

LATTICE_DATA['RHL2'] = {
    'params_func': _rhl2_params,
    'kpoints_func': _rhl2_kpoints,
    'kpath': [
        ('Γ', 'P'), ('P', 'Z'), ('Z', 'Q'), ('Q', 'Γ'),
        ('Γ', 'F'), ('F', 'P1'), ('P1', 'Q1'), ('Q1', 'L'),
        ('L', 'Z'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'F': 'F', 'L': 'L',
        'P': 'P', 'P1': r'$P_1$',
        'Q': 'Q', 'Q1': r'$Q_1$', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# MCL - Simple Monoclinic (HPKOT/seekpath mP1, Table 87)
# Crystallographic convention (unique axis b): beta is non-right angle
# -------------------------------------------------------
def _mcl_params(a, b, c, alpha):
    """
    Keep the argument name 'alpha' for API compatibility.
    Here it is interpreted as monoclinic beta (angle between a and c).
    """
    beta_rad = np.radians(alpha)
    eta = (1 + (a / c) * np.cos(beta_rad)) / (2 * np.sin(beta_rad)**2)
    nu = 1/2 + eta * c * np.cos(beta_rad) / a
    return {'eta': eta, 'nu': nu}


def _mcl_kpoints(p):
    eta, nu = p['eta'], p['nu']
    return {
        'Γ': [0, 0, 0],
        'Z': [0, 1/2, 0],
        'B': [0, 0, 1/2],
        'B2': [0, 0, -1/2],
#        'Y': [1/2, 0, 0],
        'Y2': [-1/2, 0, 0],
#        'C': [1/2, 1/2, 0],
        'C2': [-1/2, 1/2, 0],
        'D': [0, 1/2, 1/2],
        'D2': [0, 1/2, -1/2],
        'A': [-1/2, 0, 1/2],
        'E': [-1/2, 1/2, 1/2],
        'H': [-eta, 0, 1 - nu],
        'H2': [-1 + eta, 0, nu],
        'H4': [-eta, 0, -nu],
        'M': [-eta, 1/2, 1 - nu],
        'M2': [-1 + eta, 1/2, nu],
        'M4': [-eta, 1/2, -nu],
    }

LATTICE_DATA['MCL'] = {
    'params_func': _mcl_params,
    'kpoints_func': _mcl_kpoints,
    'kpath': [
        ('Γ', 'Z'), ('Z', 'D'), ('D', 'B'), ('B', 'Γ'),
        ('Γ', 'A'), ('A', 'E'), ('E', 'Z'),
        ('Z', 'C2'), ('C2', 'Y2'), ('Y2', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$',
        'A': 'A', 'B': 'B', 'B2': r'$B_2$',
        'C': 'C', 'C2': r'$C_2$',
        'D': 'D', 'D2': r'$D_2$',
        'E': 'E', 'H': 'H', 'H2': r'$H_2$', 'H4': r'$H_4$',
        'M': 'M', 'M2': r'$M_2$', 'M4': r'$M_4$',
        'Y': 'Y', 'Y2': r'$Y_2$', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# MCLC1 - C-Centered Monoclinic 1 (HPKOT Table 88)
# b < a*sin(beta)
# -------------------------------------------------------
def _mclc1_params(a, b, c, alpha):
    alpha_rad = np.radians(alpha)
    zeta = (2 + (a / c) * np.cos(alpha_rad)) / (4 * np.sin(alpha_rad)**2)
    eta = 1/2 - 2 * zeta * c * np.cos(alpha_rad) / a
    psi = 3/4 - b**2 / (4 * a**2 * np.sin(alpha_rad)**2)
    phi = psi - (3/4 - psi) * a * np.cos(alpha_rad) / c
    tau = (psi - phi) / (3/4 - psi)
    return {'zeta': zeta, 'eta': eta, 'psi': psi, 'phi': phi, 'tau': tau}


def _mclc1_kpoints(p):
    zeta, eta, psi, phi, tau = (
        p['zeta'], p['eta'], p['psi'], p['phi'], p['tau'])
    return {
        'Γ': [0, 0, 0],
        'Y2': [-1/2, 1/2, 0],
        'Y4': [1/2, -1/2, 0],
        'A': [0, 0, 1/2],
        'M2': [-1/2, 1/2, 1/2],
        'V': [1/2, 0, 0],
        'V2': [0, 1/2, 0],
        'L2': [0, 1/2, 1/2],
        'C': [1 - psi, 1 - psi, 0],
        'C2': [-1 + psi, psi, 0],
        'C4': [psi, -1 + psi, 0],
        'D': [-1 + phi, phi, 1/2],
        'D2': [1 - phi, 1 - phi, 1/2],
        'E': [-1 + zeta, 1 - zeta, 1 - eta],
        'E2': [-zeta, zeta, eta],
        'E4': [zeta, -zeta, 1 - eta],
        # Extra vertices needed to close the IRBZ polyhedron (hidden from plot).
        '_F4': [-1/2 + zeta + psi, -1/2 - zeta + psi, 1 - eta],
        '_Q1': [-1/2 - zeta + psi, -1/2 + zeta + psi, eta],                 # C2(F4) + G
        '_Q2': [1/2 + zeta - psi, 3/2 - zeta - psi, 1 - eta],               # E + C - Y2
        '_P1': [-3/2 + zeta + psi, 1/2 - zeta + psi, 1 - eta],              # C2 - Y2 + E
    }


LATTICE_DATA['MCLC1'] = {
    'params_func': _mclc1_params,
    'kpoints_func': _mclc1_kpoints,
    'kpath': [
        ('Γ', 'C'),
        ('C2', 'Y2'), ('Y2', 'Γ'), ('Γ', 'M2'), ('M2', 'D'),
        ('D2', 'A'), ('A', 'Γ'),
        ('L2', 'Γ'), ('Γ', 'V2'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$',
        'A': 'A',
        'C': 'C', 'C2': r'$C_2$', 'C4': r'$C_4$',
        'D': 'D', 'D2': r'$D_2$',
        'E': 'E', 'E2': r'$E_2$', 'E4': r'$E_4$',
        'L2': r'$L_2$',
        'M2': r'$M_2$',
        'V': 'V', 'V2': r'$V_2$',
        'Y2': r'$Y_2$', 'Y4': r'$Y_4$',
    },
}

# -------------------------------------------------------
# MCLC2 - C-Centered Monoclinic 2 (HPKOT Table 89)
# b > a*sin(beta), -a*cos(beta)/c + a^2*sin(beta)^2/b^2 < 1
# -------------------------------------------------------
def _mclc2_params(a, b, c, alpha):
    alpha_rad = np.radians(alpha)
    mu = (1 + a**2 / b**2) / 4
    delta = -a * c * np.cos(alpha_rad) / (2 * b**2)
    zeta = (a**2 / b**2 + (1 + (a / c) * np.cos(alpha_rad)) / np.sin(alpha_rad)**2) / 4
    eta = 1/2 - 2 * zeta * c * np.cos(alpha_rad) / a
    phi = 1 + zeta - 2 * mu
    psi = eta - 2 * delta
    return {'mu': mu, 'delta': delta, 'zeta': zeta, 'eta': eta, 'phi': phi, 'psi': psi}


def _mclc2_kpoints(p):
    mu, delta, zeta, eta, phi, psi = (
        p['mu'], p['delta'], p['zeta'], p['eta'], p['phi'], p['psi'])
    return {
        'Γ': [0, 0, 0],
        'Y': [1/2, 1/2, 0],
        'A': [0, 0, 1/2],
        'M': [1/2, 1/2, 1/2],
        'V2': [0, 1/2, 0],
        'L2': [0, 1/2, 1/2],
        'F': [-1 + phi, 1 - phi, 1 - psi],
        'F2': [1 - phi, phi, psi],
        'F4': [phi, 1 - phi, 1 - psi],
        'H': [-zeta, zeta, eta],
        'H2': [zeta, 1 - zeta, 1 - eta],
        'H4': [zeta, -zeta, 1 - eta],
        'G': [-mu, mu, delta],
        'G2': [mu, 1 - mu, -delta],
        'G4': [mu, -mu, -delta],
        'G6': [1 - mu, mu, delta],
    }


LATTICE_DATA['MCLC2'] = {
    'params_func': _mclc2_params,
    'kpoints_func': _mclc2_kpoints,
    'kpath': [
        ('Γ', 'Y'), ('Y', 'M'), ('M', 'A'), ('A', 'Γ'),
        ('L2', 'Γ'), ('Γ', 'V2'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$',
        'A': 'A',
        'F': 'F', 'F2': r'$F_2$', 'F4': r'$F_4$',
        'G': 'G', 'G2': r'$G_2$', 'G4': r'$G_4$', 'G6': r'$G_6$',
        'H': 'H', 'H2': r'$H_2$', 'H4': r'$H_4$',
        'L2': r'$L_2$',
        'M': 'M',
        'V2': r'$V_2$',
        'Y': 'Y',
    },
}

# -------------------------------------------------------
# MCLC3 - C-Centered Monoclinic 3 (HPKOT Table 90)
# b > a*sin(beta), -a*cos(beta)/c + a^2*sin(beta)^2/b^2 >= 1
# -------------------------------------------------------
def _mclc3_params(a, b, c, alpha):
    alpha_rad = np.radians(alpha)
    zeta = (a**2 / b**2 + (1 + (a / c) * np.cos(alpha_rad)) / np.sin(alpha_rad)**2) / 4
    rho = 1 - zeta * b**2 / a**2
    eta = 1/2 - 2 * zeta * c * np.cos(alpha_rad) / a
    mu = eta / 2 + a**2 / (4 * b**2) + a * c * np.cos(alpha_rad) / (2 * b**2)
    nu = 2 * mu - zeta
    omega = c * (1 - 4 * nu + a**2 * np.sin(alpha_rad)**2 / b**2) / (2 * a * np.cos(alpha_rad))
    delta = -1/4 + omega / 2 - zeta * c * np.cos(alpha_rad) / a
    return {'zeta': zeta, 'rho': rho, 'eta': eta, 'mu': mu, 'nu': nu,
            'omega': omega, 'delta': delta}


def _mclc3_kpoints(p):
    zeta, rho, eta = p['zeta'], p['rho'], p['eta']
    mu, nu, omega, delta = p['mu'], p['nu'], p['omega'], p['delta']
    return {
        'Γ': [0, 0, 0],
        'Y': [1/2, 1/2, 0],
        'A': [0, 0, 1/2],
        'M2': [-1/2, 1/2, 1/2],
        'V': [1/2, 0, 0],
        'V2': [0, 1/2, 0],
        'L2': [0, 1/2, 1/2],
        'I': [-1 + rho, rho, 1/2],
        'I2': [1 - rho, 1 - rho, 1/2],
        'K': [-nu, nu, omega],
        'K2': [-1 + nu, 1 - nu, 1 - omega],
        'K4': [1 - nu, nu, omega],
        'H': [-zeta, zeta, eta],
        'H2': [zeta, 1 - zeta, 1 - eta],
        'H4': [zeta, -zeta, 1 - eta],
        'N': [-mu, mu, delta],
        'N2': [mu, 1 - mu, -delta],
        'N4': [mu, -mu, -delta],
        'N6': [1 - mu, mu, delta],
    }


LATTICE_DATA['MCLC3'] = {
    'params_func': _mclc3_params,
    'kpoints_func': _mclc3_kpoints,
    'kpath': [
        ('Γ', 'A'), ('A', 'I2'),
        ('I', 'M2'), ('M2', 'Γ'), ('Γ', 'Y'),
        ('L2', 'Γ'), ('Γ', 'V2'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$',
        'A': 'A',
        'H': 'H', 'H2': r'$H_2$', 'H4': r'$H_4$',
        'I': 'I', 'I2': r'$I_2$',
        'K': 'K', 'K2': r'$K_2$', 'K4': r'$K_4$',
        'L2': r'$L_2$',
        'M2': r'$M_2$',
        'N': 'N', 'N2': r'$N_2$', 'N4': r'$N_4$', 'N6': r'$N_6$',
        'V': 'V', 'V2': r'$V_2$',
        'Y': 'Y',
    },
}

# Aliases: SC boundary cases use the same k-point tables as their neighbors
LATTICE_DATA['MCLC2_SC'] = LATTICE_DATA['MCLC1']  # SC MCLC2 (k_gamma=90) -> mC1 table
LATTICE_DATA['MCLC4_SC'] = LATTICE_DATA['MCLC2']  # SC MCLC4 (cond=1) -> mC2 table
LATTICE_DATA['MCLC4'] = LATTICE_DATA['MCLC2']
LATTICE_DATA['MCLC5'] = LATTICE_DATA['MCLC3']


# -------------------------------------------------------
# TRI1a - Triclinic 1a (Table 20)
# All reciprocal angles < 90° k_alpha* < 90, k_beta* < 90, k_gamma* < 90
# -------------------------------------------------------
LATTICE_DATA['TRI1a'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'L': [1/2, 1/2, 0],
        'M': [0, 1/2, 1/2],
        'N': [1/2, 0, 1/2],
        'R': [1/2, 1/2, 1/2],
        'X': [1/2, 0, 0],
        'Y': [0, 1/2, 0],
        'Z': [0, 0, 1/2],
    },
    'kpath': [
        ('X', 'Γ'), ('Γ', 'Y'), ('L', 'Γ'), ('Γ', 'Z'),
        ('N', 'Γ'), ('Γ', 'M'), ('R', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'L': 'L', 'M': 'M', 'N': 'N',
        'R': 'R', 'X': 'X', 'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# TRI1b - Triclinic 1b (Table 20)
# Some reciprocal angles = 90°
# Same k-points as TRI1a, different path
# -------------------------------------------------------
LATTICE_DATA['TRI1b'] = {
    'kpoints': LATTICE_DATA['TRI1a']['kpoints'],
    'kpath': [
        ('X', 'Γ'), ('Γ', 'Y'), ('L', 'Γ'), ('Γ', 'Z'),
        ('N', 'Γ'), ('Γ', 'M'), ('R', 'Γ'),
    ],
    'display_labels': LATTICE_DATA['TRI1a']['display_labels'],
}

# -------------------------------------------------------
# TRI2a - Triclinic 2a (Table 21)
# All reciprocal angles > 90° k_alpha* > 90, k_beta* > 90, k_gamma* > 90
# -------------------------------------------------------
LATTICE_DATA['TRI2a'] = {
    'kpoints': {
        'Γ': [0, 0, 0],
        'L': [1/2, -1/2, 0],
        'M': [0, 0, 1/2],
        'N': [-1/2, -1/2, 1/2],
        'R': [0, -1/2, 1/2],
        'X': [0, -1/2, 0],
        'Y': [1/2, 0, 0],
        'Z': [-1/2, 0, 1/2],
    },
    'kpath': [
        ('X', 'Γ'), ('Γ', 'Y'), ('L', 'Γ'), ('Γ', 'Z'),
        ('N', 'Γ'), ('Γ', 'M'), ('R', 'Γ'),
    ],
    'display_labels': {
        'Γ': r'$\Gamma$', 'L': 'L', 'M': 'M', 'N': 'N',
        'R': 'R', 'X': 'X', 'Y': 'Y', 'Z': 'Z',
    },
}

# -------------------------------------------------------
# TRI2b - Triclinic 2b (Table 21)
# Some reciprocal angles = 90°
# Same k-points as TRI2a, different path
# -------------------------------------------------------
LATTICE_DATA['TRI2b'] = {
    'kpoints': LATTICE_DATA['TRI2a']['kpoints'],
    'kpath': [
        ('X', 'Γ'), ('Γ', 'Y'), ('L', 'Γ'), ('Γ', 'Z'),
        ('N', 'Γ'), ('Γ', 'M'), ('R', 'Γ'),
    ],
    'display_labels': LATTICE_DATA['TRI2a']['display_labels'],
}


# ============================================================================
# Convenience function: get k-points for a given lattice type
# ============================================================================
def get_kpoints(lattice_type, a=None, b=None, c=None, alpha=None):
    """
    Return the k-points dict for the given lattice type.

    Parameters
    ----------
    lattice_type : str
        E.g. 'CUB', 'FCC', 'BCT1', 'MCLC5', ...
    a, b, c : float
        Conventional cell axial lengths (needed for parametric types).
    alpha : float
        Relevant angle in degrees (rhombohedral angle for RHL,
        conventional alpha for MCL/MCLC).

    Returns
    -------
    dict : {label: [b1, b2, b3], ...}
    """
    data = LATTICE_DATA[lattice_type]

    if 'kpoints' in data:
        return data['kpoints']

    params = data['params_func'](a, b, c, alpha)
    return data['kpoints_func'](params)


def get_kpath(lattice_type):
    """Return the k-path for the given lattice type."""
    return LATTICE_DATA[lattice_type]['kpath']


def get_display_labels(lattice_type):
    """Return the display label dict for the given lattice type."""
    return LATTICE_DATA[lattice_type]['display_labels']


def get_params(lattice_type, a=None, b=None, c=None, alpha=None):
    """Return the parameter dict (eta, nu, etc.) if parametric, else empty."""
    data = LATTICE_DATA[lattice_type]
    if 'params_func' in data:
        return data['params_func'](a, b, c, alpha)
    return {}


# ============================================================================
# List all supported lattice types
# ============================================================================
ALL_LATTICE_TYPES = list(LATTICE_DATA.keys())

FIXED_LATTICES = [k for k in ALL_LATTICE_TYPES if 'kpoints' in LATTICE_DATA[k]]
PARAMETRIC_LATTICES = [k for k in ALL_LATTICE_TYPES if 'params_func' in LATTICE_DATA[k]]

if __name__ == '__main__':
    print("=" * 60)
    print("Setyawan-Curtarolo k-point database")
    print("=" * 60)
    print(f"Total lattice types: {len(ALL_LATTICE_TYPES)}")
    print(f"  Fixed:      {FIXED_LATTICES}")
    print(f"  Parametric: {PARAMETRIC_LATTICES}")
    print()
    for lt in ALL_LATTICE_TYPES:
        kp = get_kpoints(lt, a=5.0, b=6.0, c=7.0, alpha=60.0)
        print(f"{lt:8s}: {len(kp)} k-points, {len(get_kpath(lt))} path segments")



