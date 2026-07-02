import numpy as np
import spglib
from findspingroup import (
    find_spin_group_basic,
    find_spin_group_basic_from_data,
    find_spin_group_input_ssg,
)
from findspingroup.find_spin_group import (
    Tolerances,
    _find_spin_group_input_ssg_from_parsed,
)
from findspingroup.data.MSGMPG_DB import MSG_INT_TO_BNS
from ase.io import read
import sys
import os
import re
import sympy as sp
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module=r"pymatgen\.io\.cif")


def _laue_group_from_point_group(point_group):
    pg = str(point_group).strip().replace(' ', '')
    pg = pg.replace('−', '-').replace('bar', '-')
    mapping = {
        '1': '-1', '-1': '-1',
        '2': '2/m', 'm': '2/m', '2/m': '2/m',
        '222': 'mmm', 'mm2': 'mmm', 'mmm': 'mmm',
        '4': '4/m', '-4': '4/m', '4/m': '4/m',
        '422': '4/mmm', '4mm': '4/mmm', '-42m': '4/mmm',
        '-4m2': '4/mmm', '4/mmm': '4/mmm',
        '3': '-3', '-3': '-3',
        '32': '-3m', '3m': '-3m', '-3m': '-3m',
        '6': '6/m', '-6': '6/m', '6/m': '6/m',
        '622': '6/mmm', '6mm': '6/mmm', '-6m2': '6/mmm',
        '-62m': '6/mmm', '6/mmm': '6/mmm',
        '23': 'm-3', 'm-3': 'm-3',
        '432': 'm-3m', '-43m': 'm-3m', 'm-3m': 'm-3m',
    }
    return mapping.get(pg)


# --- HELPER 1: Write FULL details for human reading ---
def write_operations_to_file(filename, rotations, translations, spin_rotations, label_info, verbose=True):
    """Writes all spin symmetry operations to a text file."""
    with open(filename, 'w') as f:
        f.write("="*40 + "\n")
        f.write("SPIN SYMMETRY LOG\n")
        f.write("="*40 + "\n\n")
        f.write(f"{label_info}\n\n")
        f.write(f"Full space-group operations: {len(rotations)}\n")
        f.write(f"Unique point operations: {count_unique_point_operations(rotations)}\n")
        f.write("-" * 40 + "\n")
        for i in range(len(rotations)):
            f.write(f"Operation {i+1}:\n")
            f.write(f"  Rotation:\n{rotations[i]}\n")
            f.write(f"  Translation:\n{translations[i]}\n")
            f.write(f"  Spin Rotation:\n{spin_rotations[i]}\n")
            f.write("-" * 20 + "\n")
    if verbose:
        print(f"[INFO] All operations written to '{filename}'")

# --- HELPER 2: Write ONLY Flip Operations for automation ---
def _spin_axis_from_moments(magmoms):
    """Return a normalized nonzero moment direction for a collinear structure."""
    for moment in np.asarray(magmoms, dtype=float):
        norm = np.linalg.norm(moment)
        if norm > 1e-10:
            return moment / norm
    raise ValueError("No nonzero magnetic moment found.")


def _parse_spin_axis(direction):
    """Parse FindSpinGroup's Cartesian collinear-axis summary when available."""
    if not direction:
        return None
    try:
        entries = [
            float(sp.N(sp.sympify(value.strip())))
            for value in str(direction).split(",")
        ]
        axis = np.asarray(entries, dtype=float)
        norm = np.linalg.norm(axis)
        return axis / norm if norm > 1e-10 else None
    except (TypeError, ValueError, sp.SympifyError):
        return None


def _operation_class_indices(spin_rotations, spin_axis, flip):
    indices = []
    for i, s_rot in enumerate(spin_rotations):
        mapped_axis = np.asarray(s_rot, dtype=float) @ spin_axis
        is_flip = np.allclose(mapped_axis, -spin_axis, atol=1e-7)
        is_preserve = np.allclose(mapped_axis, spin_axis, atol=1e-7)
        if not is_flip and not is_preserve:
            continue
        if is_flip == flip:
            indices.append(i)
    return indices


def _collect_point_ops(rotations, indices, include_inversion=False):
    ops = []
    source_indices = []
    for i in indices:
        variants = [np.array(rotations[i], dtype=int)]
        if include_inversion:
            variants.append(-np.array(rotations[i], dtype=int))
        for rot in variants:
            if not any(np.array_equal(rot, ex) for ex in ops):
                ops.append(rot)
                source_indices.append(i + 1)
    return ops, source_indices


def operation_count_summary(rotations, spin_rotations, spin_axis):
    flip_indices = _operation_class_indices(spin_rotations, spin_axis, flip=True)
    preserve_indices = _operation_class_indices(spin_rotations, spin_axis, flip=False)
    flip_point, _ = _collect_point_ops(rotations, flip_indices)
    preserve_point, _ = _collect_point_ops(rotations, preserve_indices)
    flip_ext_point, _ = _collect_point_ops(rotations, flip_indices, include_inversion=True)
    preserve_ext_point, _ = _collect_point_ops(rotations, preserve_indices, include_inversion=True)
    return {
        'full_space_group_operations': len(rotations),
        'actual_spin_flip_operations': len(flip_indices),
        'actual_spin_preserve_operations': len(preserve_indices),
        'unique_point_operations': count_unique_point_operations(rotations),
        'actual_spin_flip_point_operations': len(flip_point),
        'actual_spin_preserve_point_operations': len(preserve_point),
        'extended_spin_flip_point_operations': len(flip_ext_point),
        'extended_spin_preserve_point_operations': len(preserve_ext_point),
        'extended_spin_flip_operations': 2 * len(flip_indices),
        'extended_spin_preserve_operations': 2 * len(preserve_indices),
    }


def write_flip_ops_to_file(filename, rotations, spin_rotations, spin_axis, verbose=True):
    """
    Filters operations whose spin rotation reverses the collinear spin axis.
    For each spin-flip spatial operation R, also include the inversion-extended
    partner -R, then deduplicate. Translations do not affect reciprocal-space
    k-point mapping.
    """
    flip_indices = _operation_class_indices(spin_rotations, spin_axis, flip=True)
    flip_ops, source_indices = _collect_point_ops(
        rotations, flip_indices, include_inversion=True
    )

    if not flip_ops:
        if verbose:
            print("\n[WARNING] No spin-flipping operations found! File not created.")
        return 0

    with open(filename, 'w') as f:
        f.write(f"# Found {len(flip_ops)} inversion-extended spin-flipping point operations\n")
        f.write(f"# Original Indices: {source_indices}\n")
        for i, rot in enumerate(flip_ops):
            f.write(f"Operation_{i+1}\n")
            # Write matrix row by row
            for row in rot:
                f.write(f"{row[0]} {row[1]} {row[2]}\n")
            f.write("\n")

    if verbose:
        print(f"[INFO] {len(flip_ops)} spin-flipping matrices written to '{filename}'")
    return len(flip_ops)


def write_preserve_ops_to_file(filename, rotations, spin_rotations, spin_axis, verbose=True):
    """Write inversion-extended spin-preserving point operations."""
    preserve_indices = _operation_class_indices(spin_rotations, spin_axis, flip=False)
    preserve_ops, source_indices = _collect_point_ops(
        rotations, preserve_indices, include_inversion=True
    )

    if not preserve_ops:
        return 0

    with open(filename, 'w') as f:
        f.write(f"# Found {len(preserve_ops)} inversion-extended spin-preserving point operations\n")
        f.write(f"# Original Indices: {source_indices}\n")
        for i, rot in enumerate(preserve_ops):
            f.write(f"Operation_{i+1}\n")
            for row in rot:
                f.write(f"{row[0]} {row[1]} {row[2]}\n")
            f.write("\n")

    if verbose:
        print(f"[INFO] {len(preserve_ops)} spin-preserving matrices written to '{filename}'")
    return len(preserve_ops)


def count_unique_point_operations(rotations):
    """Count unique spatial point operations, ignoring translations."""
    unique_rots = []
    for rot in rotations:
        rot_arr = np.array(rot, dtype=int)
        if not any(np.array_equal(rot_arr, existing) for existing in unique_rots):
            unique_rots.append(rot_arr)
    return len(unique_rots)


def has_spin_flip_inversion(rotations, spin_rotations, spin_axis):
    """Return True when inversion is an actual spin-flip operation."""
    inversion = -np.eye(3, dtype=int)
    for i in _operation_class_indices(spin_rotations, spin_axis, flip=True):
        rot = rotations[i]
        if (
            np.array_equal(np.array(rot, dtype=int), inversion)
        ):
            return True
    return False


def has_spin_flip_translation(rotations, translations, spin_rotations, spin_axis):
    """Return True when a pure nonzero translation flips spin."""
    identity = np.eye(3, dtype=int)
    for i in _operation_class_indices(spin_rotations, spin_axis, flip=True):
        rot, trans = rotations[i], translations[i]
        if not np.array_equal(np.array(rot, dtype=int), identity):
            continue
        trans_mod = np.mod(np.array(trans, dtype=float), 1.0)
        trans_mod[np.isclose(trans_mod, 1.0, atol=1e-8)] = 0.0
        if not np.allclose(trans_mod, 0.0, atol=1e-8):
            return True
    return False


def altermagnetic_diagnostic(rotations, translations, spin_rotations, spin_axis, magnetic_phase):
    """Summarize whether spin-flip operations indicate AM splitting or PT."""
    if "Altermagnet" in str(magnetic_phase):
        return ""
    flip_indices = _operation_class_indices(spin_rotations, spin_axis, flip=True)
    flip_point_ops, _ = _collect_point_ops(rotations, flip_indices)
    inversion = -np.eye(3, dtype=int)
    pt = any(np.array_equal(op, inversion) for op in flip_point_ops)
    if pt:
        return "PT symmetry detected, not altermagnet."
    if has_spin_flip_translation(rotations, translations, spin_rotations, spin_axis):
        return "Ut symmetry detected, not altermagnet."
    if flip_point_ops:
        return ""
    return "No spin-flip point operation detected."


def parse_magmoms(moments_str):
    """
    Parse scalar collinear moments from manual input.

    Supports both expanded input, e.g. ``0 0 0 0 1 1``, and VASP MAGMOM-style
    shorthand, e.g. ``4*0 2*1``.  Whitespace around ``*`` is accepted.
    """
    if not moments_str:
        return []

    float_pat = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    normalized = re.sub(
        rf"(\d+)\s*\*\s*({float_pat})",
        r"\1*\2",
        moments_str.strip(),
    )

    values = []
    for token in normalized.split():
        if "*" in token:
            parts = token.split("*")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(f"invalid MAGMOM token '{token}'")
            count = int(parts[0])
            if count < 0:
                raise ValueError(f"negative repeat count in '{token}'")
            values.extend([float(parts[1])] * count)
        else:
            values.append(float(token))
    return values


def parse_cartesian_spin_axis(axis_str):
    """Parse and normalize a Cartesian collinear spin axis."""
    if axis_str is None or not str(axis_str).strip():
        return np.array([0.0, 0.0, 1.0])
    parts = str(axis_str).split()
    if len(parts) != 3:
        raise ValueError("spin axis must contain three Cartesian components")
    axis = np.asarray([float(value) for value in parts], dtype=float)
    norm = np.linalg.norm(axis)
    if norm <= 1e-10:
        raise ValueError("spin axis cannot be the zero vector")
    return axis / norm


def _display_ssg_symbol(symbol):
    """Make FindSpinGroup's database symbol readable in terminal output."""
    if not symbol:
        return "Unknown"
    return (
        str(symbol)
        .replace("\u221e", "infinity")
        .replace("\u0401\u043e", "infinity")
    )


def _fsg_operations_from_payload(payload):
    """Convert FindSpinGroup input-setting SSG operations to NumPy arrays."""
    operations = payload["ssg"]["ops"]
    rotations = np.asarray([op["real_rotation"] for op in operations], dtype=float)
    if np.allclose(rotations, np.rint(rotations), atol=1e-7):
        rotations = np.rint(rotations).astype(int)
    translations = np.asarray([op["translation"] for op in operations], dtype=float)
    spin_rotations = np.asarray([op["spin_rotation"] for op in operations], dtype=float)
    return rotations, translations, spin_rotations


def _deduplicate_collinear_operations(rotations, translations, spin_rotations, spin_axis):
    """Drop spin-only duplicates while keeping the input-setting spatial operations."""
    compact = []
    for rot, trans, spin_rot in zip(rotations, translations, spin_rotations):
        mapped_axis = np.asarray(spin_rot, dtype=float) @ spin_axis
        if np.allclose(mapped_axis, spin_axis, atol=1e-7):
            spin_class = 1
        elif np.allclose(mapped_axis, -spin_axis, atol=1e-7):
            spin_class = -1
        else:
            continue
        if any(
            old_class == spin_class
            and np.allclose(old_rot, rot, atol=1e-7)
            and np.allclose(np.mod(old_trans - trans, 1.0), 0.0, atol=1e-7)
            for old_rot, old_trans, _, old_class in compact
        ):
            continue
        compact.append((rot, trans, spin_rot, spin_class))
    return (
        np.asarray([item[0] for item in compact]),
        np.asarray([item[1] for item in compact]),
        np.asarray([item[2] for item in compact]),
    )


def _magnetic_type_label(msg_type):
    labels = {
        1: "I",
        2: "II",
        3: "III",
        4: "IV",
    }
    return labels.get(msg_type, str(msg_type))


def compute_msg_without_soc(rotations, translations, spin_rotations, spin_axis):
    """
    Compute the FindMagSym-style MSG without SOC from spin-space operations.

    Only the original FindSpinGroup operations are used here. The
    inversion-extended point operations written for k mapping are intentionally
    excluded because they are not physical space-group operations.
    """
    msg_rotations = []
    msg_translations = []
    msg_time_reversals = []

    for rot, trans, spin_rot in zip(rotations, translations, spin_rotations):
        mapped_axis = np.asarray(spin_rot, dtype=float) @ spin_axis
        if np.allclose(mapped_axis, spin_axis, atol=1e-7):
            time_reversal = False
        elif np.allclose(mapped_axis, -spin_axis, atol=1e-7):
            time_reversal = True
        else:
            continue
        rot_arr = np.asarray(rot, dtype=float)
        rounded_rot = np.rint(rot_arr)
        if not np.allclose(rot_arr, rounded_rot, atol=1e-7):
            raise ValueError("non-integer spatial rotation found in MSG without SOC operations")
        msg_rotations.append(rounded_rot.astype(int))
        msg_translations.append(np.asarray(trans, dtype=float))
        msg_time_reversals.append(time_reversal)

    if not msg_rotations:
        return None, 0, 0

    msg_type = spglib.get_magnetic_spacegroup_type_from_symmetry(
        np.asarray(msg_rotations, dtype=int),
        np.asarray(msg_translations, dtype=float),
        np.asarray(msg_time_reversals, dtype=bool),
    )
    return msg_type, len(msg_rotations), sum(bool(value) for value in msg_time_reversals)


def format_msg_without_soc(msg_type):
    if msg_type is None:
        return "Unknown"
    bns_number, bns_symbol = MSG_INT_TO_BNS.get(
        msg_type.uni_number,
        (msg_type.bns_number, None),
    )
    if bns_symbol:
        return f"{bns_symbol} (BNS {bns_number}), Type {_magnetic_type_label(msg_type.type)}"
    return f"BNS {bns_number}, Type {_magnetic_type_label(msg_type.type)}"

# ==========================================
# MAIN FUNCTION
# ==========================================
def run(structure_file, moments_str, verbose=True, spin_axis_cart=None):
    """
    Run spin-flip operations analysis.
    Called by auto-generate-general-kpath.py or used standalone.
    Returns True on success, False on failure.
    """
    # 1. Structure Loading
    if verbose:
        print("="*40)
        print("1. Structure Loading")
        print("="*40)

    try:
        if not os.path.exists(structure_file):
            raise FileNotFoundError
        is_mcif = structure_file.lower().endswith('.mcif')
        if is_mcif:
            from pymatgen.io.cif import CifParser
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                parser = CifParser(structure_file)
                pmg_struct = parser.parse_structures(primitive=False)[0]
            lattice   = np.array(pmg_struct.lattice.matrix)
            positions = np.array([site.frac_coords for site in pmg_struct])
            numbers   = np.array([site.specie.Z for site in pmg_struct])
            elements  = [str(site.specie) for site in pmg_struct]
            num_atoms = len(pmg_struct)
            structure = None   # not used for mcif path
            if verbose:
                print(f"Successfully loaded '{structure_file}' containing {num_atoms} atoms.")
        else:
            structure = read(structure_file)
            lattice   = structure.get_cell()
            positions = structure.get_scaled_positions()
            numbers   = structure.get_atomic_numbers()
            elements  = structure.get_chemical_symbols()
            num_atoms = len(structure)
            pmg_struct = None
            if verbose:
                print(f"Successfully loaded '{structure_file}' containing {num_atoms} atoms.")
    except FileNotFoundError:
        print(f"Error: File '{structure_file}' not found.")
        return False
    except Exception as e:
        print(f"Error reading file: {e}")
        return False

    # --- PART 2: Non-Magnetic Space Group (SPG) ---
    if verbose:
        print("\n" + "="*40)
        print("2. Non-Magnetic Space Group Analysis")
        print("="*40)

    cell = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell)

    non_mag_label = "Unknown"
    point_group = "Unknown"
    laue_group = "Unknown"
    if dataset:
        non_mag_label = f"{dataset.international} ({dataset.number})"
        point_group = dataset.pointgroup
        laue_group = _laue_group_from_point_group(point_group) or "Unknown"
        if verbose:
            print(f"Space Group: {non_mag_label}")
    else:
        if verbose:
            print("Non-magnetic symmetry detection failed.")

    # --- PART 3: Magnetic Configuration ---
    if verbose:
        print("\n" + "="*40)
        print("3. Magnetic Configuration")
        print("="*40)

    is_mcif = structure_file.lower().endswith('.mcif')
    magmoms = None

    if is_mcif and pmg_struct is not None:
        try:
            magmoms = np.array([
                np.array(site.properties['magmom'].moment)
                if 'magmom' in site.properties else np.zeros(3)
                for site in pmg_struct
            ])
            if verbose:
                print(f"Read moments from mcif:\n{magmoms}")
        except Exception as e:
            if verbose:
                print(f"[Warning] Could not read moments from mcif: {e}. Falling back to manual input.")

    if magmoms is None:
        if verbose:
            print(f"Moments: {moments_str}")
        try:
            manual_axis = parse_cartesian_spin_axis(spin_axis_cart)
            if not moments_str:
                user_mags = []
            else:
                user_mags = parse_magmoms(moments_str)
        except ValueError as exc:
            print(
                "Error: Invalid manual magnetic input. Enter a nonzero Cartesian "
                "axis such as '0 0 1' and moments such as '4*0 2*1'. "
                f"({exc})"
            )
            return False
        if len(user_mags) < num_atoms:
            user_mags.extend([0.0] * (num_atoms - len(user_mags)))
        elif len(user_mags) > num_atoms:
            user_mags = user_mags[:num_atoms]
        magmoms = np.asarray(user_mags, dtype=float)[:, None] * manual_axis[None, :]
        if verbose:
            print(f"Cartesian spin axis: {manual_axis}")

    if verbose:
        print(f"Using magnetic moments:\n{magmoms}")

    # --- PART 4: Spin Space Group (FindSpinGroup) ---
    if verbose:
        print("\n" + "="*40)
        print("4. Spin Space Group Analysis")
        print("="*40)

    try:
        spin_axis = _spin_axis_from_moments(magmoms)
        if is_mcif:
            fsg_basic = find_spin_group_basic(structure_file)
            fsg_input = find_spin_group_input_ssg(structure_file)
            reported_axis = _parse_spin_axis(
                fsg_input["summary"].get("input_spin_only_direction")
            )
            if reported_axis is not None:
                spin_axis = reported_axis
        else:
            occupancies = [1.0] * num_atoms
            fsg_basic = find_spin_group_basic_from_data(
                structure_file,
                lattice,
                positions,
                elements,
                occupancies,
                magmoms,
                input_spin_setting="cartesian",
            )
            # FindSpinGroup currently exposes this input-setting operation route
            # for parsed data internally; keep its use isolated here until a
            # public from-data equivalent is available.
            fsg_input = _find_spin_group_input_ssg_from_parsed(
                structure_file,
                lattice,
                positions,
                elements,
                occupancies,
                magmoms,
                Tolerances(0.02, 0.02, 0.00002, m_matrix_tol=0.01),
                input_spin_setting="cartesian",
                source_format="poscar",
            )
        rotations, translations, spin_rotations = _fsg_operations_from_payload(fsg_input)
        rotations, translations, spin_rotations = _deduplicate_collinear_operations(
            rotations, translations, spin_rotations, spin_axis
        )
    except Exception as e:
        print(f"Error running FindSpinGroup: {e}")
        return False

    sog = f"{fsg_basic.get('conf', 'Unknown')}(axis={np.array2string(spin_axis, precision=6)})"
    magnetic_phase = fsg_basic.get("magnetic_phase", "Unknown")
    msg_label = (
        f"{fsg_basic.get('msg_symbol', 'Unknown')} "
        f"(BNS {fsg_basic.get('msg_bns_number', 'Unknown')}, "
        f"OG {fsg_basic.get('msg_og_number', 'Unknown')})"
    )
    msg_without_soc, msg_without_soc_ops, msg_without_soc_tr_ops = compute_msg_without_soc(
        rotations, translations, spin_rotations, spin_axis
    )
    msg_without_soc_label = format_msg_without_soc(msg_without_soc)
    ssg_label = fsg_basic.get("index", "Unknown")
    ssg_symbol = _display_ssg_symbol(
        fsg_input["summary"].get("input_ssg_database_symbol")
    )
    g0_label = f"{fsg_basic.get('g0_symbol', 'Unknown')} ({fsg_basic.get('g0_number', 'Unknown')})"
    l0_label = f"{fsg_basic.get('l0_symbol', 'Unknown')} ({fsg_basic.get('l0_number', 'Unknown')})"
    empg = fsg_basic.get("empg", "Unknown")

    # Print info
    counts = operation_count_summary(rotations, spin_rotations, spin_axis)
    unique_point_operations = counts['unique_point_operations']
    spin_split_diagnostic = altermagnetic_diagnostic(
        rotations, translations, spin_rotations, spin_axis, magnetic_phase
    )

    if verbose:
        print(f"Magnetic Phase: {magnetic_phase}")
        print(f"Oriented SSG: {ssg_label}")
        print(f"SSG Symbol (Chen-Liu): {ssg_symbol}")
        print(f"G0: {g0_label}; L0: {l0_label}; EMPG: {empg}")
        print(f"MSG with SOC: {msg_label}")
        print(f"MSG without SOC: {msg_without_soc_label}")
        print(f"Spin-Only Group Type: {sog}")
        print(f"Full space-group operations: {counts['full_space_group_operations']}")
        print(f"Unique point operations: {unique_point_operations}")
        print(
            "Actual point operations: "
            f"{counts['actual_spin_flip_point_operations']} spin-flip, "
            f"{counts['actual_spin_preserve_point_operations']} spin-preserving"
        )
        print(
            "Inversion-extended point operations for k mapping: "
            f"{counts['extended_spin_flip_point_operations']} spin-flip, "
            f"{counts['extended_spin_preserve_point_operations']} spin-preserving"
        )
        if spin_split_diagnostic:
            print(f"\033[1mWarning! {spin_split_diagnostic}\033[0m")

    # --- PART 5: Output Files ---
    if verbose:
        print("\n" + "="*40)
        print("5. Saving Results")
        print("="*40)

    # Prepare label info for text file
    label_info_str = f"""Non-Magnetic Label: {non_mag_label}
Magnetic Phase: {magnetic_phase}
Oriented SSG: {ssg_label}
SSG Symbol (Chen-Liu): {ssg_symbol}
G0: {g0_label}
L0: {l0_label}
Effective MPG: {empg}
Spin-Only Group Type: {sog}
MSG with SOC: {msg_label}
MSG without SOC: {msg_without_soc_label}"""

    # 1. Write the full readable log with LABELS
    write_operations_to_file("spin_operations.txt", rotations, translations, spin_rotations, label_info_str, verbose=verbose)

    # 2. Write the automation file
    flip_filename = "spin_flip_operations.txt"
    preserve_filename = "spin_preserve_operations.txt"
    flip_count = write_flip_ops_to_file(flip_filename, rotations, spin_rotations, spin_axis, verbose=verbose)
    preserve_count = write_preserve_ops_to_file(preserve_filename, rotations, spin_rotations, spin_axis, verbose=verbose)
    return {
        'structure_file': structure_file,
        'num_atoms': num_atoms,
        'moments': magmoms,
        'space_group': non_mag_label,
        'point_group': point_group,
        'laue_group': laue_group,
        'magnetic_space_group': msg_label,
        'magnetic_space_group_without_soc': msg_without_soc_label,
        'msg_without_soc_bns_number': getattr(msg_without_soc, 'bns_number', None),
        'msg_without_soc_type': getattr(msg_without_soc, 'type', None),
        'msg_without_soc_operation_count': msg_without_soc_ops,
        'msg_without_soc_time_reversal_count': msg_without_soc_tr_ops,
        'spin_group': str(sog),
        'magnetic_phase': magnetic_phase,
        'ssg_index': ssg_label,
        'ssg_symbol': ssg_symbol,
        'g0_symbol': fsg_basic.get('g0_symbol'),
        'g0_number': fsg_basic.get('g0_number'),
        'l0_symbol': fsg_basic.get('l0_symbol'),
        'l0_number': fsg_basic.get('l0_number'),
        'empg': empg,
        'findspingroup_warning': fsg_input['summary'].get('warning'),
        'total_operations': len(rotations),
        'unique_point_operations': unique_point_operations,
        'actual_spin_flip_point_operations': counts['actual_spin_flip_point_operations'],
        'actual_spin_preserve_point_operations': counts['actual_spin_preserve_point_operations'],
        'extended_spin_flip_point_operations': counts['extended_spin_flip_point_operations'],
        'extended_spin_preserve_point_operations': counts['extended_spin_preserve_point_operations'],
        'extended_spin_flip_operations': counts['extended_spin_flip_operations'],
        'extended_spin_preserve_operations': counts['extended_spin_preserve_operations'],
        'pt_spin_flip': has_spin_flip_inversion(rotations, spin_rotations, spin_axis),
        'ut_spin_flip': has_spin_flip_translation(rotations, translations, spin_rotations, spin_axis),
        'spin_split_diagnostic': spin_split_diagnostic,
        'spin_flip_operations': flip_count,
        'spin_preserve_operations': preserve_count,
        'saved_files': [
            'spin_operations.txt',
            flip_filename,
            preserve_filename,
        ],
    }


# ==========================================
# STANDALONE SCRIPT (python find_sf_operations.py)
# ==========================================
if __name__ == "__main__":
    filename = input("Enter structure file name (default: POSCAR): ").strip()
    if not filename:
        filename = "POSCAR"

    spin_axis_input = None
    if filename.lower().endswith(".mcif"):
        print("Detected .mcif file -- magnetic moments will be read from file.")
        moments_input = ""
    else:
        print("Spin axis in Cartesian coordinates (default: 0 0 1):")
        spin_axis_input = input("Axis: ").strip()
        print("Enter magnetic moments along this axis (space-separated, e.g., '1 -1'):")
        moments_input = input("Moments: ").strip()

    success = run(filename, moments_input, spin_axis_cart=spin_axis_input)
    if not success:
        sys.exit(1)
