import numpy as np
import spglib
import spinspg
from ase.io import read
import sys
import os
import re
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
        f.write(f"SPIN SYMMETRY LOG\n")
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
def _operation_class_indices(spin_rotations, flip):
    indices = []
    for i, s_rot in enumerate(spin_rotations):
        is_flip = np.isclose(np.linalg.det(s_rot), -1)
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


def operation_count_summary(rotations, spin_rotations):
    flip_indices = _operation_class_indices(spin_rotations, flip=True)
    preserve_indices = _operation_class_indices(spin_rotations, flip=False)
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


def write_flip_ops_to_file(filename, rotations, spin_rotations, verbose=True):
    """
    Filters operations where Spin Rotation is a flip (det approx -1).
    For each spin-flip spatial operation R, also include the inversion-extended
    partner -R, then deduplicate. Translations do not affect reciprocal-space
    k-point mapping.
    """
    flip_indices = _operation_class_indices(spin_rotations, flip=True)
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


def write_preserve_ops_to_file(filename, rotations, spin_rotations, verbose=True):
    """Write inversion-extended spin-preserving point operations."""
    preserve_indices = _operation_class_indices(spin_rotations, flip=False)
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


def has_spin_flip_inversion(rotations, spin_rotations):
    """Return True when inversion is an actual spin-flip operation."""
    inversion = -np.eye(3, dtype=int)
    for rot, s_rot in zip(rotations, spin_rotations):
        if (
            np.array_equal(np.array(rot, dtype=int), inversion)
            and np.isclose(np.linalg.det(s_rot), -1)
        ):
            return True
    return False


def has_spin_flip_translation(rotations, translations, spin_rotations):
    """Return True when a pure nonzero translation flips spin."""
    identity = np.eye(3, dtype=int)
    for rot, trans, s_rot in zip(rotations, translations, spin_rotations):
        if not np.array_equal(np.array(rot, dtype=int), identity):
            continue
        if not np.isclose(np.linalg.det(s_rot), -1):
            continue
        trans_mod = np.mod(np.array(trans, dtype=float), 1.0)
        trans_mod[np.isclose(trans_mod, 1.0, atol=1e-8)] = 0.0
        if not np.allclose(trans_mod, 0.0, atol=1e-8):
            return True
    return False


def altermagnetic_diagnostic(rotations, translations, spin_rotations):
    """Summarize whether spin-flip operations indicate AM splitting or PT."""
    flip_indices = _operation_class_indices(spin_rotations, flip=True)
    flip_point_ops, _ = _collect_point_ops(rotations, flip_indices)
    inversion = -np.eye(3, dtype=int)
    pt = any(np.array_equal(op, inversion) for op in flip_point_ops)
    if pt:
        return "PT symmetry detected, not altermagnet."
    if has_spin_flip_translation(rotations, translations, spin_rotations):
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

# ==========================================
# MAIN FUNCTION
# ==========================================
def run(structure_file, moments_str, verbose=True):
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
            num_atoms = len(pmg_struct)
            structure = None   # not used for mcif path
            if verbose:
                print(f"Successfully loaded '{structure_file}' containing {num_atoms} atoms.")
        else:
            structure = read(structure_file)
            lattice   = structure.get_cell()
            positions = structure.get_scaled_positions()
            numbers   = structure.get_atomic_numbers()
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
            if not moments_str:
                user_mags = []
            else:
                user_mags = parse_magmoms(moments_str)
        except ValueError:
            print("Error: Invalid input. Please enter numbers or VASP-style repeats like '4*0 2*1'.")
            return False
        if len(user_mags) < num_atoms:
            user_mags.extend([0.0] * (num_atoms - len(user_mags)))
        elif len(user_mags) > num_atoms:
            user_mags = user_mags[:num_atoms]
        magmoms = np.zeros((num_atoms, 3))
        for i, m in enumerate(user_mags):
            magmoms[i] = [0, 0, m]

    if verbose:
        print(f"Using magnetic moments:\n{magmoms}")

    # Run spglib for Magnetic Space Group Label (if available).
    # spglib's magnetic API expects magnetic moments as the 4th item of the
    # cell tuple: (lattice, positions, numbers, magmoms).
    msg_label = "Not found"
    try:
        mag_cell = (lattice, positions, numbers, magmoms)
        mag_dataset = spglib.get_magnetic_symmetry_dataset(
            mag_cell, symprec=1e-5
        )
        if mag_dataset:
            uni_number = getattr(mag_dataset, 'uni_number', None)
            msg_type = getattr(mag_dataset, 'msg_type', None)
            if uni_number is None and hasattr(mag_dataset, 'get'):
                uni_number = mag_dataset.get('uni_number')
                msg_type = mag_dataset.get('msg_type')

            if uni_number is not None:
                msg_type_info = spglib.get_magnetic_spacegroup_type(int(uni_number))
                if msg_type_info:
                    bns_number = getattr(msg_type_info, 'bns_number', None)
                    og_number = getattr(msg_type_info, 'og_number', None)
                    litvin_number = getattr(msg_type_info, 'litvin_number', None)
                    msg_number = getattr(msg_type_info, 'number', None)
                    if bns_number is None and hasattr(msg_type_info, 'get'):
                        bns_number = msg_type_info.get('bns_number')
                        og_number = msg_type_info.get('og_number')
                        litvin_number = msg_type_info.get('litvin_number')
                        msg_number = msg_type_info.get('number')

                    msg_label = (
                        f"BNS {bns_number}, OG {og_number} "
                        f"(UNI {uni_number}, Litvin {litvin_number}, "
                        f"parent SG {msg_number}, type {msg_type})"
                    )
                else:
                    msg_label = f"UNI {uni_number} (type {msg_type})"
    except Exception as e:
        msg_label = f"Not found ({e})"

    if verbose:
        print(f"Magnetic Space Group: {msg_label}")

    # --- PART 4: Spin Space Group (SpinSPG) ---
    if verbose:
        print("\n" + "="*40)
        print("4. Spin Space Group Analysis")
        print("="*40)

    # Run spinspg
    sog, rotations, translations, spin_rotations = spinspg.get_spin_symmetry(
        lattice, positions, numbers, magmoms, symprec=1e-5
    )

    # Print info
    counts = operation_count_summary(rotations, spin_rotations)
    unique_point_operations = counts['unique_point_operations']
    spin_split_diagnostic = altermagnetic_diagnostic(
        rotations, translations, spin_rotations
    )

    if verbose:
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
Spin-Only Group Type: {sog}
Magnetic Space Group Label: {msg_label}"""

    # 1. Write the full readable log with LABELS
    write_operations_to_file("spin_operations.txt", rotations, translations, spin_rotations, label_info_str, verbose=verbose)

    # 2. Write the automation file
    flip_filename = "spin_flip_operations.txt"
    preserve_filename = "spin_preserve_operations.txt"
    flip_count = write_flip_ops_to_file(flip_filename, rotations, spin_rotations, verbose=verbose)
    preserve_count = write_preserve_ops_to_file(preserve_filename, rotations, spin_rotations, verbose=verbose)
    return {
        'structure_file': structure_file,
        'num_atoms': num_atoms,
        'moments': magmoms,
        'space_group': non_mag_label,
        'point_group': point_group,
        'laue_group': laue_group,
        'magnetic_space_group': msg_label,
        'spin_group': str(sog),
        'total_operations': len(rotations),
        'unique_point_operations': unique_point_operations,
        'actual_spin_flip_point_operations': counts['actual_spin_flip_point_operations'],
        'actual_spin_preserve_point_operations': counts['actual_spin_preserve_point_operations'],
        'extended_spin_flip_point_operations': counts['extended_spin_flip_point_operations'],
        'extended_spin_preserve_point_operations': counts['extended_spin_preserve_point_operations'],
        'extended_spin_flip_operations': counts['extended_spin_flip_operations'],
        'extended_spin_preserve_operations': counts['extended_spin_preserve_operations'],
        'pt_spin_flip': has_spin_flip_inversion(rotations, spin_rotations),
        'ut_spin_flip': has_spin_flip_translation(rotations, translations, spin_rotations),
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

    print("Enter magnetic moments (space-separated, e.g., '1 -1'):")
    moments_input = input("Moments: ").strip()

    success = run(filename, moments_input)
    if not success:
        sys.exit(1)
