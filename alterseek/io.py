"""Structure and file I/O: reading input structures (POSCAR/.vasp, .cif,
.mcif, and any other ASE-readable format), writing helper POSCAR/mcif files,
lattice geometry utilities, atomic text writes, and band-plot config writing.

Extracted from alterseek_path.py (restructuring phase 4, renamed from
io_vasp.py). pymatgen/ase imports are function-local (kept lazy).
"""
import itertools
import os
import uuid
import numpy as np

from .find_sf_operations import (
    fit_magmoms_to_structure,
    parse_cartesian_spin_axis,
    parse_magmoms,
)


def _atomic_write_text(path, text):
    """Write UTF-8 text beside *path*, then atomically replace the target."""
    target = os.path.abspath(os.fspath(path))
    parent = os.path.dirname(target)
    # O_CREAT|O_EXCL with mode 0666 (instead of tempfile.mkstemp) lets the
    # kernel apply the process umask itself, so a brand-new target gets the
    # same mode a plain open() would create — without mkstemp's 0600 or a
    # process-global os.umask() round-trip that could race other threads.
    temporary = os.path.join(
        parent, f".{os.path.basename(target)}.{uuid.uuid4().hex}.tmp"
    )
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    fd = os.open(temporary, flags, 0o666)
    try:
        handle = os.fdopen(fd, "w", encoding="utf-8", newline="\n")
    except Exception:
        # os.fdopen never took ownership, so the raw descriptor is still ours.
        os.close(fd)
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise
    try:
        with handle:
            handle.write(text)
        # Keep an existing target's mode (e.g. group-readable outputs in
        # shared cluster directories); a new target keeps the umask default.
        try:
            os.chmod(temporary, os.stat(target).st_mode & 0o7777)
        except OSError:
            pass
        os.replace(temporary, target)
    except Exception:
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise


def _group_poscar_sites(elements, positions, moment_keys=None):
    groups = []
    site_to_group = {}
    for idx, element in enumerate(elements):
        key = (element, moment_keys[idx]) if moment_keys is not None else (element,)
        if key not in site_to_group:
            site_to_group[key] = len(groups)
            groups.append((element, []))
        groups[site_to_group[key]][1].append(idx)
    ordered_indices = [idx for _, indices in groups for idx in indices]
    symbols = [symbol for symbol, _ in groups]
    counts = [len(indices) for _, indices in groups]
    ordered_positions = [positions[idx] for idx in ordered_indices]
    return symbols, counts, ordered_positions, ordered_indices


def _write_poscar(path, title, lattice, symbols, counts, positions):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"{title}\n")
        f.write("1.0\n")
        for row in lattice:
            f.write(f"  {row[0]: .10f}  {row[1]: .10f}  {row[2]: .10f}\n")
        f.write(" ".join(symbols) + "\n")
        f.write(" ".join(str(count) for count in counts) + "\n")
        f.write("Direct\n")
        for pos in positions:
            f.write(f"  {pos[0]: .10f}  {pos[1]: .10f}  {pos[2]: .10f}\n")


def _read_grouped_poscar(path):
    # pymatgen handles scale factors, Cartesian/Direct mode, and Selective
    # dynamics; import locally to keep module import light.
    from pymatgen.io.vasp.inputs import Poscar

    structure = Poscar.from_file(path, check_for_potcar=False).structure
    lattice = np.array(structure.lattice.matrix, dtype=float)
    expanded_symbols = [site.specie.symbol for site in structure]
    positions = [np.array(site.frac_coords, dtype=float) for site in structure]
    return lattice, expanded_symbols, positions


def _write_poscar_from_sites(path, title, lattice, elements, positions):
    symbols, counts, grouped_positions, _ = _group_poscar_sites(elements, positions)
    _write_poscar(path, title, lattice, symbols, counts, grouped_positions)


def _write_without_species(source_path, target_path, species_to_remove, title):
    lattice, elements, positions = _read_grouped_poscar(source_path)
    kept = [
        (element, position)
        for element, position in zip(elements, positions)
        if element not in set(species_to_remove)
    ]
    if not kept:
        raise RuntimeError(f"No atoms left after removing {species_to_remove} from {source_path}.")
    kept_elements, kept_positions = zip(*kept)
    _write_poscar_from_sites(target_path, title, lattice, list(kept_elements), list(kept_positions))


def _reciprocal_from_poscar(path):
    lattice, _, _ = _read_grouped_poscar(path)
    return 2 * np.pi * np.linalg.inv(np.array(lattice, dtype=float)).T


def _dedupe_frac_positions(positions, tol=1e-7):
    unique = []
    for pos in positions:
        wrapped = np.mod(np.array(pos, dtype=float), 1.0)
        duplicate = False
        for existing in unique:
            delta = wrapped - existing
            delta -= np.rint(delta)
            if np.linalg.norm(delta) < tol:
                duplicate = True
                break
        if not duplicate:
            unique.append(wrapped)
    return unique


def _min_periodic_cart_distance(frac_positions, lattice):
    if len(frac_positions) < 2:
        return float("inf")
    min_dist = float("inf")
    for i, pos_i in enumerate(frac_positions):
        for pos_j in frac_positions[i + 1:]:
            delta = np.array(pos_i) - np.array(pos_j)
            delta -= np.rint(delta)
            dist = np.linalg.norm(delta @ lattice)
            min_dist = min(min_dist, dist)
    return min_dist


def _lattice_lengths_angles(lattice):
    a_vec, b_vec, c_vec = [np.array(row, dtype=float) for row in lattice]
    lengths = [np.linalg.norm(vec) for vec in (a_vec, b_vec, c_vec)]

    def angle(u, v):
        denom = np.linalg.norm(u) * np.linalg.norm(v)
        cosang = np.dot(u, v) / denom if denom > 0 else 1.0
        return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))

    alpha = angle(b_vec, c_vec)
    beta = angle(a_vec, c_vec)
    gamma = angle(a_vec, b_vec)
    return lengths, (alpha, beta, gamma)


def _write_magnetic_mcif(path, title, lattice, elements, positions, moments_cart):
    lengths, angles = _lattice_lengths_angles(lattice)
    lattice_axes = np.array(lattice, dtype=float)
    unit_axes = lattice_axes / np.linalg.norm(lattice_axes, axis=1)[:, None]
    moments_crystal = np.array(moments_cart, dtype=float) @ np.linalg.inv(unit_axes)
    label_counts = {}

    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"data_{title.replace(' ', '_')}\n")
        f.write("_symmetry_space_group_name_H-M 'P 1'\n")
        f.write("_space_group_IT_number 1\n")
        f.write(f"_cell_length_a    {lengths[0]:.10f}\n")
        f.write(f"_cell_length_b    {lengths[1]:.10f}\n")
        f.write(f"_cell_length_c    {lengths[2]:.10f}\n")
        f.write(f"_cell_angle_alpha {angles[0]:.10f}\n")
        f.write(f"_cell_angle_beta  {angles[1]:.10f}\n")
        f.write(f"_cell_angle_gamma {angles[2]:.10f}\n")
        f.write("loop_\n")
        f.write("_space_group_symop_operation_xyz\n")
        f.write("'x,y,z'\n")
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        labels = []
        for element, pos in zip(elements, positions):
            label_counts[element] = label_counts.get(element, 0) + 1
            label = f"{element}{label_counts[element]}"
            labels.append(label)
            wrapped = np.mod(np.array(pos, dtype=float), 1.0)
            f.write(
                f"{label} {element} "
                f"{wrapped[0]:.10f} {wrapped[1]:.10f} {wrapped[2]:.10f}\n"
            )
        f.write("loop_\n")
        f.write("_atom_site_moment.label\n")
        f.write("_atom_site_moment.crystalaxis_x\n")
        f.write("_atom_site_moment.crystalaxis_y\n")
        f.write("_atom_site_moment.crystalaxis_z\n")
        for label, moment in zip(labels, moments_crystal):
            f.write(f"{label} {moment[0]:.10f} {moment[1]:.10f} {moment[2]:.10f}\n")


def _match_periodic_standard_sites(
    target_positions,
    target_types,
    candidate_positions,
    candidate_types,
    target_lattice,
):
    """Match two standardized structures up to one origin shift.

    The structures are already in the same lattice basis at this point, but
    spglib's structural and magnetic standardizations may choose different
    origins and site orders. Return candidate indices in target-site order and
    the largest Cartesian mismatch.
    """
    target_positions = np.asarray(target_positions, dtype=float)
    target_types = np.asarray(target_types, dtype=int)
    candidate_positions = np.asarray(candidate_positions, dtype=float)
    candidate_types = np.asarray(candidate_types, dtype=int)
    target_lattice = np.asarray(target_lattice, dtype=float)
    if len(target_positions) != len(candidate_positions):
        return None
    if sorted(target_types.tolist()) != sorted(candidate_types.tolist()):
        return None

    best = None
    first_type = target_types[0]
    for first_candidate in np.flatnonzero(candidate_types == first_type):
        shift = target_positions[0] - candidate_positions[first_candidate]
        shifted = np.mod(candidate_positions + shift, 1.0)
        used = set()
        mapping = []
        max_distance = 0.0
        for target_pos, atomic_number in zip(target_positions, target_types):
            choices = []
            for candidate in np.flatnonzero(candidate_types == atomic_number):
                candidate = int(candidate)
                if candidate in used:
                    continue
                delta = shifted[candidate] - target_pos
                delta -= np.rint(delta)
                distance = float(np.linalg.norm(delta @ target_lattice))
                choices.append((distance, candidate))
            if not choices:
                mapping = []
                break
            distance, candidate = min(choices)
            used.add(candidate)
            mapping.append(candidate)
            max_distance = max(max_distance, distance)
        if mapping and (best is None or max_distance < best[1]):
            best = (mapping, max_distance)
    return best


def _write_seekpath_standard_mcif(
    target_structure_path,
    output_path,
    title,
    source_lattice,
    source_positions,
    source_elements,
    source_moments_cart,
    symprec=1e-3,
):
    """Write spins in the exact standardized conventional cell SeeK-path uses.

    Structural and magnetic standardization can return the same conventional
    cell with different axis order, Cartesian orientation, origin, and site
    order. Align the magnetic standard cell to the already-written SeeK-path
    VASP before writing the MCIF; never copy moment rows by index.
    """
    import spglib
    from pymatgen.core import Element, Structure

    target = Structure.from_file(target_structure_path)
    target_lattice = np.asarray(target.lattice.matrix, dtype=float)
    target_positions = np.asarray(target.frac_coords, dtype=float)
    target_types = np.asarray([site.specie.Z for site in target], dtype=int)

    source_lattice = np.asarray(source_lattice, dtype=float)
    source_positions = np.asarray(source_positions, dtype=float)
    source_elements = [str(element) for element in source_elements]
    source_types = np.asarray(
        [Element(element).Z for element in source_elements], dtype=int)
    source_moments_cart = np.asarray(source_moments_cart, dtype=float)
    if source_moments_cart.shape != (len(source_positions), 3):
        raise ValueError("magnetic moments must contain one Cartesian vector per site")

    tolerance = 1e-3 if symprec is None else float(symprec)
    magnetic = spglib.get_magnetic_symmetry_dataset(
        (source_lattice, source_positions, source_types, source_moments_cart),
        is_axial=True,
        symprec=tolerance,
        mag_symprec=1e-5,
    )
    if magnetic is None:
        raise RuntimeError("spglib could not standardize the magnetic structure")

    magnetic_lattice = np.asarray(magnetic.std_lattice, dtype=float)
    magnetic_positions = np.asarray(magnetic.std_positions, dtype=float)
    magnetic_types = np.asarray(magnetic.std_types, dtype=int)
    magnetic_moments = np.asarray(magnetic.std_tensors, dtype=float)
    if len(magnetic_positions) != len(target_positions):
        raise RuntimeError(
            "magnetic and SeeK-path standard cells have different site counts "
            f"({len(magnetic_positions)} != {len(target_positions)})"
        )

    # The standard settings of a magnetic subgroup and its G0 parent normally
    # differ only by axis permutation/reversal plus a rigid Cartesian rotation.
    # Prefer proper transformations so axial moments rotate like ordinary
    # vectors and no artificial reflection of the spin pattern is introduced.
    best = None
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            basis_change = np.zeros((3, 3), dtype=float)
            for row, column in enumerate(permutation):
                basis_change[row, column] = signs[row]
            if np.linalg.det(basis_change) < 0.0:
                continue
            reindexed_lattice = basis_change @ magnetic_lattice
            rotation = np.linalg.solve(reindexed_lattice, target_lattice)
            if np.linalg.det(rotation) < 0.0:
                continue
            if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-5):
                continue
            reindexed_positions = np.mod(
                magnetic_positions @ np.linalg.inv(basis_change), 1.0
            )
            matched = _match_periodic_standard_sites(
                target_positions,
                target_types,
                reindexed_positions,
                magnetic_types,
                target_lattice,
            )
            if matched is None:
                continue
            mapping, max_distance = matched
            if best is None or max_distance < best[0]:
                best = (max_distance, mapping, rotation)
    if best is None:
        raise RuntimeError(
            "could not align the magnetic standard cell with the SeeK-path cell"
        )

    max_distance, mapping, rotation = best
    if max_distance > max(1e-5, 5.0 * tolerance):
        raise RuntimeError(
            "magnetic-to-SeeK-path site mismatch is too large "
            f"({max_distance:.3g} A)"
        )
    aligned_moments = (magnetic_moments @ rotation)[mapping]
    target_elements = [str(site.specie) for site in target]
    _write_magnetic_mcif(
        output_path,
        title,
        target_lattice,
        target_elements,
        target_positions,
        aligned_moments,
    )
    return {
        "path": os.path.normpath(output_path),
        "max_site_mismatch": max_distance,
    }


def _load_magnetic_input_data(structure_file, moments_str, spin_axis_cart):
    is_mcif = structure_file.lower().endswith(".mcif")
    if is_mcif:
        from pymatgen.io.cif import CifParser
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = CifParser(structure_file).parse_structures(primitive=False)[0]
        lattice = np.array(structure.lattice.matrix, dtype=float)
        positions = np.array([site.frac_coords for site in structure], dtype=float)
        elements = [str(site.specie) for site in structure]
        moments = np.array([
            np.array(site.properties["magmom"].moment, dtype=float)
            if "magmom" in site.properties else np.zeros(3)
            for site in structure
        ])
        return lattice, positions, elements, moments, "cartesian"

    from ase.io import read
    structure = read(structure_file)
    lattice = np.array(structure.get_cell(), dtype=float)
    positions = np.array(structure.get_scaled_positions(), dtype=float)
    elements = structure.get_chemical_symbols()
    axis = parse_cartesian_spin_axis(spin_axis_cart)
    scalars = parse_magmoms(moments_str) if moments_str else []
    scalars = fit_magmoms_to_structure(scalars, len(elements))
    moments = np.asarray(scalars, dtype=float)[:, None] * axis[None, :]
    return lattice, positions, elements, moments, "cartesian"


def write_bandplot_lattice_config(lattice_type, filename="alterband.toml"):
    """Record the detected lattice type for later band plotting."""
    if not lattice_type:
        return

    line = f'lattice_type = "{lattice_type}"\n'
    try:
        if os.path.exists(filename):
            with open(filename, "r", encoding="utf-8-sig") as f:
                lines = f.readlines()
            updated = False
            new_lines = []
            for existing in lines:
                stripped = existing.split("#", 1)[0].strip()
                if stripped.startswith("lattice_type") and "=" in stripped:
                    new_lines.append(line)
                    updated = True
                else:
                    new_lines.append(existing)
            if not updated:
                if new_lines and not new_lines[-1].endswith("\n"):
                    new_lines[-1] += "\n"
                new_lines.append(line)
            _atomic_write_text(filename, "".join(new_lines))
        else:
            _atomic_write_text(
                filename,
                "# AlterSeeK band-plot settings\n"
                "# Reads: KLABELS, REFORMATTED_BAND_UP.dat, "
                "REFORMATTED_BAND_DW.dat\n"
                + line,
            )
        print(f"Band plot config updated: {filename} ({line.strip()})")
    except Exception as exc:
        print(f"[Warning] Could not update band plot config '{filename}': {exc}")


def write_qe_bandplot_config(filename="alterband_qe.toml"):
    """Create the QE band-plot config file with a header comment, if missing."""
    if os.path.exists(filename):
        return
    try:
        _atomic_write_text(
            filename,
            "# AlterSeeK QE band-plot settings\n"
            "# Reads: band_up.gnu, band_down.gnu, KPOINTS_alter_qe\n",
        )
        print(f"Band plot config created: {filename}")
    except Exception as exc:
        print(f"[Warning] Could not create band plot config '{filename}': {exc}")
