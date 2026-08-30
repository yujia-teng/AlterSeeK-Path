"""Structure and file I/O: reading input structures (POSCAR/.vasp, .cif,
.mcif, and any other ASE-readable format), writing helper POSCAR/mcif files,
lattice geometry utilities, atomic text writes, and band-plot config writing.

pymatgen/ASE imports are function-local to keep module import lightweight.
"""
import os
import json
import warnings
import numpy as np

from .atomic_write import _atomic_write_text, _atomic_open_text
from .mcif import _validate_collinear_moments

from .find_sf_operations import (
    fit_magmoms_to_structure,
    parse_cartesian_spin_axis,
    parse_magmoms,
)

DEFAULT_PLOT_CONFIG = "alterseek_plot_vasp.toml"
DEFAULT_PLOT_CONFIG_QE = "alterseek_plot_qe.toml"
DEFAULT_PLOT_CONFIG_ABINIT = "alterseek_plot_abinit.toml"


def _group_poscar_sites(elements, positions, moment_keys=None):
    groups = []
    key_to_group = {}
    for idx, element in enumerate(elements):
        key = (element, moment_keys[idx]) if moment_keys is not None else (element,)
        if key not in key_to_group:
            key_to_group[key] = len(groups)
            groups.append((element, []))
        groups[key_to_group[key]][1].append(idx)
    ordered_indices = [idx for _, indices in groups for idx in indices]
    symbols = [symbol for symbol, _ in groups]
    counts = [len(indices) for _, indices in groups]
    ordered_positions = [positions[idx] for idx in ordered_indices]
    return symbols, counts, ordered_positions, ordered_indices


def _write_poscar(path, title, lattice, symbols, counts, positions):
    with _atomic_open_text(path) as f:
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
    # Import pymatgen locally to keep module import light while retaining its handling of scale factors, Cartesian/Direct coordinates, and Selective Dynamics.
    from pymatgen.io.vasp.inputs import Poscar

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="We strongly encourage explicit.*encoding",
        )
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
    removed = set(species_to_remove)
    kept = [
        (element, position)
        for element, position in zip(elements, positions)
        if element not in removed
    ]
    if not kept:
        raise RuntimeError(f"No atoms left after removing {species_to_remove} from {source_path}.")
    kept_elements, kept_positions = zip(*kept)
    _write_poscar_from_sites(target_path, title, lattice, list(kept_elements), list(kept_positions))


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

    with _atomic_open_text(path) as f:
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
        _validate_collinear_moments(moments)
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


def write_bandplot_lattice_config(lattice_type, filename=DEFAULT_PLOT_CONFIG):
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
                "# AlterSeeK VASP band-plot settings\n"
                "# Commented values are defaults; uncomment them to override.\n\n"
                + line
                + "\n"
                "# klabels = \"KLABELS\"\n"
                "# up = \"REFORMATTED_BAND_UP.dat\"\n"
                "# down = \"REFORMATTED_BAND_DW.dat\"\n"
                "# emin = -2\n"
                "# emax = 2\n"
                "# fig_width = 12\n"
                "# fig_height = 5\n"
                "# gap_width_inches = 0.05\n"
                "# split_panels = 1\n"
                "# output = \"alterband.png\"\n"
                "# save_pdf = false\n",
            )
        print(f"Band plot config updated: {filename} ({line.strip()})")
    except Exception as exc:
        print(f"[Warning] Could not update band plot config '{filename}': {exc}")


def write_qe_bandplot_config(filename=DEFAULT_PLOT_CONFIG_QE):
    """Create the QE band-plot config file with a header comment, if missing."""
    if os.path.exists(filename):
        return
    try:
        _atomic_write_text(
            filename,
            "# AlterSeeK QE band-plot settings\n"
            "# Commented values are defaults; uncomment them to override.\n\n"
            "# band_up = \"band_up.gnu\"  # replace with the spin-up bands.x .gnu file\n"
            "# band_down = \"band_down.gnu\"  # replace with the spin-down bands.x .gnu file\n"
            "# kpoints_qe = \"KPOINTS_alter_qe\"\n"
            "# fermi_ev = 0.0  # replace with the Fermi energy from the QE output\n"
            "# emin = -2\n"
            "# emax = 2\n"
            "# fig_width = 12\n"
            "# fig_height = 5\n"
            "# gap_width_inches = 0.05\n"
            "# split_panels = 1\n"
            "# output = \"alterband_qe.png\"\n"
            "# save_pdf = false\n",
        )
        print(f"Band plot config created: {filename}")
    except Exception as exc:
        print(f"[Warning] Could not create band plot config '{filename}': {exc}")


def write_abinit_bandplot_config(
    structure_file=None, filename=DEFAULT_PLOT_CONFIG_ABINIT
):
    """Create or update the ABINIT band-plot config with the structure file."""
    structure_line = (
        f"structure = {json.dumps(os.fspath(structure_file), ensure_ascii=False)}\n"
        if structure_file is not None else None
    )
    try:
        if os.path.exists(filename):
            if structure_line is None:
                return
            with open(filename, "r", encoding="utf-8-sig") as f:
                lines = f.readlines()
            updated = False
            new_lines = []
            for existing in lines:
                stripped = existing.split("#", 1)[0].strip()
                if stripped.startswith("structure") and "=" in stripped:
                    new_lines.append(structure_line)
                    updated = True
                else:
                    new_lines.append(existing)
            if not updated:
                new_lines.append(structure_line)
            _atomic_write_text(filename, "".join(new_lines))
            action = "updated"
        else:
            _atomic_write_text(
                filename,
                "# AlterSeeK ABINIT band-plot settings\n"
                "# Commented values are defaults; uncomment them to override.\n\n"
                + (structure_line or "")
                + "\n"
                "# eig = \"EIG\"  # replace with the ABINIT _EIG file\n"
                "# kpoints_abinit = \"KPOINTS_alter_abinit\"\n"
                "# abo = \"abo\"  # replace with the ABINIT .abo output\n"
                "# fermi_ev = 0.0  # use only when the .abo file is unavailable\n"
                "# emin = -2\n"
                "# emax = 2\n"
                "# fig_width = 12\n"
                "# fig_height = 5\n"
                "# gap_width_inches = 0.05\n"
                "# split_panels = 1\n"
                "# output = \"alterband_abinit.png\"\n"
                "# save_pdf = false\n",
            )
            action = "created"
        print(f"Band plot config {action}: {filename}")
    except Exception as exc:
        print(
            f"[Warning] Could not create or update band plot config "
            f"'{filename}': {exc}"
        )
