# Script for generating general k-path for band structure calculation
# Author: Yujia Teng, Rutgers University
# 06/14/2025 - Creation of the script
# 01/12/2026 - Automatically reads the spin-flip operation based on outputs of find_sf_operations.py
# 03/2026    - Integrated find_sf_operations and compute_centroid_hybrid into one workflow
import numpy as np
import os
import shutil
import sys
from typing import List, Optional

# Try to import spin-flip operations finder
try:
    from find_sf_operations import run as find_sf_run
    from find_sf_operations import parse_cartesian_spin_axis, parse_magmoms
    FIND_SF_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] find_sf_operations unavailable ({_exc}); "
          "Step 0 spin-operation detection is disabled.")
    FIND_SF_AVAILABLE = False
    parse_cartesian_spin_axis = None
    parse_magmoms = None

try:
    from findspingroup import (
        find_spin_group_acc_primitive,
        find_spin_group_acc_primitive_from_data,
    )
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] findspingroup unavailable ({_exc}); --ssg-setting is disabled. "
          "Install with: pip install \"findspingroup>=0.15.6,<0.16\"")
    find_spin_group_acc_primitive = None
    find_spin_group_acc_primitive_from_data = None
    FIND_SG_MAGNETIC_SETTING_AVAILABLE = False

# Try to import IBZ centroid calculator
try:
    from compute_centroid_hybrid import (run as compute_centroid,
                                         no_altermagnetism_reason,
                                         plot_spin_flip_figure,
                                         plot_spin_bz_figure,
                                         plot_spin_bz_top_view_figure)
    import matplotlib.pyplot as plt
    CENTROID_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] compute_centroid_hybrid/matplotlib unavailable ({_exc}); "
          "centroid and figure generation are disabled.")
    CENTROID_AVAILABLE = False
    no_altermagnetism_reason = None
    plot_spin_flip_figure = None
    plot_spin_bz_figure = None
    plot_spin_bz_top_view_figure = None
    plt = None


STEP0_VERBOSE_SUMMARY = False


_MARKER_SEEDS = [
    np.array([0.11000000, 0.12000000, 0.15000001]),
    np.array([0.13000000, 0.17000000, 0.23000001]),
    np.array([0.07100000, 0.19300000, 0.31700001]),
    np.array([0.21100000, 0.13700000, 0.29300001]),
]


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
    with open(path, "w") as f:
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


def _marker_orbit(seed, operations):
    positions = []
    for op in operations:
        rotation = np.array(op["real_rotation"], dtype=float)
        translation = np.array(op.get("translation", [0.0, 0.0, 0.0]), dtype=float)
        positions.append(seed @ rotation.T + translation)
    return _dedupe_frac_positions(positions)


def _build_marker_helper(lattice, symbols, counts, positions, operations):
    best = None
    for seed in _MARKER_SEEDS:
        markers = _marker_orbit(seed, operations)
        helper_positions = list(positions) + markers
        candidate = {
            "seed": seed,
            "symbols": list(symbols) + ["He"],
            "counts": list(counts) + [len(markers)],
            "positions": helper_positions,
            "markers": markers,
            "min_distance": _min_periodic_cart_distance(helper_positions, lattice),
        }
        if best is None or candidate["min_distance"] > best["min_distance"]:
            best = candidate
    if best is None:
        raise RuntimeError("Could not build SSG-setting marker helper.")
    return best


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

    with open(path, "w") as f:
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
        return lattice, positions, elements, moments, "cartesian"

    from ase.io import read
    structure = read(structure_file)
    lattice = np.array(structure.get_cell(), dtype=float)
    positions = np.array(structure.get_scaled_positions(), dtype=float)
    elements = structure.get_chemical_symbols()
    if parse_cartesian_spin_axis is None or parse_magmoms is None:
        raise RuntimeError("find_sf_operations magnetic input parsers are not available.")
    axis = parse_cartesian_spin_axis(spin_axis_cart)
    scalars = parse_magmoms(moments_str) if moments_str else []
    if len(scalars) < len(elements):
        scalars.extend([0.0] * (len(elements) - len(scalars)))
    elif len(scalars) > len(elements):
        scalars = scalars[:len(elements)]
    moments = np.asarray(scalars, dtype=float)[:, None] * axis[None, :]
    return lattice, positions, elements, moments, "cartesian"


def _spin_axis_from_moments(moments):
    for moment in np.asarray(moments, dtype=float):
        norm = np.linalg.norm(moment)
        if norm > 1e-8:
            return moment / norm
    raise ValueError("No nonzero magnetic moment found.")


def _operation_class_indices(operations, spin_axis, flip):
    indices = []
    for i, op in enumerate(operations):
        spin_rotation = np.array(op["spin_rotation"], dtype=float)
        mapped_axis = spin_rotation @ spin_axis
        is_flip = np.allclose(mapped_axis, -spin_axis, atol=1e-7)
        is_preserve = np.allclose(mapped_axis, spin_axis, atol=1e-7)
        if not is_flip and not is_preserve:
            continue
        if is_flip == flip:
            indices.append(i)
    return indices


def _collect_point_ops_from_payload(operations, indices, include_inversion=True):
    point_ops = []
    source_indices = []
    for idx in indices:
        rotation = np.array(operations[idx]["real_rotation"], dtype=float)
        variants = [rotation]
        if include_inversion:
            variants.append(-rotation)
        for candidate in variants:
            rounded = np.rint(candidate).astype(int)
            compare = rounded if np.allclose(candidate, rounded, atol=1e-7) else candidate
            if not any(np.allclose(compare, existing, atol=1e-7) for existing in point_ops):
                point_ops.append(compare)
                source_indices.append(int(operations[idx].get("index", idx + 1)))
    return point_ops, source_indices


def _write_operation_file(filename, rotations, source_indices, label):
    with open(filename, "w") as f:
        f.write(f"# Found {len(rotations)} inversion-extended spin-{label} point operations\n")
        f.write(f"# Original Indices: {source_indices}\n")
        for i, rotation in enumerate(rotations):
            f.write(f"Operation_{i + 1}\n")
            for row in rotation:
                f.write(" ".join(f"{value:.10g}" for value in row) + "\n")
            f.write("\n")
    return len(rotations)


def _write_magnetic_setting_operation_files(operations, spin_axis, output_dir="."):
    flip_indices = _operation_class_indices(operations, spin_axis, flip=True)
    preserve_indices = _operation_class_indices(operations, spin_axis, flip=False)
    flip_ops, flip_sources = _collect_point_ops_from_payload(operations, flip_indices)
    preserve_ops, preserve_sources = _collect_point_ops_from_payload(operations, preserve_indices)
    flip_count = _write_operation_file(
        os.path.join(output_dir, "spin_flip_operations.txt"),
        flip_ops,
        flip_sources,
        "flipping",
    )
    preserve_count = _write_operation_file(
        os.path.join(output_dir, "spin_preserve_operations.txt"),
        preserve_ops,
        preserve_sources,
        "preserving",
    )
    return flip_count, preserve_count


def _magnetic_primitive_ssg_operations(result):
    views = (
        result.get("operation_views", {})
        .get("magnetic_primitive_cartesian", {})
        .get("views", {})
    )
    nssg_view = views.get("nssg")
    if nssg_view and nssg_view.get("ops"):
        return nssg_view["ops"]
    all_view = views.get("all")
    if all_view and all_view.get("ops"):
        return all_view["ops"]
    operations = result.get("acc_primitive_ssg_ops_cartesian")
    if operations is None:
        operations = result.get("acc_primitive_ssg_operation_matrices", [])
    return operations


def prepare_magnetic_setting_files(structure_file, moments_str="", spin_axis_cart=None, output_dir="."):
    """Write real magnetic primitive POSCAR/MCIF files from FindSpinGroup."""
    if not FIND_SG_MAGNETIC_SETTING_AVAILABLE:
        raise RuntimeError("FindSpinGroup accurate primitive API is not available.")

    is_mcif = structure_file.lower().endswith(".mcif")
    if is_mcif:
        result = find_spin_group_acc_primitive(structure_file)
    else:
        lattice_in, positions_in, elements_in, moments_in, spin_setting = _load_magnetic_input_data(
            structure_file, moments_str, spin_axis_cart
        )
        result = find_spin_group_acc_primitive_from_data(
            structure_file,
            lattice_in,
            positions_in,
            elements_in,
            [1.0] * len(elements_in),
            moments_in,
            input_spin_setting=spin_setting,
        )
    cell = result["acc_primitive_cell_detail"]
    lattice = np.array(cell["lattice"], dtype=float)
    positions = [np.array(pos, dtype=float) for pos in cell["positions"]]
    elements = [str(el) for el in cell["elements"]]
    moments = np.array(cell.get("moments", np.zeros((len(elements), 3))), dtype=float)

    basename = os.path.splitext(os.path.basename(structure_file))[0]
    temp_dir = os.path.join(output_dir, ".alterseek_ssgstd_tmp")
    os.makedirs(temp_dir, exist_ok=True)
    real_path = os.path.join(temp_dir, f"{basename}_ssgprim.vasp")
    mcif_path = os.path.join(output_dir, f"{basename}_ssgprim.mcif")
    magmom_path = os.path.join(temp_dir, f"{basename}_ssgprim_MAGMOM.txt")
    helper_path = os.path.join(temp_dir, f"{basename}_ssgstd.vasp")

    symbols, counts, grouped_positions, ordered_indices = _group_poscar_sites(
        elements, positions)
    _write_poscar(
        real_path,
        f"{basename} magnetic primitive from FindSpinGroup",
        lattice,
        symbols,
        counts,
        grouped_positions,
    )

    ordered_moments = moments[ordered_indices]
    ordered_elements = [elements[idx] for idx in ordered_indices]
    _write_magnetic_mcif(
        mcif_path,
        f"{basename}_ssgprim",
        lattice,
        ordered_elements,
        grouped_positions,
        ordered_moments,
    )
    with open(magmom_path, "w") as f:
        f.write("# Magnetic primitive POSCAR atom order matches:\n")
        f.write(f"# {real_path}\n")
        f.write("# Vector moments from FindSpinGroup acc_primitive_cell_detail:\n")
        for moment in ordered_moments:
            f.write(f"{moment[0]: .10f} {moment[1]: .10f} {moment[2]: .10f}\n")
        axis = None
        for moment in ordered_moments:
            norm = np.linalg.norm(moment)
            if norm > 1e-8:
                axis = moment / norm
                break
        if axis is not None:
            scalars = [float(np.dot(moment, axis)) for moment in ordered_moments]
            f.write("# Collinear scalar MAGMOM along first nonzero moment axis:\n")
            f.write("MAGMOM = " + " ".join(f"{value:.8g}" for value in scalars) + "\n")

    operations = _magnetic_primitive_ssg_operations(result)
    marker_helper = _build_marker_helper(
        lattice,
        symbols,
        counts,
        grouped_positions,
        operations,
    )
    _write_poscar(
        helper_path,
        f"{basename} SSG-setting marker helper from FindSpinGroup",
        lattice,
        marker_helper["symbols"],
        marker_helper["counts"],
        marker_helper["positions"],
    )
    spin_axis = _spin_axis_from_moments(moments)
    flip_count, preserve_count = _write_magnetic_setting_operation_files(
        operations, spin_axis, output_dir=output_dir
    )

    return {
        "real_poscar_path": real_path,
        "helper_path": helper_path,
        "mcif_path": mcif_path,
        "magmom_path": magmom_path,
        "temp_dir": temp_dir,
        "basename": basename,
        "seekpath_type_numbers": None,
        "spin_flip_operations": flip_count,
        "spin_preserve_operations": preserve_count,
        "summary": {
            "index": result.get("index"),
            "acc_symbol": result.get("acc_symbol"),
            "setting": result.get("acc_primitive_cell_setting"),
            "marker_seed": marker_helper["seed"].tolist(),
            "marker_count": len(marker_helper["markers"]),
            "marker_min_distance": marker_helper["min_distance"],
        },
    }


def finalize_magnetic_setting_outputs(
    mag_setting,
    centroid_result,
    output_dir=".",
    verbose_output=False,
):
    helper_source = centroid_result.get("standardized_structure_path")
    if not helper_source or not os.path.exists(helper_source):
        return {}

    basename = mag_setting["basename"]
    real_final = os.path.join(
        output_dir, f"{basename}_ssgstd.vasp"
    )

    _write_without_species(
        helper_source,
        real_final,
        {"He"},
        f"{basename} magnetic setting standardized",
    )
    helper_final = None
    if verbose_output:
        helper_final = os.path.join(
            output_dir, f"{basename}_ssgstd_helper.vasp"
        )
        if os.path.abspath(helper_source) != os.path.abspath(helper_final):
            shutil.copyfile(helper_source, helper_final)

    # Remove or relocate low-level seekpath artifacts for the hidden marker helper.
    # The clean final standardized structure above is the user-facing record.
    temp_dir = mag_setting.get("temp_dir")
    for path in (
        helper_source,
        centroid_result.get("standard_mapping_path"),
    ):
        if not path or not os.path.exists(path):
            continue
        try:
            if verbose_output and temp_dir and os.path.isdir(temp_dir):
                shutil.move(path, os.path.join(temp_dir, os.path.basename(path)))
            else:
                os.remove(path)
        except OSError:
            pass

    if not verbose_output and temp_dir and os.path.isdir(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except OSError:
            pass

    result = {
        "standard_real_path": real_final,
        "b_matrix_output": _reciprocal_from_poscar(real_final),
    }
    if helper_final:
        result["standard_with_helper_path"] = helper_final
    if verbose_output:
        result["intermediate_dir"] = mag_setting.get("temp_dir")
    return result


def write_bandplot_lattice_config(lattice_type, filename="alterband.toml"):
    """Record the detected lattice type for later band plotting."""
    if not lattice_type:
        return

    line = f'lattice_type = "{lattice_type}"\n'
    try:
        if os.path.exists(filename):
            with open(filename, "r") as f:
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
            with open(filename, "w") as f:
                f.writelines(new_lines)
        else:
            with open(filename, "w") as f:
                f.write("# AlterSeeK band-plot settings\n")
                f.write(line)
        print(f"Band plot config updated: {filename} ({line.strip()})")
    except Exception as exc:
        print(f"[Warning] Could not update band plot config '{filename}': {exc}")


class KPointsModifier:
    def __init__(self, magnetic_setting: bool = False, output_verbose: bool = False):
        self.kpoints_data = []
        self.header_lines = []
        self.extra_general_points = []
        self.kpoints_basis_matrix = None
        self.output_basis_matrix = None
        self.kpoints_basis_rotation = None
        self.magnetic_setting = magnetic_setting
        self.output_verbose = output_verbose

    @staticmethod
    def _display_label(label: str) -> str:
        # Console display currently matches the VASP-safe form; delegate so the
        # two normalizations can never drift apart.
        return KPointsModifier._kpoints_label(label)

    @staticmethod
    def _kpoints_label(label: str) -> str:
        """Return labels in a VASP-safe form for KPOINTS files."""
        label = str(label)
        return 'GAMMA' if label.strip().upper() == 'GAMMA' or label == '\u0393' else label

    @classmethod
    def _format_path(cls, path_segments) -> str:
        parts = []
        prev_end = None
        for seg_start, seg_end in path_segments:
            start = cls._kpoints_label(seg_start)
            end = cls._kpoints_label(seg_end)
            if seg_start != prev_end:
                parts.append(f"| {start}-{end}" if prev_end else f"{start}-{end}")
            else:
                parts.append(end)
            prev_end = seg_end
        return "-".join(parts).replace("-| ", " | ")

    @staticmethod
    def _format_matrix(op: np.ndarray) -> str:
        rows = []
        for row in op:
            row_str = " ".join(
                f"{int(x): 3d}" if float(x) % 1 == 0 else f"{x: .2f}"
                for x in row
            )
            rows.append(f"    [ {row_str} ]")
        return "\n".join(rows)

    @staticmethod
    def _count_written_segments(kpoints) -> int:
        count = 0
        for i in range(len(kpoints) - 1):
            start_point = kpoints[i]
            end_point = kpoints[i + 1]
            if start_point is None or end_point is None:
                continue
            if start_point[3] == end_point[3]:
                continue
            if {start_point[3], end_point[3]} == {"k", "k'"}:
                continue
            count += 1
        return count
    
    def read_kpoints_file(self, filename: str = "KPOINTS") -> bool:
        """Read a line-mode KPOINTS file (e.g. generated by VASPKIT)."""
        try:
            # pymatgen validates the format; import locally to keep module
            # import light.
            from pymatgen.io.vasp.inputs import Kpoints

            kpoints = Kpoints.from_file(filename)
            if kpoints.style.name != "Line_mode":
                print(f"Error: {filename} is not a line-mode KPOINTS file "
                      f"(style: {kpoints.style.name}).")
                return False
            coord_type = kpoints.coord_type or "Reciprocal"
            if not coord_type.lower().startswith("r"):
                print(f"Error: {filename} uses {coord_type} coordinates; only "
                      "Reciprocal line-mode KPOINTS are supported.")
                return False

            with open(filename, 'r') as f:
                lines = f.readlines()
            self.header_lines = [line.strip() for line in lines[:4]]

            self.kpoints_data = []
            self.kpoints_basis_matrix = None
            self.output_basis_matrix = None
            self.kpoints_basis_rotation = None
            for coords, label in zip(kpoints.kpts, kpoints.labels or []):
                label = (label or "").strip()
                if not label:  # Skip unlabeled rows (matches prior behavior)
                    continue
                self.kpoints_data.append(
                    [float(coords[0]), float(coords[1]), float(coords[2]), label])

            print(f"Successfully read {len(self.kpoints_data)} k-points from {filename}")
            return True

        except FileNotFoundError:
            print(f"Error: File {filename} not found!")
            return False
        except Exception as e:
            print(f"Error reading file: {e}")
            return False

    def _kpoint_for_output_basis(self, point: List) -> List:
        """Convert an internal k-point to the POSCAR reciprocal basis for VASP."""
        if self.kpoints_basis_matrix is None or self.output_basis_matrix is None:
            return point

        try:
            k_frac = np.array(point[:3], dtype=float)
            b_kpoints = np.array(self.kpoints_basis_matrix, dtype=float)
            if self.kpoints_basis_rotation is not None:
                b_kpoints = b_kpoints @ np.array(self.kpoints_basis_rotation, dtype=float)
            b_output = np.array(self.output_basis_matrix, dtype=float)
            k_out = k_frac @ b_kpoints @ np.linalg.inv(b_output)
        except Exception as exc:
            # Writing unconverted coordinates would silently reproduce the
            # oI3/221-P-d conventional-cell mismatch; stop instead.
            raise RuntimeError(
                f"Output-basis conversion failed for k-point '{point[3]}': {exc}. "
                "Refusing to write unconverted coordinates into KPOINTS."
            ) from exc
        return [k_out[0], k_out[1], k_out[2], point[3]]

    def load_flip_operations(self, filename: str = "spin_flip_operations.txt") -> List[np.ndarray]:
        """Reads pre-calculated rotation matrices from file"""
        matrices = []
        unique = []
        current_matrix = []
        if not os.path.exists(filename):
            legacy = (
                "preserve_spin_operations.txt"
                if filename == "spin_preserve_operations.txt"
                else "flip_spin_operations.txt"
            )
            if filename != legacy and os.path.exists(legacy):
                filename = legacy
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Operation"):
                    continue
                if "|" in line:
                    line = line.split("|", 1)[0].strip()
                
                parts = line.split()
                if len(parts) == 3:
                    current_matrix.append([float(x) for x in parts])
                
                if len(current_matrix) == 3:
                    mat = np.array(current_matrix)
                    if not any(np.allclose(mat, ex, atol=1e-8) for ex in unique):
                        unique.append(mat)
                        matrices.append(mat)
                    current_matrix = []
            return matrices
        except FileNotFoundError:
            return []

    def load_preserve_operations(self, filename: str = "spin_preserve_operations.txt") -> List[np.ndarray]:
        """Reads pre-calculated spin-preserve rotation matrices from file."""
        return self.load_flip_operations(filename)

    def transform_kpoint(self, kpoint: List[float], transformation_matrix: np.ndarray) -> List[float]:
        """Apply transformation matrix to k-point: R^(-1)^T * k"""
        k = np.array(kpoint[:3])  # Take only x, y, z coordinates
        # Calculate R^(-1)^T
        R_inv_T = np.linalg.inv(transformation_matrix).T
        k_transformed = R_inv_T @ k
        return k_transformed.tolist()
    
    def insert_general_kpoints(self,
                               general_kpoint: List[float],
                               transformation_matrix: np.ndarray,
                               extra_general_points: Optional[List[List]] = None) -> List[List]:
        """
        Insert general k-points into every segment of the high symmetry path.
        Per-segment butterfly: for each (A,B) segment:
          - If A not yet opened: A-k | k'-A'
          - Connect A' to B'
          - Close: B'-k' | k-B
        Points already butterflied in earlier chains appear as plain segments.
        """
        if not self.kpoints_data:
            print("Error: No k-points data loaded. Please read KPOINTS file first.")
            return []

        k_prime = self.transform_kpoint(general_kpoint, transformation_matrix)
        kpt = general_kpoint
        kp  = k_prime

        def coords_eq(p, q, tol=1e-6):
            return abs(p[0]-q[0]) < tol and abs(p[1]-q[1]) < tol and abs(p[2]-q[2]) < tol

        def pt_key(p):
            return (round(p[0], 6), round(p[1], 6), round(p[2], 6))

        def is_gamma(p):
            label = str(p[3])
            return label.strip().upper() == 'GAMMA' or label == '\u0393'

        def get_prime(p):
            """Return primed version of p (Gamma stays unprimed)."""
            if is_gamma(p):
                return p.copy()
            tc = self.transform_kpoint(p, transformation_matrix)
            return [tc[0], tc[1], tc[2], f"{p[3]}'"]

        # --- Step 1: group flat kpoints_data into segment pairs ---
        raw = self.kpoints_data
        seg_pairs = [(raw[i], raw[i+1]) for i in range(0, len(raw) - 1, 2)]

        # --- Step 2: build connected chains ---
        chains = []
        current_chain = [seg_pairs[0][0], seg_pairs[0][1]]
        for sp_start, sp_end in seg_pairs[1:]:
            if coords_eq(current_chain[-1], sp_start):
                current_chain.append(sp_end)
            else:
                chains.append(current_chain)
                current_chain = [sp_start, sp_end]
        chains.append(current_chain)

        # --- Step 3: alternating plain / butterfly segments ---
        #
        # Pattern:
        #   Even-indexed segment --plain  (emit A, B)
        #   Odd-indexed segment  --butterfly (emit A, k, k', A', B', k', k, B)
        #
        # First chain: parity 0 --segment 0 is plain, segment 1 is butterfly, ...
        # Other chains: parity 1 --segment 0 is butterfly, segment 1 is plain, ...
        #
        # Consecutive segments share an endpoint; the duplicate is suppressed by
        # the same-label skip in write_kpoints_file and the display builder.
        # GAMMA is self-conjugate: get_prime(GAMMA) = GAMMA.
        path_sequence = []

        def emit_plain(A, B):
            path_sequence.append(A.copy())
            path_sequence.append(B.copy())

        def emit_butterfly(A, B):
            path_sequence.append(A.copy())
            path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
            path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
            path_sequence.append(get_prime(A))
            path_sequence.append(get_prime(B))
            # Universal rule: skip close if B already received butterfly treatment.
            # Applies to all points equally (GAMMA, X, W, etc.).
            if pt_key(B) not in butterflied:
                path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                path_sequence.append(B.copy())

        butterflied = set()  # pt_key of points that have received butterfly treatment

        for ci, chain in enumerate(chains):
            # Dedup within chain (by coordinates, not labels)
            unique = []
            for pt in chain:
                if not unique or not coords_eq(unique[-1], pt):
                    unique.append(pt)

            # Degenerate chain: all points share the same coordinates
            # (e.g. U_0 --T for certain lattice parameters). Skip silently.
            if len(unique) < 2:
                labels = [p[3] for p in chain]
                print(f"  [Note] Part {ci+1} ({' - '.join(labels)}) skipped: "
                      f"endpoints coincide in coordinates")
                continue

            if ci > 0:
                path_sequence.append(None)

            start_parity = 0 if ci == 0 else 1

            parity = start_parity  # tracks alternation independently of s
            for s, (A, B) in enumerate(zip(unique, unique[1:])):
                A_key = pt_key(A)
                B_key = pt_key(B)
                is_last = (s == len(unique) - 2)
                # If A was already butterflied, continue from the side that is
                # already active. If the previous butterfly ended at A', the next
                # plain segment must be A'->B, not A'->A->B.
                if A_key in butterflied:
                    A_start = A
                    if path_sequence:
                        prev_pt = path_sequence[-1]
                        A_prime = get_prime(A)
                        if (prev_pt is not None and prev_pt[3] == A_prime[3]
                                and A_prime[3] != A[3]):
                            A_start = A_prime
                    emit_plain(A_start, B)
                    # B is a new endpoint --partial butterfly B-k|k'-B'
                    if B_key not in butterflied and is_last:
                        path_sequence.append(B.copy())
                        path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                        path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                        path_sequence.append(get_prime(B))
                        butterflied.add(B_key)
                    # Only hold parity when butterflied-A overrides a butterfly slot;
                    # if parity is already even (plain slot), advance normally.
                    if parity % 2 == 0:
                        parity += 1
                elif parity % 2 == 1:
                    emit_butterfly(A, B)
                    butterflied.add(A_key)
                    butterflied.add(B_key)
                    parity += 1
                else:
                    emit_plain(A, B)
                    parity += 1

            # After processing all segments in this chain, check if the last
            # point still lacks butterfly treatment (happens when chain 0 has
            # only 1 segment --the single plain pair leaves B un-butterflied).
            last_key = pt_key(unique[-1])
            if last_key not in butterflied:
                last_pt = unique[-1]
                path_sequence.append(last_pt.copy())
                path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                path_sequence.append(get_prime(last_pt))
                butterflied.add(last_key)

        # Remove trailing sentinel
        if path_sequence and path_sequence[-1] is None:
            path_sequence.pop()

        # Doubled-IBZ append-only anchors. These are project-only copied
        # vertices that should be sampled through the general point without
        # adding duplicated high-symmetry edges.
        extra_general_points = extra_general_points or []
        for pt in extra_general_points:
            if path_sequence:
                path_sequence.append(None)
            path_sequence.append(pt.copy())
            path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
            path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
            path_sequence.append(get_prime(pt))

        # Print generated path as label string
        tokens = []  # list of label strings or '|' for breaks
        prev = None
        i = 0
        while i < len(path_sequence) - 1:
            cur = path_sequence[i]
            nxt = path_sequence[i + 1]
            # Chain boundary sentinel --insert break
            if cur is None or nxt is None:
                if tokens and tokens[-1] != '|':
                    tokens.append('|')
                i += 1
                continue
            # k↔k' connection --insert break
            if (cur[3] == "k" and nxt[3] == "k'") or \
               (cur[3] == "k'" and nxt[3] == "k"):
                if tokens and tokens[-1] != '|':
                    tokens.append('|')
                i += 1
                continue
            # same label --skip
            if cur[3] == nxt[3]:
                i += 1
                continue
            if prev != cur[3]:
                tokens.append(cur[3])
            tokens.append(nxt[3])
            prev = nxt[3]
            i += 1
        # Build display string: labels joined by '-', breaks become '|'.
        display_segments = []
        current_segment = []
        for t in tokens:
            if t == '|':
                if current_segment:
                    display_segments.append('-'.join(self._display_label(x) for x in current_segment))
                    current_segment = []
            else:
                current_segment.append(t)
        if current_segment:
            display_segments.append('-'.join(self._display_label(x) for x in current_segment))

        if len(display_segments) > 6:
            display = " | ".join(display_segments[:3] + ["..."] + display_segments[-3:])
        else:
            display = " | ".join(display_segments)
        print(f"Generated path: {display}")

        generated_segments = self._count_written_segments(path_sequence)
        generated_points = sum(1 for pt in path_sequence if pt is not None)
        print(f"Full path: {len(seg_pairs)} original segments -> "
              f"{generated_segments} generated segments, {generated_points} k-points")
        if extra_general_points:
            labels = ", ".join(str(pt[3]) for pt in extra_general_points)
            print(f"Added doubled-IBZ general anchors: {labels}")

        return path_sequence

    def insert_general_kpoint_anchors(self, kpoint: List[float],
                                      extra_general_points: Optional[List[List]] = None) -> List[List]:
        """Keep the ordinary path, then append compact high-symmetry/k comparisons."""
        kpt = [kpoint[0], kpoint[1], kpoint[2], "k"]
        raw = self.kpoints_data
        seg_pairs = [(raw[i], raw[i + 1]) for i in range(0, len(raw) - 1, 2)]
        path_sequence = []
        for idx, (start, end) in enumerate(seg_pairs):
            if idx:
                path_sequence.append(None)
            path_sequence.append(start.copy())
            path_sequence.append(end.copy())

        seen = set()
        anchors = []
        for pt in self.kpoints_data + (extra_general_points or []):
            if pt is None:
                continue
            key = (
                round(float(pt[0]), 8),
                round(float(pt[1]), 8),
                round(float(pt[2]), 8),
                str(pt[3]),
            )
            if key in seen:
                continue
            seen.add(key)
            anchors.append(pt)

        for idx in range(0, len(anchors), 2):
            if path_sequence:
                path_sequence.append(None)
            path_sequence.append(anchors[idx].copy())
            path_sequence.append(kpt.copy())
            if idx + 1 < len(anchors):
                path_sequence.append(anchors[idx + 1].copy())

        pair_count = len(anchors) // 2
        leftover = len(anchors) % 2
        print(
            f"Kept ordinary path and added {pair_count} A-k-B comparison segments"
            f"{' plus 1 A-k tail' if leftover else ''}."
        )
        return path_sequence
    
    def write_kpoints_file(self, new_kpoints: List[List], output_file: str = "KPOINTS_modified",
                           transformation_matrix: Optional[np.ndarray] = None,
                           transformation_label: Optional[str] = None):
        """Write modified KPOINTS file with proper Line-Mode format and discontinuity"""
        try:
            with open(output_file, 'w') as f:
                # Write header 
                if transformation_matrix is not None:
                    flat_matrix = np.array(transformation_matrix).flatten()
                    matrix_str = " ".join(f"{x:.8f}" for x in flat_matrix)
                    label = f" ({transformation_label})" if transformation_label else ""
                    f.write(f"Selected spin-flip operation{label} in input-cell fractional basis: {matrix_str}\n")
                else:
                    f.write(f"{self.header_lines[0]}\n")
                
                f.write("   30\n")    # Number of points between each segment
                f.write(f"{self.header_lines[2]}\n")   # Line-Mode
                f.write(f"{self.header_lines[3]}\n")   # Reciprocal
                
                # Write k-points as segments with special handling for discontinuity
                i = 0
                while i < len(new_kpoints) - 1:
                    start_point = new_kpoints[i]
                    end_point = new_kpoints[i + 1]

                    # Skip chain boundary sentinels
                    if start_point is None or end_point is None:
                        i += 1
                        continue

                    # Skip any k to k' connections (there should be none)
                    if start_point[3] == "k" and end_point[3] == "k'":
                        i += 1
                        continue

                    if start_point[3] == "k'" and end_point[3] == "k":
                        i += 1
                        continue

                    if start_point[3] == end_point[3]:
                        i+=1
                        continue

                    start_point_out = self._kpoint_for_output_basis(start_point)
                    end_point_out = self._kpoint_for_output_basis(end_point)
                    
                    # Write segment: start_point -> end_point
                    start_label = self._kpoints_label(start_point_out[3])
                    end_label = self._kpoints_label(end_point_out[3])
                    f.write(f"   {start_point_out[0]:.10f}   {start_point_out[1]:.10f}   {start_point_out[2]:.10f}     {start_label}\n")
                    f.write(f"   {end_point_out[0]:.10f}   {end_point_out[1]:.10f}   {end_point_out[2]:.10f}     {end_label}\n")
                    
                    # Check if this creates a discontinuity (k is dead end)
                    if end_point[3] == "k":
                        f.write("\n")
                    else:
                        if i < len(new_kpoints) - 2:
                            f.write("\n")
                    i += 1
            
            print(f"Modified KPOINTS file written to: {output_file}")
            return True
            
        except Exception as e:
            print(f"Error writing file: {e}")
            return False

    def interactive_modify(self):
        """Interactive modification of KPOINTS file"""
        BOLD  = "\033[1m"
        RESET = "\033[0m"
        print("=== Altermagnetic K-Path Generator ===")

        # Step 0: Compute spin-flip operations from structure
        print(f"\n{BOLD}>>> Step 0: Spin symmetry{RESET}")
        print("Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): ", end='', flush=True)
        struct_file = input().strip()
        if not struct_file: struct_file = "POSCAR"

        # None = Step 0 not run; True = file freshly written; False = ran but no flip ops found
        _step0_wrote_flip_file = None
        _flip_file = 'spin_flip_operations.txt'
        standard_path_reason = None
        standard_path_reason_reported = False
        centroid_result = None
        centroid_error = None
        centroid_struct_file = struct_file
        centroid_seekpath_type_numbers = None
        display_figures = []
        self.extra_general_points = []

        if not os.path.exists(struct_file):
            print(f"[Note] '{struct_file}' not found. Skipping Step 0, using existing spin_flip_operations.txt")
            struct_file = None
            centroid_struct_file = None
        elif not FIND_SF_AVAILABLE:
            print("[Note] find_sf_operations.py not found. Skipping Step 0.")
            struct_file = None
            centroid_struct_file = None
        else:
            is_mcif = struct_file.lower().endswith('.mcif')
            spin_axis_cart = None
            sf_result = None
            if is_mcif:
                print("Detected .mcif file --magnetic moments will be read from file.")
                moments_str = ""
            else:
                print("Spin axis in Cartesian coordinates (default: 0 0 1): ", end='', flush=True)
                spin_axis_cart = input().strip()
                print("Magnetic moments along this axis (atom order, trailing atoms auto-fill to 0): ", end='', flush=True)
                moments_str = input().strip()
            if not is_mcif and not moments_str:
                standard_path_reason = "No magnetic moments entered."
                standard_path_reason_reported = True
                _step0_wrote_flip_file = False
                print(f"{BOLD}[Note] {standard_path_reason}{RESET} Ordinary structural path will be written.")
            else:
                # Record mtime before so we know if find_sf_run wrote a fresh file
                _mtime_before = os.path.getmtime(_flip_file) if os.path.exists(_flip_file) else None
                sf_result = find_sf_run(
                    struct_file,
                    moments_str,
                    verbose=False,
                    spin_axis_cart=spin_axis_cart,
                )
                _mtime_after = os.path.getmtime(_flip_file) if os.path.exists(_flip_file) else None
                _step0_wrote_flip_file = (
                    _mtime_after is not None and _mtime_after != _mtime_before
                )
            if isinstance(sf_result, dict):
                magnetic_setting_counts = None
                magnetic_setting_outputs = None
                if self.magnetic_setting:
                    try:
                        mag_setting = prepare_magnetic_setting_files(
                            struct_file,
                            moments_str=moments_str,
                            spin_axis_cart=spin_axis_cart,
                            output_dir='.',
                        )
                        centroid_struct_file = mag_setting["helper_path"]
                        centroid_seekpath_type_numbers = mag_setting["seekpath_type_numbers"]
                        magnetic_setting_counts = mag_setting
                    except Exception as e:
                        print(f"[Warning] --ssg-setting failed: {e}")
                        print("[Warning] Falling back to structural SeeK-path setting.")
                        centroid_struct_file = struct_file
                        centroid_seekpath_type_numbers = None
                if CENTROID_AVAILABLE:
                    try:
                        centroid_result = compute_centroid(
                            centroid_struct_file, output_dir='.', show_plot=True,
                            defer_show=True, verbose=False,
                            seekpath_type_numbers=centroid_seekpath_type_numbers,
                        )
                        if self.magnetic_setting and magnetic_setting_counts is not None:
                            magnetic_setting_outputs = finalize_magnetic_setting_outputs(
                                magnetic_setting_counts,
                                centroid_result,
                                output_dir='.',
                                verbose_output=self.output_verbose,
                            )
                            if magnetic_setting_outputs:
                                centroid_result["b_matrix_output"] = magnetic_setting_outputs[
                                    "b_matrix_output"
                                ]
                                if self.output_verbose:
                                    print(
                                        "[SSG setting] Kept intermediates in "
                                        f"{magnetic_setting_outputs.get('intermediate_dir')}"
                                    )
                        display_figures.extend(centroid_result.get('display_figures', []))
                    except Exception as e:
                        centroid_error = e
                if magnetic_setting_counts is not None:
                    _step0_wrote_flip_file = (
                        magnetic_setting_counts.get('spin_flip_operations', 0) > 0
                    )
                else:
                    _step0_wrote_flip_file = sf_result.get('spin_flip_operations', 0) > 0
                laue_no_altermag = None
                if no_altermagnetism_reason is not None:
                    laue_no_altermag = no_altermagnetism_reason(
                        sf_result.get('point_group'),
                    )
                elif sf_result.get('laue_group') in {'-1', '-3', 'm-3'}:
                    laue_no_altermag = {
                        'laue_group': sf_result.get('laue_group'),
                        'reason': 'No altermagnetism',
                    }
                spin_split_diagnostic = sf_result.get('spin_split_diagnostic', '')

                if centroid_result is not None:
                    print(
                        "\nLattice type: "
                        f"{centroid_result.get('sc_type', centroid_result.get('seekpath_bravais', 'unknown'))}"
                    )
                print(f"Structure: {sf_result['structure_file']}, atoms: {sf_result['num_atoms']}")
                print(f"SG {sf_result['space_group']}, "
                      f"PG {sf_result['point_group']}, "
                      f"Laue {sf_result['laue_group']}")
                print(f"Phase: {sf_result['magnetic_phase']}")
                print(f"Oriented SSG: {sf_result['ssg_index']}")
                print(f"SSG Symbol (Chen-Liu): {sf_result['ssg_symbol']}")
                print(f"MSG without SOC: {sf_result['magnetic_space_group_without_soc']}")

                if laue_no_altermag:
                    laue = laue_no_altermag.get('laue_group', sf_result.get('laue_group'))
                    standard_path_reason = f"Laue group {laue}: no altermagnetism."
                    print(f"{BOLD}[Note] {standard_path_reason}{RESET} Default path will be written.")
                    standard_path_reason_reported = True
                else:
                    print("Spin operations: "
                          f"{sf_result['actual_spin_flip_point_operations']} flip, "
                          f"{sf_result['actual_spin_preserve_point_operations']} preserve")
                    if spin_split_diagnostic:
                        standard_path_reason = spin_split_diagnostic
                        print(f"{BOLD}[Note] {standard_path_reason}{RESET} Default path will be written.")
                        standard_path_reason_reported = True

                if STEP0_VERBOSE_SUMMARY:
                    print(f"Magnetic SG: {sf_result['magnetic_space_group']}")
                    print(f"G0: {sf_result['g0_symbol']} ({sf_result['g0_number']}), "
                          f"L0: {sf_result['l0_symbol']} ({sf_result['l0_number']}), "
                          f"EMPG: {sf_result['empg']}")
                    if sf_result.get('findspingroup_warning'):
                        print(f"FindSpinGroup warning: {sf_result['findspingroup_warning']}")
                    print(f"Spin axis: {sf_result['spin_group']}")
                    unique_ops = sf_result.get('unique_point_operations',
                                               sf_result['total_operations'])
                    print(f"Space-group operations: {sf_result['total_operations']} total")
                    print(f"Point operations: {unique_ops} unique")
                    print("Inversion-extended k operations: "
                          f"{sf_result['extended_spin_flip_point_operations']} spin-flip, "
                          f"{sf_result['extended_spin_preserve_point_operations']} spin-preserving "
                          f"({sf_result['extended_spin_flip_operations']} + "
                          f"{sf_result['extended_spin_preserve_operations']} with translations)")
                    print(f"Saved: {', '.join(sf_result['saved_files'])}")

        if centroid_struct_file and CENTROID_AVAILABLE:
            try:
                if centroid_result is None:
                    if centroid_error is not None:
                        raise centroid_error
                    centroid_result = compute_centroid(
                        centroid_struct_file, output_dir='.', show_plot=True,
                        defer_show=True, verbose=False,
                        seekpath_type_numbers=centroid_seekpath_type_numbers,
                    )
                    display_figures.extend(centroid_result.get('display_figures', []))
                    print(
                        "Lattice type: "
                        f"{centroid_result.get('sc_type', centroid_result.get('seekpath_bravais', 'unknown'))}"
                    )
                print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
                sp_path   = centroid_result['sp_path']
                sp_coords = centroid_result['sp_point_coords']
                displayed_path = centroid_result.get(
                    'band_kpath',
                    centroid_result.get('ibz_kpath', sp_path)
                )
                print(f"Path: {self._format_path(displayed_path)}")
                print("Press [Enter] to use this path, or type a filename to load your own: ", end='', flush=True)
                path_choice = input().strip()
                if not path_choice:
                    # Build kpoints_data in the same HPKOT/SeeK-path convention
                    # as Figure 1.  lattice_kpoints.py may include curated
                    # closure vertices for the hull, but ibz_kpath contains only
                    # the public band-path labels.  path_kpoints_frac keeps
                    # optional path-only labels such as H_2 available without
                    # adding them to the centroid hull.
                    self.kpoints_data = []
                    sc_type_auto = centroid_result.get('sc_type', '')
                    if (
                        ('band_kpath' in centroid_result and 'band_kpoints_frac' in centroid_result)
                        or ('ibz_kpath' in centroid_result and 'ibz_kpoints_frac' in centroid_result)
                    ):
                        self.header_lines = [f'K-Path generated by AlterSeeK-Path (HPKOT {sc_type_auto})',
                                             '20', 'Line-Mode', 'Reciprocal']
                        self.kpoints_basis_matrix = np.array(centroid_result['b_matrix'], dtype=float)
                        self.kpoints_basis_rotation = np.array(
                            centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                            dtype=float,
                        )
                        self.output_basis_matrix = np.array(
                            centroid_result.get(
                                'b_matrix_output',
                                centroid_result.get('b_matrix_input', centroid_result['b_matrix']),
                            ),
                            dtype=float,
                        )
                        # Prefer the selected band path when present. This
                        # keeps the prompt, Figure 1 path overlay, and KPOINTS
                        # path consistent.
                        auto_path = centroid_result.get(
                            'band_kpath',
                            centroid_result['ibz_kpath']
                        )
                        ibz_coords = centroid_result.get(
                            'band_kpoints_frac',
                            centroid_result.get(
                                'path_kpoints_frac',
                                centroid_result['ibz_kpoints_frac']
                            )
                        )
                        for seg_start, seg_end in auto_path:
                            for label in (seg_start, seg_end):
                                coords = ibz_coords[label]
                                self.kpoints_data.append([coords[0], coords[1], coords[2], label])
                        extra_vertices = centroid_result.get('extra_general_vertices', [])
                        self.extra_general_points = []
                        for label in extra_vertices:
                            if label in ibz_coords:
                                coords = ibz_coords[label]
                                self.extra_general_points.append([coords[0], coords[1], coords[2], label])
                        print(f"Using HPKOT {sc_type_auto} path ({len(auto_path)} segments, {len(self.kpoints_data)} k-points)")
                        if self.extra_general_points:
                            labels = ", ".join(str(pt[3]) for pt in self.extra_general_points)
                            print(f"Extra doubled-IBZ general anchors: {labels}")
                    else:
                        self.header_lines = ['K-Path generated by AlterSeeK-Path (seekpath)', '20', 'Line-Mode', 'Reciprocal']
                        self.kpoints_basis_matrix = np.array(centroid_result['b_matrix'], dtype=float)
                        self.kpoints_basis_rotation = np.array(
                            centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                            dtype=float,
                        )
                        self.output_basis_matrix = np.array(
                            centroid_result.get(
                                'b_matrix_output',
                                centroid_result.get('b_matrix_input', centroid_result['b_matrix']),
                            ),
                            dtype=float,
                        )
                        for seg_start, seg_end in sp_path:
                            for label in (seg_start, seg_end):
                                coords = sp_coords[label]
                                self.kpoints_data.append([coords[0], coords[1], coords[2], label])
                        print(f"Using auto-generated path ({len(sp_path)} segments, {len(self.kpoints_data)} k-points)")
                else:
                    if not self.read_kpoints_file(path_choice): return
            except Exception as e:
                print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
                print(f"[Warning] Auto path generation failed: {e}")
                print("Falling back to manual file input.")
                print("Enter KPOINTS file name (default: KPATH.in): ", end='', flush=True)
                filename = input().strip()
                if not filename: filename = "KPATH.in"
                if not self.read_kpoints_file(filename): return
        else:
            print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
            print("Enter KPOINTS file name (default: KPATH.in): ", end='', flush=True)
            filename = input().strip()
            if not filename: filename = "KPATH.in"
            if not self.read_kpoints_file(filename): return

        # Laue groups -1, -3, and m-3 do not have a one-dimensional,
        # nonidentical inversion-even irrep, so no altermagnetic splitting is
        # possible. Write the ordinary IBZ path without butterfly insertion.
        no_altermag = None
        if centroid_result is not None:
            no_altermag = centroid_result.get('no_altermagnetism')
            if no_altermag is None and no_altermagnetism_reason is not None:
                no_altermag = no_altermagnetism_reason(
                    centroid_result.get('point_group'),
                    centroid_result.get('spacegroup'),
                )
        if standard_path_reason is None and no_altermag:
            laue = no_altermag.get('laue_group', 'unknown')
            standard_path_reason = f"Laue group {laue}: no altermagnetism."

        # Step 2: Auto-compute or manually enter general k-point
        print(f"\n{BOLD}>>> Step 2: General k-point{RESET}")
        general_kpoint = None

        if centroid_result is not None:
            try:
                c = centroid_result['centroid_frac']
                general_kpoint = [c[0], c[1], c[2]]
                print(f"IBZ centroid: [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
            except Exception as e:
                print(f"[Warning] Centroid retrieval failed: {e}")
        elif struct_file and CENTROID_AVAILABLE:
            try:
                result = compute_centroid(centroid_struct_file, output_dir='.', show_plot=True,
                                          defer_show=True, verbose=False,
                                          seekpath_type_numbers=centroid_seekpath_type_numbers)
                display_figures.extend(result.get('display_figures', []))
                c = result['centroid_frac']
                general_kpoint = [c[0], c[1], c[2]]
                print(f"IBZ centroid: [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
            except Exception as e:
                print(f"[Warning] Centroid computation failed: {e}")

        # Append centroid to spin_operations.txt for reference
        if general_kpoint is not None:
            try:
                with open("spin_operations.txt", "a") as f:
                    f.write(f"\nGeneral k-point (IBZ centroid, fractional): "
                            f"[{general_kpoint[0]:.6f}, {general_kpoint[1]:.6f}, {general_kpoint[2]:.6f}]\n")
            except Exception:
                pass

        if general_kpoint is None:
            print("Format: kx ky kz (space-separated)")
            while True:
                try:
                    print("Enter k-point: ", end='', flush=True)
                    k_input = input().strip().split()
                    if len(k_input) == 3:
                        general_kpoint = [float(x) for x in k_input]
                        break
                    else:
                        print("Please enter exactly 3 coordinates.")
                except ValueError:
                    print("Invalid input. Please enter three numbers.")

        if standard_path_reason:
            if not standard_path_reason_reported:
                print(f"\n{BOLD}[Note] {standard_path_reason}{RESET}")
            print("Writing ordinary k-path with non-spin-flip general-k anchors.")
            print(f"\n{BOLD}>>> Step 5: Save{RESET}")
            print("Enter output filename (default: KPOINTS_modified): ", end='', flush=True)
            output_file = input().strip()
            if not output_file:
                output_file = "KPOINTS_modified"
            standard_general_path = self.insert_general_kpoint_anchors(
                general_kpoint,
                self.extra_general_points,
            )
            if self.write_kpoints_file(standard_general_path, output_file, None):
                if centroid_result is not None:
                    write_bandplot_lattice_config(
                        centroid_result.get('lattice_key', centroid_result.get('sc_type'))
                    )
            print("\nDone.")
            if display_figures and plt is not None:
                print('Displaying generated figure(s)...')
                plt.show()
                for fig in display_figures:
                    save_after_show = getattr(fig, '_alterseek_save_after_show', None)
                    if save_after_show is not None:
                        save_after_show()
            return

        # Step 3: Input transformation matrix
        print(f"\n{BOLD}>>> Step 3: Spin-flip operation{RESET}")
        if _step0_wrote_flip_file is False:
            # Step 0 ran but found no flip ops for this structure.
            # Don't load a stale file from a previous run on a different structure.
            print("[Note] Step 0 found no spin-flip operations --skipping any existing spin_flip_operations.txt.")
            flip_ops = []
        else:
            flip_ops = self.load_flip_operations()
        preserve_ops = self.load_preserve_operations()

        # Always include inversion-extended spatial partners. The spin-flip
        # classification comes from the spin rotation in find_sf_operations.py;
        # multiplying the spatial operation by inversion does not change that
        # spin-flip status. Deduplicate after extension.
        def _inversion_extended(ops):
            expanded = list(ops)
            for op in ops:
                neg_op = -np.array(op, dtype=float)
                if not any(np.allclose(neg_op, ex, atol=1e-8) for ex in expanded):
                    expanded.append(neg_op)
            return expanded

        flip_ops = _inversion_extended(flip_ops)
        preserve_ops = _inversion_extended(preserve_ops)

        R = None
        selected_transformation_label = None
        if flip_ops:
            print(f"Found {len(flip_ops)} spin-flip operations.")
            print("Default R: Option 1")
            while R is None:
                print("Press [Enter] to use default, type a number, 'list' to show matrices, or 'manual': ", end='', flush=True)
                choice = input().strip().lower()

                if not choice:
                    R = flip_ops[0]
                    selected_transformation_label = "Option 1"
                    print("Selected: Option 1")
                elif choice == 'list':
                    for i, op in enumerate(flip_ops):
                        print(f"\n  Option {i+1}:")
                        print(self._format_matrix(op))
                    print()
                elif choice == 'manual':
                    break
                else:
                    try:
                        idx = int(choice) - 1
                        if 0 <= idx < len(flip_ops):
                            R = flip_ops[idx]
                            selected_transformation_label = f"Option {idx+1}"
                            print(f"Selected: Option {idx+1}")
                        else:
                            print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")
                    except ValueError:
                        print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")

        if R is None:
            print("Enter custom transformation matrix R.")
            print("Enter row by row (3 numbers per row, space-separated):")
            transformation_matrix = []
            for i in range(3):
                while True:
                    try:
                        print(f"Row {i+1}: ", end='', flush=True)
                        row_input = input().strip().split()
                        if len(row_input) == 3:
                            transformation_matrix.append([float(x) for x in row_input])
                            break
                        print("Please enter exactly 3 numbers.")
                    except ValueError: pass
            R = np.array(transformation_matrix)
            selected_transformation_label = "manual"
        
        # Step 4: Process k-points
        print(f"\n{BOLD}>>> Step 4: Build altermagnetic path{RESET}")

        # find_sf_operations.py writes FindSpinGroup rotations in the input file's
        # fractional basis.  KPOINTS/IBZ coordinates are in the seekpath primitive
        # reciprocal basis, which can differ for centered lattices (e.g. BCT, RHL).
        # Convert through Cartesian k-space so k', Figure 2, and the path all use
        # the same physical spin-flip operation:
        #   R_cart_k       = b_input.T @ inv(R_input).T @ inv(b_input.T)
        #   R_prim^{-T}    = inv(b_prim.T) @ R_cart_k @ b_prim.T
        R_cart_for_plot = None
        flip_ops_for_plot = flip_ops
        preserve_ops_for_plot = preserve_ops
        if centroid_result is not None and 'b_matrix' in centroid_result:
            _b_input = np.array(centroid_result.get('b_matrix_input',
                                                    centroid_result['b_matrix_conv']), dtype=float)
            _b_prim = np.array(centroid_result['b_matrix'], dtype=float) @ np.array(
                centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                dtype=float,
            )
            def _convert_input_frac_R_to_prim(_R):
                _R_arr = np.array(_R, dtype=float)
                _R_cart = _b_input.T @ np.linalg.inv(_R_arr).T @ np.linalg.inv(_b_input.T)
                _R_prim_inv_T = np.linalg.inv(_b_prim.T) @ _R_cart @ _b_prim.T
                return np.linalg.inv(_R_prim_inv_T.T), _R_cart

            R_for_kpts, _ = _convert_input_frac_R_to_prim(R)
            # Figure 2 draws the HPKOT hull in seekpath's standardized
            # Cartesian frame.  Let it reconstruct R_cart from R_for_kpts in
            # that same frame; the Cartesian matrix above remains in the
            # orientation of the input structure (notably different for MCIF).
            R_cart_for_plot = None
            if flip_ops:
                flip_ops_for_plot = [
                    _convert_input_frac_R_to_prim(op)[0] for op in flip_ops
                ]
            if preserve_ops:
                preserve_ops_for_plot = [
                    _convert_input_frac_R_to_prim(op)[0] for op in preserve_ops
                ]
            def _annotate_ops_with_standardized_basis(filename, input_ops, standardized_ops, label):
                try:
                    with open(filename, "w") as f:
                        f.write(f"# Found {len(input_ops)} inversion-extended {label} point operations\n")
                        f.write("# Left matrix: input POSCAR fractional basis.\n")
                        f.write("# Right matrix: standardized SeeK-path/HPKOT primitive basis used by Figures 1-4 and KPOINTS.\n")
                        for i, (input_op, std_op) in enumerate(zip(input_ops, standardized_ops), 1):
                            f.write(f"Operation_{i}\n")
                            input_int = np.rint(np.array(input_op, dtype=float)).astype(int)
                            std_int = np.rint(np.array(std_op, dtype=float)).astype(int)
                            for input_row, std_row in zip(input_int, std_int):
                                left = " ".join(str(int(x)) for x in input_row)
                                right = " ".join(str(int(x)) for x in std_row)
                                f.write(f"{left}    |    {right}\n")
                            f.write("\n")
                except Exception as exc:
                    print(f"[Warning] Could not write {filename}: {exc}")

            _annotate_ops_with_standardized_basis(
                "spin_flip_operations.txt",
                flip_ops,
                flip_ops_for_plot,
                "spin-flipping",
            )
            _annotate_ops_with_standardized_basis(
                "spin_preserve_operations.txt",
                preserve_ops,
                preserve_ops_for_plot,
                "spin-preserving",
            )
            print("[Basis] Converted R from input-cell basis to primitive basis.")
            print("[Basis] Annotated spin operation files with standardized-basis matrices.")
            print("Primitive-basis R used for KPOINTS:")
            print(self._format_matrix(R_for_kpts))
        else:
            R_for_kpts = R

        # Calculate and show k'
        k_prime = self.transform_kpoint(general_kpoint, R_for_kpts)
        print(f"k' = [{k_prime[0]:.4f}, {k_prime[1]:.4f}, {k_prime[2]:.4f}]")

        try:
            new_kpoints = self.insert_general_kpoints(
                general_kpoint, R_for_kpts, self.extra_general_points
            )

            if new_kpoints:
                # Generate Figures 2-4 (spin-flip, spin-BZ, kz=0 top view).
                # One shared call scaffold; per-figure kwargs hold the
                # differences between the three plots.
                if centroid_result is not None:
                    basename = (os.path.splitext(os.path.basename(struct_file))[0]
                                if struct_file else 'output')
                    sc_type = centroid_result.get('sc_type', 'BZ')
                    flip_kwargs = dict(
                        flip_ops_frac=flip_ops_for_plot if flip_ops_for_plot else None,
                        preserve_ops_frac=preserve_ops_for_plot if preserve_ops_for_plot else None,
                    )
                    view_kwargs = dict(
                        bz_center=centroid_result.get('bz_center'),
                        bz_span=centroid_result.get('bz_span'),
                        elev=centroid_result.get('elev', 14),
                        azim=centroid_result.get('azim', 20),
                    )
                    figure_specs = []
                    if plot_spin_flip_figure is not None and 'b_matrix' in centroid_result:
                        figure_specs.append((
                            plot_spin_flip_figure, 'spin-flip',
                            f'{basename}_spinflip_{sc_type}.png',
                            dict(
                                kpoints_data=self.kpoints_data,
                                ibz_kpoints_frac=centroid_result.get('ibz_kpoints_frac', {}),
                                centroid_frac=general_kpoint,
                                R_cart=R_cart_for_plot,
                                block=False,
                                path_sequence=new_kpoints,
                                unique_ops=centroid_result.get('unique_ops'),
                                **view_kwargs,
                            ),
                        ))
                    if 'unique_ops' in centroid_result:
                        spin_bz_kwargs = dict(
                            unique_ops=centroid_result['unique_ops'],
                            centroid_cart=centroid_result.get('centroid_cart'),
                            hull_labels=centroid_result.get('hull_labels'),
                            **flip_kwargs,
                        )
                        if plot_spin_bz_figure is not None:
                            figure_specs.append((
                                plot_spin_bz_figure, 'spin-BZ',
                                f'{basename}_spinbz_{sc_type}.png',
                                dict(**spin_bz_kwargs, **view_kwargs),
                            ))
                        if plot_spin_bz_top_view_figure is not None:
                            figure_specs.append((
                                plot_spin_bz_top_view_figure, 'spin-BZ top-view',
                                f'{basename}_spinbz_top_{sc_type}.png',
                                dict(z0=0.0, **spin_bz_kwargs),
                            ))
                    for plot_fn, fig_name, fig_path, extra_kwargs in figure_specs:
                        try:
                            fig = plot_fn(
                                b_matrix=centroid_result['b_matrix'],
                                bz_loops=centroid_result['bz_loops'],
                                hull_pts=centroid_result.get('hull_pts'),
                                hull_simplices=centroid_result.get('hull_simplices'),
                                R=R_for_kpts,
                                output_path=fig_path,
                                show_plot=True,
                                defer_show=True,
                                **extra_kwargs,
                            )
                            if fig is not None:
                                display_figures.append(fig)
                        except Exception as _e:
                            print(f"[Warning] Could not generate {fig_name} figure: {_e}")
                # Step 5: Save modified file
                print(f"\n{BOLD}>>> Step 5: Save{RESET}")
                print("Enter output filename (default: KPOINTS_modified): ", end='', flush=True)
                output_file = input().strip()
                if not output_file:
                    output_file = "KPOINTS_modified"
                
                if self.write_kpoints_file(
                    new_kpoints, output_file, R, selected_transformation_label
                ):
                    if centroid_result is not None:
                        write_bandplot_lattice_config(
                            centroid_result.get('lattice_key', centroid_result.get('sc_type'))
                        )
                print("\nDone.")
                if display_figures and plt is not None:
                    print('Displaying generated figures...')
                    plt.show()
                    for fig in display_figures:
                        save_after_show = getattr(fig, '_alterseek_save_after_show', None)
                        if save_after_show is not None:
                            save_after_show()
            else:
                print("Error: Failed to process k-points.")
                
        except Exception as e:
            print(f"Error processing k-points: {e}")

def main():
    import argparse

    argv = sys.argv[1:]
    # Forward the bandplot subcommand untouched so plot_alterband's own
    # parser handles its options (including --help).
    if argv and argv[0].lower() in {"bandplot", "plot-band", "plot"}:
        from plot_alterband import main as plot_alterband_main
        plot_alterband_main(argv[1:])
        return

    parser = argparse.ArgumentParser(
        prog="alterseek-path",
        description=(
            "Generate an altermagnetic KPOINTS path interactively. "
            "Subcommand: 'bandplot' plots spin-resolved bands from KLABELS "
            "and reformatted band data (run: alterseek-path bandplot --help)."
        ),
    )
    parser.add_argument(
        "--ssg-setting",
        action="store_true",
        help="Experimental: generate Figure 1/KPOINTS from FindSpinGroup SSG setting.",
    )
    parser.add_argument(
        "--output",
        choices=["verbose"],
        help="verbose: keep SSG-setting intermediate/helper structures for debugging.",
    )
    args = parser.parse_args(argv)

    modifier = KPointsModifier(
        magnetic_setting=args.ssg_setting,
        output_verbose=args.output == "verbose",
    )
    modifier.interactive_modify()


if __name__ == "__main__":
    main()

