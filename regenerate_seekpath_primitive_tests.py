#!/usr/bin/env python3
"""Regenerate primitive_test inputs in the SeeK-path/HPKOT primitive basis.

This intentionally does not promote primitive files to active POSCAR/INCAR.
It rewrites only primitive_test/ contents for folders that already have a
primitive_test directory.
"""

from __future__ import annotations

import contextlib
import os
from pathlib import Path
import re
import shutil

import numpy as np
import seekpath
import spglib
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element

from find_sf_operations import (
    altermagnetic_diagnostic,
    count_unique_point_operations,
    has_spin_flip_inversion,
    has_spin_flip_translation,
    operation_count_summary,
    write_flip_ops_to_file,
    write_operations_to_file,
    write_preserve_ops_to_file,
)


ROOT = Path(__file__).resolve().parents[1] / "structure" / "all paper structures"


def read_poscar_species_counts(path: Path) -> tuple[list[str], list[int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    species = lines[5].split()
    counts_line_index = 6
    if all(_is_number(tok) for tok in species):
        species = [str(site.specie) for site in Structure.from_file(path).types_of_species]
        counts_line_index = 5
    counts = [int(float(tok)) for tok in lines[counts_line_index].split()]
    return species, counts


def _is_number(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def expand_magmom(text: str) -> list[float]:
    values: list[float] = []
    for token in text.replace("=", " ").split():
        if token.upper() == "MAGMOM":
            continue
        if "*" in token:
            n, value = token.split("*", 1)
            values.extend([float(value)] * int(float(n)))
        else:
            values.append(float(token))
    return values


def compress_magmom(values: list[float]) -> str:
    def fmt(value: float) -> str:
        if abs(value - round(value)) < 1e-8:
            return str(int(round(value)))
        return f"{value:.6g}"

    out: list[str] = []
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and abs(values[j] - values[i]) < 1e-8:
            j += 1
        n = j - i
        out.append(f"{n}*{fmt(values[i])}" if n > 1 else fmt(values[i]))
        i = j
    return " ".join(out)


def read_incar_magmom(path: Path) -> list[float]:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if re.match(r"^\s*MAGMOM\s*=", line, flags=re.IGNORECASE):
            return expand_magmom(line)
    raise ValueError(f"No MAGMOM line in {path}")


def update_incar_magmom(source: Path, target: Path, magmoms: list[float]) -> None:
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    replacement = f" MAGMOM =  {compress_magmom(magmoms)}"
    updated = []
    done = False
    for line in lines:
        if re.match(r"^\s*MAGMOM\s*=", line, flags=re.IGNORECASE):
            updated.append(replacement)
            done = True
        else:
            updated.append(line)
    if not done:
        updated.append(replacement)
    target.write_text("\n".join(updated) + "\n", encoding="utf-8")


def write_poscar(path: Path, title: str, lattice: np.ndarray, species_order: list[str],
                 rows: list[tuple[str, np.ndarray]]) -> None:
    grouped: list[tuple[str, np.ndarray]] = []
    counts: list[int] = []
    for element in species_order:
        items = [(sym, frac) for sym, frac in rows if sym == element]
        if items:
            counts.append(len(items))
            grouped.extend(items)

    with path.open("w", encoding="utf-8") as f:
        f.write(title + "\n")
        f.write("1.0\n")
        for vec in lattice:
            f.write("  " + " ".join(f"{x:20.16f}" for x in vec) + "\n")
        f.write(" ".join(species_order) + "\n")
        f.write(" ".join(str(x) for x in counts) + "\n")
        f.write("Direct\n")
        for sym, frac in grouped:
            frac = np.mod(frac, 1.0)
            frac[np.isclose(frac, 1.0, atol=1e-10)] = 0.0
            frac[np.isclose(frac, 0.0, atol=1e-10)] = 0.0
            f.write("  " + " ".join(f"{x:20.16f}" for x in frac) + f" {sym}\n")


def map_magmoms_to_primitive(source: Structure, conv_lattice: np.ndarray,
                             primitive_lattice: np.ndarray,
                             primitive_positions: np.ndarray, primitive_types: np.ndarray,
                             rotation_matrix: np.ndarray,
                             source_magmoms: list[float]) -> tuple[list[tuple[str, np.ndarray]], list[float]]:
    source_lattice = np.array(source.lattice.matrix, dtype=float)
    source_positions = np.array([site.frac_coords for site in source], dtype=float)
    source_numbers = np.array([site.specie.Z for site in source], dtype=int)
    source_magmoms_arr = np.array(source_magmoms, dtype=float)
    rotation_matrix = np.array(rotation_matrix, dtype=float)
    conv_lattice = np.array(conv_lattice, dtype=float)

    # SeeK-path may rotate the input cell into its standardized conventional
    # Cartesian frame before applying the primitive transformation.  Map source
    # atoms into that same conventional fractional frame before matching.
    source_conv_positions = []
    source_primitive_positions = []
    for frac in source_positions:
        source_cart = source_lattice.T @ frac
        conv_cart = rotation_matrix @ source_cart
        source_conv_positions.append(np.linalg.solve(conv_lattice.T, conv_cart) % 1.0)
        source_primitive_positions.append(np.linalg.solve(primitive_lattice.T, conv_cart) % 1.0)
    source_conv_positions = np.array(source_conv_positions)
    source_primitive_positions = np.array(source_primitive_positions)

    def best_primitive_origin_shift() -> tuple[float, np.ndarray]:
        candidates: list[np.ndarray] = []
        first_z = int(primitive_types[0])
        first_frac = np.array(primitive_positions[0], dtype=float) % 1.0
        for frac, z0 in zip(source_primitive_positions, source_numbers):
            if int(z0) == first_z:
                candidates.append((first_frac - frac) % 1.0)

        best: tuple[float, np.ndarray] | None = None
        for shift in candidates:
            worst = 0.0
            for frac_prim, z in zip(primitive_positions, primitive_types):
                frac_prim = np.array(frac_prim, dtype=float) % 1.0
                nearest = None
                for source_frac, z0 in zip(source_primitive_positions, source_numbers):
                    if int(z0) != int(z):
                        continue
                    delta = (frac_prim - (source_frac + shift) + 0.5) % 1.0 - 0.5
                    distance = float(np.linalg.norm(primitive_lattice.T @ delta))
                    if nearest is None or distance < nearest:
                        nearest = distance
                worst = max(worst, nearest if nearest is not None else 1e9)
            if best is None or worst < best[0]:
                best = (worst, shift)
        if best is None:
            return 1e9, np.zeros(3)
        return best

    primitive_shift_error, primitive_shift = best_primitive_origin_shift()

    rows: list[tuple[str, np.ndarray]] = []
    magmoms: list[float] = []
    for frac_prim, z in zip(primitive_positions, primitive_types):
        cart = primitive_lattice.T @ np.array(frac_prim, dtype=float)
        frac_source = np.linalg.solve(conv_lattice.T, cart) % 1.0
        best: tuple[float, int] | None = None
        for idx, (frac, z0) in enumerate(zip(source_conv_positions, source_numbers)):
            if int(z0) != int(z):
                continue
            delta = (frac_source - frac + 0.5) % 1.0 - 0.5
            distance = float(np.linalg.norm(conv_lattice.T @ delta))
            if best is None or distance < best[0]:
                best = (distance, idx)
        if best is None or best[0] > 1e-4:
            frac_prim_mod = np.array(frac_prim, dtype=float) % 1.0
            best = None
            for idx, (frac, z0) in enumerate(zip(source_primitive_positions, source_numbers)):
                if int(z0) != int(z):
                    continue
                delta = (frac_prim_mod - (frac + primitive_shift) + 0.5) % 1.0 - 0.5
                distance = float(np.linalg.norm(primitive_lattice.T @ delta))
                if best is None or distance < best[0]:
                    best = (distance, idx)
        if best is None or best[0] > 1e-4 or primitive_shift_error > 1e-4:
            raise ValueError(
                f"Could not map primitive atom Z={z} at {frac_prim}; "
                f"best atom distance={best[0] if best else None}, "
                f"origin-shift error={primitive_shift_error}"
            )
        rows.append((Element.from_Z(int(z)).symbol, np.array(frac_prim, dtype=float)))
        magmoms.append(float(source_magmoms_arr[best[1]]))
    return rows, magmoms


def reciprocal_basis(lattice: np.ndarray) -> np.ndarray:
    return 2 * np.pi * np.linalg.inv(lattice).T


def basis_status(lattice: np.ndarray, positions: np.ndarray, numbers: np.ndarray) -> tuple[bool, float, int]:
    sp_result = seekpath.get_path((lattice, positions, numbers), with_time_reversal=True)
    b_input = reciprocal_basis(lattice)
    b_seek = np.array(sp_result["reciprocal_primitive_lattice"], dtype=float)
    transform = np.linalg.inv(b_seek.T) @ b_input.T
    nearest = np.rint(transform)
    err = float(np.max(np.abs(transform - nearest)))
    det = int(round(np.linalg.det(nearest)))
    return err < 1e-5 and abs(det) == 1, err, det


def validate_spin(primitive_dir: Path, magmoms: list[float]) -> dict[str, object]:
    structure = Structure.from_file(primitive_dir / "POSCAR")
    lattice = np.array(structure.lattice.matrix, dtype=float)
    positions = np.array([site.frac_coords for site in structure], dtype=float)
    numbers = np.array([site.specie.Z for site in structure], dtype=int)
    magmom_vectors = np.array([[0.0, 0.0, m] for m in magmoms], dtype=float)

    dataset = spglib.get_symmetry_dataset((lattice, positions, numbers))
    sog, rotations, translations, spin_rotations = __import__("spinspg").get_spin_symmetry(
        lattice, positions, numbers, magmom_vectors, symprec=1e-5
    )
    counts = operation_count_summary(rotations, spin_rotations)
    diagnostic = altermagnetic_diagnostic(rotations, translations, spin_rotations)

    label_info = (
        f"Non-Magnetic Label: {dataset.international} ({dataset.number})\n"
        f"Spin-Only Group Type: {sog}\n"
        "Magnetic Space Group Label: not queried in primitive regeneration"
    )
    cwd = Path.cwd()
    try:
        os.chdir(primitive_dir)
        write_operations_to_file("spin_operations.txt", rotations, translations, spin_rotations, label_info, verbose=False)
        write_flip_ops_to_file("spin_flip_operations.txt", rotations, spin_rotations, verbose=False)
        write_preserve_ops_to_file("spin_preserve_operations.txt", rotations, spin_rotations, verbose=False)
    finally:
        os.chdir(cwd)

    basis_ok, basis_err, basis_det = basis_status(lattice, positions, numbers)
    return {
        "spacegroup": f"{dataset.international} ({dataset.number})",
        "pointgroup": dataset.pointgroup,
        "spin_group": str(sog),
        "unique_point_operations": count_unique_point_operations(rotations),
        "spin_flip_point_operations": counts["actual_spin_flip_point_operations"],
        "spin_preserve_point_operations": counts["actual_spin_preserve_point_operations"],
        "pt": has_spin_flip_inversion(rotations, spin_rotations),
        "ut": has_spin_flip_translation(rotations, translations, spin_rotations),
        "diagnostic": diagnostic or "AM (no PT/Ut warning)",
        "basis_ok": basis_ok,
        "basis_err": basis_err,
        "basis_det": basis_det,
    }


def regenerate_case(primitive_dir: Path) -> str:
    case_dir = primitive_dir.parent
    if "7-cubic" in case_dir.parts:
        return f"SKIP cubic {case_dir}"

    source_poscar = case_dir / "POSCAR_conventional"
    source_incar = case_dir / "INCAR_conventional"
    if not source_poscar.exists():
        source_poscar = case_dir / "POSCAR"
    if not source_incar.exists():
        source_incar = case_dir / "INCAR"

    source = Structure.from_file(source_poscar)
    species_order, source_counts = read_poscar_species_counts(source_poscar)
    source_magmoms = read_incar_magmom(source_incar)
    if len(source_magmoms) != len(source):
        raise ValueError(f"{case_dir}: MAGMOM count {len(source_magmoms)} != atom count {len(source)}")

    lattice = np.array(source.lattice.matrix, dtype=float)
    positions = np.array([site.frac_coords for site in source], dtype=float)
    numbers = np.array([site.specie.Z for site in source], dtype=int)
    sp_result = seekpath.get_path((lattice, positions, numbers), with_time_reversal=True)
    primitive_lattice = np.array(sp_result["primitive_lattice"], dtype=float)
    primitive_positions = np.array(sp_result["primitive_positions"], dtype=float) % 1.0
    primitive_types = np.array(sp_result["primitive_types"], dtype=int)

    rows, primitive_magmoms_raw = map_magmoms_to_primitive(
        source,
        np.array(sp_result["conv_lattice"], dtype=float),
        primitive_lattice,
        primitive_positions,
        primitive_types,
        np.array(sp_result["rotation_matrix"], dtype=float),
        source_magmoms,
    )

    grouped_rows: list[tuple[str, np.ndarray]] = []
    grouped_magmoms: list[float] = []
    for element in species_order:
        for (sym, frac), mag in zip(rows, primitive_magmoms_raw):
            if sym == element:
                grouped_rows.append((sym, frac))
                grouped_magmoms.append(mag)

    primitive_dir.mkdir(exist_ok=True)
    write_poscar(
        primitive_dir / "POSCAR",
        f"{source.composition.reduced_formula}; SeeK-path primitive test",
        primitive_lattice,
        species_order,
        grouped_rows,
    )
    update_incar_magmom(source_incar, primitive_dir / "INCAR", grouped_magmoms)

    validation = validate_spin(primitive_dir, grouped_magmoms)
    old_count = len(source)
    new_count = len(grouped_rows)
    pass_ok = (
        validation["basis_ok"]
        and validation["spin_flip_point_operations"] > 0
        and not validation["pt"]
        and not validation["ut"]
    )
    text = [
        f"Case: {case_dir.name}",
        f"Folder: {case_dir.relative_to(ROOT)}",
        f"Source POSCAR: {source_poscar.name}",
        f"Atom count: {old_count} -> {new_count}",
        f"MAGMOM count: {len(grouped_magmoms)}",
        f"Seekpath Bravais: {sp_result['bravais_lattice_extended']}",
        f"Space group / point group: {validation['spacegroup']} / {validation['pointgroup']}",
        f"Spin group: {validation['spin_group']}",
        f"Actual spin-flip point operations: {validation['spin_flip_point_operations']}",
        f"Actual spin-preserve point operations: {validation['spin_preserve_point_operations']}",
        f"PT spin flip: {validation['pt']}",
        f"Ut spin flip: {validation['ut']}",
        f"Diagnostic warning: {validation['diagnostic']}",
        f"Seekpath basis check: {validation['basis_ok']} (err={validation['basis_err']:.3e}, det={validation['basis_det']})",
        f"PASS: {pass_ok}",
        "",
    ]
    (primitive_dir / "primitive_validation.txt").write_text("\n".join(text), encoding="utf-8")
    return f"{'PASS' if pass_ok else 'FAIL'} {case_dir.relative_to(ROOT)} {old_count}->{new_count}"


def main() -> None:
    primitive_dirs = sorted(ROOT.rglob("primitive_test"))
    results = []
    for primitive_dir in primitive_dirs:
        if not (primitive_dir / "POSCAR").exists():
            continue
        try:
            results.append(regenerate_case(primitive_dir))
        except Exception as exc:
            rel = primitive_dir.parent.relative_to(ROOT)
            results.append(f"ERROR {rel}: {exc}")
    print("\n".join(results))


if __name__ == "__main__":
    main()
