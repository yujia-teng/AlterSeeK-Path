# API

```python
from alterseek import (
    spin_symmetry,
    brillouin_zone,
    general_k_path,
    run_workflow,
    SpinSymmetryError,
)
```

*Added in version 1.1.0.*

## `spin_symmetry`

```python
spin_symmetry(structure_file, moments, spin_axis="0 0 1", symprec=None,
              output_dir="alterseek_output", verbose=False)
```

Return the spin-symmetry analysis for a structure.

### Parameters

| Name | Description |
|------|-------------|
| `structure_file` | A `POSCAR` or a `.vasp`, `.cif`, or `.mcif` structure file. |
| `moments` | Collinear moments in atom order, for example `"5 -5 4*0"`. Ignored for MCIF input. |
| `spin_axis` | Cartesian spin axis, for example `"0 0 1"`. Ignored for MCIF input. |
| `symprec` | Symmetry tolerance passed to spglib. `None` uses the default. |
| `output_dir` | Directory for the operation files. |
| `verbose` | Print the analysis summary. |

### Returns

A dictionary containing, among other entries:

| Key | Value |
|-----|-------|
| `magnetic_phase` | Magnetic phase, such as `AFM(Altermagnet)` or `FM`. |
| `spin_split_diagnostic` | Empty when spin splitting is allowed; otherwise the reason it is forbidden. |
| `magnetic_space_group` | Magnetic space group with spin-orbit coupling. |
| `magnetic_space_group_without_soc` | Magnetic space group without spin-orbit coupling. |
| `actual_spin_flip_point_operations` | Number of unique detected spin-flip point operations before inversion extension. |
| `extended_spin_flip_point_operations` | Number of spin-flip point operations after adding inversion partners for k mapping. |
| `spin_group`, `ssg_index`, `ssg_symbol` | Spin-space-group information. |
| `g0_symbol`, `g0_number` | Spatial subgroup information. |
| `saved_files` | Paths of the three operation files. The spin-flip file is listed even when no such file was written. |

### Raises

| Exception | Condition |
|-----------|-----------|
| `SpinSymmetryError` | The structure or magnetic input cannot be read or analysed. |

### Example

```python
sf = spin_symmetry(
    "example/MnF2/VASP/POSCAR",
    "5 -5 4*0",
    spin_axis="0 0 1",
)
sf["magnetic_phase"]
# 'AFM(Altermagnet)'

sf["magnetic_space_group_without_soc"]
# "P4_2'/mn'm (BNS 136.498), Type III"
```

## `brillouin_zone`

```python
brillouin_zone(structure_file, moments=None, spin_axis="0 0 1", mode_2d=False,
               vacuum_axis="c", symprec=None, output_dir="alterseek_output",
               show_plot=True, view_elev=None, view_azim=None, save_pdf=False,
               verbose=False)
```

Return the Brillouin-zone, high-symmetry-path, and IBZ-centroid data of the submitted cell, with the magnetic structure taken into account.

### Parameters

| Name | Description |
|------|-------------|
| `structure_file` | A `POSCAR` or a `.vasp`, `.cif`, or `.mcif` structure file. |
| `moments` | Collinear moments in atom order, for example `"5 -5 4*0"`. `None` or empty treats the structure as non-magnetic. Ignored for MCIF input. |
| `spin_axis` | Cartesian spin axis. Ignored for MCIF input. |
| `mode_2d` | Use the two-dimensional Brillouin-zone implementation. |
| `vacuum_axis` | Vacuum axis of the submitted cell: `"a"`, `"b"`, or `"c"`. Affects the calculation only when `mode_2d=True`. |
| `symprec` | Symmetry tolerance passed to spglib. `None` uses the default. |
| `output_dir` | Output directory. |
| `show_plot` | Show the 3D Brillouin-zone figure. The figure file is written regardless. |
| `view_elev`, `view_azim` | 3D camera elevation and azimuth. |
| `save_pdf` | Also save the 3D figure as PDF. |
| `verbose` | Print the analysis summary. |

### Returns

A dictionary containing, among other entries:

| Key | Value |
|-----|-------|
| `sc_type` | Extended lattice label, for example `tP1`. |
| `centroid_frac`, `centroid_cart` | IBZ centroid in fractional and Cartesian reciprocal coordinates. |
| `sp_path`, `sp_point_coords` | SeeK-path high-symmetry path and point coordinates. |
| `spacegroup`, `sg_symbol`, `point_group`, `laue_group` | Symmetry information. |
| `b_matrix_input` | Reciprocal vectors of the submitted cell. |
| `ibz_volume`, `n_symmetry_ops` | Irreducible-zone volume and operation count. |

### Files

Output files are in `output_dir` (default `alterseek_output/`):

| input | files |
|---|---|
| `POSCAR`, `.vasp` or `.cif` with moments, or `.mcif` | `spin_operations.txt`, `spin_preserve_operations.txt`, `spin_flip_operations.txt` (only when a spin-flip operation exists) and `*_magnetic_primitive.mcif` |
| 3D | also `*_ibz_*.png` and `*_seekpath_basis_mapping.txt` |

### Example

```python
bz = brillouin_zone(
    "example/MnF2/VASP/POSCAR",
    moments="5 -5 4*0",
    show_plot=False,
)
bz["sc_type"]
# 'tP1'

bz["centroid_frac"]
# array([0.16666667, 0.33333333, 0.25      ])
```

## `general_k_path`

```python
general_k_path(structure_file, moments=None, spin_axis="0 0 1", flip_option=1,
               mode_2d=False, vacuum_axis="c", symprec=None,
               output_dir="alterseek_output", verbose=False)
```

Return the general-k path in the submitted cell's reciprocal basis.

### Parameters

| Name | Description |
|------|-------------|
| `structure_file` | A `POSCAR` or a `.vasp`, `.cif`, or `.mcif` structure file. |
| `moments` | Collinear moments in atom order. `None` returns the ordinary path. Ignored for MCIF input. |
| `spin_axis` | Cartesian spin axis. Ignored for MCIF input. |
| `flip_option` | One-based index of the spin-flip operation. Must be a positive integer. |
| `mode_2d` | Use the two-dimensional path implementation. |
| `vacuum_axis` | Vacuum axis of the submitted cell: `"a"`, `"b"`, or `"c"`. Affects the calculation only when `mode_2d=True`. |
| `symprec` | Symmetry tolerance passed to spglib. `None` uses the default. |
| `output_dir` | Output directory. |
| `verbose` | Print the analysis summaries. |

### Returns

A dictionary containing:

| Key | Value |
|-----|-------|
| `segments` | Path segments with `start`, `start_label`, `end`, `end_label`, and `break_before`. |
| `path` | Formatted path labels. |
| `k`, `k_prime` | General k point and its spin-flip image. |
| `spin_flip_operation` | Selected operation, or `None`. |
| `flip_option_count` | Number of selectable inversion-extended spin-flip operations. |
| `no_splitting_reason` | Reason an ordinary path was returned, or `None`. |
| `lattice` | Extended lattice label. |
| `magnetic_phase` | Magnetic phase, or `None` when moments were omitted. |
| `magnetic_space_group_without_soc` | Magnetic space group without spin-orbit coupling. |
| `spin_symmetry` | Result from `spin_symmetry`, or `None` when moments were omitted. |
| `brillouin_zone` | Result from `brillouin_zone` for the same structure. |

### Raises

| Exception | Condition |
|-----------|-----------|
| `ValueError` | `flip_option` or `vacuum_axis` is invalid, no irreducible-wedge centre or path can be built, or no current-run spin-flip operation is available when spin splitting is allowed. |
| `SpinSymmetryError` | The structure is missing, or spin-symmetry analysis cannot read the structure or analyse its moments. |
| `RuntimeError` | Submitted-cell or Brillouin-zone analysis fails. This can include unreadable non-MCIF input when `moments=None`. |

### Files

| input | files |
|---|---|
| `POSCAR`, `.vasp` or `.cif` with moments, or `.mcif` | `spin_operations.txt`, `spin_preserve_operations.txt`, `spin_flip_operations.txt` (only when a spin-flip operation exists) and `*_magnetic_primitive.mcif` |
| 3D | also `*_ibz_*.png` and `*_seekpath_basis_mapping.txt` |

If the call fails halfway, files already written are kept. Operation files left in `output_dir` by an earlier run are ignored.

### Example

```python
kp = general_k_path(
    "example/MnF2/VASP/POSCAR",
    moments="5 -5 4*0",
    flip_option=2,
)
kp["path"]
# "Γ-X-k | k'-X'-M'-k' | k-M-Γ-k | k'-Γ-Z'-k' | k-Z-R-k | k'-R'-A'-k' | k-A-Z | X-R | M-A"

[round(value, 4) for value in kp["k"]]
# [0.1667, 0.3333, 0.25]

[round(value, 4) for value in kp["k_prime"]]
# [0.3333, -0.1667, 0.25]
```

MCIF input: the moments and spin axis are read from the file.

```python
kp = general_k_path("example/python/structures/La2NiO4.mcif")
kp["lattice"]
# 'oP1'

kp["path"]
# "Γ-X-k | k'-X'-S'-k' | k-S-Y-k | k'-Y'-Γ-k' | k-Γ-Z-k | k'-Z'-U'-k' | k-U-R-k | k'-R'-T'-k' | k-T-Z | X-U | Y-T | S-R"
```

2D structure with vacuum along `c`:

```python
kp = general_k_path(
    "example/python/structures/square_4m_POSCAR",
    moments="1 -1 -1 1 6*0",
    mode_2d=True,
    vacuum_axis="c",
)
kp["lattice"]
# 'square'

kp["path"]
# "Γ-X-k | k'-X'-M'-k' | k-M | Γ-k | k'-Γ | X_A-k | k'-X_A'"

[round(value, 4) for value in kp["k"]]
# [0.25, 0.25, 0.0]
```

## `run_workflow`

```python
run_workflow(mode_2d=False, vacuum_axis=None)
```

Run the same interactive Step 0--5 workflow as the `alterseek-path` command,
but start it from Python. Inputs, prompts, and console output are the same.

### Parameters

| Name | Description |
|------|-------------|
| `mode_2d` | Run the two-dimensional workflow, as `alterseek-path --2d` does. |
| `vacuum_axis` | Vacuum axis of the submitted cell: `"a"`, `"b"`, or `"c"`. `None` leaves the choice to `vacuum_axis` in `alterseek_input.toml`. Used only when `mode_2d=True`. |

### Returns

`True` after a successful run and `False` after a failure the workflow has
already reported. Any other exception propagates.

### Raises

| Exception | Condition |
|-----------|-----------|
| `ValueError` | `vacuum_axis` is invalid. |

### Files

The path file `KPOINTS_alter`, `KPOINTS_alter_qe`, or `KPOINTS_alter_abinit` and the matching `alterseek_plot_*.toml` are written to the working directory. The operation files, figures, `*_seekpath_basis_mapping.txt`, `*_magnetic_primitive.mcif`, and `alterseek_run.log` are written to `alterseek_output/`. The run log is saved only after a successful run.

### Examples

```python
from alterseek import run_workflow

run_workflow()
```

```python
run_workflow(mode_2d=True, vacuum_axis="c")
```

See [Workflow](workflow.md) for the inputs, interaction, and generated files.
Use `general_k_path()` instead when the path is needed directly as Python
data.

## `SpinSymmetryError`

```python
SpinSymmetryError(message)
```

Raised when spin-symmetry analysis or a required spin-symmetry output step
fails.

## Example

Save the following as `quickstart.py` and run it from the repository root with
`python quickstart.py`:

```python
from alterseek import general_k_path

result = general_k_path(
    "example/MnF2/VASP/POSCAR",
    moments="5 -5 4*0",
    spin_axis="0 0 1",
    flip_option=2,
    output_dir="alterseek_output",
)

print("phase:", result["magnetic_phase"])
print("#R:", result["flip_option_count"])
print("lattice:", result["lattice"])
print("k:", [round(value, 4) for value in result["k"]])
print("k-prime:", [round(value, 4) for value in result["k_prime"]])
print("path:", result["path"])
```

Output:

```text
phase: AFM(Altermagnet)
#R: 8
lattice: tP1
k: [0.1667, 0.3333, 0.25]
k-prime: [0.3333, -0.1667, 0.25]
path: Γ-X-k | k'-X'-M'-k' | k-M-Γ-k | k'-Γ-Z'-k' | k-Z-R-k | k'-R'-A'-k' | k-A-Z | X-R | M-A
```

The same code can be tried without creating a file by pasting it into an
interactive Python session started with `python`.

## Example script

The example
[`example/python/analyze_structures.py`](https://github.com/yujia-teng/AlterSeeK-Path/blob/main/example/python/analyze_structures.py)
analyzes MnF2, GdAuGe, La2NiO4 from an MCIF, and a 2D square structure, and
prints a summary table of the returned results. From the repository root, run:

```console
python example/python/analyze_structures.py
```
