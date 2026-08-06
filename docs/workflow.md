# Workflow

The main command is interactive:

```bash
alterseek-path
```

It proceeds through structure reading, spin-symmetry analysis, path selection,
general k-point selection, spin-flip operation selection, and final `KPOINTS`
writing.

| Step | What it does | Input needed |
|------|--------------|--------------|
| **0** | Finds spin-flip symmetry operations and prints a compact symmetry summary | Structure file; Cartesian spin axis and magnetic moments for non-mcif inputs |
| **1** | Builds or reads the high-symmetry IBZ path | Press Enter for auto path, or enter a KPOINTS-style file |
| **2** | Chooses the general k point | Automatic IBZ centroid by default |
| **3** | Selects the spin-flip operation | Press Enter for default, enter a number, type `list`, or type `manual` |
| **4** | Builds the altermagnetic path | Automatic |
| **5** | Saves the output | Output code (`vasp` or `qe`); filenames are fixed |

## Step 0: Spin Symmetry

Reads the structure and magnetic moments, then calls `findspingroup` to
identify the altermagnetic phase, oriented spin space group, and spin-flip/
spin-preserving operations. Supported formats: `POSCAR`, `.vasp`, `.cif`
(moments entered manually), `.mcif` (moments read from the file).

For non-`.mcif` input, enter a Cartesian spin axis (default `0 0 1`) and
scalar moments along it in atom order, VASP `MAGMOM`-style (`1 -1` or
`5*0 2*1.0`); untyped atoms default to `0`.

## Step 1: High-Symmetry K-Path

Press Enter for the auto-detected SeeK-path, or provide a line-mode
`KPATH.in`/KPOINTS-style file to start from a custom path. The detected
lattice type is reported using SeeK-path keys (`hP2`, `oI3`, `mC2`, etc.).

A custom path is interpreted in the reciprocal fractional basis of the
structure submitted in Step 0. AlterSeeK-Path converts it to the standardized
SeeK-path basis used for internal processing and converts the final path to
the calculation cell's reciprocal basis when writing output. Custom files
must use reciprocal line mode and contain labeled endpoint pairs; one complete
segment is valid, but a lone point, an odd endpoint count, a missing label, or
a non-finite coordinate is rejected.

## Step 2: General K Point

Automatic: the IBZ centroid, chosen because it's always away from special
symmetry lines/planes. No prompt in the normal case.

## Step 3: Spin-Flip Operation

Lists available spin-flip operations with a default. Accept the default,
pick another by number, `list` the matrices, or `manual` for a custom
transformation.

## Step 4: Build The Altermagnetic Path

The selected operation maps k to k'; both are inserted into the path. Path
construction uses the standardized SeeK-path primitive reciprocal basis
internally; written `KPOINTS` coordinates are converted to the calculation
cell's reciprocal basis at output time.

## Step 5: Save Outputs

Choose VASP (the default) or Quantum ESPRESSO. The workflow writes the fixed
filename `KPOINTS_alter` or `KPOINTS_alter_qe`; it does not ask for an output
filename. It also creates or updates the corresponding plotting configuration,
`alterband.toml` or `alterband_qe.toml`. See Output Files below for the full
list, including figures.

## Output Files

Files the user must act on are written to the working directory; diagnostic and
record files go into `alterseek_output/`.

**Working directory**

| File | Description |
|------|-------------|
| `KPOINTS_alter` | Altermagnetic k-path for VASP line-mode band calculations |
| `KPOINTS_alter_qe` | Altermagnetic k-path in QE `K_POINTS crystal_b` format |
| `alterband.toml` | VASP band-plot configuration; the band plotter reads it from here |
| `alterband_qe.toml` | QE band-plot configuration; likewise |
| `*_magnetic_primitive.vasp` | Magnetic primitive calculation cell in the SSG-adapted setting, written only when the calculation setting changes |
| `*_magnetic_primitive_MAGMOM.txt` | Matching species order, Cartesian spin axis, and scalar MAGMOM line for that generated POSCAR |

**`alterseek_output/`**

| File | Description |
|------|-------------|
| `*_ibz_*.png` | IBZ/BZ figure with the selected general k point |
| `*_spinflip_*.png` | Spin-up/spin-down IBZ connection figure |
| `*_spinbz_*.png` | Spin-colored BZ figure |
| `*_spinbz_top_*.png` | Top-view spin-colored BZ figure |
| `spin_operations.txt` | Full spin-symmetry operation log |
| `spin_flip_operations.txt` | Spin-flip rotation matrices used by the main workflow |
| `spin_preserve_operations.txt` | Spin-preserving rotation matrices used for completion and diagnostics |
| `*_magnetic_primitive.mcif` | Magnetic primitive cell with its vector moments (default route) |
| `*_seekpath_standard.vasp` | SeeK-path standardized conventional diagnostic cell; it is not the calculation POSCAR or final KPOINTS basis |
| `*_seekpath_standard.mcif` | The same SeeK-path standardized cell with transformed vector moments for magnetic-structure visualization; written with explicit sites in P1 |
| `*_seekpath_marker_helper.vasp` | Marker-bearing standardized helper, retained only with verbose output |
| `*_seekpath_basis_mapping.txt` | Analysis input, internal primitive path, standardized conventional, and final KPOINTS output lattices |

For Laue groups `-1`, `-3`, and `m-3`, no altermagnetic splitting is supported.
The workflow prints a note and writes the ordinary default path. The same
happens whenever no spin-flip point operation survives for other reasons
(e.g. every candidate is trivial in-plane in 2D mode) -- it means "not
altermagnetic" for this configuration, not a request for a manual matrix.

## Skipping Repeated Prompts With `alterseek_input.toml`

If an `alterseek_input.toml` file is present in the working directory,
`alterseek-path` reads it and skips the interactive prompt for any key it
finds. Any key can be omitted; the prompt for that step then runs as usual.

```toml
structure = "POSCAR"
spin_axis = "0 0 1"
moments = "5 -5"
path = ""
flip_option = 1
output_code = "vasp"
save_pdf = true
```

| Key | Matches | Notes |
|-----|---------|-------|
| `structure` | Step 0 structure-file prompt | default `POSCAR` if omitted |
| `spin_axis` | Step 0 spin-axis prompt | ignored for `.mcif` input |
| `moments` | Step 0 magnetic-moments prompt | the one field that changes every run |
| `path` | Step 1 path-choice prompt | `""` (or omit) uses the auto-generated path; a filename loads a custom `KPATH.in`/KPOINTS-style path |
| `flip_option` | Step 3 spin-flip operation prompt | plain integer, picks that numbered option; omit for the interactive numbered menu (`list`/`manual` still available) |
| `output_code` | Step 5 output-code prompt | `"vasp"` or `"qe"` |
| `save_pdf` | BZ figure output format | `true`/`false` (default `false`); when true, also saves a vector PDF alongside each BZ figure's default PNG output |
| `symprec` | symmetry-detection tolerance | positive number in angstrom (default `1e-3`); applies to every input format |

### `symprec`

The default is `1e-3` A rather than spglib's own `1e-5`. Deposited structures
routinely carry coordinates rounded to five decimals, and at `1e-5` that
rounding noise hides symmetry that is really present -- MnSe2's genuinely cubic
`Pa-3` parent reads as orthorhombic `Pbca`, which then decides whether the
structure is treated as altermagnetic at all. `1e-3` A stays far below any
deliberate distortion (a 0.5% strain on a 4 A lattice is 0.02 A, twenty times
larger) and changes none of the reference cases.

Lower it if you are studying a structure with a real distortion smaller than
`1e-3` A and need that distortion resolved rather than averaged away.

The setting applies to every input format. For `.mcif` input that declares a
parent space group, AlterSeeK-Path additionally tries the tighter sequence
`1e-5`, `1e-4`, `1e-3` and accepts the tightest one that reproduces the
declared parent; `symprec` is the fallback when there is no declaration to
validate against. Note that this affects only the nonmagnetic parent. The
altermagnetism check in the default (magnetic primitive cell) route reads G0
from the spin space group and uses no tolerance at all.

Comment out a key (or delete the line) to make that one step interactive
again while the rest of the file still drives the run.

The file is validated before the workflow starts. Malformed TOML, unknown
keys, invalid types, an unsupported `output_code`, or an out-of-range
`flip_option` stops the run with an error instead of silently falling back to
interactive input.

## 2D / Slab Mode

For 2D materials computed as slabs (vacuum along one lattice vector), run:

```bash
alterseek-path --2d
```

The full 3D symmetry analysis runs unchanged; output is restricted to the
physical in-plane (vacuum k = 0) reciprocal plane -- k-path, general k point
(2D IBZ area centroid), and all written k-points. Step 3 additionally reports
a verdict:

```text
[2D mode] In-plane spin splitting: YES
```

Operations that are trivial in-plane (identity or k to -k) don't count. If
none remain, the verdict is `NO` and the ordinary in-plane path is written
without `k'`.

The vacuum axis is auto-detected in the standardized cell regardless of
`--vacuum-axis {a,b,c}` (default `c`), which only sanity-checks the input
cell. The workflow warns if the input has no clear vacuum gap or a tilted
vacuum axis.

Instead of the 3D figures, 2D mode saves top-down figures:

| File | Description |
|------|-------------|
| `*_2d_ibz_*.png` | 2D BZ outline, in-plane path, and the general k point |
| `*_2d_spinflip_*.png` | spin-up IBZ and its spin-flip image with path connections |
| `*_2d_spinbz_*.png` | spin-colored 2D BZ domain pattern |

## `--parent-setting`

By default, Figure 1 and the path symmetry are determined from the
**magnetic primitive cell** -- the cell that reflects the symmetry of the
magnetic state -- which AlterSeeK-Path identifies from the input magnetic
structure with FindSpinGroup's spin-space-group setting. Final KPOINTS
coordinates are expressed in the reciprocal basis of the calculation cell.

That cell coincides with the ordinary nonmagnetic structure when the magnetic
order does not lower the lattice symmetry. It differs when the order does
lower it: either a **q = 0** order selects a different magnetic primitive cell
of the same volume (e.g. a hexagonal structural lattice whose magnetic pattern
only respects orthorhombic symmetry), or a nonzero propagation vector
**q != 0** enlarges it (e.g. a magnetic cell containing three parent cells).
If the submitted calculation cell is itself an integer supercell of that
magnetic primitive cell, AlterSeeK-Path keeps it and maps the path back into
its reciprocal basis. A same-volume nontrivial basis change instead produces
the matching `*_magnetic_primitive.vasp` and MAGMOM companion for the band
calculation.
If construction of the required magnetic primitive cell fails, AlterSeeK-Path
stops without writing a parent-cell replacement path. Fix the reported input or
dependency problem, or explicitly rerun with `--parent-setting` when a
nonmagnetic reference path is actually intended.

```bash
alterseek-path --parent-setting
```

builds Figure 1/KPOINTS in the nonmagnetic parent cell instead. Displaying
bands along a common parent-cell reference path is useful for comparing
different magnetic orders while keeping the path fixed, but that reference
path does not respect the symmetry of the magnetic state. Where the parent is
not altermagnetic at all, it cannot represent the magnetic state, and the path
built in it will not show the splitting the magnetic structure actually has.

Add `--output verbose` to keep the intermediate helper files for debugging
instead of deleting them.

The altermagnetism check follows the cell in use. By default the Laue-group
exclusion is evaluated on **G0**, the spatial part of the spin space group
reported by FindSpinGroup, which is the symmetry the magnetic primitive cell
actually has; under `--parent-setting` it is evaluated on the submitted cell's
own space group. So a structure whose nonmagnetic parent forbids altermagnetism
but whose magnetic order lowers the symmetry to a Laue class that permits it is
correctly treated as altermagnetic by default.

This distinction is not cosmetic for supercell altermagnets. MnSe2 is deposited
in a cubic `Pa-3` (205) parent, Laue `m-3`, which forbids altermagnetism; its
magnetic primitive cell is a 3x1x1 supercell of that same cubic crystal, so
with the moments stripped it is *still* cubic. Only G0 (`Pbca`, 61, Laue `mmm`)
reflects the symmetry the magnetic state actually has, and Step 0 prints both
the nonmagnetic parent and G0 so the two are never confused.

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.
