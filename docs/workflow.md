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
structure submitted in Step 0. With `--ssg-setting`, the analysis works in
the magnetic primitive cell instead (e.g. a hexagonal lattice whose moments
lower it to orthorhombic), so the custom path is interpreted in that magnetic
primitive cell's reciprocal basis. AlterSeeK-Path converts it to the standardized
SeeK-path basis for internal processing and converts the final path to the
calculation cell's reciprocal basis when writing output. Custom files must use
reciprocal line mode and contain labeled endpoint pairs; one complete segment
is valid, but a lone point, an odd endpoint count, a missing label, or a
non-finite coordinate is rejected.

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

| File | Description |
|------|-------------|
| `KPOINTS_alter` | Altermagnetic k-path for VASP line-mode band calculations |
| `KPOINTS_alter_qe` | Altermagnetic k-path in QE `K_POINTS crystal_b` format |
| `alterband.toml` | VASP band-plot configuration written by the main workflow |
| `alterband_qe.toml` | QE band-plot configuration written by the main workflow |
| `spin_operations.txt` | Full spin-symmetry operation log |
| `spin_flip_operations.txt` | Spin-flip rotation matrices used by the main workflow |
| `spin_preserve_operations.txt` | Spin-preserving rotation matrices used for completion and diagnostics |
| `*_ibz_*.png` | IBZ/BZ figure with the selected general k point |
| `*_spinflip_*.png` | Spin-up/spin-down IBZ connection figure |
| `*_spinbz_*.png` | Spin-colored BZ figure |
| `*_spinbz_top_*.png` | Top-view spin-colored BZ figure |

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

## `--ssg-setting`

```bash
alterseek-path --ssg-setting
```

By default, Figure 1/KPOINTS are built from the original structural cell's
symmetry. Magnetic order can lower that symmetry -- e.g. a hexagonal
structural lattice whose magnetic pattern only respects orthorhombic
symmetry. `--ssg-setting` instead builds Figure 1/KPOINTS from the magnetic
primitive cell (the lower, magnetically-correct symmetry) via
FindSpinGroup's spin-space-group setting, rather than the higher-symmetry
structural cell. Falls back to the structural setting with a warning if this
fails. Add `--output verbose` to keep the intermediate helper files for
debugging instead of deleting them.

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.
