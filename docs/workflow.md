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
| **3** | Selects a detected spin-flip operation | Press Enter for default, enter a number, or type `list` |
| **4** | Builds the altermagnetic path | Automatic |
| **5** | Saves the output | Output code (`vasp` or `qe`); filenames are fixed |

## Step 0: Spin Symmetry

Reads the structure and magnetic moments, then calls `findspingroup` to
identify the altermagnetic phase, oriented spin space group, and spin-flip/
spin-preserving operations. Supported formats: `POSCAR`, `.vasp`, `.cif`
(moments entered manually), `.mcif` (moments read from the file).

The current workflow supports collinear magnetism only. For `.mcif` input,
zero moments are ignored and every nonzero vector moment must be parallel or
antiparallel to one common axis within an absolute transverse-moment tolerance
of `0.02` in the MCIF moment units (normally Bohr magnetons). This matches
FindSpinGroup's default moment-tolerance scale and accommodates rounded
deposited values. A noncollinear MCIF is rejected instead of projecting its
moments onto a collinear axis.

For non-`.mcif` input, enter a Cartesian spin axis (default `0 0 1`) and
scalar moments along it in atom order, VASP `MAGMOM`-style (`1 -1` or
`5*0 2*1.0`). Missing trailing values default to `0`; providing more moments
than structure sites is rejected as an atom-order/count error.

## Step 1: High-Symmetry K-Path

The path is always generated from the detected symmetry; the detected
lattice type is reported using SeeK-path keys (`hP2`, `oI3`, `mC2`, etc.).
Drop whichever generated segments you do not want from the written KPOINTS
file. The path is converted to the submitted cell's reciprocal basis when
writing output.

## Step 2: General K Point

Automatic: the IBZ centroid, chosen because it's always away from special
symmetry lines/planes. No prompt in the normal case.

## Step 3: Spin-Flip Operation

Lists the detected spin-flip operations with a default. Accept the default,
pick another by number, or `list` the matrices. Every selectable operation
comes from the detected symmetry set.

## Step 4: Build The Altermagnetic Path

The selected operation maps k to k'; both are inserted into the path. Path
construction uses the standardized SeeK-path primitive reciprocal basis
internally; written `KPOINTS` coordinates are converted to the submitted
cell's reciprocal basis at output time.

## Step 5: Save Outputs

Choose VASP (the default), Quantum ESPRESSO, or ABINIT. The workflow writes
the fixed filename `KPOINTS_alter`, `KPOINTS_alter_qe`, or
`KPOINTS_alter_abinit`; it does not ask for an output filename. It also
creates or updates the corresponding plotting configuration —
`alterseek_plot_vasp.toml`, `alterseek_plot_qe.toml`, or `alterseek_plot_abinit.toml`. See
Output Files below for the full list, including figures.

## Output Files

Files the user must act on are written to the working directory; diagnostic and
record files go into `alterseek_output/`.

**Working directory**

| File | Description |
|------|-------------|
| `KPOINTS_alter` | Altermagnetic k-path for VASP line-mode band calculations |
| `KPOINTS_alter_qe` | Altermagnetic k-path in QE `K_POINTS crystal_b` format |
| `KPOINTS_alter_abinit` | Altermagnetic k-path as an ABINIT `kptopt`/`ndivk`/`kptbounds` block |
| `alterseek_plot_vasp.toml` | VASP band-plot configuration; the band plotter reads it from here |
| `alterseek_plot_qe.toml` | QE band-plot configuration; likewise |
| `alterseek_plot_abinit.toml` | ABINIT band-plot configuration; likewise |

**`alterseek_output/`**

| File | Description |
|------|-------------|
| `*_ibz_*.png` | IBZ/BZ figure with the selected general k point |
| `*_spinflip_*.png` | Spin-up/spin-down IBZ connection figure |
| `*_spinbz_*.png` | Spin-colored BZ figure |
| `*_spinbz_top_*.png` | Top-view spin-colored BZ figure |
| `spin_operations.txt` | Full spin-symmetry operation log, written only when the current run performs spin analysis; a no-moments run neither creates nor changes it |
| `spin_flip_operations.txt` | Spin-flip rotation matrices used by the main workflow |
| `spin_preserve_operations.txt` | Spin-preserving rotation matrices used for completion and diagnostics |
| `*_magnetic_primitive.mcif` | Direct FindSpinGroup magnetic primitive reference with vector moments; it is not the calculation cell |
| `*_seekpath_basis_mapping.txt` | Submitted analysis/output lattice, internal primitive path lattice, standardized conventional lattice, and rotation mapping |

## Submitted-Cell Analysis Contract

The submitted structure remains the calculation and output cell. Its volume is
first compared with the physical primitive translation cell. An integer volume
index of one is a change of basis and keeps the ordinary physical `(R|t)`
analysis. An index above one activates the conventional/supercell mode: the
submitted edge-vector translations, without hidden centering or smaller-cell
translations, define the reciprocal lattice and calculation-cell BZ. The
submitted metric and compatible point rotations then determine the IBZ,
centroid, high-symmetry path, and Figures 1-4. A conventional cI or cF cube has
a cP calculation-cell BZ; a conventional R-centered hexagonal cell has an hP
BZ; and an axis-aligned cubic `2 x 2 x 1` supercell has a tetragonal BZ.
AlterSeeK-Path never emits a replacement calculation structure or MAGMOM file.

For magnetic input, FindSpinGroup determines physical G0, its magnetic
primitive cell, setting transforms, and spin operations. For volume index one,
AlterSeeK-Path constructs the ordinary helper from the complete physical Seitz
set, including screw and glide columns. For an index above one, it transforms
the G0 point parts into the submitted fractional basis, retains only integral
unimodular rotations that preserve the submitted metric, and verifies
multiplication closure. It applies each retained rotation to deterministic
generic seeds as `Rr` in a marker-only in-memory SeeK-path proxy. Each seed
orbit has a distinct artificial type so accidental orbit exchange cannot add
symmetry. Spglib must recover exactly the intended `(R|0)` set, no fractional
translations, and `volume_original_wrt_prim = 1` before the proxy is accepted.

The nonprimitive-cell proxy is a reciprocal-space/path construction, not a
relabeling of the physical crystal. Complete physical `(R|t)` operations
remain separate for spin transformations, nonsymmorphic phases, and diagnostics. SeeK-path may
standardize the proxy internally; all generated points are converted through
Cartesian reciprocal space into the submitted VASP, QE, or ABINIT basis. Figures use
the corresponding standardized-basis spin matrices, while output coordinates
and operation logs retain the submitted basis. Both representations describe
the same `k' = R^(-T) k` mapping. The same contract applies to magnetic,
q != 0, and no-moments inputs.

The console reports `Input cell`, `Nonmagnetic primitive cell`, and, when
magnetic moments are present, `Magnetic primitive cell`. The no-moments route
prints the same summary without the magnetic row. For a conventional/supercell
input, a detection notice directly below `Input structure` states how many
magnetic or nonmagnetic primitive cells it contains, and
`Conventional/supercell BZ` reports the new BZ lattice tag. The artificial
proxy space-group symbol is not printed.

The marker cell is never written. Because its SeeK-path standard
would not be an honest transformed real structure, the workflow also does not
publish `*_seekpath_standard.vasp`. The direct FindSpinGroup
`*_magnetic_primitive.mcif` remains as an internal-symmetry reference, and the
basis-mapping diagnostic records the submitted analysis/output lattice
explicitly.

For Laue groups `-1`, `-3`, and `m-3`, no altermagnetic splitting is supported.
The workflow prints a note and writes the ordinary default path without an
operation-selection step. The same expected fallback applies when spin-symmetry
analysis classifies the configuration as non-altermagnetic, or when every
candidate is trivial in-plane in 2D mode. If analysis instead reaches the
altermagnetic Step 3 path but no detected operation is available, the results
are inconsistent; the workflow reports an error and aborts rather than writing
a misleading ordinary path. There is no custom-matrix fallback.

## Skipping Repeated Prompts With `alterseek_input.toml`

If an `alterseek_input.toml` file is present in the working directory,
`alterseek-path` reads it and skips the interactive prompt for any key it
finds. Any key can be omitted; the prompt for that step then runs as usual.

```toml
structure = "POSCAR"
spin_axis = "0 0 1"
moments = "5 -5"
flip_option = 1
output_code = "vasp"
save_pdf = true
```

| Key | Matches | Notes |
|-----|---------|-------|
| `structure` | Step 0 structure-file prompt | default `POSCAR` if omitted |
| `spin_axis` | Step 0 spin-axis prompt | ignored for `.mcif` input |
| `moments` | Step 0 magnetic-moments prompt | the one field that changes every run |
| `flip_option` | Step 3 spin-flip operation prompt | plain integer, picks that numbered detected operation; omit for the interactive numbered menu (`list` remains available) |
| `output_code` | Step 5 output-code prompt | `"vasp"` or `"qe"` |
| `save_pdf` | BZ figure output format | `true`/`false` (default `false`); when true, also saves a vector PDF alongside each BZ figure's default PNG output |
| `symprec` | symmetry-detection tolerance | positive number in angstrom (default `1e-3`); applies to every input format |
| `vacuum_axis` | 2D-slab vacuum axis of the input cell | `"a"`, `"b"` or `"c"` (default `"c"`); 2D mode only, overridden by an explicit `--vacuum-axis` |

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
altermagnetism check reads G0 from the spin space group and uses no tolerance
at all; `symprec` still controls validation of the submitted-cell marker helper.

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
(2D IBZ area centroid), and all written k-points. Step 3 filters the detected
spin-flip operations to those that act nontrivially within that plane.

Operations that are trivial in-plane (identity or k to -k) don't count. If
none remain, the verdict is `NO` and the ordinary in-plane path is written
without `k'`.

`--vacuum-axis {a,b,c}` (default `c`) identifies the submitted cell's vacuum
vector. The other two submitted lattice vectors define one physical Cartesian
slab plane, which is carried through the submitted, SeeK-path-standardized,
and final KPOINTS bases. Spin operations are tested in
Cartesian reciprocal space, so basis standardization cannot change
the 2D verdict. Written k-points are projected only to remove numerical residue
normal to that same plane; a genuinely out-of-plane point aborts instead of
being silently flattened. SeeK-path still auto-detects its own standardized
vacuum-axis index for the 2D hull and path. The workflow warns if the input has
no clear vacuum gap or a tilted vacuum axis.

Instead of the 3D figures, 2D mode saves top-down figures:

| File | Description |
|------|-------------|
| `*_2d_ibz_*.png` | 2D BZ outline, in-plane path, and the general k point |
| `*_2d_spinflip_*.png` | spin-up IBZ and its spin-flip image with path connections |
| `*_2d_spinbz_*.png` | spin-colored 2D BZ domain pattern |

For supercell altermagnets, the submitted lattice fixes the BZ while G0 fixes
the magnetic point symmetry inside it. MnSe2 illustrates why those roles must
remain separate: the submitted magnetic cell has a cubic moment-free parent,
but G0 is orthorhombic `Pbca` (61, Laue `mmm`). Step 0 therefore reports the
nonmagnetic primitive cell, the FindSpinGroup magnetic primitive cell, and the
conventional/supercell BZ helper as distinct concepts.

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.
