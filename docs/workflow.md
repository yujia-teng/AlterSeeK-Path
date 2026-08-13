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

Press Enter for the auto-detected SeeK-path, or provide a line-mode
`KPATH.in`/KPOINTS-style file to start from a custom path. The detected
lattice type is reported using SeeK-path keys (`hP2`, `oI3`, `mC2`, etc.).

A custom path is interpreted in the reciprocal fractional basis of the
structure submitted in Step 0. AlterSeeK-Path converts it to the standardized
SeeK-path basis used for internal processing and converts the final path to
the submitted cell's reciprocal basis when writing output. Custom files
must use reciprocal line mode and contain labeled endpoint pairs; one complete
segment is valid, but a lone point, an odd endpoint count, a missing label, or
a non-finite coordinate is rejected.

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

The submitted translation lattice is authoritative. It defines the BZ, IBZ,
centroid, high-symmetry path, Figures 1-4, and the reciprocal basis of final
VASP/QE output. A true integer supercell receives its own native submitted-cell
BZ rather than a primitive path folded into supercell coordinates. A
determinant-one basis change likewise remains the calculation and output cell.
AlterSeeK-Path never emits a replacement calculation POSCAR or MAGMOM file.

For magnetic input, FindSpinGroup still determines G0 and the spin operations.
Its magnetic-primitive spatial Seitz operations `(R|t)` are transformed into
the submitted basis and applied to deterministic generic seeds as `Rr+t` in an
in-memory SeeK-path cell alongside the submitted real atoms. This retains the
old He helper's symmetry construction while replacing its chemical marker with
one positive artificial type distinct from the real atom types. Every complete
operation returned in the magnetic-primitive set is retained, including screw,
glide, and centering translations. Before the helper can be used, spglib must
recover exactly the intended full spatial operation set. A
`volume_original_wrt_prim` greater than one is valid and expected for a
centered conventional setting; for example, the GdAuGe 211 `Cmc2_1` helper has
ratio two. SeeK-path determines the HPKOT labels and fractional geometry, but
those fractions are interpreted directly in the submitted reciprocal basis.
They are not converted to preserve Cartesian k-points from SeeK-path's internal
standardized primitive cell. The same submitted-basis convention is used for
magnetic, q != 0, and no-moments inputs.

The marker cell is never written. Because a marker-only SeeK-path standard
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
| `flip_option` | Step 3 spin-flip operation prompt | plain integer, picks that numbered detected operation; omit for the interactive numbered menu (`list` remains available) |
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
(2D IBZ area centroid), and all written k-points. Step 3 additionally reports
a verdict:

```text
[2D mode] In-plane spin splitting: YES
```

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
submitted analysis cell as distinct concepts.

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.
