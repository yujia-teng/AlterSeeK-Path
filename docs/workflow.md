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
| **5** | Saves the output | Output filename |

## Step 0: Spin Symmetry

AlterSeeK-Path reads the structure and magnetic moments, then calls
`findspingroup` to identify the altermagnetic phase, report the oriented spin
space group index and Chen-Liu symbol, and extract spin-flip and
spin-preserving operations.

Supported structure inputs:

- `POSCAR`
- `.vasp`
- `.cif`
- `.mcif`

For non-magnetic structure formats, enter the collinear spin axis in Cartesian
coordinates. The default is `0 0 1`, meaning Cartesian z. Then enter scalar
moments along that axis in atom order using VASP `MAGMOM` syntax. For example:

```text
1 -1
```

or:

```text
5*0 2*1.0
```

Untyped atoms default to zero moment.

The entered axis is normalized and multiplied by each scalar moment before it
is passed to FindSpinGroup. It is a Cartesian spin vector; it is not converted
to fractional lattice coordinates. For `.mcif` files, moments are read directly
from the file and this manual-axis prompt is skipped.

For collinear structures, spin-flip operations are identified by the action
of each spin-space operation on the collinear spin axis. The written spatial
rotation matrices remain in the input-structure basis for use by subsequent
path-generation steps.

The structural lattice type determines the Brillouin-zone boundary, the Figure
1 IBZ hull, and HPKOT path. Figures 3 and 4 map that exact Figure 1 hull using
FindSpinGroup spin-preserving and spin-flipping operations. If a magnetic state
has lower spin symmetry than its structural lattice, these mapped hulls need
not fill the whole structural BZ.

## Step 1: High-Symmetry K-Path

Press Enter to use the automatically detected SeeK-path.

Alternatively, provide a line-mode `KPATH.in` or KPOINTS-style file if you want
to start from a custom path.

The detected lattice type is reported using SeeK-path lattice keys such as
`hP2`, `oI3`, or `mC2`.

## Step 2: General K Point

By default, AlterSeeK-Path uses the centroid of the irreducible Brillouin zone
as the general k point.

This point is inserted into the high-symmetry path so the band calculation
samples spin splitting away from special symmetry lines.

## Step 3: Spin-Flip Operation

The workflow lists available spin-flip operations and chooses a default when
possible.

You can accept the default, list the available matrices, select another
operation, or enter a manual transformation.

## Step 4: Build The Altermagnetic Path

The selected operation maps k to k'. AlterSeeK-Path then inserts k and k' into
the path so paired spin-splitting segments are sampled systematically.

Internally, the path construction uses the standardized SeeK-path primitive
reciprocal basis. At output time, the written `KPOINTS` coordinates
are converted into the reciprocal basis of the input structure used by VASP.

## Step 5: Save Outputs

The default output file is:

```text
KPOINTS_modified
```

The workflow also updates:

```text
alterband.toml
```

This file records band-plot settings, including the detected lattice type when
available.

Depending on the plotting settings, the workflow may also save Brillouin-zone
and spin-BZ figures.

## Output Files

| File | Description |
|------|-------------|
| `KPOINTS_modified` | Altermagnetic k-path for VASP line-mode band calculations |
| `alterband.toml` | Band-plot configuration written by the main workflow |
| `spin_operations.txt` | Full spin-symmetry operation log |
| `spin_flip_operations.txt` | Spin-flip rotation matrices used by the main workflow |
| `spin_preserve_operations.txt` | Spin-preserving rotation matrices used for completion and diagnostics |
| `*_ibz_*.png` | IBZ/BZ figure with the selected general k point |
| `*_spinflip_*.png` | Spin-up/spin-down IBZ connection figure |
| `*_spinbz_*.png` | Spin-colored BZ figure |
| `*_spinbz_top_*.png` | Top-view spin-colored BZ figure |

For Laue groups `-1`, `-3`, and `m-3`, no altermagnetic splitting is supported.
The workflow prints a note and writes the ordinary default path.

## 2D / Slab Mode

For 2D materials computed as slabs (vacuum along one lattice vector), run:

```bash
alterseek-path --2d
```

2D mode keeps the full 3D symmetry analysis unchanged and restricts the
output to the physical k-plane: the k-path, the general k point (now the 2D
IBZ area centroid), and all written k-points lie in the in-plane
(vacuum k = 0) reciprocal plane.

Step 3 additionally classifies each spin-flip operation by its in-plane
action and reports a verdict:

```text
[2D mode] In-plane spin splitting: YES
```

Operations whose in-plane action is trivial (identity or a simple k to -k
map) produce no in-plane splitting. If no valid operation remains, the
verdict is `NO` and the ordinary in-plane path is written without `k'`.

The vacuum direction is auto-detected in the standardized cell.
`--vacuum-axis {a,b,c}` (default `c`) describes the vacuum axis of the input
cell and is used for an input sanity check; the slicing axis is detected
automatically regardless of this flag. The slab should follow the usual
2D-structure convention, with the vacuum axis perpendicular to the layer
plane. The workflow warns if the input cell has no clear vacuum gap or a
tilted vacuum axis.

Instead of the 3D figures, 2D mode saves top-down figures:

| File | Description |
|------|-------------|
| `*_2d_ibz_*.png` | 2D BZ outline, in-plane path, and the general k point |
| `*_2d_spinflip_*.png` | spin-up IBZ and its spin-flip image with path connections |
| `*_2d_spinbz_*.png` | spin-colored 2D BZ domain pattern |

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.
