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
| **0** | Finds spin-flip symmetry operations and prints a compact symmetry summary | Structure file; magnetic moments for non-mcif inputs |
| **1** | Builds or reads the high-symmetry IBZ path | Press Enter for auto path, or enter a KPOINTS-style file |
| **2** | Chooses the general k point | Automatic IBZ centroid by default |
| **3** | Selects the spin-flip operation | Press Enter for default, enter a number, type `list`, or type `manual` |
| **4** | Builds the altermagnetic path | Automatic |
| **5** | Saves the output | Output filename |

## Step 0: Spin Symmetry

AlterSeeK-Path reads the structure and magnetic moments, then calls `spinspg`
to identify spin-flip and spin-preserving operations.

Supported structure inputs:

- `POSCAR`
- `.vasp`
- `.cif`
- `.mcif`

For non-magnetic structure formats, enter moments in atom order using VASP
`MAGMOM` syntax. For example:

```text
1 -1
```

or:

```text
5*0 2*1.0
```

Untyped atoms default to zero moment.

## Step 1: High-Symmetry K-Path

Press Enter to use the automatically detected SeeK-path/HPKOT path.

Alternatively, provide a line-mode `KPATH.in` or KPOINTS-style file if you want
to start from a custom path.

The detected lattice type is reported using HPKOT-style keys such as `hP2`,
`oI3`, or `mC2`.

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

Internally, the path construction uses the standardized SeeK-path/HPKOT
primitive reciprocal basis. At output time, the written `KPOINTS` coordinates
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

## Next Step

After running the VASP band calculation, see [Plotting](plotting.md) for
spin-resolved band-plot generation.

---

[Back to user guide](index.md)
