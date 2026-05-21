# Workflow

The main command is interactive:

```bash
alterseek-path
```

It proceeds through structure reading, spin-symmetry analysis, path selection,
general k-point selection, spin-flip operation selection, and final `KPOINTS`
writing.

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

## Band Plotting

After the VASP band calculation, run VASPKIT task `303`, then plot with:

```bash
alterseek-path bandplot
```

By default, the plotter reads:

```text
KLABELS
REFORMATTED_BAND_UP.dat
REFORMATTED_BAND_DW.dat
```

and writes:

```text
alterband.png
```

Command-line options can override the TOML settings:

```bash
alterseek-path bandplot -o alterband.pdf --emin -2 --emax 2
```

---

[Back to user guide](index.md)
