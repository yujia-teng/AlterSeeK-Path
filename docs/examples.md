# Examples

## Minimal Interactive Run

Run the main workflow in a directory containing a structure file:

```bash
alterseek-path
```

Example prompts:

```text
Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): GdAuGe.vasp
Magnetic moments (atom order, trailing atoms auto-fill to 0): 1 -1
```

Press Enter in Step 1 to use the automatically detected SeeK-path.

Press Enter in Step 2 to use the IBZ centroid as the general k point.

Press Enter in Step 3 to use the default spin-flip operation when one is
available.

The default output is:

```text
KPOINTS_modified
```

## Manual Moment Syntax

Magnetic moments use VASP `MAGMOM` style syntax:

```text
1 -1
```

Repeated values are supported:

```text
5*0 2*1.0
```

If fewer values are entered than the number of atoms, trailing atoms are filled
with zero moment.

## Custom Starting Path

In Step 1, instead of pressing Enter, provide a line-mode path file:

```text
KPATH.in
```

AlterSeeK-Path then inserts k and k' into that starting path rather than the
automatically generated default path.

## Band Plot Example

After VASP finishes and VASPKIT task `211` has generated the reformatted band
files, run:

```bash
alterseek-path bandplot
```

To write PDF output:

```bash
alterseek-path bandplot -o alterband.pdf
```

To adjust the energy window:

```bash
alterseek-path bandplot --emin -3 --emax 3
```

For detailed plotting settings, see [Plotting](plotting.md).
