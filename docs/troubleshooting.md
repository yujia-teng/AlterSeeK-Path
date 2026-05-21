# Troubleshooting

## `alterseek-path` Is Not Found

Install the package in editable mode from the repository root:

```bash
pip install -e .
```

Then try again:

```bash
alterseek-path
```

## Dependency Import Errors

Install the listed requirements:

```bash
pip install -r requirements.txt
```

If the environment is old or mixed with other scientific packages, using a fresh
virtual environment can avoid version conflicts.

## No Spin-Flip Operation Is Found

If no valid spin-flip operation is found, the magnetic pattern may not describe
a collinear altermagnetic case in the expected scalar up/down form.

Check:

- the structure file,
- atom order,
- magnetic moment signs,
- whether the moment pattern matches the intended magnetic symmetry.

For `.mcif` inputs, verify that the extracted moment signs match the POSCAR atom
order expected by the workflow.

## The Generated KPOINTS Look Different From SeeK-Path Coordinates

AlterSeeK-Path constructs paths internally in the standardized SeeK-path/HPKOT
primitive reciprocal basis.

VASP reads `KPOINTS` in the reciprocal basis of the actual input structure. When
the input structure basis differs from the standardized primitive basis,
AlterSeeK-Path converts the final numerical coordinates before writing the
VASP file.

This is expected and preserves the same Cartesian reciprocal-space k vectors.

## Band Plot Files Are Missing

Before running:

```bash
alterseek-path bandplot
```

make sure VASPKIT task `303` has produced:

```text
KLABELS
REFORMATTED_BAND_UP.dat
REFORMATTED_BAND_DW.dat
```

If your files have different names, pass them explicitly with command-line
options such as `--klabels`, `--up`, and `--down`.

---

[Back to user guide](index.md)
