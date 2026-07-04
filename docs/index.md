# AlterSeeK-Path User Guide

AlterSeeK-Path generates k-point paths for altermagnetic band-structure
calculations. It starts from a standard SeeK-path high-symmetry path,
inserts a general k point, maps it through a spin-flip operation to k', and
writes a VASP `KPOINTS` file for sampling spin splitting along paired path
segments.

This guide expands the quick-start information in the repository README.

## Main Commands

Run the interactive k-path workflow:

```bash
alterseek-path
```

Run the same workflow in 2D/slab mode (see [Workflow](workflow.md) for
details):

```bash
alterseek-path --2d
```

Plot a generated spin-resolved band structure:

```bash
alterseek-path bandplot
```

See [Plotting](plotting.md) for band-plot inputs, configuration, and output
options.

Standalone utilities are also available:

```bash
python -m alterseek.find_sf_operations
python -m alterseek.compute_centroid_hybrid POSCAR
```

## Coordinate Conventions

The internal path conventions follow SeeK-path labels and reciprocal
bases. When writing VASP `KPOINTS`, the final numerical coordinates are
converted into the reciprocal basis of the actual input structure used by VASP.
