# AlterSeeK-Path User Guide

AlterSeeK-Path generates general-k paths for altermagnetic band-structure
calculations. It starts from a standard SeeK-path high-symmetry path,
inserts a general k point, maps it through a spin-flip operation to k', and
writes a VASP, Quantum ESPRESSO, or ABINIT k-path file for sampling spin
splitting along paired path segments.

This guide expands the quick-start information in the repository README.

## Main Commands

Run the k-path workflow:

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
alterseek-plot
```

The generated plotting configuration selects the code automatically when it
is the only one present. It can also be specified explicitly:

```bash
alterseek-plot vasp
alterseek-plot qe
alterseek-plot abinit
```

See [Plotting](plotting.md) for their inputs and settings.

Standalone utilities are also available:

```bash
python -m alterseek.find_sf_operations
python -m alterseek.compute_centroid_3d POSCAR
```

## Coordinate Conventions

SeeK-path may construct the path in a standardized reciprocal basis.
AlterSeeK-Path converts the resulting fractional coordinates into the
submitted cell's reciprocal basis while preserving the same physical
k-vector. VASP, QE, and ABINIT outputs therefore use fractional reciprocal
coordinates in the submitted cell basis.
