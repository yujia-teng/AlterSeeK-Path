# AlterSeeK-Path User Guide

AlterSeeK-Path generates k-point paths for altermagnetic band-structure
calculations. It starts from a standard SeeK-path high-symmetry path,
inserts a general k point, maps it through a spin-flip operation to k', and
writes a VASP `KPOINTS` file for sampling spin splitting along paired path
segments.

This guide expands the quick-start information in the repository README.

## Contents

- [Installation](installation.md)
- [Workflow](workflow.md)
- [Plotting](plotting.md)
- [Examples](examples.md)

## Main Commands

Run the interactive k-path workflow:

```bash
alterseek-path
```

Plot a generated spin-resolved band structure:

```bash
alterseek-path bandplot
```

See [Plotting](plotting.md) for band-plot inputs, configuration, and output
options.

Standalone utilities are also available:

```bash
python find_sf_operations.py
python compute_centroid_hybrid.py POSCAR
```

## Current Scope

AlterSeeK-Path currently focuses on VASP workflows. Quantum ESPRESSO support is
partial.

The internal path conventions follow SeeK-path labels and reciprocal
bases. When writing VASP `KPOINTS`, the final numerical coordinates are
converted into the reciprocal basis of the actual input structure used by VASP.
