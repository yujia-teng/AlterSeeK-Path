# Installation

AlterSeeK-Path requires Python 3.11 or newer.

## Install

```bash
pip install alterseek-path
```

## Install From Source

For the latest development version, clone the repository and install it:

```bash
git clone https://github.com/yujia-teng/AlterSeeK-Path.git
cd AlterSeeK-Path
pip install .
```

The installation provides two commands:

```bash
alterseek-path
alterseek-plot
```

Use `alterseek-path` to generate a path and `alterseek-plot` to plot bands.

## Input File Expectations

Run the workflow in a directory containing a `POSCAR`, `.vasp`, `.cif`, or
`.mcif` structure. Magnetic moments are entered for the first three formats and
read from `.mcif` files. Only collinear magnetism is supported; noncollinear
MCIF input is rejected. See [Workflow](workflow.md) for the full input rules.
