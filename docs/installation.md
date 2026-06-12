# Installation

AlterSeeK-Path requires Python 3.11 or newer.

## Install From Source

Clone the repository and install the package in editable mode:

```bash
git clone https://github.com/yujia-teng/AlterSeeK-Path.git
cd AlterSeeK-Path
pip install -r requirements.txt
pip install -e .
```

After installation, the main command should be available:

```bash
alterseek-path
```

The band-plotting shortcut should also be available:

```bash
alterseek-bandplot
```

The same plotting command can be run through the main entry point:

```bash
alterseek-path bandplot
```

## Python Dependencies

The current dependency list is:

- `numpy`
- `matplotlib`
- `scipy`
- `sympy`
- `spglib`
- `findspingroup >= 0.15.6`
- `ase`
- `seekpath`
- `pymatgen`

These are installed through `requirements.txt` or through the package metadata
in `pyproject.toml`.

## Input File Expectations

For the main workflow, run AlterSeeK-Path in a directory containing a structure
file such as:

- `POSCAR`
- `.vasp`
- `.cif`
- `.mcif`

For `POSCAR`, `.vasp`, and `.cif` inputs, magnetic moments are entered manually.
For `.mcif` inputs, magnetic moments are read from the file when available.
