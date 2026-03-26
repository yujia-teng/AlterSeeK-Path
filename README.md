# AlterSeeK-Path

A tool for generating general k-point paths for band structure calculations of altermagnet.

It inserts a general k point into standard high-symmetry path like `Γ−M−K−Γ`:

```
Γ−M-k|k'−M'−K'−k'|k−K−Γ...
```

where `k` is a general k-point (the volume centroid of the irreducible Brillouin zone) and `k'` is obtained from spin-flip operation acting on `k`.

![band](./example/band.png)

**Current support:** VASP (full). Quantum ESPRESSO (partially).

---

## Installation

```bash
git clone https://github.com/yujia-teng/AlterSeeK-Path.git
cd AlterSeeK-Path
pip install -r requirements.txt
```

---

## Quick Start

Place your `POSCAR` and `KPATH.in` (line-mode KPOINTS file, e.g. from VASPKIT) in the same directory as the scripts, then run:

```bash
python generate_kpath.py
```

The script will guide you through 5 steps interactively:

| Step | What it does | Input needed |
|------|-------------|--------------|
| **0** | Finds spin-flip symmetry operations from structure | Structure file name, magnetic moments |
| **1** | Reads the high-symmetry k-path | KPOINTS file name |
| **2** | Auto-computes the general k-point (IBZ centroid) | *(automatic — no input needed)* |
| **3** | Selects the spin-flip transformation matrix | Choose from list, or press Enter for default |
| **4** | Generates the enriched k-path | *(automatic)* |
| **5** | Saves the output file | Output file name |

**Output:** `KPOINTS_modified` — ready to use directly in a VASP band structure calculation.

---

## Input Files

### Structure file (`POSCAR` or `cif` file)
Standard structure file.

### `KPATH.in`
Line-mode KPOINTS file with the high-symmetry path. Generate with [VASPKIT](https://vaspkit.com/) (task 303) or [seekpath](https://www.materialscloud.org/work/tools/seekpath).

> **Tip:** Use a **continuous** path (e.g. `Γ-M-K-Γ-A-L-H-A`) rather than disconnected segments (e.g. `Γ-M | H-K`). Disconnected segments cause duplicate k-points in the output.

---

## Output Files

| File | Description |
|------|-------------|
| `KPOINTS_modified` | Enriched k-path for VASP band structure calculation |
| `spin_operations.txt` | Full log of all spin symmetry operations |
| `flip_spin_operations.txt` | Rotation matrices of spin-flip operations (used internally) |
| `BZ_ibz_<type>.png` | Plot of the irreducible Brillouin zone with the general k-point |
| `BZ_mapped_<type>.png` | Plot of the IBZ mapped by all symmetry operations |

---

## Example usage
```
$ python generate_kpath.py
=== Altermagnetic K-Path Generator ===
Recommend to use continues high symmetry kpath as input like G-M-K-G rather than L-M|H-K, otherwise there will be duplicated paths.

>>> Step 0: Compute spin-flip symmetry operations
Enter structure file name (default: POSCAR): 
Enter magnetic moments (space-separated, e.g., '1 -1'):
Moments: 1 -1
========================================
1. Structure Loading
========================================
Successfully loaded 'POSCAR' containing 6 atoms.

========================================
2. Non-Magnetic Space Group Analysis
========================================
Space Group: P6_3mc (186)

========================================
3. Magnetic Configuration Input
========================================
Enter magnetic moments (space-separated, e.g., '1 -1'):
Moments: 1 -1
Using magnetic moments:
[[ 0.  0.  1.]
 [ 0.  0. -1.]
 [ 0.  0.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  0.]]

========================================
4. Spin Space Group Analysis
========================================
Spin-Only Group Type: COLLINEAR(axis=[0. 0. 1.])
Magnetic Space Group: Not found (spglib too old?)
Total Symmetry Operations: 12

========================================
5. Saving Results
========================================
[INFO] All operations written to 'spin_operations.txt'
[INFO] 6 spin-flipping matrices written to 'flip_spin_operations.txt'

>>> Step 1: Reading KPOINTS file...
Enter KPOINTS file name (default: KPATH.in): 
Successfully read 14 k-points from KPATH.in

>>> Step 2: Enter general k-point coordinates
Computing IBZ centroid from 'POSCAR'...
============================================================
Processing: POSCAR
============================================================

Space Group: 186 (P6_3mc)
Point Group: 6mm
Seekpath Bravais: hP2
Setyawan-Curtarolo type: HEX

High-symmetry k-points (6):
  Γ       : [  0.0000,   0.0000,   0.0000]
  A       : [  0.0000,   0.0000,   0.5000]
  H       : [  0.3333,   0.3333,   0.5000]
  K       : [  0.3333,   0.3333,   0.0000]
  L       : [  0.5000,   0.0000,   0.5000]
  M       : [  0.5000,   0.0000,   0.0000]

Symmetry operations: 12
With time-reversal: 24

==================================================
NUMERICAL VOLUME CENTROID
==================================================
Cartesian:  [0.394158, 0.409621, 0.211766]
Fractional: [0.277778, 0.111111, 0.250000]
IBZ Volume: 8.205802e-02

==================================================
SYMBOLIC VOLUME CENTROID
==================================================
  k1 = 5/18
  k2 = 1/9
  k3 = 1/4

Verification:
  Symbolic:   [0.277778, 0.111111, 0.250000]
  Numerical:  [0.277778, 0.111111, 0.250000]
==================================================

Captured view angles:
  Fig1 (IBZ):    elev=25.0, azim=-55.0
  Fig2 (Mapped): elev=25.0, azim=-55.0
Saved: ./BZ_ibz_HEX.png
Saved: ./BZ_mapped_HEX.png
General k-point (IBZ centroid): [0.277778, 0.111111, 0.250000]

>>> Step 3: Selecting Transformation Matrix R
Found 6 pre-calculated spin-flip operations:

  Option 1:
    [  1 -1  0 ]
    [  1  0  0 ]
    [  0  0  1 ]

  Option 2:
    [ -1  0  0 ]
    [  0 -1  0 ]
    [  0  0  1 ]

  Option 3:
    [  0  1  0 ]
    [ -1  1  0 ]
    [  0  0  1 ]

  Option 4:
    [  0  1  0 ]
    [  1  0  0 ]
    [  0  0  1 ]

  Option 5:
    [  1 -1  0 ]
    [  0 -1  0 ]
    [  0  0  1 ]

  Option 6:
    [ -1  0  0 ]
    [ -1  1  0 ]
    [  0  0  1 ]

Select an operation number (1-6)
Press [Enter] for default (1), or type number: 
Selected default: Option 1

>>> Step 4: Processing k-points...
Using Transformation Matrix R:
  [ 1. -1.  0.]
  [1. 0. 0.]
  [0. 0. 1.]
k' (transformed k): [-0.1111, 0.3889, 0.2500]
Found 8 unique high-symmetry points
  0: GAMMA = (0.0000, 0.0000, 0.0000)
  1: M = (0.5000, 0.0000, 0.0000)
  2: K = (0.3333, 0.3333, 0.0000)
  3: GAMMA = (0.0000, 0.0000, 0.0000)
  4: A = (0.0000, 0.0000, 0.5000)
  5: L = (0.5000, 0.0000, 0.5000)
  6: H = (0.3333, 0.3333, 0.5000)
  7: A = (0.0000, 0.0000, 0.5000)

Modified k-points (original: 14, new: 28):

>>> Step 5: Save modified file
Enter output filename (default: KPOINTS_modified): 
Modified KPOINTS file written to: KPOINTS_modified

Process completed successfully!
```
---
## Additional Utilities

### Standalone spin-flip analysis
```bash
python3 find_sf_operations.py
```
Runs Step 0 independently. Useful if you already have `flip_spin_operations.txt` and only need to re-run the k-path generation.

### IBZ centroid for a single structure
```bash
python3 compute_centroid_hybrid.py POSCAR
```
Computes the IBZ centroid for any structure file and produces 3D plots. Supports all 14 Bravais lattice types.

### Batch processing (multiple structures)
```bash
python3 batch_centroid_hybrid.py /path/to/structures/ --output summary.csv
python3 batch_centroid_hybrid.py file1.vasp file2.cif file3.vasp
```
Processes a directory of structure files and writes a CSV summary of all centroids.

### Monoclinic IRBZ vertex finder
```bash
python3 find_irbz_vertices.py
```
Advanced diagnostic for C-centered monoclinic (MCLC1 / C2/m) structures. Finds all IRBZ vertices by computing Voronoi + symmetry-plane intersections, and identifies any unlabeled vertices not covered by the standard Setyawan-Curtarolo tables.

---

## Requirements

- Python ≥ 3.9
- See `requirements.txt` for package versions

```bash
pip install -r requirements.txt
```

---

## Note
Doesn't work for triclinic system at the moment.

## Please cite
Bibtex:
```
@article{v3fg-6smc,
  title = {$G$-type antiferromagnetic ${\mathrm{BiFeO}}_{3}$ is a multiferroic $g$-wave altermagnet},
  author = {Urru, Andrea and Seleznev, Daniel and Teng, Yujia and Park, Se Young and Reyes-Lillo, Sebastian E. and Rabe, Karin M.},
  journal = {Phys. Rev. B},
  volume = {112},
  issue = {10},
  pages = {104411},
  numpages = {14},
  year = {2025},
  month = {Sep},
  publisher = {American Physical Society},
  doi = {10.1103/v3fg-6smc},
  url = {https://link.aps.org/doi/10.1103/v3fg-6smc}
}

```


