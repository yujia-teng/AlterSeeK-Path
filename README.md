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

Run:

```bash
python generate_kpath.py
```

The script guides you through 5 steps interactively. The two most important choices are your **structure file** and your **k-path source**:

---

### Workflow A — Magnetic CIF (`.mcif`)

If you have a magnetic CIF file, magnetic moments are detected automatically:

```
Enter structure file name (default: POSCAR): MnF2.mcif
# → moments extracted automatically, no manual input needed
```

---

### Workflow B — POSCAR with manual moments

If you use a POSCAR, you will be prompted to enter moments:

```
Enter structure file name (default: POSCAR): POSCAR
Enter magnetic moments in atom order from structure file (space-separated).
  Trailing non-magnetic atoms auto-fill to 0 (e.g., just '1 -1' if mag atoms come first).
  If non-magnetic atoms appear first, include 0s (e.g., '0 0 1 -1').
Moments: 1 -1
```

Moments are entered **in atom order** as they appear in the POSCAR. Only type moments for magnetic atoms — trailing non-magnetic atoms are filled with 0 automatically.

---

### K-path: auto-generate or provide file

**Option 1 — Auto-generate (recommended):** Just press Enter at the KPOINTS prompt. The path is generated automatically via seekpath.

**Option 2 — Provide your own file:** Place a line-mode `KPATH.in` in the same directory (e.g. from VASPKIT task 303), then enter the filename when prompted.

---

### Step summary

| Step | What it does | Input needed |
|------|-------------|--------------|
| **0** | Finds spin-flip symmetry operations from structure | Structure file name; magnetic moments (manual for POSCAR, auto for mcif) |
| **1** | Reads the high-symmetry k-path | KPOINTS file name (or Enter for auto) |
| **2** | Auto-computes the general k-point (IBZ centroid) | *(automatic)* |
| **3** | Selects the spin-flip transformation matrix | Choose from list, or press Enter for default — any choice gives an equivalent band structure, since k' always lands in the opposite-spin IBZ |
| **4** | Generates the enriched k-path | *(automatic)* |
| **5** | Saves the output file | Output file name |

**Output:** `KPOINTS_modified` — ready to use directly in a VASP band structure calculation.

---

## Input Files

### Structure file
- **POSCAR** — standard VASP structure file; magnetic moments entered manually at prompt
- **`.mcif` file** — magnetic CIF; moments extracted automatically
- **`.cif` file** — standard CIF; magnetic moments entered manually at prompt

### `KPATH.in` (optional)
Line-mode KPOINTS file with the high-symmetry path. If not provided, the path is auto-generated via seekpath. To generate manually: [VASPKIT](https://vaspkit.com/) (task 303) or [seekpath](https://www.materialscloud.org/work/tools/seekpath).

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

## Example

```
$ python generate_kpath.py
=== Altermagnetic K-Path Generator ===

>>> Step 0: Compute spin-flip symmetry operations
Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): 
Enter magnetic moments in atom order from structure file (space-separated).
  Trailing non-magnetic atoms auto-fill to 0 (e.g., just '1 -1' if mag atoms come first).
  If non-magnetic atoms appear first, include 0s (e.g., '0 0 1 -1').
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
3. Magnetic Configuration
========================================        
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
Total Symmetry Operations: 12

========================================
5. Saving Results
========================================
[INFO] All operations written to 'spin_operations.txt'
[INFO] 6 spin-flipping matrices written to 'flip_spin_operations.txt'

>>> Step 1: High-symmetry k-path
Computing IBZ centroid and k-path from 'POSCAR'...
============================================================
Processing: POSCAR
============================================================

Space Group: 186 (P6_3mc)
Point Group: 6mm
Seekpath Bravais: hP2

Symmetry operations: 12
With time-reversal: 24
Saved: .\POSCAR_ibz_HEX.png
Saved: .\POSCAR_mapped_HEX.png
Auto-generated path: GAMMA-M-K-GAMMA-A-L-H-A-|L-M-|H-K
Press [Enter] to use this path, or type a filename to load your own:
Using auto-generated path (9 segments, 18 k-points)

>>> Step 2: Enter general k-point coordinates
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

  Note: all options produce equivalent band structures —
  k' = R⁻ᵀk always lands in the opposite-spin IBZ for any spin-flip R.

Select an operation number (1-6)
Press [Enter] for default (1), or type number:
Selected default: Option 1

>>> Step 4: Processing k-points...
Using Transformation Matrix R:
  [ 1. -1.  0.]
  [1. 0. 0.]
  [0. 0. 1.]
k' (transformed k): [-0.1111, 0.3889, 0.2500]
High-symmetry path: GAMMA-M-K-GAMMA-A-L-H-A|L-M|H-K
Found 3 path segment(s):
  Part 1: GAMMA - M - K - GAMMA - A - L - H - A
  Part 2: L - M
  Part 3: H - K
Generated path: GAMMA-M-k|k'-M'-K'-k'|k-K-GAMMA-k|k'-GAMMA-A'-k'|k-A-L-k|k'-L'-H'-k'|k-H-A|L-M|H-K

>>> Step 5: Save modified file
Enter output filename (default: KPOINTS_modified):
Modified KPOINTS file written to: KPOINTS_modified

Process completed successfully!
```

---

## Additional Utilities

### Standalone spin-flip analysis
```bash
python find_sf_operations.py
```
Runs Step 0 independently. Useful if you already have `flip_spin_operations.txt` and only need to re-run the k-path generation.

### IBZ centroid for a single structure
```bash
python compute_centroid_hybrid.py POSCAR
```
Computes the IBZ centroid for any structure file and produces 3D plots. Supports all 14 Bravais lattice types.

### Batch processing (multiple structures)
```bash
python batch_centroid_hybrid.py /path/to/structures/ --output summary.csv
python batch_centroid_hybrid.py file1.vasp file2.cif file3.vasp
```
Processes a directory of structure files and writes a CSV summary of all centroids.

### Monoclinic IRBZ vertex finder
```bash
python find_irbz_vertices.py
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
Triclinic systems (space groups 1–2) have no altermagnetic splitting. The tool writes a plain seekpath k-path for these cases and exits.

---

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
