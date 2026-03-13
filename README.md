# DeepseeK-Path

A tool for generating enriched k-point paths for band structure calculations of altermagnetic and magnetically ordered materials.

It transforms a standard high-symmetry path like `Γ−M−K−Γ` into a symmetry-enriched path such as:

```
Γ−M | k−k'−M'−K'−k' | k−K−Γ | k−k'−Γ | ...
```

where `k` is a general k-point (the volume centroid of the irreducible Brillouin zone) and `k_t` is its spin-flip-transformed partner, revealing split bands that are degenerate along the standard path.

**Current support:** VASP (full). Quantum ESPRESSO support planned.

---

## Installation

```bash
git clone https://github.com/yujia-teng/DeepseeK-Path.git
cd DeepseeK-Path
pip install -r requirements.txt
```

---

## Quick Start

Place your `POSCAR` and `KPATH.in` (line-mode KPOINTS file, e.g. from VASPKIT) in the same directory as the scripts, then run:

```bash
python3 generate_kpath.py
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

### `POSCAR`
Standard VASP structure file for the magnetic material.

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

## Example

The `example/` directory contains a complete worked example for GdAuGe (hexagonal, altermagnetic):

```
example/
├── POSCAR              # GdAuGe structure (P6_3mc, space group 186)
├── KPATH.in            # Standard Γ-M-K-Γ-A-L-H-A path
├── KPOINTS_modified    # Enriched path (output of generate_kpath.py)
├── band_structure.png  # Band structure showing altermagnetic splitting
└── plot_altermagnet_band.py   # Script to plot the band structure
```

Run the example:
```bash
cd example
python3 ../generate_kpath.py
# Enter POSCAR, moments: 1 -1, KPATH.in, then defaults for the rest
```

---

## How It Works

1. **Spin-flip operations** (`find_sf_operations.py`): Uses [spinspg](https://github.com/spglib/spinspg) to find all spin symmetry operations of the magnetic structure, then filters for those with determinant −1 on the spin rotation (spin-flip operations).

2. **General k-point** (`compute_centroid_hybrid.py`): Determines the Bravais lattice type via [seekpath](https://seekpath.readthedocs.io/), retrieves the standard IRBZ vertices from the Setyawan-Curtarolo tables ([Comp. Mat. Sci. 49, 299, 2010](https://doi.org/10.1016/j.commatsci.2010.05.010)), and computes the volume centroid of the IRBZ. This point is in a general position — not on any symmetry element — ensuring the spin-flip partner is a distinct point.

3. **Path enrichment** (`generate_kpath.py`): For each segment X→Y in the original path, inserts the segments `X→Y | k→k_t−Y_t−Z_t−k_t | k→Z` where `k_t = R⁻ᵀ k` under the spin-flip rotation R.

---

## Requirements

- Python ≥ 3.9
- See `requirements.txt` for package versions

```bash
pip install -r requirements.txt
```

---

## Citation

If you use this tool in your research, please cite:

> Teng, Y. *DeepseeK-Path: A tool for generating enriched k-paths for altermagnetic band structure calculations.* Rutgers University (2025). https://github.com/yujia-teng/DeepseeK-Path

---

## Notes

- Spin alignment along the z-axis is assumed by default
- Spin-group analysis does not include spin-orbit coupling
- Output uses `k_t` notation instead of `k'` due to plain-text formatting constraints
- For monoclinic structures, `find_irbz_vertices.py` provides additional IRBZ analysis
