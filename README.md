# AlterSeeK-Path

AlterSeeK-Path generates k-point paths for altermagnet band-structure calculations. It inserts a general k point `k` and its spin-flip partner `k'` into a standard high-symmetry path, using the IBZ centroid as the default general point.

![AlterSeeK-Path example](./example/VASP/HEX.png)

---

## Installation

Requires Python >= 3.11.

```bash
git clone https://github.com/yujia-teng/AlterSeeK-Path.git
cd AlterSeeK-Path
pip install -r requirements.txt
pip install -e .
```

---

## Documentation

For a longer user guide, see
[yujia-teng.github.io/AlterSeeK-Path](https://yujia-teng.github.io/AlterSeeK-Path/).

---

## Quick Start

```bash
alterseek-path
```

This runs interactively, prompting for the structure file, spin axis,
moments, k-path, and so on. The default output file is:

```text
KPOINTS_alter
```

To skip the prompts (e.g. for repeated runs on the same structure), put an
`alterseek_input.toml` file in the working directory:

```toml
structure = "POSCAR"
spin_axis = "0 0 1"
moments = "5 -5"
path = ""
flip_option = 1
output_code = "vasp"
```

`alterseek-path` reads any keys it finds and only prompts for the ones that
are missing. See [Workflow](https://yujia-teng.github.io/AlterSeeK-Path/workflow/#skipping-repeated-prompts-with-alterseek_inputtoml)
for the full field reference.

---

## Inputs

- Structure file: `POSCAR` / `.vasp` / `.cif` (moments entered manually) or
  `.mcif` (moments read from the file). MCIF input must be collinear: all
  nonzero vector moments must be parallel or antiparallel within an absolute
  transverse-moment tolerance of `0.02` in the MCIF moment units (normally
  Bohr magnetons); noncollinear MCIFs are rejected.
- Spin axis + moments: Cartesian axis (default `0 0 1`) and VASP
  `MAGMOM`-style scalar moments (e.g. `5*0 2*1.0`); missing trailing values
  default to `0`, while excess values are rejected.
- K-path: press Enter for the auto-generated path, or give a line-mode
  `KPATH.in`/KPOINTS-style file. A custom path is interpreted in the reciprocal
  fractional basis of the submitted structure and must contain labeled,
  complete endpoint pairs.

These (plus the Step 3 operation choice and Step 5 output code) can also be
supplied via `alterseek_input.toml` — see Quick Start above.

---

## Example Run

```text
=== Altermagnetic K-Path Generator ===

>>> Step 0: Spin symmetry
Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): POSCAR
Spin axis in Cartesian coordinates (default: 0 0 1): 0 0 1
Magnetic moments along this axis (atom order, trailing atoms auto-fill to 0): 8 -8 4*0

Input structure: POSCAR, 6 atoms
Nonmagnetic primitive cell:   SG P6_3mc (186)  PG 6mm  Laue 6/mmm  [6 atoms, hP2]
Magnetic primitive cell:      SG P6_3mc (186)  PG 6mm  Laue 6/mmm  [6 atoms, hP2]
Submitted analysis cell:      SG P6_3mc (186)  PG 6mm  Laue 6/mmm  [6 atoms, hP2]
Phase: AFM(Altermagnet)
Oriented SSG: 186.156.1.1.L
SSG Symbol (Chen-Liu): P -1|6_{3} 1|m -1|c infinity_{001}m|1
MSG without SOC: P6_3'mc' (BNS 186.206), Type III
Spin operations: 6 flip, 6 preserve

>>> Step 1: High-symmetry k-path
Path: GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K
Press [Enter] to use this path, or type a filename to load your own:
Using HPKOT hP2 path (9 segments, 18 k-points)

>>> Step 2: General k-point
IBZ centroid (SeeK-path fractions reused in submitted-cell basis): [0.277778, 0.111111, 0.250000]
IBZ centroid (KPOINTS output basis): [0.277778, 0.111111, 0.250000]

>>> Step 3: Spin-flip operation
Found 12 spin-flip operations R.
  Note: R is in the submitted structure 'POSCAR' fractional basis;
  rotation axis/mirror plane indices are in the reciprocal (b1,b2,b3) basis.
Default R: Option 1
Press [Enter] to use default, type a number, or 'list' to show matrices: 1
Selected: Option 1  (C6+ [0 0 1])

>>> Step 4: Build altermagnetic path
k' = [-0.1111, 0.3889, 0.2500]
Generated path: GAMMA-M-k | k'-M'-K'-k' | k-K-GAMMA-k | ... | k-H-A | L-M | H-K
Full path: 9 original segments -> 21 generated segments, 36 k-points

>>> Step 5: Save
Output code ([vasp]/qe): vasp
Modified KPOINTS file written to: KPOINTS_alter
Band plot config updated: alterband.toml (lattice_type = "hP2")

Done.
Displaying generated figure(s)...
Saved: alterseek_output\POSCAR_ibz_hP2.png
Saved: alterseek_output\POSCAR_spinflip_hP2.png
Saved: alterseek_output\POSCAR_spinbz_hP2.png
Saved: alterseek_output\POSCAR_spinbz_top_hP2.png
```

---

## 2D / Slab Mode

For 2D materials computed as slabs (vacuum along one lattice vector), run:

```bash
alterseek-path --2d
```

2D mode restricts the k-path and IBZ centroid to the physical in-plane
(vacuum k = 0) reciprocal plane, and reports whether any spin-flip operation
produces in-plane spin splitting. See [Workflow](https://yujia-teng.github.io/AlterSeeK-Path/workflow/#2d-slab-mode)
for the vacuum-axis detection details and output figures.

---

## Cell Setting and Brillouin Zone

The submitted structure remains the calculation and output cell. Its complete
space group determines the translation lattice used for the BZ, IBZ, centroid,
SeeK-path labels, and figures. A true primitive-lattice supercell therefore
receives its own BZ and path, while a centered conventional setting is reduced
to its crystallographic primitive lattice by SeeK-path. A determinant-one
basis change is kept exactly as submitted. AlterSeeK-Path never replaces the
calculation POSCAR or emits a replacement MAGMOM handoff.

For magnetic input, FindSpinGroup still supplies G0, its setting transforms,
and the spin operations. AlterSeeK-Path obtains the complete standard Hall
operation set of G0, transforms every Seitz operation `(R|t)` into the
submitted basis, and generates validated generic marker orbits with `Rr+t`
alongside the submitted real atoms in an in-memory SeeK-path cell. The
magnetic-primitive representatives select the compatible Hall setting; they
are not treated as the complete conventional-cell operation set. This retains
nonsymmorphic screw/glide shifts and genuine Bravais-centering copies, while
excluding extra smaller-parent translations that are not members of G0. The
no-moments route uses the same complete-operation construction from spglib's
detected Hall setting. The helper is accepted only when spglib recovers exactly
that operation set and the expected space group. `volume_original_wrt_prim`
therefore equals two for C-centered GdAuGe 211 and three for
conventional-hexagonal R3c BiFeO3. SeeK-path constructs the BZ, IBZ, centroid,
and HPKOT path in its standardized primitive basis. The workflow then reuses
those fractional coordinates unchanged in the submitted structure's reciprocal
basis; it does not transform them to preserve the standardized Cartesian
k-points. Spin figures use the spin operations transformed into that
standardized basis, while the written KPOINTS spin partners use the original
submitted-basis operations after the fractional reuse. This keeps the mapped
IBZ inside the plotted first BZ and keeps `k` and `k'` symmetry-related in the
calculation cell. The FindSpinGroup magnetic primitive cell is retained
separately as the `*_magnetic_primitive.mcif` diagnostic.

---

## Band Plotting

After the VASP band calculation and VASPKIT task `211`, run:

```bash
alterseek-path bandplot
```

This reads `KLABELS`/`REFORMATTED_BAND_UP.dat`/`REFORMATTED_BAND_DW.dat` and
writes `alterband.png`, using settings from `alterband.toml` (written
automatically by the main workflow) if present. See
[Plotting](https://yujia-teng.github.io/AlterSeeK-Path/plotting/) for the
full settings reference (energy window, split panels, gap width, etc.).

---

## Citation

```bibtex
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
