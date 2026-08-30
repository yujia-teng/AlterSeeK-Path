# Plotting

AlterSeeK-Path includes spin-resolved band plotters for VASP, Quantum
ESPRESSO, and ABINIT.

## Basic Workflow

1. Generate the k-path and plotting configuration with `alterseek-path`.
2. Run the band calculation and the code-specific post-processing described
   below.
3. Run:

```bash
alterseek-plot
```

`alterseek-plot` detects the code from the generated
`alterseek_plot_<code>.toml` file. To choose the code explicitly, run:

```bash
alterseek-plot vasp
alterseek-plot qe
alterseek-plot abinit
```

## VASP Band Plotting

Run VASPKIT task `211` after the VASP calculation to generate the reformatted
band files, then use `alterseek-plot` or `alterseek-plot vasp`.

### Required Files

By default, the plotter reads the standard VASPKIT output filenames:

| File | Source | Purpose |
|------|--------|---------|
| `KLABELS` | VASPKIT task `211` | k-point tick labels |
| `REFORMATTED_BAND_UP.dat` | VASPKIT task `211` | spin-up band data |
| `REFORMATTED_BAND_DW.dat` | VASPKIT task `211` | spin-down band data |
| `alterseek_plot_vasp.toml` | `alterseek-path` | optional plotting settings |

The default output is:

```text
alterband.png
```

To write PDF output:

```bash
alterseek-plot vasp -o alterband.pdf
```

### Plot Settings

`alterseek-path` creates `alterseek_plot_vasp.toml` and sets or updates
`lattice_type`. The remaining settings are included as commented defaults;
uncomment and change them as needed. Only `--config` and `-o` are available as
command-line overrides.

```toml
lattice_type = "hP2"  # written by alterseek-path
klabels = "KLABELS"
up = "REFORMATTED_BAND_UP.dat"
down = "REFORMATTED_BAND_DW.dat"
emin = -2
emax = 2
fig_width = 12
fig_height = 5
gap_width_inches = 0.05
split_panels = 1
output = "alterband.png"
save_pdf = false
```

| Setting | Meaning |
|---------|---------|
| `klabels` | VASPKIT label file |
| `up`, `down` | Spin-up/down VASPKIT band files |
| `lattice_type` | Enables lattice-aware label handling (e.g. repairing VASPKIT-truncated `oI2`/`oI3` labels) |
| `emin`, `emax` | Energy window in eV |
| `fig_width`, `fig_height` | Figure size in inches |
| `gap_width_inches` | Visual width of each `k\|k'` separator gap, kept consistent across path lengths |
| `split_panels` | `1` for one panel, `2`/`3` for stacked panels (rendering only, doesn't change KPOINTS/band data) |
| `output` | Output image filename (`.png` or `.pdf`) |
| `save_pdf` | Also save a `.pdf` copy alongside `output`, if `output` isn't already a `.pdf` |

## Quantum ESPRESSO Band Plotting

The QE plotter reads user-named spin-up and spin-down `bands.x` `.gnu` files
and the generated `KPOINTS_alter_qe` k-path file:

```bash
alterseek-plot qe
```

`alterseek-path` creates `alterseek_plot_qe.toml` with commented defaults.
Uncomment the QE filenames and other settings as needed; only `--config` and
`-o` are available as command-line overrides.

```toml
band_up = "band_up.gnu"  # replace with your spin-up bands.x .gnu file
band_down = "band_down.gnu"  # replace with your spin-down bands.x .gnu file
kpoints_qe = "KPOINTS_alter_qe"
fermi_ev = 0.0  # replace with the Fermi energy from the QE output file
emin = -2
emax = 2
fig_width = 12
fig_height = 5
gap_width_inches = 0.05
split_panels = 1
output = "alterband_qe.png"
save_pdf = false
```

| Setting | Meaning |
|---------|---------|
| `band_up`, `band_down` | Spin-up/down `bands.x` `.gnu` files |
| `kpoints_qe` | The `K_POINTS crystal_b` k-path file written by `alterseek-path` |
| `fermi_ev` | Fermi energy in eV subtracted from the `.gnu` energies |
| `emin`, `emax` | Energy window in eV |
| `fig_width`, `fig_height` | Figure size in inches |
| `gap_width_inches` | Visual width of each `k\|k'` separator gap, kept consistent across path lengths |
| `split_panels` | `1` for one panel, `2`/`3` for stacked panels |
| `output` | Output image filename (`.png` or `.pdf`) |
| `save_pdf` | Also save a `.pdf` copy alongside `output`, if `output` isn't already a `.pdf` |

## ABINIT Band Plotting

The ABINIT plotter reads the plain-text `_EIG` output and generated
`KPOINTS_alter_abinit` k-path file:

```bash
alterseek-plot abinit
```

The plotter converts `_EIG` energies from Hartree to eV, obtains the Fermi
level from `abo`, and computes path distance from the submitted structure's
reciprocal lattice. Set `fermi_ev` only when the `.abo` file is unavailable.

`alterseek-path` creates `alterseek_plot_abinit.toml` and records the submitted
`structure`. The remaining settings are included as commented defaults;
uncomment and change them as needed. Only `--config` and `-o` are available as
command-line overrides.

```toml
structure = "POSCAR"  # written by alterseek-path
eig = "EIG"  # replace with your ABINIT _EIG file
kpoints_abinit = "KPOINTS_alter_abinit"
abo = "abo"  # replace with your ABINIT .abo output file
emin = -2
emax = 2
fig_width = 12
fig_height = 5
gap_width_inches = 0.05
split_panels = 1
output = "alterband_abinit.png"
save_pdf = false
```

| Setting | Meaning |
|---------|---------|
| `eig` | ABINIT `_EIG` file (plain text, Hartree, both spin channels) |
| `kpoints_abinit` | The `kptopt`/`ndivk`/`kptbounds` k-path file written by `alterseek-path` |
| `structure` | Submitted POSCAR/`.vasp`, `.cif`, or `.mcif` structure used to build the reciprocal lattice for the k-path distance |
| `abo` | ABINIT output used to read the Fermi level |
| `fermi_ev` | Manual Fermi-level override in eV; takes priority over `abo` |
| `emin`, `emax` | Energy window in eV |
| `fig_width`, `fig_height` | Figure size in inches |
| `gap_width_inches` | Visual width of each `k\|k'` separator gap, kept consistent across path lengths |
| `split_panels` | `1` for one panel, `2`/`3` for stacked panels |
| `output` | Output image filename (`.png` or `.pdf`) |
| `save_pdf` | Also save a `.pdf` copy alongside `output`, if `output` isn't already a `.pdf` |
