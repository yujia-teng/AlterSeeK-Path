# Plotting

AlterSeeK-Path includes a band plotter for spin-resolved VASP band structures.

## Basic Workflow

1. Generate `KPOINTS_alter` with `alterseek-path`.
2. Run the VASP band calculation.
3. Run VASPKIT task `211` to generate the reformatted band files.
4. Run:

```bash
alterseek-path bandplot
```

The optional shortcut command is equivalent:

```bash
alterseek-bandplot
```

## Required Files

By default, the plotter reads the standard VASPKIT output filenames:

| File | Source | Purpose |
|------|--------|---------|
| `KLABELS` | VASPKIT task `211` | k-point tick labels |
| `REFORMATTED_BAND_UP.dat` | VASPKIT task `211` | spin-up band data |
| `REFORMATTED_BAND_DW.dat` | VASPKIT task `211` | spin-down band data |
| `alterband.toml` | `alterseek-path` | optional plotting settings |

The default output is:

```text
alterband.png
```

To write PDF output:

```bash
alterseek-path bandplot -o alterband.pdf
```

## Plot Settings

If `alterband.toml` exists in the same directory, the band plotter uses it
automatically. The main `alterseek-path` workflow writes this file after
KPOINTS generation, recording the detected lattice type. A typical
configuration:

```toml
lattice_type = "hP2"
emin = -2
emax = 2
fig_width = 12
fig_height = 5
gap_width_inches = 0.05
split_panels = 1
output = "alterband.png"
save_pdf = false
```

| Setting | Meaning | Equivalent CLI flag |
|---------|---------|----------------------|
| `lattice_type` | Enables lattice-aware label handling (e.g. repairing VASPKIT-truncated `oI2`/`oI3` labels) | `--lattice-type` |
| `emin`, `emax` | Energy window in eV | `--emin`, `--emax` |
| `fig_width`, `fig_height` | Figure size in inches | -- |
| `gap_width_inches` | Visual width of each `k\|k'` separator gap, kept consistent across path lengths | `--gap-width-inches` |
| `split_panels` | `1` for one panel, `2`/`3` for stacked panels (rendering only, doesn't change KPOINTS/band data) | `--split-panels` |
| `output` | Output image filename (`.png` or `.pdf`) | `-o` |
| `save_pdf` | Also save a `.pdf` copy alongside `output`, if `output` isn't already a `.pdf` | `--save-pdf` |

The plotter repairs known VASPKIT boundary-label and truncation errors in
memory and reports any correction it uses. It does not rewrite `KLABELS`.
Both VASP and QE plotting configurations are validated; malformed TOML,
unknown settings, and invalid value types stop with an error.

Command-line flags override the TOML file, e.g.:

```bash
alterseek-path bandplot --emin -3 --emax 3 --lattice-type mC2 --split-panels 2 -o alterband.pdf
```

Use explicit filenames if your VASPKIT outputs have different names:

```bash
alterseek-path bandplot --klabels KLABELS --up REFORMATTED_BAND_UP.dat --down REFORMATTED_BAND_DW.dat
```

## Quantum ESPRESSO Band Plotting

For QE workflows, use the separate `plot_alterband_qe.py` script (or
`alterseek-bandplot-qe`) instead. It reads `bands.x` `.gnu` output and the
`KPOINTS_alter_qe` waypoint file written by `alterseek-path`:

```bash
alterseek-bandplot-qe
```

Settings are config-file only (there is no CLI-flag equivalent besides
`--config`/`-o`). If `alterband_qe.toml` exists in the working directory, it
is used automatically:

```toml
band_up = "band_up.gnu"
band_down = "band_down.gnu"
fermi_ev = 0.0
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
| `kpoints_qe` | The `K_POINTS crystal_b` waypoint file written by `alterseek-path` |
| `fermi_ev` | Energy shift applied before plotting (QE `.gnu` output is not pre-shifted to E_F, unlike VASPKIT's reformatted files) |
| `emin`, `emax` | Energy window in eV |
| `fig_width`, `fig_height` | Figure size in inches |
| `gap_width_inches` | Same meaning as the VASP plotter's setting |
| `split_panels` | `1` for one panel, `2`/`3` for stacked panels |
| `output` | Output image filename (`.png` or `.pdf`) |
| `save_pdf` | Also save a `.pdf` copy alongside `output`, if `output` isn't already a `.pdf` |

There is no `lattice_type` setting for QE plotting — it exists on the VASP
side only to repair VASPKIT's truncated labels, which doesn't apply here
since QE labels come from AlterSeeK-Path's own `KPOINTS_alter_qe` writer.
