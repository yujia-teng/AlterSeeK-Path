# Plotting

AlterSeeK-Path includes a band plotter for spin-resolved VASP band structures.

## Basic Workflow

1. Generate `KPOINTS_modified` with `alterseek-path`.
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

If a file named `alterband.toml` exists in the same directory, the band plotter
uses it automatically. The main `alterseek-path` workflow writes this file after
KPOINTS generation and records the detected lattice type when available.

A typical configuration is:

```toml
lattice_type = "hP2"
emin = -2
emax = 2
fig_width = 12
fig_height = 5
gap_width_inches = 0.05
split_panels = 0
output = "alterband.png"
```

| Setting | Meaning |
|---------|---------|
| `lattice_type` | SeeK-path lattice key used for special interval shading |
| `emin`, `emax` | energy window in eV |
| `fig_width`, `fig_height` | figure size in inches |
| `gap_width_inches` | visual width of each `k|k'` separator gap |
| `split_panels` | `0` for one panel, `2` or `3` for stacked panels |
| `output` | output image filename, usually `.png` or `.pdf` |

Then run:

```bash
alterseek-path bandplot
```

Command-line options override the TOML file:

```bash
alterseek-path bandplot -o alterband.pdf
```

## Energy Window And Files

Set the energy window from the command line:

```bash
alterseek-path bandplot --emin -3 --emax 3
```

Use explicit input filenames if your VASPKIT outputs have different names:

```bash
alterseek-path bandplot --klabels KLABELS --up REFORMATTED_BAND_UP.dat --down REFORMATTED_BAND_DW.dat
```

## Special Interval Shading

When `lattice_type` is present, special path intervals are shaded light grey in
the band plot. The main workflow writes this value automatically after KPOINTS
generation. For direct plotting, set it manually, for example:

```toml
lattice_type = "oF3"
```

or pass it on the command line:

```bash
alterseek-path bandplot --lattice-type mC2
```

## Split Panels

Use split panels for long paths that are difficult to read in one row:

```toml
split_panels = 2
```

or:

```toml
split_panels = 3
```

Use `split_panels = 0` for one panel, `2` for two stacked panels, or `3` for
three stacked panels. The split affects only rendering; it does not change
KPOINTS or band data.

The same option can be passed from the command line:

```bash
alterseek-path bandplot --split-panels 2
```

## Gap Width

`gap_width_inches` sets the full visual width of every `k|k'` separator gap.
This keeps the printed separator size consistent across figures with different
path lengths.

```bash
alterseek-path bandplot --gap-width-inches 0.04
```

---

[Back to user guide](index.md)
