# Plotting

AlterSeeK-Path includes a band plotter for spin-resolved VASP band structures.

After the VASP band calculation, run VASPKIT task `303` to generate the
reformatted band files, then run:

```bash
alterseek-path bandplot
```

The optional shortcut command is equivalent:

```bash
alterseek-bandplot
```

## Default Inputs And Output

By default, the plotter reads:

```text
KLABELS
REFORMATTED_BAND_UP.dat
REFORMATTED_BAND_DW.dat
```

and writes:

```text
alterband.png
```

To write PDF output:

```bash
alterseek-path bandplot -o alterband.pdf
```

## Plot Settings With `alterband.toml`

If a file named `alterband.toml` exists in the same directory, the band plotter
uses it automatically. The main `alterseek-path` workflow writes this file after
KPOINTS generation and records the detected lattice type when available.

A typical configuration is:

```toml
emin = -2
emax = 2
fig_width = 16
fig_height = 5
gap_width_inches = 0.05
lattice_type = "tI1"
split_panels = 0
rotate_xtick_labels = false
xtick_rotation = 45
output = "alterband.png"
```

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

## Special HPKOT Shading

When `lattice_type` is present, special HPKOT path intervals are shaded light
grey in the band plot. The main workflow writes this value automatically after
KPOINTS generation. For direct plotting, set it manually, for example:

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

Missing or `0` keeps a single panel. The split affects only rendering; it does
not change KPOINTS or band data.

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
