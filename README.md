# SPAE
Stellar Parameters, Abundances, and Errors

SPAE is a Python code that uses a state-of-the-art Bayesian method to self-consistently propagate uncertainties from the stellar atmosphere solutions in calculating individual abundances. It uses a pure-Python reimplementation of the LTE plane-parallel spectral analysis code MOOG (Sneden 1973; https://www.as.utexas.edu/~chris/moog.html) and Kurucz model atmosphere grids (http://kurucz.harvard.edu/grids.html) to derive the abundances of elements included in a user-provided linelist.


## Cloning and Installing SPAE

Clone SPAE from Github to a directory on your local machine:

```
git clone https://github.com/simon-schuler/SPAE
```

Install SPAE into a Python environment:

```
pip install ./
```


## Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| `numpy` | Array math throughout | `pip install numpy` |
| `numba` | JIT-compiled physics kernels (`voigt`, `expn2`) in `SPAE/moog` | `pip install numba` |
| `scipy` | Interpolation, optimization, integration (`SPAE/atmos.py`, `SPAE/xspect_ew`) | `pip install scipy` |
| `astropy` | FITS I/O (`SPAE/xspect_ew`, model/linelist parsing) | `pip install astropy` |
| `george` | Gaussian Process continuum fitting in `SPAE/xspect_ew` | `pip install george` |
| `matplotlib` | Plotting and interactive widget | `pip install matplotlib` |
| `emcee` | Bayesian MCMC parameter estimation | `pip install emcee` |
| `tqdm` | Per-step progress bar during `run_spae()` MCMC runs | `pip install tqdm` |
| `ipympl` | Interactive `%matplotlib widget` support in Jupyter | `pip install ipympl` |

All dependencies except `ipympl` and `tqdm` are required for core SPAE functionality. `ipympl` is only needed when using the interactive synthesis widget inside a Jupyter notebook; `tqdm` is optional but strongly recommended — without it, `run_spae()` runs silently with no progress feedback until completion.

## XSpect-EW (equivalent-width measurement)

`SPAE.xspect_ew` measures equivalent widths from a stellar spectrum and writes a MOOG-format linelist-with-EW file directly usable by `SPAE.moog.abfind`. Originally developed by George Vejar ([github.com/forgeousgeorge/XSpect](https://github.com/forgeousgeorge/XSpect)), incorporated into SPAE 2026-09-14.

`Spectrum_Data(filename)` auto-detects the spectrum's format — Keck/MAKEE, GRACES/OPERA, MAROON-X, or a FITS binary table with named wave/flux columns (including this project's own KOA HIRES output, see `Verification/Brewer2016_HIRES`) — no need to say which instrument it came from. For anything not yet recognized, pass already-extracted arrays directly (`KECK_file=False, spectx=..., specty=...`) or a `custom_reader=` callable; see `SPAE/xspect_ew/readers.py` for the detection logic and how to add a new format.

```python
from SPAE.xspect_ew import Spectrum_Data

spec = Spectrum_Data('star_blue.fits')                    # format auto-detected
spec.normalize_all()                                       # per-order GP continuum fit (needed before shifting)
spec.apply_rv_shift(verbose=True)                          # recommended: single RV from strong lines, no reference spectrum needed
spec.load_lines('linelist.txt')                           # MOOG-format: wave, species, EP, loggf, damping
spec.measure_line_ew(4779.439)                              # or measure_all_ew() / measure_ew()
spec.make_ew_doc('linelist_with_ew.txt')                   # MOOG-format output
```

`apply_rv_shift()` is the recommended way to align a spectrum before EW measurement — it measures one effective RV from a handful of strong, well-identified lines (`SPAE.xspect_ew.radial_velocity.RV_REFERENCE_LINES`) and applies it as a proper multiplicative (1+v/c) shift to every order. The older `estimate_shift()`/`clean_shift()` (per-order cross-correlation against a reference spectrum, e.g. a solar atlas) remains available and is still the right tool when you specifically need reference-spectrum registration (e.g. `combine_spectra()`), but its extrapolation for orders with no reference-spectrum overlap was measured to be substantially less reliable (RMS ~67 mA vs ~37 mA position error on a real test, worst case 119 mA) — unreliable enough to risk misidentifying a line during EW measurement. `clean_shift()` now warns when it's extrapolating outside its actual reference coverage.

`measure_ew()`/`measure_all_ew()` fit each line twice: a **global-continuum fit** (assumes the existing per-order normalization already puts the continuum at 1.0 — this is what gets reported in `lines_ew`/the output linelist) and, when `fit_continuum=True`, an additional **local-continuum diagnostic fit** (estimates a flat local continuum from the line's own wing data and refits against that — kept only as a comparison value in `lines_ew_local`, never the reported EW). Passing a high-S/N `load_reference_atlas()` spectrum (e.g. the Kurucz solar flux atlas, for solar-type targets) lets the fit additionally cross-check a line's wing against that reference to catch shallow blends a simple sigma-clip would miss.

Multiple exposures of the same star/order coverage can be combined before measurement: `combine_spectra()` aligns and co-adds two `Spectrum_Data` instances order-by-order, with cosmic-ray/bad-pixel rejection across the pair.

### Interactive widget

`ew_interactive()` launches a step-by-step widget (Load → Normalize → RV Shift → Measure EW) for visually walking through the whole measurement process instead of calling each stage from a script — inspect/tune the continuum fit per order, review or override the RV shift, and step through lines adjusting `measure_ew()`'s `ex_params` live while watching the fit update. Works standalone or embedded in Jupyter, same as the `pymoog` synthesis widget below.

```python
from SPAE.xspect_ew import ew_interactive

w = ew_interactive('star_blue.fits')                        # single exposure
# or, for multiple exposures of the same star:
w = ew_interactive(spectra=['star_blue_1.fits', 'star_blue_2.fits'])
```

**From a Jupyter notebook**, run `%matplotlib widget` in its own cell first, then call `ew_interactive()` as above — keep the returned widget assigned to a variable (`w = ...`) to prevent it from being garbage-collected and disconnecting its callbacks.

Sample data (`SPAE/xspect_ew/data/`) — a solar HIRES spectrum and Fe linelist — is included for testing.

Minimum Python version: **3.6**


---

## pymoog — standalone spectral analysis tool

`pymoog` is a command-line interface to the pure-Python MOOG engine bundled with SPAE. After installing SPAE it is available as a shell command:

```
pymoog [batch.par]
```

It reads a MOOG-style `batch.par` parameter file (defaulting to `./batch.par`), runs the analysis mode specified on the first line, and writes output files in the same format as MOOGSILENT.

### Supported modes

| Mode | Description |
|------|-------------|
| `abfind` | Equivalent-width abundance analysis |
| `blends` | Blended-line abundance fitting |
| `synth` | Synthetic spectrum calculation |
| `ewfind` | Predicted equivalent widths from known abundances |
| `cog` | Curve-of-growth mapping |
| `weedout` | Line culling by opacity ratio |
| `doflux` | Emergent continuum flux |


---

## Input file formats

### Model atmosphere

SPAE supports multiple model atmosphere types via the keyword on the first line of the model file:

| Keyword | Code | Description |
|---------|------|-------------|
| `KURUCZ` | ATLAS | Most common; columns are ρx, T, P_gas, N_e, κ_Ross |
| `BEGN` | MARCS | Columns are τ_Ross, T, P_gas, N_e, µ, κ_Ross |
| `NEWMARCS` | MARCS (new) | Columns are τ_Ross, T, N_e, P_gas, ρ, v_turb, κ_Ross |
| `KURTYPE` | ATLAS (no κ) | Like KURUCZ but without pre-computed κ |
| `KUR-PADOVA` | Padova ATLAS | Columns are τ_ref, T, κ_ref, N_e, P_gas, ρ |
| `GENERIC` | generic | τ_ref depth scale, no pre-computed κ |

All model files share the same three-line header structure: type keyword, comment line (appears in output), and `NTAU nn` giving the number of depth layers. After the depth layers come the microturbulence (km/s or cm/s), the abundance overrides (`NATOMS n [M/H]` followed by atomic-number / log ε pairs), and the molecule list (`NMOL n` followed by molecule codes).

Example KURUCZ header:
```
KURUCZ
#OVER72: T= 5763,[g]=4.42,[Fe/H]=0.01,vt=1.27e+05
NTAU            72
 5.12e-04   3691.2  13.71  2.756e+09  2.659e-04  7.862e-02  2.000e+05
 ...
```

### Linelist

The first line of every linelist is a comment (appears in output). Each data line contains the following columns in order; optional columns may be left blank in formatted reads:

```
  wavelength  species   E_low    log(gf)   [C6]   [D0]   [EW]
```

- **wavelength** — Å. In `blends` mode, the first line of each blend is positive; continuation lines are negative.
- **species** — atomic number + ionization state after decimal (26.0 = Fe I, 26.1 = Fe II). For diatomic molecules, concatenate the two two-digit atomic numbers: CO = 608.0, CN = 607.0, MgH = 112.0.
- **E_low** — lower excitation potential (eV)
- **log(gf)** — oscillator strength (MOOG detects the sign to distinguish gf from log gf)
- **C6** — van der Waals damping constant (optional; if > 10⁻¹⁰ it multiplies the Unsöld value; if smaller it replaces it)
- **D0** — molecular dissociation energy in eV (required for molecular lines; omit or zero for atomic lines)
- **EW** — measured equivalent width in mÅ (required for `abfind`; omit or leave blank for `synth`)

Example (abfind linelist with EW):
```
Y1194_linelistExtended
  5052.167    06.000     7.685   -1.304    2.2    28.14    2.846
  6767.772    26.000     4.141   -3.878    2.2   130.62    3.875
  6842.690    26.100     5.526   -1.222    2.2    63.45    2.846
```

Example (synth atomic linelist, no EW):
```
a line list in the range:  6432.00  6442.00
  6432.128      22.0     3.291   -2.70
  6437.640      63.0     0.207    0.47
```

Example (synth molecular linelist with D0):
```
A CH synthesis list
  4195.423    106.00112    1.150    4.39E-02    3.47
  4195.771    106.00113    1.145    4.58E-02    3.47
  4196.208     26.0        3.40     3.556E-01
```

### Observed spectrum (synth only)

Two-column ASCII, wavelength (Å) and normalized flux, no header required:

```
  6401.777   0.9655
  6401.807   0.9945
  6401.836   0.9701
```

Scientific notation is also accepted (`6.40177686E+03  9.65545714E-01`).


---

## batch.par — parameter file reference

The first line of `batch.par` sets the mode. All keywords are case-insensitive. Paths may be quoted or unquoted.

### Common keywords

| Keyword | Description |
|---------|-------------|
| `standard_out` | Long-form output file (line-by-line results) |
| `summary_out` | Short-form summary output |
| `model_in` | Model atmosphere file |
| `lines_in` | Linelist file |
| `observed_in` | Observed spectrum (synth only) |
| `smoothed_out` | Output file for smoothed synthetic spectrum |
| `atmosphere` | Interpolation switch (0 = linear, 1 = log) |
| `molecules` | Molecular equilibrium (1 = on, 2 = off) |
| `lines` | Line-opacity treatment (0 = all, 1 = no scattering) |
| `flux/int` | 0 = flux, 1 = intensity |
| `damping` | van der Waals damping: 0 = Unsöld (default), 1 = Unsöld × 6.3, 2 = Unsöld × Blackwell factor |
| `freeform` | Linelist read mode: 0 = formatted 7e10.3 (blanks → 0), 1 = free-format (zeros must be explicit) |
| `strong` | 0 = no strong-line file (default), 1 = read a separate strong-lines file whose opacity is added at every synthesis step |
| `plot` | Plotting flag (ignored by pymoog; use interactive widget) |

### abfind-specific keywords

None beyond the common set. All lines in the linelist that have a non-zero EW are analyzed.

### synth-specific keywords

| Keyword | Description |
|---------|-------------|
| `synlimits` | Synthesis range and step: `λ_start λ_stop step delta_lambda` |
| `plotpars` | Display range and initial smoothing (see below) |
| `abundances` | Number of elements and passes; abundance offsets (see below) |
| `isotopes` | Isotope fractions per element per pass (see below) |
| `obspectrum` | Observed spectrum format (5 = two-column ASCII) |

**`synlimits`** — line immediately following the keyword:
```
synlimits
  6432.0  6442.0  0.02  1.00
```
Fields: λ_start (Å), λ_stop (Å), step (Å), delta_lambda (internal broadening half-width in Å).

**`plotpars`** — three lines following the keyword:
```
plotpars  1
  6432.0  6442.0  0.00  1.03      # display λ_min λ_max flux_min flux_max
  0.0  -1.23  0.000  0.00         # y-shift vrad (unused) continuum vshift
   g  0.090  0.0  0.0  0.0  0.0  # smtype fwhmgauss vsini limbdark vmac fwhmloren
```
`smtype`: `g` = Gaussian, `l` = Lorentzian, `v` = rotational, `m` = macroturbulent; combinations are also accepted (e.g. `gv`). Multiple broadenings can be combined simultaneously in the interactive widget regardless of what is set here.

**`abundances`** — format: `abundances  n_elements  n_passes`; one line per element with the element number followed by one offset per pass (log ε offset from the model value):
```
abundances  1  3
  63  -0.10  0.00  0.10
```
This runs three passes for Eu (Z=63) at model − 0.10, model + 0.00, and model + 0.10 dex.

**`isotopes`** — format: `isotopes  n_isotopes  n_passes`; one line per isotope with the isotope code and one inverse abundance ratio per pass. The isotope code is: atomic number (left of decimal) + ionization digit (first right of decimal) + three-digit atomic mass. For atoms: Li6 → `3.006`, Li7 → `3.007`, Eu-151 (neutral) → `63.0151`, Eu-151 (ionized) → `63.1151`. For molecules: C13O16 → `608.01316`.
```
isotopes  2  3
  63.1151  2.092  2.092  2.092    # Eu-151, same ratio all passes
  63.1153  1.916  1.916  1.916    # Eu-153
```


---

## Tutorial: equivalent-width abundance analysis (abfind)

This walkthrough derives abundances for a solar-type star from a list of measured equivalent widths.

### 1. Prepare your files

You need three files in a working directory:

```
run/
  batch.par
  star.mod
  linelist.txt
```

**`batch.par`:**
```
abfind
standard_out   'moog_out.1'
summary_out    'moog_out.2'
model_in       'star.mod'
lines_in       'linelist.txt'
atmosphere     0
molecules      1
lines          0
flux/int       0
damping        0
plot           0
```

**`linelist.txt`** — header line, then one line per measurement:
```
My star linelist
  5052.167    06.000     7.685   -1.304    2.2    28.14    2.846
  6767.772    26.000     4.141   -3.878    2.2   130.62    3.875
  6842.690    26.100     5.526   -1.222    2.2    63.45    2.846
```

### 2. Run the analysis

```bash
cd run/
pymoog batch.par
```

Or from Python:

```python
import os
os.chdir('run/')
from SPAE.moog import abfind_from_files
results = abfind_from_files('batch.par')
```

### 3. Read the output

`pymoog` writes two files. **`moog_out.2`** (summary) is most useful — one block per species:

```
Abundance Results for Species Fe 1/ 26.0
  average abundance = 7.451   std. deviation = 0.082   # lines = 55
  E.P. slope = -0.003   R.W. slope = 0.112
```

From Python, `results` is a dictionary:

```python
results['species']['26.0']
# {'mean': 7.451, 'sigma': 0.082, 'n': 55,
#  'ep_slope': -0.003, 'rw_slope': 0.112, ...}

results['lines'][0]
# {'wave': 5052.167, 'species': 6.0, 'abund': 8.43, 'ew': 28.14, ...}
```

### 4. Iterate on stellar parameters

Converged stellar parameters satisfy:
- No trend of abundance with excitation potential (EP slope ≈ 0) → constrains T_eff
- No trend with reduced equivalent width log(EW/λ) (RW slope ≈ 0) → constrains v_micro
- Fe I and Fe II give the same mean abundance (ionization balance) → constrains log g

Edit the microturbulence on the last line of `star.mod`, or swap in a new model, then re-run until these conditions are met.


---

## Tutorial: spectrum synthesis and abundance fitting (synth + interactive widget)

This walkthrough fits the abundance of a single element by overplotting synthetic spectra on an observed normalized spectrum.

### 1. Prepare your files

```
run/
  batch.par
  star.mod
  linelist.txt
  spectrum_norm.dat
```

**`spectrum_norm.dat`** — two-column ASCII, wavelength (Å) and normalized flux:
```
  6401.78   0.9655
  6401.81   0.9945
```

**`batch.par`** — three-pass Eu synthesis as an example:
```
synth
standard_out   'moog_out.1'
summary_out    'moog_out.2'
smoothed_out   'results.3'
model_in       'star.mod'
lines_in       'linelist.txt'
observed_in    'spectrum_norm.dat'
atmosphere     1
molecules      2
lines          1
flux/int       0
plot           2
abundances     1   3
  63   -0.10   0.00   0.10
synlimits
  6432.0  6442.0  0.02  1.00
plotpars  1
  6432.0  6442.0  0.00  1.03
  0.0  0.0  0.000  0.00
   g  0.090  0.0  0.0  0.0  0.0
obspectrum  5
```

To synthesize a single spectrum (no abundance grid), use `abundances 0 1` and omit the element line.

### 2. Launch the interactive widget

**From a terminal:**
```bash
cd run/
python -c "
from SPAE.moog.interactive import synth_interactive
w = synth_interactive('batch.par')
"
```

**From a Jupyter notebook:**
```python
%matplotlib widget          # must be in its own cell, run first

import os
os.chdir('run/')
from SPAE.moog.interactive import synth_interactive
w = synth_interactive('batch.par')   # keep w = to prevent garbage collection
```

### 3. Widget controls reference

| Control | Location | Function |
|---------|----------|----------|
| **Synthesize** | bottom-left button | Run synthesis with current `batch.par` abundances |
| **Reset** | bottom button | Reset all sliders to their initial values |
| **Save spectra** | bottom button | Write output files (see §4 below) |
| **Save batch.par** | bottom button | Write current abundances back to `batch.par` |
| **[+ Add element]** | bottom-right button | Add an element to the abundance slider list |
| **Pass 1 / 2 / 3 …** | pass buttons | Select which synthesis pass the abundance sliders edit |
| **◀ ▶** | pagination buttons | Cycle through elements when more than six are fitted (hidden when ≤ 6) |
| **Element slider** | element panel | Δ log ε from the model value; absolute log ε shown to the right |
| **FWHM Gauss** | smoothing panel | Gaussian broadening (Å) — instrumental + atmospheric |
| **vmac** | smoothing panel | Radial-tangential macroturbulence (km/s) |
| **vsini** | smoothing panel | Projected rotation (km/s) |
| **limb dark** | smoothing panel | Limb-darkening coefficient for rotation profile |
| **Δλ shift** | smoothing panel | Wavelength zero-point offset (Å); applied without re-synthesis |
| **FWHM Loren** | smoothing panel | Lorentzian broadening (Å) |
| **continuum** | smoothing panel | Multiplicative continuum scale (0.80–1.20); display only |
| **λ min / max** | zoom row | Set wavelength axis limits (press Enter to apply) |
| **F min / max** | zoom row | Set flux axis limits (press Enter to apply) |
| **Reset zoom** | zoom row | Restore automatic axis scaling |
| Cursor readout | spectrum panel | Shows λ (Å) and F at the current mouse position |

Multiple broadening types can be combined simultaneously (e.g. Gaussian + rotation).

### 4. Fitting workflow

1. Click **Synthesize** — all passes are computed and plotted over the observed spectrum.
2. Zoom in on the absorption line of interest using the λ/F boxes or the zoom tool.
3. Select **Pass 1** and adjust the element abundance slider until that synthetic profile matches the observed line. Repeat for each pass to bracket the best-fit abundance.
4. Tune broadening sliders until the synthetic line widths match. Start with FWHM Gauss (instrumental profile), then add vmac/vsini if needed.
5. Adjust **Δλ shift** if the synthetic and observed lines are not aligned in wavelength.
6. Use the **continuum** slider to match the continuum level if the observed spectrum was not perfectly normalized.
7. Click **Save batch.par** to record the fitted abundances.
8. Click **Save spectra** to write all output files.

### 5. Output files

Clicking **Save spectra** writes four files to the run directory. The base name is taken from `smoothed_out` in `batch.par` (or the `batch.par` stem if not set):

| File | Contents |
|------|----------|
| `<base>_abund.txt` | Pass number, element, Δ log ε, and absolute log ε for each fitted element |
| `<base>_spec.dat` | Wavelength (Å) + one smoothed flux column per pass |
| `<base>_obs.dat` | Observed wavelength and flux (written only if an observed spectrum was loaded) |
| `<base>_plot.py` | Self-contained Python script that reproduces the publication figure |

### 6. Making the publication figure

The saved `_plot.py` script requires only NumPy and matplotlib and can be run directly:

```bash
python <base>_plot.py
```

It produces a two-panel PDF (`<base>_figure.pdf`) with:
- **Top panel** — observed spectrum (black) + synthetic passes (color-coded), legend showing log ε for each pass
- **Bottom panel** — residuals (Obs − Syn) for each pass

Customize the figure by editing the configuration block at the top of the script:

```python
XLIM    = (6434.0, 6441.0)   # wavelength range (Å)
YLIM    = (0.910, 1.030)     # flux range
RLIM    = (-0.040, 0.040)    # residual range
COLORS  = ['#377eb8', '#e41a1c', '#4daf4a', ...]
OUTFILE = 'mystar_Eu.pdf'    # set to None to display interactively
```

The pass labels in the legend (`PASS_LABELS`) are auto-generated from the saved abundances and can also be edited freely.


---

## Running pymoog from Python

All modes are accessible programmatically without going through the CLI:

```python
from SPAE.moog import (
    abfind_from_files,
    synth_from_files,
    blends_from_files,
    ewfind_from_files,
    cog_from_files,
    weedout_from_files,
    doflux_from_files,
)

# Each function takes a path to batch.par and returns a results dict
results = abfind_from_files('batch.par')
results = synth_from_files('batch.par')
```

The interactive widget can also be instantiated directly for scripting:

```python
from SPAE.moog.interactive import SynthWidget
import matplotlib.pyplot as plt

w = SynthWidget('batch.par')
# Access w._wave, w._flux_smooth, w._obs_wave, w._obs_flux, etc.
plt.show()
```


---

## End-to-end pipeline: spectrum in, stellar parameters out (`SPAE.pipeline`)

`SPAE.pipeline` chains `SPAE.xspect_ew` (EW measurement) and `SPAE.spae`/`SPAE.analysis` (MCMC stellar-parameter fitting) into a single call, so a raw spectrum plus a linelist can go straight to a posterior on Teff/logg/[Fe/H]/micro without hand-wiring the two packages together.

```python
from SPAE.pipeline import run_full_pipeline

out = run_full_pipeline('star_blue.fits', 'linelist.txt', output_dir='run_star')
print(out['summary'])
```

The two stages are also available separately — `measure_star_ew()` (EW measurement only, writing a MOOG-format linelist) and `fit_stellar_params()` (MCMC fit from an already-measured linelist) — useful when you want to inspect or hand-edit the EW linelist between stages, or re-fit the same measurements with different `run_spae()` settings. `measure_star_ew(review=True)` pauses after measurement to show each flagged line's QC plot and let a reviewer keep it anyway (`review=False`, the default, is fully unattended). `solar_reference()` is a convenience wrapper that runs the same EW-measurement stage for the Sun and returns `sun_el, sun_abs` — the differential-abundance reference point `run_spae()`/`fit_stellar_params()` expect.

### `run_pipeline_cli.py` — differential (relative-to-solar) abundances for one star

A ready-to-edit command-line script that measures a target star and one or more solar reference exposures, filters both linelists down to their common lines (avoiding the array-position misalignment `SPAE.abunds.rel_abunds()` is otherwise exposed to when the two linelists disagree on which lines exist), and runs the MCMC fit differentially against the Sun:

```bash
python run_pipeline_cli.py /path/to/target.fits [/path/to/solar_dir_or_file]
```

Omitting the solar argument uses the bundled solar sample spectra (`SPAE/xspect_ew/data/spectra_sample/Solar/`). Edit the `CONFIG` block at the top of the script first — in particular `LINELIST_PATH`. Outputs (QC plots, EW linelists, MCMC diagnostics, trace/corner plots) are written inside the target spectrum's own directory; see the script's module docstring for the full file list. Uses multiprocessing via `run_spae()`, so run it as a script (`python run_pipeline_cli.py ...`) rather than importing and calling its internals directly.
