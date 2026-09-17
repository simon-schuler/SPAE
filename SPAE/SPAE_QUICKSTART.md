# SPAE Quick Start

This is a quick start guide for running **SPAE** itself — the Bayesian MCMC
layer that derives stellar parameters (Teff, log g, [Fe/H], microturbulence)
and their uncertainties. SPAE is built on **pymoog**, a pure-Python
reimplementation of the MOOG spectral analysis engine; see `README.md` for
`pymoog`/`batch.par`/linelist-format details if you need to work with that
layer directly. This guide assumes SPAE is already installed (`pip install
./` from the project root — see `README.md`).

This is the first piece of a fuller SPAE/MOOG manual; more sections
(interactive synthesis widget, `pymoog` CLI, etc.) are covered separately in
`README.md` for now.


---

## Part 1 — Running SPAE

### 1. Prepare your linelists

SPAE derives **differential** abundances (star minus Sun) line by line, so
the star and solar linelists must be matched — every line used for the star
must have a corresponding line at the same wavelength/species in the solar
linelist. Linelist file format (wavelength, species, EP, log gf, ...) is the
same as `pymoog`'s — see the "Linelist" section of `README.md`.

If you have independent raw equivalent-width linelists for the star and Sun
(not yet matched), build the matched pair with:

```python
from SPAE.read_write import linelist_create

linelist_create('ocsn49star_ew.txt', 'ocsn49sun_ew.txt', direc_path)
# writes matched linelist_star.txt / linelist_sun.txt into direc_path
```

If you already have matched, line-for-line linelists in MOOG format, skip
this step and use them directly.

### 2. Get solar reference abundances

```python
from SPAE.abunds import abunds_func

x_sun = 5777, 4.44, 0.00, 1.38   # Teff, logg, [Fe/H], microturbulence
sun_el, sun_abs = abunds_func(x_sun, linelist_sun)
```

(`SPAE.abunds` also has a `sun_abs(sun_linelist, teff_sun=5777, logg_sun=4.44,
feh_sun=0.00, micro_sun=1.38)` convenience wrapper that does the same thing
in one call, if you don't need to override the solar linelist path
separately.)

### 3. Run the MCMC

```python
from SPAE.spae import run_spae

sampler, flat_blob, log = run_spae(
    linelist_star,
    sun_el=sun_el, sun_abs=sun_abs,
    n_steps=1000,
    n_cores=8,
    ep_slope_scale=0.010,
    rew_slope_scale=0.010,
)
```

| Parameter | Default | Meaning |
|---|---|---|
| `linelist` | required | Path to the star's matched MOOG linelist |
| `sun_el`, `sun_abs` | `None` | Solar reference from step 2; enables differential mode |
| `x_0` | `(5777, 4.44, 0.01, 1.38)` | Initial guess: (Teff, logg, [Fe/H], micro) |
| `n_walkers` | `40` | Number of emcee walkers |
| `n_steps` | `1000` | MCMC steps per walker |
| `n_cores` | `None` (all cores) | `1` disables multiprocessing |
| `ep_slope_scale` | `0.010` dex/eV | Penalty scale on the Fe I excitation-potential slope; smaller = tighter excitation-balance constraint |
| `rew_slope_scale` | `0.010` dex/dex | Penalty scale on the Fe I reduced-EW slope; smaller = tighter ionization/curve-of-growth constraint |
| `include_prior` | `False` | Add a spectroscopic prior on Teff/logg |

Install `tqdm` to get a live per-step progress bar (`pip install tqdm`) —
without it `run_spae()` runs silently until it finishes. A 1000-step/40-walker
run typically takes on the order of two hours on 8 cores, depending on
linelist size; scale roughly linearly with `n_steps`.

### 4. Save your results

```python
import pickle
import numpy as np

with open(f'{star}_sampler_{run}.pkl', 'wb') as f:
    pickle.dump(sampler, f)
with open(f'{star}_{run}.log', 'w') as f:
    f.write(log)
np.save(f'{star}_flat_blob_{run}.npy', flat_blob)
```

Save the **sampler itself** (via pickle) — that's what every post-run
analysis function below actually needs. `flat_blob` is only a convenience
cache for quick manual inspection; it isn't used by any `SPAE.analysis`
function.


---

## Part 2 — Post-run tasks

All of the following use `SPAE.analysis` (`import SPAE.analysis as az`)
and start from the pickled `sampler` object saved above.

### Step 1 — quick provisional look

It's natural to want an immediate look right after the MCMC finishes. But
the correct burn-in is a **visual judgment call** from the trace plot (Step
2) — so treat anything computed before that as provisional:

```python
result = az.analyze_run(sampler, linelist_star, sun_el=sun_el, sun_abs=sun_abs,
                         burn_in=100, a_frac=0.40, max_cor=0.01)
print(az.summarize(result, provisional=True))
```

`provisional=True` prepends a banner to the summary flagging `burn_in=100`
as a hardcoded default, not a real convergence-based answer, and pointing
back to this recipe.

### Step 2 — pick the real burn-in from the trace plot

```python
fig, axes = az.plot_trace(sampler.get_chain())
```

`plot_trace()` deliberately takes the **full, unclipped** chain — burn-in
is a visual judgment call. Look for the step number where all four
parameter panels (Teff, log g, [Fe/H], micro) stop drifting and settle into
a steady band across all walkers; that step number is your `burn_in`.

(Automating this via emcee's `get_autocorr_time()` was considered, but it's
unreliable on runs not much longer than the autocorrelation time — exactly
the case where visual inspection matters most — so it isn't used here.)

### Step 3 — final analysis with the chosen burn-in

```python
result = az.analyze_run(sampler, linelist_star, sun_el=sun_el, sun_abs=sun_abs,
                         burn_in=<your value>, a_frac=0.40, max_cor=0.01)
print(az.summarize(result))
```

What `analyze_run()` does, in order:
1. **`clip_run`** — discards the first `burn_in` steps, and drops entire
   walkers whose acceptance fraction is `<= a_frac` (not just truncates them).
2. **`median_solution`** — posterior median Teff/logg/[Fe/H]/micro (with
   16th/84th-percentile offsets) plus a per-species `[X/H]` table.
3. **`balance_stats`** (at the median) — EP/REW slope diagnostics.
4. **`traditional_solution`** — the excitation/ionization-balance solution:
   the mean of unique posterior samples satisfying `|EP_r|, |REW_r| <
   max_cor` and Fe I = Fe II, i.e. the classical "balance" result rather
   than the full posterior median.

Per-species table columns: `mean` ([X/H]), `std` (scatter across posterior
samples), `sigma_mean` (line-to-line scatter within a single sample), and
`tot_uncert` (`sqrt(std**2 + sigma_mean**2)`, combining both).

### Step 4 (optional) — refine the traditional solution

```python
result = az.analyze_run(sampler, linelist_star, sun_el=sun_el, sun_abs=sun_abs,
                         burn_in=<your value>, refine_p_threshold=0.94)
```

Costs one extra MOOG call per unique traditional-solution sample — only use
it if you need the tightest possible traditional solution (keeps only
samples whose EP/REW slopes are consistent with zero at the given p-value).

### Step 5 (optional) — diagnostic plots

```python
fig, ax    = az.plot_balance(result.median_balance)                       # EP/REW correlation at the median
fig, axes  = az.plot_corner(result.clipped)                               # posterior corner plot (needs `corner`)
fig, ax    = az.plot_species_histogram(result.clipped.flat_blob, 'Fe I ', 'Fe II ')
```

### Step 6 — write the text summary to a file

```python
with open(f'{star}_{run}_summary.txt', 'w') as f:
    f.write(az.summarize(result))
```

### Step 7 (optional) — export for a paper table

```python
az.write_mrt(result, name=star, path=f'{star}_{run}.mrt')
```

Writes a single-row CDS/MRT-format file (median **and** traditional
solutions side by side, per-species `[X/H]`/uncertainty/line-count columns)
built for easy transfer into a paper table.

**Call this deliberately, not after every run.** Most runs are exploratory
(tuning `ep_slope_scale`/`rew_slope_scale`, testing initial guesses, etc.);
an MRT file represents "this is the paper-table row for this star," so
generate it once you've settled on a star's final run.

### Step 8 (optional) — a multi-star batch

`write_mrt()` is a per-star, per-run call — there's no bulk API, since
different stars in a sample often need different `burn_in` (see Step 3)
even when everything else about the run is identical. For a batch of
already-completed runs (samplers already saved to disk, no need to re-run
the MCMC), loop over stars and just vary `burn_in` per star:

```python
import pickle
from SPAE.abunds import abunds_func
from SPAE import analysis as az

# per-star: (directory, star id, linelist dir, confirmed burn_in)
stars = {
    'Kepler100': ('kep100', 100),
    'Kepler20':  ('kep20',  200),   # this one needed a longer burn-in
    # ... etc
}
x_sun = 5777, 4.44, 0.00, 1.38

for star_dir, (star, burn_in) in stars.items():
    d = f'{star_dir}/'
    sun_el, sun_abs = abunds_func(x_sun, d + 'linelist_sun.txt')

    with open(d + f'{star}_sampler_1.pkl', 'rb') as f:
        sampler = pickle.load(f)

    result = az.analyze_run(sampler, d + 'linelist_star.txt',
                             sun_el=sun_el, sun_abs=sun_abs, burn_in=burn_in)
    az.write_mrt(result, name=star_dir, path=d + f'{star}_1.mrt')
```

Each star still gets its own MRT file (Step 7's one-row-per-file design is
unchanged) — this is just the loop that avoids re-running MCMC for stars
whose sampler is already on disk, while still letting each star use its
own trace-inspected `burn_in`. Combining the resulting per-star files into
one paper table (`astropy.table.vstack` after reading each back with
`format='ascii.cds'`) is, as noted in Step 7, the caller's own job.

Validated this way on the 7-star Kepler host-star sample (Schuler et al.
2015) — see `SPAE/moog/DEVELOPMENT_LOG.md` §1 for the full validation
writeup and the resulting comparison against that paper's published
parameters.


---

## Reference example

`Verification/SPAE_test/star.py` is a complete, working example combining
all of the steps above (linelist matching → MCMC → save → provisional
analysis). It hardcodes `burn_in=100` for its automatic post-run summary,
consistent with Step 1 above — re-run Steps 2–3 manually once you've
inspected the trace plot for a real result.
