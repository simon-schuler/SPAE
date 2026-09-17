# SPAE.moog / pymoog: development record

Durable, third-person technical record of real development/validation work on
`SPAE/moog` (the pure-Python MOOG conversion) and its SPAE integration
(`SPAE/spae.py`, `SPAE/abunds.py`, `SPAE/analysis.py`) -- kept for paper
support, the same role `SPAE/xspect_ew/DEVELOPMENT_LOG.md` plays for that
subsystem. This file starts with the first real-star validation round;
earlier phases of the F77-to-Python conversion itself (state.py/physics
layer, abfind/synth/all-modes drivers, the numba performance pass, the
`SPAE.analysis` post-processing module) predate this file and are not yet
backfilled here.

## 1. Seven-star Kepler host-star validation (2026-09-16/17)

**Purpose**: validate the pymoog implementation (as opposed to unit-level
Fortran-reference checks, e.g. the existing `Verification/Claude_MOOG/`
164-line test) against a real, previously-published multi-star analysis,
run end-to-end through SPAE's own MCMC pipeline (`run_spae()`).

**Sample**: the seven Kepler planet-host stars from Schuler et al. (2015,
ApJ 815, 5; arXiv:1511.00934) -- Kepler-20, -21, -22, -37, -68, -100, -130
-- each already present in this project's `Verification/SWPs/Schuler/`
tree with a matched `linelist_star.txt`/`linelist_sun.txt` pair (144-192
lines, 19-21 species) from that prior (Fortran-MOOG-era) analysis.

**Setup**: new `Verification/KepTest/` directory, one subdirectory per
star, each with its own copy of `linelist_star.txt`/`linelist_sun.txt`
and a driver script (`kepXXX_run1.py`) following the pattern already
validated for Kepler-22 alone in an earlier round (`kep22_run2.py`,
confirming the numba speedup at full-run scale) -- independent-per-star
`run_spae()` call, 1000 steps/40 walkers/8 cores, `ep_slope_scale=
rew_slope_scale=0.010`. `x_0` for six of the seven stars was seeded from
that star's own prior (Fortran-MOOG-era) solution, read directly from the
first line of the original `kepXXX_abunds.csv` files in
`Verification/SWPs/Schuler/`. All 7 pre-flight-checked with a single fast
`abunds_func()` call per star/Sun before committing to the full runs.

Run sequentially (not in parallel) inside one `screen` session with
`caffeinate -s` tied to its lifetime, per this project's established
unattended-run pattern -- this machine (Apple M3 Pro, 12 cores) cannot
support 7 concurrent 8-core runs without the same kind of resource
contention documented earlier in `[[project-moog-conversion]]`'s Dropbox-
contention investigation. All 7 completed cleanly (exit code 0), 34-40
min each, ~4h16m total.

### 1.1 A stale `x_0`, not a pymoog bug (Kepler-20)

Comparing each star's new `analyze_run()` traditional (excitation/
ionization balance) solution against its `x_0` reference showed six
stars in excellent agreement but Kepler-20 off by 269 K / 0.45 dex (log
g) / 0.28 dex ([Fe/H]) / 0.64 km/s (vt) -- initially looked like a real
implementation problem. Root cause, found via `Verification/SWPs/
Schuler/Kepler20/kep20_notes.txt`: Kepler-20 was "the first star run...
and the first for the new SWP test of SPAE, [requiring] five runs...
adjusting the code to work on helium" -- its `kepXXX_abunds.csv` "Run 1"
entry (the source of `x_0` for every other star) is an early debugging
attempt on a different code path, not the adopted solution. The real
adopted model atmosphere (`Kepler20/kep20.mod` header: `T=5555,
[g]=4.53,[Fe/H]=0.08,vt=1.27e+05` cm/s) is close to the actual published
value (Table 2 below) and nowhere near "Run 1"'s stale numbers. The
MCMC's own posterior still recovered a good answer despite the bad
starting point (see below) -- a mild positive sign for the sampler's own
robustness, not just luck.

### 1.2 Burn-in: 100 for six stars, 200 for Kepler-20

Every run's first-pass `analyze_run()` used the provisional
`burn_in=100` default. After visual trace inspection (not detailed
further here -- a user judgment call, not an automated one), `burn_in=
100` was confirmed adequate for six of the seven stars; Kepler-20 needed
`burn_in=200`. Re-analysis used the ALREADY-SAVED sampler pickle (no new
MCMC run needed) -- `az.analyze_run(sampler, ..., burn_in=200)` +
`az.summarize(result, provisional=False)`. The traditional-balance
solution barely moved (Teff 5520.0->5520.1 K) since it only draws from
posterior samples independently satisfying strict `|EP_r|,|REW_r|<0.01`
constraints, largely insensitive to burn-in past actual convergence; the
median solution's confidence interval tightened somewhat (Teff
-14.7/+16.2 -> -14.3/+13.9), consistent with discarding more of the
initial transient. Final, non-provisional summaries for all 7 stars
saved as `<star>_1_final_summary.txt` in each `KepTest/<Star>/` directory
(`kep20_1_final_summary.txt` uses the burn_in=200 re-analysis).

### 1.3 Comparison against Schuler et al. (2015)'s actual published Table 2

An initial comparison against each star's OWN internal reference (the
`kepXXX_abunds.csv` "Run 1" x_0, or `kep20.mod`'s header for Kepler-20)
looked good but was not the right ground truth -- confirmed by fetching
the real paper (arXiv:1511.00934, Table 2) and comparing directly.

**Traditional (excitation/ionization balance) solution vs. published:**

| Star | Published Teff/logg/[Fe/H]/vt | New pymoog (traditional) | dTeff (K) | dlogg | d[Fe/H] | dvt (km/s) |
|---|---|---|---|---|---|---|
| Kepler-20 | 5514/4.44/+0.059/1.15 | 5520.1(+/-2.0)/4.436(+/-0.005)/+0.060(+/-0.004)/1.167(+/-0.004) | +6.0 | -0.004 | +0.001 | +0.017 |
| Kepler-21 | 6177/3.99/-0.079/1.98 | 6182.3(+/-1.8)/4.030(+/-0.007)/-0.083(+/-0.004)/2.081(+/-0.006) | +5.3 | +0.040 | -0.004 | +0.101 |
| Kepler-22 | 5622/4.57/-0.251/1.47 | 5621.5(+/-1.3)/4.572(+/-0.005)/-0.252(+/-0.004)/1.459(+/-0.005) | -0.5 | +0.002 | -0.001 | -0.011 |
| Kepler-37 | 5406/4.49/-0.323/1.36 | 5399.7(+/-1.3)/4.447(+/-0.006)/-0.328(+/-0.003)/1.321(+/-0.004) | -6.3 | -0.043 | -0.005 | -0.039 |
| Kepler-68 | 5887/4.45/+0.126/1.57 | 5885.4(+/-1.6)/4.432(+/-0.006)/+0.123(+/-0.006)/1.572(+/-0.003) | -1.6 | -0.018 | -0.003 | +0.002 |
| Kepler-100 | 5855/3.98/+0.061/1.49 | 5861.3(+/-2.6)/3.986(+/-0.006)/+0.071(+/-0.005)/1.479(+/-0.002) | +6.3 | +0.006 | +0.010 | -0.011 |
| Kepler-130 | 5958/4.41/-0.199/1.78 | 5964.1(+/-0.9)/4.398(+/-0.015)/-0.196(+/-0.003)/1.742(+/-0.005) | +6.1 | -0.012 | +0.003 | -0.038 |

**MCMC median solution vs. published:**

| Star | Published Teff/logg/[Fe/H]/vt | New pymoog (MCMC median) | dTeff (K) | dlogg | d[Fe/H] | dvt (km/s) |
|---|---|---|---|---|---|---|
| Kepler-20 | 5514/4.44/+0.059/1.15 | 5523.4/4.449/+0.063/1.179 | +9.4 | +0.009 | +0.004 | +0.029 |
| Kepler-21 | 6177/3.99/-0.079/1.98 | 6194.9/4.024/-0.071/2.005 | +17.9 | +0.034 | +0.008 | +0.025 |
| Kepler-22 | 5622/4.57/-0.251/1.47 | 5622.6/4.566/-0.247/1.405 | +0.6 | -0.004 | +0.004 | -0.065 |
| Kepler-37 | 5406/4.49/-0.323/1.36 | 5412.5/4.484/-0.316/1.316 | +6.5 | -0.006 | +0.007 | -0.044 |
| Kepler-68 | 5887/4.45/+0.126/1.57 | 5892.9/4.453/+0.130/1.566 | +5.9 | +0.003 | +0.004 | -0.004 |
| Kepler-100 | 5855/3.98/+0.061/1.49 | 5865.8/3.997/+0.067/1.506 | +10.8 | +0.017 | +0.006 | +0.016 |
| Kepler-130 | 5958/4.41/-0.199/1.78 | 5962.1/4.393/-0.196/1.729 | +4.1 | -0.017 | +0.003 | -0.051 |

**Conclusion**: every star, on both solution types, agrees with Schuler et
al. (2015)'s published values well within that paper's own quoted
uncertainties (Teff: 25-45 K; [Fe/H]: 0.04-0.08 dex) -- Teff to within
~18 K, log g to ~0.04 dex, [Fe/H] to ~0.01 dex, vt to ~0.1 km/s, across
every star. As expected, the traditional (excitation/ionization balance)
solution agrees slightly more tightly than the broader MCMC median
posterior summary, but both validate the pymoog implementation cleanly
against a real, independent, previously-published multi-star analysis --
a stronger test than the single-star/single-model-atmosphere Fortran-
reference check alone.

Full run outputs (samplers, logs, trace plots, final summaries) live in
`Verification/KepTest/<Star>/` -- not tracked in this repository (see
`[[feedback-experiment-tracking]]`), this file is the durable record of
the run configuration and results.

### 1.4 MRT paper-table exports

With all 7 stars validated and finalized (§1.2-1.3), generated a
per-star CDS/MRT export (`SPAE.analysis.write_mrt()`) for each --
`Verification/KepTest/<Star>/<star>_1.mrt`, using each star's own
confirmed `burn_in` (100 for six stars, 200 for Kepler-20). All 7
round-trip cleanly through astropy's `Table.read(..., format=
'ascii.cds')` reader. This is the first real multi-star use of
`write_mrt()` (previously only smoke-tested on a single star, Kepler-22)
-- since the function is deliberately per-star/per-run with no bulk API
(different stars can need different `burn_in`), the natural pattern is a
loop varying `burn_in` per star while reusing each star's already-saved
sampler pickle (no MCMC re-run needed). This pattern is now documented
in `SPAE_QUICKSTART.md`'s "Step 8 (optional) -- a multi-star batch".
