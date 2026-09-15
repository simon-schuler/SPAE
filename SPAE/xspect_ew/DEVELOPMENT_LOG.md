# SPAE.xspect_ew: Development and Validation Record

This document is a running technical record of the design, bugs found,
fixes made, and validation evidence for `SPAE.xspect_ew`, intended as
source material for the Methods/Software section of a future paper. It
complements (does not replace) the in-code docstrings, which describe
the *current* design rationale in detail; this document preserves the
*chronology* — what was tried, what failed and why, what fixed it, and
the quantitative evidence for each step — which the docstrings
deliberately do not carry.

Git history referenced below is on branch `pymoog-conversion`,
`git@github.com:simon-schuler/SPAE.git`. Commit hashes are given in the
order they were made.

---

## 1. Provenance and incorporation

`XSpect-EW` originated as an independent tool by George Vejar
(github.com/forgeousgeorge/XSpect, published on PyPI). It was
incorporated into SPAE as `SPAE.xspect_ew` (commit `465ea88`), carrying
forward two small unreleased fixes present in the original author's
local development copy but not yet in the published PyPI release (a
`verbose` flag, and a crash guard in the wavelength-shift estimator for
orders with no reference-spectrum match). Incorporation was verified
end-to-end against the bundled solar HIRES sample spectrum and Fe
linelist: load → normalize → measure Fe I 4779.439's equivalent width
(40.84 ± 0.50 mÅ, consistent with the sample linelist's reference value
of 39.5 mÅ) → confirm the MOOG-format output parses cleanly with
`SPAE.moog.inlines.parse_linelist()`.

The single 1324-line source file was then restructured (commit
`e92e5c7`) into topic modules mirroring `SPAE.moog`'s one-file-per-
concern convention: `constants.py`, `line_profile.py`, `gp_utils.py`,
`continuum.py`, `combine.py`, `plotting.py`, `io_utils.py`, and the
`Spectrum_Data` driver class in `spectrum_data.py`.

## 2. Multi-instrument format support

`Spectrum_Data` originally read only classic Keck/MAKEE FITS output.
`readers.py` (commit `9f439b7`) generalized this to an explicit,
extensible registry of (format-detector, extractor) pairs, verified
against one real file per instrument:

- **Keck/MAKEE**: the pre-existing `CRVL1_NN`/`WAT2_NNN` header logic,
  unchanged, moved out of `Spectrum_Data`.
- **GRACES/OPERA**: self-documenting via `COL1`..`COLn` FITS header
  keywords; orders are split on the `Order` data column rather than a
  fragile "wavelength decreasing" heuristic.
- **MAROON-X**: HDF5/pandas-HDFStore format (`spec_blue`/`spec_red`
  keys, MultiIndex `(Fiber, Order)`, `optimal_extraction` column). Which
  fiber index is the science target is not recoverable from file
  metadata (the header's `FIBERn` fields use different numbering than
  the DataFrame's own `Fiber` index); defaults to fiber 6, overridable.
- **Generic FITS binary table** with named wavelength/flux columns.

Design principle, applied throughout: prefer raw/uncalibrated flux and
uncorrected wavelength over any pre-normalized or pre-corrected
quantity a format may also offer, since `Spectrum_Data` performs its
own continuum normalization and wavelength-shift correction.

## 3. Wavelength shift and radial velocity

Two issues were found and fixed in the coarse-grid wavelength-shift
search (`estimate_shift()`, commit `99c7207`): its per-trial-shift chi²
was unnormalized by the number of overlapping points (biasing toward
shifts that minimized window overlap), and it had no sub-grid
refinement. Fixing both (`np.mean` instead of `np.sum`; a free 3-point
parabolic refinement) improved recovery of a known synthetic +0.30 Å
offset from 0.024 Å error to 0.0003 Å error — roughly 80× — at the same
grid resolution/compute cost.

This preceded a larger change: rather than continue hardening the
per-order reference-spectrum matching heuristic, a single effective
radial velocity is now measured from a handful of strong, well-
identified lines (Ca II H&K, Balmer series, Mg b, Na D;
`radial_velocity.py`, commit `8c3b32c`) and applied to every order via
correct multiplicative `(1 + v/c)` scaling
(`Spectrum_Data.apply_rv_shift()`). This was motivated by measuring the
actual EW-measurement impact of the two methods on real data: running
the full `measure_all_ew()` pipeline on a real GRACES spectrum with the
bundled 78-line Fe linelist, 56/78 lines (72%) showed a >20 mÅ line-
center disagreement between the two shift methods, with EW differences
up to ~11× on the same nominal line. Comparing both methods against
each line's independently-predicted position (from the cross-validated
RV) showed the RV-based method tracked substantially tighter (RMS 37
mÅ, max 63 mÅ) than the older extrapolation-based method (RMS 67 mÅ,
max 119 mÅ) for orders lacking direct reference-spectrum overlap. The
RV-based method is now the documented default; the original method
remains available and correct for its original use (reference-spectrum
cross-correlation).

## 4. EW flagging infrastructure

Rather than silently including a suspect measurement, `Spectrum_Data`
flags per-line issues and routes flagged lines to a companion file
(`check_for_flags()`, `make_ew_doc()`; commit `44a33e6`). Checks
include: EW fractional error > 10%, EW < 2 mÅ (too shallow to trust), a
poor Gaussian fit (χ² above threshold), and a found line center more
than `position_thresh` (0.07 Å, calibrated against the real 78-line
GRACES test in §3) from the line's rest wavelength — the last check
specifically targets likely misidentification (a silently wrong line,
not just an imprecise measurement), which is the failure mode most
worth catching automatically. `make_ew_doc()` writes two files: the
main linelist (only lines that passed every check) and a companion
`..._flagged` file (every flagged line, annotated with its flag
reason(s), for manual review). This infrastructure was later extended
to include cross-order overlap disagreements (§8).

## 5. Continuum-fitting algorithm: abandoned approaches

Four approaches were tried and abandoned before the current algorithm,
each because it reproduced (or introduced a new form of) the same
underlying failure: a stretch of an order left insufficiently
constrained by whatever hard selection rule was in use.

1. **`Continuum_scan`'s local-window percentile selection** (the
   original inherited method): a narrow window sitting mostly or
   entirely inside a broad/strong line has no true-continuum points to
   select at all. Measured 20-70% continuum underestimate near a
   synthetic 3-Å-wide line.
2. **A globally-refit Gaussian Process**, iteratively sigma-clipped
   against the same kind of local selection (commit `e14e912`'s
   history, superseded within the same commit): did not fix the
   broad-line bug at all (-26%, max 78.6% error on the synthetic test)
   — a GP's local flexibility bends into a contaminated region just as
   easily as the window selection did, and one-sided sigma-clip
   rejection cannot self-correct once the fit is already biased low
   across a contiguous stretch.
3. **A single global low-order polynomial** (commit `e14e912`, an
   earlier state, since superseded): fixed the broad-line bug
   (-0.7%, max 0.8% on the synthetic test — every point constrains the
   whole curve, so it cannot be dragged into one local gap) but was too
   rigid to track real order-wide continuum curvature. Raising the
   polynomial degree to compensate made this worse (degree 5: -5.6%,
   max +7% broad-line error, up to 14% error in nominally clean
   regions) via Runge's-phenomenon-style oscillation through gaps left
   by rejected points.
4. **A piecewise cubic spline** with knots at fixed spacing: tracked
   curvature well on the synthetic test (~0.2% error) but reintroduced
   selection-mask fragility in three distinct real-data forms as its
   initial-guess and edge-handling were patched: (a) a milder version
   of the original broad-line bias, (b) catastrophic divergence
   (-1.2 million counts against a ~50,000-100,000-count real order) at
   a real Keck order edge where too few points constrained the
   outermost spline segment, corrupting two real lines' EW
   measurements, and (c) a real GRACES order's genuine ~5× large-scale
   continuum decline being excluded from the selection mask entirely
   (the same mechanism that correctly excludes absorption lines,
   wrongly applied to a real trend). An "edge guard" (holding the fit
   flat beyond an adaptively-grown margin from each domain edge) was
   built and iterated through three versions, each fixing the case that
   motivated it while breaking a different real order, before being
   abandoned as a fundamentally unsound strategy rather than patched a
   fourth time.

Each of these four failures shares a root cause: a **hard** point
selection (definitely-continuum vs. definitely-not) that, whenever it
leaves a stretch of an order with too few or zero trusted points,
leaves whatever is fit to that stretch either unconstrained or
constrained by unrepresentative data.

## 6. Current continuum-fitting algorithm

Replaced entirely (commit `dba6464`) with Asymmetric Least Squares
(AsLS) smoothing (Eilers & Boelens 2005; standard in Raman/chromatography
baseline correction, here inverted to fit an upper envelope rather than
a baseline). Implemented in `continuum.py`, function
`fit_als_continuum()`. Every point receives a **soft**, iteratively
re-estimated weight — points above the current fit (likely continuum)
are upweighted, points below (likely absorption) are downweighted, but
never fully excluded — so no stretch of an order can become completely
unconstrained by construction, closing off the entire failure family in
§5 at once rather than patching around it.

The algorithm solves, each iteration,

```
(W + Dᵀ diag(λ_vec) D) pred = W · flux
```

where `D` is the discrete second-difference operator (the smoothness
penalty, replacing knot placement/polynomial degree as a *continuous*
regularizer) and `W` holds each point's weight: base inverse-variance
weight (`1/err²`), further modulated by which side of the current fit
the point falls on. As of the current state (below), several
refinements sit on top of this base scheme, each added in response to a
specific, confirmed real-data failure:

### 6.1 Adaptive per-point stiffening (two signals)

A single global smoothness penalty λ cannot simultaneously resist a
densely-blended real stretch (needs stiff) and track a genuine
large-scale continuum decline (needs flexible) — confirmed on a real
Keck order (6198-6309 Å) with a ~30 Å densely-blended Fe I stretch
alongside a real GRACES order (5698-5898 Å) with a genuine ~5× decline
toward its edges: no single λ served both.

`fit_als_continuum()` computes a per-point *severity* signal from the
raw data (independent of the fit itself, so it cannot inherit the
fit's own bias) and multiplies λ up to `stiffen_factor`× wherever
severity is high, combining two distinct sub-signals via `max()`:

- **Signal 1 (broad/deep troughs)**: compares each point's local peak
  (max flux in `local_window`, default 3 Å) against a wider window's
  peak (`wide_window`, default 25 Å). A genuine large-scale slope does
  not trigger this (both windows' peaks decline together, since
  `wide_window` is much narrower than any slope of interest); a broad
  line or saturated telluric band does (the local peak stays well below
  the wide peak over an extended stretch).
- **Signal 2 (widespread, individually-modest blending)**, added later
  (commit `0a08008`): Signal 1 only asks whether *some* point nearby
  reaches true continuum — satisfied by a single bright pixel even when
  95-99% of surrounding points sit measurably below it. Confirmed
  invisible on real MAROON-X data (two independent "clean" regions,
  §9.3-9.4). Computed by dividing flux by a wide-window MEDIAN trend
  (removing a genuine large-scale slope's own contribution) and then
  comparing the LOCAL median to the LOCAL max of that ratio — both over
  the same narrow window, so a genuine slope's variation across just
  that window cannot masquerade as blending. Uses its own, much smaller
  threshold/scale (`severity_threshold=0.02`, `severity_scale=0.06`)
  than Signal 1's fixed 0.08/0.92, since this shortfall saturates far
  below 1.0 (typically 0.03-0.08) even for an obviously real forest,
  unlike a genuine broad trough (0.3-1.0).

### 6.2 Noise-scaled below-fit weighting

A hard weight step at the current fit (flat weight `p` for anything
below it) biases the fit toward the upper envelope of the noise rather
than its mean, since real Poisson noise scatters roughly half of
genuinely unabsorbed points below the current estimate at any moment.
Confirmed on a real Keck order (5698-5800 Å): the fit tracked at or
above the local peak everywhere checked. Replaced with a smooth,
noise-scaled decay (points within `low_reject_sigma`, default 2.5σ, of
the fit are trusted close to fully; further below, trust decays toward
`p`).

A **symmetric** version (also downweighting above-fit outliers) was
tried and reverted: it caused a self-reinforcing lockout, where genuine
continuum points adjacent to a still-recovering, transiently-low fit
looked like implausible upward outliers and were wrongly downweighted
too, preventing recovery near any strong/broad line (confirmed
directly: a real continuum point at 49,667 counts, next to a synthetic
broad line, was assigned z ≥ 20 against a fit still at ~45,000, and
rejected as if it were a cosmic ray).

### 6.3 Cumulative low-`p` bias correction

Even with 6.1-6.2, the below-fit weight floor `p` (default 0.01) never
reaches zero; in densely-lined regions the cumulative pull of many
small non-zero weights, coupled through the smoothness penalty across
the whole order, produces a measurable systematic UNDERSHOOT even in
stretches with no single point locally flagged as an outlier. Confirmed
via a synthetic test varying only line density (0/30/100/250 scattered
lines against a known-flat continuum): error grew from 0% to -4.4% with
density, and was confirmed not a convergence artifact (identical result
at n_iter = 15 vs. 120). Fixed (commit `4e6c5f4`) by reusing the same
severity signal from §6.1 to also shrink `p` by up to
`p_reduction_factor` (default 30) in high-severity regions — `p` cannot
simply be made small everywhere, since that degrades tracking of a
genuine slope the fit has not yet caught up to (confirmed: naively
shrinking `p` globally worsened the GRACES-slope synthetic test's
clean-region error from ~4% to ~15%).

### 6.4 Percentile targeting (upper-envelope offset)

Even with 6.1-6.3, `low_reject_sigma`'s width (2.5σ) barely discounts
anything within ~2σ of the fit — the bulk of a Gaussian — so despite
nominally asymmetric weighting, the converged fit sits near the noise
**mean**, not a genuine upper envelope. This under-corrects for real
pervasive weak-line blanketing (the true, line-free continuum sits
above the observed centroid nearly everywhere in a real spectrum, not
just where a resolved line exists), and was the root cause of both a
"segments of spectrum sit above the fitted continuum" symptom and,
compounded with real line density, much of an observed Keck red-order-
edge droop.

Added (commit `e2be2e2`) a post-hoc `target_percentile` parameter
(default 80 ⇒ mean + 0.84σ via `scipy.stats.norm.ppf`), expressed as a
multiple of photon-noise σ (not a fraction of flux — 20% of flux would
correspond to ~45σ at typical Keck S/N, physically nonsensical) and
rescaled from the observed (possibly absorbed) flux level to the
FITTED continuum level (`err × sqrt(pred/flux)`) before use, so the
offset reflects true local continuum brightness rather than being
artificially small at line cores. Verified on the synthetic test: bias
moved from ~0% to +0.32-0.73%, matching the predicted
`0.84σ × relative-noise-level` closely. `target_percentile=50`
reproduces the pre-fix (mean-tracking) behavior.

This offset's noise scale is itself empirically calibrated (commit
`0a08008`) rather than trusted at face value — see §9.4 for why, and
§9.4 for the fix.

### 6.5 Base smoothness-penalty recalibration

`lam` (base λ, used wherever 6.1's adaptive stiffening is not
triggered) was never revisited after adaptive stiffening (§6.1) was
introduced to take over resisting broad/blended regions, leaving it
needlessly rigid elsewhere. Confirmed on a real Keck order (5760-5800
Å, a red-order-edge transition with no genuinely broad feature, so
adaptive stiffening correctly does not engage): the achievable local
envelope (real data's peaks) declined measurably faster than `pred`
could track, undershooting true continuum by 3-4%. Forcing MORE
adaptive stiffening across the same zone made this dramatically worse
(envelope collapsed to 0.52-0.93), confirming the problem was base
rigidity, not an under-triggered adaptive mechanism. Scanning `lam`
directly against both synthetic ground-truth tests and the real Keck
envelope found `lam = 2×10³` (recalibrated down from `2×10⁴`, commit
`0a08008`) fixes the Keck envelope (0.96-0.97 → ~1.00) while also
*improving* both synthetic tests' own error metrics independently.

### 6.6 Weight-scale invariance across instruments

`lam`/`p` are absolute numbers, implicitly calibrated against one
particular data scale (Poisson noise on a ~50,000-count continuum, per
the synthetic ground-truth test: err ≈ 224, weight ≈ 1/224²).
Response-corrected MAROON-X flux sits many orders of magnitude away
(millions of counts, weight ~10⁻⁷ to 10⁻¹⁰) — not merely a bigger/
smaller version of the same problem, but a qualitatively different,
badly-conditioned regime once the response-correction error-propagation
fix (§7) made a MAROON-X order's error correctly small in a wide, near-
total-absorption region. `fit_als_continuum()` now rescales `lam` by
the ratio of the order's own typical weight to the calibration
reference (commit `e2be2e2`) — a no-op for Keck/GRACES (ratio ≈ 1,
confirmed bit-identical smoke/synthetic-test output), restoring the
same λ:weight balance regardless of an instrument's absolute flux
units.

### 6.7 Empirical noise-scale calibration

The §6.4 offset (and, implicitly, the base inverse-variance weighting)
assumes `err` correctly predicts real point-to-point scatter
(Poisson-style `err ∝ √flux`). Confirmed FALSE for MAROON-X
specifically: measuring the actual observed scatter in the flattest
short (~50-point) windows of real MAROON-X data and comparing to what
`err` predicted there, real scatter was only 30-40% of theoretical,
while the identical test on Keck data matched almost exactly (~100%).
Interpreted as a consequence of MAROON-X's "optimal extraction"
pipeline, which combines multiple raw CCD pixels per output point via
inverse-variance-weighted PSF fitting — correlating adjacent extracted
points and yielding less independent per-pixel scatter than naive
photon-counting statistics predict. `err`'s relative SHAPE (how it
scales with brightness/response across an order) remains correct; only
its absolute scale was inflated. Since the §6.4 offset is sized
directly from `err`, this caused a real, visible, spectrum-wide
continuum overshoot specifically on MAROON-X (user-identified: "the
observed spectrum is pretty consistently depressed compared to the
[fitted] continuum").

Fixed (commit `0a08008`) by calibrating empirically each time, using
the fit's own above-fit residuals (least likely to be real absorption,
by construction): the ratio of `median(residual | residual > 0)` to
what a half-normal distribution with `err`'s σ predicts for that same
statistic (`0.6745σ`) rescales the offset's noise term before use. This
is self-calibrating rather than an instrument-specific correction
factor — confirmed near 1.0 (no-op) on Keck and the synthetic tests,
0.3-0.6 on MAROON-X, matching the direct empirical measurement.

## 7. Instrument response/blaze correction

`response_correction.py` (introduced in commit `e14e912`, MAROON-X arm
bug fixed in `4e6c5f4`, error propagation fixed in `e2be2e2`) divides a
per-order response/blaze curve into raw flux before continuum fitting,
matching response chunks to science orders by wavelength overlap.

Two real bugs were found and fixed:

1. **Near-zero-response edge amplification**: a response curve tapers
   toward zero at order edges by construction; dividing raw counts (and
   their noise) by a near-zero response amplifies a handful of edge
   pixels far more than their neighbors (confirmed: one real MAROON-X
   order's lowest-throughput edge pixel, response 0.0117% of that
   chunk's peak, amplified 570 raw counts into ~4.9×10⁶ "counts" —
   1000×+ its already-corrected neighbors). Falling back to raw flux at
   these points is worse (raw counts then sit ~1000× BELOW corrected
   neighbors, mimicking a sharp absorption line). Fixed by interpolating
   across below-floor points from neighboring corrected values.
2. **MAROON-X blue/red dichroic-arm ambiguity**: MAROON-X's two arms
   physically overlap in wavelength near their dichroic split (both
   independently cover echelle orders 91-94), so wavelength-overlap-
   only matching can silently pick the wrong arm's response chunk,
   producing a spurious ~4× monotonic trend across an entire science
   order — not just an edge artifact, since the mismatch is in the
   whole chunk's shape. Found from comparing two science orders known
   to cover the same physical line (Hα, in orders 2 and 60): order 2
   looked correct, order 60 did not. Fixed by labeling response chunks
   and science orders by arm (`readers.get_maroonx_bands()`,
   `load_maroonx_response()` now returns a 3-tuple including band
   labels) and preferring same-arm matches. Verified: order 60's median
   normalized flux went from 0.829 to 0.980 (matching order 2's 0.977);
   orders 59/61 corrected the same way; order 48 (a genuinely different
   echelle order, unaffected by the arm ambiguity) was unaffected as
   expected.
3. **Error propagation** (see §6.7's lead-in and §9.1): `obs_err` was
   not scaled by the same response division as `flux`, silently
   dropping the response-division noise-amplification factor and
   biasing the continuum fit's confidence low in low-response
   stretches specifically. This also silently corrupted MAROON-X's EW
   uncertainty estimates (`obs_err` feeds line-fitting error
   propagation elsewhere in `spectrum_data.py`), independent of the
   continuum-shape effect. Fixed by scaling `obs_err` identically to
   `flux` during response correction, and having `normalize()` use that
   maintained value instead of recomputing `sqrt(flux)` from the
   already-corrected flux.

## 8. Cross-order overlap consistency check

`overlap_check.py` (commit `e6d016d`) compares `normalized_flux`
between every pair of orders whose wavelength ranges overlap — a
ground-truth-free consistency check, since two overlapping orders are
independent measurements of the same wavelengths of the same star,
near-simultaneously; if both orders' continuum placement is correct,
they should agree without needing to know the true continuum level.
Runtime is negligible (measured: 5 ms for MAROON-X's 62 orders, a post-
processing step over already-computed `normalized_flux`, not part of
the per-order fitting loop).

Immediately validated on real data: correctly re-identified MAROON-X's
independently-known-bad regions (order 0's low S/N, orders 58-61's
dichroic-arm-boundary wavelength duplication) as the largest
disagreements (2.2-6.2%) without being told anything about them in
advance, and surfaced one new, real, previously-unknown issue: GRACES
order 5 vs. order 6, +5.3% median disagreement in their ~14 Å overlap,
well above the ~0.5-2% baseline residual disagreement seen everywhere
else on all three instruments (see §9.5 for the investigation and
disposition).

Wired into the existing EW-flagging pipeline (§4) rather than building
a parallel mechanism: `Spectrum_Data.flag_order_overlaps
(threshold_pct=2.0)` keeps only pairs disagreeing at or above this
threshold (chosen to sit above the common ~0.5-1.5% baseline while
catching both the MAROON-X and GRACES cases above) and stores them;
`check_for_flags()` then flags any line whose rest wavelength falls
inside a disputed range, without attempting to determine which of the
two disagreeing orders is at fault (often genuinely ambiguous) —
deliberately conservative. Fully opt-in/backward-compatible: nothing
changes for existing pipelines unless `flag_order_overlaps()` is called
first (confirmed: unrelated smoke-test output unchanged).

## 9. Chronological bug log with quantitative verification (real-data phase)

This section is a condensed, dated timeline of the real-data
verification phase (2026-09-14 through 2026-09-15), cross-referencing
§§6-8 above by commit. Verification throughout used two complementary
methods: (a) a synthetic ground-truth test (`test_normalize.py`, not
under version control — a flat continuum plus narrow lines and one
deliberately broad/strong line, with a `--slope` variant reproducing a
real GRACES-like large-scale decline, all with Poisson noise added to a
KNOWN true continuum, so recovery error is directly measurable), and
(b) full-order-set visual review of all real orders for all three
instruments (Keck, GRACES, MAROON-X — 16, 35, and 62 orders
respectively, 113 pages total), not just cherry-picked single orders,
which caught several issues invisible in spot-checks.

1. **(commit `dba6464`)** AsLS replaces the spline+edge-guard patch
   stack (§5-6). Verified on the synthetic test and an initial 3-real-
   spectrum spot check; full 113-page review then caught:
2. **(commit `4e6c5f4`)** MAROON-X arm mismatch (§7.2) and cumulative
   low-`p` bias (§6.3), both described above.
3. **(commit `e2be2e2`)** Response-corrected error propagation (§7.3),
   the weight-scale/λ-conditioning fix it exposed (§6.6, see §9.1), and
   percentile targeting (§6.4).
4. **(commit `0a08008`)** Base λ recalibration (§6.5, a separate
   parameter from §6.6's weight-scale ratio), Signal 2 (§6.1), and
   empirical noise calibration (§6.7).
5. **(commit `e6d016d`)** Cross-order overlap check (§8) and its
   integration into EW flagging (§4, §8).

### 9.1 MAROON-X order 48 vs. order 56 (motivated §7.3, §6.6)

User observation: order 56 has comparable telluric (O₂ B-band)
absorption depth to order 48 (O₂ A-band), but only order 48's *clean*
(non-telluric) region showed a continuum-placement gradient (fitted
continuum 0.80 at one end of an 83 Å clean stretch, rising to 1.02 at
the other). Traced to §7.3's error-propagation bug: order 48's clean
stretch spans response 37%→100%→76% of its chunk's peak (a large
relative-response swing across the stretch), while order 56's telluric
onset happens to sit right at its own chunk's response peak, masking
the identical underlying bug there. Fixing §7.3 alone caused the
previously-working AsLS iteration to diverge catastrophically on order
48 specifically (predicted continuum ~35× actual flux) because the
corrected, now-much-smaller error fell far outside the weight regime
`lam`/`p` were calibrated for, inside a ~26 Å near-total-absorption
telluric trough wider than `wide_window` (so Signal 1, §6.1, under-
detects it — severity capped at 0.27, never approached 1). Two false
starts were tried and rejected before the real fix (§6.6): rescaling
error back toward the old absolute level (defeats the correctness fix);
clipping any single point's weight to within 100× the order's median
(still diverged — the disproportionate-weight cluster was ~150+
contiguous points, not one outlier).

### 9.2 Recalibrating the synthetic slope test

The `--slope` synthetic test's decline rate was originally ~3× steeper
than any real spectrum in hand (5× decline over a 30 Å half-width, vs.
real GRACES's 5× over ~100 Å), making it a falsely pessimistic
regression check mid-investigation. Recalibrated to
`0.55 + 0.45 sin(πx)^1.5` (matching the real rate), after which it
passed cleanly.

### 9.3 MAROON-X orders 1-25 ("sits below the fitted continuum")

Distinguished from a real bug by checking the local 5 Å rolling-MAX
envelope (the actual achievable peak level) rather than the median: the
envelope already sat at ~1.008, comfortably at or above 1.0. The
"sits below the line" visual impression is the expected, correct
appearance of a densely-lined spectrum — most individual pixels carry
real weak absorption, so the bulk/median legitimately sits below 1.0
while genuine continuum-reaching peaks are unaffected. General lesson
recorded for future continuum-quality assessment: judge placement by
the local envelope/peak metric, not the median, which is dominated by
line density rather than fit quality.

### 9.4 MAROON-X, spectrum-wide overshoot

A distinct, more serious follow-up complaint after §6.4 was deployed:
"the observed spectrum is pretty consistently depressed compared to the
[fitted] continuum," specifically and only on MAROON-X. This motivated
and is resolved by §6.7.

### 9.5 GRACES order 6 red edge (found via §8, NOT fixed — see disposition)

Order 5 vs. order 6 overlap disagreement (+5.3%) traced to order 6's
true red edge: raw flux there declines genuinely and steeply (-28.5%
over the last 25 Å, a real echelle blaze rolloff) while simultaneously
containing several real, moderately deep absorption features in the
same stretch — i.e., the region genuinely needs LESS smoothness-penalty
rigidity (to track the fast real decline) and MORE rigidity (Signal 2,
§6.1, correctly measures severity 0.37-0.58 there — real blending, not
a false positive) at the same time. Confirmed directly that this is a
genuine architectural limit, not a mistunable threshold: disabling
adaptive stiffening entirely (§6.1) IMPROVES the envelope there (0.867
→ 0.946), proving the current mechanism (which can only ever add
rigidity, never reduce it below the base value) cannot resolve a
region needing both properties simultaneously. Also confirmed this is
not a fixable detrending-window artifact (shrinking Signal 2's
detrending window from 25 Å to 6 Å does not change its severity
reading).

**Disposition, decided explicitly with the user**: not fixed at the
algorithm level. A real fix would require a new mechanism able to
*reduce* λ in response to a genuinely fast large-scale trend — a change
to the shared fitting routine, requiring the same full re-verification
burden (synthetic tests + all three instruments) as every other change
in §6-7, for the benefit of a single order's ~14 Å edge on one
instrument. Instead, handled via the flagging path (§8): any line
landing in a disputed overlap region is routed to the flagged-lines
file rather than the main linelist, the same treatment already applied
(by convention, not yet by this automatic mechanism) to MAROON-X's
independently-known-bad edge regions.

## 10. Known, deliberately out-of-scope limitations

- **MAROON-X order 48's telluric O₂ A-band** (and by the same
  reasoning, order 56's O₂ B-band, and analogous deep telluric bands
  elsewhere): the fit oscillates within the saturated telluric trough
  itself (not the surrounding clean continuum, which is correctly
  fitted — see §9.1). Explicitly descoped by the user: "I am not
  worried about fitting orders with such strong telluric absorption. I
  would rather not engineer a 'fix' for that, if it could negatively
  impact the general continuum fitting routine."
- **MAROON-X order 0 and orders 58-61**: low S/N (order 0) and
  wavelength-range duplication at the blue/red dichroic-arm boundary
  (58-61). By user instruction, not intended for EW measurement; now
  additionally caught automatically by the overlap check (§8) for
  orders 58-61.
- **GRACES order 6's red edge**: see §9.5 disposition above.
- **MAROON-X order 32, ~180 points at the array start pinned to one
  constant value**: a response-correction edge-bridging artifact
  (constant extrapolation beyond the response curve's valid range, per
  §7.1), predates the current session's changes, cosmetic/low-impact
  only (well outside where any line would be measured). Not yet fixed.
- **Calibration coverage**: all quantitative tuning (§6.5-6.7's
  specific constants in particular) is validated against one synthetic
  ground-truth test plus real data from two stars (a solar Keck
  spectrum, one GRACES target, one MAROON-X M-dwarf-like target). A
  meaningfully different spectral regime (much lower S/N, a hot star
  with very few lines, a cool star dominated by molecular bands) has
  not been tested and could plausibly surface a new failure mode, the
  same way each real bug in §9 was found via a new real spectrum rather
  than by reasoning about the algorithm in the abstract.

## 11. Test infrastructure

- `test_normalize.py` (`/tmp`, not under version control): the
  synthetic ground-truth test referenced throughout §6 and §9 — a flat
  continuum (plus `--slope` variant) with narrow lines, one
  deliberately broad/strong line, and Poisson noise added to a KNOWN
  true continuum level, so recovery error is directly measurable
  in several diagnostic regions (broad-line core/wings, clean
  continuum, a narrow-line control, order edges).
- `test_keck_smoke.py` (`/tmp`): full real-pipeline regression check
  (load → normalize → shift → measure EWs → flag → write) against the
  bundled Keck sample data.
- `make_order_pdfs.py` (`/tmp`): generates one multi-page PDF per
  instrument with every order's normalized spectrum plotted, used for
  full-order-set visual review (§9's "113-page review").
- `Verification/xspect/*.pdf`: the rendered output of the above,
  regenerated after each algorithm change described in §6-8.

## 12. Commit reference

| Commit | Summary |
|---|---|
| `465ea88` | Incorporate XSpect-EW into SPAE |
| `e92e5c7` | Restructure into topic modules |
| `9f439b7` | Universal spectrum-format auto-detection |
| `99c7207` | Fix wavelength-shift grid-search precision |
| `8c3b32c` | RV-based wavelength shift (recommended default) |
| `44a33e6` | Line-position tracking and EW-flagging file |
| `e14e912` | Polynomial continuum fit; response/blaze correction |
| `dba6464` | Replace hard-selection fitting with AsLS |
| `4e6c5f4` | Fix MAROON-X arm mismatch; fix cumulative low-`p` bias |
| `e2be2e2` | Fix response-corrected error propagation; percentile targeting |
| `0a08008` | Recalibrate `lam`; add density-severity signal; empirical noise calibration |
| `e6d016d` | Cross-order overlap check; EW-flagging integration |
