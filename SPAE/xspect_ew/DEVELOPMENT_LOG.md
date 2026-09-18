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

### 9.6 Line identification and the wavelength-shift sign bug

Once the continuum-fitting work in §6-9.5 reached a stable state, work
shifted to the NEXT pipeline stage: locating each linelist entry in the
spectrum (line identification), as a distinct step prior to and
decoupled from EW measurement (which was, at this point, being
developed independently/in parallel by a collaborator) — see §12 for
the identification algorithm itself and §13 for what this work
surfaced in wavelength-shift correction, including a significant,
pre-existing sign bug in `apply_rv_shift()` found via direct,
absolute-wavelength validation of the new identification step (§13.3).

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
- **MAROON-X's much lower detection rate than Keck/GRACES (12/78 vs
  78/78)**: see §13.4-13.6 — investigated via five real, independent
  issues (a coarse-search misidentification bias in the RV estimate; a
  units bug applying §6.7's noise calibration to a new context; a
  missing minimum-prominence check in `identify_line()` that let a
  monotonic-slope artifact register as a false detection, user-caught;
  cross-order arbitration trusting raw significance over position,
  user-caught; a missing hard `position_tolerance` acceptance ceiling
  for uncontested single-order candidates). With all five fixed,
  confirmed the remaining gap is genuine, low S/N for this specific
  exposure (0.29%/0.56%/2.14% median relative photon error for
  Keck/GRACES/MAROON-X respectively -- MAROON-X is ~4-7x noisier), not
  a remaining algorithm shortfall: several of the 78 linelist lines are
  intrinsically weak enough that they are genuinely undetectable at
  MAROON-X's S/N even though they are 10-50σ detections at Keck's.
- **Coarse-wavelength-grid identification precision** (e.g. GRACES):
  see §12's validation discussion of Fe I 6716.222 Å — a discrete
  significance-profile peak can land one grid point away from the true
  (sub-pixel) line center on coarsely-sampled data. Within this
  package's established real-data precision envelope and explicitly
  not addressed (a sub-pixel refinement was proposed and declined by
  the user, since it belongs to EW measurement's own centering, not
  this prior identification step).

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
- `make_line_id_pdf.py` (`/tmp`): one page per linelist entry, showing
  a window around each line with the rest wavelength and the
  identify_lines()-identified center both marked, plus detection
  status/significance/blend flag -- the visual-review counterpart to
  §12's identification work, same role for identification as
  make_order_pdfs.py serves for continuum fitting.
- `Verification/xspect/*.pdf`: the rendered output of the above,
  regenerated after each algorithm change described in §6-8 and §12-13.

## 12. Line identification (before EW measurement)

New module `line_identification.py` locates each linelist entry in a
normalized, wavelength/RV-corrected spectrum as a distinct step
**before** any EW measurement is attempted — deliberately separate from
`line_profile.py`'s existing `get_line_window()`/`measure_ew()`
machinery (used for actual EW measurement, being developed
independently — see this document's introduction), so the two can
evolve without conflicting.

The previous (and still-used, by `get_line_window()`) approach was
"whatever the single lowest flux point within ±0.1 Å happens to be,
call that the line," with no test for whether a real feature is even
present there. `identify_line()` instead runs an explicit detection
test: it converts the local search window to a per-point "how many
local-noise-sigma below continuum" significance profile (using each
point's own already-instrument/response-corrected propagated error —
see §6.7 for why this can't be assumed uniform), lightly smooths it,
and finds candidate local minima via `scipy.signal.find_peaks` with a
minimum-significance threshold (default 3σ). A line with no candidate
clearing that threshold anywhere in the search window is reported as
**not detected** — distinct from, and never silently converted into, a
low-confidence position the way the older approach would. Among
candidates that do clear the threshold, the one nearest the rest
wavelength is preferred over a more significant but more distant one
(scored as significance discounted by squared distance from the rest
wavelength, in units of a real, previously-established position-
precision scale, 0.07 Å); a second, competitive candidate nearby is
recorded as a blend flag. `identify_lines_in_spectrum()` runs this
across a full linelist and every order whose coverage could contain
each line, keeping whichever candidate order gives the highest
significance when more than one covers it (naturally preferring
whichever order's local data/continuum quality is better in an overlap
region, the same property §8's overlap check exploits).

**Known, deliberately unaddressed limitation**: two lines close enough
together that their combined significance profile has no resolvable
valley between them are reported as one, correctly-centered detection
with no blend flag — not wrong, but understates that the region isn't
a single isolated line. Confirmed with a synthetic test (two
comparable-depth lines 0.09 Å apart, each FWHM 0.15 Å): reported as a
single un-blended detection. Flagging this class would need a real
per-instrument "typical single-line width" reference this package does
not establish anywhere yet; left as a known gap rather than an
unjustified absolute threshold.

**Validation** (bundled Keck solar sample, all 3 files combined for
full linelist coverage; real GRACES target; the same 78-line Fe
linelist used throughout this package's history):
- Keck: 78/78 lines detected, 0 blend flags (after the wavelength-
  shift sign fix in §13.3 — see there for the dramatically worse
  numbers beforehand, which is how that bug was found).
- GRACES: 78/78 lines detected, 0 blend flags.
- Position precision (Keck, wide-search residual after RV correction):
  RMS 14.6-26 mA depending on which named-line RV was used (§13.2),
  matching or beating this package's previously-established real-data
  precision benchmark (37 mA RMS, §3).
- A real, spot-checked case (GRACES, Fe I 6716.222 Å, order 12) showed
  a 36.6 mA offset traced to the true line center falling almost
  exactly between two adjacent, nearly-tied-significance sample points
  (18.51 vs 18.32σ) on GRACES's coarser wavelength grid — a genuine
  pixel-grid precision limit, not a misidentification (the profile has
  exactly one real, correctly-resolved peak; a separate, much deeper
  feature ~0.57 Å away has no numerical influence at this search
  radius). Confirmed as within the established real-data precision
  envelope, not a new failure mode. A parabolic/sub-pixel refinement
  around the winning peak (mirroring `combine.py`'s existing
  `parabolic_refine()`) would recover this last bit of precision;
  explicitly deferred by the user, since sub-pixel refinement belongs
  to EW measurement's own centering, not this prior identification
  step ("we just need to identify the line, so the EW fitting routine
  can measure the line strength").
- MAROON-X: 12/78 detected at final, strict settings, far below
  Keck/GRACES — see §13.4-13.6 for the full investigation (five real,
  independent issues found and fixed along the way, including three
  genuine misidentifications/false-negatives the user caught visually).
  Confirmed the remaining gap is real, low S/N for this specific
  exposure (~4-7x worse than Keck/GRACES), not an algorithm shortfall.

## 13. Radial-velocity / wavelength-shift robustness

### 13.1 Motivation

`apply_rv_shift()` (§3) depends on a small, fixed set of named
reference lines (`RV_REFERENCE_LINES`: Ca II H&K, Balmer series, Mg b,
Na D) being both present in an instrument's coverage and individually
well-behaved. Validating the new line-identification step (§12)
surfaced two real, separate problems with this: individual reference
lines can disagree with each other by an amount too large to be
measurement noise, and a spectrum with different wavelength coverage
may have none of them at all. Per explicit user direction, addressed
by keeping the named-line approach as the primary/fast path, but adding
a linelist-based measurement as a fallback (when no named lines are
usable) and cross-check (when they are, comparing the two and
preferring the more robust one on disagreement) — not replacing the
named-line approach outright, since "the overall goal is to identify
the lines, so the code can measure EWs," not to build a general-purpose
RV pipeline.

### 13.2 Linelist-based RV and Balmer-line deprioritization

New `radial_velocity.measure_rv_from_linelist()` reuses
`identify_lines_in_spectrum()` (§12) — the same detection-based
centering used for EW-measurement line identification itself, rather
than introducing a third line-centering method into the package — with
a much wider search radius (default 1.0 Å vs. identify_lines()'s own
0.15 Å default for EW identification), since this runs BEFORE any
wavelength correction and must tolerate however large the spectrum's
real, uncorrected offset is. Velocities from all detected lines are
combined via the same sigma-clipped mean `measure_effective_rv()`
already used, but over dozens of lines instead of 3-4, so the clip is
far more effective.

Applying this to the bundled Keck sample directly motivated a second,
smaller fix: the named-line RV varied by up to ~2.6 km/s across the
three sample files (sunb/sunr/suni.fits) depending on which reference
lines happened to be available, while the linelist-based RV was
consistent to ~0.65 km/s across all three. Tracing this down: whichever
file happened to lean on Balmer lines (sunb.fits: Ca II H&K + H-delta +
H-gamma) showed internal disagreement up to ~11 km/s between individual
reference lines — a real difference in line-formation physics between
H and metal lines (broader, more pressure/NLTE-sensitive Balmer wings
centering less precisely), not noise a 3-4-line sigma clip can reliably
separate from a genuine measurement. `measure_effective_rv()` now tries
metal lines (Ca II/Mg b/Na D) alone first when using the default
reference set, falling back to include Balmer lines only when fewer
than `min_metal_lines` (default 2) metal lines are available. An
explicitly-passed custom `lines` dict bypasses this split entirely.
Verified: the three Keck sample files' final RVs tightened from
(-5.191, -2.628, -3.742) to (-3.338, -2.628, -3.742) km/s, all now
within ~1.2 km/s of the linelist-based estimate, down from up to 2.6 km/s.

`Spectrum_Data.apply_rv_shift()` gained `cross_check=True` (default):
if the named-line RV is unavailable, the linelist RV is used as the
sole estimate; if both are available and disagree by more than
`disagreement_kms` (default 2.0 km/s), the linelist RV is preferred (as
the more robust, larger-N estimate) with a clear warning; otherwise the
(cheaper, already-computed) named-line RV is kept. Fully backward
compatible: if no linelist is loaded (`load_lines()` not yet called),
`cross_check` is a silent no-op and behavior is identical to before
this work.

### 13.3 A significant, pre-existing sign bug in `apply_rv_shift()`

Building the cross-check above required computing an ABSOLUTE
comparison (identified position vs. a real linelist's rest
wavelengths) for the first time — every previous validation of
`apply_rv_shift()` (§3) had been DIFFERENTIAL (one spectrum's RV minus
another's). This absolute comparison immediately surfaced a serious,
pre-existing bug: `apply_rv_shift()` (introduced in an earlier session,
commit `8c3b32c`) applied `shifted_wavelength = wavelength * (1.0 +
rv/c)`, but `measure_line_velocity()`'s velocity convention is `v =
c*(observed-rest)/rest` (standard, positive = redshifted), meaning
`observed = rest*(1+v/c)` is the FORWARD relation from rest to
observed frame. Correcting an observed spectrum back to rest frame
(the actual purpose of this method) needs the INVERSE, `rest =
observed/(1+v/c) ≈ observed*(1-v/c)` for `v << c` — the sign was
backwards.

Confirmed directly and unambiguously on real data (Keck, sunr.fits): a
raw, pre-shift residual (identified position vs. rest wavelength) of
-73.1 mA became **-146.1 mA (doubled, same sign) under the existing
`(1+v/c)` formula**, and **-0.0 mA under the corrected `(1-v/c)`
formula**. This explains why the bug went undetected for as long as it
did: the original validation (§3) compared `apply_rv_shift()`'s
DIFFERENTIAL RV between two independently-shifted spectra against
`estimate_shift()`'s own differential measurement of the same pair,
and a consistent sign error applied to both sides of a difference
partially cancels rather than clearly failing.

Impact, measured directly with the full 78-line Keck+GRACES validation
from §12: applying the (buggy) `(1+v/c)` formula collapsed Keck's
identify_lines() default-window detection rate from 76/78 (with NO
shift applied at all, relying on the search window's own margin to
tolerate the raw uncorrected offset) to 14/78. After the fix, Keck and
GRACES both reach 78/78. Fixed in `apply_rv_shift()`
(`spectrum_data.py`); `estimated_shift`'s sign was flipped to match for
consistency with `wave_shift()`'s convention (`shifted_wavelength =
wavelength + shift`).

**This bug affected every spectrum previously processed with
`apply_rv_shift()`** (documented as the package's "RECOMMENDED
default") prior to this fix. Anything downstream of an
`apply_rv_shift()`-corrected wavelength solution from before this
point should be treated as having a wavelength/RV error of
approximately double the star's true RV-implied shift, not a small
correction.

### 13.4 MAROON-X's residual scatter (investigated and resolved -- see also §13.5)

Applying the same pipeline to the real MAROON-X target: the named-line
RV (Na D-based, -30.995 km/s) disagreed with the linelist RV (1.316 ±
3.258 km/s, 45 lines) by 32 km/s — the cross-check correctly identified
this as untrustworthy and fell back to the linelist estimate. Na D is
itself a well-known problem line for stellar RV work independent of
this package (frequently contaminated by interstellar-medium or
telluric/geocoronal absorption near the line core), consistent with
the named-line estimate being the wrong one here, not the linelist.

However, even using the (then-preferred) linelist RV, MAROON-X's
default-window (0.15 Å) detection rate was only 13/78, far below
Keck/GRACES's 78/78. Investigated and found TWO real, independent,
compounding causes (both now fixed):

1. **Wide-search-radius misidentification biasing the RV estimate
   itself, not just adding scatter.** Printing individual per-line
   velocities from the coarse (1.0 Å) linelist RV pass showed no smooth
   trend with wavelength or order -- including wide swings WITHIN a
   single order (e.g. one order's own lines implying -30, +21, 0, -11,
   and +27 km/s) -- inconsistent with a real RV or calibration drift,
   consistent with the search radius being wide enough to lock onto an
   unrelated real neighboring line in this densely-lined spectrum's
   crowded regions. Confirmed directly: narrowing the search radius
   alone (with nothing else changed) collapsed the scatter smoothly
   (std 21.9/12.0/7.5/4.3/3.5 km/s at 1.0/0.5/0.3/0.2/0.15 Å) while the
   MEDIAN stayed stable near +0.2 to +0.4 km/s throughout -- meaning the
   true RV was already small and well-determined, but the WIDE pass's
   OWN reported value (1.316 km/s, itself computed at 1.0 Å) was
   measurably biased by the same misidentification, not merely noisier.
   Fixed (§13.2's `measure_rv_from_linelist()`): now runs COARSE-then-
   FINE, the way a standard cross-correlation RV search does -- the
   coarse pass's own median velocity is applied as a trial correction,
   then a second, much narrower pass (`refine_radius`, default 0.15 Å)
   refines it, largely immune to the wide pass's own bias since most of
   the true offset is already removed before the narrow, less-
   ambiguous search runs.
2. **A units bug in applying this document's own §6.7 noise-calibration
   fix to a NEW context.** Attempting to apply the same empirical
   noise-recalibration idea from continuum fitting (`err`'s theoretical,
   Poisson-style scale doesn't match real MAROON-X data -- confirmed
   independently, again, here: computed calibration factor median 0.34,
   matching §6.7's continuum-fitting-context value almost exactly)
   directly inside `identify_lines_in_spectrum()` at first produced NO
   improvement, tracked down to a units mismatch: that function receives
   NORMALIZED flux (~1.0 baseline), but the calibration formula
   (borrowed verbatim from continuum.py, where `flux` means RAW counts)
   computed `resid = flux - pred` directly, which is nonsensical when
   `flux` is dimensionless and `pred` is in raw-count units. Fixed using
   the algebraic equivalent for normalized flux, `resid = pred*(flux -
   1.0)` (since `normalized_flux = raw_flux/pred` by definition) --
   confirmed exactly reproduces the correct, unit-consistent calibration
   once fixed.

Combined, these two fixes took MAROON-X's detection rate from 13/78 to
28/78, with a well-converged linelist RV (0.503 ± 0.677 km/s from 28
lines, consistent with the ~0.2-0.5 km/s value found stable across
every search radius during the investigation). Verified no regression
on Keck/GRACES from either fix (both stay at 78/78; Keck's own
linelist RV cross-check moved by <0.04 km/s).

### 13.5 A real misidentification, and a minimum-prominence fix (user-caught)

Visually reviewing the 28-detection MAROON-X result (`identify_lines()`
diagnostic PDF), the user caught a genuine misidentification: Fe I
5661.346 Å (order 17) was reported DETECTED at 3.0σ, but the "center"
sat in the extreme blue wing of a much deeper, unrelated line, with no
real feature of its own. Diagnosed precisely: the winning candidate's
significance profile showed a smooth, monotonic RISE from the window's
start right up to that point (values 0.60, 0.65, 0.59, 0.55, 1.19,
0.83, 0.82, 1.42, 2.24, 3.03), heading into a real, much stronger
feature that the search window (±0.15 Å) cut off before reaching --
the "peak" registered as one only because the very next (and last)
point in the window happened to be marginally lower (3.01 vs 3.03).
Its `scipy.signal.find_peaks` prominence was 0.02 -- essentially zero,
confirming it was not a real local feature at all, just the second-to-
last sample of an unresolved slope, sitting only 0.023 Å from the
search window's own edge.

Fixed: `identify_line()` gained `min_prominence` (default:
`min_significance` itself, so no new unjustified constant), passed to
`find_peaks(..., prominence=min_prominence)`. A real, isolated
detection should stand out from ITS OWN local surroundings by roughly
as much as its absolute height; a point on a monotonic slope into a
different, cut-off feature has near-zero prominence by construction,
regardless of its absolute significance. Confirmed this specific case
is now correctly rejected (`detected: False`).

Effect across all three instruments: Keck unchanged (78/78, 0 blends —
its real detections all have genuine, isolated prominence already).
GRACES lost exactly one previously-marginal detection (77/78) -- Fe I
7114.549 Å, ref EW 8.0 mA (the single weakest line in the entire
linelist), sitting in a noisy, ambiguous stretch with no clear isolated
dip at rest wavelength; a legitimate rejection, not a loss. **MAROON-X
dropped sharply, from 28/78 to 10/78, with all 4 blend flags also
disappearing** -- i.e., most of the fixed round's apparent recovery
(§13.4) was itself low-confidence, edge-driven, or shoulder-of-another-
line detections that a real prominence check correctly discards. The
user anticipated this drop explicitly before it was measured.

**This resolved what "the remaining gap" (§13.4) actually was**: not
primarily the fiber-selection ambiguity or a subtler remaining bug, but
real, quantitatively confirmed low S/N. Measuring each instrument's
typical relative photon error directly: Keck 0.29% (S/N ~349), GRACES
0.56% (S/N ~180), MAROON-X 2.14% (S/N ~47) -- roughly 4-7x worse than
the other two. Since several of the 78 linelist lines are intrinsically
weak (multiple under 20 mA, one as low as 8.0 mA), a line that is a
confident 10-50σ detection at Keck's S/N would only be ~1.5-7σ at
MAROON-X's for the identical true depth -- genuinely too marginal to
trust, not an identification-algorithm shortfall. 10/78 at strict,
uniform significance/prominence standards is the honest, correct
answer for this specific exposure's real data quality, not a remaining
bug to chase further. The fiber-selection ambiguity (`readers.py`)
remains a separate, valid, but now lower-priority open question --
low S/N alone is sufficient to explain the detection count without
invoking it.

### 13.6 Three more real misidentifications (user-caught), and two structural fixes

Visually reviewing the (now 78/78) GRACES PDF, the user flagged four
more questionable identifications. Direct diagnosis of each against the
raw per-order data:

- **Fe I 6392.535 Å** ("near order edge"): the winning candidate came
  from order 14, whose own boundary sits only ~2.4 Å past this
  wavelength. Confirmed order 14's flux there is a continuum-
  normalization artifact, not a real line: it declines monotonically,
  never recovering, all the way to that order's literal last data
  point (0.93 five Å out, down to 0.48 at the very edge). Order 13 (57
  Å interior), covering the same true wavelength via the overlap,
  independently found a correctly-centered, appropriately shallow
  (~6%, 21.6σ) dip essentially exactly at rest wavelength -- the real
  answer, but LOWER raw significance (21.6 vs. order 14's 76.8) than
  the artifact, so the existing highest-significance-wins cross-order
  arbitration in `identify_lines_in_spectrum()` picked the wrong one.
- **Fe I 6745.090 Å** ("strong absorption to the blue"): same root
  cause, different flavor. Order 12 found a correctly-centered,
  appropriately weak (~4%, 8.1σ) dip essentially exactly at rest
  wavelength. Order 11 found a real, well-formed, but 0.033 Å-offset,
  much deeper (~10%, 18.9σ) feature -- too deep for this 8.1 mÅ line,
  more likely a genuinely different absorption nearby. Cross-order
  arbitration again picked the deeper, mispositioned answer.
- **Fe I 7114.549 Å** ("strong absorption surrounding line"): the
  opposite failure -- a real line wrongly REJECTED, not misidentified.
  There IS a genuine ~5-6% dip essentially at rest wavelength (raw
  significance 10.5σ, comfortably above the 3.0 floor), but this
  specific part of order 10 is coarsely sampled (6 points across the
  ±0.15 Å window), and a separate, shallower dip just outside the
  window to the left kept that window's own edge value elevated
  (8.8σ), artificially capping the correctly-centered peak's prominence
  at 1.7 -- under the 3.0 floor -- purely because the narrow window
  never showed the point further out where the profile truly returns
  toward baseline.
- **Fe II 6113.222 Å**: no such entry exists; nearest is Fe I 6113.322
  Å (12.2 mÅ), which the linelist itself annotates `#close to another
  line at base`. The identified center is 0.031 Å off, consistent with
  that documented neighbor pulling the apparent center without
  producing a resolvable second peak -- the blend-without-a-resolvable-
  valley limitation already noted in the module's own docstring (§12),
  not a new bug.

Two of these four are the SAME underlying defect (cross-order
arbitration trusting raw significance with no position or plausibility
check); one is a distinct, opposite defect (prominence miscomputed from
too narrow a window); one is a already-documented, accepted limitation.
Fixed the two real bugs in `line_identification.py`:

1. **`identify_line()` now computes the significance profile (and thus
   prominence) over a WIDER `context_radius` window** (default
   `2*search_radius`) than the window used to accept a candidate's
   POSITION (still `search_radius`, further tightened below).
   Candidates found outside the acceptance window but inside the
   context window inform prominence only, never a reported position.
   This directly fixes the Fe I 7114.549 false negative (its true
   local valley only becomes visible with the wider context) while
   simultaneously making order 14's Fe I 6392.535 artifact WORSE, not
   better disguised -- confirmed it never returns toward baseline even
   ±0.6 Å out, i.e. it's a real, extended, non-line-shaped artifact, not
   just a narrow-window sampling accident.
2. **Cross-order arbitration in `identify_lines_in_spectrum()` now
   prefers the candidate positioned CLOSER to rest_wave, not the one
   with higher raw significance** (ties within 0.005 Å broken by
   significance). Two overlapping orders are two independent
   measurements of the same true spectrum; if they disagree on WHERE
   the line is by more than noise can explain, that disagreement is
   itself the signal that one of them is looking at something else.
   Confirmed this alone flips both the 6392.535 and 6745.090 cases to
   the correct order.
3. A widened context window alone does not stop a single, uncontested
   candidate (no competing order) from winning purely because nothing
   else was in its window. Checking the offset distribution directly
   confirmed this WAS a real, separate problem specifically for
   MAROON-X: every confirmed-good Keck/GRACES detection has
   `|offset| <= 0.058` Å, but before this fix MAROON-X had several
   accepted detections at 0.08-0.14 Å with statistically absurd
   significance for their catalogued strength (e.g. 115σ for a 10.5 mÅ
   line -- clearly a different, real, much stronger feature, not the
   catalogued one). Fixed: the final accepted candidate must now fall
   within `min(search_radius, position_tolerance)` (0.07 Å default),
   not the full, wider `search_radius` (0.15 Å) -- `search_radius`
   still governs how far the search for CANDIDATES extends, but
   `position_tolerance`, this package's own measured RV-corrected
   real-data precision, now gates final acceptance. Zero effect on
   Keck/GRACES (neither ever had an offset past 0.058 Å to begin with).

Effect across all three instruments, verified directly: **Keck 78/78,
GRACES 78/78** (both 0 blends -- GRACES gained back Fe I 7114.549 via
fix 1, its only previous loss). **MAROON-X: 12/78** (0 blends) -- up
from 10/78 (net +2 genuine recoveries, both near-zero offset: Fe I
6726.666 and Fe I 7114.549), after an intermediate check (fix 1 alone,
before fix 3) had spiked it to a false-looking 22/78 that the offset-
distribution check above caught and fix 3 corrected back down. The
original user-caught misidentification from §13.5 (Fe I 5661.346)
remains correctly rejected throughout. All final MAROON-X detections
now have plausible significance-to-catalogued-strength ratios and
`|offset| <= 0.067` Å.

**One remaining MAROON-X case flagged for later, not fixed here**: Fe I
5546.991 Å (order 19, 6.0σ, offset +0.034 Å) looks like a real, clean,
isolated dip -- returns fully to baseline on both sides, no cross-order
disagreement, no order-edge or window-truncation signature, and the
AsLS continuum fit through this region is smooth with no ringing near
the huge unrelated line 0.5 Å redward. But it's shallower (~5-6%) than
naively expected for a 27.7 mA line at this resolution, and the user
noted the region is noisy enough that a cosmic ray may have hit the
center pixel, distorting the profile. Deliberately left as detected
(not enough here to justify rejecting a statistically legitimate,
isolated 6σ candidate): the user's assessment is that this class of
problem -- a single bad pixel inside an otherwise real line -- is
better addressed with an explicit cosmic-ray/outlier check at the EW
MEASUREMENT step, not by making identification itself more skeptical
near the position-tolerance boundary. Left as a known open item for the
EW-measurement phase (independently developed -- see §12's decoupling
note), not a bug in this module.

## 14. Line-identification S/N floor (synthetic sweep)

Motivation: MAROON-X's much lower real detection rate (§13.4-13.6)
raised a broader, forward-looking question -- what is this package's
own FUNCTIONAL S/N limit for line identification, independent of any
one real spectrum's particular quirks? Split deliberately into two
separate studies matching this module's own decoupling from EW
measurement (§12): this section characterizes IDENTIFICATION only; an
EW-measurement-stage characterization is deferred until that
(independently developed) code is available.

**Method**: start from real Keck data (`sunr.fits`, 53/78 linelist
lines covered), not a fully synthetic spectrum, so real line
shapes/blends/continuum shape are preserved and only the noise is
synthetic. For each target relative photon error, add zero-mean
Gaussian noise to the RAW (pre-normalization) counts with variance
`max(0, target^2*flux^2 - flux)` -- i.e. exactly enough, in quadrature
with the existing native Poisson noise, to bring the TOTAL relative
error to the target value -- and override `obs_err` to the same target
exactly (`target * flux`), rather than re-deriving it from the noisy
flux. This isolates "how does identification degrade with true S/N"
from "does the code correctly know its own noise level", which was
already investigated separately (§6.7, §13.4's units bug). The full
realistic pipeline (`normalize_all -> load_lines -> apply_rv_shift ->
identify_lines`) is re-run at each level, so continuum-fitting quality
degrades too, matching how a genuinely fainter real exposure behaves.
5 independent noise seeds per level; detection rate measured against
only the 53 lines this file can ever cover (not all 78).

**Result** (detection rate = mean +/- std over 5 seeds; see
`Verification/xspect/sn_sweep_identification.png`):

| target rel. error | S/N (~1/err) | detection rate | med \|offset\| |
|---|---|---|---|
| 0.29% (Keck native) | 349 | 100.0% | 22.1 mA |
| 0.56% (GRACES-like) | 180 | 100.0% | 20.9 mA |
| 2.14% (MAROON-X's old, order-median figure) | 47 | 99.6 +/- 0.8% | 21.7 mA |
| 3% | 33 | 96.2 +/- 2.4% | 21.8 mA |
| 4% | 25 | 93.2 +/- 0.9% | 22.8 mA |
| 6% | 17 | 87.2 +/- 1.4% | 22.8 mA |
| 8% | 12.5 | 81.5 +/- 2.2% | 23.2 mA |
| 10% | 10 | 78.1 +/- 3.3% | 24.3 mA |
| 15% | 6.7 | 69.1 +/- 2.6% | 26.5 mA |
| 20% | 5 | 59.2 +/- 1.9% | 28.5 mA |
| 25% | 4 | 50.6 +/- 6.6% | 28.7 mA |
| 30% | 3.3 | 43.0 +/- 6.5% | 25.9 mA |
| 40% | 2.5 | 37.0 +/- 4.4% | 24.9 mA |
| 50% | 2 | 35.5 +/- 6.2% | 34.7 mA |

**Conclusion**: the identification ALGORITHM's own degradation is
smooth and graceful, not a cliff -- detection rate is still ~99.6% at a
relative error matching MAROON-X's real (old, order-median) noise
figure, and only crosses 50% around a ~25% relative error (S/N ~4),
flattening to an asymptotic floor around 35% (the population of very
strong lines that survive almost any noise level). This threshold is
more than an order of magnitude worse than the noisiest real instrument
used anywhere in this project (MAROON-X's real per-line local relative
error is ~0.59%, confirmed in §13.6 -- essentially indistinguishable
from GRACES's 0.56%). **This confirms directly, not just by inference,
that MAROON-X's real 12/78 detection rate is NOT an identification-
algorithm S/N limitation**: a spectrum with MAROON-X's actual per-line
noise level, if it had Keck-like (narrow, un-broadened) line profiles,
would be expected to detect ~99.6% of lines, not ~15%. Combined with
§13.6's rotational-broadening hypothesis (this is a moving-group target
likely to be a faster rotator than the Sun even with otherwise
solar-like Teff/logg/[Fe/H], per the user), the remaining MAROON-X gap
is best attributed to a genuine property of that specific star, not a
package limitation -- consistent with the user's decision to keep this
spectrum in the analysis sample rather than exclude it, since its
independently-derived stellar parameters will allow a direct check once
the full SPAE pipeline (EW measurement through abundances) is run on it.

A parallel EW-measurement-stage S/N characterization is deferred until
that (separately developed) code is available, per the user's original
two-step framing of this question.

## 15. New GRACES format, a negative-flux bug, and whole-order outlier flagging

Motivated by testing/improving `combine_spectra()` (multi-exposure
co-addition) -- the user provided two real GRACES exposures of the same
star (Theia 456-6, `N20220115G0040i.fits`/`N20220118G0038i.fits`) as test
data, which turned out to use a GRACES/OPERA output variant this package
didn't support yet.

### 15.1 New reader: GRACES/OPERA merged "intensity" (.i) format

Different from the already-supported per-order "spectrum" (.m) format
(§7's original GRACES reader): no `Order` column at all -- every echelle
order pre-merged into one continuous ~194000-point array. The header is
fully self-documenting (verified via its own source, not assumed): 4
wavelength/intensity/errorbar column triplets (`COL1-3`, `4-6`, `7-9`,
`10-12`), Normalized/UnNormalized crossed with autowave-corrected
(telluric-line-based)/uncorrected, distinguished only by each `COLn`'s
header COMMENT (the VALUES are identical repeated 'Wavelength'/
'Intensity'/'ErrorBar' labels). Confirmed directly against the OPERA
pipeline's own C++ source (`operaGenerateLEFormats.cpp`,
github.com/CFHT/OPERA): "UnNormalized" = `RawFluxInElectronsPerElement`
(a literal copy of the raw extracted flux, no background subtraction),
"no autowave correction" = `ThArCalibratedInNM` (comparison-lamp
wavelength solution only, no telluric-line fine correction) -- both the
least-processed options, matching this module's established "prefer
raw/uncorrected" design principle exactly the way GRACES's `.m` format's
`RawFlux` column did.

New `_read_graces_intensity()`/`_detect_graces_intensity()` in
`readers.py` (registered as `'graces_intensity'`). Since there's no Order
column, orders are recovered by detecting wavelength BOUNDARIES -- a
less robust heuristic than the `.m` reader's explicit column, used only
because this format gives us nothing better. A first attempt (negative
jumps only, i.e. order OVERLAP) silently merged the reddest ~5 orders
into one 39101-point, 2270-A-wide "order" at ~22x coarser effective
sampling than everywhere else, because real order overlap stops past
~820 nm and adjacent orders there have a small wavelength GAP instead
(a forward jump, not a negative one). Fixed: flag a boundary at either a
negative jump OR an abnormally large positive jump relative to a ROLLING
local median spacing (needed either way, since real point spacing itself
grows gradually and smoothly from blue to red). Confirmed against both
real files: recovers exactly 35 orders, matching this instrument's
known-good order count from the `.m` format, with smoothly increasing
order sizes throughout. No regression confirmed on the existing `.m`
GRACES file, Keck, or MAROON-X.

### 15.2 Negative raw flux crashing an entire order's continuum fit

GRACES's real raw ("UnNormalized") flux can be genuinely negative --
confirmed via the OPERA source above (ordinary CCD read noise
symmetrically scattering a near-zero-signal pixel below zero after
bias/dark subtraction, not a data error -- more common here simply
because these two exposures are fainter, SNR ~20-70 per the header, than
anything tested against so far). `Spectrum_Data.__init__`'s initial
`obs_err = sqrt(flux)` produced NaN for those points, and even a modest
fraction of NaN weights (539/4055, ~13%, in one real order) was confirmed
to corrupt the ENTIRE order's AsLS continuum fit output, not just the
affected points. Fixed: `sqrt(abs(flux))` -- keeps the correct order of
magnitude of Poisson-like noise at that pixel either way, and is a
general fix (not format-specific), since any sufficiently faint spectrum
in any format could hit the same failure mode.

### 15.3 Whole-order outlier flagging (new `outlier_check.py`)

Fixing 15.2 let `normalize_all()` complete without crashing, but order 0
of both real files still came back with an absurd median normalized flux
of ~0.055 (should be ~0.9-1.0). Root cause: a single raw-flux point at
~28,000-60,000x the order's own typical scale (154.16 million vs. a
~5,475 median) was dragging the AsLS fit up across the WHOLE order via
its smoothness-penalty coupling. Confirmed NOT a random cosmic ray: the
same defect appears at the same wavelength (~400.56 nm) and nearly the
same pixel index (540/4055 vs. 540/4054) in BOTH independent exposures,
taken 3 days apart -- a fixed, reproducible detector/pipeline defect
(plausibly related to GRACES's dual-amplifier readout boundary), not
noise. Per the user's explicit preference and this package's established
flag-and-exclude precedent (§8's overlap check) rather than trying to
make the shared AsLS routine itself robust to every possible corruption:
new `outlier_check.py` / `Spectrum_Data.flag_bad_orders()`, wired into
`check_for_flags()` exactly like `overlap_flag_ranges` (new
`self.bad_order_ranges`, same `wave_lo`/`wave_hi` key convention, same
silent-no-op-until-called behavior).

**Getting the detection metric right took one real iteration**: a first
attempt (median absolute deviation, symmetric around the order's median)
false-flagged 11 of 35 real orders in the first test file alone, whose
only "outlier" was a genuine, deep, real spectral feature (Na D at
5889.95 A, H-alpha at 6562.63 A, the O2 B telluric band at 6869.59 A)
sitting in an otherwise very low-scatter stretch. Root cause: deep
ABSORPTION already has purpose-built, asymmetric handling in the
continuum-fitting routine itself (§6) -- this check only needs to catch
the failure mode that mechanism can't, an extreme point ABOVE the
order's own natural peak level, so checking symmetrically was both
unnecessary and actively wrong. Fixed: compare the order's single
highest point against its own 99th-percentile flux (max/p99), not a
median-centered measure. Confirmed directly against every real spectrum
in this project's test set (Keck, both GRACES formats, MAROON-X
pre-/post-response-correction): every genuine order sits at max/p99 <=
1.51, while the real defect sits at 96.7-713 in the two test files --
a >60x gap between the two clusters, with the default `outlier_factor`
(20.0) sitting at that gap's log-space midpoint.

## 16. `combine_spectra()` redesign (multi-exposure co-addition)

User request: two or more exposures of the same star are often co-added
for higher combined S/N; XSpect-EW already had a `combine_spectra()` for
this, asked to make it more robust. Test data: the two real GRACES
exposures of Theia 456-6 from §15 (`N20220115G0040i.fits`/
`N20220118G0038i.fits`), the only reason that section's reader/NaN/
outlier work happened at all.

**Problems found in the original implementation** (code review, before
any real-data testing):
1. Orders matched by nearest *median* wavelength only, no check that A
   and B's candidate orders actually overlap.
2. Wavelength alignment reused `estimate_shift()`/`clean_shift()` --
   built for cross-correlating a star against a REFERENCE spectrum (e.g.
   a solar atlas), with a documented extrapolation failure mode for
   orders with no reference overlap (§3/§13). Wrong tool for this case:
   A and B are two exposures of the *same* star, so the natural, already-
   more-robust approach is to rest-frame each one independently via
   `apply_rv_shift()` (§13's work) rather than cross-correlate them
   against each other.
3. The actual flux combination did not interpolate at all: for each of
   A's pixels, it grabbed B's single nearest sample within a hardcoded
   +/-5-pixel INDEX window around the SAME pixel index -- silently
   assuming A and B share an identical pixel grid/length per order, not
   a safe assumption for two independently-reduced exposures.
4. Error was recomputed as `sqrt(combined counts)`, discarding each
   exposure's own already-computed `obs_err` rather than propagating it.
5. A blocking `plt.show()` per order made it impossible to run non-
   interactively.

**Redesign** (`Spectrum_Data.combine_spectra()`, `spectrum_data.py`):
rest-frames both spectra independently via `apply_rv_shift()`; matches
each of A's orders to whichever B order has the MOST real wavelength-
range overlap (minimum overlap fraction required, default 0.5, else that
order is left as A alone rather than paired with something unrelated);
properly interpolates B's flux/error onto A's rest-frame grid
(`scipy.interpolate.interp1d`) over the real overlap only; sums flux and
propagates error as a quadrature sum of each spectrum's own (interpolated)
`obs_err`; skips combining any order flagged by `flag_bad_orders()` (§15)
in EITHER spectrum -- combining a good order with a known-corrupted one
would spread the corruption, not average it out -- falling back to that
spectrum's own order alone instead, same flag-and-fall-back philosophy as
everywhere else rather than trying to salvage a known-bad order; plotting
is now opt-in (`plot=True`) and non-blocking (one PNG per combined order,
not a per-order `plt.show()`). `update_combined()` now also restores
`self.obs_err` (previously only `self.flux`), matching that the combined
error is now a real, propagated quantity worth keeping.

**Verified end-to-end on the two real Theia 456-6 exposures** (both
already `flag_bad_orders()`'d and `normalize_all()`'d first, per the new
docstring's stated requirement -- `apply_rv_shift()` needs
`normalized_flux` to measure each spectrum's own RV):
- Each spectrum's own RV measured independently and reasonably (-3.97
  and -4.45 km/s, from the same 4 named reference lines in both) --
  confirms two independent measurements of the same real target agree at
  the ~0.5 km/s level, consistent with real per-line scatter already
  characterized in §13.
- Order 0 (flagged bad in BOTH exposures, §15) correctly skipped and left
  as A's own data, not combined -- combined count 34/35, matching "1
  flagged bad in A, 1 flagged bad in B" (the same order in both, as
  expected: it's a fixed detector defect, not per-exposure noise).
- All 34 other orders matched 1:1 (both files share the same GRACES order
  structure) at overlap=1.00 and combined successfully -- finite,
  sane flux/error everywhere (no NaN/inf), median relative error on a
  representative order improved from 1.87% (A alone) to 1.25% (combined)
  -- a 1.51x improvement, close to (and, since the two exposures are not
  perfectly equal quality, reasonably above) the sqrt(2)~1.41x expected
  for two equal-quality co-added exposures.
- Fallback path confirmed by deliberately displacing one B order out of
  range: that order correctly falls back to "using A alone" (best
  overlap 0.26 < the 0.5 minimum), combined count correctly drops to
  33/35.
- `plot=True` confirmed to write one non-blocking PNG per combined order
  (34 files, no `plt.show()` blocking); visually, A and B track each
  other closely (same real absorption features, well RV-aligned) and
  A+B shows the expected ~2x continuum level with all line depths
  preserved.

## 17. Precise alignment and cosmic-ray/bad-pixel handling for `combine_spectra()`

After reviewing the §16 before/after plot, the user raised two requirements
that §16 did not yet address: (1) "ensuring the spectra are perfectly
aligned before adding is of paramount importance" -- a residual mismatch
between A and B injects noise/smearing into the co-add; (2) a general way
to detect AND CORRECT cosmic rays/bad pixels -- both when a second
exposure is available for cross-checking, and when only a SINGLE
spectrum exists (explicit requirement: "not every star we analyze will
have more than one spectrum"); and (3) any line whose measurement window
contains a corrected pixel should be excluded from the linelist
entirely, not just flagged.

### 17.1 Alignment was NOT actually verified before this

Direct answer to the user's question: no. §16's independent
`apply_rv_shift()` on each spectrum gets each one close to its OWN rest
frame, but nothing verified the two ended up mutually aligned. Confirmed
this matters: two independent RV measurements of the same real star
disagreed by ~0.5 km/s (§16, from real per-line measurement noise alone)
-- a real several-mA residual mismatch at optical wavelengths.

New `combine.measure_order_alignment()`: a fine, direct cross-correlation
between A and B's own overlapping NORMALIZED flux (not raw -- removes
each exposure's own absolute throughput level), reusing the existing
`parabolic_refine()` sub-pixel machinery. Distinct from (and more precise
than) trusting either spectrum's own independent RV, since it compares
the two spectra directly over many points at once rather than inheriting
either one's own per-line noise. Measured directly on the real Theia
456-6 pair: residuals of 5-29 mA across different orders (order 15's own
cross-correlation chi^2 improved 8x, 0.00524->0.00065, confirming this is
a real, measurable improvement, not fitting noise). `combine_spectra()`
now measures this per matched order pair and applies it to B's wavelength
query before interpolating (`align=True` by default; warns above
`align_warn_threshold`, default 0.05 A, comfortably above every residual
confirmed normal on real data).

### 17.2 Single-spectrum spike detection: `detect_spikes()`/`correct_spikes()`

Needed to work with NO second exposure, per the user's explicit
requirement. Two sigma-based designs were tried and rejected on real
data before landing on a factor-based one -- both failure modes are
documented in `outlier_check.py`'s docstring in detail, summarized here:
(1) comparing each point to a local window's median-absolute-deviation
badly UNDER-estimated the true noise inside any real large-scale slope,
since a short window straddling real curvature has an artificially small
point-to-point MAD -- false-flagged ordinary, smooth Keck continuum
points at up to ~250x their true significance. (2) Comparing to the
spectrum's own theoretical Poisson error (`obs_err`) fared even worse:
real, genuinely RESOLVED spectral structure (a real continuum peak
between two blended lines, or a real sky-emission line) routinely varies
pixel-to-pixel by many sigma of pure photon noise over just a few
pixels -- confirmed directly on real Keck data (a smooth, symmetric,
several-pixel-wide bump, visually indistinguishable in shape from a real
line) and on the real GRACES 5577 A sky-emission-line case (a genuinely
resolved feature, not a defect) -- neither test can tell a modest real
feature from a modest defect using shape/statistics alone.

A simple **ratio to a tight (5-point) local median** can, because on
every real spectrum checked (Keck, both GRACES formats, MAROON-X, the
real GRACES sky-emission line) this ratio never exceeds ~1.6x, while the
real defect motivating `check_bad_orders()` (§15) reaches 5.5-126x at
the same window size -- confirmed via full-spectrum regression, zero
false positives anywhere in this project's validated test set.
`detect_spikes()` (`outlier_check.py`) implements this; deliberately
only catches EXTREME, unambiguous cases (the user's own framing:
"strong spikes") -- a moderate, ambiguous excursion (like the real
5577 A line, ratio ~1.2-1.25x) is deliberately left for the two-exposure
comparison below, which can resolve the ambiguity with a second,
independent measurement instead of guessing from shape alone.
`Spectrum_Data.correct_spikes()` runs this per order, REPLACES (not just
flags) each affected point with its local median, recomputes `obs_err`
there, and records the wavelength in the new
`self.corrected_pixel_wavelengths` (same silent-no-op-until-called
convention as every other `flag_*()`/`correct_*()` method; `skip_orders=`
lets a caller exclude orders already handled by `flag_bad_orders()`).

### 17.3 Cross-exposure spike detection: `check_cross_exposure_spikes()`

When a second exposure IS available, far more sensitive detection is
possible: compare what SHOULD be (for real absorption) two independent
measurements of the same true signal directly. Two designs were tried
and rejected before the final one, both documented in detail in
`outlier_check.py`:

1. A raw-flux difference test (`flux_A - flux_B_interp`, scaled by
   propagated error) flagged nearly every point in the bluest orders --
   root cause: A and B have genuinely different absolute continuum
   levels (different nights' throughput), which a raw difference doesn't
   account for at all.
2. Switching to a NORMALIZED-flux difference (removes the throughput
   issue) still over-flagged low-S/N blue orders by hundreds of points --
   root cause: naive theoretical error underestimates true noise there
   by up to 3x (the same effect `_empirical_noise_calibration()`, §12,
   already exists to correct) -- and even after applying that
   calibration, still over-flagged red orders by dozens of points each,
   this time because two independently-noisy real spectra routinely
   disagree by many sigma of pure photon noise at deep, narrow line
   CORES from ordinary residual sub-pixel sampling differences (confirmed
   directly: order 29's flagged points were all real line cores, with a
   consistent, one-directional bias -- A always deeper than B at dozens
   of unrelated lines, the signature of a real alignment/sampling effect,
   not contamination).

The fix: stop comparing the two spectra's DIFFERENCE at all. Instead,
test each spectrum's own normalized flux against the CONTINUUM (1.0)
independently -- a real absorption line can only ever push flux DOWN, so
"significantly above 1.0" (using each spectrum's own calibrated relative
error) is an unambiguous, physically-motivated excess signature
regardless of local line density or alignment noise. Confirmed on real
data: the known line-core mismatches show excess significance NEGATIVE
(below continuum) in BOTH spectra and are correctly never flagged, while
the real 5577 A and 6300 A [OI] auroral sky-emission lines show a clear,
positive excess.

An asymmetry requirement (flag only if one side is significant AND the
other is not) was tried first, motivated by exactly this real-vs-defect
distinction, but this LEFT THE WORST PART of real contamination
uncorrected: a real sky line present in both exposures at different
strength (confirmed on real data, order 27) shows comparably significant
excess in BOTH spectra right at its peak, only becoming "asymmetric"
toward its weaker edges -- the asymmetry gate correctly avoided real
line-core noise but also skipped the CENTER of genuine contamination,
catching only its edges (confirmed directly: the real 5577 A peak,
~13,300 combined counts, was untouched by the asymmetry-gated version,
only the flanking points at ~3,100-3,800 counts were corrected). Fixed:
flag whenever EITHER side's excess significance clears the threshold (no
asymmetry gate), always using the LOWER (less-contaminated) side's value
-- safe because real absorption can never trigger a positive-excess test
in either spectrum, so no real protection is lost, while the full extent
of real contamination (not just its edges) is now corrected. Verified
directly: order 15's combined peak dropped from ~13,300 (uncorrected sum)
to a sane 7,156 counts (matching the higher of the two individual
exposures' own contaminated values, not their sum) after this fix, with
zero new false positives in the order-29 line-core region.

Wired into `combine_spectra()` as `cross_exposure_check=True` (default),
using the SAME per-order alignment residual from 17.1 so both checks
operate on identically-defined points. A flagged point uses only the
less-contaminated side's raw flux/error (not the sum) and is appended to
`self.corrected_pixel_wavelengths`.

### 17.4 Line exclusion for corrected pixels

`check_for_flags()` extended with a new check: excludes (not just flags)
any line whose rest wavelength falls within 0.1 A (matching
`get_line_window()`'s own default search radius) of any wavelength in
`self.corrected_pixel_wavelengths` -- populated by both
`correct_spikes()` and `combine_spectra()`'s cross-exposure correction.
Deliberately EXCLUDES rather than just flags, since the measurement at
that point reflects a corrected, not originally observed, value.

### 17.5 Final verification on the real Theia 456-6 pair

- 693 total points cross-exposure-corrected across 35 orders (up from
  103 with the rejected asymmetry-gated version) -- concentrated
  overwhelmingly in the red orders (dozens to ~100 per order from order
  25 onward), consistent with the well-known, dense real population of
  OH airglow emission lines longward of ~7000 A; minimal in blue/mid
  orders (0-4 points each) where only the two confirmed [OI] auroral
  lines and a few similar features appear.
- All combined flux/error values finite (no NaN/inf) across every order;
  `update_combined()` correctly restores both `self.flux` and
  `self.obs_err`.
- `check_for_flags()` integration confirmed: a synthetic line placed
  within 0.1 A of a real corrected pixel is excluded with the correct
  reason string; a clean line elsewhere is untouched.

## 18. Interactive widget (`interactive.py`)

Prior to this work, running xspect-ew meant driving `Spectrum_Data`
directly from a Jupyter notebook -- no way to visually confirm each
pipeline stage before moving to the next. Added `EWWidget`
(`interactive.py`), a 4-stage step wizard -- Load, Normalize, RV Shift,
Measure EW -- built entirely from `matplotlib.widgets`
(Slider/Button/TextBox), mirroring the architecture already established
in `SPAE/moog/interactive.py`'s `SynthWidget`: manual `fig.add_axes()`
figure-fraction layout, no separate GUI framework, so the identical code
opens as a standalone window from a terminal script or embeds in a
Jupyter cell via `%matplotlib widget`. Launch via `ew_interactive(...)`.

Every stage calls the same `Spectrum_Data` methods an automated script
would call -- `normalize()`/`normalize_all()`, `apply_rv_shift()`,
`measure_ew()` -- so the widget is a thin, optional layer on the one
shared pipeline API, not a second code path:

- **Load**: confirms order count, wavelength coverage, and linelist size;
  plots every order's raw flux for a first visual sanity check.
- **Normalize**: per-order continuum overlay (raw flux + AsLS fit) with
  live `lam`/`p` sliders (`normalize()` re-fits just the current order on
  each change) and an "Apply to ALL orders" button (`normalize_all()`).
- **RV Shift**: runs `apply_rv_shift()` on first entry, displays the
  measured RV (with its named-line/linelist-cross-check warning text, if
  any -- see §13), and lets the user toggle the current order's display
  between observed and rest-frame wavelength, remeasure, or type in a
  manual RV override.
- **Measure EW**: pages through loaded lines; for each, locates the
  covering order and calls `measure_ew(..., axes=(fit_view, data_view))`
  to draw directly into the widget's persistent axes. Sliders expose
  `measure_ew()`'s own `ex_params` (continuum shift, left/right bound,
  center) for live manual adjustment.

**Refactor required first**: `measure_ew()`'s ~50-line inline plotting
block was extracted into a new pure-drawing function, `plot_ew_fit()`
(`plotting.py`), which either creates a new figure (`axes=None`,
preserving `measure_ew()`'s original `plot=True` behavior exactly) or
clears and redraws into a caller-supplied `(fit_view, data_view)` axes
pair (the new capability the widget needs). This was done deliberately
as a black-box wrapper -- the widget's EW stage does not depend on any
of `measure_ew()`'s internal fitting logic, only on the `axes=` parameter
and the array shapes it already produces, since that fitting logic is
expected to keep changing on the colleague's own branch.

**Degenerate-fit guard**: some `ex_params` slider combinations (e.g. a
center/bound shift narrow enough to leave too few points in the fit
window) make `get_line_window()`'s downstream `np.gradient()` call raise
outright (`ValueError: Shape of array too small...`), confirmed by
smoke-testing every stage against `sunb.fits`. Since a slider drag must
never crash the whole session, `_remeasure_ew()` catches any exception
from `measure_ew()` and reports it in-place in the fit axes instead;
resetting the sliders recovers cleanly (verified: EW measurement for Fe I
4779.439/order 22 returns to ~40.3 +/- 0.5 mA, matching the established
regression value, both before triggering the guard and after).

Regression-tested all three `measure_ew()` call paths (`plot=False`
no-op, `plot=True` original standalone-figure behavior, `axes=(...)`
embedded behavior) against that same Fe I 4779.439/order 22/`sunb.fits`
case -- all three agree on EW within GP-sampling noise.

## 19. Merging `fix_ew_fitting`: new global/local-continuum `measure_ew()`

Merged a colleague's parallel branch (`fix_ew_fitting`) into
`pymoog-conversion` after committing the interactive-widget work above.
Their rewrite replaced the old GP-based (`SEKernel`/`Pred_GP`) EW fit
entirely with a direct weighted Gaussian fit (`gfit_direct`) plus an
independent reference-atlas cross-check (`reference_atlas.py`,
`load_reference_atlas()`), and introduced two parallel continuum
treatments per line: a **global-continuum fit** (assumes the existing
normalization already puts this window's continuum at 1.0 -- this is
what `self.lines_ew` actually reports) and an optional **local-continuum
diagnostic fit** (`fit_continuum=True`: estimate a flat local continuum
from the line's own wing data and refit against that instead, kept only
as a comparison value in `lines_ew_local`). `measure_ew()`'s plotting
changed to match: two side-by-side panels, "Global continuum
(REPORTED)" and "Local continuum (diagnostic only)".

`spectrum_data.py` had 7 conflicting blocks against the interactive
widget's own recent changes (imports, `obs_err`/`pred_all` init,
`combine_spectra()`'s redesigned body, `measure_ew()`'s signature and
its plotting tail). Resolved by taking main's side for the EW-fitting
rewrite itself (imports, init bugfixes, `measure_ew()`, the new
`load_reference_atlas()` method) while keeping this branch's
`combine_spectra()` redesign (§16-17) and `bad_order_ranges`/
`corrected_pixel_wavelengths` state intact -- confirmed via `grep` that
main's own `combine_spectra()` on that side was the pre-redesign
nearest-neighbor version this branch had already superseded, not a
competing rewrite.

**`axes=` re-integration**: the merge dropped `measure_ew()`'s `axes=`
parameter (and `plotting.py`'s `plot_ew_fit()` helper) that
`interactive.EWWidget`'s EW stage (§18) depends on to draw into its own
persistent figure rather than popping up a new one every remeasurement.
Rather than resolve this inside the merge itself, took main's plotting
code as-is first to get a clean, working `measure_ew()`, then reapplied
the same black-box wrapper strategy as §18 as a separate follow-up:
`plot_ew_fit()` was rewritten to draw the new Global/Local two-panel
layout, still either into a new figure (`axes=None`, unchanged
`plot=True` behavior) or into a caller-supplied `(fit_view, local_view)`
pair (cleared first) -- folding the order/flagged title into the left
panel when embedded, since an embedded figure has no `fig.suptitle` of
its own to use. `measure_ew()` regained `axes=None`, forwarded through
the `auto_widen` retry call; plotting now runs whenever `plot=True` OR
`axes is not None`, and `show_plot`/`plt.show()`/`plt.close()` are
skipped whenever `axes` is supplied, matching §18's original contract.
Smoke-tested `plot_ew_fit()` directly (synthetic line profile, `Agg`
backend) in all three combinations (`axes=None`/flagged,
`axes=(ax1,ax2)`/`fit_continuum=False` placeholder panel,
`axes=(ax1,ax2)`/flagged) -- all draw without error.

**Full widget regression, `test_ew_widget_sunb.py` (`/tmp`, not under
version control)**: end-to-end drive of `EWWidget` against the bundled
`sunb.fits` + `Sun_fe_sample.txt` (`Agg` backend, no real display) --
Stage 1's `normalize_all()`/per-order refit/Apply-to-ALL, Stage 2's
`apply_rv_shift()` (correctly preferred the cross-correlation RV over
the disagreeing named-line RV here) and manual override, and Stage 3
stepping through all 78 linelist lines (75 fall outside this
spectrum's 3643-4795 A coverage, exercising the "no order covers this
line" path; the 3 in-range lines -- Fe I 4779.439/4788.757, Fe II
4620.521 -- measured real EWs of 40.4/66.6/54.3 mA) plus slider-driven
`ex_params` changes, all run with zero exceptions. Rendered Stage 4 to
PNG and visually confirmed both the default (`fit_continuum=False`,
widget's own path) and `fit_continuum=True` (called directly, since
the widget itself doesn't expose that flag) panel layouts populate
correctly, including the local-continuum estimate line on the right
panel.

## 20. Commit reference

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
| `0dd202c` | Line identification (§12); RV linelist cross-check and wavelength-shift sign fix (§13) |
| `802cf66` | Add interactive `EWWidget` (§18) |
| `5905744` | Merge `fix_ew_fitting` (new global/local-continuum `measure_ew()`) (§19) |
| `a00c5d5` | Re-integrate `axes=` support into the merged `measure_ew()`/`plot_ew_fit()` (§19) |
