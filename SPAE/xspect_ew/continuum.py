"""Continuum fitting for Spectrum_Data.normalize().

Every approach tried before this one shared the same structural
weakness: a HARD point-selection mask deciding "this point is definitely
continuum" / "this point is definitely a line," computed from some
threshold that isn't locally aware. Whenever that mask left a stretch of
an order with too few (or zero) trusted points, whatever was fit to it
next either couldn't be constrained there or got it badly wrong:

- Continuum_scan's local-window percentile selection: a ~1.5 A window
  sitting mostly or entirely inside a broad/strong line has no true-
  continuum points to select at all -- measured 20-70% continuum
  underestimate near a synthetic 3-A-wide line.
- A globally-refit Gaussian Process, iteratively sigma-clipped against
  that same kind of hard selection: did not fix the broad-line bug at
  all (still -26%/-19% error) -- a GP's local flexibility bends into a
  contaminated region just as easily as the local-window selection did,
  and low-side sigma-clip rejection can't self-correct once the first
  fit is already biased low across a whole contiguous stretch.
- A single low-order global polynomial: fixed the broad-line bug (every
  point influences the whole curve, so it can't be dragged into a local
  gap) but was too rigid to track real order-wide continuum curvature.
- A piecewise cubic spline with sparse knots: tracked curvature well in
  isolated synthetic tests, but reintroduced the same hard-mask fragility
  in three different real-data forms as the initial-guess/edge-handling
  around it was patched again and again: (1) a milder version of the
  original broad-line bias with a bad initial guess, (2) catastrophic
  divergence (-1.2 MILLION counts) right at real order edges where a
  knot's span had too few good points, (3) a real GRACES order's own
  large-scale continuum decline getting excluded from the selection
  mask entirely (same mechanism that correctly excludes absorption
  lines, wrongly applied to a real trend), and (4) an outright negative
  continuum in the INTERIOR of a real Keck order, where a 20-A stretch
  of ordinary-looking continuum left only 12 of 733 points passing a
  global percentile threshold, starving that spline segment completely.
  Patching each of these individually kept surfacing a new failure on
  the next real spectrum -- evidence the underlying strategy (hard
  selection mask + whatever-fits-it) was the problem, not any one
  threshold or margin choice.

fit_als_continuum() replaces the hard mask with Asymmetric Least Squares
(AsLS) smoothing (Eilers & Boelens 2005; standard in Raman/chromatography
baseline correction, applied here to fit continuum as an upper envelope
rather than a baseline under emission peaks). Every point gets a SOFT
weight, iteratively re-estimated: points above the current fit (likely
real continuum) are upweighted, points below (likely absorption) are
downweighted, but never fully excluded. A smoothness penalty (`lam`),
not knot placement, controls rigidity -- continuously, not via discrete
segments -- so there is no way for a stretch of an order to end up
completely unconstrained the way a starved spline segment could: the
fit there is still tied to its neighbors through the same smoothness
penalty that applies everywhere else. This isn't one more patch to the
old approach; it's a different algorithm chosen specifically because it
can't reproduce this whole family of failures by construction.

A single GLOBAL `lam` still has a real limit, though: a real Keck order
(6198-6309 A) has a ~30 A stretch (6272-6300 A) so densely packed with
blended Fe I lines that local peaks never fully recover to true
continuum -- a `lam` loose enough to track a real large-scale continuum
slope (needed for a GRACES order with a genuine ~5x decline toward its
edges) isn't stiff enough to resist being dragged down through that
stretch, and a `lam` stiff enough to resist it starts overshooting past
GRACES's real edges. No single global value serves both well.

fit_als_continuum() therefore stiffens the penalty ADAPTIVELY, per
point, based on a measurement that's independent of the fit itself (so
it can't inherit a bias from a fit that's already wrong). Two DIFFERENT
things can make a stretch untrustworthy for the base (lam, p) to handle
well, and one metric can't catch both without breaking on a real slope
-- so there are two signals, combined via max():

1. A genuinely BROAD/deep trough (one strong or wide line, a saturated
telluric band): compare each point's local peak (max flux in
`local_window` Angstroms) against a wider window's peak (`wide_window`
Angstroms). In a normal region, isolated lines still leave points
nearby where the local peak nearly reaches the wider peak; across a
broad trough, the local peak stays well below the wider peak over an
EXTENDED stretch. A real large-scale slope (GRACES) doesn't trigger
this: `wide_window` is much narrower than the slope's own scale, so
local and wide PEAKS stay close together even while both decline
together across the order (both are governed by the same nearby-noise
maximum, largely independent of the slope's own rate).

2. WIDESPREAD but individually modest weak/blended lines -- confirmed a
real, common case this session, and invisible to signal 1: a single
bright pixel anywhere in `local_window` satisfies "local peak reaches
wide peak" even when 95-99% of the surrounding points sit measurably
below it. Detected instead by comparing the local MEDIAN to the local
MAX (same narrow window, not local-vs-wide) of flux first divided by a
wide-window MEDIAN trend -- the division removes a real large-scale
slope's own contribution (which would otherwise make local median vs.
local max differ across just `local_window` purely from the slope's own
rate, a second, smaller version of the same false-positive signal 1
was designed to avoid), while a genuine weak-line forest's dips survive
the division untouched, since they're far narrower than `wide_window`.
Confirmed on real MAROON-X data: two different orders' "clean" regions
(one flagged by the user, one this module's own earlier verification
had called fine) both showed the ORIGINAL (signal-1-only) severity
staying at exactly 0 throughout that stretch, silently leaving it
uncorrected.

Both shortfalls are computed ONCE up front from the raw data (a few
cheap 1-D filter passes), not re-derived every reweighting iteration,
so this adds negligible runtime. Signal 2 uses its own, much smaller
threshold/scale (`severity_threshold`/`severity_scale`) than signal 1's
fixed 0.08/0.92 -- a widespread-but-moderate forest saturates nowhere
near the 0.3-1.0 shortfall a genuinely broad trough does.

A second, more subtle bias remained even with adaptive stiffening: the
weight update used a HARD step at the current fit -- any point below it
got flat weight p, no matter how close to the fit it actually was. Real
Poisson noise scatters points both above AND below the true continuum
by chance, so roughly half of genuinely-unabsorbed points sit slightly
below the current estimate at any moment -- treating all of them as
"probably absorption" biases the fit toward the upper envelope of the
noise rather than its mean. Confirmed on a real Keck order (5698-5800
A): the fit tracked at or even slightly ABOVE the local peak everywhere
checked, when it should sit at the mean continuum level. The below-fit
side of the weight is now a smooth, noise-scaled decay (using `err`,
the same per-point noise already used for the base inverse-variance
weight) instead of a hard step: points within roughly `low_reject_sigma`
worth of noise below the fit are trusted close to fully (consistent
with ordinary scatter, not absorption), so the fit is pulled toward
their AVERAGE rather than only ever being told "no" by anything sitting
below it; points further below decay smoothly toward the rejection
floor `p`, same as real absorption always needed.

The above-fit side deliberately keeps its ORIGINAL flat (1-p) trust,
not a symmetric sigma-based decay -- tried that first and it was worse:
whenever the fit is transiently biased low (which it always is, early
in the iteration, near any broad/strong line, before the fit has
recovered), genuinely-unabsorbed points nearby look like implausible
upward outliers relative to that still-wrong fit and get wrongly
downweighted too, which only reinforces the bias instead of letting the
fit recover from it (confirmed directly: a real continuum point at
49,667 counts, next to a synthetic broad line, got assigned z=20+
against a fit still sitting at ~45,000 from the line's influence, and
was rejected as if it were a cosmic ray). The smoothness penalty
(`lam`/adaptive stiffening) is what keeps the fit from literally
chasing the single highest point on the upper side -- it doesn't need
its own sigma cap, and adding one reintroduces the same kind of
self-reinforcing-bias failure this whole module was built to avoid.

A third bias, smaller per-point but systematic, showed up specifically
in densely-lined orders: the fit sat measurably (2-10%, worse with more
lines) below the true continuum even in stretches with NO line nearby
at all -- confirmed with a controlled synthetic test (known-flat true
continuum, only the NUMBER of scattered lines varied): 0% error with no
lines, growing to -4.4% with very dense lines, and NOT a convergence
issue (identical at 15 vs 120 iterations -- a genuinely different
equilibrium, not an under-run one). Root cause: the below-fit weight
floor `p` (default 0.01) is small per point but never zero, and a
densely-lined order has MANY absorbed points -- their cumulative pull,
summed over the whole order via the smoothness penalty that ties
everything together, measurably drags down the fit even where no single
absorbed point is nearby. Confirmed directly: shrinking `p` toward zero
on the same synthetic test shrinks the bias toward zero too (0.01 to
0.00001 took -5.8% down to -0.08%). But `p` can't just be made small
everywhere: doing that on the SLOPED synthetic test (a real large-scale
decline, no dense blending) made it much WORSE (clean-region error rose
from ~4% to ~15%) -- a small `p` rejects below-fit points so fast that a
genuine downward slope the fit hasn't caught up to yet gets treated as
absorption and locked out, the same self-reinforcing-bias shape as the
above-fit case, just triggered from the other side. The fix reuses the
SAME severity signal already computed for adaptive stiffening (dense
blending and genuine slopes are exactly what it already tells apart):
`p` is shrunk by up to `p_reduction_factor` in high-severity regions,
left at its lenient base value elsewhere.

A fourth bias is really a mismatch between what the fit targets and
what "continuum" means for real, not-perfectly-clean stellar spectra:
the low_reject_sigma-based decay above was tuned to bring the fit down
to the MEAN of the noise (undoing the earlier "tracks above the local
peak" overshoot), but `low_reject_sigma` (2.5) is wide enough that it
barely discounts anything within about 2 sigma either side -- the bulk
of a Gaussian -- so in practice the fit settles close to the raw
CENTROID of nearby points, not a true upper envelope, even though the
weighting is nominally asymmetric (full trust above, decay below).
Confirmed by inspection: real Poisson noise scatters roughly
symmetrically around that centroid, so about half of any clean stretch
visibly pokes up above the fitted line -- consistent with fitting a
mean, not an envelope. For real echelle data this matters even away
from any resolved line: pervasive weak/blended absorption depresses the
observed centroid below the TRUE (line-free) continuum almost
everywhere, so a mean-tracking fit inherits that depression, worst in
densely-lined stretches (e.g. Keck's red order edges) but present even
in nominally "clean" regions. `target_percentile` (default 80) restores
a controlled upward bias, expressed in units of each point's own
photon-noise sigma rather than an arbitrary fraction of flux (which
would be a wildly different amount of correction depending on
brightness) -- applied as a POST-HOC shift of the already-converged
mean-tracking fit, not folded into the iteration itself, so none of the
above convergence/stability behavior changes; only the final reported
level does.

A fifth issue is specific to real order EDGES, not interior structure:
an echelle order's outermost pixels are where detector/blaze-rolloff
artifacts concentrate (a spurious upward spike immediately followed by a
cliff-like drop, etc.) -- confirmed on a real Keck order (HD_10383,
~4977-4982 A): a spike-then-drop right at the boundary dragged the fit up
across the whole nearby region (curving from ~52000 to ~53500 toward the
edge) even though the rest of the order's true local peaks sat flat
around ~49700-50000. `edge_ignore_aa` (default 2.0) excludes the outer
this-many-Angstroms at EACH end from influencing the fit (zero data-
weight there, same mechanism as `p`'s below-fit decay, just permanent and
unconditional for that stretch) while still returning a continuum value
at every point, including the excluded edges -- the smoothness penalty
extrapolates the curve into them from the trusted interior, rather than
leaving them unfit. This is NOT free: a real, physically genuine blaze
decline that's only captured by data in the last couple of Angstroms
(confirmed on a real solar order, sunb.fits order 8: true continuum
keeps falling all the way to the edge) gets extrapolated instead of
tracked, overshooting there by >100%. There's no universally-correct
choice -- 2 A is a compromise default, small enough to rarely eat into
genuine order-edge curvature on typically-wide (tens of A) orders, large
enough to exclude the kind of narrow edge artifact seen on HD_10383.
Set to 0 to disable and fit every point as before."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import maximum_filter1d, median_filter, uniform_filter1d
from scipy.stats import norm


def fit_als_continuum(wave, flux, err, lam=2e3, p=0.01, n_iter=15, adaptive=True,
                       stiffen_factor=15.0, local_window=3.0, wide_window=25.0,
                       severity_threshold=0.02, severity_scale=0.06, low_reject_sigma=2.5,
                       p_reduction_factor=30.0, target_percentile=80.0, edge_ignore_aa=2.0):
    """
    Fit the continuum as the (noise-aware) upper envelope of flux via
    Asymmetric Least Squares (AsLS) smoothing: iteratively solve the
    penalized weighted least-squares problem

        (W + D^T diag(lam_vec) D) z = W * flux

    where D is the discrete second-difference operator (the penalty
    resists curvature, the same role knot_spacing/polynomial degree
    played before, but as a continuous regularizer instead of discrete
    segments) and W holds each point's weight: the inverse-variance
    measurement weight (1/err^2), further scaled each iteration by
    which side of the current fit the point is on and, for below-fit
    points, how consistent the residual is with pure noise (see module
    docstring) -- full trust above the fit, smoothly decaying trust
    below it as the shortfall grows past ordinary noise. A few
    iterations let the fit settle at the noise-consistent continuum
    level, not chase absorption lines OR treat ordinary downward noise
    scatter as absorption -- without ever needing a hard yes/no
    selection mask.

    Parameters
    ----------
    wave, flux, err : arrays for one order (err used as the base
        per-point weight, 1/err**2 -- e.g. sqrt(flux) for Poisson noise
        -- AND as the noise scale the below-fit decay weighting is
        measured in, so it directly reflects this spectrum's own S/N)
    lam : BASE smoothness penalty (used everywhere adaptive stiffening
        isn't triggered). Larger = more rigid; smaller = more flexible.
        Calibrated against the same synthetic ground-truth test used
        throughout this module's history, now WITH adaptive stiffening
        doing the work of resisting broad/blended regions -- lam itself
        can be smaller than it needed to be without adaptive stiffening,
        since it's no longer solely responsible for that. Recalibrated
        this session (2e4 -> 2e3): the OLD value was left over from
        before adaptive stiffening existed and was never revisited once
        the adaptive mechanism took over that job, leaving the base
        penalty needlessly rigid everywhere else. Confirmed on a real
        Keck order (5760-5800 A, a red-order-edge transition with no
        genuinely broad/deep feature, so adaptive stiffening correctly
        never engages there): even the ACHIEVABLE local envelope
        (between real lines) undershot true continuum by 3-4%, because
        `lam=2e4` was too rigid to bend down fast enough as line
        crowding gradually thickened -- confirmed NOT fixable by adding
        more adaptive stiffening (forcing severity=1 there made it
        dramatically WORSE, down to 50-65%, since over-rigidifying an
        already-too-stiff base penalty just resists tracking real local
        structure even harder). Scanning lam over the same synthetic
        ground-truth test used throughout this module's history found
        2e3 recovers that envelope (~1.00) while still improving (not
        regressing) both the flat and sloped synthetic tests' own error
        metrics -- the old value wasn't apparently just cautious, it
        was too rigid full stop, verified independent of the Keck case
        that motivated re-examining it.
    p : BASE weight floor for points far below the current fit (as a
        fraction of full trust), 0 < p < 0.5, used everywhere adaptive
        stiffening isn't triggered -- see module docstring for why this
        can't just be made small everywhere.
    n_iter : number of reweighting iterations.
    adaptive : if True (default), compute a per-point stiffness
        multiplier from local-vs-wide peak shortfall (see module
        docstring) and use `lam * multiplier` instead of a flat `lam`
        (and shrink `p` by up to `p_reduction_factor` in the same
        regions -- see module docstring for why both are needed).
    stiffen_factor : maximum multiple of `lam` applied where either
        severity signal (see module docstring) reaches 1.
    local_window, wide_window : Angstrom widths of the two windows used
        by BOTH severity signals (see module docstring for what each is
        for).
    severity_threshold, severity_scale : threshold and (linear) span,
        in units of the SIGNAL-2 (widespread weak-line forest) shortfall
        (see module docstring), below/over which severity ramps from 0
        to 1. Signal 1 (broad troughs) uses its own fixed 0.08/0.92,
        unchanged from this module's original calibration -- it doesn't
        need retuning since a genuinely broad/deep trough's shortfall
        saturates close to 1 regardless. Signal 2's shortfall saturates
        far short of 1 for even an obviously-real widespread forest
        (0.03-0.08 is typical, not 0.3-1.0), so it needs its own, much
        smaller default scale to respond to it at all.
    low_reject_sigma : noise-width (in units of `err`) over which trust
        decays for points below the current fit -- points within
        roughly this many sigma of the fit are trusted close to fully
        (ordinary scatter, not absorption); points further below decay
        smoothly toward the rejection floor `p`. There is no equivalent
        cap above the fit -- see module docstring for why adding one
        made things worse, not better.
    p_reduction_factor : maximum factor `p` is divided by in the same
        high-severity regions `stiffen_factor` targets -- see module
        docstring for the cumulative-bias bug this fixes.
    target_percentile : where the reported continuum should sit within
        the LOCAL photon-noise distribution, not within flux itself --
        80 means roughly "mean + 0.84 sigma" (norm.ppf(0.80)), not "80%
        of the flux value." Converted once to a sigma multiplier and
        applied as a flat additive shift, in units of the noise scale AT
        the fitted continuum level (err rescaled by sqrt(pred/flux), so
        it reflects true continuum brightness rather than being
        artificially small inside an absorption line where flux itself
        is depressed) -- see module docstring for why the base iteration
        settles near the noise MEAN and needs this correction on top.
        50 reproduces the old (uncorrected) mean-tracking behavior. The
        noise scale itself is further calibrated against the fit's own
        above-fit residuals before use (see the code just before the
        `return`) -- real extracted spectra don't always match `err`'s
        theoretical Poisson scaling (confirmed on MAROON-X), and this
        keeps the offset's actual SIZE correct even when they don't.
    edge_ignore_aa : Angstroms of each order's outer edge (both ends) to
        exclude from influencing the fit -- see module docstring for the
        real-data motivation (a boundary artifact vs. genuine edge
        curvature) and the tradeoff. 0 disables this (every point fits as
        before). A continuum value is still returned at every point,
        including excluded edges -- the smoothness penalty extrapolates
        into them from the trusted interior.

    Returns
    -------
    pred : the fitted continuum, evaluated at every point in wave
    pred_var : formal variance of the final iteration's residuals,
        broadcast to every point (a single order-wide number, not a
        per-point predictive variance; nothing downstream currently
        consumes more than that)
    """
    # median_filter (unlike maximum_filter1d/uniform_filter1d) rejects
    # non-native-byte-order input, which raw Keck/MAKEE FITS data is
    # (big-endian float32) -- normalize dtype once up front rather than
    # requiring every caller to know that.
    wave = np.asarray(wave, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)
    err = np.asarray(err, dtype=np.float64)

    L = len(flux)
    d2 = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(L - 2, L))
    base_weight = 1.0 / err**2

    # edge_ignore_aa: zero the data-weight within this many Angstroms of
    # EACH end -- see module docstring ("A fifth issue...") for why.
    # trusted/edge_ignore_mask are reused below (typical_weight's median,
    # the below-fit reweighting each iteration, and the target_percentile
    # calibration all need to exclude these points the same way).
    edge_ignore_mask = np.zeros(L, dtype=bool)
    if edge_ignore_aa > 0 and L > 2:
        edge_ignore_mask = ((wave - wave[0] < edge_ignore_aa) |
                             (wave[-1] - wave < edge_ignore_aa))
        base_weight = base_weight.copy()
        base_weight[edge_ignore_mask] = 0.0
    trusted = ~edge_ignore_mask

    # lam/p are absolute numbers, calibrated (via this module's synthetic
    # ground-truth test, test_normalize.py: Poisson noise on a ~50,000-
    # count continuum, giving err~224, weight~4e-10... no -- weight =
    # 1/err^2 ~ 1/224^2 ~ 2e-5) against ONE particular data scale. lam
    # itself doesn't know what units flux/err are in, so (W + penalty)'s
    # balance between "trust the data" and "stay smooth" only behaves the
    # way lam/p were tuned for when the data's typical weight stays near
    # that same ~2e-5 scale. Response-corrected flux (MAROON-X) can sit
    # many orders of magnitude away from it (millions of counts, weight
    # ~1e-7 to 1e-10) -- not just a uniformly bigger/smaller problem, but
    # a qualitatively different, badly-conditioned regime where lam is
    # effectively far stiffer, relative to the data, than it was ever
    # tuned to be. Confirmed on a real MAROON-X order: correcting its
    # error to properly reflect response-division noise (see
    # response_correction.py) pushed typical weight far below this scale
    # and made the AsLS reweighting iteration converge to a non-physical
    # fixed point (predicted continuum ~35x the actual flux) instead of
    # the small, stable correction seen on every other order -- even
    # though the per-point error itself was now correct. Rescaling lam by
    # the order's own typical weight relative to the calibration
    # reference restores the SAME relative balance regardless of what
    # absolute units flux/err happen to be in -- a no-op for Keck/GRACES-
    # scale data (where this ratio is already ~1).
    _CALIBRATION_WEIGHT = 1.0 / 224.0**2
    typical_weight = np.median(base_weight[trusted]) if trusted.any() else np.median(base_weight)
    if typical_weight > 0:
        lam = lam * typical_weight / _CALIBRATION_WEIGHT

    lam_vec = np.full(L - 2, lam)
    p_vec = np.full(L, p)
    if adaptive and L > 4:
        dx = (wave[-1] - wave[0]) / (L - 1)
        local_pts = max(int(round(local_window / dx)), 1)
        wide_pts = max(int(round(wide_window / dx)), local_pts + 1)

        # Signal 1: a genuinely BROAD/deep trough (one strong/wide line,
        # a saturated telluric band) -- local peak vs a WIDER window's
        # peak. Unchanged from the original design; deliberately still
        # max-vs-max (see module docstring for why median-vs-max here
        # would falsely fire on an ordinary large-scale slope).
        local_max = maximum_filter1d(flux, size=local_pts, mode='nearest')
        wide_max = maximum_filter1d(flux, size=wide_pts, mode='nearest')
        shortfall_broad = uniform_filter1d(1.0 - local_max / wide_max, size=local_pts, mode='nearest')
        severity_broad = np.clip((shortfall_broad - 0.08) / 0.92, 0.0, 1.0)

        # Signal 2: WIDESPREAD but individually modest weak/blended lines
        # -- real, common, and invisible to signal 1 (see module
        # docstring). Detrend by a wide-window median first (tracks a
        # real large-scale slope without a real line's much narrower
        # dip), then compare the local median to the local max of that
        # detrended ratio -- both computed over the SAME narrow window,
        # so a genuine slope's own variation across just that window
        # can't masquerade as blending the way a raw local-vs-wide
        # comparison would.
        wide_trend = median_filter(flux, size=wide_pts, mode='nearest')
        detrended = flux / wide_trend
        local_max_dt = maximum_filter1d(detrended, size=local_pts, mode='nearest')
        local_med_dt = median_filter(detrended, size=local_pts, mode='nearest')
        shortfall_density = uniform_filter1d(1.0 - local_med_dt / local_max_dt, size=local_pts, mode='nearest')
        severity_density = np.clip((shortfall_density - severity_threshold) / severity_scale, 0.0, 1.0)

        severity = np.maximum(severity_broad, severity_density)
        mult = 1.0 + (stiffen_factor - 1.0) * severity
        lam_vec = lam * mult[1:-1]  # align to d2's L-2 interior rows
        p_vec = p / (1.0 + (p_reduction_factor - 1.0) * severity)

    penalty = d2.T @ sparse.diags(lam_vec, 0) @ d2
    w = base_weight.copy()
    pred = flux.copy()
    for _ in range(n_iter):
        W = sparse.diags(w, 0)
        pred = spsolve((W + penalty).tocsc(), w * flux)
        z = (flux - pred) / err
        below_decay = np.exp(-0.5 * (z / low_reject_sigma)**2)
        w = np.where(z >= 0, base_weight * (1.0 - p_vec),
                     base_weight * (p_vec + (1.0 - 2.0 * p_vec) * below_decay))
        w[edge_ignore_mask] = 0.0  # keep excluded edges at zero trust every iteration

    resid = flux - pred

    if target_percentile != 50.0:
        sigma_offset = norm.ppf(target_percentile / 100.0)
        flux_floor = np.maximum(flux, 1e-3 * np.median(flux[flux > 0]) if np.any(flux > 0) else 1.0)
        continuum_scale = np.clip(pred / flux_floor, 0.0, 100.0)
        err_at_continuum = err * np.sqrt(continuum_scale)

        # `err` is a THEORETICAL noise scale (Poisson-style sqrt
        # scaling); real extracted spectra don't always match it. An
        # "optimal extraction" algorithm (MAROON-X) combines several raw
        # CCD pixels per output point via inverse-variance-weighted PSF
        # fitting, correlating adjacent output points and leaving LESS
        # actual point-to-point scatter than sqrt(counts) alone would
        # predict -- confirmed directly: real MAROON-X data's observed
        # scatter in genuinely flat stretches was only 30-40% of what
        # err implied, while Keck's matched almost exactly. Trusting
        # err's absolute scale blindly overshot MAROON-X's continuum by
        # several percent (a real, user-caught, spectrum-wide bug), even
        # though the SAME formula was correct for Keck. Rather than a
        # per-instrument fudge factor, calibrate against the fit's OWN
        # above-fit residuals (least likely to be real absorption, by
        # construction, so their spread is a genuine empirical noise
        # estimate) each time: compare the median of positive residuals
        # to what a Gaussian with err's sigma would predict for that same
        # statistic (median of a positive half-normal = 0.6745*sigma),
        # and rescale err_at_continuum by that ratio before applying the
        # offset. A no-op when err already matches reality (Keck/GRACES);
        # self-corrects when it doesn't, without needing to know why.
        above = (resid > 0) & trusted  # exclude edge_ignore_aa's edges from this too
        calib = 1.0
        if above.sum() > 10:
            empirical = np.median(resid[above])
            theoretical = np.median(err_at_continuum[above]) * 0.6744897501960817
            if theoretical > 0:
                calib = np.clip(empirical / theoretical, 0.05, 3.0)

        pred = pred + sigma_offset * calib * err_at_continuum
        resid = flux - pred

    pred_var = np.full_like(wave, np.var(resid[trusted]) if trusted.any() else np.var(resid))
    return pred, pred_var
