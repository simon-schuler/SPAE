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
it can't inherit a bias from a fit that's already wrong): compare each
point's local peak (max flux in a narrow window, `local_window`
Angstroms) against a wider window's peak (`wide_window` Angstroms). In
a normal region, isolated lines still leave points nearby where the
local peak nearly reaches the wider peak. In a densely blended stretch
(or, just as well, inside one broad/strong line -- the same criterion
naturally covers the original broad-line stress case too), the local
peak stays well below the wider peak over an EXTENDED stretch, and that
shortfall is what triggers extra stiffness there. A real large-scale
slope (GRACES) doesn't trigger this: `wide_window` is much narrower
than the slope's own scale, so the local and wide peaks stay close
together even while both decline together across the order. This is
computed ONCE up front from the raw data (two cheap 1-D filter passes),
not re-derived every reweighting iteration, so it adds negligible
runtime.

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
self-reinforcing-bias failure this whole module was built to avoid."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import maximum_filter1d, uniform_filter1d


def fit_als_continuum(wave, flux, err, lam=2e4, p=0.01, n_iter=15, adaptive=True,
                       stiffen_factor=15.0, local_window=3.0, wide_window=25.0,
                       severity_threshold=0.08, low_reject_sigma=2.5):
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
        since it's no longer solely responsible for that.
    p : weight floor for points far below the current fit (as a fraction
        of full trust), 0 < p < 0.5.
    n_iter : number of reweighting iterations.
    adaptive : if True (default), compute a per-point stiffness
        multiplier from local-vs-wide peak shortfall (see module
        docstring) and use `lam * multiplier` instead of a flat `lam`.
    stiffen_factor : maximum multiple of `lam` applied where local peaks
        fall furthest short of the wider window's peak.
    local_window, wide_window : Angstrom widths of the two peak-finding
        windows (see module docstring for what each is for).
    severity_threshold : fractional local-vs-wide peak shortfall below
        which a point is treated as normal (no stiffening) -- guards
        against stiffening every ordinary line-to-line gap, not just
        genuinely extended blended/broad stretches.
    low_reject_sigma : noise-width (in units of `err`) over which trust
        decays for points below the current fit -- points within
        roughly this many sigma of the fit are trusted close to fully
        (ordinary scatter, not absorption); points further below decay
        smoothly toward the rejection floor `p`. There is no equivalent
        cap above the fit -- see module docstring for why adding one
        made things worse, not better.

    Returns
    -------
    pred : the fitted continuum, evaluated at every point in wave
    pred_var : formal variance of the final iteration's residuals,
        broadcast to every point (a single order-wide number, not a
        per-point predictive variance; nothing downstream currently
        consumes more than that)
    """
    L = len(flux)
    d2 = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(L - 2, L))
    base_weight = 1.0 / err**2

    lam_vec = np.full(L - 2, lam)
    if adaptive and L > 4:
        dx = (wave[-1] - wave[0]) / (L - 1)
        local_pts = max(int(round(local_window / dx)), 1)
        wide_pts = max(int(round(wide_window / dx)), local_pts + 1)
        local_max = maximum_filter1d(flux, size=local_pts, mode='nearest')
        wide_max = maximum_filter1d(flux, size=wide_pts, mode='nearest')
        shortfall = uniform_filter1d(1.0 - local_max / wide_max, size=local_pts, mode='nearest')
        severity = np.clip((shortfall - severity_threshold) / (1.0 - severity_threshold), 0.0, 1.0)
        mult = 1.0 + (stiffen_factor - 1.0) * severity
        lam_vec = lam * mult[1:-1]  # align to d2's L-2 interior rows

    penalty = d2.T @ sparse.diags(lam_vec, 0) @ d2
    w = base_weight.copy()
    pred = flux.copy()
    for _ in range(n_iter):
        W = sparse.diags(w, 0)
        pred = spsolve((W + penalty).tocsc(), w * flux)
        z = (flux - pred) / err
        below_decay = np.exp(-0.5 * (z / low_reject_sigma)**2)
        w = np.where(z >= 0, base_weight * (1.0 - p),
                     base_weight * (p + (1.0 - 2.0 * p) * below_decay))

    resid = flux - pred
    pred_var = np.full_like(wave, np.var(resid))
    return pred, pred_var
