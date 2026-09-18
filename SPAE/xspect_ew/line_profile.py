"""Line-window finding and Gaussian profile fitting for EW measurement."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def get_line_window(line, wave, flux, left_bound, right_bound,
                     line_input, window_size=1.5):
    boundaries = [0, 0]
    #if no line is specified, auto fine best line guess
    if line_input == 0.0:
        #find line tip
        left_look = np.where((wave <= line)&(wave >= line - 0.1))
        right_look = np.where((wave >= line)&(wave <= line + 0.1))
        #find min
        mins = [flux[left_look].min(),flux[right_look].min()]
        best_line_guess = wave[np.where(flux == np.min(mins))][0]
    else:
        best_line_guess = line_input

    #get_window around line
    window = np.where((wave >= best_line_guess-window_size/2.0)&(wave <= best_line_guess+window_size/2.0))

    #calc derivative
    dy = np.gradient(flux[window])
    dy_std = dy.std()

    #if no left or right bound given set using std
    auto_bound_l = False
    auto_bound_r = False
    if left_bound == 0:
        dy_l = dy_std/2.0
        auto_bound_l = True
    if right_bound == 0:
        dy_r = dy_std/2.0
        auto_bound_r = True

    #if no line boundaries specified auto find boundaries
    if auto_bound_l:
        left_look = np.where(wave[window] <= best_line_guess - 0.05)
        dy1_left = np.where((dy[left_look] < dy_l)&(dy[left_look] > (-1)*dy_l))
        if len(wave[window][left_look][dy1_left]) ==0:
            print('line ',line,' very close to edge or dy selection value too small')
            print('will attempt to remeasure, if not possible, add line to exclude lines list in .measure_all_ew() function')
            plt.clf()
            plt.plot(wave[window],flux[window])
            plt.plot([line,line],[0.95,1.0], 'k')
            plt.annotate(str(line), xy=[line,1.01])
            plt.plot([best_line_guess,best_line_guess],[0.95,1.0], 'k--')
            plt.annotate(str(best_line_guess), xy=[best_line_guess,1.01])
            plt.show()
            return 1,1,0,0

        else:
            boundaries[0] = wave[window][left_look][dy1_left][-1]
    else:
        boundaries[0] = left_bound
    if auto_bound_r:
        right_look = np.where(wave[window] >= best_line_guess + 0.05)
        dy1_right = np.where((dy[right_look] < dy_r)&(dy[right_look] > (-1)*dy_r))
        if len(wave[window][right_look][dy1_right]) ==0:
            print('line ',line,' very close to edge or dy selection value too small')
            print('will attempt to remeasure, if not possible, add line to exclude lines list in .measure_all_ew() function')
            plt.clf()
            plt.plot(wave[window],flux[window])
            plt.plot([line,line],[0.95,1.0], 'k')
            plt.annotate(str(line), xy=[line,1.01])
            plt.plot([best_line_guess,best_line_guess],[0.95,1.0], 'k--')
            plt.annotate(str(best_line_guess), xy=[best_line_guess,1.01])
            plt.show()
            return 0,0,1,1

        else:
            boundaries[1] = wave[window][right_look][dy1_right][0]
    else:
        boundaries[1] = right_bound

    return window,best_line_guess, boundaries,dy


def gauss_model(x,A,mu,sigma, baseline):
    return A*np.exp(-(x-mu)**2/2/sigma**2) + baseline


def gfit_simple(x_array, y_array, mu, sigma, baseline):
    """Single unweighted Gaussian fit. Used by radial_velocity.py for line-
    center/RV fitting; measure_ew()'s EW fitting uses gfit_direct() below,
    which adds per-pixel error weighting."""
    A = y_array.max()
    p0 = [A, mu, sigma, baseline]
    try:
        bf, cov = curve_fit(gauss_model, x_array, y_array, p0)
        return bf, np.sqrt(np.diag(cov)), p0
    except:
        bf, cov = [0,0,0,0],None
        return bf, cov, p0


def gfit_direct(x_array, y_array, y_err, mu, sigma, baseline):
    """Single weighted least-squares Gaussian fit to real (unsmoothed) data.

    Replaces the old GP-smoothing + 500x Monte-Carlo-resampled curve_fit
    scheme: one scipy.optimize.curve_fit call, weighted by the actual
    per-pixel measurement errors (absolute_sigma=True so the returned
    covariance is on the same absolute scale as y_err, not just relative).
    Parameter uncertainties -- and, via gauss_ew_err(), the EW uncertainty
    -- come directly from that covariance matrix instead of from the
    scatter of repeated fits.

    Returns
    -------
    bf : ndarray [A, mu, sigma, baseline], or None if the fit failed
    pcov : full 4x4 parameter covariance matrix, or None if the fit failed
    p0 : initial guess used
    """
    A = y_array.max()
    p0 = [A, mu, sigma, baseline]
    try:
        bf, pcov = curve_fit(gauss_model, x_array, y_array, p0=p0,
                              sigma=y_err, absolute_sigma=True)
        return bf, pcov, p0
    except (RuntimeError, ValueError):
        return None, None, p0


def estimate_local_continuum(x, y, y_err, min_points=5, clip_sigma=3.0, ref_keep=None):
    """Robust local (flat) continuum level from a set of presumed-
    continuum points -- e.g. a line's wing, outside its own detected
    boundary -- instead of assuming the global normalization already put
    this window's continuum at exactly norm.

    Deliberately a separate ESTIMATION step, not a parameter jointly
    fit alongside the line: letting continuum and line amplitude trade
    off against each other in one fit, seeded from only the line's own
    handful of core points, was confirmed in practice to overfit badly on
    weaker lines (many >50% EW swings, some the wrong direction, on a
    real test spectrum) -- a local continuum should be set by the many
    nearby continuum points, not by the line's own few points fighting an
    optimizer for it.

    Two-pass: first a straight median/MAD pass rejects points that are
    themselves part of a different, deeper absorption feature (not just
    noise -- a plain sigma-clip breaks down if that feature occupies a
    large fraction of the wing, so this is a coarse defense, not a
    substitute for the caller excluding an obviously separate line's own
    core); second, a photon-noise-weighted mean (a pure vertical bias, no
    slope) of the surviving points. A fitted slope was tried and dropped:
    it gave a one-sided contamination (a neighboring line's wing entering
    only one side of the window, not caught by the clip above -- e.g. Fe I
    5587.574 in the bundled sunr.fits sample) direct leverage to tilt the
    whole local continuum, visibly biasing it low on the contaminated
    side. A flat bias can't be tilted that way; it can still be pulled
    off-level if contamination survives the clip, but not systematically
    worse on one side of the window than the other.

    ref_keep : optional boolean mask into x/y, e.g. from
        reference_atlas.reference_continuum_mask() -- an independent,
        externally-sourced exclusion (a real feature confirmed against a
        high-S/N reference spectrum, too shallow for THIS spectrum's own
        noise to catch via the median/MAD clip above) ANDed into the clip
        below rather than replacing it. None (default) leaves behavior
        identical to not having a reference atlas at all. If the atlas
        flags nearly everything in the candidate set as non-continuum
        (fewer than min_points survive ref_keep alone), that's respected
        as a real, atlas-driven contamination verdict -- see the comment
        at its use below for why this can (correctly) return a hard
        failure (0., inf) rather than quietly falling back to the atlas-
        unaware clip.

    Returns
    -------
    c0 : the fitted local continuum level (flat, i.e. no slope)
    c0_err : standard error on c0 -- large when few/noisy wing points
        actually constrain it, so a correction from a poorly-sampled wing
        doesn't get treated as confidently as one from a clean,
        well-sampled one (see measure_ew()'s use of this in its EW error
        budget)
    keep : boolean mask into x/y of points actually used (after clipping)
    """
    keep = np.zeros(len(x), dtype=bool)
    if len(x) < min_points:
        return 0., np.inf, keep

    med = np.median(y)
    mad = np.median(np.abs(y-med)) * 1.4826
    clip = max(mad, np.median(y_err))
    keep = y > (med - clip_sigma*clip)
    if ref_keep is not None:
        if ref_keep.sum() < min_points:
            #the atlas itself found fewer than min_points points that look
            #like continuum ANYWHERE in this candidate set -- a real,
            #atlas-driven contamination verdict (confirmed on Fe I
            #6220.776's blue wing and Fe II 5234.625's WHOLE window: both
            #matched a window riddled with real, broad absorption -- the
            #atlas's own flux tracks our data's shape closely there, it's
            #not an atlas artifact), not noise from a marginal disagreement
            #with the median/MAD clip below -- respect it rather than
            #silently reverting to the atlas-unaware clip, which would
            #override a unanimous, correct contamination verdict with
            #contaminated data (confirmed: this was silently producing a
            #biased-low continuum on 5234.625 before this fix). A caller
            #with no usable continuum left here isn't being deprived of a
            #fallback that would have helped -- it's being told the truth:
            #this window doesn't have one. measure_ew()'s auto_widen retry
            #already exists to try a wider window when that happens, and
            #check_for_flags() already excludes a line whose error comes
            #back too large as a result -- an honest large error beats a
            #confident wrong number.
            keep = keep & ref_keep
        else:
            #the atlas DID find enough usable continuum somewhere in this
            #candidate set -- a low overlap specifically with THIS clip's
            #own surviving points is more likely a marginal/edge
            #disagreement (e.g. poor atlas coverage, a bad resolving-power
            #estimate) than "no continuum exists here", so fall back to
            #the median/MAD-only clip rather than failing outright
            keep_with_ref = keep & ref_keep
            if keep_with_ref.sum() >= min_points:
                keep = keep_with_ref
    if keep.sum() < min_points:
        return 0., np.inf, keep

    yc, ec = y[keep], np.clip(y_err[keep], 1e-6, None)
    w = 1./ec**2
    c0 = np.sum(w*yc) / np.sum(w)
    c0_err = 1./np.sqrt(np.sum(w))
    return c0, c0_err, keep


def _weighted_mean_err(y, y_err):
    ec = np.clip(y_err, 1e-6, None)
    w = 1./ec**2
    return np.sum(w*y)/np.sum(w), 1./np.sqrt(np.sum(w))


def _internal_trend_significant(xs, ys, errs, x0, thresh):
    """Near-half vs far-half (relative to x0) comparison within ONE side's
    own already-clipped points. A real, benign gradient's own side should
    look flat internally; a still-recovering contamination tail (e.g. Fe I
    5587.574's truncated neighbor) forms its own smooth trend even among
    points that individually survived that side's median/MAD + reference-
    atlas clip -- this is what catches it. Returns (trend_found, sig), with
    sig=None (trend_found=False) when there aren't enough points on this
    side to test either way -- too little data to detect a trend isn't
    evidence there isn't one, but it also shouldn't count against a side
    already past the caller's own min_points floor.
    """
    if len(xs) < 4:
        return False, None
    order = np.argsort(np.abs(xs - x0))
    half = len(xs)//2
    near_idx, far_idx = order[:half], order[half:]
    if len(near_idx) < 2 or len(far_idx) < 2:
        return False, None
    mean_near, err_near = _weighted_mean_err(ys[near_idx], errs[near_idx])
    mean_far, err_far = _weighted_mean_err(ys[far_idx], errs[far_idx])
    combined = np.sqrt(err_near**2 + err_far**2)
    sig = abs(mean_far - mean_near)/combined if combined > 0 else 0.
    return sig > thresh, sig


def estimate_local_continuum_sloped(x, y, y_err, x0, min_points=5, clip_sigma=3.0,
                                     ref_keep=None, sig_thresh=3.0, internal_sig_thresh=3.0):
    """DIAGNOSTIC / under evaluation -- not yet used for the reported EW.

    Conditionally allow a linear (sloped) local continuum instead of
    estimate_local_continuum()'s flat bias, WITHOUT reopening the original
    failure mode that got the slope term removed in the first place (see
    estimate_local_continuum()'s docstring): an unconditional linear fit
    over the raw wing points lets one-sided contamination masquerade as a
    trend, since a monotonic decline on just one side is indistinguishable
    from a real slope to an ordinary least-squares fit.

    Instead of fitting the raw points directly, this runs
    estimate_local_continuum() -- the SAME robust median/MAD + optional
    reference-atlas clip, unchanged -- independently on the blue (x < x0)
    and red (x > x0) halves of the wing. A slope is only used if ALL of
    these hold:
      1. each side independently retains at least min_points after its
         own clip (not the combined total) -- a side starved down to a
         handful of points by contamination can't anchor a slope endpoint
      2. NEITHER side shows a significant internal near-vs-far trend of
         its own (see _internal_trend_significant()) -- condition 1 alone
         was confirmed INSUFFICIENT on real data: Fe I 5587.574's
         contaminated side kept 21 points that individually survived the
         per-side clip (they're mutually consistent with EACH OTHER, just
         all part of the same still-recovering contamination tail), so it
         passed condition 1 and even a huge disagreement significance
         (30.5 sigma) in condition 3 below, producing a visibly wrong,
         upward-tilting continuum that overshot the real red-wing data.
         This condition catches exactly that: a benign gradient's own
         side should look flat internally; a recovering tail won't.
      3. the two sides' independently-fitted levels disagree by more than
         sig_thresh times their combined uncertainty -- so the slope only
         engages when the data actually demands a trend, not whenever
         ordinary noise happens to differ between the two sides
    Fe I 5587.574 (truncated neighbor) and Fe I 5579.335 (6-point-starved
    wing) both fail conditions 1 or 2 and fall back to the flat estimate;
    Fe I 5522.447 (a genuine ~1-2% real level difference between wings,
    both independently well-populated, uncontaminated, and internally
    flat) is the motivating case where all three pass.

    sig_thresh, internal_sig_thresh, and min_points are exposed because
    there isn't yet a principled a priori choice for any of them -- this
    function exists to let all three be tuned empirically against real
    lines (see measure_ew()'s parallel flat-vs-sloped diagnostic output)
    before anything here becomes a default.

    Returns
    -------
    c0, c1 : continuum level (at x0) and slope -- c1 is exactly 0. when
        the flat fallback was used
    c0_err : uncertainty on c0 -- from estimate_local_continuum() directly
        for the flat fallback, or propagated as a linear interpolation
        between the two independently-fitted side levels otherwise
    keep : boolean mask into x/y of points used across BOTH sides
    used_slope : bool, whether the slope was actually engaged
    diagnostics : dict with n_blue, n_red, c0_blue, c0_blue_err, c0_red,
        c0_red_err, significance (the disagreement/combined-error ratio
        actually compared against sig_thresh), and internal_trend_blue/
        internal_trend_red/internal_sig_blue/internal_sig_red (condition
        2's per-side result) -- for inspecting near misses, not just
        yes/no, while tuning thresholds
    """
    c0_flat, c0_err_flat, keep_flat = estimate_local_continuum(
        x, y, y_err, min_points=min_points, clip_sigma=clip_sigma, ref_keep=ref_keep)
    diagnostics = {'n_blue': 0, 'n_red': 0, 'c0_blue': None, 'c0_blue_err': None,
                   'c0_red': None, 'c0_red_err': None, 'significance': None,
                   'internal_trend_blue': None, 'internal_sig_blue': None,
                   'internal_trend_red': None, 'internal_sig_red': None}

    blue = x < x0
    red = ~blue
    if blue.sum() < min_points or red.sum() < min_points:
        return c0_flat, 0., c0_err_flat, keep_flat, False, diagnostics

    ref_keep_blue = ref_keep[blue] if ref_keep is not None else None
    ref_keep_red = ref_keep[red] if ref_keep is not None else None
    c0_b, err_b, keep_b = estimate_local_continuum(
        x[blue], y[blue], y_err[blue], min_points=min_points, clip_sigma=clip_sigma,
        ref_keep=ref_keep_blue)
    c0_r, err_r, keep_r = estimate_local_continuum(
        x[red], y[red], y_err[red], min_points=min_points, clip_sigma=clip_sigma,
        ref_keep=ref_keep_red)
    diagnostics.update(n_blue=int(keep_b.sum()), n_red=int(keep_r.sum()),
                        c0_blue=c0_b, c0_blue_err=err_b, c0_red=c0_r, c0_red_err=err_r)

    if keep_b.sum() < min_points or keep_r.sum() < min_points:
        return c0_flat, 0., c0_err_flat, keep_flat, False, diagnostics

    trend_b, sig_b = _internal_trend_significant(
        x[blue][keep_b], y[blue][keep_b], y_err[blue][keep_b], x0, internal_sig_thresh)
    trend_r, sig_r = _internal_trend_significant(
        x[red][keep_r], y[red][keep_r], y_err[red][keep_r], x0, internal_sig_thresh)
    diagnostics.update(internal_trend_blue=trend_b, internal_sig_blue=sig_b,
                        internal_trend_red=trend_r, internal_sig_red=sig_r)
    if trend_b or trend_r:
        return c0_flat, 0., c0_err_flat, keep_flat, False, diagnostics

    combined_err = np.sqrt(err_b**2 + err_r**2)
    significance = abs(c0_r - c0_b)/combined_err if combined_err > 0 else 0.
    diagnostics['significance'] = significance
    if significance <= sig_thresh:
        return c0_flat, 0., c0_err_flat, keep_flat, False, diagnostics

    x_b_mean = np.mean(x[blue][keep_b])
    x_r_mean = np.mean(x[red][keep_r])
    c1 = (c0_r - c0_b) / (x_r_mean - x_b_mean)
    c0 = c0_b + c1*(x0 - x_b_mean)

    #linear-interpolation error propagation between the two independently
    #-fitted side levels, evaluated at x0
    w_b = (x_r_mean - x0) / (x_r_mean - x_b_mean)
    w_r = (x0 - x_b_mean) / (x_r_mean - x_b_mean)
    c0_err = np.sqrt((w_b*err_b)**2 + (w_r*err_r)**2)

    keep = np.zeros(len(x), dtype=bool)
    keep[blue] = keep_b
    keep[red] = keep_r
    return c0, c1, c0_err, keep, True, diagnostics


def gauss_ew(a, fwhm):
    if a == 0 or fwhm == 0:
        return 0
    else:
        return 500.*a*np.sqrt(np.pi/np.log(2))*fwhm #From Adamow pyMOOG ew measure


def gauss_model_err(x, bf, pcov):
    """1-sigma uncertainty envelope of gauss_model(x, *bf) for plotting a
    shaded fit-quality band, propagated from the fit's parameter covariance
    matrix via the model's Jacobian: var(x) = J(x) . pcov . J(x)^T, the
    standard linear error propagation used for a fitted curve's confidence
    band. Returns an array the same shape as x (all zeros if pcov is None
    or the fit is degenerate).
    """
    A, mu, sigma_, baseline = bf
    if pcov is None or sigma_ == 0:
        return np.zeros_like(x, dtype=float)
    g = np.exp(-(x-mu)**2/(2*sigma_**2))
    dA = g
    dmu = A*g*(x-mu)/sigma_**2
    dsigma = A*g*(x-mu)**2/sigma_**3
    dbase = np.ones_like(x, dtype=float)
    J = np.stack([dA, dmu, dsigma, dbase], axis=-1)  # (N, 4)
    var = np.einsum('ni,ij,nj->n', J, pcov, J)
    return np.sqrt(np.abs(var))


#EW = 500*A*sqrt(pi/ln2)*FWHM = 500*A*sqrt(pi/ln2)*(sigma*2.355) = EW_K*A*sigma
#(see gauss_ew) -- exposed so callers propagating an extra, independently-
#estimated source of amplitude uncertainty (e.g. measure_ew()'s local
#continuum uncertainty) into an EW uncertainty can reuse the same constant.
EW_K = 500. * np.sqrt(np.pi/np.log(2)) * 2.355


def gauss_ew_err(a, sigma, pcov):
    """Propagate a Gaussian fit's parameter covariance to an EW uncertainty.

    EW = EW_K*a*sigma, so this is the standard first-order error
    propagation for a product of two correlated fit parameters:
    var(EW) = EW_K^2 * (sigma^2*var(a) + a^2*var(sigma) +
    2*a*sigma*cov(a,sigma)). a is gauss_model's amplitude (index 0) and
    sigma its stddev (index 2), matching gfit_direct's parameter order.
    """
    if a == 0 or sigma == 0 or pcov is None:
        return 0.
    var_a, var_sigma, cov_a_sigma = pcov[0, 0], pcov[2, 2], pcov[0, 2]
    var_ew = EW_K**2 * (sigma**2*var_a + a**2*var_sigma + 2*a*sigma*cov_a_sigma)
    return np.sqrt(abs(var_ew))
