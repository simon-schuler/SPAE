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
        identical to not having a reference atlas at all.

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
        #a real reference-atlas exclusion should only ever SHRINK how much
        #wing survives -- if it shrinks it below min_points, that's a sign
        #the reference cross-check isn't well-conditioned here (e.g. poor
        #atlas coverage, a bad resolving-power estimate), not that this
        #line's continuum can't be estimated at all -- fall back to the
        #median/MAD-only clip rather than failing outright
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
