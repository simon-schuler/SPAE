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


def gauss_ew_err(a, sigma, pcov):
    """Propagate a Gaussian fit's parameter covariance to an EW uncertainty.

    EW = 500*a*sqrt(pi/ln2)*(sigma*2.355) = K*a*sigma (see gauss_ew), so this
    is the standard first-order error propagation for a product of two
    correlated fit parameters: var(EW) = K^2 * (sigma^2*var(a) +
    a^2*var(sigma) + 2*a*sigma*cov(a,sigma)). a is gauss_model's amplitude
    (index 0) and sigma its stddev (index 2), matching gfit_direct's
    parameter order.
    """
    if a == 0 or sigma == 0 or pcov is None:
        return 0.
    K = 500. * np.sqrt(np.pi/np.log(2)) * 2.355
    var_a, var_sigma, cov_a_sigma = pcov[0, 0], pcov[2, 2], pcov[0, 2]
    var_ew = K**2 * (sigma**2*var_a + a**2*var_sigma + 2*a*sigma*cov_a_sigma)
    return np.sqrt(abs(var_ew))
