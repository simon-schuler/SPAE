"""
Effective radial-velocity measurement from a handful of strong, well-
identified spectral lines, and its application as a single, physically
correct multiplicative (v/c) wavelength shift applied to every order.

This is an alternative to estimate_shift()'s per-order grid search against
a reference spectrum. Where estimate_shift() needs each order to overlap a
reference spectrum's order (and treats every order's shift as an
independent free parameter, in Angstrom, regardless of wavelength), this
needs only a few of RV_REFERENCE_LINES to fall somewhere in the spectrum's
coverage, and applies one consistent shift, scaled correctly by each
order's own wavelength -- structurally immune to the reference-spectrum
coverage-gap problem, and confirmed against real data (see conversation/
session notes) that a single velocity is in fact what's really going on:
per-order shifts recovered by estimate_shift() on a real cross-instrument
test showed Delta(lambda)/lambda consistent to ~3% across widely separated
orders, exactly as expected for a true velocity shift and not independent
per-order noise.

No reference (solar) spectrum needed at all -- just the target spectrum's
own normalized flux.
"""

import numpy as np

from .line_profile import gauss_model, gfit_simple

C_KMS = 299792.458  # speed of light, km/s (exact by definition)

# name: (rest wavelength [Angstrom, air], suggested fit window [Angstrom])
# Balmer lines get much wider windows -- their wings are intrinsically much
# broader than metal lines like Ca II/Mg b/Na D.
RV_REFERENCE_LINES = {
    'Ca II K':  (3933.66, 1.5),
    'Ca II H':  (3968.47, 1.5),
    'H-delta':  (4101.73, 4.0),
    'H-gamma':  (4340.46, 4.0),
    'H-beta':   (4861.35, 4.0),
    'Mg b2':    (5172.68, 1.5),
    'Mg b1':    (5183.60, 1.5),
    'Na D2':    (5889.95, 1.5),
    'Na D1':    (5895.92, 1.5),
    'H-alpha':  (6562.79, 4.0),
}


def measure_line_velocity(wave, flux, rest_wavelength, window, min_depth=0.02):
    """
    Fit one absorption line's center and return the velocity (km/s)
    implied relative to rest_wavelength, or None if the line isn't usable
    here (out of range, too shallow, or the fit didn't converge to
    something reasonable).

    Follows the same convention used elsewhere in this package
    (Spectrum_Data.measure_ew()): the continuum-normalized profile is
    depth-inverted (absorption dip -> positive peak at ~0 baseline) before
    fitting gauss_model via curve_fit, since gauss_model is written for a
    positive peak.
    """
    if rest_wavelength - window < wave.min() or rest_wavelength + window > wave.max():
        return None  # not covered here, with margin for the fit window

    mask = (wave >= rest_wavelength - window) & (wave <= rest_wavelength + window)
    if mask.sum() < 5:
        return None

    x = wave[mask]
    y = 1.0 - flux[mask]  # invert: absorption dip -> positive peak, continuum -> ~0

    if y.max() < min_depth:
        return None  # line not really there (too shallow to trust)

    bf, err, p0 = gfit_simple(x, y, rest_wavelength, window / 4.0, 0.0)
    if bf[0] == 0 and bf[1] == 0:
        return None  # gfit_simple's failure sentinel

    fitted_center = bf[1]
    if abs(fitted_center - rest_wavelength) > window / 2.0:
        return None  # fit wandered too far to trust as this line

    return C_KMS * (fitted_center - rest_wavelength) / rest_wavelength


def measure_effective_rv(spectrum, lines=None, min_depth=0.02, sigma_clip=3.0):
    """
    Measure one effective RV (km/s) from whichever of `lines` fall within
    `spectrum`'s wavelength coverage. Uses spectrum.normalized_flux, so
    normalize()/normalize_all() must be run first.

    Parameters
    ----------
    spectrum : Spectrum_Data
    lines : {name: (rest_wavelength, window)}, optional -- defaults to
        RV_REFERENCE_LINES
    min_depth : float -- minimum line depth (in normalized flux) to trust
    sigma_clip : float -- reject individual line velocities more than this
        many standard deviations from the median (guards against a
        misidentified or blended line), only applied when >2 lines matched

    Returns
    -------
    rv : float or None (km/s) -- None if no usable line was found
    rv_err : float (km/s) -- standard error on the mean of the lines used
    used : list of (name, rest_wavelength, order, velocity) -- diagnostics
    """
    if lines is None:
        lines = RV_REFERENCE_LINES

    velocities = []
    used = []
    for name, (rest_wavelength, window) in lines.items():
        for order in range(len(spectrum.wavelength)):
            wave = spectrum.wavelength[order]
            if wave.min() <= rest_wavelength <= wave.max():
                v = measure_line_velocity(wave, spectrum.normalized_flux[order],
                                           rest_wavelength, window, min_depth)
                if v is not None:
                    velocities.append(v)
                    used.append((name, rest_wavelength, order, v))
                break  # only need the first covering order for this line

    if not velocities:
        return None, None, []

    velocities = np.array(velocities)
    if len(velocities) > 2:
        med = np.median(velocities)
        std = velocities.std()
        if std > 0:
            good = np.abs(velocities - med) < sigma_clip * std
            if good.sum() > 0:
                velocities = velocities[good]
                used = [u for u, g in zip(used, good) if g]

    rv = float(np.mean(velocities))
    rv_err = float(velocities.std() / np.sqrt(len(velocities))) if len(velocities) > 1 else 0.0
    return rv, rv_err, used
