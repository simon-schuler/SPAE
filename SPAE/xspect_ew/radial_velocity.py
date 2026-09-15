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
from .line_identification import identify_lines_in_spectrum

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


def _measure_rv_over_lines(spectrum, lines, min_depth, sigma_clip):
    """Shared measurement loop for measure_effective_rv() -- see there."""
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


def measure_effective_rv(spectrum, lines=None, min_depth=0.02, sigma_clip=3.0, min_metal_lines=2):
    """
    Measure one effective RV (km/s) from whichever of `lines` fall within
    `spectrum`'s wavelength coverage. Uses spectrum.normalized_flux, so
    normalize()/normalize_all() must be run first.

    When `lines` is left at its default (RV_REFERENCE_LINES), Balmer
    lines (H-alpha/beta/gamma/delta) are deliberately tried only as a
    fallback, not averaged in with the metal lines (Ca II/Mg b/Na D) from
    the start: confirmed on real data (the bundled Keck solar sample) that
    the two families can disagree by ~10 km/s, a real difference in line
    formation between H and metal lines (broader, more pressure/NLTE-
    sensitive Balmer wings centering less precisely), not measurement
    noise that sigma-clipping over only 3-4 mixed-species lines can
    reliably separate out. Metal lines are used alone whenever at least
    `min_metal_lines` of them are available; Balmer lines are added back
    in only when there aren't enough metal lines to trust alone. An
    explicitly-passed `lines` dict is used exactly as given, with no such
    splitting -- this only changes the DEFAULT set's own behavior.

    Parameters
    ----------
    spectrum : Spectrum_Data
    lines : {name: (rest_wavelength, window)}, optional -- defaults to
        RV_REFERENCE_LINES (with the metal-first behavior described
        above; pass an explicit dict to bypass it entirely).
    min_depth : float -- minimum line depth (in normalized flux) to trust
    sigma_clip : float -- reject individual line velocities more than this
        many standard deviations from the median (guards against a
        misidentified or blended line), only applied when >2 lines matched
    min_metal_lines : minimum number of non-Balmer lines required before
        Balmer lines are excluded from the default set entirely (only
        relevant when `lines` is left at its default).

    Returns
    -------
    rv : float or None (km/s) -- None if no usable line was found
    rv_err : float (km/s) -- standard error on the mean of the lines used
    used : list of (name, rest_wavelength, order, velocity) -- diagnostics
    """
    if lines is not None:
        return _measure_rv_over_lines(spectrum, lines, min_depth, sigma_clip)

    metal_lines = {name: v for name, v in RV_REFERENCE_LINES.items() if not name.startswith('H-')}
    rv, rv_err, used = _measure_rv_over_lines(spectrum, metal_lines, min_depth, sigma_clip)
    if len(used) >= min_metal_lines:
        return rv, rv_err, used
    # not enough metal lines available/usable -- fall back to the full
    # set, Balmer lines included, same as always trying everything at once
    return _measure_rv_over_lines(spectrum, RV_REFERENCE_LINES, min_depth, sigma_clip)


def _rv_from_detections(lines, wavelength, flux, err, pred, radius, min_significance):
    """One identify-and-convert-to-velocity pass -- shared by both
    stages of measure_rv_from_linelist() below."""
    results = identify_lines_in_spectrum(lines, wavelength, flux, err, pred,
                                          search_radius=radius,
                                          min_significance=min_significance)
    velocities = []
    for rest_wavelength, r in zip(lines, results):
        if r['detected']:
            velocities.append(C_KMS * (r['center'] - rest_wavelength) / rest_wavelength)
    return np.array(velocities)


def measure_rv_from_linelist(lines, wavelength, flux, err, pred, search_radius=1.0,
                              refine_radius=0.15, min_significance=3.0, sigma_clip=3.0,
                              min_lines=3):
    """
    Measure one effective RV (km/s) from a full science linelist, as a
    generalization of measure_effective_rv() that isn't tied to a small
    fixed set of named lines -- see Spectrum_Data.apply_rv_shift()'s
    docstring for why this matters (RV_REFERENCE_LINES can be absent
    from a given spectrum's coverage entirely, or individually behave
    inconsistently: confirmed on real data that Balmer lines specifically
    can disagree with metal lines by ~10 km/s, a real astrophysical
    difference in line formation, not just measurement noise, which a
    small, mixed-species reference set has no way to average past).

    Reuses line_identification.identify_lines_in_spectrum() -- the same
    detection-based centering used for EW-measurement line identification
    itself, rather than introducing a third line-centering method into
    the package.

    Runs in TWO passes, the way a standard coarse-then-fine cross-
    correlation RV search does, rather than one single wide search:
    1. COARSE: `search_radius` (default 1.0 A, deliberately much wider
       than identify_lines()'s own 0.15 A EW-identification default),
       since this runs before any wavelength correction and the true
       center can be offset by however large the spectrum's real,
       uncorrected RV is.
    2. FINE: the coarse pass's own median velocity is applied as a
       trial correction, then the search is repeated with a much
       narrower `refine_radius` (default 0.15 A, matching
       identify_lines()'s own default).

    The second pass is not just extra precision -- confirmed on real,
    densely-lined MAROON-X data that it's necessary for correctness: a
    wide search radius in a densely-lined spectrum can lock onto an
    unrelated NEARBY real line instead of the intended one often enough
    to measurably BIAS the coarse pass's own median velocity, not just
    add scatter around the right answer (confirmed directly: individual
    per-line velocities from the coarse pass showed no smooth trend with
    wavelength or order, including wide swings WITHIN a single order --
    inconsistent with a real RV or calibration drift, consistent with
    scattered misidentification in a dense line forest). Once a first-
    pass correction removes most of the true offset, a narrow window is
    both sufficient and far less likely to contain more than one real
    candidate.

    Parameters
    ----------
    lines : array of rest wavelengths (e.g. Spectrum_Data.lines).
    wavelength, flux, err, pred : lists of per-order arrays, UNSHIFTED
        (e.g. Spectrum_Data.wavelength, not shifted_wavelength -- this
        measures the shift that hasn't been applied yet).
    search_radius, refine_radius : Angstrom half-widths for the coarse
        and fine passes respectively -- see above.
    min_significance : passed through to identify_line() -- how many
        local-noise-sigma below continuum a candidate must clear to
        count, same meaning as there (identify_lines_in_spectrum()'s own
        empirical noise calibration applies here too, so this threshold
        means the same thing regardless of an instrument's own error
        propagation quirks -- see its docstring).
    sigma_clip : reject individual lines' implied velocity more than
        this many standard deviations from the median, same convention
        as measure_effective_rv() -- much more effective here than on a
        small reference set, since a real linelist has dozens of lines
        to average over instead of 3-4.
    min_lines : refuse to report an RV from fewer than this many
        detected lines at EACH pass (a handful of detections isn't
        enough to trust over a single bad blend/misidentification); if
        the fine pass comes up short, falls back to the coarse-only
        estimate rather than failing outright.

    Returns
    -------
    rv : float or None (km/s) -- None if the coarse pass found fewer
        than min_lines usable lines.
    rv_err : float (km/s) -- standard error on the mean of the lines the
        FINAL estimate is based on (fine pass if it had enough lines,
        else the coarse pass).
    n_used : int -- number of lines the reported RV is averaged over.
    """
    def _clip(v):
        if len(v) == 0:
            return v
        med, std = np.median(v), v.std()
        if std > 0:
            good = np.abs(v - med) < sigma_clip * std
            if good.sum() > 0:
                return v[good]
        return v

    coarse_v = _clip(_rv_from_detections(lines, wavelength, flux, err, pred,
                                          search_radius, min_significance))
    if len(coarse_v) < min_lines:
        return None, None, len(coarse_v)
    rv0 = float(np.median(coarse_v))  # median, not mean: robust to a coarse
                                       # pass skewed by misidentification, not
                                       # just noisy around the right answer

    shifted = [w * (1.0 - rv0 / C_KMS) for w in wavelength]
    fine_v = _clip(_rv_from_detections(lines, shifted, flux, err, pred,
                                        refine_radius, min_significance))
    if len(fine_v) < min_lines:
        rv_err = float(coarse_v.std() / np.sqrt(len(coarse_v))) if len(coarse_v) > 1 else 0.0
        return rv0, rv_err, len(coarse_v)

    rv = rv0 + float(np.mean(fine_v))
    rv_err = float(fine_v.std() / np.sqrt(len(fine_v))) if len(fine_v) > 1 else 0.0
    return rv, rv_err, len(fine_v)
