"""Line identification: locate each linelist entry in a normalized,
wavelength/RV-corrected spectrum, BEFORE any equivalent-width
measurement is attempted.

Deliberately kept separate from line_profile.py's existing
get_line_window()/measure_ew() machinery (used for actual EW
measurement, which may be replaced independently -- see
DEVELOPMENT_LOG.md). This module answers a narrower, prior question:
given a line's rest wavelength, IS there a real absorption feature near
it in this spectrum at all, and if so, where exactly is its center --
not how deep/wide it is or what its EW is.

The previous (and still-used, by get_line_window()) approach was
"whatever the single lowest flux point within +/-0.1 A happens to be,
call that the line" -- with no test for whether a real feature is even
present. That always reports SOMETHING, including pure noise, and a
downstream check (Spectrum_Data.check_for_flags()'s position-offset
test) can only catch a resulting misidentification after the fact, by
noticing the reported center strayed too far from the rest wavelength.
This module instead tests for a real detection up front: it converts
the local window to a per-point "how many local-noise-sigma below
continuum" significance profile (using each point's own propagated
error, already correctly instrument/response-corrected -- see
continuum.py's module docstring for why this can't just be assumed
uniform), finds significant local minima in THAT profile via
scipy.signal.find_peaks, and only reports a detection if at least one
candidate clears a minimum significance. A line with no real,
significant absorption feature nearby (too weak for this spectrum's
S/N, or genuinely absent) is reported as NOT DETECTED rather than
silently handed noise as if it were a real position.

Candidates are scored by significance discounted by distance from the
rest wavelength (not simply "deepest wins"): a moderately significant
dip right at the expected position is more likely the intended line
than a more significant but more distant one, which is more likely a
different, neighboring real line. A second candidate close enough to
compete with the winner is recorded as a blend flag, since its presence
means the measured region isn't a clean, isolated absorption feature
even if a center was still confidently identified.

Known limitation, not addressed here: this only flags a blend when
find_peaks() resolves TWO separate local maxima in the significance
profile. Two lines close enough together that their combined profile
has no resolvable valley between them (confirmed with a synthetic
test: two comparable-depth lines 0.09 A apart, each FWHM 0.15 A) are
reported as a single, correctly-centered detection with no blend flag
-- not wrong (the reported center is still accurate), but understates
that the region isn't a single isolated line. Flagging this class would
need a width-based check (an unusually wide detected feature relative
to what a single line normally looks like at this resolution), which
needs a real per-instrument "normal single-line width" reference this
package does not yet establish anywhere -- left as a known gap rather
than guessing at an unjustified absolute threshold.
"""

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks


def identify_line(rest_wave, wave, flux, err, pred, search_radius=0.15,
                   min_significance=3.0, position_tolerance=0.07,
                   smooth_points=3, blend_significance_fraction=0.5):
    """
    Locate one line's real center in one order, or report it as not
    detected.

    Parameters
    ----------
    rest_wave : the line's rest (lab) wavelength, Angstrom.
    wave, flux : this order's shifted_wavelength, normalized_flux.
    err, pred : this order's obs_err, pred_all (continuum) -- used
        together as err/pred, the same normalized-flux-scale noise
        estimate already used elsewhere (e.g. measure_ew()'s error
        bars), so the significance profile below is correctly informed
        by each point's own noise, not a single assumed scale.
    search_radius : Angstrom half-width of the window searched around
        rest_wave. 0.15 sits comfortably above position_tolerance
        (0.07 A, the RV-corrected real-data position precision this
        package has measured elsewhere -- see DEVELOPMENT_LOG.md), so a
        correctly-shifted spectrum's true line essentially always falls
        well inside the search window, while still being tight enough
        that a moderately-dense real spectrum's NEXT line over usually
        falls outside it.
    min_significance : minimum "sigma below continuum" for a candidate
        to count as a real detection, not noise. 3.0 is a standard
        detection threshold; real, usable absorption lines in a
        decent-S/N spectrum clear this by a wide margin (see
        DEVELOPMENT_LOG.md for real-data significance values found
        during validation).
    position_tolerance : Angstrom scale over which a candidate's score
        is discounted with distance from rest_wave (see module
        docstring -- NOT a hard cutoff, a soft preference).
    smooth_points : light smoothing applied to the significance profile
        before peak-finding, so ordinary point-to-point noise near the
        bottom of a real line doesn't get counted as multiple separate
        candidates.
    blend_significance_fraction : a second candidate is flagged as a
        blend risk if its significance is at least this fraction of the
        winning candidate's.

    Returns
    -------
    dict with keys:
        detected : bool
        center : float (Angstrom) or nan if not detected
        significance : float, the winning candidate's sigma-below-
            continuum (0 if not detected)
        blended : bool -- a competitive second candidate exists nearby
        n_candidates : int -- number of candidates clearing
            min_significance in the search window (0 if none)
    """
    mask = (wave >= rest_wave - search_radius) & (wave <= rest_wave + search_radius)
    if mask.sum() < max(smooth_points, 3):
        return {'detected': False, 'center': np.nan, 'significance': 0.0,
                'blended': False, 'n_candidates': 0}

    w = wave[mask]
    depth = 1.0 - flux[mask]
    local_err = err[mask] / pred[mask]
    local_err = np.where(local_err > 0, local_err, np.inf)
    significance = uniform_filter1d(depth / local_err, size=smooth_points, mode='nearest')

    peak_idx, _ = find_peaks(significance, height=min_significance)
    if len(peak_idx) == 0:
        return {'detected': False, 'center': np.nan, 'significance': 0.0,
                'blended': False, 'n_candidates': 0}

    peak_wave = w[peak_idx]
    peak_sig = significance[peak_idx]
    score = peak_sig / (1.0 + ((peak_wave - rest_wave) / position_tolerance)**2)

    best = np.argmax(score)
    center = float(peak_wave[best])
    best_sig = float(peak_sig[best])

    others = np.delete(peak_sig, best)
    blended = bool(len(others) > 0 and others.max() >= blend_significance_fraction * best_sig)

    return {'detected': True, 'center': center, 'significance': best_sig,
            'blended': blended, 'n_candidates': int(len(peak_idx))}


def identify_lines_in_spectrum(lines, wavelength, flux, err, pred, **kwargs):
    """
    Run identify_line() for every rest wavelength in `lines`, against
    every order whose wavelength range could contain it -- i.e. every
    order for which the search window (rest_wave +/- search_radius)
    overlaps that order's own coverage, not just orders whose exact
    range contains rest_wave itself, so a line right at one order's
    edge still gets a fair chance via a neighboring, overlapping order.
    When more than one order is a candidate, keeps whichever gives the
    HIGHEST significance -- naturally prefers whichever order's local
    data/continuum quality is better in an overlap region, the same
    property overlap_check.py's cross-order comparison exploits.

    Parameters
    ----------
    lines : array of rest wavelengths.
    wavelength, flux, err, pred : lists of per-order arrays (e.g.
        Spectrum_Data's shifted_wavelength, normalized_flux, obs_err,
        pred_all).
    **kwargs : passed through to identify_line() (search_radius,
        min_significance, etc.)

    Returns
    -------
    list of dicts, one per line (same order as `lines`), each an
    identify_line() result dict plus 'order' (the order index the
    result came from, or None if never detected in any candidate
    order).
    """
    search_radius = kwargs.get('search_radius', 0.15)
    results = []
    for rest_wave in lines:
        best_result = {'detected': False, 'center': np.nan, 'significance': 0.0,
                       'blended': False, 'n_candidates': 0}
        best_order = None
        for order in range(len(wavelength)):
            w = wavelength[order]
            if w[-1] < rest_wave - search_radius or w[0] > rest_wave + search_radius:
                continue
            result = identify_line(rest_wave, w, flux[order], err[order], pred[order], **kwargs)
            if result['detected'] and result['significance'] > best_result['significance']:
                best_result = result
                best_order = order
        best_result['order'] = best_order
        results.append(best_result)
    return results
