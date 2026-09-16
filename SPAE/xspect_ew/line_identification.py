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
                   min_significance=3.0, min_prominence=None, position_tolerance=0.07,
                   smooth_points=3, blend_significance_fraction=0.5, context_radius=None):
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
    search_radius : Angstrom half-width of the window a candidate's
        POSITION must fall within to be accepted. 0.15 sits comfortably
        above position_tolerance (0.07 A, the RV-corrected real-data
        position precision this package has measured elsewhere -- see
        DEVELOPMENT_LOG.md), so a correctly-shifted spectrum's true line
        essentially always falls well inside it, while still being
        tight enough that a moderately-dense real spectrum's NEXT line
        over usually falls outside it.
    context_radius : Angstrom half-width of the (wider) window the
        significance profile/prominence are actually COMPUTED over.
        Defaults to 2*search_radius. Kept separate from search_radius
        because the two need different sizes: candidate POSITIONS must
        stay tightly restricted to search_radius (so a real nearby-but-
        different line is never mistaken for the intended one), but
        computing prominence from that same narrow window is unreliable
        -- confirmed on two distinct real GRACES failure modes needing
        opposite-looking fixes that this one change resolves together.
        (1) A real, correctly-centered line can be UNDER-counted: a
        separate, shallower dip just outside a narrow window keeps that
        window's own edge significance elevated, artificially capping
        the correctly-centered peak's prominence below threshold even
        though its raw significance clears min_significance easily
        (confirmed: Fe I 7114.549, GRACES order 10 -- edge significance
        8.8 next to a legitimate interior peak at 10.5, prominence only
        1.7 against a 3.0 floor purely because the narrow window never
        shows the point further out where it truly returns to near
        baseline). (2) A fake, non-isolated "line" can be OVER-counted:
        an order-edge continuum-normalization artifact or a genuinely
        different, much stronger nearby feature can fill an entire
        narrow window with sustained elevated significance, giving high
        apparent prominence simply because the window is too narrow to
        show it never actually returns toward baseline (confirmed: Fe I
        6392.535 and 6745.090, GRACES -- an order sitting near ITS OWN
        edge, or containing a genuinely different stronger absorption a
        fraction of an Angstrom away, won cross-order arbitration in
        identify_lines_in_spectrum() purely because its narrow-window
        prominence looked fine; with prominence computed over the wider
        context instead, both no longer clear the threshold at all, so
        the correct order's candidate wins by default -- no separate
        cross-order fix needed). Candidate POSITIONS found outside
        search_radius (but inside context_radius) are computed for
        prominence purposes only and never returned as a detection.
    min_significance : minimum "sigma below continuum" for a candidate
        to count as a real detection, not noise. 3.0 is a standard
        detection threshold; real, usable absorption lines in a
        decent-S/N spectrum clear this by a wide margin (see
        DEVELOPMENT_LOG.md for real-data significance values found
        during validation).
    min_prominence : minimum topographic prominence (scipy.signal
        convention: height above the higher of the two valleys
        flanking a peak) the winning candidate must have, in the same
        sigma units as min_significance. Defaults to min_significance
        itself if not given -- a real, isolated detection should stand
        out from ITS OWN local surroundings by roughly as much as its
        absolute height, not just clear a fixed floor. Without this,
        a point sitting on the monotonic wing of a much deeper,
        DIFFERENT nearby line can register as a "peak" purely because
        it is (very marginally) higher than its immediate neighbor,
        with near-zero real prominence. Confirmed on real MAROON-X
        data (user-caught): a candidate at 3.03 sigma had prominence
        0.02 -- essentially the second-to-last point of a smoothly
        rising slope into a real, much stronger line nearby, not a
        genuine local feature at all. See context_radius above for why
        this alone isn't sufficient without also widening the window
        the prominence itself is computed over.
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
        n_candidates : int -- number of candidates clearing both
            min_significance and min_prominence AND positioned within
            search_radius (0 if none)
    """
    if min_prominence is None:
        min_prominence = min_significance
    if context_radius is None:
        context_radius = 2.0 * search_radius

    mask = (wave >= rest_wave - context_radius) & (wave <= rest_wave + context_radius)
    if mask.sum() < max(smooth_points, 3):
        return {'detected': False, 'center': np.nan, 'significance': 0.0,
                'blended': False, 'n_candidates': 0}

    w = wave[mask]
    depth = 1.0 - flux[mask]
    local_err = err[mask] / pred[mask]
    local_err = np.where(local_err > 0, local_err, np.inf)
    significance = uniform_filter1d(depth / local_err, size=smooth_points, mode='nearest')

    peak_idx, _ = find_peaks(significance, height=min_significance, prominence=min_prominence)
    if len(peak_idx) == 0:
        return {'detected': False, 'center': np.nan, 'significance': 0.0,
                'blended': False, 'n_candidates': 0}

    peak_wave = w[peak_idx]
    peak_sig = significance[peak_idx]

    # prominence/significance were computed over the wider context window
    # (see context_radius above), but only a candidate actually POSITIONED
    # within position_tolerance may be reported as the line itself -- NOT
    # the full (wider) search_radius. search_radius exists to find
    # candidates despite some real residual uncorrected offset; but
    # position_tolerance is this package's own measured RV-corrected
    # real-data precision (max |offset| 0.054-0.058 A across every
    # confirmed-good Keck/GRACES detection), so a lone candidate beyond
    # it -- with no competing, better-positioned alternative to weigh it
    # against -- is more likely a different, coincidentally nearby real
    # feature than the catalogued line itself. Confirmed on real
    # MAROON-X data: without this, several implausible detections passed
    # with offsets of 0.08-0.14 A and statistically absurd significance
    # for their catalogued strength (e.g. 115 sigma for a 10.5 mA line).
    in_window = np.abs(peak_wave - rest_wave) <= min(search_radius, position_tolerance)
    if not in_window.any():
        return {'detected': False, 'center': np.nan, 'significance': 0.0,
                'blended': False, 'n_candidates': 0}
    peak_wave = peak_wave[in_window]
    peak_sig = peak_sig[in_window]

    score = peak_sig / (1.0 + ((peak_wave - rest_wave) / position_tolerance)**2)

    best = np.argmax(score)
    center = float(peak_wave[best])
    best_sig = float(peak_sig[best])

    others = np.delete(peak_sig, best)
    blended = bool(len(others) > 0 and others.max() >= blend_significance_fraction * best_sig)

    return {'detected': True, 'center': center, 'significance': best_sig,
            'blended': blended, 'n_candidates': int(len(peak_wave))}


def _empirical_noise_calibration(flux, pred, err):
    """
    Rescale one order's err to match its ACTUAL above-fit residual
    spread, rather than trusting err's theoretical (Poisson-style)
    scale blindly -- the same fix as continuum.py's fit_als_continuum()
    (target_percentile's noise calibration), independently re-derived
    here because it matters for line DETECTION too, not just continuum
    placement. Confirmed on real MAROON-X data: without this, real
    lines' computed significance is deflated by the same ~0.25-0.4
    factor found there ("optimal extraction" pipelines that combine
    multiple raw CCD pixels per output point correlate adjacent points,
    giving less real point-to-point scatter than err's Poisson-style
    scaling predicts), pushing many real lines below the detection
    threshold that would otherwise clear it easily -- confirmed
    directly: this alone took one real MAROON-X order set's detection
    count from 13/78 to 27/78 lines, with the calibration factor
    (median ~0.34) matching the continuum-fitting context's
    independently-derived value almost exactly. A no-op on Keck/GRACES
    (confirmed: full detection preserved, calibration factor near 1).

    `flux` here is NORMALIZED flux (~1.0 baseline), this module's
    convention throughout -- NOT raw flux the way continuum.py's
    fit_als_continuum() uses it (resid = flux - pred there). The
    equivalent raw-scale residual from normalized flux is
    pred*(flux-1.0), since normalized_flux = raw_flux/pred by
    definition; `err` is still expected in raw/absolute units
    (Spectrum_Data.obs_err), matching `pred`'s scale, since that's what
    identify_line() itself expects (it divides by pred internally).
    """
    resid = pred * (flux - 1.0)
    above = resid > 0
    if above.sum() < 10:
        return 1.0
    empirical = np.median(resid[above])
    theoretical = np.median(err[above]) * 0.6744897501960817  # median of a positive half-normal
    if theoretical <= 0:
        return 1.0
    return float(np.clip(empirical / theoretical, 0.05, 3.0))


# Cross-order arbitration below prefers the candidate positioned CLOSER
# to rest_wave over one with merely higher significance; two offsets
# within this many Angstrom of each other are treated as a tie (broken
# by significance instead). Set well below position_tolerance's own
# default (0.07 A) -- this only needs to separate genuine sub-noise
# centering differences (typically <0.005 A between two good detections
# of the same real line -- see the two order-13/-14, -11/-12 examples in
# identify_lines_in_spectrum()'s docstring) from a candidate that's
# measuring something else entirely (offsets of several hundredths of
# an Angstrom in the confirmed real cases below).
_CROSS_ORDER_POSITION_TIE = 0.005


def identify_lines_in_spectrum(lines, wavelength, flux, err, pred, calibrate_noise=True, **kwargs):
    """
    Run identify_line() for every rest wavelength in `lines`, against
    every order whose wavelength range could contain it -- i.e. every
    order for which the search window (rest_wave +/- search_radius)
    overlaps that order's own coverage, not just orders whose exact
    range contains rest_wave itself, so a line right at one order's
    edge still gets a fair chance via a neighboring, overlapping order.

    When more than one order detects the line, keeps whichever is
    positioned CLOSER to rest_wave (ties within _CROSS_ORDER_POSITION_TIE
    broken by higher significance) -- NOT whichever has the highest raw
    significance. Two overlapping orders are two independent
    measurements of the same true spectrum; if they disagree on WHERE
    the line is by more than noise can explain, that disagreement itself
    is the signal that one of them is looking at something else, and
    precise centering (not depth) is what distinguishes "the correct,
    catalogued line" from "a different, coincidentally nearby feature".

    Confirmed on two real GRACES cases where raw-significance selection
    picked the wrong order outright: (1) Fe I 6392.535 -- order 13 (57 A
    interior) finds a correctly-centered, appropriately-shallow (~6%)
    dip essentially exactly at rest_wave (21.6 sigma), while order 14,
    whose own order boundary sits only ~2.4 A past this wavelength,
    shows a much deeper (76.8 sigma) but 0.046 A offset dip that turns
    out to be a continuum-normalization artifact from that order's own
    edge (confirmed separately: order 14's flux declines monotonically,
    never recovering, all the way to its literal last data point).
    (2) Fe I 6745.090 -- order 12 finds a correctly-centered, appropriately
    weak (8.1 sigma) dip essentially exactly at rest_wave, while order 11
    shows a real, well-formed but 0.033 A offset, much deeper (18.9
    sigma) feature -- too deep for this 8.1 mA line, more likely a
    genuinely different absorption feature nearby. Raw-significance
    selection picked the wrong (deeper, mispositioned) order in both
    cases; position-based selection picks correctly in both.

    Parameters
    ----------
    lines : array of rest wavelengths.
    wavelength, flux, err, pred : lists of per-order arrays (e.g.
        Spectrum_Data's shifted_wavelength, normalized_flux, obs_err,
        pred_all).
    calibrate_noise : if True (default), rescale each order's err by
        _empirical_noise_calibration() before computing significance --
        see there for why. Set False to use err exactly as given.
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
    if calibrate_noise:
        err = [err[o] * _empirical_noise_calibration(flux[o], pred[o], err[o])
               for o in range(len(wavelength))]
    results = []
    for rest_wave in lines:
        best_result = {'detected': False, 'center': np.nan, 'significance': 0.0,
                       'blended': False, 'n_candidates': 0}
        best_order = None
        best_offset = np.inf
        for order in range(len(wavelength)):
            w = wavelength[order]
            if w[-1] < rest_wave - search_radius or w[0] > rest_wave + search_radius:
                continue
            result = identify_line(rest_wave, w, flux[order], err[order], pred[order], **kwargs)
            if not result['detected']:
                continue
            offset = abs(result['center'] - rest_wave)
            is_tie = abs(offset - best_offset) <= _CROSS_ORDER_POSITION_TIE
            better = (best_order is None or
                      offset < best_offset - _CROSS_ORDER_POSITION_TIE or
                      (is_tie and result['significance'] > best_result['significance']))
            if better:
                best_result = result
                best_order = order
                best_offset = offset
        best_result['order'] = best_order
        results.append(best_result)
    return results
