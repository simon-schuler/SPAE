"""Data-quality flagging for a spectrum's own raw data (as opposed to
overlap_check.py, which compares two DIFFERENT orders' agreement): two
checks at two different scales.

check_bad_orders() flags a whole order whose raw flux contains an
extreme, non-astrophysical outlier (e.g. a cosmic ray or a fixed
detector/amplifier-boundary defect) that would corrupt any continuum fit
run on it -- rather than trying to make the continuum-fitting routine
itself robust to every possible corruption. Same flag-and-exclude
philosophy as overlap_check.py, just for a different failure mode.

detect_spikes() flags/locates individual extreme pixels within an
otherwise-good order (e.g. for Spectrum_Data.correct_spikes()/
combine_spectra() to actually replace, not just flag out the whole
order) -- see its own docstring for why this needs a different,
factor-based (not sigma-based) test than a naive per-pixel noise check.

Confirmed on a real GRACES order (reproducible at the same wavelength/
pixel index across two independent exposures of the same star taken 3
days apart -- a fixed defect, not random cosmic-ray noise, though this
check doesn't need to tell the two apart): a single point at ~700x the
order's own 99th-percentile flux corrupted the ENTIRE order's AsLS
continuum fit via its smoothness-penalty coupling, not just the affected
point (normalized_flux median crashed to ~0.055 instead of a reasonable
~0.9-1.0).

The check compares the order's single highest point against its OWN
99th-percentile flux, not a median-centered/symmetric measure (e.g.
MAD): a first attempt using median absolute deviation false-flagged 11
of 35 real orders whose only "outlier" was a genuine, deep, real
spectral feature (Na D, H-alpha, an O2 telluric band) sitting in an
otherwise very low-scatter stretch -- deep ABSORPTION already has a
purpose-built, asymmetric handling mechanism in the continuum-fitting
routine itself (see continuum.py), so this check only needs to catch
the failure mode that mechanism can't: an extreme point ABOVE the
order's own natural peak level. Confirmed directly against every real
spectrum in this project's test set: Keck, GRACES (both formats),
MAROON-X all sit at max/p99 <= 1.51, while the real defect sits at
96.7-713 -- a >60x gap between the two clusters.
"""

import numpy as np
from scipy.ndimage import median_filter
from scipy.interpolate import interp1d

from .line_identification import _empirical_noise_calibration


def check_bad_orders(spectrum, outlier_factor=20.0):
    """
    Flag whole orders whose RAW flux (spectrum.flux, not the fitted
    continuum -- so this can run before normalize_all() too, catching the
    problem at its source rather than only after it has already corrupted
    a fit) has a maximum value more than `outlier_factor` times its own
    99th-percentile flux.

    Parameters
    ----------
    spectrum : Spectrum_Data
    outlier_factor : default 20.0 sits in the middle (in log space) of
        the gap between every real order in this project's test set
        (max/p99 <= 1.51) and the confirmed real defect (max/p99 =
        96.7-713) -- wide margin either way. See module docstring for
        why this is max-vs-p99, not a symmetric/median-centered measure.

    Returns
    -------
    list of dicts, one per flagged order: {'order', 'wave_lo', 'wave_hi',
    'outlier_factor'} -- same wave_lo/wave_hi key convention as
    overlap_check.py's flagged ranges, so check_for_flags() can share one
    loop shape across both checks.
    """
    flagged = []
    for i in range(len(spectrum.flux)):
        flux = np.asarray(spectrum.flux[i], dtype=float)
        p99 = np.percentile(flux, 99)
        if p99 <= 0:
            continue
        worst = flux.max() / p99
        if worst > outlier_factor:
            wave = spectrum.wavelength[i]
            flagged.append({'order': i, 'wave_lo': float(wave.min()),
                             'wave_hi': float(wave.max()),
                             'outlier_factor': float(worst)})
    return flagged


def detect_spikes(flux, window=5, factor=5.0):
    """
    Detect isolated, extreme single/few-pixel POSITIVE spikes (cosmic
    rays, hot pixels, uncorrected sky-emission-line contamination) in one
    order's raw flux, usable on a SINGLE spectrum with no second exposure
    to cross-check against (see Spectrum_Data.combine_spectra() for a
    more sensitive two-exposure comparison when a second exposure IS
    available).

    Compares each point to a TIGHT local median (default 5 points) by
    RATIO, not a sigma/noise-based measure. Two sigma-based attempts were
    tried and rejected on real data before this: (1) comparing to a local
    median-absolute-deviation badly UNDER-estimated the true local noise
    inside any real large-scale slope/trend, since a short window
    straddling real curvature has artificially small point-to-point MAD
    -- false-flagged smooth, ordinary continuum points at ~15-25x their
    true significance; (2) comparing to the spectrum's own theoretical
    Poisson error (obs_err) fared even worse: real, genuinely RESOLVED
    spectral structure (e.g. the local continuum peak between two
    blended absorption lines, or a real sky-emission line) routinely
    varies pixel-to-pixel by MANY sigma of pure photon noise over just a
    few pixels -- confirmed on real Keck data (a smooth, symmetric,
    several-pixel-wide real continuum bump, indistinguishable in
    amplitude from noise-based statistics) and on the real GRACES sky-
    emission-line case from combine_spectra()'s testing (a genuinely
    resolved emission feature, not a defect). A sigma-based test cannot
    tell these apart from a genuine defect.

    A simple factor-based ratio can, because the two populations turned
    out to occupy completely different regimes on real data: every
    genuine resolved feature checked (real Keck spectral structure, the
    real GRACES sky-emission line) stays under ~1.6x its own tight local
    median, while a confirmed real defect (the same one motivating
    check_bad_orders()) reaches 5.5-126x. `factor` sits well inside that
    gap. This deliberately only catches EXTREME, unambiguous cases (the
    user's own framing: "strong spikes") -- a moderate, ambiguous
    excursion that could plausibly be real, resolved structure (like
    that same sky-emission line) is left for the two-exposure comparison
    in combine_spectra(), which can resolve the ambiguity with a second,
    independent measurement instead of guessing from shape alone.

    Parameters
    ----------
    flux : one order's raw flux array.
    window : local median window (points). 5 is tight enough to track
        real curvature closely (confirmed: real Keck structure varies
        smoothly enough that a 5-point window's residual is small) while
        still being disrupted only by a genuinely narrow defect.
    factor : flux/local_median ratio a point must exceed to be flagged.

    Returns
    -------
    bad : boolean array, same shape as flux.
    local_median : the array used as the reference/replacement value
        (so callers needing to correct the spectrum don't recompute it).
    """
    flux = np.asarray(flux, dtype=float)
    local_median = median_filter(flux, size=window, mode='reflect')
    ratio = np.divide(flux, local_median, out=np.ones_like(flux), where=local_median > 0)
    return ratio > factor, local_median


def check_cross_exposure_spikes(wave_A, nflux_A, pred_A, err_A,
                                 wave_B, nflux_B, pred_B, err_B,
                                 residual=0.0, significance=8.0):
    """
    Detect an extreme, non-astrophysical EXCESS (cosmic ray, hot pixel,
    uncorrected sky-emission-line contamination) in EITHER of two
    exposures of the same star, at a given aligned wavelength --
    substantially more sensitive than detect_spikes() alone, since
    comparing two independent measurements of the same true signal can
    resolve cases a single spectrum's shape can't (see
    Spectrum_Data.combine_spectra()).

    Uses normalized flux exceeding the CONTINUUM (1.0), not a comparison
    between the two spectra's raw values or a local window -- a real
    absorption line can only ever push flux DOWN, so "significantly
    above 1.0" is an unambiguous, physically-motivated excess signature
    regardless of local line density. Both spectra's significance is
    computed independently (each spectrum's own empirically-calibrated
    relative error, see line_identification._empirical_noise_calibration()
    -- needed here as much as there: confirmed a real GRACES order's
    naive theoretical error underestimated the true noise by up to the
    calibration function's own 3x clip ceiling).

    A point is flagged whenever EITHER side's excess significance clears
    the threshold -- NOT only when one side is significant and the other
    isn't. An asymmetry-gated version was tried first and rejected on
    real data: a real, resolved sky-emission line (confirmed: the
    classic 5577 A and 6300 A [OI] auroral lines) is present in BOTH
    exposures, just at different strength (real night-to-night sky
    variability) -- requiring the OTHER side be negligible correctly
    avoided real line-core noise, but also left the CENTER of the
    contamination (where both sides happen to show comparable excess)
    uncorrected, catching only its weaker edges. Since real absorption
    can never trigger a positive-excess test in EITHER exposure (line
    cores show excess significance NEGATIVE, below continuum, in both),
    dropping the asymmetry requirement loses no real protection while
    fixing this: whenever either side is significant, the lower of the
    two is always the safer value, regardless of what the other side is
    doing.

    Parameters
    ----------
    wave_A, nflux_A, pred_A, err_A : spectrum A's FULL shifted_wavelength/
        normalized_flux/pred_all/obs_err for one order (the full order,
        not a pre-restricted subset -- _empirical_noise_calibration()
        needs the whole order's above-fit residuals to calibrate
        correctly).
    wave_B, nflux_B, pred_B, err_B : same for spectrum B's matched order.
    residual : Angstrom, the SAME fine-alignment residual (see
        combine.measure_order_alignment()) combine_spectra() applies to
        B's raw flux/error before combining -- passed here too so this
        check's own A/B pairing lines up with exactly the same points
        combine_spectra() actually combines, not a subtly different one.
    significance : how many (calibrated) sigma above continuum EITHER
        side must clear to flag that point.

    Returns
    -------
    in_range : boolean mask into wave_A -- which of A's points had a
        valid B counterpart to compare against (identical formula to
        combine_spectra()'s own in_range, given the same residual).
    bad : boolean array, same length as wave_A[in_range] -- flagged points.
    a_is_higher : boolean array, same length -- True where A's excess
        significance is the larger of the two, for flagged points only
        (meaningless where bad is False) -- the caller should use B's
        (lower, less-contaminated) value there, and vice versa.
    """
    cal_A = _empirical_noise_calibration(nflux_A, pred_A, err_A)
    cal_B = _empirical_noise_calibration(nflux_B, pred_B, err_B)
    relerr_A = np.where(pred_A > 0, (err_A * cal_A) / pred_A, np.inf)
    relerr_B = np.where(pred_B > 0, (err_B * cal_B) / pred_B, np.inf)
    sig_A_full = (nflux_A - 1.0) / relerr_A
    sig_B_full = (nflux_B - 1.0) / relerr_B

    in_range = (wave_A >= wave_B.min() + abs(residual)) & (wave_A <= wave_B.max() - abs(residual))
    query = wave_A[in_range] + residual
    sig_B = interp1d(wave_B, sig_B_full, kind='linear', bounds_error=False, fill_value=0.0)(query)
    sig_A = sig_A_full[in_range]

    bad = (sig_A > significance) | (sig_B > significance)
    a_is_higher = sig_A >= sig_B
    return in_range, bad, a_is_higher
