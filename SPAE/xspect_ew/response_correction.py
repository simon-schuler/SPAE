"""
General instrument response/blaze correction.

Not tied to any one instrument's file format -- takes a response curve as
plain per-chunk (wave, response) arrays, from wherever it came from (an
embedded blaze table, a separate calibration file like MAROON-X's PHOENIX-
based response correction -- see readers.load_maroonx_response() for that
one -- or anything else), and divides it into a Spectrum_Data's flux.

Response chunks are matched to each of the spectrum's own orders by
wavelength overlap, not by any order-numbering convention -- those aren't
guaranteed to agree between the science spectrum and wherever the response
curve was sourced from (confirmed necessary: MAROON-X's response file uses
real echelle order numbers as column labels, e.g. 92, 93, ..., which have
no reason to line up with array positions in an arbitrarily-read science
spectrum). Matching one response chunk per order (rather than
concatenating all chunks into one global curve and interpolating across
chunk boundaries) avoids introducing spurious discontinuities at order
edges, since per-order blaze shape can differ in absolute throughput even
between adjacent orders.

Wavelength overlap alone is NOT always enough to disambiguate, though:
confirmed on real MAROON-X data, its two arms (blue/red) physically
overlap near their dichroic split (both arms independently cover echelle
orders 91-94), so a science order there can have TWO response chunks
that both "contain" its mean wavelength, describing two different light
paths with nearly identical wavelength ranges but very different
absolute response curves. Picking the first match found (this function's
original behavior) silently applies the wrong arm's calibration whenever
that happens -- confirmed to produce a spurious, smoothly WRONG ~4x
monotonic trend across an entire order (not just an edge artifact,
because the mismatch is in the whole chunk's shape, not one bad pixel).
When `science_bands`/`response_bands` are supplied, candidates are
filtered to same-arm matches before picking one, resolving the ambiguity
directly instead of guessing; see readers.get_maroonx_bands() for how to
obtain `science_bands` for a MAROON-X file.
"""

import numpy as np
from scipy.interpolate import interp1d


def apply_response_correction(spectrum, response_wave, response, min_overlap_fraction=0.5,
                               min_response_fraction=0.1, response_bands=None, science_bands=None):
    """
    Divide a response/blaze correction curve into spectrum.flux, per order.

    Run this BEFORE normalize() -- it corrects the raw counts (flux), not
    normalized_flux, and normalize()'s continuum fit will be far more
    robust on an already-flattened spectrum.

    Parameters
    ----------
    spectrum : Spectrum_Data
    response_wave, response : lists of per-chunk arrays (Angstrom, unitless
        response/throughput). Any number of chunks; doesn't need to match
        this spectrum's own order count or boundaries.
    min_overlap_fraction : float -- if the matched response chunk covers
        less than this fraction of an order's wavelength range, that order
        is left uncorrected (with a message) rather than extrapolating the
        response curve into territory it doesn't actually describe.
    min_response_fraction : float -- points where the response curve has
        fallen below this fraction of ITS OWN chunk's peak are not divided
        directly; their corrected value is instead INTERPOLATED from the
        nearest above-floor points elsewhere in the same order. Needed
        because a blaze/response curve tapers toward zero at order edges
        (and sometimes elsewhere) by construction -- dividing raw counts
        (and their noise) by a near-zero response inflates a handful of
        pixels far more than their neighbors (confirmed: one real
        MAROON-X order's lowest-throughput edge pixel, response 0.0117%
        of that chunk's median, amplified 570 raw counts into ~4.9e6
        "counts" -- 1000x+ its already-corrected neighbors). Falling back
        to RAW flux at those pixels (an earlier version of this function)
        is actually worse: raw counts then sit ~1000x BELOW the
        surrounding corrected values, which looks exactly like a sharp
        absorption line at that pixel -- precisely the line-
        misidentification risk this whole package exists to avoid.
        Interpolating from neighboring above-floor points bridges the gap
        smoothly instead, at the cost of a few pixels' worth of
        independent information (an order typically overlaps its
        neighbors in wavelength anyway).
    response_bands, science_bands : optional, parallel labels (any
        hashable, e.g. 'blue'/'red') for response_wave/response and for
        spectrum's own orders respectively. When both are given, a
        response chunk is only considered a candidate match for an order
        if their labels are equal -- see module docstring for why this
        matters on real MAROON-X data. When either is None (default),
        matching is wavelength-overlap-only, unchanged from before.

    Returns
    -------
    corrected_orders : list of order indices that were actually corrected
    """
    corrected_orders = []
    for i in range(len(spectrum.wavelength)):
        w = spectrum.wavelength[i]
        flux = spectrum.flux[i]
        order_mean = w.mean()

        if spectrum.continuum[i].any():
            print(f'order {i}: normalize() already ran on this order -- its continuum fit and '
                  f'obs_err were computed from the PRE-correction flux and are now stale. '
                  f're-run normalize() for this order after response correction.')

        # find the response chunk(s) whose range contains this order's mean
        # wavelength (same "mean falls inside" matching principle used
        # elsewhere in this package, e.g. estimate_shift()), then narrow to
        # a same-arm match if band labels were supplied (see module
        # docstring -- wavelength overlap alone is ambiguous for MAROON-X's
        # two physically-overlapping arms)
        band_iter = response_bands if response_bands is not None else [None] * len(response_wave)
        candidates = [(rw, rf, b) for rw, rf, b in zip(response_wave, response, band_iter)
                      if rw.min() <= order_mean <= rw.max()]
        if science_bands is not None and response_bands is not None:
            same_band = [c for c in candidates if c[2] == science_bands[i]]
            if same_band:
                candidates = same_band
        candidates = [(rw, rf) for rw, rf, b in candidates]
        best_chunk = candidates[0] if candidates else None
        if best_chunk is None:
            print(f'order {i} (mean {order_mean:.1f} A): no response-curve chunk covers this '
                  f'wavelength, left uncorrected')
            continue

        rw, rf = best_chunk
        overlap_lo = max(rw.min(), w.min())
        overlap_hi = min(rw.max(), w.max())
        covered_frac = (overlap_hi - overlap_lo) / (w.max() - w.min())
        if covered_frac < min_overlap_fraction:
            print(f'order {i}: matched response chunk only covers {covered_frac:.0%} of this '
                  f'order, left uncorrected (pass a lower min_overlap_fraction to force)')
            continue

        interp_resp = interp1d(rw, rf, kind='linear', bounds_error=False, fill_value=np.nan)
        resp_on_grid = interp_resp(w)
        response_floor = min_response_fraction * rf.max()
        valid = ~np.isnan(resp_on_grid) & (resp_on_grid > response_floor)

        if valid.sum() < 2:
            print(f'order {i}: fewer than 2 points above the response floor, left uncorrected')
            continue

        corrected = flux.copy()
        corrected[valid] = flux[valid] / resp_on_grid[valid]

        # obs_err must be divided by the SAME response curve as flux, not
        # just recomputed as sqrt(corrected) downstream -- confirmed real
        # bug: sqrt(raw_counts)/response is the correct propagated Poisson
        # error on response-corrected flux, but sqrt(corrected_flux) alone
        # silently drops the 1/response factor. Since response tapers
        # non-uniformly across an order (blaze-like shape, not flat), this
        # under-estimates the true error MORE in low-response stretches
        # than high-response ones, handing fit_als_continuum's below-fit
        # weighting (1/err^2) artificially high confidence there and
        # letting ordinary noise excursions get chased as if they were
        # trustworthy signal. Confirmed on a real MAROON-X order (48):
        # this alone produced a smooth ~19-percentage-point fitted-
        # continuum overshoot across an 83-A "clean" stretch with no
        # absorption at all, worst at the order's true edge (lowest
        # relative response) and fading toward the order's response peak
        # -- not a telluric- or line-density-driven effect, since a
        # neighboring order (56) with comparably deep telluric absorption
        # but a flatter relative-response profile in ITS clean region
        # showed no such gradient.
        obs_err = spectrum.obs_err[i]
        corrected_err = obs_err.copy()
        corrected_err[valid] = obs_err[valid] / resp_on_grid[valid]

        n_bridged = (~valid).sum()
        if n_bridged:
            # w is wavelength-sorted, so w[valid]/corrected[valid] are too --
            # the first/last valid entries are the correct edge fill values
            bridge = interp1d(w[valid], corrected[valid], kind='linear', bounds_error=False,
                               fill_value=(corrected[valid][0], corrected[valid][-1]))
            corrected[~valid] = bridge(w[~valid])
            err_bridge = interp1d(w[valid], corrected_err[valid], kind='linear', bounds_error=False,
                                   fill_value=(corrected_err[valid][0], corrected_err[valid][-1]))
            corrected_err[~valid] = err_bridge(w[~valid])
            print(f'order {i}: {n_bridged}/{len(w)} points below the response floor '
                  f'({min_response_fraction:.0%} of chunk peak), bridged via interpolation '
                  f'from neighboring corrected points')

        spectrum.flux[i] = corrected
        spectrum.obs_err[i] = corrected_err
        corrected_orders.append(i)

    return corrected_orders
