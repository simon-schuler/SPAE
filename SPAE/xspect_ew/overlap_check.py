"""Cross-order overlap consistency check for Spectrum_Data.

Real echelle orders commonly overlap in wavelength with a neighbor (by
physical design), and -- confirmed on real MAROON-X data -- sometimes
with a non-neighbor too (its two arms both cover echelle orders 91-94
near their dichroic split). Two overlapping orders are independent
measurements of the SAME wavelengths of the SAME star, near-
simultaneously: if both orders' continuum placement is correct, their
normalized flux should agree throughout the overlap WITHOUT needing to
know the true continuum level at all. This is exactly how the MAROON-X
arm-mismatch bug was originally found -- comparing order 2 against
order 60, both covering H-alpha, because the user happened to notice
they shared a line -- just done here systematically across every
overlapping pair instead of relying on a human noticing two orders
share a feature.

This can't catch a bias that affects every order the same way (two
equally-biased orders still agree with each other), and has nothing to
say about an order with no overlapping neighbor (typically the
reddest/bluest order in an instrument's coverage). It's a real,
complementary diagnostic to synthetic ground-truth testing, not a
replacement for it -- it validates INTERNAL CONSISTENCY between orders,
not absolute correctness against a known truth.
"""

import numpy as np
from scipy.interpolate import interp1d


def check_order_overlaps(spectrum, min_overlap_points=10):
    """
    Compare normalized_flux between every pair of orders whose
    wavelength ranges overlap.

    Parameters
    ----------
    spectrum : Spectrum_Data -- normalize_all() (and apply_rv_shift(),
        if used) should already have been run.
    min_overlap_points : skip pairs whose overlap has fewer than this
        many valid points (too few for the median comparison to mean
        much).

    Returns
    -------
    list of dicts, one per overlapping pair with enough points, sorted
    by |median_diff| (worst first). Each dict:
        order_i, order_j : order indices (i < j)
        wave_lo, wave_hi : overlap wavelength range (Angstrom)
        n_points : points compared
        median_diff : median(normalized_flux_i - normalized_flux_j),
            order_j's normalized_flux interpolated onto order i's grid
        median_pct : median_diff expressed as a percent
    """
    wave = spectrum.shifted_wavelength
    flux = spectrum.normalized_flux
    n = len(wave)
    bounds = [(w.min(), w.max()) for w in wave]

    results = []
    for i in range(n):
        lo_i, hi_i = bounds[i]
        for j in range(i + 1, n):
            lo_j, hi_j = bounds[j]
            overlap_lo, overlap_hi = max(lo_i, lo_j), min(hi_i, hi_j)
            if overlap_hi <= overlap_lo:
                continue

            wi, fi = wave[i], flux[i]
            mask_i = (wi >= overlap_lo) & (wi <= overlap_hi)
            if mask_i.sum() < min_overlap_points:
                continue

            wj, fj = wave[j], flux[j]
            interp_j = interp1d(wj, fj, bounds_error=False, fill_value=np.nan)
            fj_on_i = interp_j(wi[mask_i])
            valid = ~np.isnan(fj_on_i)
            if valid.sum() < min_overlap_points:
                continue

            diff = fi[mask_i][valid] - fj_on_i[valid]
            median_diff = float(np.median(diff))
            results.append({
                'order_i': i, 'order_j': j,
                'wave_lo': float(overlap_lo), 'wave_hi': float(overlap_hi),
                'n_points': int(valid.sum()),
                'median_diff': median_diff,
                'median_pct': median_diff * 100.0,
            })

    results.sort(key=lambda r: -abs(r['median_diff']))
    return results


def flagged_overlap_ranges(results, threshold_pct=2.0):
    """
    Filter check_order_overlaps() results down to pairs that disagree
    badly enough to matter, keeping just what's needed to flag a LINE
    landing inside one of them (see Spectrum_Data.flag_order_overlaps()
    and check_for_flags()) -- the wavelength range and a human-readable
    reason, not the full diagnostic record.

    threshold_pct : flag pairs whose |median_pct| meets or exceeds this.
        2.0 sits above the ~0.5-1.5% typical residual disagreement seen
        on real, otherwise-unremarkable order pairs (Keck, GRACES,
        MAROON-X all show this baseline level), while still catching
        genuine problems -- confirmed on real data: MAROON-X's known-bad
        edge regions (order 0's noise, orders 58-61's dichroic-arm
        overlap) measured 2.2-6.2%, and a real, previously-unnoticed
        GRACES order-edge continuum bug measured 5.3%.
    """
    return [
        {
            'wave_lo': r['wave_lo'], 'wave_hi': r['wave_hi'],
            'order_i': r['order_i'], 'order_j': r['order_j'],
            'median_pct': r['median_pct'],
        }
        for r in results if abs(r['median_pct']) >= threshold_pct
    ]


def print_overlap_report(results, flag_threshold_pct=1.0, top=None):
    """
    Print check_order_overlaps() results, worst first.

    flag_threshold_pct : pairs whose |median_pct| meets or exceeds this
        get a trailing "<-- FLAG" marker.
    top : only print this many (worst) pairs; None prints all.
    """
    shown = results[:top] if top else results
    for r in shown:
        flag = ' <-- FLAG' if abs(r['median_pct']) >= flag_threshold_pct else ''
        print(f"order {r['order_i']:2d} vs order {r['order_j']:2d}  "
              f"({r['wave_lo']:.1f}-{r['wave_hi']:.1f} A, {r['n_points']} pts)  "
              f"median diff = {r['median_pct']:+.2f}%{flag}")
