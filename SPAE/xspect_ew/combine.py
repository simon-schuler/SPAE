"""Combining multiple Spectrum_Data objects, and helpers used by wavelength-
shift estimation/cleaning (Spectrum_Data.estimate_shift()/clean_shift())."""

import numpy as np
from scipy.interpolate import interp1d


def make_line(x,m,b):
    return m*x+b


def parabolic_refine(shifts, chi, k_min):
    """
    Sub-grid-resolution refinement of a coarse grid-search minimum via a
    3-point parabolic fit around shifts[k_min]/chi[k_min].

    Assumes a uniform grid (as produced by np.linspace, which is what
    estimate_shift() uses). Falls back to the unrefined grid value at the
    edges of the search range, where a symmetric 3-point fit isn't possible.

    Returns
    -------
    refined_shift : float
    """
    if k_min == 0 or k_min == len(shifts) - 1:
        return shifts[k_min]

    y_lo, y_mid, y_hi = chi[k_min - 1], chi[k_min], chi[k_min + 1]
    denom = (y_lo - 2.0 * y_mid + y_hi)
    if denom == 0.0:
        return shifts[k_min]

    d = shifts[k_min] - shifts[k_min - 1]  # grid spacing
    delta = 0.5 * (y_lo - y_hi) / denom
    return shifts[k_min] + delta * d


def measure_order_alignment(wave_A, flux_A, wave_B, flux_B, search_radius=0.1, n_shifts=201):
    """
    Fine-grained DIRECT cross-correlation between two spectra's own
    overlapping data, to measure any REMAINING residual wavelength
    offset between them -- for Spectrum_Data.combine_spectra(), used
    AFTER both spectra have already been independently rest-framed via
    apply_rv_shift(). Distinct from (and much more precise than) either
    spectrum's own independent absolute RV measurement: two independent
    RV measurements of the same real star routinely disagree by several
    tenths of a km/s from real per-line measurement noise alone
    (confirmed on real data: -3.97 vs -4.45 km/s from the same 4
    reference lines in both spectra) -- a real residual misalignment of
    several mA at optical wavelengths, easily enough to measurably smear
    a naively-co-added spectrum. This function instead compares the two
    spectra directly, over many points at once, which averages out each
    spectrum's own per-line noise rather than inheriting it.

    Pass NORMALIZED flux (continuum-divided) for the cleanest cross-
    correlation signal, not raw flux -- confirmed on real data: raw flux
    carries each exposure's own absolute throughput/continuum level,
    which normalized_flux removes and raw flux does not.

    Parameters
    ----------
    wave_A, flux_A : spectrum A's (already rest-framed) wavelength/
        normalized flux for one order.
    wave_B, flux_B : spectrum B's, for the matched order.
    search_radius : Angstrom half-width searched. 0.1 A comfortably
        covers the residual sizes confirmed on real data (5-19 mA
        typical, one order up to ~29 mA) after both spectra are already
        independently rest-framed -- this is hunting for a small leftover
        mismatch, not a real uncorrected RV.
    n_shifts : grid resolution before parabolic_refine()'s sub-grid
        interpolation.

    Returns
    -------
    residual : float, Angstrom -- the amount to ADD to a wavelength
        query into B's flux array so it aligns with A (i.e., evaluate
        B's flux at `wave_A_point + residual` to get the value that
        truly corresponds to `wave_A_point`). None if the overlap is too
        small/noisy to measure (fewer than 20 usable points, or every
        trial shift left fewer than 20 valid overlapping points).
    """
    lo = max(wave_A.min(), wave_B.min()) + search_radius
    hi = min(wave_A.max(), wave_B.max()) - search_radius
    if hi <= lo:
        return None
    ref_mask = (wave_A >= lo) & (wave_A <= hi)
    if ref_mask.sum() < 20:
        return None
    ref_wave, ref_flux = wave_A[ref_mask], flux_A[ref_mask]

    shifts = np.linspace(-search_radius, search_radius, n_shifts)
    chi = np.empty(n_shifts)
    interp_B = interp1d(wave_B, flux_B, bounds_error=False, fill_value=np.nan)
    for k, s in enumerate(shifts):
        test = interp_B(ref_wave + s)
        valid = np.isfinite(test)
        chi[k] = np.mean((ref_flux[valid] - test[valid]) ** 2) if valid.sum() >= 20 else np.inf

    k_min = np.argmin(chi)
    if not np.isfinite(chi[k_min]):
        return None
    return parabolic_refine(shifts, chi, k_min)


def combine_files(empty_obj,objects = []):
    final_wavelength = []
    final_flux = []
    final_norm_flux = []
    final_shifted_wavelength = []
    final_estimated_shift = []
    final_continuum = []
    final_obs_err = []
    final_pred_all = []
    final_pred_var_all = []
    final_gain = []

    for j in objects:

        for i in range(len(j.flux)):
            final_norm_flux.append(j.normalized_flux[i])
            final_shifted_wavelength.append(j.shifted_wavelength[i])
            final_wavelength.append(j.wavelength[i])
            final_flux.append(j.flux[i])
            final_estimated_shift.append(j.estimated_shift[i])
            final_continuum.append(j.continuum[i])
            final_obs_err.append(j.obs_err[i])
            final_pred_all.append(j.pred_all[i])
            final_pred_var_all.append(j.pred_var_all[i])
            final_gain.append(j.gain[i])

    empty_obj.wavelength = np.array(final_wavelength)
    empty_obj.flux = np.array(final_flux)
    empty_obj.shifted_wavelength = np.array(final_shifted_wavelength)
    empty_obj.normalized_flux = np.array(final_norm_flux)
    empty_obj.estimated_shift = np.array(final_estimated_shift)
    empty_obj.continuum = np.array(final_continuum)
    empty_obj.obs_err = np.array(final_obs_err)
    empty_obj.pred_all = np.array(final_pred_all)
    empty_obj.pred_var_all = np.array(final_pred_var_all)
    empty_obj.gain = np.array(final_gain)
    del final_wavelength
    del final_flux
    del final_norm_flux
    del final_shifted_wavelength
    del final_estimated_shift
    del final_continuum
    del final_obs_err
    del final_pred_all
    del final_pred_var_all
    del final_gain

    return empty_obj


def reduce_cc(x,y,lines,lines_removed,limit=0.12):
    #check correlation before going further
    cc = np.corrcoef(x,y)
    print('starting cc', cc[0,1])

    if abs(cc[0,1]) < limit:
        print('cc good enough')
        return lines,x,y,lines_removed

    check_ccs = np.zeros(len(x))

    #remove largest cc difference
    for i in range(len(x)):
        new_x = np.delete(x,i)
        new_y = np.delete(y,i)
        new_cc = np.corrcoef(new_x,new_y)
        check_ccs[i] = new_cc[0,1]

    #Calculate differences
    diffs = [abs(cc[0,1])- abs(j) for j in check_ccs]
    #which gives largest difference?
    biggest_diff = np.where(diffs == max(diffs))[0][0]
    #remove that one line
    lines_removed.append([lines[biggest_diff],x[biggest_diff],y[biggest_diff]])
    x = np.delete(x,biggest_diff)
    y = np.delete(y,biggest_diff)
    lines = np.delete(lines,biggest_diff)
    print('line removed: ', lines_removed)

    #recalculate cc
    cc = np.corrcoef(x,y)
    print('ending cc', cc[0,1])

    #Call function again to remove lines until 0.12 is passed
    lines,x,y,lines_removed = reduce_cc(x,y,lines,lines_removed)

    return lines,x,y,lines_removed
