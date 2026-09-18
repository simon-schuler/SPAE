"""Cross-check a line's local wing against an independent, high-S/N
reference spectrum (e.g. the Kurucz solar flux atlas) to flag wing points
that sit on a real -- if too shallow for this spectrum's own noise to
reject -- absorption feature, rather than genuine continuum.

This is deliberately independent of Spectrum_Data's own flux: a feature
too shallow to clear estimate_local_continuum()'s median/MAD clip (a
gradual, few-percent-deep blend, not a sharp outlier) is exactly the case
that clip can't catch on its own -- see line_profile.py's
estimate_local_continuum() docstring. A reference atlas settles it
directly instead of trying to infer it from the noisier working spectrum.

Only useful when a suitable reference actually exists for the target (so
far: the Sun). Everything here is opt-in -- Spectrum_Data.measure_ew()
falls back to its existing behavior when no reference atlas is loaded.
"""

import os

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .line_profile import gfit_simple


def load_reference_atlas(path):
    """Load a Kurucz-style fluxspliced.2005 reference atlas: whitespace-
    delimited text, one header line, columns (wavelength_nm, flux,
    wavenumber_cm-1). Returns (wave_ang, flux), both sorted ascending in
    wavelength.

    Caches the parsed result as `<path>.npy` next to the source file --
    the source text file is tens of MB and re-parsing it with np.loadtxt
    on every run is wasteful (same reasoning as this package's existing
    kurucz_atmosphere_data_4D.npy cache). Re-parses if the source file is
    newer than the cache.
    """
    cache_path = path + '.npy'
    if os.path.exists(cache_path) and os.path.getmtime(cache_path) >= os.path.getmtime(path):
        wave_ang, flux = np.load(cache_path)
        return wave_ang, flux

    wave_nm, flux = np.loadtxt(path, skiprows=1, usecols=(0, 1), unpack=True)
    wave_ang = wave_nm * 10.0  # nm -> Angstrom; empirically matches this
    # package's air-wavelength linelists directly (confirmed against Fe I
    # 5054.643: atlas minimum at 5054.635 A, 8 mA agreement, no vacuum/air
    # correction needed)
    order = np.argsort(wave_ang)
    wave_ang, flux = wave_ang[order], flux[order]

    np.save(cache_path, np.vstack([wave_ang, flux]))
    return wave_ang, flux


def estimate_resolving_power(spectrum, lines=None, window=0.75,
                              min_depth=0.05, max_depth=0.4, percentile=25):
    """Empirically estimate the spectrograph's resolving power R = lambda
    / FWHM by fitting a Gaussian to a sample of weak-to-moderate lines and
    taking a robust LOW percentile of the resulting FWHM -- not the
    median.

    Uses the spectrum's own already-loaded linelist by default (`lines`
    is spectrum.lines, i.e. the same lines measure_all_ew() will measure)
    rather than a fixed reference set: EARLIER VERSION used a hardcoded
    set of famously strong lines (Mg b1/b2, Na D1/D2, Ca II H/K), which
    was WRONG -- confirmed on real data, those lines are pressure/damping-
    broadened far beyond the instrumental profile (fitted FWHM 0.48-1.38 A,
    implying R as low as ~3800 on a spectrograph actually around R~50-70k),
    so their width reflects line physics, not resolution.

    min_depth/max_depth restrict the sample to weak-to-moderate,
    unsaturated lines (max_depth=0.4 specifically excludes strong,
    damping-broadened lines like the ones above); a LOW percentile
    (default 25th) of the surviving FWHM is used rather than the median
    because blending or any remaining intrinsic broadening can only
    WIDEN a line, never narrow it below the true instrumental floor -- so
    the narrowest lines in the sample are the best available proxy for
    that floor, and a percentile is more robust to any single
    anomalously-narrow bad fit than a strict minimum would be.

    Returns
    -------
    R : resolving power from the chosen percentile of fitted FWHM across
        the lines that fit successfully, or None if fewer than 5 lines
        were usable (too little to trust a percentile from)
    """
    if lines is None:
        if spectrum.lines is None:
            raise ValueError(
                "estimate_resolving_power() needs a line sample: either call "
                "spectrum.load_lines() before load_reference_atlas(), or pass "
                "lines= explicitly.")
        lines = spectrum.lines

    FWHM_over_lambda = []
    for rest_wavelength in lines:
        for order in range(len(spectrum.shifted_wavelength)):
            wave = spectrum.shifted_wavelength[order]
            if wave.min() <= rest_wavelength <= wave.max():
                flux = spectrum.normalized_flux[order]
                mask = (wave >= rest_wavelength - window) & (wave <= rest_wavelength + window)
                if mask.sum() < 5:
                    break
                x = wave[mask]
                y = 1.0 - flux[mask]
                depth = y.max()
                if not (min_depth <= depth <= max_depth):
                    break  # too shallow to trust, or too strong to trust as unsaturated
                bf, err, p0 = gfit_simple(x, y, rest_wavelength, window / 4.0, 0.0)
                if bf[0] == 0 and bf[1] == 0:
                    break  # gfit_simple's failure sentinel
                fitted_center, fitted_sigma = bf[1], bf[2]
                if abs(fitted_center - rest_wavelength) > window / 2.0 or fitted_sigma <= 0:
                    break  # fit wandered too far, or degenerate width, to trust
                fwhm = fitted_sigma * 2.3548
                FWHM_over_lambda.append(fwhm / rest_wavelength)
                break  # only need the first covering order for this line

    if len(FWHM_over_lambda) < 5:
        return None
    return float(1.0 / np.percentile(FWHM_over_lambda, percentile))


def smooth_to_resolution(wave, flux, R):
    """Gaussian-convolve a reference-atlas slice down to the target
    resolving power R = lambda/FWHM. R is treated as constant in velocity
    space (the standard echelle-spectrograph approximation), so the
    smoothing FWHM scales with wavelength -- evaluated once at this slice's
    own median wavelength, which is accurate enough for the ~1.5 A windows
    this is actually called with (FWHM(lambda) barely changes over that
    span).

    wave must be ~uniformly sampled (true of the Kurucz atlas over a small
    slice); dwave is estimated from the median point spacing.
    """
    if len(wave) < 3:
        return flux

    fwhm_ang = np.median(wave) / R
    sigma_ang = fwhm_ang / 2.3548
    dwave = np.median(np.diff(wave))
    sigma_pix = sigma_ang / dwave
    if sigma_pix <= 0:
        return flux
    return gaussian_filter1d(flux, sigma_pix, mode='nearest')


def reference_continuum_mask(x, y_err, ref_wave, ref_flux, R, clip_sigma=3.0):
    """For wavelengths `x` (e.g. a line's wing points), flag which ones
    sit on a real absorption feature in the reference atlas rather than
    genuine continuum.

    Slices the reference atlas to (a margin around) x's own range, smooths
    that slice to the target spectrograph's resolving power R (see
    smooth_to_resolution()), then flags x as "not continuum" wherever the
    smoothed reference flux there is more than clip_sigma*y_err below 1.0
    -- the same significance convention estimate_local_continuum()'s own
    median/MAD clip uses, just sourced from the atlas instead of guessed
    from the noisy wing itself: a real feature only needs excluding if
    it's big enough that this spectrum's own noise wouldn't already hide
    it.

    Returns
    -------
    ok : boolean mask into x, True where the point looks like genuine
        continuum (safe to use), False where the reference atlas shows a
        real feature there. All-True (nothing flagged) if x falls outside
        the reference atlas's coverage.
    """
    lo, hi = x.min() - 0.5, x.max() + 0.5
    sel = (ref_wave >= lo) & (ref_wave <= hi)
    if sel.sum() < 3:
        return np.ones(len(x), dtype=bool)

    ref_wave_local, ref_flux_local = ref_wave[sel], ref_flux[sel]
    smoothed = smooth_to_resolution(ref_wave_local, ref_flux_local, R)
    ref_on_x = np.interp(x, ref_wave_local, smoothed)

    threshold = 1.0 - clip_sigma * np.clip(y_err, 1e-6, None)
    return ref_on_x >= threshold
