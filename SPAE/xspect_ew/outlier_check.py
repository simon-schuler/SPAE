"""Whole-order data-quality flagging: detect an order whose raw flux
contains an extreme, non-astrophysical outlier (e.g. a cosmic ray or a
fixed detector/amplifier-boundary defect) that would corrupt any
continuum fit run on it -- rather than trying to make the continuum-
fitting routine itself robust to every possible corruption. Same flag-
and-exclude philosophy as overlap_check.py, just for a different failure
mode (one order's own raw data, not a disagreement between two orders).

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
