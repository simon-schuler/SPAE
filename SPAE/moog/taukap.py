"""
Line absorption coefficient and optical depth at wavelength state.wave.
Translated from Taukap.f.
"""
import numpy as np

from .math_utils import rinteg, voigt, rinteg_batch_cumsum


def taukap(state) -> None:
    """
    Compute total line opacity kapnu[:ntau] and line optical depth taunu[:ntau].

    Uses lines in range [lim1, lim2] (0-indexed, inclusive) plus any strong
    lines.  state.wave must be set to the current wavelength [Å] before calling.
    """
    ntau = state.ntau
    wave = state.wave          # current wavelength [Å]
    c    = 2.997929e10         # speed of light [cm/s]

    kapnu = np.zeros(ntau)

    # Regular lines — vectorize over all lines at once (nlines × ntau)
    if state.lim1 <= state.lim2:
        js   = np.arange(state.lim1, state.lim2 + 1)          # (nl,)
        w1   = state.wave1[js]                                  # (nl,)
        v2d  = (c * np.abs(wave - w1[:, None])
                / (w1[:, None] * state.dopp[js, :ntau]))        # (nl, ntau)
        k0   = state.kapnu0[js, :ntau]                          # (nl, ntau)
        a2d  = state.a[js, :ntau]                               # (nl, ntau)
        kapnu = (k0 * voigt(a2d, v2d)).sum(axis=0)             # (ntau,)

    # Strong lines (typically few, keep loop)
    if state.dostrong > 0:
        for j in range(state.nlines, state.nlines + state.nstrong):
            v_arr = c * np.abs(wave - state.wave1[j]) / (state.wave1[j] * state.dopp[j, :ntau])
            kapnu += state.kapnu0[j, :ntau] * voigt(state.a[j, :ntau], v_arr)

    state.kapnu[:ntau] = kapnu

    # Integrate taunu — Fortran uses start=0 then overwrites taunu[0]; using
    # start=first directly produces the same cumulative sum.
    first = state.tauref[0] * kapnu[0] / state.kapref[0]
    integrand = state.tauref[:ntau] * kapnu / (0.4343 * state.kapref[:ntau])
    _, fint = rinteg(state.xref[:ntau], integrand, ntau, first)
    state.taunu[:ntau] = np.cumsum(fint)


def taukap_batch(state, waves):
    """
    Compute kapnu and taunu for multiple wavelengths in one vectorized pass.

    Parameters
    ----------
    state : State
    waves : ndarray, shape (nwave,) — wavelengths [Å]

    Returns
    -------
    kapnu_batch : (nwave, ntau)
    taunu_batch : (nwave, ntau)
    """
    ntau  = state.ntau
    nwave = len(waves)
    c     = 2.997929e10

    kapnu_batch = np.zeros((nwave, ntau))

    if state.lim1 <= state.lim2:
        js  = np.arange(state.lim1, state.lim2 + 1)   # (nl,)
        w1  = state.wave1[js]                          # (nl,)
        k0  = state.kapnu0[js, :ntau]                  # (nl, ntau)
        a3d = state.a[js, :ntau]                       # (nl, ntau)

        # v3d: (nwave, nl, ntau)
        v3d = (c * np.abs(waves[:, None, None] - w1[None, :, None])
               / (w1[None, :, None] * state.dopp[None, js, :ntau]))

        kapnu_batch = (k0[None, :, :] * voigt(a3d[None, :, :], v3d)).sum(axis=1)

    if state.dostrong > 0:
        for j in range(state.nlines, state.nlines + state.nstrong):
            v2d = (c * np.abs(waves[:, None] - state.wave1[j])
                   / (state.wave1[j] * state.dopp[None, j, :ntau]))
            kapnu_batch += (state.kapnu0[None, j, :ntau]
                            * voigt(state.a[None, j, :ntau], v2d))

    tauref = state.tauref[:ntau]
    kapref = state.kapref[:ntau]
    xref   = state.xref[:ntau]

    integrand_b = tauref[None, :] * kapnu_batch / (0.4343 * kapref[None, :])
    first_b     = tauref[0] * kapnu_batch[:, 0] / kapref[0]

    taunu_batch = rinteg_batch_cumsum(xref, integrand_b, first_b)

    return kapnu_batch, taunu_batch
