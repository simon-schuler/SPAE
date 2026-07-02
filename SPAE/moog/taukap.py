"""
Line absorption coefficient and optical depth at wavelength state.wave.
Translated from Taukap.f.
"""
import numpy as np

from .math_utils import rinteg, voigt


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

    # Regular lines
    for j in range(state.lim1, state.lim2 + 1):
        v_arr = c * np.abs(wave - state.wave1[j]) / (state.wave1[j] * state.dopp[j, :ntau])
        kapnu += state.kapnu0[j, :ntau] * voigt(state.a[j, :ntau], v_arr)

    # Strong lines
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
