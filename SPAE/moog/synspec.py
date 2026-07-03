"""
Synthetic spectrum computation.
Translated from Synspec.f (synthesis loop) and the mode-3 section of Linlimit.f.

Entry point: synspec(state) — returns (wave_arr, depth_arr).
Depths are line depths (0 = no absorption, 1 = fully absorbed).
Conversion to flux (1 - depth) is done by the caller (synth.py).
"""
import numpy as np

from .opacit   import opacit
from .taukap   import taukap
from .cdcalc   import cdcalc
from .math_utils import rinteg
from .inlines  import inlines
from .nearly   import nearly


# ---------------------------------------------------------------------------
# Line range finder for mode 3 (spectrum synthesis)
# ---------------------------------------------------------------------------

def _linlimit_synth(state) -> None:
    """
    Set state.lim1line, state.lim2line, state.lineflag for the current
    synthesis wavelength (state.wave).

    Matches Linlimit.f mode-3 logic exactly.  The sentinel lim2line = -1
    is never returned here because in Python we read all lines at once
    (inlines mode=1 reads the full file); the block-read branch is omitted.
    """
    nlines  = state.nlines
    wave    = state.wave
    delta   = state.delta
    wave1   = state.wave1

    state.lineflag = 0

    wavelo = wave - delta
    wavehi = wave + delta

    # synthesis range too far outside linelist
    if wavehi < wave1[0] - 10.0:
        raise ValueError(
            f"Synthesis begins >10 Å blueward of linelist start "
            f"({wave1[0]:.2f} Å); wavehi={wavehi:.2f} Å")
    if wavelo > wave1[nlines - 1] + 10.0:
        raise ValueError(
            f"Synthesis ends >10 Å redward of linelist end "
            f"({wave1[nlines-1]:.2f} Å); wavelo={wavelo:.2f} Å")

    # blank synthesis (no lines in window)
    if wavehi < wave1[0]:
        state.lim1line = 0
        state.lim2line = 0
        state.lineflag = -1
        return
    if wavelo > wave1[nlines - 1]:
        state.lim1line = nlines - 1
        state.lim2line = nlines - 1
        state.lineflag = -1
        return

    # advance lower limit
    lim1 = state.lim1line
    if lim1 == 0 or wavelo < wave1[0]:
        lim1 = 0
    for j in range(lim1, nlines):
        if wavelo < wave1[j]:
            lim1 = j
            break

    # advance upper limit
    lim2 = nlines - 1   # default: end of list
    for j in range(lim1, nlines):
        if wavehi < wave1[j]:
            lim2 = j - 1
            if lim1 == lim2 and wavelo > wave1[lim1]:
                state.lineflag = -1
            state.lim1line = lim1
            state.lim2line = lim2
            return

    # fell off end of list — all remaining lines are in window
    state.lim1line = lim1
    state.lim2line = nlines - 1


# ---------------------------------------------------------------------------
# Core synthesis loop
# ---------------------------------------------------------------------------

def synspec(state) -> tuple:
    """
    Compute a synthetic spectrum for the range [state.start, state.sstop].

    Requires state to be fully populated:
      inlines(state, 1), eqlib(state), nearly(state, 1) must already have run.

    Returns
    -------
    wave_arr : np.ndarray, wavelength at each synthesis point [Å]
    d_arr    : np.ndarray, line depth at each point (0 = continuum, 1 = black)
    """
    ntau  = state.ntau
    start = state.oldstart   # Fortran uses oldstart in the wave= assignment
    sstop = state.sstop
    step  = state.step

    kount = int(round((sstop - start + step / 4.0) / step)) + 1

    wave_arr = np.empty(kount, dtype=np.float64)
    d_arr    = np.empty(kount, dtype=np.float64)

    wavl = 0.0   # wavelength at which continuum was last computed
    state.lim1line = 0

    for n in range(kount):
        wave = start + n * step
        state.wave = wave
        wave_arr[n] = wave

        # Recompute continuum when wavelength shifts by ≥ 0.1%
        if abs(wave - wavl) / wave >= 0.001:
            wavl = wave
            opacit(state, 2, wave)
            cdcalc(state, 1)
            first = 0.4343 * state.cd[0]
            state.flux, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

        # Find lines in window around this wavelength (mode 3 linlimit)
        _linlimit_synth(state)
        state.lim1 = state.lim1line
        state.lim2 = state.lim2line

        # Compute line depth
        if state.lineflag < 0:
            d_arr[n] = 0.0
        else:
            taukap(state)
            cdcalc(state, 2)
            first = 0.4343 * state.cd[0]
            d_arr[n], _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

    return wave_arr, d_arr
