"""
Compute a single line profile and equivalent width.
Translated from Oneline.f.

oneline() fills state.w[state.ncurve] with the EW [Å] but does NOT
increment state.ncurve — that is the caller's responsibility.
"""
import numpy as np

from .opacit  import opacit
from .taukap  import taukap
from .cdcalc  import cdcalc
from .math_utils import rinteg


def oneline(state, imode: int) -> None:
    """
    Compute the EW of line state.lim1 using state.gf1[state.ncurve].

    imode == 0 : initialise gf1[ncurve] from gf[lim1] then compute
    imode == 1 : normal iteration pass
    imode == 2 : silent pass (COG or final adjustment)

    Result stored in state.w[state.ncurve].  Caller must increment ncurve.
    """
    ntau = state.ntau
    lim1 = state.lim1

    if imode == 0:
        state.gf1[state.ncurve] = state.gf[lim1]

    state.dellam[0] = 0.0

    # Wavelength step in Å (rounded to 4 decimal places, matching Fortran idnint)
    if state.wavestep == 0.0:
        st1 = round(state.wave1[lim1] * state.dopp[lim1, state.jtau5]
                    / (2.997929e10 * 5.0), 4)
    else:
        st1 = state.wavestep

    # Continuum: recompute only when wavelength shifts by more than 30 Å
    if abs(state.wave1[lim1] - state.waveold) > 30.0:
        state.waveold = state.wave1[lim1]
        opacit(state, 2, state.wave1[lim1])
        cdcalc(state, 1)
        first = 0.4343 * state.cd[0]
        state.flux, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

    # Adapt step size so that the line depth at 5×st1 is 60–80 % of line centre
    if state.wavestep == 0.0:
        state.wave = state.wave1[lim1]
        taukap(state)
        cdcalc(state, 2)
        first = 0.4343 * state.cd[0]
        d0, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

        for _ in range(30):
            state.wave = state.wave1[lim1] + 5.0 * st1
            taukap(state)
            cdcalc(state, 2)
            first = 0.4343 * state.cd[0]
            d1, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)
            d2d1 = d1 / d0 if d0 != 0.0 else 1.0
            if d2d1 <= 0.2:
                st1 /= 1.5
            elif d2d1 <= 0.6:
                st1 /= 1.2
            elif d2d1 < 0.8:
                break          # target 60–80 % depth ratio
            else:
                st1 *= 1.6     # d2d1 >= 0.80 (Fortran: >= 0.80 before >= 0.90 branch)
            st1 = round(st1, 4)

    state.st1 = st1

    # Line profile: compute depth at each wavelength step until wings recovered
    maxsteps = 100
    d_prof = np.empty(maxsteps)
    ndepths = maxsteps

    for n in range(maxsteps):
        state.dellam[n] = n * st1
        state.wave = state.wave1[lim1] + state.dellam[n]
        taukap(state)
        cdcalc(state, 2)
        first = 0.4343 * state.cd[0]
        d_prof[n], _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)
        if n > 0 and d_prof[0] > 0.0 and d_prof[n] / d_prof[0] < 0.005:
            ndepths = n + 1
            break

    state.ndepths = ndepths

    # Symmetrize profile: [wing_far .. wing_1, centre, wing_1 .. wing_far]
    d_wing = d_prof[:ndepths]
    d_sym  = np.concatenate([d_wing[ndepths-1:0:-1], d_wing])
    ndep   = 2 * ndepths - 1
    dellam_sym = np.arange(-(ndepths - 1), ndepths, dtype=float) * st1

    # Boundary condition: linear extrapolation from the outermost wing point
    first_ew = 2.0 * dellam_sym[-1] * d_sym[-1]
    state.w[state.ncurve], _ = rinteg(dellam_sym, d_sym, ndep, first_ew)
