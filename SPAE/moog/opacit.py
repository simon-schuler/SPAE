"""
Continuous opacity driver — computes kaplam and taulam at a given wavelength.
Translated from Opacit.f.
"""
import numpy as np

from .math_utils import rinteg
from .opac_hydrogen import opac_h1, opac_hminus
from .opac_helium   import opac_heminus
from .opac_scatter  import opac_escat, opac_hscat, opac_h2scat, opac_hescat
from .opac_metals   import (opac_c1, opac_mg1, opac_mg2, opac_al1,
                             opac_si1, opac_si2, opac_fe1)

_H = 6.6256e-27    # Planck constant [erg·s]
_K = 1.38065e-16   # Boltzmann constant [erg/K]
_C_ANG = 2.997925e18   # speed of light [Å/s]


def opacit(state, modeop: int, waveop: float) -> None:
    """
    Compute continuous opacity at wavelength *waveop* [Å].

    modeop == 1 : set kapref = kaplam at the reference wavelength (no taulam)
    modeop != 1 : compute kaplam and integrate taulam
    """
    ntau = state.ntau

    # Frequency and stimulated-emission factor
    state.freq   = _C_ANG / waveop
    state.freqlg = np.log(state.freq)
    hkt_arr = _H / (_K * state.t[:ntau])
    state.evhkt[:ntau] = np.exp(-state.freq * hkt_arr)

    # Reset all opacity arrays to a negligible floor
    floor = 1.0e-99
    state.aH1[:ntau]      = floor
    state.aHminus[:ntau]  = floor
    state.aHeminus[:ntau] = floor
    state.aC1[:ntau]      = floor
    state.aMg1[:ntau]     = floor
    state.aMg2[:ntau]     = floor
    state.aAl1[:ntau]     = floor
    state.aSi1[:ntau]     = floor
    state.aSi2[:ntau]     = floor
    state.aFe1[:ntau]     = floor

    # Compute individual opacity contributions
    opac_h1(state)
    opac_hminus(state)
    opac_hscat(state)
    opac_h2scat(state)
    opac_heminus(state)
    opac_hescat(state)
    opac_escat(state)
    opac_c1(state)
    opac_mg1(state)
    opac_mg2(state)
    opac_al1(state)
    opac_si1(state)
    opac_si2(state)
    opac_fe1(state)

    # Sum opacities
    state.kaplamabs[:ntau] = (state.aH1[:ntau]     + state.aHminus[:ntau] +
                               state.aHeminus[:ntau] + state.aC1[:ntau]    +
                               state.aMg1[:ntau]    + state.aMg2[:ntau]   +
                               state.aAl1[:ntau]    + state.aSi1[:ntau]   +
                               state.aSi2[:ntau]    + state.aFe1[:ntau])
    state.kaplamsca[:ntau] = (state.sigH[:ntau]  + state.sigH2[:ntau] +
                               state.sigHe[:ntau] + state.sigel[:ntau])
    state.kaplam[:ntau]    = state.kaplamabs[:ntau] + state.kaplamsca[:ntau]

    # Optional fudge factor
    if state.fudge > 0.0:
        state.kaplam[:ntau] *= (state.fudge * 10000.0) / state.t[:ntau]

    if modeop == 1:
        # Reference opacity assignment
        state.kapref[:ntau] = state.kaplam[:ntau]
        return

    # Integrate taulam
    integrand = (state.tauref[:ntau] * state.kaplam[:ntau] /
                 (0.4343 * state.kapref[:ntau]))
    first = state.tauref[0] * state.kaplam[0] / state.kapref[0]
    _, fint = rinteg(state.xref[:ntau], integrand, ntau, first)
    state.taulam[:ntau] = np.cumsum(fint)
