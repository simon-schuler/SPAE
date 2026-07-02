"""
Continuous opacity from scattering: electron, H I, H2, He I Rayleigh.
Translated from Opacscat.f.
"""
import numpy as np

_C_ANG = 2.997925e18    # speed of light [Å/s]


def opac_escat(state) -> None:
    """Electron scattering: σ_el = 6.653e-25 × ne [cm²/e⁻ × cm⁻³ = cm⁻¹]."""
    ntau = state.ntau
    state.sigel[:ntau] = 6.653e-25 * state.ne[:ntau]


def opac_hscat(state) -> None:
    """H I Rayleigh scattering opacity into state.sigH."""
    ntau = state.ntau
    freq = state.freq
    wavetemp = _C_ANG / min(freq, 2.463e15)   # λ in Å, capped at Lyα
    ww = wavetemp ** 2
    sig = (5.799e-13 + 1.422e-6 / ww + 2.784 / (ww * ww)) / (ww * ww)
    state.sigH[:ntau] = sig * 2.0 * state.numdens[0, 0, :ntau] / state.u[0, 0, :ntau]


def opac_h2scat(state) -> None:
    """H2 Rayleigh scattering opacity into state.sigH2."""
    ntau = state.ntau
    freq = state.freq
    wavetemp = _C_ANG / min(freq, 2.463e15)
    ww = wavetemp ** 2
    sig = (8.14e-13 + 1.28e-6 / ww + 1.61 / (ww * ww)) / (ww * ww)
    state.sigH2[:ntau] = sig * state.numdens[7, 0, :ntau]


def opac_hescat(state) -> None:
    """He I Rayleigh scattering opacity into state.sigHe."""
    ntau = state.ntau
    freq = state.freq
    wavetemp = _C_ANG / min(freq, 5.15e15)    # λ in Å, capped at He ionization
    ww = wavetemp ** 2
    sig = (5.484e-14 / (ww * ww)
           * (1.0 + (2.44e5 + 5.94e10 / (ww - 2.90e5)) / ww) ** 2)
    state.sigHe[:ntau] = sig * state.numdens[1, 0, :ntau] / state.u[1, 0, :ntau]
