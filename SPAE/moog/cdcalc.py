"""
Continuum and line+continuum contribution functions.
Translated from Cdcalc.f and Jexpint.f.
"""
import numpy as np
from scipy.special import expn as _scipy_expn


def expint(x: float, n: int) -> float:
    """n-th order exponential integral E_n(x) via scipy.special.expn."""
    if abs(x) >= 100.0:
        return 0.0
    if x == 0.0:
        return 0.0 if n <= 1 else 1.0 / float(n - 1)
    return float(_scipy_expn(n, float(x)))


def cdcalc(state, number: int) -> None:
    """
    Compute contribution functions into state.cd[:ntau].

    number == 1 : continuum contribution function
    number != 1 : line + continuum contribution function

    Uses state.fluxintopt (0=flux integral, 1=disk-center intensity).
    """
    ntau = state.ntau
    wave = state.wave   # current wavelength [Å]

    # Planck function Bλ [erg/cm²/s/Å/sr]
    # cdcalc.f: (1.19089e25/wave²)*1e10 / (wave³*(exp(1.43879e8/(wave*T))-1))
    scont = ((1.19089e25 / wave**2) * 1.0e10) / (
        wave**3 * (np.exp(1.43879e8 / (wave * state.t[:ntau])) - 1.0))
    state.scont[:ntau] = scont

    taulam = state.taulam[:ntau]
    kap    = state.kaplam[:ntau]
    kref   = state.kapref[:ntau]
    tref   = state.tauref[:ntau]

    if number == 1:
        if state.fluxintopt == 1:
            cd = kap * tref * scont * np.exp(-taulam) / (0.4343 * kref)
        else:
            e2 = np.array([expint(taulam[i], 2) for i in range(ntau)])
            cd = 2.0 * kap * tref * scont * e2 / (0.4343 * kref)
        state.cd[:ntau] = cd

    else:
        taunu  = state.taunu[:ntau]
        kapnu  = state.kapnu[:ntau]
        sline  = scont.copy()   # LTE: sline = scont
        state.sline[:ntau] = sline
        flux   = state.flux

        if state.fluxintopt == 1:
            tau_tot = taulam + taunu
            exptau  = np.where(tau_tot <= 50.0, np.exp(-tau_tot), 0.0)
            cd = (tref * kap / (0.4343 * flux * kref) *
                  (scont * np.exp(-taulam) -
                   (1.0 + kapnu / kap) * sline * exptau))
        else:
            e2_lam = np.array([expint(taulam[i],          2) for i in range(ntau)])
            e2_tot = np.array([expint(taulam[i]+taunu[i], 2) for i in range(ntau)])
            cd = (2.0 * tref * kap / (0.4343 * flux * kref) *
                  (scont * e2_lam -
                   (1.0 + kapnu / kap) * sline * e2_tot))
        state.cd[:ntau] = cd
