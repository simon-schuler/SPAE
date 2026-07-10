"""
Continuum and line+continuum contribution functions.
Translated from Cdcalc.f and Jexpint.f.
"""
import numpy as np
from scipy.special import expn as _scipy_expn


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
            e2 = _scipy_expn(2, taulam)
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
            e2_lam = _scipy_expn(2, taulam)
            e2_tot = _scipy_expn(2, taulam + taunu)
            cd = (2.0 * tref * kap / (0.4343 * flux * kref) *
                  (scont * e2_lam -
                   (1.0 + kapnu / kap) * sline * e2_tot))
        state.cd[:ntau] = cd


def cdcalc_batch(state, kapnu_batch, taunu_batch, waves):
    """
    Batched cdcalc (line+continuum) for multiple wavelengths.

    taulam and kaplam must already be set in state (from opacit at line centre).

    Parameters
    ----------
    state        : State
    kapnu_batch  : (nwave, ntau)
    taunu_batch  : (nwave, ntau)
    waves        : (nwave,) — wavelengths [Å] for Planck function

    Returns
    -------
    cd_batch : (nwave, ntau)
    """
    ntau   = state.ntau
    taulam = state.taulam[:ntau]   # (ntau,) — fixed
    kap    = state.kaplam[:ntau]   # (ntau,)
    kref   = state.kapref[:ntau]   # (ntau,)
    tref   = state.tauref[:ntau]   # (ntau,)
    flux   = state.flux

    # Planck function per wavelength: (nwave, ntau)
    scont_b = ((1.19089e25 / waves[:, None] ** 2) * 1.0e10) / (
        waves[:, None] ** 3 * (
            np.exp(1.43879e8 / (waves[:, None] * state.t[:ntau][None, :])) - 1.0))

    factor = 2.0 * tref[None, :] * kap[None, :] / (0.4343 * flux * kref[None, :])

    if state.fluxintopt == 1:
        tau_tot  = taulam[None, :] + taunu_batch      # (nwave, ntau)
        exptau   = np.where(tau_tot <= 50.0, np.exp(-tau_tot), 0.0)
        cd_batch = (factor
                    * (scont_b * np.exp(-taulam)[None, :]
                       - (1.0 + kapnu_batch / kap[None, :]) * scont_b * exptau))
    else:
        e2_lam = _scipy_expn(2, taulam)                            # (ntau,)
        e2_tot = _scipy_expn(2, taulam[None, :] + taunu_batch)    # (nwave, ntau)
        cd_batch = (factor
                    * (scont_b * e2_lam[None, :]
                       - (1.0 + kapnu_batch / kap[None, :]) * scont_b * e2_tot))

    return cd_batch
