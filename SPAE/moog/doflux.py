"""
Continuum flux curve computation.
Translated from Doflux.f.

Entry points
------------
doflux(state)            — I/O-free; returns result dict
doflux_from_files(state) — reads params/model then calls doflux()
"""
import numpy as np

from .params     import params
from .inmodel    import inmodel
from .opacit     import opacit
from .cdcalc     import cdcalc
from .math_utils import rinteg


def doflux(state) -> dict:
    """
    Compute the emergent continuum flux at each wavelength from start to sstop.

    Wavelength range and step are taken from state.start, state.sstop,
    state.step (set by the 'fluxlimits' keyword in batch.par).

    No linelist is required; eqlib is not called.

    Returns
    -------
    dict with key 'flux': list of per-wavelength result dicts, each containing:
      wave    : wavelength [Å]
      flux    : emergent continuum flux [erg/cm²/s/Å]
      waveinv : inverse wavelength [1/μm] = 1e4/wave
      fluxlog : log10(flux), or -1.0 if flux == 0
    """
    ntau = state.ntau
    wave = state.start
    results = []

    while wave <= state.sstop:
        state.wave = wave         # cdcalc uses state.wave for the Planck function
        opacit(state, 2, wave)    # compute kaplam and taulam at this wavelength
        cdcalc(state, 1)          # compute continuum contribution function → state.cd

        first = 0.4343 * state.cd[0]
        flux, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

        if flux <= 0.1:
            flux = 0.0

        results.append({
            'wave':    float(wave),
            'flux':    float(flux),
            'waveinv': float(1.0e4 / wave),
            'fluxlog': float(np.log10(flux)) if flux > 0.0 else -1.0,
        })

        wave += state.step

    return {'flux': results}


def doflux_from_files(state) -> dict:
    """
    Full pipeline: read params → model → doflux.

    state.fparam must point to a valid batch.par (or equivalent) before calling.
    The batch.par must contain a 'fluxlimits' entry to set start/sstop/step.
    """
    params(state, state.fparam)
    inmodel(state)
    return doflux(state)
