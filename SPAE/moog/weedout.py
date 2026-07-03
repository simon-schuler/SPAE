"""
Line list culling by line/continuum opacity ratio.
Translated from Weedout.f.

Entry points
------------
weedout(state, xratio)            — I/O-free; returns result dict
weedout_from_files(state, xratio) — reads params/model/lines then calls weedout()
"""
import numpy as np

from .params  import params
from .inmodel import inmodel
from .inlines import inlines
from .eqlib   import eqlib
from .nearly  import nearly
from .opacit  import opacit


def weedout(state, xratio: float) -> dict:
    """
    Classify lines as 'kept' or 'discarded' based on line/continuum opacity ratio.

    A line is kept when  strength[j] / kaplam[jtau5] >= xratio, where
      strength[j]  = kapnu0[j, jtau5]  (line-centre opacity at tau_5000 ≈ 0.5)
      kaplam[jtau5] = continuum opacity at the linelist midpoint wavelength

    nearly(state, 1) is called internally (matches Weedout.f structure).

    Parameters
    ----------
    state  : State, fully populated — inlines and eqlib must already have run.
    xratio : minimum line/continuum opacity ratio; lines below this are discarded.

    Returns
    -------
    dict with keys:
      'kept'        : list of kept line info dicts
      'discarded'   : list of discarded line info dicts
      'xratio'      : threshold used
      'kaplam_tau5' : continuum opacity at jtau5 (the reference value)
      'wave_mid'    : midpoint wavelength [Å] used for continuum computation

    Each line info dict contains:
      wave        : wavelength [Å]
      species     : atom1 species code
      ep          : lower excitation potential [eV]
      loggf       : log10(gf)
      dampnum     : raw damping value as stored in state (van der Waals C6 or
                    Barklem gamma, depending on the original linelist)
      d0          : dissociation energy [eV] (molecules) or 0 (atoms)
      ew_obs      : observed EW [mÅ] (0 for linelists without EW column)
      logstrength : log10(kapnu0 at jtau5), or None if strength ≤ 0
      ratio       : strength / kaplam_tau5
    """
    state.mode    = 1
    state.waveold = 0.0   # force continuum recompute at first line

    # Compute line-centre opacities for all lines; sets state.strength[j]
    nearly(state, 1)

    # Recompute continuum at linelist midpoint wavelength (Weedout.f line 81)
    wave_mid = (state.wave1[0] + state.wave1[state.nlines - 1]) / 2.0
    opacit(state, 2, wave_mid)
    kaplam_tau5 = float(state.kaplam[state.jtau5])

    kept      = []
    discarded = []

    for j in range(state.nlines):
        ratio = state.strength[j] / kaplam_tau5
        s     = float(state.strength[j])

        line_info = {
            'wave':        float(state.wave1[j]),
            'species':     float(state.atom1[j]),
            'ep':          float(state.e[j, 0]),
            'loggf':       float(np.log10(state.gf[j])),
            'dampnum':     float(state.dampnum[j]),
            'd0':          float(state.d0[j]),
            'ew_obs':      float(state.width[j]) * 1000.0,
            'logstrength': float(np.log10(s)) if s > 0.0 else None,
            'ratio':       float(ratio),
        }

        if ratio >= xratio:
            kept.append(line_info)
        else:
            discarded.append(line_info)

    return {
        'kept':        kept,
        'discarded':   discarded,
        'xratio':      xratio,
        'kaplam_tau5': kaplam_tau5,
        'wave_mid':    float(wave_mid),
    }


def weedout_from_files(state, xratio: float) -> dict:
    """
    Full pipeline: read params → model → eqlib → lines → weedout.

    Note: eqlib is called before inlines here, matching Weedout.f.
    state.fparam must point to a valid batch.par (or equivalent) before calling.
    xratio is the minimum line/continuum opacity ratio to keep a line.
    """
    params(state, state.fparam)
    inmodel(state)
    eqlib(state)        # called before inlines — matches Weedout.f ordering
    inlines(state, 1)
    return weedout(state, xratio)
