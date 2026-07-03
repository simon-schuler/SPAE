"""
Curve-of-growth computation for each line in the linelist.
Translated from Cog.f and Curve.f.

Entry points
------------
cog(state)             — I/O-free; returns result dict
cog_from_files(state)  — reads params/model/lines then calls cog()
"""
import numpy as np

from .params   import params
from .inmodel  import inmodel
from .inlines  import inlines
from .eqlib    import eqlib
from .nearly   import nearly
from .fakeline import _curve


def cog(state) -> dict:
    """
    Compute a curve of growth for every line in the linelist.

    Uses state.rwlow, state.rwhigh, state.rwstep to control the COG range
    and step (set by 'coglimits' in batch.par).  If rwstep == 0 (not set),
    sensible defaults of -6.5, -3.5, 0.15 are used.

    Requires state to be fully populated:
      inlines(state, 1), eqlib(state), nearly(state, 1) must already have run.

    Returns
    -------
    dict with key 'lines': list of per-line result dicts, each containing:
      wave    : wavelength [Å]
      species : atom1 species code (e.g. 26.0 = Fe I)
      ep      : excitation potential [eV]
      abund   : log-epsilon abundance from state.xabund
      loggf   : list of log10(gf) values along the COG
      logrw   : list of log10(EW/λ) values along the COG
      n       : number of COG points
    """
    state.mode    = 1
    state.waveold = 0.0   # force continuum recompute at first line

    # Apply defaults if coglimits was not set in the parameter file
    if state.rwstep == 0.0:
        state.rwlow  = -6.5
        state.rwhigh = -3.5
        state.rwstep =  0.15

    lines_out = []

    for lim1 in range(state.nlines):
        state.lim1 = lim1
        state.lim2 = lim1

        # _curve() computes the full COG for state.lim1 and stores:
        #   state.gf1[0..ncurve]  = log10(gf) at each COG point
        #   state.w[0..ncurve]    = log10(RW) at each COG point
        #   state.ncurve          = index of the last point
        _curve(state)

        n_pts = state.ncurve + 1

        iatom = int(state.atom1[lim1] + 0.0001)
        if iatom >= 100:
            iatom = 1   # molecules: use H abundance (matches Fortran Curve.f)
        abund = np.log10(state.xabund[iatom - 1]) + 12.0

        lines_out.append({
            'wave':    float(state.wave1[lim1]),
            'species': float(state.atom1[lim1]),
            'ep':      float(state.e[lim1, 0]),
            'abund':   float(abund),
            'loggf':   state.gf1[:n_pts].tolist(),
            'logrw':   state.w[:n_pts].tolist(),
            'n':       n_pts,
        })

    return {'lines': lines_out}


def cog_from_files(state) -> dict:
    """
    Full pipeline: read params → model → lines → eqlib → nearly → cog.

    state.fparam must point to a valid batch.par (or equivalent) before calling.
    """
    params(state, state.fparam)
    inmodel(state)
    inlines(state, 1)
    eqlib(state)
    nearly(state, 1)
    return cog(state)
