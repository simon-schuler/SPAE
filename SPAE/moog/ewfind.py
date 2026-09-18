"""
Equivalent width prediction from known abundances.
Translated from Ewfind.f.

Entry points
------------
ewfind(state)             — I/O-free; returns result dict
ewfind_from_files(state)  — reads params/model/lines then calls ewfind()
"""
import numpy as np

from .params   import params
from .inmodel  import inmodel
from .inlines  import inlines
from .eqlib    import eqlib
from .nearly   import nearly
from .oneline  import oneline
from .cdcalc   import cdcalc
from .math_utils import rinteg


def ewfind(state) -> dict:
    """
    Predict equivalent widths for all lines in the linelist using the
    model abundances already stored in state.xabund.

    Requires state to be fully populated:
      inlines(state, 1), eqlib(state), nearly(state, 1) must already have run.

    Returns
    -------
    dict with key 'lines': list of per-line result dicts, each containing:
      wave        : wavelength [Å]
      species     : atom1 species code (e.g. 26.0 = Fe I)
      ep          : excitation potential [eV]
      loggf       : log10(gf)
      abund       : log-epsilon abundance used
      ew_pred     : predicted equivalent width [mÅ]
      taunu0      : cumulative line-centre optical depth at each depth layer
      cd          : contribution function at each depth layer
      depth_taulam1  : geometric depth [km] where tau_continuum = 1 (None if not reached)
      depth_taunu1   : geometric depth [km] where tau_line_centre = 1 (None if not reached)
      depth_tautot1  : geometric depth [km] where tau_cont + tau_line = 1 (None if not reached)
      depth_mean     : mean line-centre formation depth [km] weighted by |C_d|
    """
    state.mode    = 1
    state.waveold = 0.0   # force continuum recompute at first line

    ntau     = state.ntau
    lines_out = []

    for lim1 in range(state.nlines):
        state.lim1   = lim1
        state.lim2   = lim1
        state.ncurve = 0
        state.gf1[0] = state.gf[lim1]

        # Single oneline pass at the input abundance (no iteration)
        oneline(state, 1)
        ew_pred = state.w[0]   # EW [Å]

        iatom = int(state.atom1[lim1] + 0.0001)
        abund = np.log10(state.xabund[iatom - 1]) + 12.0

        # ------------------------------------------------------------------
        # Formation-depth diagnostics (Ewfind.f lines 103–182)
        # ------------------------------------------------------------------
        # Line-centre kapnu from the pre-computed kapnu0 array (set by nearly)
        kapnu_c            = state.kapnu0[lim1, :ntau].copy()
        state.kapnu[:ntau] = kapnu_c

        integrand = state.tauref[:ntau] * kapnu_c / (0.4343 * state.kapref[:ntau])
        first     = state.tauref[0] * kapnu_c[0] / state.kapref[0]

        # rinteg fills fint with incremental contributions (fint[0] = 0 from start=0)
        _, fint = rinteg(state.xref[:ntau], integrand, ntau, 0.0)
        fint[0] = first   # override boundary (matches Fortran taunu0(1)=first)

        # Fortran sets taunu = taunu0 (incremental) then calls cdcalc(2) before
        # computing the cumulative sum — faithfully translated here.
        state.taunu[:ntau] = fint
        cdcalc(state, 2)
        cd = state.cd[:ntau].copy()

        # Cumulative line-centre optical depth (for depth search)
        taunu0 = np.cumsum(fint)

        taulam = state.taulam[:ntau]
        xdepth = state.xdepth[:ntau]
        xref   = state.xref[:ntau]

        def _interp_depth(tau_arr, tau_thresh=1.0):
            """Linear interpolation to find depth where tau_arr crosses tau_thresh."""
            for i in range(1, ntau):
                if tau_arr[i] >= tau_thresh:
                    t0, t1 = tau_arr[i - 1], tau_arr[i]
                    x0, x1 = xdepth[i - 1], xdepth[i]
                    return x0 + (tau_thresh - t0) * (x1 - x0) / (t1 - t0)
            return None

        depth_taulam1 = _interp_depth(taulam)
        depth_taunu1  = _interp_depth(taunu0) if taunu0[-1] >= 1.0 else None
        depth_tautot1 = _interp_depth(taulam + taunu0)

        # Mean formation depth weighted by |C_d|
        abs_cd = np.abs(cd)
        _, xref_cd_fint = rinteg(xref, xref * abs_cd, ntau, 0.0)
        _, cd_fint      = rinteg(xref, abs_cd,        ntau, 0.0)
        xrefcdinteg = float(np.sum(xref_cd_fint))
        cdinteg     = float(np.sum(cd_fint))
        if cdinteg != 0.0:
            xrefmean   = xrefcdinteg / cdinteg
            depth_mean = None
            for i in range(1, ntau):
                if xrefmean <= xref[i]:
                    x0, x1 = xdepth[i - 1], xdepth[i]
                    r0, r1 = xref[i - 1], xref[i]
                    depth_mean = x0 + (xrefmean - r0) * (x1 - x0) / (r1 - r0)
                    break
        else:
            depth_mean = None

        lines_out.append({
            'wave':          float(state.wave1[lim1]),
            'species':       float(state.atom1[lim1]),
            'ep':            float(state.e[lim1, 0]),
            'loggf':         float(np.log10(state.gf[lim1])),
            'abund':         float(abund),
            'ew_pred':       float(ew_pred) * 1000.0,   # mÅ
            'taunu0':        taunu0.tolist(),
            'cd':            cd.tolist(),
            'depth_taulam1': depth_taulam1,
            'depth_taunu1':  depth_taunu1,
            'depth_tautot1': depth_tautot1,
            'depth_mean':    depth_mean,
        })

    return {'lines': lines_out}


def ewfind_from_files(state) -> dict:
    """
    Full pipeline: read params → model → lines → eqlib → nearly → ewfind.

    state.fparam must point to a valid batch.par (or equivalent) before calling.
    """
    params(state, state.fparam)
    inmodel(state)
    inlines(state, 1)
    eqlib(state)
    nearly(state, 1)
    return ewfind(state)
