"""
Synthetic spectrum driver.
Translated from Synth.f and Getsyns.f.

Entry points:
  synth(state)           — I/O-free; returns result dict
  synth_from_files(state) — reads params/model/lines from disk then calls synth()
"""
import numpy as np

from .params   import params
from .inmodel  import inmodel
from .inlines  import inlines
from .eqlib    import eqlib
from .nearly   import nearly
from .synspec  import synspec
from .smooth   import smooth as _smooth


def synth(state,
          smtype: str = None,
          vsini: float = None,
          limbdark: float = None,
          vmac: float = None,
          fwhmgauss: float = None,
          fwhmloren: float = None,
          addflux: float = None) -> dict:
    """
    Synthesise a spectrum and return wavelength, raw flux, and smoothed flux.

    Smoothing parameters default to whatever is already in state (set by
    params() from batch.par).  Keyword arguments override state values.

    Parameters
    ----------
    state       : State, fully populated (model + lines loaded, eqlib run)
    smtype      : smoothing type ('n','g','l','v','c','m','d','r'); None → state.smtype
    vsini       : rotation [km/s]
    limbdark    : limb darkening coefficient (with vsini)
    vmac        : macroturbulence [km/s]
    fwhmgauss   : Gaussian FWHM [Å]
    fwhmloren   : Lorentzian FWHM [Å]
    addflux     : veiling fraction

    Returns
    -------
    dict with keys:
      'wave'        : np.ndarray, wavelengths [Å]
      'flux_raw'    : np.ndarray, unsmoothed flux (= 1 - depth)
      'flux_smooth' : np.ndarray, smoothed flux (same as flux_raw if smtype='n')
    """
    # Apply kwarg overrides to state
    if smtype    is not None: state.smtype    = smtype
    if vsini     is not None: state.vsini     = vsini
    if limbdark  is not None: state.limbdark  = limbdark
    if vmac      is not None: state.vmac      = vmac
    if fwhmgauss is not None: state.fwhmgauss = fwhmgauss
    if fwhmloren is not None: state.fwhmloren = fwhmloren
    if addflux   is not None: state.addflux   = addflux

    state.mode = 3

    # If multiple abundance variations, run the first one only (the common case)
    if state.numpecatom == 0 or state.numatomsyn == 0:
        state.isynth = 1
        state.isorun = 1
        state.nlines = 0
        state.waveold = 0.0    # force continuum recompute at first step
        inlines(state, 1)
        eqlib(state)
        nearly(state, 1)
        wave_arr, d_arr = synspec(state)
    else:
        # Multiple synthesis passes (abundance grid)
        wave_arr_list = []
        d_arr_list = []
        for n in range(state.numatomsyn):
            state.isynth = n + 1
            state.isorun = n + 1
            state.start  = state.oldstart
            state.sstop  = state.oldstop
            state.mode   = 3
            state.waveold = 0.0
            inlines(state, 1)
            eqlib(state)
            nearly(state, 1)
            w, d = synspec(state)
            wave_arr_list.append(w)
            d_arr_list.append(d)
        wave_arr = wave_arr_list[0]
        d_arr    = d_arr_list[0]    # primary synthesis (first abundance variation)

    # Convert depths → flux (Smooth.f line 347: y(i) = 1 - y(i))
    flux_raw = 1.0 - d_arr

    # Apply smoothing
    sm = state.smtype
    step = state.step

    do_rot  = sm in ('v', 'c', 'r')
    do_mac  = sm in ('m', 'd', 'r')
    do_gau  = sm in ('g', 'c', 'd', 'r')
    do_lor  = sm in ('l',)

    flux_smooth = _smooth(
        wave_arr, flux_raw, step,
        vsini     = state.vsini     if do_rot else 0.0,
        limbdark  = state.limbdark  if do_rot else 0.0,
        vmac      = state.vmac      if do_mac else 0.0,
        fwhmgauss = state.fwhmgauss if do_gau else 0.0,
        fwhmloren = state.fwhmloren if do_lor else 0.0,
        addflux   = state.addflux,
    )

    return {
        'wave':        wave_arr,
        'flux_raw':    flux_raw,
        'flux_smooth': flux_smooth,
    }


def synth_from_files(state, **kwargs) -> dict:
    """
    Full pipeline: read params → model → lines → eqlib → synth.

    state.fparam must point to a valid batch.par before calling.
    Smoothing kwargs are forwarded to synth().
    """
    params(state, state.fparam)
    inmodel(state)
    return synth(state, **kwargs)
