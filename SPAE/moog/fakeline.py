"""
Build a curve-of-growth lookup table used by lineabund for abundance fitting.
Translated from Fakeline.f and Curve.f.

Uses a real Fe I 5006.126 Å line (which has Barklem damping data) as the
representative line.  The raw COG is computed with _curve(), then interpolated
to a fine 0.005-dex log(gf) grid stored in state.gftab / state.rwtab.
"""
import numpy as np

from .partition import partfn
from .nearly    import nearly
from .oneline   import oneline
from .atomic_data import XAM, XCHI1, XCHI2, XCHI3


# ---------------------------------------------------------------------------
# Catmull-Rom cubic interpolation (matches Fortran fakeline exactly)
# ---------------------------------------------------------------------------

def _catmull_rom(w0, w1, w2, w3, pp):
    """Cubic Lagrange interpolant at fractional position pp in [w1, w2]."""
    return (w0 * (-pp) * (pp - 1.) * (pp - 2.) / 6.0 +
            w1 * (pp * pp - 1.) * (pp - 2.) / 2.0 +
            w2 * (-pp) * (pp + 1.) * (pp - 2.) / 2.0 +
            w3 *  pp  * (pp * pp - 1.) / 6.0)


# ---------------------------------------------------------------------------
# Curve-of-growth computation  (Curve.f)
# ---------------------------------------------------------------------------

def _curve(state) -> None:
    """
    Build a raw COG for state.lim1 in the range [state.rwlow, state.rwhigh].

    After return:
      state.gf1[0..state.ncurve]  = log10(gf) at each COG point
      state.w[0..state.ncurve]    = log10(RW) at each COG point
      state.ncurve                = index of last point (total = ncurve+1)
    """
    ntau = state.ntau
    lim1 = state.lim1
    dec  = 10.0 ** state.rwstep

    wstart = 10.0 ** state.rwlow  * state.wave1[lim1]
    wstop  = 10.0 ** state.rwhigh * state.wave1[lim1]

    state.ncurve        = 0
    state.gf1[state.ncurve] = state.gf[lim1]

    # Phase 1: reduce gf until computed EW drops below wstart
    while True:
        oneline(state, 2)
        if state.w[state.ncurve] <= wstart:
            break
        state.gf1[state.ncurve]      /= dec
        state.kapnu0[lim1, :ntau]    /= dec

    # Phase 2: increase gf one step at a time, collecting COG points
    while state.w[state.ncurve] < wstop:
        state.ncurve += 1
        state.gf1[state.ncurve]   = state.gf1[state.ncurve - 1] * dec
        state.kapnu0[lim1, :ntau] *= dec
        oneline(state, 2)

    # Convert to log space (matches Fortran post-loop conversion)
    n_pts = state.ncurve + 1
    for i in range(n_pts):
        state.w[i]   = np.log10(state.w[i] / state.wave1[lim1])
        state.gf1[i] = np.log10(state.gf1[i])


# ---------------------------------------------------------------------------
# fakeline  (Fakeline.f)
# ---------------------------------------------------------------------------

def fakeline(state) -> None:
    """
    Set up a synthetic Fe I line, compute its curve of growth, and build
    the fine gftab/rwtab lookup table used by lineabund.

    Populates:
      state.gftab[0..ntabtot-1]  — log10(gf) at 0.005-dex steps
      state.rwtab[0..ntabtot-1]  — corresponding log10(RW)
      state.ntabtot              — number of table entries
    """
    ntau = state.ntau

    # Fake Fe I line parameters (5006.126 Å, same line as Fortran Fakeline.f)
    state.wave1[0]  = 5006.126
    state.atom1[0]  = 26.0       # Fe I
    state.e[0, 0]   = 2.833      # lower excitation potential [eV]
    state.e[0, 1]   = 5.308      # upper excitation potential [eV]
    state.gf[0]     = 1.0e-3
    iatom           = 26         # Fe, Z=26

    state.charge[0] = 1.0
    state.amass[0]  = XAM[iatom - 1]
    state.chi[0, 0] = XCHI1[iatom - 1]
    state.chi[0, 1] = XCHI2[iatom - 1]
    state.chi[0, 2] = XCHI3[iatom - 1]

    # Barklem damping for this Fe I line
    gammabk          = -7.280
    alphabk          =  0.238
    state.gambark[0] = 10.0 ** gammabk
    state.alpbark[0] = (1.0 - alphabk) / 2.0
    state.gamrad[0]  =  0.0

    opt = state.dampingopt
    if   opt == 0: state.damptype[0] = 'UNSLDc6'
    elif opt == 1: state.damptype[0] = 'BKgamma'
    elif opt == 2: state.damptype[0] = 'BLKWLc6'
    else:          state.damptype[0] = 'NEXTGEN '

    # Partition functions for all elements
    partfn(state)

    # Doppler widths, damping, kapnu0 for the fake line (nearly pass 3)
    state.lim1line = 0
    state.lim2line = 0
    state.nlines   = 1
    old_nf2out = getattr(state, '_nf2out_saved', None)  # suppress output
    nearly(state, 3)

    # COG parameters (rwstep=0.15, range [-6.7, -3.7])
    state.lim1   = 0
    state.lim2   = 0
    state.rwlow  = -6.7
    state.rwhigh = -3.7
    state.rwstep =  0.15

    _curve(state)

    # Build fine interpolation table at 0.005-dex steps in log(gf)
    # Covers the interior COG points using Catmull-Rom cubic interpolation,
    # matching the exact scheme in Fortran Fakeline.f.
    n_pts  = state.ncurve + 1  # total raw COG points
    gf_log = state.gf1[:n_pts]
    rw_log = state.w[:n_pts]

    m = 0
    for i in range(1, n_pts - 2):     # Fortran: do i=2,ncurve-2 (1-based)
        for l in range(30):
            pp = (1.0 / 30.0) * (l - 1)
            state.gftab[m] = gf_log[1] + 0.005 * m
            state.rwtab[m] = _catmull_rom(
                rw_log[i - 1], rw_log[i], rw_log[i + 1], rw_log[i + 2], pp)
            m += 1

    state.ntabtot = m
