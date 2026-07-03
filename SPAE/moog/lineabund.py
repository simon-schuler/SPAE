"""
Iteratively determine the abundance for a single spectral line.
Translated from Lineabund.f.

The algorithm works by adjusting log(gf) — equivalently kapnu0 — until the
synthetic equivalent width matches the observed width, then converts the gf
shift to an abundance correction.
"""
import numpy as np

from .oneline import oneline


# ---------------------------------------------------------------------------
# COG lookup helper
# ---------------------------------------------------------------------------

def _interp_cog(state, rwlg: float) -> float:
    """
    Interpolate the fine COG table (gftab/rwtab) to find log(gf) for a
    given log(RW).  Returns the last table entry if rwlg is above the range.
    """
    ntot  = state.ntabtot
    rwtab = state.rwtab[:ntot]
    gftab = state.gftab[:ntot]
    # Linear interpolation: first crossing of rwtab > rwlg
    for i in range(1, ntot):
        if rwtab[i] > rwlg:
            frac = (rwlg - rwtab[i - 1]) / (rwtab[i] - rwtab[i - 1])
            return gftab[i - 1] + (gftab[i] - gftab[i - 1]) * frac
    return float(gftab[ntot - 1])


# ---------------------------------------------------------------------------
# lineabund
# ---------------------------------------------------------------------------

def lineabund(state, abundin: float) -> None:
    """
    Fit the abundance for line state.lim1 to match state.width[lim1].

    Reads:
      state.lim1      — 0-based index of the line to fit
      state.width[]   — observed equivalent widths [Å]
      state.gftab/rwtab/ntabtot — COG lookup table from fakeline()

    Writes:
      state.abundout[lim1]  — derived abundance (log epsilon scale)
      state.widout[lim1]    — final computed EW [Å]
      state.wid1comp[lim1]  — first-iteration computed EW [Å]
    """
    lim1 = state.lim1
    ntau = state.ntau

    state.lim2   = lim1
    state.ncurve = 0   # 0-based; Fortran starts at ncurve=1

    # Working gf starts from the input gf
    state.gf1[state.ncurve] = state.gf[lim1]

    # Observed RW → log(gf) from COG
    rwlgobs = np.log10(state.width[lim1] / state.wave1[lim1])
    gfobs   = _interp_cog(state, rwlgobs)

    ratio = 1.0   # initialise (used in final adjustment even on first exit)

    while True:
        # Compute the line with the current gf
        oneline(state, 1)

        rwlgcal = np.log10(state.w[state.ncurve] / state.wave1[lim1])
        gfcal   = _interp_cog(state, rwlgcal)

        error = (state.w[state.ncurve] - state.width[lim1]) / state.width[lim1]
        ratio = 10.0 ** (gfobs - gfcal)

        state.ncurve += 1

        if abs(error) >= 0.0015 and state.ncurve < 20:
            # Not yet converged: adjust gf (and proportionally kapnu0)
            rwlcomp = np.log10(state.w[state.ncurve - 1] / state.wave1[lim1])
            if rwlcomp > -4.7:
                # Saturated regime: slow convergence with sqrt of proposed shift
                adj = np.sqrt(ratio)
            else:
                adj = ratio
            state.gf1[state.ncurve]   = state.gf1[state.ncurve - 1] * adj
            state.kapnu0[lim1, :ntau] *= adj
        else:
            break

    # Final small gf adjustment and one last profile computation
    state.gf1[state.ncurve]   = state.gf1[state.ncurve - 1] * ratio
    state.kapnu0[lim1, :ntau] *= ratio
    oneline(state, 2)

    state.widout[lim1]   = state.w[state.ncurve]
    state.wid1comp[lim1] = state.w[0]   # first-try EW (Fortran w(1))
    diff = np.log10(state.gf1[state.ncurve] / state.gf[lim1])
    state.abundout[lim1] = abundin + diff
