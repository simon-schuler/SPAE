"""
Abundance derivation from blended spectral features.
Translated from Blends.f and Total.f.

Entry points
------------
blends(state)            — I/O-free; returns result dict
blends_from_files(state) — reads params/model/lines then calls blends()
"""
import numpy as np

from .params     import params
from .inmodel    import inmodel
from .inlines    import inlines, _sunder
from .eqlib      import eqlib
from .nearly     import nearly
from .opacit     import opacit
from .cdcalc     import cdcalc
from .taukap     import taukap
from .stats      import stats
from .math_utils import rinteg


# ---------------------------------------------------------------------------
# Blend-group linlimit  (Linlimit.f mode 4, 0-indexed)
# ---------------------------------------------------------------------------

def _linlimit_blends(state) -> None:
    """
    Advance lim1line/lim2line to the next blend group (mode=4 linlimit).

    Groups are delimited by state.group[j]:
      0 = start of a new blend group (positive wave1 in the original linelist)
      1 = member of the preceding group (negative wave1 in the original linelist)

    Uses the same 0-indexed convention as abfind._linlimit.
    Raises RuntimeError if the first line of a group is incorrectly marked.
    """
    nlines = state.nlines + state.nstrong
    lim1   = state.lim1line

    if lim1 < 0:       # "not started" sentinel → first line
        lim1 = 0

    if state.group[lim1] != 0:
        raise RuntimeError(
            f"Line {lim1} (wave={state.wave1[lim1]:.3f} Å) is marked as a "
            "blend member but is expected to start a new group.  "
            "Check that the linelist has the correct negative-wave encoding.")

    state.lim1line = lim1

    if lim1 == nlines - 1:
        state.lim2line = lim1
        return

    for j in range(lim1 + 1, nlines):
        if state.group[j] != 1:
            state.lim2line = j - 1
            return

    state.lim2line = nlines - 1


# ---------------------------------------------------------------------------
# Synthesis loop for a fixed line range  (inner part of Synspec.f, mode=4)
# ---------------------------------------------------------------------------

def _synspec_blends(state, lim1: int, lim2: int) -> tuple:
    """
    Compute a synthetic spectrum over [state.start, state.sstop] using
    ALL lines in [lim1, lim2] at every wavelength step (no per-step linlimit).

    Returns (wave_arr, depth_arr) — same convention as synspec().
    """
    ntau  = state.ntau
    start = state.start
    sstop = state.sstop
    step  = state.step

    kount = int(round((sstop - start + step / 4.0) / step)) + 1

    wave_arr = np.empty(kount, dtype=np.float64)
    d_arr    = np.empty(kount, dtype=np.float64)

    state.lim1 = lim1
    state.lim2 = lim2
    wavl = 0.0

    for n in range(kount):
        wave = start + n * step
        state.wave = wave
        wave_arr[n] = wave

        # Recompute continuum when wavelength shifts by ≥ 0.1%
        if abs(wave - wavl) / wave >= 0.001:
            wavl = wave
            opacit(state, 2, wave)
            cdcalc(state, 1)
            first = 0.4343 * state.cd[0]
            state.flux, _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

        taukap(state)
        cdcalc(state, 2)
        first = 0.4343 * state.cd[0]
        d_arr[n], _ = rinteg(state.xref[:ntau], state.cd[:ntau], ntau, first)

    return wave_arr, d_arr


# ---------------------------------------------------------------------------
# kapnu0 scaling  (shared by the iteration and the final step)
# ---------------------------------------------------------------------------

def _scale_kapnu0(state, iatom: int, lim1: int, lim2: int, ratio: float) -> None:
    """Scale kapnu0 by ratio for every line of iatom in [lim1, lim2]."""
    ntau = state.ntau
    for j in range(lim1, lim2 + 1):
        if state.atom1[j] >= 100.0:
            ia, ib = _sunder(state.atom1[j])
            if ia == iatom or ib == iatom:
                state.kapnu0[j, :ntau] *= ratio
        elif int(state.atom1[j] + 0.0001) == iatom:
            state.kapnu0[j, :ntau] *= ratio


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def blends(state) -> dict:
    """
    Derive the abundance of one element from blended spectral features.

    The element to fit is state.cogatom (set by 'coglimits' in batch.par).
    Blend groups are encoded in the linelist by negative wave1 values for
    group members (positive wave1 = group leader).

    Requires state to be fully populated:
      inlines(state, 1), eqlib(state), nearly(state, 1) must already have run.

    Returns
    -------
    dict with keys:
      'lines'   : list of per-feature result dicts
      'species' : statistics dict for the fitted element
    """
    iatom = int(state.cogatom + 0.0001)
    if iatom == 0:
        raise ValueError(
            "state.cogatom == 0: set 'coglimits' in batch.par to specify "
            "the element (atomic number) whose abundance to fit.")

    state.mode = 4

    # Initialise blend-group iteration
    state.lim1line = -1    # "not started" sentinel
    state.lim2line = -1

    lines_out = []

    for _feature in range(1000):
        # Exit when all lines have been processed
        if state.lim2line >= state.nlines + state.nstrong - 1:
            break

        _linlimit_blends(state)
        lim1 = state.lim1line    # 0-indexed, inclusive
        lim2 = state.lim2line    # 0-indexed, inclusive

        # Check: the target element must have a line in this group
        ifind = False
        for j in range(lim1, lim2 + 1):
            if state.atom1[j] >= 100.0:
                ia, ib = _sunder(state.atom1[j])
                if ia == iatom or ib == iatom:
                    ifind = True
                    break
            elif int(state.atom1[j] + 0.0001) == iatom:
                ifind = True
                break

        if not ifind:
            for j in range(lim1, lim2 + 1):
                state.abundout[j] = 999.99
            state.lim1line = state.lim2line + 1
            if state.lim1line >= state.nlines + state.nstrong:
                break
            continue

        # Set the synthesis window for this blend feature
        state.start    = state.wave1[lim1] - state.delwave
        state.sstop    = state.wave1[lim2] + state.delwave
        state.oldstart = state.start
        state.oldstop  = state.sstop

        ew_obs  = state.width[lim1]     # observed EW [Å]
        rwlgobs = np.log10(ew_obs / state.wave1[lim1])

        # gf1 accumulates the cumulative kapnu0 scaling factor
        ncurve = 1
        state.gf1[ncurve] = 1.0
        ratio = 1.0

        for k in range(30):
            wave_arr, d_arr = _synspec_blends(state, lim1, lim2)
            ew_pred, _ = rinteg(wave_arr, d_arr, len(wave_arr), 0.0)
            state.w[ncurve] = ew_pred

            if ew_pred <= 0.0:
                break   # degenerate feature; skip

            error = (ew_pred - ew_obs) / ew_obs
            ratio = ew_obs / ew_pred
            ncurve += 1

            if abs(error) >= 0.0075:
                rwlcomp = np.log10(state.w[ncurve - 1] / state.wave1[lim1])
                if   rwlcomp < -5.2 and rwlgobs < -5.2:
                    pass                # linear regime: ratio unchanged
                elif rwlcomp >= -5.2 and rwlgobs >= -5.2:
                    ratio = ratio ** 2.0  # flat part: stronger correction
                else:
                    ratio = ratio ** 1.5  # mixed regime

                state.gf1[ncurve] = state.gf1[ncurve - 1] * ratio
                _scale_kapnu0(state, iatom, lim1, lim2, ratio)

                if k == 19:
                    raise RuntimeError(
                        f"blends: 20 iterations without convergence for "
                        f"feature at {state.wave1[lim1]:.3f} Å")
            else:
                break

        # Final synthesis at converged kapnu0
        state.gf1[ncurve] = state.gf1[ncurve - 1] * ratio
        _scale_kapnu0(state, iatom, lim1, lim2, ratio)
        wave_arr, d_arr = _synspec_blends(state, lim1, lim2)
        ew_pred, _      = rinteg(wave_arr, d_arr, len(wave_arr), 0.0)
        state.w[ncurve]      = ew_pred
        state.widout[lim1]   = ew_pred

        diff = np.log10(state.gf1[ncurve])
        state.abundout[lim1] = np.log10(state.xabund[iatom - 1]) + 12.0 + diff

        # For multi-line blends: assign abundance to the strongest target line
        if lim2 > lim1:
            abunblend    = state.abundout[lim1]
            widblend     = state.widout[lim1]
            for j in range(lim1, lim2 + 1):
                state.abundout[j] = 999.99
            strongest    = 0.0
            linstrongest = lim1
            for j in range(lim1, lim2 + 1):
                if state.atom1[j] >= 100.0:
                    ia, ib = _sunder(state.atom1[j])
                    if ia == iatom or ib == iatom:
                        if state.kapnu0[j, state.jtau5] > strongest:
                            strongest    = state.kapnu0[j, state.jtau5]
                            linstrongest = j
                elif int(state.atom1[j] + 0.0001) == iatom:
                    if state.kapnu0[j, state.jtau5] > strongest:
                        strongest    = state.kapnu0[j, state.jtau5]
                        linstrongest = j
            state.abundout[linstrongest] = abunblend
            state.widout[linstrongest]   = widblend
        else:
            linstrongest = lim1

        lines_out.append({
            'wave':     float(state.wave1[lim1]),
            'wave_hi':  float(state.wave1[lim2]),
            'species':  float(state.atom1[linstrongest]),
            'ep':       float(state.e[linstrongest, 0]),
            'loggf':    float(np.log10(state.gf[linstrongest])),
            'ew_obs':   float(ew_obs) * 1000.0,
            'ew_calc':  float(state.widout[linstrongest]) * 1000.0,
            'abund':    float(state.abundout[linstrongest]),
            'n_lines':  lim2 - lim1 + 1,
            'n_iter':   ncurve,
        })

        state.lim1line = state.lim2line + 1
        if state.lim1line >= state.nlines + state.nstrong:
            break

    # Statistics over all lines (999.99 excluded by stats())
    state.lim1obs = 0
    state.lim2obs = state.nlines + state.nstrong - 1
    sp_stats = stats(state)

    return {
        'lines':   lines_out,
        'species': {
            float(iatom): {
                'average': sp_stats['average'],
                'deviate': sp_stats['deviate'],
                'n':       sp_stats['kount'],
                'ep_slope':     sp_stats['ep_slope'],
                'ep_intercept': sp_stats['ep_intercept'],
                'ep_r':         sp_stats['ep_r'],
                'rw_slope':     sp_stats['rw_slope'],
                'rw_intercept': sp_stats['rw_intercept'],
                'rw_r':         sp_stats['rw_r'],
                'wv_slope':     sp_stats['wv_slope'],
                'wv_intercept': sp_stats['wv_intercept'],
                'wv_r':         sp_stats['wv_r'],
            }
        },
    }


# ---------------------------------------------------------------------------
# Convenience wrapper: read files then run
# ---------------------------------------------------------------------------

def blends_from_files(state) -> dict:
    """
    Full pipeline: read params → model → lines → eqlib → nearly → blends.

    state.fparam must point to a valid batch.par (or equivalent) before calling.
    The batch.par must specify:
      - 'coglimits' with the target element's atomic number in the 5th field
      - 'blenlimits' with three values: delwave, step, cogatom
        (cogatom is the atomic number Z of the element to fit)
    """
    params(state, state.fparam)
    inmodel(state)
    inlines(state, 1)
    eqlib(state)
    nearly(state, 1)
    return blends(state)
