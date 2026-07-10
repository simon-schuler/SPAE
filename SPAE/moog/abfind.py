"""
Abundance determination from equivalent widths (abfind mode).
Translated from Abfind.f, Linlimit.f, and Molquery.f.

Entry point: abfind(state) — I/O-free; returns a structured results dict.

Caller is responsible for populating state before calling:
  state.fmodel, state.flines, state.fparam  — file paths
  inmodel(state), inlines(state, 1), eqlib(state) must already have run, OR
  call abfind() which handles the full pipeline.

The public function abfind_from_files() is a convenience wrapper that reads
the parameter file, model atmosphere, and line list from disk.
"""
import numpy as np

from .params    import params
from .inmodel   import inmodel, inmodel_from_array
from .inlines   import inlines, parse_linelist, apply_parsed_lines
from .eqlib     import eqlib
from .nearly    import nearly
from .fakeline  import fakeline
from .lineabund import lineabund
from .stats     import stats
from .eqlib     import _sunder


# ---------------------------------------------------------------------------
# Species-range iterator (Linlimit.f, mode 2)
# ---------------------------------------------------------------------------

def _linlimit(state) -> None:
    """
    Advance lim1line/lim2line to the next group of same-species lines.

    Uses lim2line == -1 as the "not started" sentinel (set by the driver
    before the first call).  Sets lim1line = -1 when all species are done.
    """
    nlines = state.nlines

    # All species exhausted after previous call
    if state.lim2line == nlines - 1:
        state.lim1line = -1
        return

    # Advance pointer (lim2line < 0 = first call)
    if state.lim2line < 0:
        state.lim1line = 0
    else:
        state.lim1line = state.lim2line + 1

    # Single last line
    if state.lim1line == nlines - 1:
        state.lim2line = state.lim1line
        return

    # Find end of this species group (same atom1 value)
    old_atom = state.atom1[state.lim1line]
    for j in range(state.lim1line + 1, nlines):
        if state.atom1[j] != old_atom:
            state.lim2line = j - 1
            return

    state.lim2line = nlines - 1


# ---------------------------------------------------------------------------
# Molecular equilibrium query (Molquery.f)
# ---------------------------------------------------------------------------

def _molquery(state) -> None:
    """
    Determine whether the current species is involved in molecular
    equilibrium and set state.molflag and state.iabatom accordingly.

    Translates Molquery.f faithfully, including the known Fortran quirk
    where the atomic-species loop returns after checking only iorder[0].

    For non-hydride molecules the user would normally be prompted for which
    atom to vary; here we default to the heavier atom (non-H partner).
    """
    state.molflag = 0
    state.iabatom = int(state.atom1[state.lim1obs] + 0.0001)

    if state.atom1[state.lim1line] < 100.0:
        # Atomic species: check molecular equilibrium list
        if state.neq == 0:
            return
        # Fortran has a `return` inside the do-loop, so only iorder[0] is checked
        if state.iabatom == state.iorder[0]:
            state.molflag = 1
        return

    # Molecular species
    if state.neq == 0:
        raise RuntimeError(
            f"Molecular equilibrium not computed but species {state.iabatom} is a molecule.")

    iaa, ibb = _sunder(state.atom1[state.lim1obs])
    state.iaa = int(iaa)
    state.ibb = int(ibb) if ibb else 0

    found = False
    for n in range(state.neq):
        if state.iaa == state.iorder[n] or state.ibb == state.iorder[n]:
            found = True
            break
    if not found:
        raise RuntimeError(
            f"Molecular equilibrium does not include atoms for species {state.iabatom}.")

    state.molflag = 1

    # Identify which atom's abundance to vary
    if state.iaa == 1:
        state.iabatom = state.ibb   # hydride: vary non-H partner
    elif state.ibb == 1:
        state.iabatom = state.iaa   # hydride: vary non-H partner
    else:
        # Non-hydride: default to first (heavier) atom; interactive choice skipped
        state.iabatom = state.iaa


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def abfind(state, parsed_lines=None) -> dict:
    """
    Derive elemental abundances from equivalent widths.

    Assumes state is fully populated (model atmosphere, line list, and
    molecular equilibrium already loaded).  Calls fakeline() → nearly(1) →
    per-species lineabund() → stats().

    Parameters
    ----------
    parsed_lines : dict or None
        Pre-parsed linelist snapshot from parse_linelist().  If provided,
        line data is restored from the dict instead of re-reading the file
        after fakeline() — eliminates one file read per call.

    Returns
    -------
    dict with keys:
      'lines'   : list of per-line result dicts
      'species' : dict keyed by atom1 value → species statistics dict
    """
    state.mode = 2

    # Build curve-of-growth lookup table (clobbers nlines=1 and wave1[0])
    fakeline(state)

    # Restore real line data after fakeline — from memory if pre-parsed,
    # otherwise re-read from disk (original behaviour)
    if parsed_lines is not None:
        apply_parsed_lines(state, parsed_lines)
    else:
        inlines(state, 1)

    # Doppler widths, damping, line-centre opacities for all lines
    state.waveold = 0.0   # force continuum recompute on first line
    nearly(state, 1)

    lines_out   = []
    species_out = {}

    # Iterate over species groups
    state.lim2line = -1   # "not started" sentinel for _linlimit
    while True:
        _linlimit(state)
        if state.lim1line < 0:
            break

        state.lim1obs = state.lim1line
        state.lim2obs = state.lim2line

        _molquery(state)

        if state.molflag == 0:
            # Atomic species not in molecular equilibrium
            abundin = np.log10(state.xabund[state.iabatom - 1]) + 12.0
            for lim1 in range(state.lim1line, state.lim2line + 1):
                state.lim1 = lim1
                lineabund(state, abundin)
            sp_stats = stats(state)

        else:
            # Species involved in molecular equilibrium: iterate until convergence
            abundin = np.log10(state.xabund[state.iabatom - 1]) + 12.0
            for _iter in range(6):
                for lim1 in range(state.lim1line, state.lim2line + 1):
                    state.lim1 = lim1
                    lineabund(state, abundin)
                sp_stats = stats(state)

                T_tau5 = state.t[state.jtau5]
                iatom  = int(state.atom1[state.lim1line] + 0.0001)
                converged = (
                    T_tau5 >= 3800.0 and
                    iatom not in (6, 8) and
                    iatom < 100
                )
                if converged or abs(sp_stats['average'] - abundin) <= 0.02:
                    break
                # Update abundance and redo eqlib + nearly
                state.xabund[state.iabatom - 1] = 10.0 ** (sp_stats['average'] - 12.0)
                abundin = sp_stats['average']
                eqlib(state)
                nearly(state, 2)

        # Collect per-line results
        for l in range(state.lim1obs, state.lim2obs + 1):
            lines_out.append({
                'wave':    float(state.wave1[l]),
                'species': float(state.atom1[l]),
                'ep':      float(state.e[l, 0]),
                'loggf':   float(np.log10(state.gf[l])),
                'ew_obs':  float(state.width[l]) * 1000.0,    # mÅ
                'ew_calc': float(state.widout[l]) * 1000.0,   # mÅ
                'abund':   float(state.abundout[l]),
                'delavg':  float(state.abundout[l] - sp_stats['average']),
            })

        species_key = float(state.atom1[state.lim1obs])
        species_out[species_key] = {
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

    return {'lines': lines_out, 'species': species_out}


# ---------------------------------------------------------------------------
# Convenience wrapper: read files then run
# ---------------------------------------------------------------------------

def abfind_from_files(state) -> dict:
    """
    Full pipeline: read params → model → lines → eqlib → abfind.

    state.fparam must point to a valid batch.par before calling.
    """
    params(state, state.fparam)
    inmodel(state)
    inlines(state, 1)
    eqlib(state)
    return abfind(state)


def abfind_direct(state, atmos_array, feh, vt_kms, linelist,
                  extra_overrides=None) -> dict:
    """
    Full abfind pipeline without star.mod or batch.par file I/O.

    Parameters
    ----------
    state          : State
    atmos_array    : ndarray, shape (N, 7) returned by atmos.atmos()
    feh            : float  — [Fe/H]
    vt_kms         : float  — microturbulence in km/s
    linelist       : str or dict
        Either a path to the MOOG line list file (str), or a pre-parsed
        snapshot dict from parse_linelist() for zero file I/O per call.
    extra_overrides : dict {Z: logeps} or None
        Passed through to inmodel_from_array(); defaults to {3: 3.30} (Li).
    """
    # Reset abundance overrides (mirrors params._init_defaults)
    state.numpecatom     = 0
    state.numatomsyn     = 0
    state.ninetynineflag = 0
    state.pec[:]         = 0
    state.pecabund[:]    = 0.0
    state.abfactor[:]    = 0.0

    # Accept either a filename (str) or a pre-parsed dict
    if isinstance(linelist, dict):
        parsed = linelist
    else:
        state.flines = linelist
        parsed = parse_linelist(linelist)

    inmodel_from_array(state, atmos_array, feh, vt_kms, extra_overrides)
    apply_parsed_lines(state, parsed)
    eqlib(state)
    return abfind(state, parsed_lines=parsed)
