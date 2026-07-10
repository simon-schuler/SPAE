"""
Line list reader translated from Inlines.f.

Public API
----------
inlines(state, num=1)
    Read spectral lines from state.flines into the state Linex.com arrays.
"""

import numpy as np
from .atomic_data import XAM, XCHI1, XCHI2, XCHI3


def _sunder(amol: float):
    """
    Split a MOOG molecule species code into two constituent atomic numbers.

    Translated from Sunder.f.  The encoding packs the two Z values as
    integers: e.g. 608 → C(6) + O(8), 106 → H(1) + C(6).

    Returns
    -------
    (ia, ib) : (int, int)
        Left and right atomic numbers (ia < ib guaranteed by inlines check).
    """
    im = round(amol)
    for divisor in (100_000_000, 1_000_000, 10_000, 100, 1):
        i1 = im // divisor
        if i1:
            i2 = im - i1 * divisor
            return i1, i2
    return 0, 0


def _read_fixed(line: str) -> list:
    """
    Parse one line in Fortran 7e10.3 format (7 × 10-char E fields).
    Blank fields yield 0.0; handles Fortran D/d exponent notation.
    """
    line = line.rstrip('\n')
    vals = []
    for start in range(0, 70, 10):
        field = line[start:start + 10] if start < len(line) else ''
        field = field.replace('d', 'e').replace('D', 'E').strip()
        try:
            vals.append(float(field) if field else 0.0)
        except ValueError:
            vals.append(0.0)
    while len(vals) < 7:
        vals.append(0.0)
    return vals


def _read_free(line: str) -> list:
    """Parse one free-format line into up to 7 floats (missing → 0.0)."""
    toks = line.split()
    vals = [float(t) for t in toks[:7]]
    while len(vals) < 7:
        vals.append(0.0)
    return vals


def _parse_line(line: str, linfileopt: int) -> list:
    return _read_fixed(line) if linfileopt == 0 else _read_free(line)


def inlines(state, num: int = 1) -> None:
    """
    Read spectral line data from state.flines into state.

    After this call the following State Linex.com fields are populated:
        wave1, atom1, e[:,0/1], gf, dampnum, d0, width, charge,
        amass, rdmass, chi, group, nlines, nstrong.

    Parameters
    ----------
    state : State
        MOOG State dataclass (flines, fslines, linfileopt, dostrong,
        gfstyle, iunits, and control must already be set by params).
    num : int
        Read mode: 1=first read, 5=synth (same path), 6=COG (skip title).
    """
    linfileopt = state.linfileopt
    is_blends  = state.control.strip() == 'blends'

    # Apply peculiar element abundance overrides (Inlines.f lines 27-32).
    # num==2 and num==6 skip this block in the Fortran (goto 4 / goto 340).
    # pecabund[z, k] is a log offset from xabu; applied before eqlib/nearly.
    if state.numpecatom > 0 and num not in (2, 6):
        k = max(0, state.isynth - 1)   # isynth is 1-indexed
        for z in range(3, 96):         # elements 3–95 (1-indexed), as in Fortran
            state.xabund[z - 1] = (
                10.0 ** state.pecabund[z - 1, k]
                * 10.0 ** state.abfactor[k]
                * state.xabu[z - 1]
            )

    # ------------------------------------------------------------------ #
    # 1. Strong lines (max 40)                                             #
    # ------------------------------------------------------------------ #
    swave1   = []; satom1  = []; se      = []; sgf      = []
    sdampnum = []; sd0     = []; swidth  = []; scharge  = []

    if state.dostrong > 0:
        with open(state.fslines) as sf:
            for j in range(41):
                raw = sf.readline()
                if not raw:
                    break
                vals = _parse_line(raw, linfileopt)
                iatom = int(vals[1])
                chg   = 1.0 + int(10.0 * (vals[1] - iatom) + 0.0001)
                if chg > 3:
                    raise ValueError(
                        f"Strong line λ={vals[0]:.3f} atom={vals[1]}: "
                        "triple ion or higher not supported"
                    )
                swave1.append(vals[0]);   satom1.append(vals[1])
                se.append(vals[2]);       sgf.append(vals[3])
                sdampnum.append(vals[4]); sd0.append(vals[5])
                swidth.append(vals[6]);   scharge.append(chg)
        if len(swave1) > 40:
            raise ValueError("Strong line list has more than 40 lines.")

    nstrong = len(swave1)

    # ------------------------------------------------------------------ #
    # 2. Main line list                                                    #
    # ------------------------------------------------------------------ #
    wave1_l = []; atom1_l = []; e_l   = []; gf_l   = []
    damp_l  = []; d0_l    = []; wid_l = []; chg_l  = []

    with open(state.flines) as fh:
        # Title line (num==6 skips it; COG mode doesn't rewind and re-read)
        if num != 6:
            state.linitle = fh.readline().rstrip('\n')

        max_main = 2500 - nstrong
        for _ in range(max_main):
            raw = fh.readline()
            if not raw:
                break
            vals  = _parse_line(raw, linfileopt)
            iatom = int(vals[1])
            chg   = 1.0 + int(10.0 * (vals[1] - iatom) + 0.0001)
            if chg > 3:
                raise ValueError(
                    f"Line λ={vals[0]:.3f} atom={vals[1]}: "
                    "triple ion or higher not supported"
                )
            # Skip lines whose EW < 0 (log-RW format) in non-blends modes
            if vals[6] < 0.0 and not is_blends:
                continue

            if state.iunits == 1:
                vals[0] *= 1.0e4   # microns → Angstroms

            wave1_l.append(vals[0]); atom1_l.append(vals[1])
            e_l.append(vals[2]);     gf_l.append(vals[3])
            damp_l.append(vals[4]); d0_l.append(vals[5])
            wid_l.append(vals[6]);   chg_l.append(chg)

    nlines = len(wave1_l)
    total  = nlines + nstrong

    state.nlines  = nlines
    state.nstrong = nstrong

    # Concatenate main + strong
    all_wave1   = wave1_l + swave1
    all_atom1   = atom1_l + satom1
    all_e       = e_l     + se
    all_gf      = gf_l    + sgf
    all_dampnum = damp_l  + sdampnum
    all_d0      = d0_l    + sd0
    all_width   = wid_l   + swidth
    all_charge  = chg_l   + scharge

    for j in range(total):
        state.wave1[j]   = all_wave1[j]
        state.atom1[j]   = all_atom1[j]
        state.e[j, 0]    = all_e[j]
        state.gf[j]      = all_gf[j]
        state.dampnum[j] = all_dampnum[j]
        state.d0[j]      = all_d0[j]
        state.width[j]   = all_width[j]
        state.charge[j]  = all_charge[j]

    # ------------------------------------------------------------------ #
    # 3. Post-processing                                                   #
    # ------------------------------------------------------------------ #
    _postprocess(state, total)


def _snapshot_lines(state) -> dict:
    """Return a copy of all line arrays currently loaded in state."""
    total = state.nlines + state.nstrong
    return {
        'nlines':  state.nlines,
        'nstrong': state.nstrong,
        'linitle': state.linitle,
        'wave1':   state.wave1[:total].copy(),
        'atom1':   state.atom1[:total].copy(),
        'e':       state.e[:total, :].copy(),
        'gf':      state.gf[:total].copy(),
        'dampnum': state.dampnum[:total].copy(),
        'd0':      state.d0[:total].copy(),
        'width':   state.width[:total].copy(),
        'charge':  state.charge[:total].copy(),
        'group':   state.group[:total].copy(),
        'amass':   state.amass[:total].copy(),
        'rdmass':  state.rdmass[:total].copy(),
        'chi':     state.chi[:total, :].copy(),
    }


def apply_parsed_lines(state, parsed: dict) -> None:
    """Copy a pre-parsed line snapshot back into state (no file I/O)."""
    total = parsed['nlines'] + parsed['nstrong']
    state.nlines   = parsed['nlines']
    state.nstrong  = parsed['nstrong']
    state.linitle  = parsed['linitle']
    state.wave1[:total]   = parsed['wave1']
    state.atom1[:total]   = parsed['atom1']
    state.e[:total, :]    = parsed['e']
    state.gf[:total]      = parsed['gf']
    state.dampnum[:total] = parsed['dampnum']
    state.d0[:total]      = parsed['d0']
    state.width[:total]   = parsed['width']
    state.charge[:total]  = parsed['charge']
    state.group[:total]   = parsed['group']
    state.amass[:total]   = parsed['amass']
    state.rdmass[:total]  = parsed['rdmass']
    state.chi[:total, :]  = parsed['chi']


def parse_linelist(filename: str, linfileopt: int = 0, dostrong: int = 0,
                   iunits: int = 0, gfstyle: int = 0) -> dict:
    """
    Parse a MOOG linelist file once and return a snapshot dict.

    The dict can be passed to apply_parsed_lines() to restore line data into
    any State without reading the file again.  Intended for use in MCMC loops
    where the linelist is fixed across thousands of likelihood evaluations.

    Parameters
    ----------
    filename    : str  — path to the MOOG line list file
    linfileopt  : int  — 0=fixed format (default), 1=free format
    dostrong    : int  — 0=no strong lines (default)
    iunits      : int  — 0=Angstroms (default), 1=microns
    gfstyle     : int  — 0=log gf in file (default)
    """
    from .state import State
    s = State()
    s.flines     = filename
    s.linfileopt = linfileopt
    s.dostrong   = dostrong
    s.iunits     = iunits
    s.gfstyle    = gfstyle
    s.control    = 'abfind '
    s.numpecatom = 0
    inlines(s, 1)
    return _snapshot_lines(s)


def _postprocess(state, total: int) -> None:
    """Apply all Inlines.f post-read conversions and derived quantities."""

    # ---- group detection: negative wave1 marks a blend member ----
    for j in range(total):
        if state.wave1[j] < 0.0:
            state.group[j] = 1
            state.wave1[j] = abs(state.wave1[j])
            if j > 0:
                state.width[j] = state.width[j - 1]
        else:
            state.group[j] = 0

    # ---- excitation potential: cm^-1 → eV if any e > 50 ----
    for j in range(total):
        if state.e[j, 0] > 50.0:
            for jj in range(total):
                state.e[jj, 0] *= 1.2389e-4
            break

    # ---- gf: log(gf) → gf if any gf < 0 or gfstyle == 0 ----
    for j in range(total):
        if state.gfstyle == 0 or state.gf[j] < 0.0:
            for jj in range(total):
                state.gf[jj] = 10.0 ** state.gf[jj]
            break

    # ---- EW units: log(RW) → Å, or mÅ → Å ----
    for j in range(total):
        if state.width[j] < 0.0:
            state.width[j] = (10.0 ** state.width[j]) * state.wave1[j]
        else:
            state.width[j] /= 1000.0

    # ---- derived quantities per line ----
    for j in range(total):
        iatom  = int(state.atom1[j])
        atom10 = 10.0 * state.atom1[j]
        # Upper excitation potential from lower + photon energy
        state.e[j, 1] = state.e[j, 0] + 1.239e4 / state.wave1[j]

        if iatom >= 100:
            # ---- molecular line ----
            ia, ib = _sunder(state.atom1[j])
            if ia > ib:
                raise ValueError(
                    f"Molecular code {iatom}: constituent Z values in wrong "
                    f"order (Z1={ia} > Z2={ib})"
                )
            # Fractional part of atom10 signals isotopic mass encoding
            frac = atom10 - int(atom10)
            if frac <= 1e-10:
                # Standard masses from atomic data table
                mas1 = XAM[ia - 1]
                mas2 = XAM[ib - 1]
                state.amass[j] = mas1 + mas2
            else:
                # Isotopic masses packed into atom10 decimal digits
                jat100   = int(100.0  * (atom10 + 1e-5))
                mas1     = float(jat100 - 100 * int(atom10))
                jat10000 = int(10000.0 * (atom10 + 1e-5))
                mas2     = float(jat10000 - 100 * jat100)
                state.amass[j] = mas1 + mas2

            if state.d0[j] == 0.0:
                # Look up dissociation energy from datmol table
                found = False
                for k in range(110):
                    if round(state.datmol[0, k]) == round(state.atom1[j]):
                        state.d0[j] = state.datmol[1, k]
                        found = True
                        break
                if not found:
                    raise ValueError(
                        f"Unknown molecule {state.atom1[j]:.1f}: "
                        "no dissociation energy in datmol"
                    )

            state.rdmass[j]  = mas1 * mas2 / state.amass[j]
            state.chi[j, 0]  = 0.0
            state.chi[j, 1]  = 0.0
            state.chi[j, 2]  = 0.0

        else:
            # ---- atomic line ----
            frac = atom10 - int(atom10)
            if frac <= 1e-10:
                # Standard atomic mass (ionisation state is encoded, not mass)
                state.amass[j] = XAM[iatom - 1]
            else:
                # Isotopic mass encoded in fractional digits
                atom10_adj     = atom10 + 1e-5
                state.amass[j] = float(int(1000.0 * (atom10_adj - int(atom10_adj))))

            state.rdmass[j] = 0.0
            state.chi[j, 0] = XCHI1[iatom - 1]
            state.chi[j, 1] = XCHI2[iatom - 1]
            state.chi[j, 2] = XCHI3[iatom - 1]
