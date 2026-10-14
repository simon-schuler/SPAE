"""EW linelist output formatting for Spectrum_Data.make_ew_doc().

Kept separate from spectrum_data.py -- and passed into make_ew_doc() as a
plain callable (its own line_formatter= parameter) rather than hardcoded
there -- so a different output format can be added (here or anywhere
else) without touching make_ew_doc() itself. Write a function with the
same signature as format_moog_ew_line() below and pass it as
line_formatter= to use it instead.
"""

import numpy as np


def format_moog_ew_line(wave, species, ep, gf, rad, ew, ew_err):
    """
    Default make_ew_doc() line formatter: the classic fixed-width
    MOOG-format linelist-with-EW line (wavelength, species, EP, log(gf),
    C6/D0 damping, EW, EW error) -- see moog/inlines.py's parse_linelist()
    for the matching reader, and README.md's "Linelist" section for the
    full column reference.

    Parameters
    ----------
    wave : float -- rest wavelength (A)
    species : float -- MOOG species code (e.g. 26.0 = Fe I)
    ep : float -- lower excitation potential (eV)
    gf : float -- log(gf) (or gf; MOOG detects the sign)
    rad : float -- van der Waals damping constant (C6) or molecular
        dissociation energy (D0)
    ew : float -- measured equivalent width (mA)
    ew_err : float -- EW uncertainty (mA)

    Returns
    -------
    str -- one formatted line, WITHOUT a trailing newline (make_ew_doc()
    adds it, along with any trailing comment).
    """
    wave_s = "  " + str(wave)
    elmnt = str(species)
    if species < 10:
        elmnt = '0' + elmnt
    while len(elmnt) < 6:
        elmnt = elmnt + '0'
    ep_s = str(ep)
    gf_s = str(gf)
    rad_s = str(rad)
    # No caller currently tracks a real per-line dissociation energy;
    # always written as 0.0, matching the pre-refactor inline formatter.
    d0_s = '0.0'
    ew_s = str(np.round(ew, 3))
    err_s = str(np.round(ew_err, 3))
    # 8 fixed 10-char columns: wave1,atom1,e,gf,dampnum,d0,width,width_err
    # (the first 7 match MOOG's native 7e10.3 linelist format; the trailing
    # err column is a pymoog-only extension read by moog/inlines.py -- real
    # MOOG ignores anything past column 70). Must stay 8x10 -- a prior
    # 7-field/variable-width version shifted every column past `rad` and
    # silently dropped the err column outside abfind_direct's read window.
    return "{0:10s}{1:10s}{2:10s}{3:10s}{4:10s}{5:10s}{6:10s}{7:10s}".format(
        wave_s, elmnt, ep_s, gf_s, rad_s, d0_s, ew_s, err_s)
