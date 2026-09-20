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
    ew_s = str(np.round(ew, 3))
    err_s = str(np.round(ew_err, 3))
    return "{0:14s}{1:11s}{2:8s}{3:15s}{4:17s}{5:10s}{6:5s}".format(
        wave_s, elmnt, ep_s, gf_s, rad_s, ew_s, err_s)
