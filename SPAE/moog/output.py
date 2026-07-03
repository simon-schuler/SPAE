"""
MOOG-compatible output writers for the pymoog CLI.

Standard output (moog_out.1 / standard_out) receives a run header.
Summary output  (moog_out.2 / summary_out) receives per-mode data tables
in the same format as MOOGSILENT (Lineinfo.f formats 3001-3011, etc.).
"""
import math
from .atomic_data import ELEMENT_NAMES

_NO_FILE = 'no_filename_given'


def _open_out(path: str):
    """Open an output file; return None for 'no_filename_given'."""
    if not path or path == _NO_FILE:
        return None
    return open(path, 'w', encoding='utf-8')


def _species_parts(species: float):
    """Return (symbol, roman_numeral, iatom) for a MOOG species float code."""
    iatom = int(species + 0.0001)
    ion   = round((species - iatom) * 10)
    name  = ELEMENT_NAMES[iatom - 1].strip() if 1 <= iatom <= len(ELEMENT_NAMES) else '??'
    roman = ('I', 'II', 'III', 'IV')[min(ion, 3)]
    return name, roman, iatom


def _write_abund_block(f2, state, lines, species_dict, match_all_ions=False):
    """
    Abundance table in Lineinfo.f format (3001-3011).
    Used by abfind and blends.
    match_all_ions=True: match lines by atomic number (blends, cogatom mode).
    match_all_ions=False: match by exact species code (abfind).
    """
    # format 3001: a80 title lines
    f2.write(f"{(state.moditle or '')[:80]:<80}\n")
    f2.write(f"{(state.linitle or '')[:80]:<80}\n")

    if state.abscale != 0.0:
        # Inlines.f format 1006
        f2.write(
            f"ALL abundances NOT listed below differ from solar by"
            f" {state.abscale:6.2f} dex\n"
        )

    for sp_key in sorted(species_dict.keys()):
        sp = species_dict[sp_key]
        name, roman, iatom = _species_parts(sp_key)
        input_abund = math.log10(state.xabund[iatom - 1]) + 12.0

        # format 3002: a2 element + a4 ionization (" I  ", " II ", etc.)
        f2.write(
            f"Abundance Results for Species {name:2s} {roman:<3s}"
            f"       (input abundance = {input_abund:7.3f})\n"
        )
        # format 3003
        f2.write(
            "wavelength         ID      EP   logGF"
            "     EWin   logRWin     abund   delavg\n"
        )

        if match_all_ions:
            sp_lines = [l for l in lines if int(l['species'] + 0.0001) == iatom]
        else:
            sp_lines = [l for l in lines if abs(l['species'] - sp_key) < 0.001]
        for ln in sp_lines:
            ew = ln.get('ew_obs', 0.0)
            logrw = math.log10(ew / ln['wave'] / 1000.0) if ew > 0.0 else 0.0
            delavg = ln.get('delavg', ln['abund'] - sp['average'])
            # format 3007: f10.3 f11.5 f8.3 f8.3 f9.2 f10.3 f10.3 f9.3
            f2.write(
                f"{ln['wave']:10.3f}{ln['species']:11.5f}{ln['ep']:8.3f}"
                f"{ln['loggf']:8.3f}{ew:9.2f}{logrw:10.3f}"
                f"{ln['abund']:10.3f}{delavg:9.3f}\n"
            )

        # format 3008
        f2.write(
            f"average abundance = {sp['average']:6.3f}     "
            f"std. deviation = {sp['deviate']:6.3f}     "
            f"#lines = {sp['n']:3d}\n"
        )

        if sp.get('ep_slope') is not None:
            # format 3009
            f2.write(
                f"E.P. correlation:  slope = {sp['ep_slope']:7.3f}"
                f"  intercept = {sp['ep_intercept']:7.3f}"
                f"  corr. coeff. = {sp['ep_r']:7.3f}\n"
            )
            # format 3010
            f2.write(
                f"R.W. correlation:  slope = {sp['rw_slope']:7.3f}"
                f"  intercept = {sp['rw_intercept']:7.3f}"
                f"  corr. coeff. = {sp['rw_r']:7.3f}\n"
            )
            # format 3011: 1pd11.3 → scientific notation 11 wide, 3 decimal
            f2.write(
                f"wav. correl.:  slope = {sp['wv_slope']:11.3e}"
                f"  intercept = {sp['wv_intercept']:7.3f}"
                f"  corr. coeff. = {sp['wv_r']:7.3f}\n"
            )


def write_abfind_output(state, result, f1path: str, f2path: str) -> None:
    """Write abfind results to standard and summary output files."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog abfind\n")
            f1.write(f"model:    {state.fmodel}\n")
            f1.write(f"linelist: {state.flines}\n")
            f1.write(f"[M/H] = {state.abscale:+.2f}\n")
        if f2:
            _write_abund_block(f2, state, result['lines'], result['species'])
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_blends_output(state, result, f1path: str, f2path: str) -> None:
    """Write blends results — same abundance table format as abfind."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog blends\n")
            f1.write(f"model:    {state.fmodel}\n")
            f1.write(f"linelist: {state.flines}\n")
            f1.write(f"[M/H] = {state.abscale:+.2f}\n")
        if f2:
            _write_abund_block(f2, state, result['lines'], result['species'],
                               match_all_ions=True)
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_ewfind_output(state, result, f1path: str, f2path: str) -> None:
    """Write ewfind results — predicted EW table."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog ewfind\n")
            f1.write(f"model:    {state.fmodel}\n")
            f1.write(f"linelist: {state.flines}\n")
            f1.write(f"[M/H] = {state.abscale:+.2f}\n")
        if f2:
            f2.write(f"{(state.moditle or '')[:80]:<80}\n")
            f2.write(f"{(state.linitle or '')[:80]:<80}\n")
            f2.write(
                "wavelength         ID      EP   logGF"
                "   EW(pred)   logRW     abund\n"
            )
            for ln in result['lines']:
                ew = ln['ew_pred']
                logrw = math.log10(ew / ln['wave'] / 1000.0) if ew > 0.0 else 0.0
                f2.write(
                    f"{ln['wave']:10.3f}{ln['species']:11.5f}{ln['ep']:8.3f}"
                    f"{ln['loggf']:8.3f}{ew:9.2f}{logrw:10.3f}"
                    f"{ln['abund']:10.3f}\n"
                )
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_synth_output(state, result, f1path: str, f2path: str) -> None:
    """Write synth spectrum — Smooth.f format 1008: f10.3 f10.5."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog synth\n")
            f1.write(f"model:    {state.fmodel}\n")
            f1.write(f"linelist: {state.flines}\n")
            f1.write(f"[M/H] = {state.abscale:+.2f}\n")
        if f2:
            wave = result['wave']
            flux = result['flux_smooth']
            for w, f in zip(wave, flux):
                f2.write(f"{float(w):10.3f}{float(f):10.5f}\n")
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_cog_output(state, result, f1path: str, f2path: str) -> None:
    """Write COG results — one block per line."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog cog\n")
            f1.write(f"model:    {state.fmodel}\n")
            f1.write(f"linelist: {state.flines}\n")
            f1.write(f"[M/H] = {state.abscale:+.2f}\n")
        if f2:
            f2.write(f"{(state.moditle or '')[:80]:<80}\n")
            f2.write(f"{(state.linitle or '')[:80]:<80}\n")
            for ln in result['lines']:
                f2.write(
                    f"# {ln['wave']:.3f}  {ln['species']:.5f}"
                    f"  EP={ln['ep']:.3f}  abund={ln['abund']:.3f}\n"
                )
                for loggf, logrw in zip(ln['loggf'], ln['logrw']):
                    f2.write(f"{loggf:10.3f}{logrw:10.3f}\n")
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_weedout_output(state, result, f1path: str, f2path: str) -> None:
    """Write weedout results — kept lines to f1, discarded to f2."""
    hdr = (
        "wavelength         ID      EP   logGF  dampnum        d0"
        "    EW(obs)  logstrength    ratio\n"
    )

    def _fmt_line(ln):
        logstr = f"{ln['logstrength']:10.3f}" if ln['logstrength'] is not None else "          "
        return (
            f"{ln['wave']:10.3f}{ln['species']:11.5f}{ln['ep']:8.3f}"
            f"{ln['loggf']:8.3f}{ln['dampnum']:9.3f}{ln['d0']:10.3f}"
            f"{ln['ew_obs']:10.3f}{logstr}{ln['ratio']:10.3f}\n"
        )

    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog weedout   xratio = {result['xratio']:.4f}\n")
            f1.write(f"KEPT ({len(result['kept'])} lines)\n")
            f1.write(hdr)
            for ln in result['kept']:
                f1.write(_fmt_line(ln))
        if f2:
            f2.write(f"DISCARDED ({len(result['discarded'])} lines)\n")
            f2.write(hdr)
            for ln in result['discarded']:
                f2.write(_fmt_line(ln))
    finally:
        if f1: f1.close()
        if f2: f2.close()


def write_doflux_output(state, result, f1path: str, f2path: str) -> None:
    """Write doflux results. Doflux.f format 1001: 1p2d12.4, 0p2f10.4."""
    f1 = _open_out(f1path)
    f2 = _open_out(f2path)
    try:
        if f1:
            f1.write(f"pymoog doflux\n")
            f1.write(f"model: {state.fmodel}\n")
        if f2:
            for pt in result['flux']:
                f2.write(
                    f"{pt['wave']:12.4e}{pt['flux']:12.4e}"
                    f"{pt['waveinv']:10.4f}{pt['fluxlog']:10.4f}\n"
                )
    finally:
        if f1: f1.close()
        if f2: f2.close()
