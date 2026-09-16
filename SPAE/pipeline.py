"""End-to-end driver: measured spectrum -> equivalent widths -> stellar
parameters, in one of two modes.

- Fully automated (review=False, the default): measure_star_ew() writes
  the final EW linelist unattended, using check_for_flags()'s existing
  data-only criteria to route anything untrustworthy to a separate
  flagged file rather than the main one -- no line ever silently ends up
  in the science linelist just because nothing stopped it.
- Human-in-the-loop (review=True): the same measurement run, but pauses
  after check_for_flags() to show each flagged line's QC plot and its
  reason(s), and lets a reviewer keep it anyway. That decision is
  recorded on Spectrum_Data.lines_human_keep, which survives every later
  check_for_flags() call (including the one make_ew_doc() makes
  internally) -- see spectrum_data.py's check_for_flags()/make_ew_doc()
  for exactly how.

Both modes produce the same MOOG-format linelist that SPAE.spae.run_spae()
and SPAE.analysis.analyze_run() already consume directly -- no format
translation needed between the EW-measurement side of this package and
the stellar-parameter-fitting side.
"""

import os

from .xspect_ew import Spectrum_Data
from .moog.inlines import parse_linelist
from .abunds import sun_abs as _sun_abs
from .spae import run_spae
from .analysis import analyze_run, summarize


def measure_star_ew(spectrum_path, linelist_path, output_dir, reference_atlas_path=None,
                     resolving_power=None, window_size=1.5, review=False,
                     doc_title='STARNAME, PROJECT, YEAR; ', **measure_kwargs):
    """Stage 1: measure every line's EW for one star and write the final
    MOOG-format linelist, either unattended or with a human review pass.

    Parameters
    ----------
    spectrum_path : str -- passed straight to Spectrum_Data(); any format
        readers.py auto-detects
    linelist_path : str -- MOOG-format linelist of lines to measure
    output_dir : str -- created if it doesn't exist; holds
        linelist_with_ew.txt (main output), linelist_with_ew_flagged.txt
        (untrustworthy lines, for review, not automatic use), and --
        only when review=True -- line_plots/ (this function chdir()s
        there for the run so measure_all_ew(save_all=True)'s relative
        line_plots/ output lands inside output_dir rather than wherever
        the caller happened to be)
    reference_atlas_path : str or None -- see Spectrum_Data.
        load_reference_atlas(); only meaningful for a target this
        actually has a matching high-S/N reference for (so far: the Sun)
    review : bool -- False (default): fully automated, no prompts.
        True: after measuring, pause and show each flagged line's plot
        and reason(s) one at a time, and ask whether to keep it anyway
        (see _review_flagged_lines())
    **measure_kwargs : forwarded to measure_all_ew() (e.g. exclude_lines,
        fit_continuum, auto_widen, widen_window_size, slope_sig_thresh)

    Returns
    -------
    ew_path, flagged_path : str -- the two output file paths
    spec : the Spectrum_Data instance, for any further inspection
    """
    os.makedirs(output_dir, exist_ok=True)
    ew_path = os.path.join(output_dir, 'linelist_with_ew.txt')
    flagged_path = os.path.join(output_dir, 'linelist_with_ew_flagged.txt')

    spec = Spectrum_Data(spectrum_path)
    spec.normalize_all()
    spec.apply_rv_shift()
    spec.load_lines(linelist_path)
    if reference_atlas_path is not None:
        spec.load_reference_atlas(reference_atlas_path, resolving_power=resolving_power)

    if review:
        #save_all=True is what produces line_plots/<element>_<wave>_
        #<order>.pdf for every line -- needed so the reviewer has
        #something to look at below, not just the printed reason
        cwd = os.getcwd()
        os.chdir(output_dir)
        try:
            spec.measure_all_ew(window_size=window_size, save_all=True, **measure_kwargs)
        finally:
            os.chdir(cwd)
    else:
        spec.measure_all_ew(window_size=window_size, save_all=False, **measure_kwargs)

    spec.check_for_flags()
    if review:
        _review_flagged_lines(spec, plots_dir=os.path.join(output_dir, 'line_plots'))

    spec.make_ew_doc(ew_path, doc_title=doc_title, flagged_name=flagged_path)
    return ew_path, flagged_path, spec


def _review_flagged_lines(spec, plots_dir):
    """Human-in-the-loop step: one prompt per flagged line. Blocking
    input() by design -- this is meant to be run from an interactive
    terminal or notebook, not unattended (that's review=False's job).

    Sets spec.lines_human_keep[i] = True for every line the reviewer
    chooses to keep; leaves it False (the default) for everything else,
    which is exactly equivalent to not reviewing at all.
    """
    import numpy as np

    flagged_idx = np.where(spec.lines_check_flag)[0]
    if len(flagged_idx) == 0:
        print('No flagged lines -- nothing to review.')
        return

    print(f'\n{len(flagged_idx)} line(s) flagged for review. QC plots are in {plots_dir}/\n')
    for i in flagged_idx:
        print(f'  {spec.lines[i]}: {spec.lines_flag_reasons[i]}')
        resp = input('    Keep this line anyway? [y/N]: ').strip().lower()
        if resp == 'y':
            spec.lines_human_keep[i] = True
            print('    -> kept (will appear in the main linelist, annotated with this reason)')
        else:
            print('    -> excluded (stays in the flagged file)')


def solar_reference(sun_spectrum_path, sun_linelist_path, output_dir, reference_atlas_path=None,
                     resolving_power=None, teff_sun=5777, logg_sun=4.44, feh_sun=0.0,
                     micro_sun=1.38, **measure_kwargs):
    """Convenience: measure EWs for the Sun with the same pipeline used for
    any other star, then derive the solar abundances SPAE.spae.run_spae()'s
    differential mode (sun_el, sun_abs) needs as its reference point.

    Runs unattended (review is not exposed here) -- this is meant to be a
    one-time, low-stakes setup step reusing a well-characterized reference
    spectrum, not a per-star measurement that needs a human looking at it.
    """
    ew_path, flagged_path, _ = measure_star_ew(
        sun_spectrum_path, sun_linelist_path, output_dir,
        reference_atlas_path=reference_atlas_path, resolving_power=resolving_power,
        review=False, **measure_kwargs)
    parsed = parse_linelist(ew_path)
    sun_el, sun_abunds = _sun_abs(parsed, teff_sun=teff_sun, logg_sun=logg_sun,
                                   feh_sun=feh_sun, micro_sun=micro_sun)
    return sun_el, sun_abunds


def fit_stellar_params(ew_linelist_path, sun_el=None, sun_abs=None,
                        x_0=(5777, 4.44, 0.01, 1.38), summarize_result=True,
                        analyze_kwargs=None, **run_spae_kwargs):
    """Stage 2: MCMC stellar-parameter fit from an already-measured EW
    linelist (e.g. measure_star_ew()'s ew_path).

    Thin wrapper around SPAE.spae.run_spae() -- parses the linelist once
    and reuses the same parsed dict for both the MCMC run and
    SPAE.analysis.analyze_run(), rather than parsing it twice.

    Returns
    -------
    sampler, flat_blob, log : run_spae()'s own return values
    result : SPAE.analysis.AnalysisResult, or None if summarize_result=False
    summary : the human-readable summarize(result) string, or None
    """
    parsed = parse_linelist(ew_linelist_path)
    sampler, flat_blob, log = run_spae(parsed, sun_el=sun_el, sun_abs=sun_abs,
                                        x_0=x_0, **run_spae_kwargs)

    result, summary = None, None
    if summarize_result:
        result = analyze_run(sampler, parsed, sun_el=sun_el, sun_abs=sun_abs,
                              **(analyze_kwargs or {}))
        summary = summarize(result)
        print(summary)

    return sampler, flat_blob, log, result, summary


def run_full_pipeline(spectrum_path, linelist_path, output_dir, reference_atlas_path=None,
                       resolving_power=None, window_size=1.5, review=False,
                       sun_el=None, sun_abs=None, x_0=(5777, 4.44, 0.01, 1.38),
                       measure_kwargs=None, run_spae_kwargs=None):
    """One-call convenience: measure_star_ew() -> fit_stellar_params().

    review=False end to end is the fully-automated path; review=True
    pauses once, after EW measurement, for the human review step -- the
    MCMC stage itself has no human-in-the-loop equivalent (a stellar-
    parameter posterior isn't something to eyeball line by line the way
    an individual EW fit is).
    """
    ew_path, flagged_path, spec = measure_star_ew(
        spectrum_path, linelist_path, output_dir,
        reference_atlas_path=reference_atlas_path, resolving_power=resolving_power,
        window_size=window_size, review=review, **(measure_kwargs or {}))

    sampler, flat_blob, log, result, summary = fit_stellar_params(
        ew_path, sun_el=sun_el, sun_abs=sun_abs, x_0=x_0, **(run_spae_kwargs or {}))

    return {
        'ew_path': ew_path, 'flagged_path': flagged_path, 'spec': spec,
        'sampler': sampler, 'flat_blob': flat_blob, 'log': log,
        'result': result, 'summary': summary,
    }
