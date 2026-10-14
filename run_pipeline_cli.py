#!/usr/bin/env python
"""General CLI driver: full SPAE.xspect_ew EW-measurement + SPAE.spae
stellar-parameter MCMC pipeline, for one target star, using DIFFERENTIAL
(relative-to-solar) abundances.

Usage
-----
    python run_pipeline_cli.py /path/to/target.fits [target2.fits ...] [--solar /path/to/solar]

One or more `target.fits` files are given -- all of them the SAME star
(e.g. several exposures/epochs, or several bands/orders that together
cover more of the linelist than any one alone does; see "Why multiple
files at once" below). `--solar` is optional and may be:
  - omitted entirely: defaults to the bundled solar sample directory,
    SPAE/xspect_ew/data/spectra_sample/Solar/ (three FITS files covering
    different, mostly non-overlapping wavelength ranges/orders -- see
    below for why all three get used together)
  - a single FITS file: one solar spectrum
  - a directory: every *.fits file directly inside it is used

Why multiple files at once (both target and solar): a single echelle
exposure only covers part of the optical range, so a linelist spanning
the whole range (like the bundled Sun_fe_sample.txt) needs several
exposures to have every line measured somewhere. Each file is measured
independently (same linelist, same routine as any other star), then the
individual EW linelists are merged into one combined linelist by
wavelength -- safe because each file's own `measure_star_ew()` output
already only contains lines that fell inside ITS wavelength coverage and
were successfully measured, so the pieces contribute disjoint (or, at
worst, redundant-but-consistent) sets of lines, never conflicting values
for the same line. (This assumes every file given really is the same
star -- picking that out from a directory that might also hold OTHER
stars' spectra is the calling script's job, e.g. run_batch_pipeline.py's
header-name matching, not this script's.)

Only lines measured in BOTH the target star and the (combined) solar
reference are used for the MCMC fit: SPAE.abunds.rel_abunds() matches
the two per-line abundance arrays by ARRAY POSITION, not by searching
for the same wavelength, so it silently misaligns (and, for a species
with few lines like Fe II, can come out entirely as NaN) the moment the
two linelists disagree on which lines exist -- MOOG's abfind can
legitimately return a different subset of lines at different (Teff,
logg, [Fe/H], micro), so this isn't a hypothetical edge case. Explicitly
filtering the target's AND the solar reference's EW linelists down to
their common wavelengths (same order, both files) before the MCMC ever
starts removes the most common source of that misalignment (though a
single MCMC trial can still transiently drop a line MOOG itself finds
unphysical at that specific trial point -- a separate, deeper issue in
rel_abunds() itself, not fixed by this filtering).

Writes, inside the TARGET spectrum's own directory (if more than one
target file is given, all must live in the SAME directory -- that
directory is where every output below lands):

    line_plots/                     -- one QC plot per measured line (only
                                        target file case; see below for
                                        the multi-file layout)
    linelist_with_ew.txt            -- main MOOG-format EW linelist (ALL
                                        lines measured in the target,
                                        unfiltered; the by-wavelength
                                        merge of every target file's own
                                        linelist_with_ew.txt if more than
                                        one was given)
    linelist_with_ew_common.txt     -- the subset of the above also
                                        measured in the solar reference --
                                        what the MCMC fit actually uses
    linelist_with_ew_flagged.txt    -- untrustworthy lines, for manual
                                        review (only target file case;
                                        with multiple target files, each
                                        one's own flagged file stays in
                                        its own subdirectory, unmerged --
                                        see below)
    ew_measurement.log              -- Stage 2 diagnostics (per-line fit
                                        prints, slope diagnostics, flags)
    mcmc_fit.log                    -- Stage 4 diagnostics + final summary
    sampler.pickle                  -- the raw emcee sampler, reloadable via
                                        SPAE.analysis.load_sampler() for a
                                        later re-analysis (e.g. a different
                                        burn_in after inspecting a trace plot)
    flat_blob.npy                   -- convenience cache of the unclipped
                                        per-walker abundance blobs
    trace_plot.png                  -- chain value vs. step, FULL unclipped
                                        chain -- inspect this to pick a real
                                        BURN_IN, then re-analyze (see
                                        summarize()'s own printed recipe)
    corner_plot.png                 -- posterior corner plot of the
                                        (BURN_IN-clipped) chain; skipped
                                        with a note if the optional `corner`
                                        package isn't installed
    ep_vs_abundance.png             -- Fe I/Fe II per-line differential
                                        abundance vs. excitation potential,
                                        at the MCMC median solution, with
                                        the same Fe I linear fit/slope
                                        analyze_run()'s EP-balance check
                                        uses -- only for the common-lines
                                        subset actually fit (see Stage 3
                                        below)

If more than one target file was given, each one is measured into its
own subdirectory of the target directory first (named after that file's
own basename, minus extension), exactly like the solar files below --
each with its own line_plots/, linelist_with_ew.txt,
linelist_with_ew_flagged.txt -- and only the merged linelist_with_ew.txt
described above is written at the top level.

...and, inside a solar_reference/ subdirectory of that same directory
(kept separate from the target's own outputs, and never written into any
solar FITS file's own directory, which may be shared/reused across many
runs -- e.g. the bundled package data):

    <solar file's basename>/        -- one subdirectory per solar FITS
                                        file measured, each with its own
                                        line_plots/, linelist_with_ew.txt,
                                        linelist_with_ew_flagged.txt
    linelist_with_ew_merged.txt     -- all solar files' linelist_with_ew.txt
                                        combined by wavelength
    linelist_with_ew_common.txt     -- the subset also measured in the
                                        target -- what sun_el/sun_abs
                                        actually get derived from
    solar_measurement.log           -- Stage 1 diagnostics, for every
                                        solar file

CUSTOMIZE the CONFIG block below before running -- in particular
LINELIST_PATH (must be a linelist BOTH the target and the Sun will
plausibly show all these lines for) and REFERENCE_ATLAS_PATH (a
per-line wing cross-check against an independent high-S/N reference --
only meaningful for a target that actually matches one; left at None by
default since it's essentially never valid except for the Sun itself,
handled separately below via SOLAR_REFERENCE_ATLAS_PATH).

IMPORTANT: this script uses multiprocessing (via SPAE.spae.run_spae()),
so the `if __name__ == '__main__':` guard below is REQUIRED on macOS/
Windows (the default multiprocessing start method there re-imports this
file in every worker process) -- not optional boilerplate. If you copy
logic out of this file into your own script, keep that guard.
"""
import argparse
import glob
import os
import pickle
import sys
import time

import matplotlib
matplotlib.use('Agg')  # no interactive windows -- this is meant to run unattended
import numpy as np


# ============================== CONFIG ==============================
# Edit these for your target before running.

# Repository root, resolved from this file's own location -- NOT the
# process's current working directory, so every default path below stays
# correct no matter where this script is invoked from (e.g. a relative
# `python run_pipeline_cli.py ...` from inside a star's own data
# directory, or a subprocess launched by run_batch_pipeline.py without an
# explicit `cwd=`). Plain os.path.dirname(__file__) is NOT enough for
# this: if the script was invoked with a relative path (the common case),
# __file__ itself is relative, and joining it with anything still
# resolves against the CWD, not the repo -- os.path.abspath() first is
# what actually pins it down.
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

# MOOG-format linelist of lines to measure (wave, species, EP, loggf, ...)
# -- used for BOTH the target and every solar measurement (see module
# docstring for why they must match)
LINELIST_PATH = os.path.join(
    REPO_ROOT, 'SPAE/xspect_ew/data/Line_list_sample/Sun_fe_sample.txt')

# Default solar reference when no solar argument is given on the command
# line -- see module docstring for why all three bundled files get used
# together.
DEFAULT_SOLAR_DIR = os.path.join(
    REPO_ROOT, 'SPAE/xspect_ew/data/spectra_sample/Solar')

# Reference atlas for cross-checking a line's wing against an independent
# high-S/N spectrum (see SPAE/xspect_ew/reference_atlas.py) -- only valid
# for a target this actually matches. None (default) is correct for
# almost any real target; every solar measurement below uses
# SOLAR_REFERENCE_ATLAS_PATH instead, since it IS a valid match there.
REFERENCE_ATLAS_PATH = None
SOLAR_REFERENCE_ATLAS_PATH = os.path.join(
    REPO_ROOT, 'data/kurucz_solar_atlas/fluxspliced.2005')
RESOLVING_POWER = None  # None: estimate empirically from the linelist

WINDOW_SIZE = 1.5
REVIEW = False  # True: pause after measurement for interactive human review
                # of flagged lines (run this script from a real terminal,
                # not a notebook, if you turn this on -- it blocks on input());
                # applies to the TARGET measurement only -- every solar
                # measurement always runs unattended
SAVE_PLOTS = True  # save every line's QC plot to line_plots/ -- confirmed
                    # to cost ~0.4 s/line (matplotlib figure + vector PDF
                    # savefig), dwarfing the actual fit itself (a few ms/
                    # line); set False to skip plotting entirely for a much
                    # faster run when you don't need the full audit trail

# Solar parameters assumed when deriving per-line solar abundances from
# the measured (merged, common-lines-filtered) solar EWs (see
# SPAE.abunds.sun_abs) -- standard literature values, not fitted
TEFF_SUN = 5777
LOGG_SUN = 4.44
FEH_SUN = 0.0
MICRO_SUN = 1.38

# MCMC (see SPAE.spae.run_spae())
X_0 = (5777, 4.44, 0.01, 1.38)   # (teff, logg, feh, micro) initial guess
N_WALKERS = 40
N_STEPS = 1000
N_CORES = None                   # None: use every available core

# Discarded from the start of the chain during analysis -- must be well
# below N_STEPS (analyze_run()'s own default, 100, assumes something close
# to N_STEPS=1000; e.g. for a short N_STEPS=100 test run, drop this to ~20,
# or analysis will silently see an empty post-burn-in chain and crash).
BURN_IN = 100
# =====================================================================


class Tee:
    """Duplicate writes to every stream given (e.g. the real terminal AND
    a log file), so a long run's output is both visible live and saved."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def _resolve_solar_paths(solar_arg):
    """CLI argument (None, a file, or a directory) -> sorted list of solar
    FITS file paths. See module docstring for the default/file/directory
    behavior."""
    if solar_arg is None:
        solar_arg = DEFAULT_SOLAR_DIR
    solar_arg = os.path.abspath(solar_arg)
    if os.path.isdir(solar_arg):
        paths = sorted(glob.glob(os.path.join(solar_arg, '*.fits')))
        if not paths:
            raise FileNotFoundError(f'No .fits files found in {solar_arg}')
        return paths
    return [solar_arg]


def _read_ew_lines(path):
    """MOOG-format EW file -> {wavelength (str, as written): full line
    (str, incl. newline)}, skipping the header line. Keyed by the
    wavelength field exactly as make_ew_doc() wrote it (both files derive
    it from the same underlying linelist entries, so string equality is
    exact and avoids any float round-trip concern) -- lets the merge/
    intersection logic below work on plain text without needing to
    understand parse_linelist()'s internal (MOOG-Fortran-derived) dict
    format at all.
    """
    lines = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            if not line.strip():
                continue
            wave = line.split()[0]
            lines[wave] = line
    return lines


def _write_ew_lines(path, header, wave_to_line, waves):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(header)
        for wave in sorted(waves, key=float):
            f.write(wave_to_line[wave])


def measure_and_merge(paths, output_dir, reference_atlas_path, review, doc_title_fn,
                       merged_filename, label):
    """Measure EWs for every one of `paths` (each in its own subdirectory
    of output_dir, via the same measure_star_ew() any single-file target
    or solar spectrum uses), then merge them by wavelength into one
    combined EW linelist. See module docstring for why this merge is
    safe -- ONLY when every path really is the same star (this function
    has no way to check that itself).

    doc_title_fn : str (piece basename, no extension) -> str, the
        doc_title passed to measure_star_ew() for that piece.
    label : str -- used only in progress prints (e.g. 'solar file',
        'target file').

    Returns
    -------
    merged_path : str -- output_dir/<merged_filename>
    """
    from SPAE.pipeline import measure_star_ew

    merged = {}
    header = None
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        piece_dir = os.path.join(output_dir, name)
        print(f'--- {label}: {path} ---')
        ew_path, flagged_path, _ = measure_star_ew(
            path, LINELIST_PATH, piece_dir,
            reference_atlas_path=reference_atlas_path,
            resolving_power=RESOLVING_POWER, window_size=WINDOW_SIZE,
            review=review, save_plots=SAVE_PLOTS,
            doc_title=doc_title_fn(name))
        piece_lines = _read_ew_lines(ew_path)
        print(f'    {len(piece_lines)} line(s) measured')
        if header is None:
            with open(ew_path) as f:
                header = f.readline()
        overlap = set(piece_lines) & set(merged)
        if overlap:
            #expected to be rare (pieces are mostly disjoint wavelength
            #ranges), and harmless when it happens -- just keep whichever
            #copy was merged first, both are the same star's EW for the
            #same line measured independently, not a real conflict to
            #resolve carefully
            print(f'    {len(overlap)} line(s) already covered by an earlier '
                  f'{label} -- keeping the first measurement')
        merged.update({k: v for k, v in piece_lines.items() if k not in merged})

    merged_path = os.path.join(output_dir, merged_filename)
    _write_ew_lines(merged_path, header, merged, merged.keys())
    print(f'merged {label} linelist: {merged_path} ({len(merged)} line(s) total)')
    return merged_path


def measure_solar_reference(solar_paths, output_dir):
    """Thin measure_and_merge() wrapper for the solar reference -- see
    that function and the module docstring for details.

    Returns
    -------
    merged_path : str -- output_dir/linelist_with_ew_merged.txt
    """
    return measure_and_merge(
        solar_paths, output_dir, reference_atlas_path=SOLAR_REFERENCE_ATLAS_PATH,
        review=False, doc_title_fn=lambda name: f'Solar reference ({name}); ',
        merged_filename='linelist_with_ew_merged.txt', label='solar file')


def measure_target(target_paths, output_dir):
    """Measure EWs for the target star, from one or more FITS files that
    are all the SAME star -- see measure_and_merge() and the module
    docstring for the multi-file merge behavior.

    Returns
    -------
    ew_path, flagged_path : str -- flagged_path is None when
        target_paths has more than one entry (each piece's own flagged
        file stays in its own subdirectory, unmerged; see module
        docstring)
    """
    if len(target_paths) == 1:
        from SPAE.pipeline import measure_star_ew
        ew_path, flagged_path, _ = measure_star_ew(
            target_paths[0], LINELIST_PATH, output_dir,
            reference_atlas_path=REFERENCE_ATLAS_PATH, resolving_power=RESOLVING_POWER,
            window_size=WINDOW_SIZE, review=REVIEW, save_plots=SAVE_PLOTS)
        return ew_path, flagged_path

    ew_path = measure_and_merge(
        target_paths, output_dir, reference_atlas_path=REFERENCE_ATLAS_PATH,
        review=REVIEW, doc_title_fn=lambda name: f'{name}; ',
        merged_filename='linelist_with_ew.txt', label='target file')
    return ew_path, None


def restrict_to_common_lines(target_ew_path, solar_ew_path, target_common_path,
                              solar_common_path):
    """Filter both EW linelists down to the wavelengths present in BOTH
    (same sorted order in both output files) -- see module docstring for
    why this matters for rel_abunds()'s positional matching."""
    target_lines = _read_ew_lines(target_ew_path)
    solar_lines = _read_ew_lines(solar_ew_path)
    common = set(target_lines) & set(solar_lines)
    print(f'{len(common)} line(s) measured in both the target ({len(target_lines)} total) '
          f'and the solar reference ({len(solar_lines)} total)')

    with open(target_ew_path) as f:
        target_header = f.readline()
    with open(solar_ew_path) as f:
        solar_header = f.readline()

    _write_ew_lines(target_common_path, target_header, target_lines, common)
    _write_ew_lines(solar_common_path, solar_header, solar_lines, common)
    return target_common_path, solar_common_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('spectrum', nargs='+',
                         help='Target star FITS spectrum file(s) -- all the SAME star '
                              '(e.g. several exposures/bands); see module docstring for the '
                              'multi-file merge behavior. All must live in the same directory.')
    parser.add_argument('--solar', dest='solar_spectrum', default=None,
                         help='Observed solar FITS spectrum file or directory of FITS files, '
                              'for the differential abundance reference. Defaults to the '
                              'bundled SPAE/xspect_ew/data/spectra_sample/Solar/ (all three '
                              'files) if omitted.')
    args = parser.parse_args()

    target_paths = [os.path.abspath(p) for p in args.spectrum]
    directories = {os.path.dirname(p) for p in target_paths}
    if len(directories) > 1:
        parser.error(f'all target spectrum files must live in the same directory, got: '
                      f'{sorted(directories)}')
    directory = directories.pop()
    solar_paths = _resolve_solar_paths(args.solar_spectrum)
    solar_directory = os.path.join(directory, 'solar_reference')
    os.makedirs(directory, exist_ok=True)
    print(f'Target spectra ({len(target_paths)}):')
    for p in target_paths:
        print(f'  {p}')
    print(f'Solar spectra ({len(solar_paths)}):')
    for p in solar_paths:
        print(f'  {p}')

    stdout_orig, stderr_orig = sys.stdout, sys.stderr

    # ---------------- Stage 1: solar EW measurement (all files, merged) ----------------
    solar_log_path = os.path.join(directory, 'solar_measurement.log')
    with open(solar_log_path, 'w') as solar_log:
        sys.stdout = Tee(stdout_orig, solar_log)
        sys.stderr = Tee(stderr_orig, solar_log)
        try:
            print(f'=== Solar EW measurement (reference for differential abundances) ===')
            t0 = time.time()
            solar_merged_path = measure_solar_reference(solar_paths, solar_directory)
            print(f'Solar EW measurement time: {time.time()-t0:.2f} s')
        finally:
            sys.stdout, sys.stderr = stdout_orig, stderr_orig
    print(f'Solar measurement log: {solar_log_path}')

    # ---------------- Stage 2: target EW measurement ----------------
    ew_log_path = os.path.join(directory, 'ew_measurement.log')
    with open(ew_log_path, 'w') as ew_log:
        sys.stdout = Tee(stdout_orig, ew_log)
        sys.stderr = Tee(stderr_orig, ew_log)
        try:
            print(f'=== EW measurement: {target_paths} ===')
            t0 = time.time()
            ew_path, flagged_path = measure_target(target_paths, directory)
            print(f'EW measurement time: {time.time()-t0:.2f} s')
            print(f'main linelist:    {ew_path}')
            print(f'flagged linelist: {flagged_path}' if flagged_path is not None else
                  'flagged linelist: none merged -- see each target file\'s own subdirectory')
            if SAVE_PLOTS:
                print(f'line plots:       {os.path.join(directory, "line_plots")}')
        finally:
            sys.stdout, sys.stderr = stdout_orig, stderr_orig
    print(f'EW measurement log: {ew_log_path}')

    # ---------------- Stage 3: restrict to lines measured in both ----------------
    target_common_path = os.path.join(directory, 'linelist_with_ew_common.txt')
    solar_common_path = os.path.join(solar_directory, 'linelist_with_ew_common.txt')
    restrict_to_common_lines(ew_path, solar_merged_path, target_common_path, solar_common_path)

    from SPAE.abunds import sun_abs as compute_sun_abs
    from SPAE.moog.inlines import parse_linelist
    sun_el, sun_abs = compute_sun_abs(parse_linelist(solar_common_path), teff_sun=TEFF_SUN,
                                       logg_sun=LOGG_SUN, feh_sun=FEH_SUN, micro_sun=MICRO_SUN)

    # ---------------- Stage 4: MCMC stellar-parameter fit (differential) ----------------
    from SPAE.pipeline import fit_stellar_params
    from SPAE.analysis import analyze_run, summarize, plot_trace, plot_corner

    mcmc_log_path = os.path.join(directory, 'mcmc_fit.log')
    with open(mcmc_log_path, 'w') as mcmc_log:
        sys.stdout = Tee(stdout_orig, mcmc_log)
        sys.stderr = Tee(stderr_orig, mcmc_log)
        try:
            print(f'=== MCMC fit (relative to solar, common lines only): '
                  f'{N_WALKERS} walkers x {N_STEPS} steps ===')
            t0 = time.time()
            #summarize_result=False: save the (expensive) raw sampler to
            #disk BEFORE attempting analysis, so a post-processing problem
            #(e.g. a burn_in that doesn't fit N_STEPS) can never cost
            #re-running the sampler itself -- see BURN_IN's comment above
            sampler, flat_blob, log, _, _ = fit_stellar_params(
                target_common_path, sun_el=sun_el, sun_abs=sun_abs, x_0=X_0,
                n_walkers=N_WALKERS, n_steps=N_STEPS, n_cores=N_CORES,
                summarize_result=False)
            print(f'MCMC wall time: {time.time()-t0:.1f} s')
            print(log)

            sampler_path = os.path.join(directory, 'sampler.pickle')
            with open(sampler_path, 'wb') as f:
                pickle.dump(sampler, f)
            np.save(os.path.join(directory, 'flat_blob.npy'), flat_blob)
            print(f'saved {sampler_path}')

            parsed = parse_linelist(target_common_path)
            result = analyze_run(sampler, parsed, sun_el=sun_el, sun_abs=sun_abs, burn_in=BURN_IN)
            #provisional=True: BURN_IN above is a hardcoded config value,
            #not chosen by inspecting this run's own trace plot -- treat
            #this as a first look, then see summarize()'s own banner for
            #the recipe to get a final result
            summary = summarize(result, provisional=True)
            print('=== SUMMARY ===')
            print(summary)

            #trace plot: FULL unclipped chain, on purpose -- burn_in is a
            #visual judgment call made by looking at this, so pre-clipping
            #it here would defeat the point (see plot_trace()'s docstring)
            trace_fig, _ = plot_trace(sampler.get_chain())
            trace_path = os.path.join(directory, 'trace_plot.png')
            trace_fig.savefig(trace_path, dpi=150)
            print(f'saved {trace_path}')

            #corner plot: the CLIPPED posterior (result.clipped), i.e.
            #after BURN_IN and the acceptance-fraction cut already applied
            try:
                corner_fig, _ = plot_corner(result.clipped)
                corner_path = os.path.join(directory, 'corner_plot.png')
                corner_fig.savefig(corner_path, dpi=150)
                print(f'saved {corner_path}')
            except ImportError as e:
                print(f'skipped corner plot: {e}')

            #excitation-potential vs. Fe abundance -- Fe I and Fe II
            #separately, at the MCMC median solution (result.median_balance
            #is already computed by analyze_run() above, no extra MOOG
            #call needed). This is the "no flat continuum"-style physical
            #check for the WHOLE fit, not a per-line QC plot -- a real
            #trend here means the adopted Teff is off (see the Fe I fit
            #line/slope annotated on the plot itself).
            import matplotlib.pyplot as plt
            balance = result.median_balance
            fe1, fe2 = balance.fe1, balance.fe2
            ep_fig, ep_ax = plt.subplots(figsize=(9, 6))
            #y-error bars: each line's own abundance uncertainty, propagated
            #from its measured EW error (SPAE.moog.lineabund.line_abund_err(),
            #via the abund_err field in abunds._RESULT_DTYPE) -- not the
            #line-to-line scatter, which is a separate, already-plotted trend
            ep_ax.errorbar(fe1['EP'], fe1['abund'], yerr=fe1['abund_err'], fmt='o',
                           color='#377eb8', ecolor='#377eb8', elinewidth=1, capsize=3,
                           markersize=5, label=f'Fe I (n={len(fe1)})')
            ep_ax.errorbar(fe2['EP'], fe2['abund'], yerr=fe2['abund_err'], fmt='^',
                           color='#e41a1c', ecolor='#e41a1c', elinewidth=1, capsize=3,
                           markersize=6, label=f'Fe II (n={len(fe2)})')
            ep_line_x = np.array([fe1['EP'].min(), fe1['EP'].max()])
            ep_ax.plot(ep_line_x, balance.ep.intercept + balance.ep.slope * ep_line_x, 'k--',
                       label=f"Fe I fit: slope={balance.ep.slope:+.4f} "
                             f"(r={balance.ep.rvalue:+.3f}, p={balance.ep.pvalue:.3f})")
            ep_ax.set_xlabel('Excitation Potential (eV)', fontsize=12)
            #per-line differential abundance (this line's abundance minus
            #the Sun's own, from rel_abunds() -- see balance_stats()) --
            #NOT the same number as the bulk MCMC [Fe/H] in the title
            #below, which is the fitted model parameter, not a per-line
            #average of this plot's points
            ep_ax.set_ylabel('[Fe/H] (per line)', fontsize=12)
            ep_ax.set_title(
                f'{os.path.basename(directory)} -- Fe abundance vs. excitation potential '
                f'(MCMC median, common lines only)\nTeff={result.median.teff.median:.0f} K, '
                f'logg={result.median.logg.median:.2f}, [Fe/H]={result.median.feh.median:.2f}, '
                f'micro={result.median.micro.median:.2f} km/s')
            ep_ax.legend(fontsize=10)
            ep_fig.tight_layout()
            ep_path = os.path.join(directory, 'ep_vs_abundance.png')
            ep_fig.savefig(ep_path, dpi=150)
            print(f'saved {ep_path}')
        finally:
            sys.stdout, sys.stderr = stdout_orig, stderr_orig
    print(f'MCMC fit log: {mcmc_log_path}')

    print('Done.')


if __name__ == '__main__':
    main()
