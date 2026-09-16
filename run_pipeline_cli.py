#!/usr/bin/env python
"""Template CLI driver: full SPAE.xspect_ew EW-measurement + SPAE.spae
stellar-parameter MCMC pipeline, for one FITS spectrum.

Usage
-----
    python run_pipeline_cli.py /path/to/directory/containing/one/fits/file

Expects exactly one *.fits file in that directory. Writes, all inside
that same directory:

    line_plots/                     -- one QC plot per measured line
    linelist_with_ew.txt            -- main MOOG-format EW linelist
    linelist_with_ew_flagged.txt    -- untrustworthy lines, for manual review
    ew_measurement.log              -- Stage 1 diagnostics (per-line fit
                                        prints, slope diagnostics, flags)
    mcmc_fit.log                    -- Stage 2 diagnostics + final summary
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

CUSTOMIZE the CONFIG block below before running: it points at the bundled
solar sample data by default, which is almost certainly NOT what you want
for a real target -- in particular, REFERENCE_ATLAS_PATH is only valid
for a star with an actual matching high-S/N reference spectrum (so far:
the Sun); leave it as None for anything else.

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

# MOOG-format linelist of lines to measure (wave, species, EP, loggf, ...)
LINELIST_PATH = os.path.join(
    os.path.dirname(__file__),
    'SPAE/xspect_ew/data/Line_list_sample/Sun_fe_sample.txt')

# Reference atlas for cross-checking a line's wing against an independent
# high-S/N spectrum (see SPAE/xspect_ew/reference_atlas.py) -- only valid
# for a target this actually matches. Set to None for anything but the Sun.
REFERENCE_ATLAS_PATH = os.path.join(
    os.path.dirname(__file__), 'data/kurucz_solar_atlas/fluxspliced.2005')
RESOLVING_POWER = None  # None: estimate empirically from the linelist

WINDOW_SIZE = 1.5
REVIEW = False  # True: pause after measurement for interactive human review
                # of flagged lines (run this script from a real terminal,
                # not a notebook, if you turn this on -- it blocks on input())

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


def find_fits_file(directory):
    matches = sorted(glob.glob(os.path.join(directory, '*.fits')))
    if len(matches) == 0:
        raise FileNotFoundError(f'No .fits file found in {directory}')
    if len(matches) > 1:
        raise ValueError(
            f'Expected exactly one .fits file in {directory}, found {len(matches)}: '
            f'{matches} -- point this script at a directory with just the one '
            f'spectrum you want measured.')
    return matches[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('directory', help='Directory containing exactly one .fits spectrum file')
    args = parser.parse_args()

    directory = os.path.abspath(args.directory)
    spectrum_path = find_fits_file(directory)
    print(f'Found spectrum: {spectrum_path}')

    stdout_orig, stderr_orig = sys.stdout, sys.stderr

    # ---------------- Stage 1: EW measurement ----------------
    from SPAE.pipeline import measure_star_ew

    ew_log_path = os.path.join(directory, 'ew_measurement.log')
    with open(ew_log_path, 'w') as ew_log:
        sys.stdout = Tee(stdout_orig, ew_log)
        sys.stderr = Tee(stderr_orig, ew_log)
        try:
            print(f'=== EW measurement: {spectrum_path} ===')
            t0 = time.time()
            ew_path, flagged_path, spec = measure_star_ew(
                spectrum_path, LINELIST_PATH, directory,
                reference_atlas_path=REFERENCE_ATLAS_PATH,
                resolving_power=RESOLVING_POWER,
                window_size=WINDOW_SIZE, review=REVIEW, save_plots=True)
            print(f'EW measurement time: {time.time()-t0:.2f} s')
            print(f'main linelist:    {ew_path}')
            print(f'flagged linelist: {flagged_path}')
            print(f'line plots:       {os.path.join(directory, "line_plots")}')
        finally:
            sys.stdout, sys.stderr = stdout_orig, stderr_orig
    print(f'EW measurement log: {ew_log_path}')

    # ---------------- Stage 2: MCMC stellar-parameter fit ----------------
    from SPAE.pipeline import fit_stellar_params
    from SPAE.moog.inlines import parse_linelist
    from SPAE.analysis import analyze_run, summarize, plot_trace, plot_corner

    mcmc_log_path = os.path.join(directory, 'mcmc_fit.log')
    with open(mcmc_log_path, 'w') as mcmc_log:
        sys.stdout = Tee(stdout_orig, mcmc_log)
        sys.stderr = Tee(stderr_orig, mcmc_log)
        try:
            print(f'=== MCMC fit: {N_WALKERS} walkers x {N_STEPS} steps ===')
            t0 = time.time()
            #summarize_result=False: save the (expensive) raw sampler to
            #disk BEFORE attempting analysis, so a post-processing problem
            #(e.g. a burn_in that doesn't fit N_STEPS) can never cost
            #re-running the sampler itself -- see BURN_IN's comment above
            sampler, flat_blob, log, _, _ = fit_stellar_params(
                ew_path, x_0=X_0, n_walkers=N_WALKERS, n_steps=N_STEPS,
                n_cores=N_CORES, summarize_result=False)
            print(f'MCMC wall time: {time.time()-t0:.1f} s')
            print(log)

            sampler_path = os.path.join(directory, 'sampler.pickle')
            with open(sampler_path, 'wb') as f:
                pickle.dump(sampler, f)
            np.save(os.path.join(directory, 'flat_blob.npy'), flat_blob)
            print(f'saved {sampler_path}')

            parsed = parse_linelist(ew_path)
            result = analyze_run(sampler, parsed, burn_in=BURN_IN)
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
        finally:
            sys.stdout, sys.stderr = stdout_orig, stderr_orig
    print(f'MCMC fit log: {mcmc_log_path}')

    print('Done.')


if __name__ == '__main__':
    main()
