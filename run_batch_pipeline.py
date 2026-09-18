#!/usr/bin/env python
"""Batch driver: run run_pipeline_cli.py end-to-end, in series, over many
star directories.

Usage
-----
    python run_batch_pipeline.py /path/to/parent_dir [/path/to/solar]

`parent_dir` is a directory whose immediate subdirectories are each one
star (e.g. parent_dir/HD_10383/, parent_dir/HD_102071/, ...). `solar` is
optional and is passed straight through to run_pipeline_cli.py for every
star (see that script's docstring for what it accepts) -- omit it to use
run_pipeline_cli.py's own bundled default for every star.

Picking the target FITS file(s) in each star directory
--------------------------------------------------------
Each star directory may hold any number of .fits files, not necessarily
all the target star: some may be OTHER stars' spectra that ended up in
the same folder (as opposed to several files that ARE the target star,
e.g. different exposures/epochs or bands/orders -- run_pipeline_cli.py's
target argument now accepts any number of files and measures+merges them
exactly like it already does for its `--solar` reference, see that
script's docstring).

So for each star directory:
  - 0 .fits files: skipped, logged as `skip_no_fits`.
  - 1 or more .fits files: every file's header is checked for a
    star-name keyword (OBJECT, TARGNAME, TARGET, STARNAME, OBJNAME --
    checked in that order, skipping placeholder-like values such as
    OBJECT='object', since e.g. Keck/MAKEE solar-sample files have
    exactly that combined with TARGNAME='Ganymede [ L1]', the actual
    observed target) and compared, alphanumeric-only and
    case-insensitive, against the directory's own name. EVERY matching
    file is passed to run_pipeline_cli.py together (one match: the usual
    single-file case; several matches: all of them, merged); zero
    matches skips the directory as `skip_no_match` for manual review
    (logged with every candidate's header name -- or lack of one -- so
    the review is quick).

Output
------
Writes one line per star directory to a summary CSV (default:
batch_summary.csv in the current directory; override with --summary),
with columns: star_dir, fits_file, status, elapsed_s, note. `status` is
one of: ok, failed, skip_no_fits, skip_no_match. `fits_file` is a
semicolon-separated list when more than one file was used.

A star whose run_pipeline_cli.py invocation fails (non-zero exit, e.g. an
unhandled exception partway through EW measurement or the MCMC fit) is
logged as `failed` and the batch continues with the next star -- it does
NOT abort the whole batch. Each star's full subprocess stdout/stderr
(everything run_pipeline_cli.py printed, including anything before its
own per-stage log files were opened) is saved to
<star_dir>/batch_pipeline_stdout.log for post-mortem review; on success
this is largely a duplicate of run_pipeline_cli.py's own per-stage logs
(ew_measurement.log, mcmc_fit.log, etc.) written inside the same
directory.

Each star runs in its own subprocess (rather than importing and calling
run_pipeline_cli.py's main() in-process) so that one star's crash --
including inside SPAE.spae.run_spae()'s multiprocessing pool -- can never
take down the batch or leave stray worker processes behind for the next
star.
"""
import argparse
import csv
import os
import subprocess
import sys
import time

from astropy.io import fits

# Star-name header keywords to check, in priority order. OBJECT is often
# a generic instrument placeholder rather than the actual target name
# (e.g. Keck/MAKEE solar-sample files: OBJECT='object', TARGNAME=
# 'Ganymede [ L1]' -- reflected sunlight off Ganymede, the actual target),
# so every plausible alternative is checked before giving up.
HEADER_NAME_KEYS = ('OBJECT', 'TARGNAME', 'TARGET', 'STARNAME', 'OBJNAME')
PLACEHOLDER_VALUES = {'object', 'unknown', 'none', 'target', 'star', ''}

RUN_PIPELINE_CLI = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'run_pipeline_cli.py')


def _normalize(name):
    """Case- and punctuation-insensitive form for comparing a header
    star name against a directory name (e.g. 'HD 10383', 'hd_10383', and
    'HD10383' should all compare equal)."""
    return ''.join(ch for ch in name.upper() if ch.isalnum())


def _header_star_name(fits_path):
    """Best-effort star name from a FITS primary header, or None if no
    usable (non-placeholder) keyword is present."""
    try:
        with fits.open(fits_path) as hdul:
            header = hdul[0].header
    except Exception:
        return None
    for key in HEADER_NAME_KEYS:
        if key in header:
            value = str(header[key]).strip()
            if value and value.lower() not in PLACEHOLDER_VALUES:
                return value
    return None


def pick_target_fits(star_dir):
    """One star directory -> (list of fits_path, status, note).

    status is one of: 'ok', 'skip_no_fits', 'skip_no_match'. See the
    module docstring for the selection rules -- every .fits file present
    (even just one) is header-checked against the directory name, and
    EVERY matching file is returned (run_pipeline_cli.py merges however
    many target files it's given, the same way it already merges
    multiple solar files).
    """
    fits_files = sorted(
        f for f in os.listdir(star_dir) if f.lower().endswith('.fits'))
    dir_name = _normalize(os.path.basename(star_dir))

    if not fits_files:
        return [], 'skip_no_fits', 'no .fits files in directory'

    candidates = {f: _header_star_name(os.path.join(star_dir, f)) for f in fits_files}
    matches = [f for f, name in candidates.items()
               if name and (dir_name in _normalize(name) or _normalize(name) in dir_name)]

    if matches:
        others = [f for f in fits_files if f not in matches]
        note = (f'picked {len(matches)}/{len(fits_files)} file(s) matching directory: '
                 f'{[(f, candidates[f]) for f in matches]}')
        if others:
            note += f' -- ignored: {others}'
        return [os.path.join(star_dir, f) for f in matches], 'ok', note

    note = (f'{len(fits_files)} fits file(s) present, 0 header-name matches to directory '
             f'{os.path.basename(star_dir)!r}: '
             + ', '.join(f'{f}={candidates[f]!r}' for f in fits_files))
    return [], 'skip_no_match', note


def run_one_star(star_dir, solar_arg, dry_run=False):
    """Pick the target FITS file(s) for star_dir and, if any matched, run
    run_pipeline_cli.py on all of them (as one merged target) via a
    single subprocess -- or, if dry_run, just report which file(s) WOULD
    be picked without measuring anything or calling run_pipeline_cli.py
    at all (see --dry-run in main()).

    Returns a dict with keys star_dir, fits_file, status, elapsed_s, note
    -- one row of the batch summary (see module docstring)."""
    target_paths, status, note = pick_target_fits(star_dir)
    print(f'    {note}')
    if status != 'ok':
        return dict(star_dir=star_dir, fits_file='', status=status,
                     elapsed_s='', note=note)

    fits_file_field = ';'.join(os.path.basename(p) for p in target_paths)

    if dry_run:
        return dict(star_dir=star_dir, fits_file=fits_file_field,
                     status='dry_run_ok', elapsed_s='',
                     note=note + ' -- dry run, not executed')

    cmd = [sys.executable, RUN_PIPELINE_CLI, *target_paths]
    if solar_arg is not None:
        cmd += ['--solar', solar_arg]

    stdout_log_path = os.path.join(star_dir, 'batch_pipeline_stdout.log')
    t0 = time.time()
    with open(stdout_log_path, 'w') as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0

    if proc.returncode == 0:
        status = 'ok'
        note = f'completed in {elapsed:.1f} s'
        print(f'    OK ({elapsed:.1f} s)')
    else:
        status = 'failed'
        note = (f'run_pipeline_cli.py exited {proc.returncode} after {elapsed:.1f} s '
                 f'-- see {stdout_log_path}')
        print(f'    FAILED (exit {proc.returncode}) -- see {stdout_log_path}')

    return dict(star_dir=star_dir, fits_file=fits_file_field,
                status=status, elapsed_s=f'{elapsed:.1f}', note=note)


def run_batch(star_dirs, solar_arg, summary_path, dry_run=False):
    rows = []
    for i, star_dir in enumerate(star_dirs, 1):
        star_dir = os.path.abspath(star_dir)
        print(f'[{i}/{len(star_dirs)}] {os.path.basename(star_dir)}')
        rows.append(run_one_star(star_dir, solar_arg, dry_run=dry_run))

    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f, fieldnames=['star_dir', 'fits_file', 'status', 'elapsed_s', 'note'])
        writer.writeheader()
        writer.writerows(rows)

    counts = {}
    for row in rows:
        counts[row['status']] = counts.get(row['status'], 0) + 1
    n_ok = counts.get('ok', 0) + counts.get('dry_run_ok', 0)
    verb = 'would run' if dry_run else 'succeeded'
    other_counts = ', '.join(f'{v} {k}' for k, v in counts.items()
                              if k not in ('ok', 'dry_run_ok'))
    print(f'\nBatch summary written to {summary_path}')
    print(f'{n_ok}/{len(rows)} {verb}' + (f' ({other_counts})' if other_counts else ''))
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'root', help='Parent directory whose immediate subdirectories are each one star')
    parser.add_argument(
        'solar_spectrum', nargs='?', default=None,
        help="Solar reference passed straight through to run_pipeline_cli.py for every "
             "star (file or directory; see its own docstring). Omit to use its bundled "
             "default for every star.")
    parser.add_argument(
        '--summary', default='batch_summary.csv',
        help='Path to write the batch summary CSV (default: ./batch_summary.csv)')
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Only run the target-FITS-file selection logic for every star directory and '
             'report what WOULD be picked/skipped -- never calls run_pipeline_cli.py, so no '
             'EW measurement or MCMC actually runs. Useful as a smoke test of the directory '
             'tree and header-matching before committing to a real (slow) batch run.')
    args = parser.parse_args()

    root = os.path.abspath(args.root)
    star_dirs = sorted(
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d)))
    print(f'Found {len(star_dirs)} star director{"y" if len(star_dirs) == 1 else "ies"} '
          f'under {root}')

    run_batch(star_dirs, args.solar_spectrum, os.path.abspath(args.summary),
              dry_run=args.dry_run)


if __name__ == '__main__':
    main()
