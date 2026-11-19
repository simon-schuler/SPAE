#!/usr/bin/env python
"""Batch driver: run run_pipeline_cli.py end-to-end, in series, over many
star directories.

Usage
-----
    python run_batch_pipeline.py /path/to/parent_dir [/path/to/solar]
    python run_batch_pipeline.py /path/to/parent_dir --smoke-test-only

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

Smoke test and preflight (run automatically before a real batch)
------------------------------------------------------------------
1. Preflight: the required Python packages import, the linelist and the
   solar reference exist, and the parent directory exists, is a
   directory, is readable, and has at least one subdirectory. Any
   failure is reported in one line and the script exits non-zero
   without touching any star.
2. Smoke test: the first star(s) that pass the file-selection step are
   run end-to-end with a very short MCMC chain (20 walkers x 30 steps),
   and the expected outputs (EW linelists, mcmc_fit.log, sampler.pickle,
   trace_plot.png, ep_vs_abundance.png) are checked. Up to
   --smoke-candidates stars are tried (one star's bad data should not
   look like a broken install); if none passes the batch aborts before
   the slow runs. The smoke test writes into that star's own directory;
   the real run overwrites those files.
Resuming: a star is skipped (`skip_done`) if a previous batch run left
batch_done.json in its directory, or, for runs made before that marker existed,
if it has linelist_with_ew_common.txt, sampler.pickle and ep_vs_abundance.png
and no smoke-test log. So after an interruption just run the same command again.
--redo ignores this; --start-at N additionally skips the first N-1 stars. The
smoke test never uses a finished star. Delete a star's batch_done.json to force
a re-fit of just that star.

--smoke-test-only stops after the smoke test; --skip-smoke-test omits it;
--dry-run only reports which files would be picked (no measurement).

Missing directory / missing data
--------------------------------
None of these raise a traceback; each becomes a summary row:
  skip_no_fits     directory has no .fits files
  skip_unreadable  every .fits file is empty or cannot be opened
  skip_no_match    readable files, but none names the directory's star
  skip_missing_dir the directory vanished or cannot be listed
A mix of readable and unreadable files runs on the readable ones and
notes the rest. The summary CSV is rewritten after every star, so an
interrupted run (Ctrl-C included) keeps what finished.

Output
------
Writes one line per star directory to a summary CSV (default:
batch_summary.csv in the current directory; override with --summary),
with columns: star_dir, fits_file, status, elapsed_s, note. `status` is
one of: ok, failed, skip_no_fits, skip_unreadable, skip_no_match,
skip_missing_dir, interrupted. `fits_file` is a
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
import json
import os
import subprocess
import sys
import time

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
    from astropy.io import fits  # lazy: preflight() reports a missing install cleanly
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


def _fits_name_or_error(fits_path):
    """-> (header star name or None, error string or None). An empty or
    corrupt file is an error, which is different from a readable file
    that simply has no usable name keyword."""
    try:
        if os.path.getsize(fits_path) == 0:
            return None, 'empty file (0 bytes)'
    except OSError as e:
        return None, f'cannot stat: {e}'
    from astropy.io import fits
    try:
        with fits.open(fits_path) as hdul:
            header = hdul[0].header
    except Exception as e:
        return None, f'cannot open: {type(e).__name__}: {e}'
    return _header_star_name(fits_path), None


def pick_target_fits(star_dir):
    """One star directory -> (list of fits_path, status, note).

    status is one of: 'ok', 'skip_no_fits', 'skip_unreadable',
    'skip_no_match', 'skip_missing_dir'. See the module docstring for the
    selection rules -- every readable .fits file is header-checked against
    the directory name, and EVERY matching file is returned
    (run_pipeline_cli.py merges however many target files it is given).
    """
    try:
        fits_files = sorted(
            f for f in os.listdir(star_dir) if f.lower().endswith('.fits'))
    except OSError as e:
        return [], 'skip_missing_dir', f'cannot list directory: {e}'
    dir_name = _normalize(os.path.basename(star_dir))

    if not fits_files:
        return [], 'skip_no_fits', 'no .fits files in directory'

    candidates, unreadable = {}, {}
    for f in fits_files:
        name, err = _fits_name_or_error(os.path.join(star_dir, f))
        if err:
            unreadable[f] = err
        else:
            candidates[f] = name
    bad_note = (' -- unreadable: ' + '; '.join(f'{f} ({e})' for f, e in unreadable.items())
                if unreadable else '')

    if not candidates:
        return [], 'skip_unreadable', (
            f'all {len(fits_files)} fits file(s) unreadable' + bad_note)

    matches = [f for f, name in candidates.items()
               if name and (dir_name in _normalize(name) or _normalize(name) in dir_name)]

    if matches:
        others = [f for f in candidates if f not in matches]
        note = (f'picked {len(matches)}/{len(fits_files)} file(s) matching directory: '
                 f'{[(f, candidates[f]) for f in matches]}')
        if others:
            note += f' -- ignored: {others}'
        return [os.path.join(star_dir, f) for f in matches], 'ok', note + bad_note

    note = (f'{len(candidates)} readable fits file(s), 0 header-name matches to directory '
             f'{os.path.basename(star_dir)!r}: '
             + ', '.join(f'{f}={n!r}' for f, n in candidates.items()))
    return [], 'skip_no_match', note + bad_note


def run_one_star(star_dir, solar_arg, dry_run=False, extra_args=(), log_name='batch_pipeline_stdout.log'):
    """Pick the target FITS file(s) for star_dir and, if any matched, run
    run_pipeline_cli.py on all of them (as one merged target) via a
    single subprocess -- or, if dry_run, just report which file(s) WOULD
    be picked. `extra_args` are appended to the run_pipeline_cli.py
    command (the smoke test uses them to shorten the MCMC).

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

    cmd = [sys.executable, RUN_PIPELINE_CLI, *target_paths, *extra_args]
    if solar_arg is not None:
        cmd += ['--solar', solar_arg]

    stdout_log_path = os.path.join(star_dir, log_name)
    t0 = time.time()
    try:
        with open(stdout_log_path, 'w') as f:
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    except OSError as e:
        # e.g. directory read-only: record it and keep going
        note = f'could not launch/log run: {e}'
        print(f'    FAILED -- {note}')
        return dict(star_dir=star_dir, fits_file=fits_file_field, status='failed',
                     elapsed_s='', note=note)
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


DONE_MARKER = 'batch_done.json'
# Outputs of a finished full run, for stars fit before the marker existed.
LEGACY_DONE_OUTPUTS = ['linelist_with_ew_common.txt', 'sampler.pickle', 'ep_vs_abundance.png']


def is_already_fit(star_dir):
    """-> (True/False, reason). A star counts as fit if a previous batch run left
    DONE_MARKER, or (older runs) it has every full-run output and no smoke-test
    log -- a smoke test also writes those outputs, from a deliberately tiny chain,
    so a directory with a smoke log and no marker is re-run to be safe."""
    if os.path.isfile(os.path.join(star_dir, DONE_MARKER)):
        return True, f'{DONE_MARKER} present'
    if os.path.isfile(os.path.join(star_dir, 'batch_smoke_test_stdout.log')):
        return False, ''
    if all(os.path.isfile(os.path.join(star_dir, f)) for f in LEGACY_DONE_OUTPUTS):
        return True, 'full-run outputs present (no marker; fit by an earlier run)'
    return False, ''


SUMMARY_FIELDS = ['star_dir', 'fits_file', 'status', 'elapsed_s', 'note']


def _write_summary(summary_path, rows):
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def run_batch(star_dirs, solar_arg, summary_path, dry_run=False, start_at=1, redo=False):
    rows = []
    star_dir = ''
    try:
        for i, star_dir in enumerate(star_dirs, 1):
            star_dir = os.path.abspath(star_dir)
            print(f'[{i}/{len(star_dirs)}] {os.path.basename(star_dir)}')
            if i < start_at:
                print(f'    skipped: before --start-at {start_at}')
                rows.append(dict(star_dir=star_dir, fits_file='', status='skip_before_start',
                                  elapsed_s='', note=f'index {i} < --start-at {start_at}'))
                _write_summary(summary_path, rows)
                continue
            done, why = (False, '') if redo else is_already_fit(star_dir)
            if done:
                print(f'    skipped: already fit ({why})')
                rows.append(dict(star_dir=star_dir, fits_file='', status='skip_done',
                                  elapsed_s='', note=f'already fit: {why}'))
                _write_summary(summary_path, rows)
                continue
            row = run_one_star(star_dir, solar_arg, dry_run=dry_run)
            if row['status'] == 'ok':
                try:
                    with open(os.path.join(star_dir, DONE_MARKER), 'w') as f:
                        json.dump(dict(fits_file=row['fits_file'], elapsed_s=row['elapsed_s'],
                                       finished=time.strftime('%Y-%m-%d %H:%M:%S')), f)
                except OSError as e:
                    row['note'] += f' (could not write {DONE_MARKER}: {e})'
            rows.append(row)
            _write_summary(summary_path, rows)  # keep progress if the run is cut short
    except KeyboardInterrupt:
        print('\nInterrupted -- summary so far is saved.')
        rows.append(dict(star_dir=star_dir, fits_file='', status='interrupted',
                          elapsed_s='', note='batch interrupted during this star'))
        _write_summary(summary_path, rows)

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


# ------------------------- preflight and smoke test -------------------------

SMOKE_ARGS = ['--n-walkers', '20', '--n-steps', '30', '--burn-in', '5', '--a-frac', '0.0']
SMOKE_EXPECTED_OUTPUTS = ['linelist_with_ew.txt', 'linelist_with_ew_common.txt',
                           'ew_measurement.log', 'mcmc_fit.log', 'sampler.pickle',
                           'trace_plot.png', 'ep_vs_abundance.png']


def list_star_dirs(root):
    """-> (sorted star directories, error message or None)."""
    if not os.path.exists(root):
        return [], f'parent directory does not exist: {root}'
    if not os.path.isdir(root):
        return [], f'not a directory: {root}'
    try:
        entries = os.listdir(root)
    except OSError as e:
        return [], f'cannot read parent directory {root}: {e}'
    star_dirs = sorted(os.path.join(root, d) for d in entries
                        if not d.startswith('.') and os.path.isdir(os.path.join(root, d)))
    if not star_dirs:
        return [], f'parent directory has no star subdirectories: {root}'
    return star_dirs, None


def preflight(solar_arg):
    """Cheap environment checks. -> list of problem strings (empty = OK)."""
    problems = []
    for mod in ('numpy', 'scipy', 'astropy', 'matplotlib', 'emcee', 'SPAE'):
        try:
            __import__(mod)
        except Exception as e:
            problems.append(f'cannot import {mod} in {sys.executable}: {type(e).__name__}: {e}')
    if problems:  # run_pipeline_cli imports the same packages; report those first
        return problems
    import run_pipeline_cli as cli
    # Version check: this driver needs the SPAE package to be the same
    # checkout as run_pipeline_cli.py. The per-line EW-error field was added
    # to the abundance records later than main; an older or separately
    # installed SPAE lacks it and the run fails only at the very end
    # (analyze_mcmc_stage, "no field of name abund_err").
    import numpy as np
    import SPAE
    from SPAE import abunds
    spae_dir = os.path.dirname(os.path.abspath(SPAE.__file__))
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    if 'abund_err' not in getattr(abunds, '_RESULT_DTYPE', np.dtype([])).names:
        problems.append(
            f'SPAE imported from {spae_dir} is too old: abunds._RESULT_DTYPE has no '
            f'abund_err (needs the abundance_errors changes). This script is in {repo_dir}; '
            f'check out the same branch there and `pip install -e {repo_dir}`.')
    elif os.path.dirname(spae_dir) != repo_dir:
        print(f'WARNING: SPAE is imported from {spae_dir}, not this checkout ({repo_dir}).')
    if not os.path.isfile(cli.LINELIST_PATH):
        problems.append(f'linelist missing: {cli.LINELIST_PATH}')
    solar = solar_arg if solar_arg is not None else cli.DEFAULT_SOLAR_DIR
    if not os.path.exists(solar):
        problems.append(f'solar reference missing: {solar}')
    elif os.path.isdir(solar) and not any(
            f.lower().endswith('.fits') for f in os.listdir(solar)):
        problems.append(f'solar reference directory has no .fits files: {solar}')
    return problems


def smoke_test(star_dirs, solar_arg, max_candidates):
    """Run up to `max_candidates` runnable stars with a very short chain.
    Returns True on the first one that finishes AND writes every expected
    output. Skipped/unreadable stars do not count against the limit."""
    print(f'=== Smoke test: short chain ({" ".join(SMOKE_ARGS)}) ===')
    tried = 0
    for star_dir in star_dirs:
        if tried >= max_candidates:
            break
        star_dir = os.path.abspath(star_dir)
        target_paths, status, note = pick_target_fits(star_dir)
        if status != 'ok' or is_already_fit(star_dir)[0]:
            continue  # never overwrite a finished star's outputs with a smoke run
        tried += 1
        print(f'  smoke candidate {tried}/{max_candidates}: {os.path.basename(star_dir)}')
        row = run_one_star(star_dir, solar_arg, extra_args=SMOKE_ARGS,
                            log_name='batch_smoke_test_stdout.log')
        if row['status'] != 'ok':
            continue
        missing = [f for f in SMOKE_EXPECTED_OUTPUTS
                   if not os.path.isfile(os.path.join(star_dir, f))]
        if missing:
            print(f'    exited 0 but missing outputs: {missing}')
            continue
        print('=== Smoke test PASSED ===')
        return True
    if tried == 0:
        print('Smoke test could not run: no star directory has readable, name-matching '
              'FITS data (see --dry-run for the per-directory reasons).')
    else:
        print(f'=== Smoke test FAILED on {tried} star(s) -- see batch_smoke_test_stdout.log '
              f'in the directories above ===')
    return False


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
             'EW measurement or MCMC actually runs.')
    parser.add_argument(
        '--smoke-test-only', action='store_true',
        help='Run preflight and the short-chain smoke test, then exit.')
    parser.add_argument(
        '--skip-smoke-test', action='store_true',
        help='Go straight to the full batch (preflight still runs).')
    parser.add_argument(
        '--smoke-candidates', type=int, default=3,
        help='Try up to this many stars in the smoke test before giving up (default 3).')
    parser.add_argument(
        '--start-at', type=int, default=1, metavar='N',
        help='Skip every star before the Nth one (1-based, in the order shown by the '
             '[i/N] progress lines). Rarely needed: finished stars are skipped '
             'automatically.')
    parser.add_argument(
        '--redo', action='store_true',
        help='Ignore the already-fit check and re-run every star.')
    args = parser.parse_args()

    root = os.path.abspath(args.root)
    star_dirs, err = list_star_dirs(root)
    if err:
        print(f'ERROR: {err}', file=sys.stderr)
        return 2
    print(f'Found {len(star_dirs)} star director{"y" if len(star_dirs) == 1 else "ies"} '
          f'under {root}')

    if args.dry_run:
        run_batch(star_dirs, args.solar_spectrum, os.path.abspath(args.summary), dry_run=True,
                  start_at=args.start_at, redo=args.redo)
        return 0

    problems = preflight(args.solar_spectrum)
    if problems:
        print('ERROR: preflight failed:', file=sys.stderr)
        for p in problems:
            print(f'  - {p}', file=sys.stderr)
        return 2
    print('Preflight OK.')

    pending = [d for d in star_dirs[max(args.start_at - 1, 0):]
               if args.redo or not is_already_fit(os.path.abspath(d))[0]]
    if not pending:
        print('Nothing to do: every star is already fit (use --redo to re-run).')
        return 0

    if not args.skip_smoke_test or args.smoke_test_only:
        if not smoke_test(star_dirs, args.solar_spectrum, args.smoke_candidates):
            return 3
        if args.smoke_test_only:
            return 0

    rows = run_batch(star_dirs, args.solar_spectrum, os.path.abspath(args.summary),
                     start_at=args.start_at, redo=args.redo)
    return 1 if any(r['status'] == 'failed' for r in rows) else 0


if __name__ == '__main__':
    sys.exit(main())
