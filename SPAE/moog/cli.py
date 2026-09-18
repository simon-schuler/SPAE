#!/usr/bin/env python3
"""
pymoog — command-line interface for the pymoog spectral analysis package.

Usage:
    pymoog [batch.par]

Reads a MOOG-format batch.par (defaults to 'batch.par' in the current
directory, matching MOOGSILENT behaviour), runs the requested mode, and
writes output to moog_out.1 / moog_out.2 — or to whatever filenames are
specified by standard_out / summary_out in the parameter file.

Supported modes: abfind, blends, cog, doflux, ewfind, synth, weedout
"""
import os
import sys

from .state   import State
from .abfind  import abfind_from_files
from .blends  import blends_from_files
from .cog     import cog_from_files
from .doflux  import doflux_from_files
from .ewfind  import ewfind_from_files
from .synth   import synth_from_files
from .weedout import weedout_from_files
from .output  import (
    write_abfind_output,
    write_blends_output,
    write_cog_output,
    write_doflux_output,
    write_ewfind_output,
    write_synth_output,
    write_weedout_output,
)

_NO_FILE = 'no_filename_given'

_MODES = {
    'abfind':  (abfind_from_files,  write_abfind_output),
    'blends':  (blends_from_files,  write_blends_output),
    'cog':     (cog_from_files,     write_cog_output),
    'doflux':  (doflux_from_files,  write_doflux_output),
    'ewfind':  (ewfind_from_files,  write_ewfind_output),
    'synth':   (synth_from_files,   write_synth_output),
    'weedout': (weedout_from_files, write_weedout_output),
}


def main(argv=None):
    """Entry point for the pymoog command."""
    if argv is None:
        argv = sys.argv[1:]

    fparam = argv[0] if argv else 'batch.par'
    fparam = os.path.abspath(fparam)
    if not os.path.isfile(fparam):
        sys.exit(f"pymoog: cannot find parameter file '{fparam}'")

    # cd to batch.par directory so relative model/linelist paths resolve
    work_dir = os.path.dirname(fparam)
    os.chdir(work_dir)
    fparam_local = os.path.basename(fparam)

    # Peek at the first line to determine the mode (Fortran a7 control keyword)
    with open(fparam_local) as fh:
        mode_key = fh.readline().strip().lower()

    if mode_key not in _MODES:
        sys.exit(
            f"pymoog: unrecognized mode '{mode_key}'\n"
            f"supported: {', '.join(sorted(_MODES))}"
        )

    run_fn, write_fn = _MODES[mode_key]

    state = State()
    state.fparam = fparam_local

    try:
        result = run_fn(state)
    except Exception as exc:
        sys.exit(f"pymoog {mode_key} failed: {exc}")

    f1path = state.f1out if state.f1out != _NO_FILE else 'moog_out.1'
    f2path = state.f2out if state.f2out != _NO_FILE else 'moog_out.2'

    write_fn(state, result, f1path, f2path)

    written = []
    if f1path != _NO_FILE: written.append(f1path)
    if f2path != _NO_FILE: written.append(f2path)
    if written:
        print(f"pymoog: wrote {', '.join(written)}")


if __name__ == '__main__':
    main()
