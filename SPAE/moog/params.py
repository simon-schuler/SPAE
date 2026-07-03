"""
Parameter file reader translated from Begin.f / Params.f.

Reads a MOOG batch.par file into *state*, setting the control mode,
file paths, option flags, and synthesis limits.

Public API
----------
params(state, param_file='batch.par')
    Parse the parameter file and populate state fields.
"""


def _parse_token(s: str) -> str:
    """
    Strip outer single or double quotes from a Fortran list-directed string.
    Handles forms like 'star.mod' or "star.mod" or bare star.mod.
    """
    s = s.strip()
    if not s:
        return ''
    if len(s) >= 2 and s[0] in ("'", '"') and s[-1] == s[0]:
        return s[1:-1]
    if s[0] in ("'", '"'):
        q = s[0]
        end = s.find(q, 1)
        return s[1:end] if end > 0 else s[1:]
    return s.split()[0]


def _init_defaults(state) -> None:
    """Apply Params.f initialization block to state."""
    # File paths
    state.f1out   = 'no_filename_given'
    state.f2out   = 'no_filename_given'
    state.f3out   = 'no_filename_given'
    state.f4out   = 'no_filename_given'
    state.f5out   = 'no_filename_given'
    state.fmodel  = 'no_filename_given'
    state.flines  = 'no_filename_given'
    state.fslines = 'no_filename_given'
    state.fobs    = 'no_filename_given'

    # Option flags (Params.f defaults, some differ from State dataclass defaults)
    state.modprintopt = 1
    state.molopt      = 1
    state.linprintopt = 1
    state.linprintalt = 1
    state.fluxintopt  = 0
    state.plotopt     = 0
    state.dampingopt  = 0
    state.specfileopt = 0
    state.linfileopt  = 0
    state.iunits      = 0
    state.itru        = 0
    state.iraf        = 0
    state.scatopt     = 0
    state.gfstyle     = 0
    state.dostrong    = 0
    state.molset      = 0
    state.fudge       = -1.0

    # Synthesis range
    state.start    = 0.0
    state.sstop    = 0.0
    state.step     = 0.0
    state.delta    = 0.0
    state.cogatom  = 0.0
    state.contnorm = 1.0
    state.oldstart = 0.0
    state.oldstop  = 0.0
    state.oldstep  = 0.0
    state.olddelta = 0.0
    state.delwave  = 0.0
    state.wavestep = 0.0
    state.rwlow    = 0.0
    state.rwhigh   = 0.0
    state.rwstep   = 0.0

    # Line index limits
    state.ncurve   = 0
    state.lim1line = 0
    state.lim2line = 0
    state.lim1obs  = 0
    state.lim2obs  = 0
    state.lim1     = 0
    state.lim2     = 0

    # Abundance overrides (reset each call, as in Params.f goto-4 block)
    state.numpecatom    = 0
    state.numatomsyn    = 0
    state.ninetynineflag = 0
    state.pec[:]        = 0
    state.pecabund[:]   = 0.0
    state.abfactor[:]   = 0.0
    state.numiso        = 0
    state.numisosyn     = 0
    state.isotope[:]    = 0.0
    state.isoabund[:]   = 0.0


def params(state, param_file: str = 'batch.par') -> None:
    """
    Read a MOOG parameter file into *state*.

    The first line of the file is the 7-character control keyword (the run
    driver, e.g. 'abfind', 'synth  ').  All remaining lines are keyword/value
    pairs as defined in Params.f.

    Parameters
    ----------
    state : State
        MOOG State dataclass.
    param_file : str
        Path to the batch parameter file (default: ``'batch.par'``).
    """
    _init_defaults(state)

    with open(param_file) as fh:
        raw_lines = fh.readlines()

    if not raw_lines:
        raise ValueError(f"Empty parameter file: {param_file}")

    # First line: control keyword, padded/truncated to 7 chars (Fortran a7 format)
    first = raw_lines[0].rstrip('\n')
    state.control = (first + '       ')[:7]

    i = 1
    while i < len(raw_lines):
        raw = raw_lines[i]
        i += 1
        stripped = raw.strip()
        if not stripped:
            continue

        parts   = stripped.split(None, 1)
        keyword = parts[0]
        rest    = parts[1] if len(parts) > 1 else ''

        # ---- output files ----
        if keyword == 'standard_out':
            state.f1out = _parse_token(rest)
        elif keyword == 'summary_out':
            state.f2out = _parse_token(rest)
        elif keyword == 'smoothed_out':
            state.f3out = _parse_token(rest)
        elif keyword == 'iraf_out':
            state.f4out = _parse_token(rest)
        elif keyword == 'hardpost_out':
            state.f5out = _parse_token(rest)

        # ---- input files ----
        elif keyword == 'model_in':
            state.fmodel  = _parse_token(rest)
        elif keyword == 'lines_in':
            state.flines  = _parse_token(rest)
        elif keyword == 'stronglines_in':
            state.fslines = _parse_token(rest)
        elif keyword == 'observed_in':
            state.fobs    = _parse_token(rest)

        # ---- simple integer/float options ----
        elif keyword == 'atmosphere':
            state.modprintopt = int(rest.split()[0])
        elif keyword == 'molecules':
            state.molopt = max(1, int(rest.split()[0]))   # 0 is converted to 1
        elif keyword == 'molset':
            state.molset = int(rest.split()[0])
        elif keyword == 'lines':
            v = int(rest.split()[0])
            state.linprintopt = v
            state.linprintalt = v
        elif keyword == 'freeform':
            state.linfileopt  = int(rest.split()[0])
        elif keyword == 'gfstyle':
            state.gfstyle     = int(rest.split()[0])
        elif keyword == 'flux/int':
            state.fluxintopt  = int(rest.split()[0])
        elif keyword == 'damping':
            state.dampingopt  = int(rest.split()[0])
        elif keyword == 'plot':
            state.plotopt     = int(rest.split()[0])
        elif keyword == 'units':
            state.iunits      = int(rest.split()[0])
        elif keyword == 'iraf':
            state.iraf        = int(rest.split()[0])
        elif keyword == 'scat':
            state.scatopt     = int(rest.split()[0])
        elif keyword == 'contnorm':
            state.contnorm    = float(rest.split()[0])
        elif keyword == 'strong':
            state.dostrong    = int(rest.split()[0])
        elif keyword == 'opacit':
            state.fudge       = float(rest.split()[0])
        elif keyword == 'trudamp':
            state.itru        = int(rest.split()[0])
        elif keyword == 'obspectrum':
            v = int(rest.split()[0])
            state.specfileopt = abs(v)

        # ---- multi-value keywords: values on the NEXT line (Params.f reads
        #      with read(nfparam,*) after parsing the keyword line) ----
        elif keyword == 'synlimits':
            vals = raw_lines[i].split()
            i += 1
            state.start    = float(vals[0])
            state.sstop    = float(vals[1])
            state.step     = float(vals[2])
            state.delta    = float(vals[3])
            state.oldstart = state.start
            state.oldstop  = state.sstop
            state.oldstep  = state.step
            state.olddelta = state.delta

        elif keyword == 'fluxlimits':
            vals = raw_lines[i].split()
            i += 1
            state.start = float(vals[0])
            state.sstop = float(vals[1])
            state.step  = float(vals[2])

        elif keyword == 'blenlimits':
            vals = raw_lines[i].split()
            i += 1
            state.delwave = float(vals[0])
            state.step    = float(vals[1])
            state.cogatom = float(vals[2])

        elif keyword == 'coglimits':
            vals = raw_lines[i].split()
            i += 1
            state.rwlow    = float(vals[0])
            state.rwhigh   = float(vals[1])
            state.rwstep   = float(vals[2])
            state.wavestep = float(vals[3])
            state.cogatom  = float(vals[4])

        elif keyword == 'weedlimits':
            vals = raw_lines[i].split()
            i += 1
            state.xratio = float(vals[0])

        # ---- abundance overrides ----
        elif keyword == 'abundances':
            vals = rest.split()
            numpec  = int(vals[0])
            numsyn  = int(vals[1])
            state.numpecatom = numpec
            state.numatomsyn = numsyn
            state.pec[:]     = 0
            state.pecabund[:] = 0.0
            state.abfactor[:] = 0.0
            state.ninetynineflag = 0
            jatom = 0
            for _ in range(numpec):
                av = raw_lines[i].split()
                i += 1
                jatom = int(float(av[0]))
                factors = [float(av[k + 1]) for k in range(numsyn)]
                if jatom == 99:
                    for k, fval in enumerate(factors):
                        state.abfactor[k] = fval
                else:
                    for k, fval in enumerate(factors):
                        state.pecabund[jatom - 1, k] = fval
                    state.pec[jatom - 1] = 1
            if numpec == 1 and jatom == 99:
                state.ninetynineflag = 1

        elif keyword == 'isotopes':
            vals = rest.split()
            numiso    = int(vals[0])
            numisosyn = int(vals[1])
            state.numiso    = numiso
            state.numisosyn = numisosyn
            state.isotope[:]   = 0.0
            state.isoabund[:]  = 0.0
            for j in range(numiso):
                av = raw_lines[i].split()
                i += 1
                state.isotope[j] = float(av[0])
                for k in range(numisosyn):
                    state.isoabund[j, k] = float(av[k + 1])

        elif keyword == 'plotpars':
            iscale = int(rest.split()[0])
            if iscale != 0:
                i += 1          # xlo, xhi, ylo, yhi line
                i += 1          # veladd, xadd, yadd, ymult line
                sp = raw_lines[i].split()
                i += 1
                # smtype fwhmgauss vsini limbdark vmac fwhmloren
                if len(sp) > 1:
                    state.fwhmgauss = float(sp[1])
                if len(sp) > 2:
                    state.vsini     = float(sp[2])
                if len(sp) > 3:
                    state.limbdark  = float(sp[3])
                if len(sp) > 4:
                    state.vmac      = float(sp[4])
                if len(sp) > 5:
                    state.fwhmloren = float(sp[5])

        # ---- keywords not needed for the core physics layer ----
        elif keyword in (
            'terminal', 'histogram', 'veladjust', 'deviations',
            'lumratio', 'deltaradvel',
            'bin_raw_out', 'bin_smo_out', 'summary_in',
            'keeplines_out', 'tosslines_out', 'speccomp_out',
            'popsyn_out', 'rawbin_out', 'smoobin_out',
            'table_in', 'table_out', 'RUN',
        ):
            pass

        else:
            pass   # unknown keyword: silently ignore (MOOG would stop)
