"""
Model atmosphere reader translated from Inmodel.f.

Supports all MOOG model types: KURUCZ, BEGN, NEWMARCS, WEBMARCS,
WEB2MARC, NEXTGEN, KURTYPE, KUR-PADOVA, GENERIC.

SPAE uses the KURUCZ type; the others are provided for completeness.
"""

import numpy as np
from .atomic_data import XSOLAR, XAM
from .partition import partfn
from .math_utils import rinteg

_BOLTZMANN = 1.38054e-16   # erg K⁻¹
_AMU       = 1.6606e-24    # g


def _read_label_then_value(f):
    """
    Read one line; skip the first 10 characters (label field);
    return what remains as a stripped string.  Matches Fortran:
        read(nfmodel,2002) list
        list2 = list(11:)
    """
    line = f.readline()
    if not line:
        return ''
    return line[10:].strip()


def _parse_vturb_line(f):
    """
    Read the vturb line written in Fortran fixed format (6d13.0).
    Returns a list of up to 6 float values; missing fields are 0.
    """
    line = f.readline().rstrip('\n')
    # Fortran d13.0: each field is exactly 13 characters wide
    vals = []
    for start in range(0, 13 * 6, 13):
        field = line[start:start + 13].strip()
        try:
            vals.append(float(field.replace('d', 'e').replace('D', 'E')))
        except ValueError:
            vals.append(0.0)
    return vals


def inmodel(state, eqlib_func=None) -> None:
    """
    Read a MOOG model atmosphere file into *state*.

    The filename is taken from ``state.infile`` (set by the params reader).
    After this call the following State fields are populated:

        Atmosphere: t, theta, tkev, tlog, pgas, ne, rhox/tauref, rho,
                    vturb, molweight, nhtot, kapref, xref
        Abundances: xabund, xabu
        Partition:  u[Z-1, ion, k]

    Parameters
    ----------
    state : State
        MOOG State dataclass.
    eqlib_func : callable, optional
        Molecular equilibrium function ``eqlib(state)``.  If *None*, a
        no-op stub is used (molecules all set to zero).
    """
    if eqlib_func is None:
        def eqlib_func(s):
            pass

    state.modelnum += 1

    with open(state.fmodel) as f:
        # -----------------------------------------------------------------
        # 1. Model type keyword (first 10 chars, right-padded with spaces)
        # -----------------------------------------------------------------
        raw = f.readline()
        modtype = (raw.rstrip('\n') + '          ')[:10]
        state.modtype = modtype

        # -----------------------------------------------------------------
        # 2. Title comment line
        # -----------------------------------------------------------------
        state.moditle = f.readline().rstrip('\n')

        # -----------------------------------------------------------------
        # 3. NTAU line  (label in chars 1-10, number from char 11 onward)
        # -----------------------------------------------------------------
        ntau = int(_read_label_then_value(f))
        if ntau > 100:
            raise ValueError(f"HOUSTON, WE HAVE MORE THAN 100 DEPTH POINTS! (ntau={ntau})")
        state.ntau = ntau

        # -----------------------------------------------------------------
        # 4. Model data lines  (format depends on modtype)
        # -----------------------------------------------------------------
        kaprefmass = np.zeros(ntau)
        wavref = 5000.0

        mtype = modtype.strip()

        if mtype == 'NEWMARCS':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                state.tauref[i] = float(parts[0])
                state.t[i]      = float(parts[1])
                state.ne[i]     = float(parts[2])
                state.pgas[i]   = float(parts[3])
                state.rho[i]    = float(parts[4])
                state.vturb[0]  = float(parts[5])
                kaprefmass[i]   = float(parts[6])

        elif mtype == 'WEBMARCS':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                # layer, log(tauRoss), log(tau5000), depth, t, pe, pgas
                state.tauref[i] = float(parts[2])
                state.t[i]      = float(parts[4])
                state.ne[i]     = float(parts[5])
                state.pgas[i]   = float(parts[6])

        elif mtype == 'WEB2MARC':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                # layer, log(tau5000), t, log(Pe), log(Pgas), rhox
                state.tauref[i] = float(parts[1])
                state.t[i]      = float(parts[2])
                state.ne[i]     = float(parts[3])
                state.pgas[i]   = float(parts[4])
                state.rhox[i]   = float(parts[5])

        elif mtype == 'KURUCZ':
            for i in range(ntau):
                parts = f.readline().split()
                state.rhox[i]   = float(parts[0])
                state.t[i]      = float(parts[1])
                state.pgas[i]   = float(parts[2])
                state.ne[i]     = float(parts[3])
                kaprefmass[i]   = float(parts[4])

        elif mtype == 'NEXTGEN':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                state.tauref[i]    = float(parts[0])
                state.t[i]         = float(parts[1])
                state.pgas[i]      = float(parts[2])
                state.ne[i]        = float(parts[3])
                state.rho[i]       = float(parts[4])
                state.molweight[i] = float(parts[5])
                kaprefmass[i]      = float(parts[8])

        elif mtype == 'BEGN':
            for i in range(ntau):
                parts = f.readline().split()
                state.tauref[i]    = float(parts[0])
                state.t[i]         = float(parts[1])
                state.pgas[i]      = float(parts[2])
                state.ne[i]        = float(parts[3])
                state.molweight[i] = float(parts[4])
                kaprefmass[i]      = float(parts[5])

        elif mtype == 'KURTYPE':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                state.rhox[i]  = float(parts[0])
                state.t[i]     = float(parts[1])
                state.pgas[i]  = float(parts[2])
                state.ne[i]    = float(parts[3])

        elif mtype == 'KUR-PADOVA':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                state.tauref[i]  = float(parts[0])
                state.t[i]       = float(parts[1])
                kaprefmass[i]    = float(parts[2])
                state.ne[i]      = float(parts[3])
                state.pgas[i]    = float(parts[4])
                state.rho[i]     = float(parts[5])

        elif mtype == 'GENERIC':
            wavref = float(f.readline().split()[0])
            for i in range(ntau):
                parts = f.readline().split()
                state.tauref[i] = float(parts[0])
                state.t[i]      = float(parts[1])
                state.pgas[i]   = float(parts[2])
                state.ne[i]     = float(parts[3])

        else:
            raise ValueError(
                f"Unknown MOOG model type '{mtype}'.  "
                "Permitted: KURUCZ, BEGN, KURTYPE, KUR-PADOVA, "
                "NEWMARCS, WEBMARCS, NEXTGEN, WEB2MARC, GENERIC"
            )

        state.wavref = wavref

        # -----------------------------------------------------------------
        # 5. Derived temperature quantities
        # -----------------------------------------------------------------
        for i in range(ntau):
            state.theta[i] = 5040.0 / state.t[i]
            state.tkev[i]  = 8.6171e-5 * state.t[i]
            state.tlog[i]  = np.log(state.t[i])

        # -----------------------------------------------------------------
        # 6. Convert log Pgas scale if values span less than a factor of 10
        # -----------------------------------------------------------------
        if state.pgas[ntau - 1] / state.pgas[0] < 10.0:
            state.pgas[:ntau] = 10.0 ** state.pgas[:ntau]

        # -----------------------------------------------------------------
        # 7. Convert log Ne scale if values span less than a factor of 20
        # -----------------------------------------------------------------
        if state.ne[ntau - 1] / state.ne[0] < 20.0:
            state.ne[:ntau] = 10.0 ** state.ne[:ntau]

        # -----------------------------------------------------------------
        # 8. Convert Pe → Ne  (Ne is actually Pe when ne < 1e7)
        # -----------------------------------------------------------------
        if state.ne[ntau - 1] < 1.0e7:
            for i in range(ntau):
                state.ne[i] /= _BOLTZMANN * state.t[i]

        # -----------------------------------------------------------------
        # 9. Partition functions for all elements at all depths
        # -----------------------------------------------------------------
        partfn(state)

        # -----------------------------------------------------------------
        # 10. Microturbulence  (Fortran format: 6d13.0 per line)
        # -----------------------------------------------------------------
        vtvals = _parse_vturb_line(f)
        state.vturb[0] = vtvals[0]
        if vtvals[1] != 0.0:
            # Per-layer vturb: read additional values up to ntau
            for i in range(1, 6):
                state.vturb[i] = vtvals[i]
            # Read remaining lines as needed
            idx = 6
            while idx < ntau:
                more = _parse_vturb_line(f)
                for v in more:
                    if idx < ntau:
                        state.vturb[idx] = v
                        idx += 1
        else:
            # Constant vturb for all layers
            state.vturb[:ntau] = state.vturb[0]

        # Convert km/s → cm/s if needed
        if state.vturb[0] < 100.0:
            state.vturb[:ntau] *= 1.0e5

        # -----------------------------------------------------------------
        # 11. Abundance overrides  (NATOMS + abscale line, then element data)
        # -----------------------------------------------------------------
        natoms_line = _read_label_then_value(f)
        parts = natoms_line.split()
        natoms  = int(parts[0])
        abscale = float(parts[1])

        overrides = {}
        if natoms != 0:
            for _ in range(natoms):
                ap = f.readline().split()
                Z_flt    = float(ap[0])
                logeps   = float(ap[1])
                overrides[int(Z_flt)] = logeps

        # Set default xabund from solar values + metallicity offset
        xhyd = 10.0 ** XSOLAR[0]          # 10^12 (log H = 12)
        state.xabund[0] = 1.0             # hydrogen by definition
        state.xabund[1] = 10.0 ** XSOLAR[1] / xhyd   # He
        for i in range(2, 95):
            state.xabund[i] = 10.0 ** (XSOLAR[i] + abscale) / xhyd
            state.xabu[i]   = state.xabund[i]

        # Apply per-element overrides
        for Z, logeps in overrides.items():
            idx = Z - 1
            state.xabund[idx] = 10.0 ** logeps / xhyd
            state.xabu[idx]   = state.xabund[idx]

        # -----------------------------------------------------------------
        # 12. Mean molecular weight (ignoring molecules)
        # -----------------------------------------------------------------
        wtnum = sum(state.xabund[i] * XAM[i] for i in range(95))
        wtden = sum(state.xabund[i]           for i in range(95))
        wtmol = wtnum / (XAM[0] * wtden)

        nomolweight = mtype in ('BEGN', 'NEXTGEN')
        if not nomolweight:
            state.molweight[:ntau] = wtmol

        # -----------------------------------------------------------------
        # 13. Density  (all types except NEXTGEN which provides rho directly)
        # -----------------------------------------------------------------
        if mtype != 'NEXTGEN':
            for i in range(ntau):
                state.rho[i] = (state.pgas[i] * state.molweight[i]
                                * _AMU / (_BOLTZMANN * state.t[i]))

        # -----------------------------------------------------------------
        # 14. Fictitious H number density (quadratic formula)
        # -----------------------------------------------------------------
        for i in range(ntau):
            th  = 5040.0 / state.t[i]
            ah2 = 10.0 ** (-(12.7422 + (-5.1137 + (0.1145 - 0.0091 * th) * th) * th))
            a1  = (1.0 + 2.0 * state.xabund[1]) * ah2
            b1  = 1.0 + state.xabund[1]
            c1  = -state.pgas[i]
            ph  = (-b1 / (2.0 * a1)
                   + np.sqrt((b1 / (2.0 * a1)) ** 2 - c1 / a1))
            state.nhtot[i] = (ph + 2.0 * ph * ph * ah2) / (_BOLTZMANN * state.t[i])

        # -----------------------------------------------------------------
        # 15. Molecule list  (NMOL header, then species codes)
        # -----------------------------------------------------------------
        moremol_line = _read_label_then_value(f)
        moremol = int(moremol_line.split()[0]) if moremol_line else 0

        if moremol > 0:
            extra_codes = []
            while len(extra_codes) < moremol:
                line = f.readline()
                extra_codes.extend(float(x) for x in line.split())

            # Use molset to choose base list
            if state.molset == 0:
                for i in range(110):
                    state.amol[i] = state.smallmollist[i]
                state.nmol = 30
            else:
                for i in range(110):
                    state.amol[i] = state.largemollist[i]
                state.nmol = 59

            # Append any extra molecules not already in the list
            existing = {round(state.amol[k]) for k in range(state.nmol)}
            for code in extra_codes:
                if round(code) not in existing:
                    state.amol[state.nmol] = code
                    existing.add(round(code))
                    state.nmol += 1

        else:
            if state.molset == 0:
                for i in range(110):
                    state.amol[i] = state.smallmollist[i]
                state.nmol = 30
            else:
                for i in range(110):
                    state.amol[i] = state.largemollist[i]
                state.nmol = 59

    # ---------------------------------------------------------------------
    # 16. Molecular equilibrium
    # ---------------------------------------------------------------------
    eqlib_func(state)

    # ---------------------------------------------------------------------
    # 17. Post-processing: tauref, kapref, xref  (by model type)
    # ---------------------------------------------------------------------
    if mtype == 'NEWMARCS':
        state.kapref[:ntau] = kaprefmass[:ntau] * state.rho[:ntau]

    elif mtype == 'KURUCZ':
        first = state.rhox[0] * kaprefmass[0]
        _, fint = rinteg(state.rhox, kaprefmass, ntau, first)
        state.tauref[:ntau] = np.cumsum(fint)
        state.kapref[:ntau] = kaprefmass[:ntau] * state.rho[:ntau]

    elif mtype == 'NEXTGEN':
        state.kapref[:ntau] = kaprefmass[:ntau] * state.rho[:ntau]

    elif mtype == 'BEGN':
        state.kapref[:ntau] = kaprefmass[:ntau] * state.rho[:ntau]

    elif mtype in ('KURTYPE',):
        raise NotImplementedError(
            "KURTYPE model requires opacit() — implement opacit.py first"
        )

    elif mtype == 'KUR-PADOVA':
        state.kapref[:ntau] = kaprefmass[:ntau] * state.rho[:ntau]

    elif mtype in ('GENERIC', 'WEBMARCS', 'WEB2MARC'):
        raise NotImplementedError(
            f"{mtype} model requires opacit() — implement opacit.py first"
        )

    # ---------------------------------------------------------------------
    # 18. Optical depth scale: convert linear ↔ log
    # ---------------------------------------------------------------------
    if state.tauref[0] < 0.0:
        # tauref is already in log scale → convert to linear
        for i in range(ntau):
            state.xref[i]   = state.tauref[i]
            state.tauref[i] = 10.0 ** state.xref[i]
    else:
        for i in range(ntau):
            state.xref[i] = np.log10(state.tauref[i])
