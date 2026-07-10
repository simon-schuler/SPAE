"""
Line opacity at line center (kapnu0), Doppler widths (dopp), Voigt a-parameter,
and atmosphere level jtau5 where taulam ≈ 0.5.

Translated from Nearly.f, Damping.f, Gammabark.f, Trudamp.f.

Bug fixes vs Fortran Trudamp.f:
  - `iatom` is never defined in Fortran Trudamp.f (implicit real*8 → undefined).
    Replaced all Trudamp.f `iatom` references with `iatom10` (which IS defined).
  - `anature` on line 32 of Trudamp.f is undefined; should be `gnature`. Fixed.
"""
import os
import numpy as np

from .opacit import opacit
from .eqlib import _sunder


# ---------------------------------------------------------------------------
# Barklem damping data loader
# ---------------------------------------------------------------------------

def gammabark(state) -> None:
    """
    Read Barklem (2000, 2005) damping data and match to line list.

    Populates state.gambark[j], state.alpbark[j], state.gamrad[j].
    Unmatched lines get gambark[j] = -1 (signals use Unsold fallback).
    """
    # Pick file based on wavelength coverage (UV vs optical/IR)
    if state.nlines > 0 and state.wave1[state.nlines - 1] > 3000.0:
        barkfile = state.fbarklem
    else:
        barkfile = state.fbarklemUV

    if not barkfile or not os.path.isfile(barkfile):
        # No Barklem file → flag all lines as unmatched
        state.gambark[:state.nlines + state.nstrong] = -1.0
        state.alpbark[:state.nlines + state.nstrong] = -1.0
        state.gamrad[:state.nlines + state.nstrong]  =  0.0
        return

    wavebk = []
    idbk   = []
    gammabk  = []
    alphabk  = []
    gammarad = []
    with open(barkfile) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            wavebk.append(float(parts[0]))
            idbk.append(float(parts[1]))
            gammabk.append(float(parts[2]))
            alphabk.append(float(parts[3]))
            gammarad.append(float(parts[4]) if len(parts) >= 5 else 0.0)

    wavebk  = np.array(wavebk)
    idbk    = np.array(idbk)
    gammabk = np.array(gammabk)
    alphabk = np.array(alphabk)
    gammarad = np.array(gammarad)
    numbark = len(wavebk)

    # Wavelength range of line list
    ntotal = state.nlines + state.nstrong
    wavemin = np.min(state.wave1[:ntotal])
    wavemax = np.max(state.wave1[:ntotal])

    # Index limits in Barklem table
    nummin = 0
    for k in range(numbark):
        if wavemin - wavebk[k] < 1.0:
            nummin = k
            break
    nummax = numbark
    for k in range(nummin, numbark):
        if wavebk[k] - wavemax > 1.0:
            nummax = k
            break

    # Match each line
    for j in range(ntotal):
        state.gambark[j] = -1.0
        state.alpbark[j] = -1.0
        state.gamrad[j]  =  0.0
        if state.atom1[j] > 100.0:
            continue
        iatom10_j = round(10.0 * state.atom1[j])
        for k in range(nummin, nummax):
            waveerror = -(state.wave1[j] - wavebk[k]) / wavebk[k]
            if abs(waveerror) < 5.0e-6 and round(10.0 * idbk[k]) == iatom10_j:
                state.gamrad[j]  = gammarad[k]
                state.gambark[j] = 10.0 ** gammabk[k]
                state.alpbark[j] = (1.0 - alphabk[k]) / 2.0
                break
            if waveerror > 5.0e-6:
                break


# ---------------------------------------------------------------------------
# Accurate damping for specific well-studied lines (Trudamp.f)
# Bug fixes: iatom → iatom10 throughout; anature → gnature on Ca II K line.
# ---------------------------------------------------------------------------

def _trudamp(state, j: int) -> None:
    """
    Compute a(j,:) for lines with accurately known laboratory damping.
    Overwrites state.a[j, :ntau].
    """
    ntau    = state.ntau
    wave_j  = state.wave1[j]
    iwave   = int(wave_j)
    ich     = round(state.charge[j])
    iatom10 = round(10.0 * state.atom1[j])

    unsold = abs(
        1.61e-33 * (13.5 * state.charge[j] / (state.chi[j, ich-1] - state.e[j, 0]))**2
        - 1.61e-33 * (13.5 * state.charge[j] / (state.chi[j, ich-1] - state.e[j, 1]))**2
    )

    t   = state.t[:ntau]
    ne  = state.ne[:ntau]
    nhtot = state.nhtot[:ntau]
    nhe = state.xabund[1] * nhtot   # He number density

    if iatom10 == 201 and iwave == 3933:
        # Ca II K line
        gnature  = 1.45e8
        gvander  = 1.6e-8 * (t / 5000.0)**0.3 * nhtot
        gstark   = 3.0e-6 * ne
        gnature_arr = np.full(ntau, gnature)
        gammadamp   = gnature_arr + gvander + gstark   # fixed: anature → gnature
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])

    elif iatom10 == 201 and iwave in (8498, 8542, 8662):
        # Ca II IR triplet
        gnature = 1.5e8
        gstark  = 1.5e-6 * ne * (t / 5000.0)**0.1666
        ghelium = 3.0e-9 * nhe * (t / 5000.0)**0.4
        ghydro  = 1.0e-8 * nhtot * (t / 5000.0)**0.4
        gammadamp = gnature / 2.0 + gstark + ghelium + ghydro
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])

    elif iatom10 == 200 and iwave == 6717:
        # Ca I 6717 Å
        gnature = 0.4e-8
        ghelium = 1.0e-9 * nhe   * (t / 5000.0)**0.4
        ghydro  = 2.0e-8 * nhtot * (t / 5000.0)**0.4
        gammadamp = (gnature / 2.0 + ghelium + ghydro) * 2.0
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])

    elif iatom10 == 200 and iwave in (6318, 6343, 6361):
        # Ca I autoionization lines
        gnature = (state.dampnum[j] if state.dampnum[j] != 0 else 1.0) * 1.5e12
        gammadamp = np.full(ntau, gnature)
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])

    elif iatom10 == 110 and iwave != 0:
        # Na I lines
        gnature = 2.21e15 / wave_j**2
        v1      = np.sqrt(2.1175e8 * t * (1.0 / state.amass[j] + 1.008))
        gvander = 17.0 * unsold**0.4 * v1**0.6 * nhtot
        gcoll   = gvander * 2.1
        gammadamp = gnature / 2.0 + gcoll
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])

    elif iatom10 == 1060 and iwave == 3693:
        # CH autoionization at 3693 Å
        gnature = (state.dampnum[j] if state.dampnum[j] != 0 else 1.0) * 4.0e11
        gammadamp = np.full(ntau, gnature)
        state.a[j, :ntau] = gammadamp * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])


# ---------------------------------------------------------------------------
# Per-line damping (Damping.f)
# ---------------------------------------------------------------------------

def _damping(state, j: int) -> None:
    """
    Compute Voigt a-parameter for line j at all depth points.
    Fills state.a[j, :ntau] and sets state.gammar/gammas/gammav/gammatot
    at each depth.
    """
    ntau    = state.ntau
    wave_j  = state.wave1[j]
    iwave   = int(wave_j)
    iatom10 = round(10.0 * state.atom1[j])
    ich     = round(state.charge[j])

    # Convert dampnum from log if negative (Fortran: if dampnum < 0 → 10^dampnum)
    if state.dampnum[j] < 0.0:
        state.dampnum[j] = 10.0 ** state.dampnum[j]

    # Delegate special lines to trudamp (only if itru == 0)
    if state.itru == 0:
        is_caII_ir  = (iatom10 == 201 and iwave in (8498, 8542, 8662, 3933))
        is_ch_auto  = (iatom10 == 1060 and iwave == 3693)
        is_ca1_6717 = (iatom10 == 200 and iwave == 6717)
        is_ca1_auto = (iatom10 == 200 and iwave in (6318, 6343, 6361))
        if is_caII_ir or is_ch_auto or is_ca1_6717 or is_ca1_auto:
            _trudamp(state, j)
            return

    t    = state.t[:ntau]
    ne   = state.ne[:ntau]
    nh1  = state.numdens[0, 0, :ntau]   # H I number density
    nhe1 = state.numdens[1, 0, :ntau]   # He I number density
    nh2  = state.numdens[7, 0, :ntau]   # H2 number density

    v1 = np.sqrt(2.1175e8 * t * (1.0 / state.amass[j] + 1.008))

    # Unsold C6 for van der Waals damping
    ebreakup = 7.0 if state.atom1[j] > 100.0 else state.chi[j, ich - 1]
    e_lo = state.e[j, 0]
    e_up = state.e[j, 1]
    if e_lo >= ebreakup or e_up >= ebreakup:
        unsold = 1.0e-33
    else:
        unsold = abs(
            1.61e-33 * (13.598 * state.charge[j] / (ebreakup - e_lo))**2
            - 1.61e-33 * (13.598 * state.charge[j] / (ebreakup - e_up))**2
        )

    opt  = state.dampingopt
    gbar = state.gambark[j]

    if opt == 0 or (opt == 1 and gbar < 0):
        dn = state.dampnum[j]
        if dn == 0.0:
            gammav = 17.0 * unsold**0.4 * v1**0.6 * nh1
        elif dn < 1.0e-15:
            gammav = 17.0 * dn**0.4 * v1**0.6 * nh1
        elif dn < 1.0e-4:
            gammav = dn * (t / 10000.0)**0.3 * nh1
        else:
            gammav = 17.0 * (unsold * dn)**0.4 * v1**0.6 * nh1

    elif opt == 1 and gbar > 0.0:
        gammav = gbar * (t / 10000.0)**state.alpbark[j] * nh1

    elif opt == 2:
        gammav = 17.0 * ((1.0 + 0.67 * e_lo) * unsold)**0.4 * v1**0.6 * nh1

    elif opt == 3:
        dn = state.dampnum[j]
        if dn <= 1.0e-10:
            dn = 1.0
        c6h  = abs(1.01e-32 * state.charge[j]**2 * (13.598 / (ebreakup - e_lo))**2
                   - 1.61e-33 * (13.598 / (ebreakup - e_up))**2)
        c6he = abs((0.204956 / 0.666793) * 1.01e-32 * state.charge[j]**2
                   * (13.598 / (ebreakup - e_lo))**2
                   - 1.61e-33 * (13.598 / (ebreakup - e_up))**2)
        c6ht = abs((0.806 / 0.666793) * 1.01e-32 * state.charge[j]**2
                   * (13.598 / (ebreakup - e_lo))**2
                   - 1.61e-33 * (13.598 / (ebreakup - e_up))**2)
        gammav = (17.0 * v1**0.6 *
                  (c6h**0.4 * nh1 + c6he**0.4 * nhe1 + c6ht**0.4 * nh2)
                  * dn**0.4)
    else:
        gammav = 17.0 * unsold**0.4 * v1**0.6 * nh1

    # Radiative broadening
    if state.gamrad[j] != 0.0 and opt == 1:
        gammar = np.full(ntau, state.gamrad[j])
    else:
        gammar = np.full(ntau, 2.223e15 / wave_j**2)

    # Stark broadening
    excdiff = state.chi[j, ich - 1] - e_up
    if excdiff > 0.0 and state.atom1[j] < 100.0:
        effn2 = 13.6 * state.charge[j]**2 / excdiff
    else:
        effn2 = 25.0
    gammas = 1.0e-8 * ne * effn2**2.5

    gammatot = gammar + gammas + gammav
    state.a[j, :ntau] = gammatot * wave_j * 1.0e-8 / (12.56636 * state.dopp[j, :ntau])


# ---------------------------------------------------------------------------
# Main routine
# ---------------------------------------------------------------------------

def nearly(state, numpass: int) -> None:
    """
    Compute kapnu0, dopp, a for lines; find jtau5.

    numpass == 1 : all lines, recompute dopp/damping, find jtau5
    numpass == 2 : only lines [lim1line..lim2line], don't touch dopp/damping
    numpass == 3 : only line 0, recompute dopp/damping, find jtau5
    """
    ntau   = state.ntau

    # Load Barklem data on first full pass
    if numpass == 1 and state.dampingopt == 1:
        gammabark(state)

    # Find jtau5 (depth where taulam ≥ 0.5 at the first line's wavelength)
    if numpass in (1, 3):
        opacit(state, 2, state.wave1[0])
        jtau5 = ntau - 1   # default: last layer
        for i in range(ntau):
            if state.taulam[i] >= 0.5:
                jtau5 = i
                break
        state.jtau5 = jtau5

    # Line index range
    if numpass == 1:
        j1, j2 = 0, state.nlines + state.nstrong - 1
    elif numpass == 2:
        j1, j2 = state.lim1line, state.lim2line
    else:  # numpass == 3
        j1, j2 = 0, 0

    for j in range(j1, j2 + 1):
        ich   = int(state.charge[j] + 0.1)
        iatom = int(state.atom1[j] + 0.0001)

        factoriso = 1.0

        # Doppler widths and damping (only on first/fake-line passes)
        if numpass in (1, 3):
            t    = state.t[:ntau]
            dopp = np.sqrt(1.6631e8 * t / state.amass[j] + state.vturb[:ntau]**2)
            state.dopp[j, :ntau] = dopp
            _damping(state, j)

        # Lower state number densities
        tkev = state.tkev[:ntau]
        t    = state.t[:ntau]

        if iatom < 100:
            # Atomic line: Saha + Boltzmann
            u1  = state.u[iatom - 1, 0, :ntau]
            u2  = state.u[iatom - 1, 1, :ntau]
            u3  = state.u[iatom - 1, 2, :ntau]
            u4  = state.u[iatom - 1, 3, :ntau]

            q21 = 4.825e15 * u2 / (u1 * state.ne[:ntau]) * t**1.5 * np.exp(-state.chi[j, 0] / tkev)
            q32 = 4.825e15 * u3 / (u2 * state.ne[:ntau]) * t**1.5 * np.exp(-state.chi[j, 1] / tkev)
            q43 = 4.825e15 * u4 / (u3 * state.ne[:ntau]) * t**1.5 * np.exp(-state.chi[j, 2] / tkev)

            if ich == 1:
                q = 1.0 + q21 + q32 * q21 + q43 * q32 * q21
            elif ich == 2:
                q = 1.0 / q21 + 1.0 + q32 + q43 * q32
            else:
                q = 1.0 / (q21 * q32) + 1.0 / q32 + 1.0 + q43

            if state.control.strip() == 'abandy':
                xxab = state.xabund[iatom - 1] * 10.0**state.deltaabund
            else:
                xxab = state.xabund[iatom - 1]

            xnum = xxab * state.nhtot[:ntau] / q * np.exp(-state.e[j, 0] / tkev) / state.u[iatom - 1, ich - 1, :ntau]

        elif iatom < 10000:
            # Diatomic molecular line: Bates-Damgaard dissociation equilibrium
            iaa, ibb_code = _sunder(state.atom1[j])
            ibb = int(ibb_code) if ibb_code else 0

            ia = ib = -1
            for n in range(state.neq):
                if state.iorder[n] == iaa:
                    ia = n
                if state.iorder[n] == ibb:
                    ib = n

            uaa = state.u[iaa - 1, 0, :ntau] if iaa > 0 else np.ones(ntau)
            ubb = state.u[ibb - 1, 0, :ntau] if ibb > 0 else np.ones(ntau)

            psipri = (1.38065e-16 * t *
                      10.0**(state.d0[j] * state.theta[:ntau] - 13.670) *
                      state.theta[:ntau]**2.5 /
                      (state.rdmass[j]**1.5 * uaa * ubb))

            xa = state.xamol[ia, :ntau] if ia >= 0 else np.zeros(ntau)
            xb = state.xamol[ib, :ntau] if ib >= 0 else np.zeros(ntau)
            xnum = np.exp(-state.e[j, 0] / tkev) * psipri * xa * xb

        else:
            # Triatomic: H2O (10108) or CO2 (60808)
            if iatom == 10108:
                xnum = state.xnh2o[:ntau] / state.uh2o[:ntau] * np.exp(-state.e[j, 0] / tkev)
            elif iatom == 60808:
                xnum = state.xnco2[:ntau] / state.uco2[:ntau] * np.exp(-state.e[j, 0] / tkev)
            else:
                raise ValueError(f"Unsupported triatomic species code: {iatom}")

        # Isotope abundance factor
        if abs(state.atom1[j] - float(iatom)) >= 0.0:
            for n in range(state.numiso):
                if abs(state.atom1[j] - state.isotope[n]) < 1.0e-9:
                    factoriso = state.isoabund[n, state.isorun - 1]
                    break

        # Line opacity at line center
        stim = 1.0 - np.exp(-1.43879e8 / (state.wave1[j] * t))
        state.kapnu0[j, :ntau] = (
            2.65386e-2 * xnum * state.gf[j] * state.wave1[j] * 1.0e-8 /
            state.dopp[j, :ntau] * stim / factoriso
        )

        state.strength[j] = state.kapnu0[j, state.jtau5]
