"""
Partition function calculations translated from Partfn.f / Ucalc.f / Partnew.f.

Public API
----------
partfn(state)
    Populate state.u[Z-1, ion, k] for all elements, ions, and depth layers.

partnew(iatom, ion, temp)
    Irwin (1981) polynomial partition function for updated species.

ucalc(j, temp)
    ATLAS9-encoded partition function for standard species.
"""

import numpy as np
from .atomic_data import NUDATA, PARTFLAG, NEWPARTDATA, XCHI1, XCHI2, XCHI3

_SCALE = np.array([0.001, 0.01, 0.1, 1.0])   # ucalc decoding scale table


def ucalc(j: int, temp: float) -> float:
    """
    Return partition function for species j at temperature temp.

    Parameters
    ----------
    j : int
        0-based species index into NUDATA  (= 4*(Z-1) + ion).
    temp : float
        Temperature [K].

    Notes
    -----
    Translated from Ucalc.f.  The NUDATA entries store two packed partition
    function values per 9-digit integer; the integer encoding is:

        k1    = val // 100000          (upper 4 digits, first PF value)
        k2    = val  % 100000          (lower 5 digits)
        k3    = k2  // 10              (lower 4 digits, second PF value)
        kscale= k2   % 10             (scale index 1..4 → _SCALE[kscale-1])

    Fortran nudata is 1-based → NUDATA[j, i-1] in Python (i = 1..5 used).

    The temperature grid is scaled by chix * 2000/11 where chix is the
    relevant ionisation potential.  it = 1..9 selects the temperature bin;
    the partition function is linearly interpolated between adjacent bins.
    """
    # Map species index j (0-based) → element Z and ion
    Z   = j // 4          # 0-based element index (= Z_atomic - 1)
    ion = j  % 4          # 0-based ion index

    # Ionisation potential for the relevant ion
    if ion == 0:
        chix = XCHI1[Z]
    elif ion == 1:
        chix = XCHI2[Z]
    elif ion == 2:
        chix = XCHI3[Z]
    else:
        return 1.0         # no data beyond 3rd ion → single ground state

    # For atoms with very large or unknown IPs, avoid division by zero
    if chix <= 0.0 or chix >= 90.0:
        return 1.0

    t2000 = chix * 2000.0 / 11.0
    if t2000 <= 0.0:
        return 1.0

    # Temperature bin index 1..9
    it = int(temp / t2000 - 0.5)
    it = max(1, min(9, it))
    dt = temp / t2000 - it - 0.5

    # Row into NUDATA (1-based → 0-based)
    i = (it + 1) // 2       # Fortran i (1-based); Python index i-1 below

    val = int(NUDATA[j, i - 1])    # nudata(i, j+1) in Fortran (0-based here)
    k1     = val // 100000
    k2     = val  % 100000
    k3     = k2  // 10
    kscale = k2   % 10
    if kscale < 1 or kscale > 4:
        kscale = 1

    pmin = 1.0    # Fortran: pmin = 1.  (minimum returned value)

    if it % 2 == 0:
        # Even it: p1 from k3 of row i; p2 from k1 of row i+1
        p1 = float(k3) * _SCALE[kscale - 1]
        val2    = int(NUDATA[j, i])          # nudata(i+1, j+1)
        kscale2 = val2 % 10
        if kscale2 < 1 or kscale2 > 4:
            kscale2 = 1
        k1_2 = val2 // 100000
        p2 = float(k1_2) * _SCALE[kscale2 - 1]
    else:
        # Odd it: both p1 and p2 from the same row i
        p1 = float(k1) * _SCALE[kscale - 1]
        p2 = float(k3) * _SCALE[kscale - 1]
        # Special case (Ucalc.f lines 49-53): at low T with kscale=1
        # and integer-valued p1==p2, clamp pmin to p1 to avoid drifting below it
        if dt < 0.0 and kscale == 1:
            kp1 = int(p1)
            if kp1 == int(p2 + 0.5):
                pmin = float(kp1)

    return max(pmin, p1 + (p2 - p1) * dt)


def partnew(iatom: int, ion: int, temp: float) -> float:
    """
    Partition function from Irwin (1981) polynomial fit.

    Parameters
    ----------
    iatom : int
        1-based element number (= Z_atomic).
    ion : int
        1-based ionisation state (1=neutral, 2=singly ionised, …).
    temp : float
        Temperature [K].

    Returns
    -------
    float
        Partition function U(T).

    Notes
    -----
    Translated from Partnew.f.  The polynomial is:
        log10(U) = sum_{j=1}^{6}  C_j * log10(T)^(j-1)
    where C_j come from NEWPARTDATA[row-1, j-1] and *row* is
    PARTFLAG[iatom-1, ion-1].
    """
    row = PARTFLAG[iatom - 1, ion - 1] - 1    # 0-based index into NEWPARTDATA
    if row < 0 or row >= len(NEWPARTDATA):
        return 1.0

    # MOOG stores these coefficients in natural-log form despite the Irwin
    # (1981) paper using log10.  Partnew.f uses tlog=log(T) and dexp().
    tlog = np.log(temp)
    coeffs = NEWPARTDATA[row]
    ulog = 0.0
    tpow = 1.0
    for c in coeffs:
        ulog += c * tpow
        tpow *= tlog

    return max(1.0, np.exp(ulog))


def partfn(state) -> None:
    """
    Compute partition functions for all 95 elements, 4 ionisation states,
    and all ntau depth layers; store in state.u[Z-1, ion, k].

    Translated from Partfn.f.  For each species, uses *partnew* if
    PARTFLAG[Z-1, ion] > 0, else *ucalc*.

    Parameters
    ----------
    state : State
        The moog State dataclass (must have ntau and t populated).
    """
    ntau = state.ntau
    for Z in range(1, 96):          # 1-based element number
        for ion in range(1, 5):     # 1-based ionisation state (1=I, 2=II, …)
            j = 4 * (Z - 1) + (ion - 1)   # 0-based NUDATA species index
            use_new = (PARTFLAG[Z - 1, ion - 1] > 0)
            for k in range(ntau):
                temp = state.t[k]
                if use_new:
                    state.u[Z - 1, ion - 1, k] = partnew(Z, ion, temp)
                else:
                    state.u[Z - 1, ion - 1, k] = ucalc(j, temp)
