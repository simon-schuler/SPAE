"""
Molecular equilibrium solver translated from Eqlib.f / Setmols.f / Bmolec.f.

Public API
----------
init_mol_data(state)
    Copy BLOCK DATA molecular constants into state (call once before inmodel).

eqlib(state)
    Solve the molecular dissociation equilibrium for all tau layers and store
    number densities in state.xmol, state.xamol, state.pmol, state.numdens.
"""

import numpy as np
from .atomic_data import XCHI1, XCHI2, XCHI3, PARTFLAG
from .partition import ucalc, partnew

# ---------------------------------------------------------------------------
# Module-level constants  (from Bmolec.f BLOCK DATA)
# ---------------------------------------------------------------------------

# Small default molecule list (30 species): atomic ions + molecules
# Codes: Z.1 = first ion of element Z; integer codes = molecules
_SMALLMOLLIST_RAW = [
    1.1,   2.1,   6.1,   7.1,   8.1,   9.1,
   12.1,  13.1,  14.1,  22.1,  26.1, 101.0,
  106.0, 107.0, 108.0, 109.0, 112.0, 113.0,
  114.0, 126.0, 606.0, 607.0, 608.0, 707.0,
  708.0, 808.0, 10108.0, 60808.0, 812.0, 822.0,
]
SMALLMOLLIST = np.zeros(110)
SMALLMOLLIST[:len(_SMALLMOLLIST_RAW)] = _SMALLMOLLIST_RAW

# Large default molecule list (59 species)
_LARGEMOLLIST_RAW = [
    1.1,   2.1,   6.1,   7.1,   8.1,   9.1,
   12.1,  13.1,  14.1,  15.1,  16.1,  17.1,
   20.1,  22.1,  23.1,  24.1,  26.1,
  101.0, 106.0, 107.0, 108.0, 109.0, 112.0,
  113.0, 114.0, 115.0, 116.0, 117.0, 120.0,
  124.0, 126.0, 10106.0, 10107.0, 10108.0, 10115.0,
  10116.0, 10608.0, 10812.0, 10813.0, 10820.0, 606.0,
  607.0, 608.0, 616.0, 60808.0, 707.0, 708.0,
  714.0, 715.0, 716.0, 808.0, 812.0, 814.0,
  815.0, 816.0, 822.0, 823.0, 826.0, 1416.0,
]
LARGEMOLLIST = np.zeros(110)
LARGEMOLLIST[:len(_LARGEMOLLIST_RAW)] = _LARGEMOLLIST_RAW

# H2O and CO2 partition function polynomial coefficients (HITRAN fits)
# u(T) = sum_{j=0}^{4} c[j] * T^j
H2OCOEFF = np.array([2.01999e+00,  1.36468e-03, -3.70508e-07,  6.74255e-11, -4.94970e-15])
CO2COEFF = np.array([2.00071e+00,  1.79927e-03, -4.06396e-07,  6.00864e-11, -3.84741e-15])

# Molecular equilibrium constant data (from Bmolec.f).
# Each row: [mol_id, D0, c1, c2, c3, c4, c5]
# Kp reconstruction: log10(Kp) = c1 + c2*lth + c3*lth^2 + c4*lth^3 + c5*lth^4 - D0*theta
# where lth = log10(theta), theta = 5040/T
_DATMOL_ROWS = [
  [101.0, 4.4781,12.1174,-1.0476, 1.6851,-5.5831, 4.0060],
  [106.0, 3.4650,11.5335,-0.5211,-0.7475, 0.1494,-0.1967],
  [107.0, 3.4700,11.4657,-0.7265,-0.6439, 0.0004, 0.1269],
  [108.0, 4.3920,11.8018,-0.8525,-0.5525, 0.1625,-0.1935],
  [109.0, 5.8690,12.2897,-0.9174,-0.6416, 0.1616,-0.1222],
  [111.0, 1.8800,10.7189,-0.8053, 4.1674,-12.8321,10.8590],
  [112.0, 1.3400,10.2878,-0.3455, 0.1677,-3.8628, 4.7348],
  [113.0, 3.0600,11.4876,-0.4024,-0.4809,-1.6283, 2.5415],
  [114.0, 3.0600,11.2586,-0.6758,-0.5870, 0.0669, 0.3139],
  [115.0, 3.3000,11.3387,-0.2112, 0.5964, 0.2027, 0.2323],
  [116.0, 3.5500,11.4380,-0.7731,-0.4785, 0.1716,-0.2326],
  [117.0, 4.4336,11.9042,-0.8250,-0.6309, 0.1545,-0.1999],
  [120.0, 1.7000,10.1987,-0.9426, 1.8085,-4.6629, 3.4971],
  [124.0, 2.1700,10.4501,-3.4047,-2.5032, 1.6933,-2.1073],
  [125.0, 1.3100, 9.7219,-3.9379,-3.4116, 0.6378,-3.0173],
  [126.0, 2.4100,12.1214, 0.9531, 2.3351,-0.2231, 3.0718],
  [128.0, 2.7000,11.9592,-0.9476,-0.4685, 0.8228, 0.2487],
  [129.0, 2.8400,11.3419,-1.3372,-0.6389, 1.5957,-0.4408],
  [10106.0, 7.9400,23.8688,-1.7944, 4.4565,-10.8615, 6.6375],
  [10107.0, 7.4400,23.7463,-1.7687, 4.2349,-12.2375, 8.6009],
  [10108.0, 9.6221,24.6063,-1.8370, 3.9590,-10.9331, 7.4896],
  [10115.0, 6.4895,23.0957,-2.0802, 5.0222,-10.7703, 5.8343],
  [10116.0, 7.5946,23.8619,-1.7009, 4.4792,-11.4475, 7.2724],
  [10508.0,12.7425,25.2365,-1.2673, 5.1472,-12.0672, 7.6325],
  [10607.0,13.2363,25.1400,-1.3548, 5.4650,-12.6262, 7.5607],
  [10608.0,11.8560,24.6494,-1.6665, 4.8174,-10.9079, 6.7311],
  [10708.0, 8.6140,24.4465,-1.3261,-0.5924, 0.0156,-0.7950],
  [10811.0, 8.0150,23.3475,-1.4238, 7.1752,-17.8186,13.4471],
  [10812.0, 8.0735,23.3316,-1.3691, 6.1267,-15.3175,11.1568],
  [10813.0,10.1252,25.2641,-1.4342, 5.1091,-13.2682, 9.2414],
  [10819.0, 8.1892,23.3235,-1.9306, 8.6260,-17.6081,10.5153],
  [10820.0, 8.7035,23.2006,-1.9644, 8.3448,-17.4669,10.8898],
  [10856.0, 9.0621,23.3508,-2.9133, 7.9864,-14.7824,10.4060],
  [508.0,  8.2800,12.6247,-0.6958,-0.4147, 0.2800,-0.4535],
  [606.0,  6.2100,12.4677,-0.4434,-0.0516,-0.1304,-0.0555],
  [607.0,  7.7600,12.4439,-0.4823,-0.4724,-1.1721, 1.3124],
  [608.0, 11.0920,13.2412,-0.8502,-0.0724,-0.2098,-0.2096],
  [614.0,  4.6400,11.8943,-1.2037, 2.5280,-5.4728, 3.5183],
  [615.0,  6.8950,13.0550, 1.1581, 2.9371, 0.6616, 1.9298],
  [616.0,  7.3550,12.8508,-0.7689,-0.4254, 1.9738,-2.7599],
  [60606.0,13.8610,26.3157,-1.9852, 4.7933,-10.1191, 6.4350],
  [60614.0,13.1966,25.7482,-0.8673, 6.1588,-14.8609,10.8245],
  [60717.0,12.2076,25.3546,-0.9498, 6.1809,-13.9272, 8.7748],
  [60808.0,16.5382,26.9665,-1.5802, 5.2319,-12.8507, 7.2799],
  [61616.0,11.9993,26.2469,-1.0515, 6.3688,-12.9191, 7.5525],
  [707.0,  9.7594,12.8868,-0.8861, 0.2644,-1.4001,  0.9645],
  [708.0,  6.4968,11.9347,-0.7631, 0.0848,-0.8149,  0.4636],
  [709.0,  2.8190,11.4441,-1.2536,-0.9028,-0.5407, -0.8453],
  [714.0,  4.5100,11.9190,-0.7171,-0.7889,-1.6147,  1.2567],
  [715.0,  7.1110,12.1020,-1.2837, 0.4612,-1.4609,  0.1464],
  [716.0,  4.8000,11.9264,-1.0597, 2.7095,-6.0355,  3.6858],
  [70708.0,11.4400,25.8228,-1.8022, 5.4142,-13.7572, 7.8651],
  [70808.0, 9.6210,25.4890,-2.1013, 4.7855,-12.0992, 6.7146],
  [808.0,  5.1156,12.8763,-0.4923,-0.5474, 0.2097, -0.3337],
  [811.0,  3.0790,11.1443, 0.1006, 1.3407,-0.8726,  1.1604],
  [812.0,  3.5300,10.7965,-0.4336, 4.5543,-9.6555,  6.8587],
  [813.0,  5.2700,12.2111,-0.5018,-0.0903,-1.3613,  1.7974],
  [814.0,  8.2600,12.9276,-0.7698,-0.5315, 2.1774, -2.8482],
  [815.0,  6.0710,11.9149,-1.0730, 0.4169,-0.8595, -0.0207],
  [816.0,  5.3590,12.3424,-0.8964, 2.5742,-6.2111,  3.7515],
  [817.0,  2.7450,11.8129,-1.0573, 2.2239,-5.7356,  3.2844],
  [820.0,  4.5310,11.7806, 1.6875, 3.1813,-5.1106,  4.8826],
  [821.0,  6.9600,12.5239,-1.2329, 1.5239,-1.9903,  0.4735],
  [822.0,  6.8700,12.3189,-1.8737, 4.2408,-7.4936,  3.3890],
  [823.0,  6.4100,12.8103,-0.5642,-0.5910,-2.5261,  3.9417],
  [826.0,  4.2000,12.5333,-1.0582, 1.0158,-1.5253,  0.7221],
  [839.0,  7.2900,12.4455,-1.3319, 1.0692, 0.0949, -2.0131],
  [840.0,  7.8500,12.4688,-1.0832, 0.0935,-0.2133,  0.5480],
  [856.0,  5.4410,11.5981,-2.2831,-1.3418, 4.2093, -4.0635],
  [857.0,  8.2300,12.1926, 0.1654,-0.8084,-1.0001,  1.1994],
  [80814.0,13.0355,26.5705,-1.1245, 6.0979,-12.6131, 7.8414],
  [80816.0,11.1405,25.9338,-1.3799, 5.5465,-11.7480, 7.0670],
  [80822.0,13.2915,25.9435,-2.1804, 7.1450,-13.0296, 7.9278],
  [80839.0,15.2000,25.8616,-1.4066,-0.5348, 2.4526, -1.0943],
  [80840.0,14.4650,25.6438,-2.2794, 6.4569,-11.9176, 7.3456],
  [80857.0,21.1510,31.0796,10.7083,13.0309, 9.1626,10.4251],
  [81313.0,10.9653,24.8867,-0.8364, 6.6875,-15.6083,11.6178],
  [909.0,  1.5920,12.6197,-0.4363,-0.5840, 0.4907, -0.6261],
  [911.0,  4.9530,11.4755,-0.4941, 0.6264,-1.1793,  0.5988],
  [912.0,  3.2000, 9.4953,-4.2045,-4.0612,-3.2308, -3.0520],
  [913.0,  6.8900,12.2405,-0.4662,-0.3499,-0.6568,  1.0680],
  [914.0,  5.5700,12.0156,-0.5029,-0.1059, 0.3137,  0.0030],
  [916.0,  3.3380,11.6730,-0.9068,-0.8140,-0.0380, -0.8173],
  [917.0,  2.6160,12.2059, 0.4295,-3.9612, 6.8690, -4.3459],
  [1111.0, 0.7300,10.2142,-0.5498, 1.0576,-2.8013,  1.4272],
  [1117.0, 4.2300,11.0824,-0.4354, 1.4015,-6.7895,  7.3949],
  [1216.0, 2.4000,11.3114,-0.5242, 0.5655,-3.1955,  3.8100],
  [1217.0, 2.7010,10.2253,-1.7707,-1.1903,-1.3832, -0.9560],
  [1313.0, 1.5500,11.3645,-0.1608, 0.1321,-3.2955,  4.3483],
  [1316.0, 3.8400,11.9118,-0.4891,-0.0108,-1.2859,  1.7739],
  [1317.0, 5.1200,11.8263,-0.3140,-1.0244, 0.8982, -0.1197],
  [1414.0, 3.2100,12.1806,-0.6731,-0.1732, 0.0349,  0.4326],
  [1416.0, 6.4200,12.6372,-0.7285,-0.0244, 0.8566, -1.1574],
  [1417.0, 4.0020,11.7232,-0.2951, 0.1824, 0.8534, -0.3405],
  [1515.0, 5.0330,12.2546,-1.1717, 1.4643,-1.9459,  0.8396],
  [1516.0, 5.6370,12.5482, 1.4333, 3.0823, 0.9032,  2.0434],
  [1616.0, 4.3693,12.3238,-0.9114, 2.7995,-6.1180,  3.5118],
  [1617.0, 2.7490,11.8411,-0.0001, 0.5709, 0.7753,  0.0580],
  [1622.0, 4.7500,11.6639,-1.3994, 1.3016,-1.2659,  0.5035],
  [1717.0, 2.4760,12.2664,-0.5351,-0.6211, 0.7092, -0.7675],
  [2016.0, 3.4600,10.8548, 0.2371, 3.1884,-8.1647,  5.1342],
  [2616.0, 3.1000,11.3327,-1.5519, 0.5150,-1.5570, -0.0107],
]
# Build DATMOL: shape (7, 110), DATMOL[0,:] = mol_ids, DATMOL[1:7,:] = constants
_DATMOL_ARR = np.array(_DATMOL_ROWS)                # shape (N, 7)
DATMOL = np.zeros((7, 110))
DATMOL[:, :len(_DATMOL_ROWS)] = _DATMOL_ARR.T       # shape (7, N) → (7, 110)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _sunder(amol: float):
    """Split molecule species code into leftmost atom and remainder (Sunder.f)."""
    im = round(amol)
    for divisor in (100_000_000, 1_000_000, 10_000, 100, 1):
        i1 = im // divisor
        if i1:
            i2 = im - i1 * divisor
            return i1, i2
    return 0, 0


def _discov(amol: float, i1: int) -> int:
    """Return the number of times atom i1 appears in molecule amol (Discov.f)."""
    im = round(amol)
    count = 0
    for divisor in (100_000_000, 1_000_000, 10_000, 100, 1):
        i3 = im // divisor
        if i3 == i1:
            count += 1
        im -= i3 * divisor
    return count


# ---------------------------------------------------------------------------
# Public initialization
# ---------------------------------------------------------------------------

def init_mol_data(state) -> None:
    """
    Copy Bmolec.f BLOCK DATA molecular constants into *state*.

    Must be called once before the first call to inmodel / eqlib.
    Populates: smallmollist, largemollist, datmol, h2ocoeff, co2coeff.
    """
    state.smallmollist[:] = SMALLMOLLIST
    state.largemollist[:] = LARGEMOLLIST
    state.datmol[:, :]    = DATMOL
    state.h2ocoeff[:]     = H2OCOEFF
    state.co2coeff[:]     = CO2COEFF


# ---------------------------------------------------------------------------
# Molecular equilibrium solver
# ---------------------------------------------------------------------------

def eqlib(state) -> None:
    """
    Solve molecular dissociation equilibrium for all depth layers.

    Translated from Eqlib.f + Setmols.f.  Uses Newton-Raphson iteration at
    each depth layer (innermost to outermost) until convergence.

    After this call the following State fields are populated:
        xmol[jmol, i]   — number density of each molecular/ionic species
        xamol[k, i]     — neutral atom number density for each equilibrium atom
        pmol[jmol]      — log10 partial pressure of each species (last layer)
        patom[k]        — log10 partial pressure of each atom (last layer)
        numdens[j,0,i]  — neutral densities of H, He, C, Mg, Al, Si, Fe
        numdens[j,1,i]  — ion densities of same 7 elements
        numdens[7,0,i]  — H2 number density
        xnh2o[i], xnco2[i] — H2O and CO2 densities
        uh2o[i], uco2[i]   — H2O and CO2 partition functions

    Parameters
    ----------
    state : State
        MOOG State dataclass.  Requires ntau, t, ne, nhtot, xabund, molopt.
    """
    # Fortran BLOCK DATA constants (analogous to compile-time initialisation)
    init_mol_data(state)

    ntau = state.ntau
    nmol = state.nmol

    # Temperature step for precomputing atomic partition function ratios
    tdel = (state.t[ntau - 1] - state.t[0]) / 4.0 + 1.0

    # ------------------------------------------------------------------ #
    # 1. Molecular/ionic equilibrium constants  (const[0..5, jmol])       #
    #    const[0]     = D0 (dissociation energy) or chi1 (ionization pot) #
    #    const[1..5]  = polynomial coefficients (molecules)                #
    #                or PF ratios at 5 T points (ions)                    #
    # ------------------------------------------------------------------ #
    const = np.zeros((6, 110))

    for jmol in range(nmol):
        amol_val = state.amol[jmol]
        if amol_val >= 100.0:
            # Molecular: look up datmol table by species code
            found = False
            for k in range(110):
                if state.datmol[0, k] == amol_val:
                    const[:, jmol] = state.datmol[1:7, k]
                    found = True
                    break
            if not found:
                raise ValueError(f"eqlib: unknown molecule {amol_val:.1f} "
                                 "not found in datmol")
        else:
            # Atomic/ionic: compute PF ratios at 5 representative temperatures
            iatom1 = int(amol_val + 0.0001)   # 1-based atomic number Z
            const[0, jmol] = XCHI1[iatom1 - 1]
            for kk in range(1, 6):
                ti = state.t[0] + (kk - 1) * tdel
                uu = np.zeros(2)
                for jj in range(1, 3):          # jj=1=neutral, jj=2=singly ionized
                    z_idx = iatom1 - 1           # 0-based for PARTFLAG
                    j0    = 4 * z_idx + (jj - 1) # 0-based species index for ucalc
                    if PARTFLAG[z_idx, jj - 1] > 0:
                        uu[jj - 1] = partnew(iatom1, jj, ti)
                    else:
                        uu[jj - 1] = ucalc(j0, ti)
                const[kk, jmol] = uu[1] / uu[0]  # ratio U_ion / U_neutral

    # ------------------------------------------------------------------ #
    # 2. Build iorder and ident arrays                                     #
    #    iorder[k]       = atomic number Z of k-th equilibrium element     #
    #    ident[k, 0..n]  = 1-indexed molecule numbers containing atom k    #
    #                      (0 = empty slot)                                #
    # ------------------------------------------------------------------ #
    iorder = np.zeros(30, dtype=int)
    ident  = np.zeros((30, 110), dtype=int)   # 1-indexed jmol values, 0=empty
    neq    = 0
    nmax   = 1    # current number of columns used in ident

    for jmol in range(nmol):
        atom = float(state.amol[jmol])
        while True:
            ia, ib = _sunder(atom)
            iatom1 = ia

            # Search for iatom1 in iorder[0..neq-1]
            k_found = -1
            for k in range(neq):
                if iorder[k] == iatom1:
                    k_found = k
                    break

            jmol1 = jmol + 1   # 1-indexed molecule number (0 = empty in ident)

            if k_found < 0:
                # New element: append to iorder, start its ident list
                k_found = neq
                iorder[neq] = iatom1
                ident[neq, 0] = jmol1
                neq += 1
            else:
                # Existing element: find empty slot or existing jmol1 in ident
                inserted = False
                for kk in range(nmax):
                    if ident[k_found, kk] == 0 or ident[k_found, kk] == jmol1:
                        ident[k_found, kk] = jmol1
                        inserted = True
                        break
                if not inserted:
                    ident[k_found, nmax] = jmol1
                    nmax += 1

            if ib == 0:
                break
            atom = float(ib)   # recurse on remaining atoms in molecule

    state.neq = neq
    state.iorder[:neq] = iorder[:neq]

    # ------------------------------------------------------------------ #
    # Precompute vectorization tables (outside depth/NR loops)            #
    # ------------------------------------------------------------------ #
    disc_table = np.array([[_discov(state.amol[jm], int(iorder[k]))
                             for k in range(neq)] for jm in range(nmol)], dtype=int)

    mol_mask    = np.array([state.amol[jm] >= 100.0 for jm in range(nmol)])
    mol_indices = np.where(mol_mask)[0]
    ion_indices = np.where(~mol_mask)[0]

    # True molecules: composition matrix, atom count, ionisation, equilibrium polynomial
    mol_comp   = disc_table[mol_indices, :]        # (n_mol, neq) atom counts
    mol_count  = mol_comp.sum(axis=1).astype(float)
    mol_hion   = np.array([10.0 * (state.amol[jm] - int(state.amol[jm]))
                            for jm in mol_indices])
    poly_coeff = const[1:6, mol_indices].T         # (n_mol, 5)
    d0_mol     = const[0, mol_indices]             # (n_mol,)

    # Ionic species: index into xatom, ionisation potential, PF ratio table
    ion_k_idx  = np.array([next(k for k in range(neq)
                                if iorder[k] == int(state.amol[jm]))
                            for jm in ion_indices], dtype=int)
    ion_chi    = const[0, ion_indices]             # (n_ion,) chi1 values
    ion_pf     = const[1:6, ion_indices]           # (5, n_ion) PF ratios at 5 T points

    # ------------------------------------------------------------------ #
    # 3. Main loop: iterate over depth layers from deep to shallow         #
    # ------------------------------------------------------------------ #
    xatom = np.zeros(30)   # neutral number densities of equilibrium atoms
    xfic  = np.zeros(30)   # total number density (upper bound) per atom
    xcorr = np.zeros(30)   # Newton-Raphson correction vector

    for kev in range(1, ntau + 1):
        i = ntau - kev     # Python 0-based (Fortran: i = ntau+1-kev, 1-based)

        tk = 1.38065e-16 * state.t[i]   # kT in ergs

        # Initial guess: xfic = total available number density
        for k in range(neq):
            korder = iorder[k]
            xfic[k] = state.xabund[korder - 1] * state.nhtot[i]

        # Propagate solution from previous (deeper) layer, or initialize
        if i < ntau - 1:
            for k in range(neq):
                xatom[k] = xatom[k] * state.nhtot[i] / state.nhtot[i + 1]
        else:
            xatom[:neq] = xfic[:neq]

        # ---- Newton-Raphson iteration ----
        while True:
            # True molecules (vectorized)
            if state.t[i] > 12000.0:
                state.xmol[mol_indices, i] = 1.0e-20
            else:
                th       = 5040.0 / state.t[i]
                lth      = np.log10(th)
                poly_vec = np.array([1.0, lth, lth**2, lth**3, lth**4])
                kp_mol   = poly_coeff @ poly_vec - d0_mol * th
                log_xa   = np.log(np.maximum(xatom[:neq], 1e-300))
                log_xm   = mol_comp @ log_xa                 # (n_mol,)
                xmol_mol = (np.exp(log_xm) * tk**(mol_count - 1.0)
                             / 10.0**kp_mol / state.ne[i]**mol_hion)
                state.xmol[mol_indices, i] = xmol_mol

            # Ionic species: Saha equation — must update every NR iteration (vectorized)
            if len(ion_indices) > 0:
                delt      = (state.t[i] - state.t[0]) / tdel
                m         = min(int(delt) + 1, 4)
                delt_frac = delt - int(delt)
                u1_ion    = ion_pf[m - 1] + (ion_pf[m] - ion_pf[m - 1]) * delt_frac
                xmol_ion  = (4.825e15 * u1_ion * state.t[i]**1.5 / state.ne[i]
                              * np.exp(-1.1605e4 * ion_chi / state.t[i])
                              * xatom[ion_k_idx])
                state.xmol[ion_indices, i] = xmol_ion

            # Jacobian and residual (vectorized)
            xmol_i   = state.xmol[:nmol, i]
            weighted = disc_table * xmol_i[:, None]          # (nmol, neq)
            deltax   = -xfic[:neq] + xatom[:neq] + weighted.sum(axis=0)
            xa_safe  = np.maximum(xatom[:neq], 1e-200)
            c_mat    = np.eye(neq) + (disc_table.T @ weighted) / xa_safe[None, :]

            # Invert and compute corrections
            try:
                c_inv = np.linalg.inv(c_mat)
            except np.linalg.LinAlgError:
                break

            # First loop: compute xcorr, saving previous xcorr in x1
            x1 = 0.0
            for k in range(neq):
                x1       = xcorr[k]     # save old xcorr[k]
                xcorr[k] = float(c_inv[k, :] @ deltax)

            # Second loop: check convergence, damp, apply correction
            iflag = False
            for k in range(neq):
                # Oscillation damping using x1 from previous k (Fortran logic)
                if x1 * xcorr[k] < -0.5 * x1 ** 2:
                    xcorr[k] *= 0.5
                x1 = xatom[k]   # update x1 for next k's damping check

                if abs(xcorr[k] / xatom[k]) > 0.005:
                    iflag = True

                x_before  = xatom[k]
                xatom[k] -= xcorr[k]

                if xatom[k] <= 0.0 or xatom[k] >= 1.001 * xfic[k]:
                    iflag = True
                    xatom[k] = 1.0e-2 * abs(x_before)

            if not iflag:
                break
        # ---- end Newton-Raphson ----

        # Store converged number densities and log partial pressures
        for k in range(neq):
            state.xamol[k, i] = xatom[k]
            state.patom[k]    = np.log10(xatom[k] * tk)
        for jmol in range(nmol):
            state.pmol[jmol]  = np.log10(state.xmol[jmol, i] * tk)

    # ------------------------------------------------------------------ #
    # 4. Post-processing: setmols (Setmols.f)                             #
    # ------------------------------------------------------------------ #
    _setmols(state, iorder, neq)


def _setmols(state, iorder: np.ndarray, neq: int) -> None:
    """
    Transfer molecular equilibrium output to arrays needed elsewhere
    (opacity calculations, etc.).  Translated from Setmols.f.
    """
    ntau = state.ntau
    nmol = state.nmol

    # Elements of interest: H(1), He(2), C(6), Mg(12), Al(13), Si(14), Fe(26)
    nel = [1, 2, 6, 12, 13, 14, 26]

    # Neutral number densities
    for j, z in enumerate(nel):
        for k in range(neq):
            if iorder[k] == z:
                state.numdens[j, 0, :ntau] = state.xamol[k, :ntau]
                break

    # Ion number densities (species code = Z.1 → 10*Z+1 as integer)
    for j, z in enumerate(nel):
        ispec10 = 10 * z + 1
        for k in range(nmol):
            kmol10 = round(10.0 * state.amol[k])
            if ispec10 == kmol10:
                state.numdens[j, 1, :ntau] = state.xmol[k, :ntau]
                break

    # H2 number density (species code 101)
    for k in range(nmol):
        if round(state.amol[k]) == 101:
            state.numdens[7, 0, :ntau] = state.xmol[k, :ntau]
            break

    # H2O and CO2 partition functions
    for i in range(ntau):
        if state.t[i] > 5000.0:
            state.uh2o[i] = 1.0e8
            state.uco2[i] = 1.0e8
        else:
            h2olog = sum(state.h2ocoeff[j] * state.t[i] ** j for j in range(5))
            co2log = sum(state.co2coeff[j] * state.t[i] ** j for j in range(5))
            state.uh2o[i] = 10.0 ** h2olog
            state.uco2[i] = 10.0 ** co2log

    # Transfer H2O and CO2 number densities
    ih2o = -1
    ico2 = -1
    for j in range(nmol):
        if round(state.amol[j]) == 10108:
            ih2o = j
        if round(state.amol[j]) == 60808:
            ico2 = j

    if ih2o >= 0:
        state.xnh2o[:ntau] = state.xmol[ih2o, :ntau]
    if ico2 >= 0:
        state.xnco2[:ntau] = state.xmol[ico2, :ntau]
