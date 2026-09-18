"""
H I bound-free/free-free and H- bound-free/free-free continuous opacity.
Translated from OpacHydrogen.f and Opaccouls.f (ATLAS9 routines).

Bug fix: Fortran opacH1 references undeclared `cont8` where `cont(8)` was
intended. Python uses cont[7] (physically correct n=8 Coulomb cross-section).
"""
import numpy as np

# ---------------------------------------------------------------------------
# coulbf1s: 1s Gaunt factor table (151 points, log10 step 0.02)
# ---------------------------------------------------------------------------
_GAUNT1S = np.array([
    0.7973,0.8094,0.8212,0.8328,0.8439,0.8548,0.8653,0.8754,0.8852,
    0.8946,0.9035,0.9120,0.9201,0.9278,0.9351,0.9420,0.9484,0.9544,
    0.9601,0.9653,0.9702,0.9745,0.9785,0.9820,0.9852,0.9879,0.9903,
    0.9922,0.9938,0.9949,0.9957,0.9960,0.9960,0.9957,0.9949,0.9938,
    0.9923,0.9905,0.9884,0.9859,0.9832,0.9801,0.9767,0.9730,0.9688,
    0.9645,0.9598,0.9550,0.9499,0.9445,0.9389,0.9330,0.9269,0.9206,
    0.9140,0.9071,0.9001,0.8930,0.8856,0.8781,0.8705,0.8627,0.8546,
    0.8464,0.8381,0.8298,0.8213,0.8128,0.8042,0.7954,0.7866,0.7777,
    0.7685,0.7593,0.7502,0.7410,0.7318,0.7226,0.7134,0.7042,0.6951,
    0.6859,0.6767,0.6675,0.6584,0.6492,0.6401,0.6310,0.6219,0.6129,
    0.6039,0.5948,0.5859,0.5769,0.5680,0.5590,0.5502,0.5413,0.5324,
    0.5236,0.5148,0.5063,0.4979,0.4896,0.4814,0.4733,0.4652,0.4572,
    0.4493,0.4415,0.4337,0.4261,0.4185,0.4110,0.4035,0.3962,0.3889,
    0.3818,0.3749,0.3680,0.3611,0.3544,0.3478,0.3413,0.3348,0.3285,
    0.3222,0.3160,0.3099,0.3039,0.2980,0.2923,0.2866,0.2810,0.2755,
    0.2701,0.2648,0.2595,0.2544,0.2493,0.2443,0.2394,0.2345,0.2298,
    0.2251,0.2205,0.2160,0.2115,0.2072,0.2029,0.1987,
])  # 151 entries; index 0 corresponds to log10(freq/z²/3.28805e15)=0

# coulx polynomial coefficients for Gaunt factors at n=2..6
_COULX_A = np.array([0.9916, 1.105, 1.101, 1.101, 1.102, 1.0986])
_COULX_B = np.array([2.719e13, -2.375e14, -9.863e13, -5.765e13, -3.909e13, -2.704e13])
_COULX_C = np.array([-2.268e30, 4.077e28, 1.035e28, 4.593e27, 2.371e27, 1.229e27])

# coulff bilinear table a(igam, ihvkt), shape (11, 12), Fortran column-major
# Each block of 11 values fills one column (ihvkt dimension)
_COULFF_A = np.array([
    [5.53,5.49,5.46,5.43,5.40,5.25,5.00,4.69,4.48,4.16,3.85],
    [4.91,4.87,4.84,4.80,4.77,4.63,4.40,4.13,3.87,3.52,3.27],
    [4.29,4.25,4.22,4.18,4.15,4.02,3.80,3.57,3.27,2.98,2.70],
    [3.64,3.61,3.59,3.56,3.54,3.41,3.22,2.97,2.70,2.45,2.20],
    [3.00,2.98,2.97,2.95,2.94,2.81,2.65,2.44,2.21,2.01,1.81],
    [2.41,2.41,2.41,2.41,2.41,2.32,2.19,2.02,1.84,1.67,1.50],
    [1.87,1.89,1.91,1.93,1.95,1.90,1.80,1.68,1.52,1.41,1.30],
    [1.33,1.39,1.44,1.49,1.55,1.56,1.51,1.42,1.33,1.25,1.17],
    [0.90,0.95,1.00,1.08,1.17,1.30,1.32,1.30,1.20,1.15,1.11],
    [0.55,0.58,0.62,0.70,0.85,1.01,1.15,1.18,1.15,1.11,1.08],
    [0.33,0.36,0.39,0.46,0.59,0.76,0.97,1.09,1.13,1.10,1.08],
    [0.19,0.21,0.24,0.28,0.38,0.53,0.76,0.96,1.08,1.09,1.09],
]).T   # shape (11, 12): row=igam (0..10), col=ihvkt (0..11)

_COULFF_Z4LOG = np.array([0.0, 1.20412, 1.90849, 2.40824, 2.79588, 3.11261])

# ---------------------------------------------------------------------------
# H- opacity tables (Bell & Berrington 1987; Wishart 1979; Broad & Reinhardt)
# ---------------------------------------------------------------------------
_WBF = np.array([
     18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,  40.40,
     49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05, 113.10,
    113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00, 175.00,
    200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00, 375.00,
    400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00, 575.00,
    600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00, 775.00,
    800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00, 975.00,
   1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,1175.00,
   1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,1375.00,
   1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,1575.00,
   1600.00,1610.00,1620.00,1630.00,1643.91,
])  # 85 wavelengths [nm]

_BF = np.array([
     0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,  0.703,
      1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,  56.60,
     20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,  7.918,
     9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,  21.35,
     23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,  33.17,
     34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,  39.48,
     39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,  38.53,
     37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,  30.34,
     28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,  16.74,
     14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,  2.575,
      1.302, 0.8697, 0.4974, 0.1989,    0.0,
])  # 85 cross-sections [10^-18 cm^2]

_WAVEK = np.array([.50,.40,.35,.30,.25,.20,.18,.16,.14,.12,.10,
                   .09,.08,.07,.06,.05,.04,.03,.02,.01,.008,.006])  # 22

_THETAFF = np.array([0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6])  # 11

# ff(itheta, iwave) Fortran column-major → each group of 11 fills one iwave column
# Stored here row=iwave (0..21), col=itheta (0..10); then transposed → shape (11,22)
_FF_RAW = np.array([
    [.0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896,.131,.172],   # iwave 1
    [.0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972,.110,.160,.211],    # iwave 2
    [.0277,.0342,.0476,.0615,.0760,.0908,.105,.121,.136,.199,.262],      # iwave 3
    [.0364,.0447,.0616,.0789,.0966,.114,.132,.150,.169,.243,.318],       # iwave 4
    [.0520,.0633,.0859,.108,.131,.154,.178,.201,.225,.321,.418],         # iwave 5
    [.0791,.0959,.129,.161,.194,.227,.260,.293,.327,.463,.602],          # iwave 6
    [.0965,.117,.157,.195,.234,.272,.311,.351,.390,.549,.711],           # iwave 7
    [.121,.146,.195,.241,.288,.334,.381,.428,.475,.667,.861],            # iwave 8
    [.154,.188,.249,.309,.367,.424,.482,.539,.597,.830,1.07],            # iwave 9
    [.208,.250,.332,.409,.484,.557,.630,.702,.774,1.06,1.36],            # iwave 10
    [.293,.354,.468,.576,.677,.777,.874,.969,1.06,1.45,1.83],           # iwave 11
    [.358,.432,.572,.702,.825,.943,1.06,1.17,1.28,1.73,2.17],          # iwave 12
    [.448,.539,.711,.871,1.02,1.16,1.29,1.43,1.57,2.09,2.60],         # iwave 13
    [.579,.699,.924,1.13,1.33,1.51,1.69,1.86,2.02,2.67,3.31],         # iwave 14
    [.781,.940,1.24,1.52,1.78,2.02,2.26,2.48,2.69,3.52,4.31],         # iwave 15
    [1.11,1.34,1.77,2.17,2.53,2.87,3.20,3.51,3.80,4.92,5.97],         # iwave 16
    [1.73,2.08,2.74,3.37,3.90,4.50,5.01,5.50,5.95,7.59,9.06],         # iwave 17
    [3.04,3.65,4.80,5.86,6.86,7.79,8.67,9.50,10.3,13.2,15.6],         # iwave 18
    [6.79,8.16,10.7,13.1,15.3,17.4,19.4,21.2,23.0,29.5,35.0],         # iwave 19
    [27.0,32.4,42.6,51.9,60.7,68.9,76.8,84.2,91.4,117.,140.],          # iwave 20
    [42.3,50.6,66.4,80.8,94.5,107.,120.,131.,142.,183.,219.],           # iwave 21
    [75.1,90.0,118.,144.,168.,191.,212.,234.,253.,325.,388.],            # iwave 22
])
# _FF[itheta, iwave], shape (11, 22)
_FF = _FF_RAW.T

# Precomputed log arrays (computed once at import)
_WFFLOG = np.log(91.134 / _WAVEK)                              # shape (22,)
_FFLOG  = np.log(_FF.T * 1.0e-26)                             # shape (22, 11): [iwave, itheta]


# ---------------------------------------------------------------------------
# Helper: scalar linear interpolation (linter from ATLAS, nnew=1 case)
# ---------------------------------------------------------------------------
def _linter(xold, yold, xnew):
    """Linear interpolation at scalar xnew; extrapolates at both ends."""
    n = len(xold)
    i = 1   # 0-indexed; Fortran starts at iold=2 (1-indexed)
    while True:
        if xnew < xold[i] or i == n - 1:
            return (yold[i-1] + (yold[i]-yold[i-1]) /
                    (xold[i]-xold[i-1]) * (xnew-xold[i-1]))
        i += 1


# ---------------------------------------------------------------------------
# Helper: ATLAS quadratic interpolation (map1 from ATLAS, nnew=1 case)
# ---------------------------------------------------------------------------
def _map1_scalar(xold, fold, xnew):
    """ATLAS quadratic interpolation at scalar xnew (faithful translation of map1)."""
    nold = len(xold)
    # Find l (Fortran 1-indexed): smallest l>=2 with xnew < xold[l-1], or l>nold
    l = 2
    while xnew >= xold[l-1]:
        l += 1
        if l > nold:
            break

    # Linear cases: first two intervals or beyond-end extrapolation
    if l <= 3 or l > nold:
        l = min(l, nold)
        b = (fold[l-1]-fold[l-2]) / (xold[l-1]-xold[l-2])
        a = fold[l-1] - xold[l-1]*b
        return a + b*xnew

    # Quadratic: backward parabola through (l-2, l-1, l)
    l1, l2 = l-1, l-2
    d = (fold[l1-1]-fold[l2-1]) / (xold[l1-1]-xold[l2-1])
    cbac = (fold[l-1]/((xold[l-1]-xold[l1-1])*(xold[l-1]-xold[l2-1])) +
            (fold[l2-1]/(xold[l-1]-xold[l2-1]) -
             fold[l1-1]/(xold[l-1]-xold[l1-1])) / (xold[l1-1]-xold[l2-1]))
    bbac = d - (xold[l1-1]+xold[l2-1])*cbac
    abac = fold[l2-1] - xold[l2-1]*d + xold[l1-1]*xold[l2-1]*cbac

    if l >= nold:
        return abac + (bbac + cbac*xnew)*xnew

    # Blend with forward parabola through (l-1, l, l+1)
    d = (fold[l-1]-fold[l1-1]) / (xold[l-1]-xold[l1-1])
    cfor = (fold[l]/((xold[l]-xold[l-1])*(xold[l]-xold[l1-1])) +
            (fold[l1-1]/(xold[l]-xold[l1-1]) -
             fold[l-1]/(xold[l]-xold[l-1])) / (xold[l-1]-xold[l1-1]))
    bfor = d - (xold[l]+xold[l1-1])*cfor
    afor = fold[l1-1] - xold[l1-1]*d + xold[l]*xold[l1-1]*cfor
    wt = abs(cfor)/(abs(cfor)+abs(cbac)) if cfor != 0.0 else 0.0
    a = afor + wt*(abac-afor)
    b = bfor + wt*(bbac-bfor)
    c = cfor + wt*(cbac-cfor)
    return a + (b + c*xnew)*xnew


# ---------------------------------------------------------------------------
# coulbf1s: 1s bound-free Gaunt factor via tabulated interpolation
# ---------------------------------------------------------------------------
def _coulbf1s(freq, z):
    if freq < z*z*3.28805e15:
        return 0.0
    elog = np.log10(freq/(z*z)/3.28805e15)
    i = int(elog/0.02)
    i = max(0, min(149, i))   # Fortran: max(min(i+1,150),1) → 0-indexed [0..149]
    return _GAUNT1S[i] + (_GAUNT1S[i+1]-_GAUNT1S[i])/0.02 * (elog - i*0.02)


# ---------------------------------------------------------------------------
# coulx: Coulomb bound-free cross-section for level n
# ---------------------------------------------------------------------------
def _coulx(n, freq, z):
    if freq < z*z*3.28805e15/float(n*n):
        return 0.0
    result = 2.815e29 / freq**3 / float(n**5) * z**4
    if n > 6:
        return result
    ni = n - 1  # 0-indexed into polynomial arrays
    if n == 1:
        return result * _coulbf1s(freq, z)
    zf = z*z/freq
    return result * (_COULX_A[ni] + (_COULX_B[ni] + _COULX_C[ni]*zf)*zf)


# ---------------------------------------------------------------------------
# coulff: free-free Gaunt factor via bilinear table interpolation
# (vectorized over depth: tlog_arr is shape (ntau,), returns shape (ntau,))
# ---------------------------------------------------------------------------
def _coulff_vec(nz, tlog_arr, freq):
    gamlog  = 10.39638 - tlog_arr/1.15129 + _COULFF_Z4LOG[nz-1]
    hvktlg  = (np.log(freq) - tlog_arr)/1.15129 - 20.63764
    igam_f  = np.clip(np.floor(gamlog  + 7.0).astype(int), 1, 10)
    ihvkt_f = np.clip(np.floor(hvktlg  + 9.0).astype(int), 1, 11)
    p = gamlog  - (igam_f  - 7).astype(float)
    q = hvktlg  - (ihvkt_f - 9).astype(float)
    ia, ib = igam_f - 1, ihvkt_f - 1     # 0-indexed
    return ((1-p)*((1-q)*_COULFF_A[ia,ib]   + q*_COULFF_A[ia,ib+1]) +
              p *((1-q)*_COULFF_A[ia+1,ib]  + q*_COULFF_A[ia+1,ib+1]))


# ---------------------------------------------------------------------------
# opac_h1: H I bound-free + free-free opacity → state.aH1
# ---------------------------------------------------------------------------
def opac_h1(state) -> None:
    """H I b-f (Lyman through Paschen series, n=1..8) + f-f opacity."""
    ntau = state.ntau
    freq = state.freq

    nh1   = state.numdens[0, 0, :ntau]   # H neutral
    uh1   = state.u[0, 0, :ntau]
    nh2   = state.numdens[0, 1, :ntau]   # H+
    uh2   = state.u[0, 1, :ntau]
    tkev  = state.tkev[:ntau]
    t     = state.t[:ntau]
    ne    = state.ne[:ntau]
    tlog  = state.tlog[:ntau]
    evhkt = state.evhkt[:ntau]

    # Boltzmann populations for levels n=1..8 relative to ground
    nh1_u = nh1 / uh1
    bolt = np.empty((8, ntau))
    for n in range(1, 9):
        xn2 = float(n*n)
        bolt[n-1] = np.exp(-13.595*(1.0 - 1.0/xn2)/tkev) * 2.0*xn2 * nh1_u

    freet  = ne * nh2/uh2 / np.sqrt(t)
    xr     = nh1_u * (1.0/13.595) * tkev   # = nh1/u * (2/2/13.595) * tkev
    boltex = np.exp(-13.427/tkev) * xr      # level n=9 population proxy
    exlim  = np.exp(-13.595/tkev) * xr      # ionization limit population

    # Frequency-dependent Coulomb cross-sections for n=1..8
    cont = np.array([_coulx(n, freq, 1.0) for n in range(1, 9)])

    freq3    = freq**3
    cfree    = 3.6919e8  / freq3
    c_const  = 2.815e29  / freq3

    # Vectorized coulff over depth
    cff = _coulff_vec(1, tlog, freq)

    # High-n (>=9) population: switch at 4.05933e13 Hz (n=9 series limit)
    if freq < 4.05933e13:
        ex_arr = exlim / evhkt
    else:
        ex_arr = boltex

    stim = 1.0 - evhkt

    # n=7 and n=8 terms + high-n continuum + f-f (bug fix: cont[7] not cont8=0)
    h = ((cont[6]*bolt[6] + cont[7]*bolt[7] +
          (ex_arr - exlim)*c_const +
          cff*freet*cfree) * stim)

    # n=1..6 explicit levels
    for n in range(6):
        h += cont[n]*bolt[n]*stim

    state.aH1[:ntau] = h


# ---------------------------------------------------------------------------
# opac_hminus: H- bound-free + free-free opacity → state.aHminus
# ---------------------------------------------------------------------------
def opac_hminus(state) -> None:
    """H- b-f (Wishart/Broad&Reinhardt table) + f-f (Bell&Berrington table) opacity."""
    ntau = state.ntau
    freq = state.freq

    nh1   = state.numdens[0, 0, :ntau]
    uh1   = state.u[0, 0, :ntau]
    ne    = state.ne[:ntau]
    tkev  = state.tkev[:ntau]
    t     = state.t[:ntau]
    evhkt = state.evhkt[:ntau]

    # H- equilibrium number density (Saha-like; .754209 eV = electron affinity of H-)
    xhmin = (np.exp(0.754209/tkev) / (2.0*2.4148e15*t**1.5) * nh1/uh1 * ne)

    # wave in nm (c in nm/s = 2.99792458e17 nm/s)
    wave    = 2.99792458e17 / freq
    wavelog = np.log(wave)

    # Interpolate f-f table at this wavelength for each theta value
    # _FFLOG[:, itheta] is the log cross-section vs log(wave) curve for that theta
    fftt = np.empty(11)
    for itheta in range(11):
        fftlog = _linter(_WFFLOG, _FFLOG[:, itheta], wavelog)
        fftt[itheta] = (np.exp(fftlog) / _THETAFF[itheta] * 5040.0 * 1.380658e-16)

    # H- bound-free cross-section at this frequency (b-f threshold: 1643.9 nm)
    hminbf = 0.0
    if freq > 1.82365e14:
        hminbf = _map1_scalar(_WBF, _BF, wave)

    # Per-depth interpolation in theta and final opacity
    theta   = 5040.0 / t
    nh1_u   = nh1 / uh1

    for i in range(ntau):
        fftheta = _linter(_THETAFF, fftt, theta[i])
        hminff  = fftheta * nh1_u[i] * 2.0 * ne[i]
        h       = hminbf * 1.0e-18 * (1.0-evhkt[i]) * xhmin[i]
        state.aHminus[i] = h + hminff
