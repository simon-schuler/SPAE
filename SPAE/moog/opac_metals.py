"""
Metal bound-free continuous opacity: C I, Mg I/II, Al I, Si I/II, Fe I.
Translated from Opacmetals.f.
"""
import numpy as np

# ---------------------------------------------------------------------------
# Mg I Peach (1970) cross-section table
# peach[it, n] — it = temperature index (0..6, T=4000..10000K),
#               n = frequency segment pair index (0..14, 15 columns)
# ---------------------------------------------------------------------------
_MG1_PEACH = np.array([
    [-42.474,-42.350,-42.109,-41.795,-41.467,-41.159,-40.883],
    [-41.808,-41.735,-41.582,-41.363,-41.115,-40.866,-40.631],
    [-41.273,-41.223,-41.114,-40.951,-40.755,-40.549,-40.347],
    [-45.583,-44.008,-42.957,-42.205,-41.639,-41.198,-40.841],
    [-44.324,-42.747,-41.694,-40.939,-40.370,-39.925,-39.566],
    [-50.969,-48.388,-46.630,-45.344,-44.355,-43.568,-42.924],
    [-50.633,-48.026,-46.220,-44.859,-43.803,-42.957,-42.264],
    [-53.028,-49.643,-47.367,-45.729,-44.491,-43.520,-42.736],
    [-51.785,-48.352,-46.050,-44.393,-43.140,-42.157,-41.363],
    [-52.285,-48.797,-46.453,-44.765,-43.486,-42.480,-41.668],
    [-52.028,-48.540,-46.196,-44.507,-43.227,-42.222,-41.408],
    [-52.384,-48.876,-46.513,-44.806,-43.509,-42.488,-41.660],
    [-52.363,-48.856,-46.493,-44.786,-43.489,-42.467,-41.639],
    [-54.704,-50.772,-48.107,-46.176,-44.707,-43.549,-42.611],
    [-54.359,-50.349,-47.643,-45.685,-44.198,-43.027,-42.418],
]).T     # shape (7, 15): [temperature_idx, freq_segment_col]

_MG1_FREQMG = np.array([1.9341452e15, 1.8488510e15, 1.1925797e15,
                          7.9804046e14, 4.5772110e14, 4.1440977e14,
                          4.1113514e14])
_MG1_FLOG   = np.array([35.23123, 35.19844, 35.15334, 34.71490, 34.31318,
                          33.75728, 33.65788, 33.64994, 33.43947])
_MG1_TLG    = np.array([8.29405, 8.51719, 8.69951, 8.85367, 8.98720,
                          9.10498, 9.21034])

# ---------------------------------------------------------------------------
# Si I Peach (1970) cross-section table  shape (9, 19)
# ---------------------------------------------------------------------------
_SI1_PEACH = np.array([
    [38.136,38.138,38.140,38.141,38.143,38.144,38.144,38.145,38.145],
    [37.834,37.839,37.843,37.847,37.850,37.853,37.855,37.857,37.858],
    [37.898,37.898,37.897,37.897,37.897,37.896,37.895,37.895,37.894],
    [40.737,40.319,40.047,39.855,39.714,39.604,39.517,39.445,39.385],
    [40.581,40.164,39.893,39.702,39.561,39.452,39.366,39.295,39.235],
    [45.521,44.456,43.753,43.254,42.878,42.580,42.332,42.119,41.930],
    [45.520,44.455,43.752,43.251,42.871,42.569,42.315,42.094,41.896],
    [55.068,51.783,49.553,47.942,46.723,45.768,44.997,44.360,43.823],
    [53.868,50.369,48.031,46.355,45.092,44.104,43.308,42.652,42.100],
    [54.133,50.597,48.233,46.539,45.261,44.262,43.456,42.790,42.230],
    [54.051,50.514,48.150,46.454,45.176,44.175,43.368,42.702,42.141],
    [54.442,50.854,48.455,46.733,45.433,44.415,43.592,42.912,42.340],
    [54.320,50.722,48.313,46.583,45.277,44.251,43.423,42.738,42.160],
    [55.691,51.965,49.444,47.615,46.221,45.119,44.223,43.478,42.848],
    [55.661,51.933,49.412,47.582,46.188,45.085,44.189,43.445,42.813],
    [55.973,52.193,49.630,47.769,46.349,45.226,44.314,43.555,42.913],
    [55.922,52.141,49.577,47.715,46.295,45.172,44.259,43.500,42.858],
    [56.828,52.821,50.110,48.146,46.654,45.477,44.522,43.730,43.061],
    [56.657,52.653,49.944,47.983,46.491,45.315,44.360,43.569,42.901],
]).T     # shape (9, 19): [temperature_idx, freq_segment_col]

_SI1_FREQSI = np.array([2.1413750e15, 1.9723165e15, 1.7879689e15,
                          1.5152920e15, 5.5723927e14, 5.3295914e14,
                          4.7886458e14, 4.7216422e14, 4.6185133e14])
_SI1_FLOG   = np.array([35.45438, 35.30022, 35.21799, 35.11986, 34.95438,
                          33.95402, 33.90947, 33.80244, 33.78835,
                          33.76626, 33.70518])
_SI1_TLG    = np.array([8.29405, 8.51719, 8.69951, 8.85367, 8.98720,
                          9.10498, 9.21034, 9.30565, 9.39266])

# ---------------------------------------------------------------------------
# Si II Peach (1970) cross-section table  shape (6, 14)
# ---------------------------------------------------------------------------
_SI2_PEACH = np.array([
    [-43.8941,-43.8941,-43.8941,-43.8941,-43.8941,-43.8941],
    [-42.2444,-42.2444,-42.2444,-42.2444,-42.2444,-42.2444],
    [-40.6054,-40.6054,-40.6054,-40.6054,-40.6054,-40.6054],
    [-54.2389,-52.2906,-50.8799,-49.8033,-48.9485,-48.2490],
    [-50.4108,-48.4892,-47.1090,-46.0672,-45.2510,-44.5933],
    [-52.0936,-50.0741,-48.5999,-47.4676,-46.5649,-45.8246],
    [-51.9548,-49.9371,-48.4647,-47.3340,-46.4333,-45.6947],
    [-54.2407,-51.7319,-49.9178,-48.5395,-47.4529,-46.5709],
    [-52.7355,-50.2218,-48.4059,-47.0267,-45.9402,-45.0592],
    [-53.5387,-50.9189,-49.0200,-47.5750,-46.4341,-45.5082],
    [-53.2417,-50.6234,-48.7252,-47.2810,-46.1410,-45.2153],
    [-53.5097,-50.8535,-48.9263,-47.4586,-46.2994,-45.3581],
    [-54.0561,-51.2365,-49.1980,-47.6497,-46.4302,-45.4414],
    [-53.8469,-51.0256,-48.9860,-47.4368,-46.2162,-45.2266],
]).T     # shape (6, 14): [temperature_idx, freq_segment_col]

_SI2_FREQSI = np.array([4.9965417e15, 3.9466738e15, 1.5736321e15,
                          1.5171539e15, 9.2378947e14, 8.3825004e14,
                          7.6869872e14])
_SI2_FLOG   = np.array([36.32984, 36.14752, 35.91165, 34.99216, 34.95561,
                          34.45951, 34.36234, 34.27572, 34.20161])
_SI2_TLG    = np.array([9.21034, 9.39266, 9.54681, 9.68034, 9.79813, 9.90349])

# ---------------------------------------------------------------------------
# Fe I level data for b-f opacity
# ---------------------------------------------------------------------------
_FE1_GG  = np.array([25.,35.,21.,15., 9.,35.,33.,21.,27.,49., 9.,21.,27., 9., 9.,
                       25.,33.,15.,35., 3., 5.,11.,15.,13.,15., 9.,21.,15.,21.,25.,35.,
                        9., 5.,45.,27.,21.,15.,21.,15.,25.,21.,35., 5.,15.,45.,35.,55.,25.])
_FE1_EE  = np.array([500.,7500.,12500.,17500.,19000.,19500.,19500.,21000.,
                       22000.,23000.,23000.,24000.,24000.,24500.,24500.,26000.,26500.,
                       26500.,27000.,27500.,28500.,29000.,29500.,29500.,29500.,30000.,
                       31500.,31500.,33500.,33500.,34000.,34500.,34500.,35000.,35500.,
                       37000.,37000.,37000.,38500.,40000.,40000.,41000.,41000.,43000.,
                       43000.,43000.,43000.,44000.])
_FE1_WNO = np.array([63500.,58500.,53500.,59500.,45000.,44500.,44500.,43000.,
                       58000.,41000.,54000.,40000.,40000.,57500.,55500.,38000.,57500.,
                       57500.,37000.,54500.,53500.,55000.,34500.,34500.,34500.,34000.,
                       32500.,32500.,32500.,32500.,32000.,29500.,29500.,31000.,30500.,
                       29000.,27000.,54000.,27500.,24000.,47000.,23000.,44000.,42000.,
                       42000.,21000.,42000.,42000.])

_H = 6.6256e-27    # Planck constant [erg·s]
_K = 1.38065e-16   # Boltzmann constant [erg/K]


def _seaton(freq, freq0, xsect, power, a):
    """Seaton (1958) hydrogenic b-f cross-section approximation."""
    freqratio = freq0 / freq
    n = round(2.0 * power + 0.01)
    return xsect * (a + freqratio * (1.0 - a)) * np.sqrt(freqratio ** n)


def _peach_interp(freq, freqlg, freq_edges, flog_grid, tlg_grid, peach):
    """
    Interpolate the Peach (1970) table at (freq, T) for each depth layer.

    freq_edges : sorted list of edge frequencies (decreasing); we find which
                 frequency interval contains the current freq.
    flog_grid  : log-frequency grid for the peach table columns.
    tlg_grid   : log(T) grid for the peach table rows.
    peach      : 2-D array, shape (n_T, n_col).

    Returns xx[i] = interpolated log cross-section at each temperature in tlg_grid,
    as a (n_T,) array.
    """
    n_edge = len(freq_edges)
    # Find interval n such that freq > freq_edges[n]
    n = n_edge   # default: beyond all edges
    for k in range(n_edge):
        if freq > freq_edges[k]:
            n = k
            break

    # Map interval index to Peach table column (same as Fortran: n > 2 → n = 2*n - 2)
    n_flog = len(flog_grid)
    dd  = (freqlg - flog_grid[n]) / (flog_grid[n + 1] - flog_grid[n])
    n_col = n if n <= 1 else 2 * n - 2   # Fortran mapping
    dd1 = 1.0 - dd
    xx = peach[:, n_col + 1] * dd + peach[:, n_col] * dd1
    return xx


def opac_c1(state) -> None:
    """C I bound-free opacity (Luo & Pradhan 1989; Burke & Taylor 1979)."""
    ntau = state.ntau
    freq = state.freq

    # Frequency-independent cross-sections (depend on which edge we're above)
    ryd    = 109732.298
    waveno = freq / 2.99792458e10   # [cm⁻¹]
    xs0 = xs1 = xd0 = xd1 = xd2 = x1444 = x1240 = x1100 = 0.0

    if freq >= 2.7254e15:
        x1100 = (10.0 ** (-16.80 - (waveno - 90777.0) / 3.0 / ryd)
                 * _seaton(freq, 2.7254e15, 1.219e-17, 2.0, 3.317))

    if freq >= 2.4196e15:
        xd0 = 10.0 ** (-16.80 - (waveno - 80627.760) / 3.0 / ryd)
        eeps = (waveno - 93917.0) * 2.0 / 9230.0
        xd1 = (22.0e-18 * eeps + 26.0e-18) / (eeps ** 2 + 1.0)
        eeps = (waveno - 111130.0) * 2.0 / 2743.0
        xd2 = (-10.5e-18 * eeps + 46.0e-18) / (eeps ** 2 + 1.0)
        x1240 = xd0 + xd1 + xd2

    if freq >= 2.0761e15:
        xs0 = 10.0 ** (-16.80 - (waveno - 69172.400) / 3.0 / ryd)
        eeps = (waveno - 97700.0) * 2.0 / 2743.0
        xs1 = (68.0e-18 * eeps + 118.0e-18) / (eeps ** 2 + 1.0)
        x1444 = xs0 + xs1

    if freq >= 2.0761e15:
        tkev = state.tkev[:ntau]
        c1240 = 5.0 * np.exp(-1.264 / tkev)
        c1444 = np.exp(-2.683 / tkev)
        nc = state.numdens[2, 0, :ntau]
        uc = state.u[5, 0, :ntau]     # C neutral: Z=6 → u[5, 0, :]
        state.aC1[:ntau] = (x1100 * 9.0 + x1240 * c1240 + x1444 * c1444) * nc / uc


def opac_mg1(state) -> None:
    """Mg I bound-free opacity (Peach 1970 tables)."""
    ntau = state.ntau
    freq = state.freq
    if freq < 2.997925e14:
        return
    freqlg = state.freqlg

    xx = _peach_interp(freq, freqlg, _MG1_FREQMG, _MG1_FLOG, _MG1_TLG, _MG1_PEACH)

    t   = state.t[:ntau]
    tlg = state.tlog[:ntau]
    nmg = state.numdens[3, 0, :ntau]
    umg = state.u[11, 0, :ntau]    # Mg neutral: Z=12 → u[11, 0, :]

    n_T = len(_MG1_TLG)
    for i in range(ntau):
        n = max(0, min(5, round(t[i] / 1000.0) - 4))   # Fortran: max(1,min(6,nint(T/1000)-3))-1
        dt = (tlg[i] - _MG1_TLG[n]) / (_MG1_TLG[n + 1] - _MG1_TLG[n])
        state.aMg1[i] = np.exp(xx[n] * (1.0 - dt) + xx[n + 1] * dt) * nmg[i] / umg[i]


def opac_mg2(state) -> None:
    """Mg II bound-free opacity (edges at 824 Å and 1169 Å)."""
    ntau = state.ntau
    freq = state.freq

    x824  = _seaton(freq, 3.635492e15, 1.40e-19, 4.0, 6.7) if freq >= 3.635492e15 else 1.0e-99
    x1169 = 5.11e-19 * (2.564306e15 / freq) ** 3 if freq >= 2.564306e15 else 1.0e-99

    if x1169 < 1.0e-90:
        return

    tkev = state.tkev[:ntau]
    c1169 = 6.0 * np.exp(-4.43 / tkev)
    nmg2  = state.numdens[3, 1, :ntau]
    umg2  = state.u[11, 1, :ntau]   # Mg+: Z=12, ion=1 → u[11, 1, :]
    state.aMg2[:ntau] = (x824 * 2.0 + x1169 * c1169) * nmg2 / umg2


def opac_al1(state) -> None:
    """Al I bound-free opacity (edge at 1443 Å ~ 2.077e15 Hz)."""
    ntau = state.ntau
    freq = state.freq
    if freq < 1.443e15:
        return
    nal = state.numdens[4, 0, :ntau]
    ual = state.u[12, 0, :ntau]    # Al neutral: Z=13 → u[12, 0, :]
    state.aAl1[:ntau] = 6.5e-17 * (1.443e15 / freq) ** 5 * 6.0 * nal / ual


def opac_si1(state) -> None:
    """Si I bound-free opacity (Peach 1970 tables). Note sign: peach stores +log|σ|."""
    ntau = state.ntau
    freq = state.freq
    if freq < 2.997925e14:
        return
    freqlg = state.freqlg

    # Find frequency interval and map to Peach table columns
    n_edge = len(_SI1_FREQSI)
    n = n_edge
    for k in range(n_edge):
        if freq > _SI1_FREQSI[k]:
            n = k
            break
    n_flog = len(_SI1_FLOG)
    dd  = (freqlg - _SI1_FLOG[n]) / (_SI1_FLOG[n + 1] - _SI1_FLOG[n])
    n_col = n if n <= 1 else 2 * n - 2
    dd1 = 1.0 - dd
    xx  = _SI1_PEACH[:, n_col + 1] * dd + _SI1_PEACH[:, n_col] * dd1

    t   = state.t[:ntau]
    tlg = state.tlog[:ntau]
    nsi = state.numdens[5, 0, :ntau]
    usi = state.u[13, 0, :ntau]    # Si neutral: Z=14 → u[13, 0, :]

    for i in range(ntau):
        n = max(0, min(7, round(t[i] / 1000.0) - 4))
        dt = (tlg[i] - _SI1_TLG[n]) / (_SI1_TLG[n + 1] - _SI1_TLG[n])
        # Note: Si I stores positive log|σ| and we take exp(-xx) unlike Mg I
        state.aSi1[i] = (np.exp(-(xx[n] * (1.0 - dt) + xx[n + 1] * dt)) * 9.0
                         * nsi[i] / usi[i])


def opac_si2(state) -> None:
    """Si II bound-free opacity (Peach 1970 tables)."""
    ntau = state.ntau
    freq = state.freq
    if freq < 7.6869872e14:
        return
    freqlg = state.freqlg

    n_edge = len(_SI2_FREQSI)
    n = n_edge
    for k in range(n_edge):
        if freq > _SI2_FREQSI[k]:
            n = k
            break
    n_flog = len(_SI2_FLOG)
    dd  = (freqlg - _SI2_FLOG[n]) / (_SI2_FLOG[n + 1] - _SI2_FLOG[n])
    n_col = n if n <= 1 else 2 * n - 2
    if n_col >= 13:
        n_col = 12
    dd1 = 1.0 - dd
    xx  = _SI2_PEACH[:, n_col + 1] * dd + _SI2_PEACH[:, n_col] * dd1

    t   = state.t[:ntau]
    tlg = state.tlog[:ntau]
    nsi2 = state.numdens[5, 1, :ntau]
    usi2 = state.u[13, 1, :ntau]   # Si+: Z=14, ion=1 → u[13, 1, :]

    for i in range(ntau):
        n = max(0, min(4, round(t[i] / 2000.0) - 5))
        dt = (tlg[i] - _SI2_TLG[n]) / (_SI2_TLG[n + 1] - _SI2_TLG[n])
        state.aSi2[i] = (np.exp(xx[n] * (1.0 - dt) + xx[n + 1] * dt) * 6.0
                         * nsi2[i] / usi2[i])


def opac_fe1(state) -> None:
    """Fe I bound-free opacity (Bautista 1997 levels)."""
    ntau  = state.ntau
    freq  = state.freq
    waveno = freq / 2.99792458e10   # cm⁻¹

    if waveno < 21000.0:
        return

    # Cross-sections for each Fe I level at current frequency
    xsect = np.where(
        _FE1_WNO < waveno,
        3.0e-18 / (1.0 + ((_FE1_WNO + 3000.0 - waveno) / _FE1_WNO / 0.1) ** 4),
        0.0,
    )   # shape (48,)

    nfe = state.numdens[6, 0, :ntau]
    ufe = state.u[25, 0, :ntau]    # Fe neutral: Z=26 → u[25, 0, :]
    hkt_arr = _H / (_K * state.t[:ntau])

    for i in range(ntau):
        bolt = _FE1_GG * np.exp(-_FE1_EE * 2.99792458e10 * hkt_arr[i])
        state.aFe1[i] = float(np.dot(xsect, bolt)) * nfe[i] / ufe[i]
