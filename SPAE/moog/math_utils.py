"""
Math primitives translated from MOOG Fortran.

    rinteg  — piecewise-quadratic integration (ATLAS6 scheme), from Rinteg.f
    voigt   — Voigt profile approximation (Landolt-Börnstein), from Voigt.f
    invert  — square matrix inversion, from Invert.f
"""

import numpy as np


# ---------------------------------------------------------------------------
# rinteg  (Rinteg.f)
# ---------------------------------------------------------------------------

def _parcoe(f, x):
    """
    Compute piecewise-quadratic coefficients a, b, c such that
    f(x) ≈ a[i] + b[i]*x + c[i]*x²  on the interval [x[i], x[i+1]].

    This is the ATLAS6 scheme used in MOOG.  The boundary treatment is:
      • segments 0 and 1   : forced linear  (c = 0)
      • segment  n-2       : copied from n-1 (linear at upper end)
      • interior segments  : 3-point quadratic
    """
    n = len(f)
    a = np.empty(n)
    b = np.empty(n)
    c = np.zeros(n)

    # ---- boundary: segment 0 (linear) ----
    b[0] = (f[1] - f[0]) / (x[1] - x[0])
    a[0] = f[0] - x[0] * b[0]

    # ---- boundary: segment n-1 (linear) ----
    b[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2])
    a[n - 1] = f[n - 1] - x[n - 1] * b[n - 1]

    if n == 2:
        return a, b, c

    # ---- interior: 3-point quadratic (vectorized) ----
    j   = np.arange(1, n - 1)
    jm1 = j - 1
    jp1 = j + 1
    dx_lo = x[j] - x[jm1]
    dx_hi = x[jp1] - x[j]
    dx_sp = x[jp1] - x[jm1]
    d     = (f[j] - f[jm1]) / dx_lo
    c[j]  = (f[jp1] / (dx_hi * dx_sp)
             - f[j]   / (dx_lo * dx_hi)
             + f[jm1] / (dx_lo * dx_sp))
    b[j]  = d - (x[j] + x[jm1]) * c[j]
    a[j]  = f[jm1] - x[jm1] * d + x[j] * x[jm1] * c[j]

    # ---- override segments 1 and 2 with linear (ATLAS6 edge treatment) ----
    if n > 2:
        c[1] = 0.0
        b[1] = (f[2] - f[1]) / (x[2] - x[1])
        a[1] = f[1] - x[1] * b[1]
    if n > 3:
        c[2] = 0.0
        b[2] = (f[3] - f[2]) / (x[3] - x[2])
        a[2] = f[2] - x[2] * b[2]

    # ---- blend or copy last interior segment (ATLAS6 upper-end treatment) ----
    j = n - 2   # last interior index
    if c[j] != 0.0:
        j1 = j + 1   # = n-1
        wt = abs(c[j1]) / (abs(c[j1]) + abs(c[j]))
        a[j] = a[j1] + wt * (a[j] - a[j1])
        b[j] = b[j1] + wt * (b[j] - b[j1])
        c[j] = c[j1] + wt * (c[j] - c[j1])
    # copy n-1 coefficients to n-2 (linear at the upper boundary)
    a[n - 2] = a[n - 1]
    b[n - 2] = b[n - 1]
    c[n - 2] = c[n - 1]

    return a, b, c


def rinteg(x, f, n, start):
    """
    Integrate f(x) over x[0:n] using the ATLAS6 piecewise-quadratic scheme.

    Parameters
    ----------
    x : array-like, length >= n
        Independent variable (e.g. rhox or tauref).
    f : array-like, length >= n
        Integrand values at each x.
    n : int
        Number of points to use.
    start : float
        Value of the integral at x[0] (boundary condition).

    Returns
    -------
    total : float
        Cumulative integral from x[0] to x[n-1].
    fint : ndarray, shape (n,)
        Increment contributed by each interval; fint[0] = start,
        fint[i] = integral over [x[i-1], x[i]] for i >= 1.
    """
    xv = np.asarray(x, dtype=float)[:n]
    fv = np.asarray(f, dtype=float)[:n]

    a, b, c = _parcoe(fv, xv)

    # vectorized integration (replaces per-interval Python loop)
    dx = xv[1:] - xv[:-1]
    sx = xv[1:] + xv[:-1]
    contrib = (a[:-1] + b[:-1] / 2.0 * sx
               + c[:-1] / 3.0 * (sx * xv[1:] + xv[:-1] * xv[:-1])) * dx
    fint = np.empty(n)
    fint[0] = start
    fint[1:] = contrib
    total = start + contrib.sum()

    return total, fint


# ---------------------------------------------------------------------------
# voigt  (Voigt.f)
# ---------------------------------------------------------------------------

def voigt(a, v):
    """
    Voigt profile H(a, v) / sqrt(π), normalised so that H(0,0)=1.

    Fully vectorized numpy implementation of the Landolt-Börnstein approximation
    used in MOOG (Voigt.f). Handles both scalar and array inputs.

    Parameters
    ----------
    a : float or ndarray
        Damping parameter (ratio of Lorentzian to Doppler HWHM).
    v : float or ndarray
        Frequency offset in Doppler units.

    Returns
    -------
    float or ndarray
        Voigt function value.
    """
    scalar = np.ndim(a) == 0 and np.ndim(v) == 0
    a = np.atleast_1d(np.asarray(a, dtype=float))
    v = np.atleast_1d(np.asarray(v, dtype=float))

    a2 = a * a
    v2 = v * v
    h0 = np.exp(-v2)

    # ---- case 1: large-a approximation ----
    u1   = 1.4142136 * (a2 + v2)
    u1s  = np.where(u1 == 0.0, 1.0, u1)   # guard divide; only used where mask_c1
    val1 = (0.7978847 * a / u1s
            * (1.0 + (3.0 * v2 - a2) / u1s**2
               + (15.0 * v2**2 - 30.0 * a2 * v2 + 3.0 * a2**2) / u1s**4)
            ) / 1.772454

    # ---- cases 2 & 3: polynomial in v (h1) ----
    h1_lo  = (-1.12470432  - 0.15516677 * v + 3.28867591 * v2
              - 2.34357915 * v * v2 + 0.42139162 * v2 * v2)
    h1_mid = (-4.48480194  + 9.39456063 * v - 6.61487486 * v2
              + 1.98919585 * v * v2 - 0.2204165  * v2 * v2)
    safe_v2 = np.where(np.abs(v2 - 1.5) < 1e-30, 1.0, v2 - 1.5)
    h1_hi  = ((0.554153432 + 0.278711796 * v - 0.188325687 * v2
               + 0.042991293 * v * v2 - 0.003278278 * v2 * v2) / safe_v2)
    h1_234 = np.where(v < 1.3, h1_lo, np.where(v < 2.4, h1_mid, h1_hi))
    h2     = (1.0 - 2.0 * v2) * h0
    u234   = (0.979895023 - 0.962846325 * a + 0.532770573 * a2
              - 0.122727278 * a * a2)

    # case 3 (a <= 0.2, v < 5)
    val3 = (h0 + h1_234 * a + h2 * a2) / 1.772454

    # case 4 (a <= 0.2, v >= 5)
    safe_v2_45 = np.where(v2 == 0.0, 1.0, v2)
    val4 = (a / (1.772454 * safe_v2_45)
            * (1.0 + 1.5 / safe_v2_45 + 3.75 / (safe_v2_45 * safe_v2_45))
            / 1.772454)

    # case 2 (0.2 < a <= 1.4, a+v <= 3.2)
    h1_c2 = h1_234 + 1.1283790 * h0
    h2_c2 = h2 - h0 + 1.1283790 * h1_c2
    h3_c2 = 0.37612635 * (1.0 - h2) - 0.6666667 * v2 * h1_c2 + 1.1283790 * h2_c2
    h4_c2 = 0.6666667 * v2 * v2 * h0 - 0.37612635 * h1_c2 + 1.1283790 * h3_c2
    val2  = (u234 * (h0 + h1_c2 * a + h2_c2 * a2
                     + h3_c2 * a * a2 + h4_c2 * a2 * a2) / 1.772454)

    # case 5 (a == 0)
    val5 = h0 / 1.772454

    # ---- select by case masks ----
    mask_c1 = (a > 1.4) | ((a > 0.2) & ((a + v) > 3.2))
    mask_c5 = (a == 0.0)
    mask_c4 = (~mask_c1) & (~mask_c5) & (a <= 0.2) & (v >= 5.0)
    mask_c3 = (~mask_c1) & (~mask_c5) & (~mask_c4) & (a <= 0.2)

    result = np.where(mask_c5, val5,
             np.where(mask_c1, val1,
             np.where(mask_c4, val4,
             np.where(mask_c3, val3,
                      val2))))

    return float(result[0]) if scalar else result


# ---------------------------------------------------------------------------
# Batched rinteg helpers  (vectorized over a batch dimension)
# ---------------------------------------------------------------------------

def _parcoe_batch(f_batch, x):
    """
    Batched _parcoe: same ATLAS6 scheme as _parcoe but for f_batch (B, n).
    x is shared across all batch elements.
    Returns a_b, b_b, c_b each shape (B, n).
    """
    B, n = f_batch.shape
    xv = x

    a_b = np.empty((B, n))
    b_b = np.empty((B, n))
    c_b = np.zeros((B, n))

    # segment 0: linear boundary
    dx01 = xv[1] - xv[0]
    b_b[:, 0] = (f_batch[:, 1] - f_batch[:, 0]) / dx01
    a_b[:, 0] = f_batch[:, 0] - xv[0] * b_b[:, 0]

    # segment n-1: linear boundary
    dx_last = xv[n - 1] - xv[n - 2]
    b_b[:, n - 1] = (f_batch[:, n - 1] - f_batch[:, n - 2]) / dx_last
    a_b[:, n - 1] = f_batch[:, n - 1] - xv[n - 1] * b_b[:, n - 1]

    if n > 2:
        j   = np.arange(1, n - 1)           # interior indices (n-2,)
        jm1 = j - 1
        jp1 = j + 1
        dx_lo = xv[j] - xv[jm1]
        dx_hi = xv[jp1] - xv[j]
        dx_sp = xv[jp1] - xv[jm1]

        c_b[:, j] = (f_batch[:, jp1] / (dx_hi * dx_sp)[None, :]
                     - f_batch[:, j]  / (dx_lo * dx_hi)[None, :]
                     + f_batch[:, jm1] / (dx_lo * dx_sp)[None, :])
        d_int = (f_batch[:, j] - f_batch[:, jm1]) / dx_lo[None, :]
        b_b[:, j] = d_int - (xv[j] + xv[jm1])[None, :] * c_b[:, j]
        a_b[:, j] = (f_batch[:, jm1] - xv[jm1][None, :] * d_int
                     + (xv[j] * xv[jm1])[None, :] * c_b[:, j])

        # ATLAS6: force linear at segments 1 and 2
        b_b[:, 1] = (f_batch[:, 2] - f_batch[:, 1]) / (xv[2] - xv[1])
        a_b[:, 1] = f_batch[:, 1] - xv[1] * b_b[:, 1]
        c_b[:, 1] = 0.0
        if n > 3:
            b_b[:, 2] = (f_batch[:, 3] - f_batch[:, 2]) / (xv[3] - xv[2])
            a_b[:, 2] = f_batch[:, 2] - xv[2] * b_b[:, 2]
            c_b[:, 2] = 0.0

        # ATLAS6 upper end: copy n-1 coefficients to n-2
        a_b[:, n - 2] = a_b[:, n - 1]
        b_b[:, n - 2] = b_b[:, n - 1]
        c_b[:, n - 2] = c_b[:, n - 1]

    return a_b, b_b, c_b


def _rinteg_batch_contrib(x, f_batch):
    """
    Return per-interval contrib array of shape (B, n-1) using _parcoe_batch.
    """
    xv = np.asarray(x, dtype=float)
    a_b, b_b, c_b = _parcoe_batch(f_batch, xv)
    dx   = xv[1:] - xv[:-1]       # (n-1,)
    sx   = xv[1:] + xv[:-1]       # (n-1,)
    sx_xp1 = sx * xv[1:]          # (n-1,)
    xi2    = xv[:-1] ** 2          # (n-1,)
    return ((a_b[:, :-1]
             + b_b[:, :-1] / 2.0 * sx[None, :]
             + c_b[:, :-1] / 3.0 * (sx_xp1[None, :] + xi2[None, :]))
            * dx[None, :])         # (B, n-1)


def rinteg_batch_total(x, f_batch, start_batch):
    """
    Batched rinteg returning only the total integral.

    Parameters
    ----------
    x           : (n,)  shared abscissa
    f_batch     : (B, n)  integrand values per batch element
    start_batch : (B,)   boundary value at x[0]

    Returns
    -------
    total : (B,)  integral from x[0] to x[n-1]
    """
    contrib = _rinteg_batch_contrib(x, f_batch)   # (B, n-1)
    return start_batch + contrib.sum(axis=1)


def rinteg_batch_cumsum(x, f_batch, start_batch):
    """
    Batched rinteg returning cumulative integrals — matches np.cumsum(fint).

    Parameters
    ----------
    x           : (n,)
    f_batch     : (B, n)
    start_batch : (B,)

    Returns
    -------
    result : (B, n)  cumulative integral at each x[i]
    """
    B = f_batch.shape[0]
    n = len(x)
    contrib = _rinteg_batch_contrib(x, f_batch)   # (B, n-1)
    fint = np.empty((B, n))
    fint[:, 0] = start_batch
    fint[:, 1:] = contrib
    return np.cumsum(fint, axis=1)


# ---------------------------------------------------------------------------
# invert  (Invert.f)
# ---------------------------------------------------------------------------

def invert(a):
    """
    Invert a square matrix.

    Wraps numpy.linalg.inv, which is more numerically stable than the
    Gaussian elimination in the original Fortran.  Returns the zero matrix
    with a warning if the matrix is singular (matching MOOG's behavior of
    printing a warning and returning without raising).

    Parameters
    ----------
    a : ndarray, shape (n, n)
        Matrix to invert.

    Returns
    -------
    ndarray, shape (n, n)
        Inverse of a, or zeros if singular.
    """
    try:
        return np.linalg.inv(a)
    except np.linalg.LinAlgError:
        import warnings
        warnings.warn("WARNING: AN UN-INVERTABLE ARRAY HAS BEEN ENCOUNTERED!")
        return np.zeros_like(a)
