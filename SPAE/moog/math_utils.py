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

    # ---- interior: 3-point quadratic ----
    for j in range(1, n - 1):
        d = (f[j] - f[j - 1]) / (x[j] - x[j - 1])
        c[j] = (f[j + 1] / ((x[j + 1] - x[j]) * (x[j + 1] - x[j - 1]))
                - f[j]   / ((x[j] - x[j - 1]) * (x[j + 1] - x[j]))
                + f[j - 1] / ((x[j] - x[j - 1]) * (x[j + 1] - x[j - 1])))
        b[j] = d - (x[j] + x[j - 1]) * c[j]
        a[j] = f[j - 1] - x[j - 1] * d + x[j] * x[j - 1] * c[j]

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

    fint = np.empty(n)
    fint[0] = start
    total = start

    for i in range(n - 1):
        dx = xv[i + 1] - xv[i]
        sx = xv[i + 1] + xv[i]
        contrib = (a[i] + b[i] / 2.0 * sx
                   + c[i] / 3.0 * (sx * xv[i + 1] + xv[i] * xv[i])) * dx
        fint[i + 1] = contrib
        total += contrib

    return total, fint


# ---------------------------------------------------------------------------
# voigt  (Voigt.f)
# ---------------------------------------------------------------------------

def _voigt_scalar(a, v):
    """Voigt scalar core — direct translation of Voigt.f."""
    a2 = a * a
    v2 = v * v

    # case 5: a == 0
    if a == 0.0:
        return np.exp(-v2) / 1.772454

    # case 1: a > 1.4  OR  (a > 0.2 and a+v > 3.2)
    if a > 1.4 or (a > 0.2 and a + v > 3.2):
        u = 1.4142136 * (a2 + v2)
        r = (0.7978847 * a / u
             * (1.0 + (3.0 * v2 - a2) / u**2
                + (15.0 * v2**2 - 30.0 * a2 * v2 + 3.0 * a2**2) / u**4))
        return r / 1.772454

    # cases 2 & 3: a <= 1.4 and a+v <= 3.2 (enter via Fortran label 30)
    u = (0.979895023 - 0.962846325 * a + 0.532770573 * a2
         - 0.122727278 * a * a2)
    h0 = np.exp(-v2)

    if v < 1.3:
        h1 = (-1.12470432 - 0.15516677 * v + 3.28867591 * v2
              - 2.34357915 * v * v2 + 0.42139162 * v2 * v2)
    elif v < 2.4:
        h1 = (-4.48480194 + 9.39456063 * v - 6.61487486 * v2
              + 1.98919585 * v * v2 - 0.2204165 * v2 * v2)
    else:
        h1 = ((0.554153432 + 0.278711796 * v - 0.188325687 * v2
               + 0.042991293 * v * v2 - 0.003278278 * v2 * v2)
              / (v2 - 1.5))

    h2 = (1.0 - 2.0 * v2) * h0

    # case 3: a <= 0.2  (Fortran label 52: just h0 + h1*a + h2*a2)
    if a <= 0.2:
        # case 4: a <= 0.2 and v >= 5
        if v >= 5.0:
            return (a / (1.772454 * v2)
                    * (1.0 + 1.5 / v2 + 3.75 / (v2 * v2)) / 1.772454)
        return (h0 + h1 * a + h2 * a2) / 1.772454

    # case 2: 0.2 < a <= 1.4 and a+v <= 3.2
    h1 = h1 + 1.1283790 * h0
    h2p = h2
    h2 = h2 - h0 + 1.1283790 * h1
    h3 = 0.37612635 * (1.0 - h2p) - 0.6666667 * v2 * h1 + 1.1283790 * h2
    h4 = 0.6666667 * v2 * v2 * h0 - 0.37612635 * h1 + 1.1283790 * h3
    return u * (h0 + h1 * a + h2 * a2 + h3 * a * a2 + h4 * a2 * a2) / 1.772454


_voigt_vec = np.vectorize(_voigt_scalar, otypes=[float])


def voigt(a, v):
    """
    Voigt profile H(a, v) / sqrt(π), normalised so that H(0,0)=1.

    Approximation from Landolt-Börnstein New Series — Astronomy &
    Astrophysics, as used in MOOG.  Handles both scalar and array inputs.

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
    result = _voigt_vec(a, v)
    return result.item() if scalar else result


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
