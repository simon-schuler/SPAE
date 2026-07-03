"""
Abundance statistics for a single species.
Translated from Stats.f.

Computes mean, standard deviation, and linear regression of derived
abundances against excitation potential, reduced equivalent width, and
wavelength for lines lim1obs..lim2obs.

Results are written to state and also returned as a dict.
"""
import math
import numpy as np

_SENTINEL = 999.99   # abundout value marking a failed/skipped line


def stats(state) -> dict:
    """
    Compute statistics over state.abundout[lim1obs..lim2obs].

    Writes to state:
      average, deviate, kount
      xxm1/xxb1/xxr1  (excitation potential trend)
      xxm2/xxb2/xxr2  (reduced equivalent width trend)
      xxm3/xxb3/xxr3  (wavelength trend)
      deltaep, deltarw, deltawv

    Returns a dict with the same keys for the caller's convenience.
    """
    lo = state.lim1obs
    hi = state.lim2obs

    # Mean
    total = 0.0
    kount = 0
    for l in range(lo, hi + 1):
        if state.abundout[l] != _SENTINEL:
            total += state.abundout[l]
            kount += 1
    average = total / kount if kount > 0 else 0.0

    # Standard deviation
    deviate = 0.0
    if kount > 1:
        ss = sum((state.abundout[l] - average) ** 2
                 for l in range(lo, hi + 1)
                 if state.abundout[l] != _SENTINEL)
        deviate = math.sqrt(ss / (kount - 1))

    state.average = average
    state.deviate = deviate
    state.kount   = kount

    # Linear regressions (only meaningful with > 2 lines)
    result = {
        'average': average, 'deviate': deviate, 'kount': kount,
        'ep_slope': None, 'ep_intercept': None, 'ep_r': None,
        'rw_slope': None, 'rw_intercept': None, 'rw_r': None,
        'wv_slope': None, 'wv_intercept': None, 'wv_r': None,
        'deltaep': 0.0,   'deltarw': 0.0,       'deltawv': 0.0,
    }

    if kount <= 2:
        return result

    epmin = rwmin = wvmin =  1e30
    epmax = rwmax = wvmax = -1e30
    x1=x2=x3=x4=x5=x6 = 0.0
    y1=y2 = xy = yz = za = 0.0

    for l in range(lo, hi + 1):
        if state.abundout[l] == _SENTINEL:
            continue
        rw = math.log10(state.width[l] / state.wave1[l])
        ep = state.e[l, 0]
        wv = state.wave1[l]
        ab = state.abundout[l]

        x1 += ep;    x2 += ep * ep
        x3 += rw;    x4 += rw * rw
        x5 += wv;    x6 += wv * wv
        y1 += ab;    y2 += ab * ab
        xy += ep * ab
        yz += rw * ab
        za += wv * ab

        epmin = min(epmin, ep);  epmax = max(epmax, ep)
        rwmin = min(rwmin, rw);  rwmax = max(rwmax, rw)
        wvmin = min(wvmin, wv);  wvmax = max(wvmax, wv)

    deltaep = epmax - epmin
    deltarw = rwmax - rwmin
    deltawv = wvmax - wvmin

    def _lsq(sx, sx2, sxy, n, sy, sy2):
        denom = n * sx2 - sx * sx
        if denom == 0.0:
            return 0.0, 0.0, 0.0
        m = (n * sxy - sx * sy) / denom
        b = (sy * sx2 - sxy * sx) / denom
        denom_r = math.sqrt(abs(denom) * abs(n * sy2 - sy * sy))
        r = (n * sxy - sx * sy) / denom_r if denom_r > 0 else 0.0
        return m, b, r

    xxm1, xxb1, xxr1 = _lsq(x1, x2, xy, kount, y1, y2)
    xxm2, xxb2, xxr2 = _lsq(x3, x4, yz, kount, y1, y2)
    xxm3, xxb3, xxr3 = _lsq(x5, x6, za, kount, y1, y2)

    # Write trend coefficients to state
    state.xxm1 = xxm1;  state.xxb1 = xxb1;  state.xxr1 = xxr1
    state.xxm2 = xxm2;  state.xxb2 = xxb2;  state.xxr2 = xxr2
    state.xxm3 = xxm3;  state.xxb3 = xxb3;  state.xxr3 = xxr3
    state.deltaep = deltaep
    state.deltarw = deltarw
    state.deltawv = deltawv

    result.update({
        'ep_slope': xxm1, 'ep_intercept': xxb1, 'ep_r': xxr1,
        'rw_slope': xxm2, 'rw_intercept': xxb2, 'rw_r': xxr2,
        'wv_slope': xxm3, 'wv_intercept': xxb3, 'wv_r': xxr3,
        'deltaep': deltaep, 'deltarw': deltarw, 'deltawv': deltawv,
    })
    return result
