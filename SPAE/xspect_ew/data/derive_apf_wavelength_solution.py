"""Derive the edge-corrected APF static wavelength solution shipped as
apf_wavelength_solution.npz (see readers._read_apf()).

The APF reduction's own static solution (recovered from text exports, see
the npz's 'note') is wrong near order edges -- mostly the first ~20% of each
order, up to ~5-10 km/s -- which is fixed in (order, pixel) because the
solution is static. This script measures that error per exposure, fits a
smooth per-order correction, and writes the corrected solution plus the
pixel range of each order that could actually be calibrated.

Usage:
    python derive_apf_wavelength_solution.py measure  <apf_dir> <work_dir>
    python derive_apf_wavelength_solution.py fit      <work_dir>
    python derive_apf_wavelength_solution.py validate <apf_dir> <work_dir>

measure: for every r*.fits in apf_dir, normalized single-exposure flux in
    the OBSERVED frame (the error is fixed on the detector, while the star
    moves ~50 px between exposures, so averaging over exposures averages out
    template mismatch). Two independent offset measurements per chunk
    (CHUNK_PX wide, every STRIDE_PX):
      'atlas': star (shifted to rest by that exposure's RV) vs the Kurucz
               solar flux atlas.
      'overlap': order o+1 vs order o of the SAME exposure, wherever their
               wavelength ranges overlap -- a perfect template, but only
               available near the blue end of o+1 (and the red end of o).
"""
import glob, io, contextlib, os, pickle, sys, warnings
import numpy as np
from scipy.ndimage import gaussian_filter1d
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

C_KMS = 299792.458
N_ORDERS, N_PX = 79, 4608
CHUNK_PX, STRIDE_PX = 115, 57          # ~2.5% of an order, half-overlapping
VGRID = np.arange(-15.0, 15.001, 0.2)   # km/s, searched around the expected position
R_APF = 60000.0                         # atlas broadening (only shapes the template)
ATLAS = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))), 'data', 'kurucz_solar_atlas', 'fluxspliced.2005')
LINELIST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'Line_list_sample', 'Sun_fe_sample.txt')


def chunk_starts():
    return np.arange(0, N_PX - CHUNK_PX + 1, STRIDE_PX)


def ccf_peak(w, y, tw, ty, vgrid=VGRID):
    """Velocity v (km/s) maximizing corr(y(w), template(w / (1 + v/c))),
    parabola-refined; returns (v, peak correlation) or (nan, nan)."""
    q = (w[None, :] / (1.0 + vgrid[:, None] / C_KMS)).ravel()
    m = np.interp(q, tw, ty, left=np.nan, right=np.nan).reshape(len(vgrid), len(w))
    if not np.all(np.isfinite(m)):
        return np.nan, np.nan
    yc = y - y.mean()
    mc = m - m.mean(axis=1, keepdims=True)
    denom = np.sqrt((mc ** 2).sum(axis=1) * (yc ** 2).sum())
    cc = (mc @ yc) / np.where(denom > 0, denom, np.nan)
    if not np.any(np.isfinite(cc)):
        return np.nan, np.nan
    j = int(np.nanargmax(cc))
    if j == 0 or j == len(vgrid) - 1:
        return np.nan, cc[j]
    d = cc[j - 1] - 2 * cc[j] + cc[j + 1]
    step = vgrid[1] - vgrid[0]
    return vgrid[j] + (0.5 * (cc[j - 1] - cc[j + 1]) / d * step if d != 0 else 0.0), cc[j]


def load_exposure(path, wave_override=None):
    """Spectrum_Data in the observed frame, normalized, with its measured RV."""
    from xspect_ew.spectrum_data import Spectrum_Data
    with contextlib.redirect_stdout(io.StringIO()):
        # the pipeline's own solution on the full detector -- this script
        # is what derives the correction and the trimming
        s = Spectrum_Data(path, wavelength_solution='static', trim_edges=False)
        if wave_override is not None:
            for o in range(N_ORDERS):
                s.wavelength[o] = wave_override[o].copy()
                s.shifted_wavelength[o] = wave_override[o].copy()
        s.normalize_all()
        s.load_lines(LINELIST)
        s.apply_rv_shift()
    return s


def measure(apf_dir, work_dir, wave_override=None, tag='static'):
    from xspect_ew.reference_atlas import load_reference_atlas
    from astropy.io import fits
    aw, af = load_reference_atlas(ATLAS)
    files = sorted(glob.glob(os.path.join(apf_dir, '*.fits')))
    starts = chunk_starts()
    centers = starts + CHUNK_PX / 2.0
    n_e, n_k = len(files), len(starts)
    atlas_v = np.full((n_e, N_ORDERS, n_k), np.nan); atlas_cc = np.full_like(atlas_v, np.nan)
    ov_v = np.full_like(atlas_v, np.nan); ov_cc = np.full_like(atlas_v, np.nan)
    meta = []
    for e, f in enumerate(files):
        s = load_exposure(f, wave_override)
        rv = s.rv[0]
        W = [np.asarray(w, float) for w in s.wavelength]          # observed frame
        F = [np.clip(1.0 - np.asarray(n, float), -0.5, 1.5) for n in s.normalized_flux]
        for o in range(N_ORDERS):
            w_rest = W[o] / (1.0 + rv / C_KMS)
            near = (aw > w_rest[0] - 30) & (aw < w_rest[-1] + 30)
            if near.sum() < 100:
                continue
            ta_w = aw[near]
            ta_y = 1.0 - af[near] / np.percentile(af[near], 98)
            ta_y = gaussian_filter1d(ta_y, (np.median(w_rest) / R_APF / 2.355) / np.median(np.diff(ta_w)))
            for k, p0 in enumerate(starts):
                sl = slice(p0, p0 + CHUNK_PX)
                y = F[o][sl]
                if not np.all(np.isfinite(y)) or np.std(y) < 0.01:
                    continue
                atlas_v[e, o, k], atlas_cc[e, o, k] = ccf_peak(w_rest[sl], y, ta_w, ta_y)
                if o > 0:   # same wavelengths in the previous (redder-end) order of this exposure
                    wp = W[o - 1]
                    if W[o][sl][0] > wp[0] + 2 and W[o][sl][-1] < wp[-1] - 2:
                        ov_v[e, o, k], ov_cc[e, o, k] = ccf_peak(W[o][sl], y, wp, F[o - 1],
                                                                 np.arange(-10, 10.001, 0.1))
        h = fits.getheader(f)
        meta.append(dict(file=os.path.basename(f), date=h.get('DATE-OBS'), rv=rv,
                         median_counts=float(np.median(np.concatenate([np.asarray(x) for x in s.flux])))))
        print(f'  [{e + 1}/{n_e}] {os.path.basename(f)} rv={rv:+.2f}', flush=True)
    out = dict(files=files, meta=meta, centers=centers, atlas_v=atlas_v, atlas_cc=atlas_cc,
               ov_v=ov_v, ov_cc=ov_cc, chunk_px=CHUNK_PX)
    os.makedirs(work_dir, exist_ok=True)
    with open(os.path.join(work_dir, f'offsets_{tag}.pkl'), 'wb') as fh:
        pickle.dump(out, fh)
    return out


MIN_CC = 0.6          # chunk CCF peak correlation to count as a measurement
MID = (0.25, 0.75)    # order "middle" (fraction of N_PX) -- the reference level


def relative_maps(m):
    """Per-exposure offsets relative to each order's own middle (km/s), from
    the atlas and from the overlap with the previous order. The atlas map
    carries the star's zero point in each order's middle; subtracting it
    leaves e(order, px) - e(order, middle). The overlap map measures
    e(o, px) - e(o-1, same wavelength) directly; adding the atlas middle
    difference (o-1 minus o) puts it on the same footing, assuming o-1's
    red end sits at its middle's level (true to ~0.4 km/s, see fit())."""
    x = m['centers'] / N_PX
    mid = (x > MID[0]) & (x < MID[1])
    av = np.where(m['atlas_cc'] > MIN_CC, m['atlas_v'], np.nan)
    ov = np.where(m['ov_cc'] > MIN_CC, m['ov_v'], np.nan)
    mids = np.nanmedian(av[:, :, mid], axis=2)                    # (exposure, order)
    atlas_rel = av - mids[:, :, None]
    prev_mid = np.concatenate([np.full((av.shape[0], 1), np.nan), mids[:, :-1]], axis=1)
    ov_rel = ov + (prev_mid - mids)[:, :, None]
    return x, atlas_rel, ov_rel, mids


def epoch_report(m, work_dir):
    """Does the offset map change between observations? Compare each
    exposure's blue-edge offsets (overlap measurement: same-star template,
    so nothing shifts with the exposure's RV) to the all-exposure median."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x, atlas_rel, ov_rel, mids = relative_maps(m)
    lines = []
    for name, rel, raw in (('overlap', ov_rel, np.where(m['ov_cc'] > MIN_CC, m['ov_v'], np.nan)),
                           ('atlas', atlas_rel, None)):
        med = np.nanmedian(rel, axis=0)                           # (order, chunk) consensus
        resid = rel - med[None]
        edge = x < 0.15
        lines.append(f'\n[{name}] per-exposure deviation from the all-exposure median map, '
                     f'blue 15% of each order (km/s):')
        lines.append(f'  {"file":15s} {"date":10s} {"RV":>7s} {"mean":>7s} {"sem":>6s} {"rms":>6s} {"n":>5s}')
        stats = []
        for e, meta in enumerate(m['meta']):
            r = resid[e][:, edge].ravel(); r = r[np.isfinite(r)]
            if len(r) < 10:
                continue
            sem = 1.2533 * np.std(r) / np.sqrt(len(r))
            stats.append((meta['date'], np.mean(r), sem))
            lines.append(f'  {meta["file"]:15s} {meta["date"][:10]:10s} {meta["rv"]:+7.2f} '
                         f'{np.mean(r):+7.2f} {sem:6.2f} {np.std(r):6.2f} {len(r):5d}')
        means = np.array([s[1] for s in stats]); sems = np.array([s[2] for s in stats])
        chi2 = np.sum(((means - np.average(means, weights=1 / sems**2)) / sems) ** 2)
        lines.append(f'  scatter of exposure means {np.std(means):.2f} km/s; chi2 = {chi2:.1f} '
                     f'for {len(means) - 1} dof (consistent with noise if ~{len(means) - 1})')
        if raw is not None:
            # also the raw overlap offset itself (e_o - e_o-1), which needs no atlas at all
            r_raw = np.nanmedian(raw[:, :, edge], axis=2)                     # (exposure, order)
            dev = r_raw - np.nanmedian(r_raw, axis=0)[None]
            lines.append('  raw overlap offset per exposure, deviation from all-exposure median '
                         '(median over orders): ' + ', '.join(
                             f'{meta["file"][:4]} {np.nanmedian(dev[e]):+.2f}' for e, meta in enumerate(m['meta'])))
    fig, axs = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for ax, (name, rel) in zip(axs, (('overlap (same-star template)', ov_rel), ('solar atlas', atlas_rel))):
        for e in range(rel.shape[0]):
            prof = np.nanmedian(rel[e][16:63], axis=0)
            ax.plot(x, prof, lw=0.8, alpha=0.8, label=m['meta'][e]['date'][:10])
        ax.axhline(0, color='k', lw=0.5); ax.set_ylim(-8, 3)
        ax.set_ylabel('offset rel. to order middle (km/s)'); ax.set_title(f'{name}: median over orders 16-62, one line per exposure')
    axs[1].set_xlabel('position along order (fraction)'); axs[0].legend(fontsize=6, ncol=3)
    plt.tight_layout(); plt.savefig(os.path.join(work_dir, 'epoch_comparison.png'), dpi=100)
    report = '\n'.join(lines)
    open(os.path.join(work_dir, 'epoch_report.txt'), 'w').write(report)
    print(report)


KNOTS = [0.06, 0.13, 0.22, 0.35, 0.60, 0.85]   # spline knots, fraction of order
# observed-frame telluric bands (A) -- chunks touching these are not used:
# the bands sit still on the detector while the star moves, so they pull a
# star-vs-template CCF toward the wrong answer
TELLURIC_BANDS = [(5885, 5905), (6275, 6330), (6865, 6935), (7155, 7345), (7585, 7705),
                  (8100, 8360), (8900, 9900)]
CLIP_SIGMA, CLIP_FLOOR = 3.0, 1.5        # outlier rejection: max(3 x robust sigma, 1.5 km/s)
MAX_FIT_RMS = 1.0                        # km/s, after clipping, to accept an order's fit
MIN_EDGE_CHUNKS, EDGE_X = 3, 0.12        # the steep blue edge must actually be measured
MAX_ABS_CORR = 15.0                      # km/s, sanity bound on the correction
MIN_EXPOSURES = 6         # exposures with a good chunk measurement for a consensus point
MIN_GOOD_CHUNKS = 30      # of ~79 per order, for the order to be fitted at all
FALLBACK_TRIM = (0.15, 0.025)   # blue/red fraction trimmed from orders that can't be calibrated
# ... redward of this the uncalibrated orders' blue-edge offset persists further
# (measured median -1.7 km/s at 15-20% of the order vs -0.8 km/s in the blue)
FALLBACK_TRIM_RED, FALLBACK_RED_FROM_AA = (0.20, 0.025), 7100.0
NOTE_STATIC = ("APF (Lick Automated Planet Finder) static wavelength solution, Angstrom (air), "
               "observed (topocentric) frame, shape (79 orders, 4608 px), row 0 = bluest order. "
               "Derived 2026-10-12 from Gaia_DR3_1125337614021152256_epoch_1.txt (K. Taggart, "
               "ziggy.ucolick.org) by dividing out (1 + BC/c), BC = astropy barycentric correction "
               "at Lick for DATE-OBS + EXPTIME/2 and header RA/DEC. All 18 epochs (2024-08 to 2026-05) "
               "give the same solution to 0.2 m/s. Source values were float32-rounded before the BC "
               "shift, so this is within ~18 m/s of the pipeline's own array.")


def consensus(rel):
    """Median over exposures per (order, chunk), with its error, where at
    least MIN_EXPOSURES exposures measured it."""
    n = np.sum(np.isfinite(rel), axis=0)
    med = np.nanmedian(rel, axis=0)
    err = 1.2533 * np.nanstd(rel, axis=0) / np.sqrt(np.maximum(n, 1))
    bad = n < MIN_EXPOSURES
    med[bad] = np.nan; err[bad] = np.nan
    return med, err, n


def fit(work_dir, wave_static):
    """Fit the per-order correction; returns (wave_corrected, valid_lo,
    valid_hi, report lines, per-order info)."""
    from scipy.interpolate import LSQUnivariateSpline
    with open(os.path.join(work_dir, 'offsets_static.pkl'), 'rb') as fh:
        m = pickle.load(fh)
    x, atlas_rel, ov_rel, _ = relative_maps(m)
    a_med, a_err, _ = consensus(atlas_rel)
    o_med, o_err, _ = consensus(ov_rel)
    px = np.arange(N_PX); xp = px / N_PX
    wave = wave_static.copy()
    valid_lo = np.zeros(N_ORDERS, int); valid_hi = np.full(N_ORDERS, N_PX - 1)
    info, rep = [], []
    for o in range(N_ORDERS):
        cw = np.interp(x * N_PX, px, wave_static[o])                   # chunk centre wavelengths
        half = CHUNK_PX / 2 * np.median(np.diff(wave_static[o]))
        tell = np.zeros(len(x), bool)
        for lo_b, hi_b in TELLURIC_BANDS:
            tell |= (cw + half > lo_b) & (cw - half < hi_b)
        good = np.isfinite(a_med[o]) & (a_err[o] > 0) & ~tell
        # merge the overlap measurement where available (independent of the
        # atlas template at that chunk); record their agreement
        ys_all = a_med[o].copy(); ws_all = 1.0 / np.maximum(a_err[o], 0.1)
        both = good & np.isfinite(o_med[o])
        agree = np.nanmedian(np.abs(a_med[o][both] - o_med[o][both])) if both.any() else np.nan
        wo = 1.0 / np.maximum(o_err[o], 0.1)
        ys_all[both] = (ys_all[both] * ws_all[both] ** 2 + o_med[o][both] * wo[both] ** 2) / (ws_all[both] ** 2 + wo[both] ** 2)
        ws_all[both] = np.sqrt(ws_all[both] ** 2 + wo[both] ** 2)

        def fallback(why):
            trim = FALLBACK_TRIM_RED if wave_static[o][0] > FALLBACK_RED_FROM_AA else FALLBACK_TRIM
            lo, hi = int(trim[0] * N_PX), int((1 - trim[1]) * N_PX) - 1
            valid_lo[o], valid_hi[o] = lo, hi
            info.append(dict(order=o, fitted=False, n_good=int(good.sum()), why=why))
            rep.append(f'order {o:2d} {wave_static[o][0]:7.1f}-{wave_static[o][-1]:7.1f}: NOT calibrated '
                       f'({why}) -- trimmed to px {lo}-{hi}')

        if good.sum() < MIN_GOOD_CHUNKS or not np.any(good & (x > MID[0]) & (x < MID[1])):
            fallback(f'{good.sum()} usable chunks'); continue
        use = good.copy()
        for _ in range(4):                                             # iterative outlier rejection
            xs, ys, ws = x[use], ys_all[use], ws_all[use]
            knots, prev = [], xs[0]
            for t in KNOTS:
                if xs[0] < t < xs[-1] and np.sum((xs > prev) & (xs <= t)) >= 3:
                    knots.append(t); prev = t
            spl = LSQUnivariateSpline(xs, ys, knots, w=ws, k=3)
            r = ys_all - spl(np.clip(x, xs[0], xs[-1]))
            sig = 1.4826 * np.nanmedian(np.abs(r[use]))
            new_use = good & (np.abs(r) <= max(CLIP_SIGMA * sig, CLIP_FLOOR))
            if np.array_equal(new_use, use):
                break
            use = new_use
        xs, ys = x[use], ys_all[use]
        zero = np.median(spl(np.linspace(MID[0], MID[1], 50)))
        resid = ys - (spl(xs) - zero) - zero
        rms = float(np.std(resid))
        n_edge = int(np.sum(use & (x < EDGE_X)))
        lo = max(int(round(xs[0] * N_PX - CHUNK_PX / 2)), 0)
        hi = min(int(round(xs[-1] * N_PX + CHUNK_PX / 2)) - 1, N_PX - 1)
        corr = spl(np.clip(xp, xs[0], xs[-1])) - zero
        inside = corr[lo:hi + 1]
        if n_edge < MIN_EDGE_CHUNKS:
            fallback(f'blue edge unmeasured: {n_edge} usable chunks below {EDGE_X:.0%}'); continue
        if rms > MAX_FIT_RMS:
            fallback(f'fit rms {rms:.2f} km/s'); continue
        if (np.max(np.abs(inside)) > MAX_ABS_CORR or inside.max() > ys.max() - zero + 1.5
                or inside.min() < ys.min() - zero - 1.5):
            fallback(f'spline leaves its data range (max |{np.max(np.abs(inside)):.1f}| km/s)'); continue
        wave[o] = wave_static[o] / (1.0 + corr / C_KMS)
        valid_lo[o], valid_hi[o] = lo, hi
        info.append(dict(order=o, fitted=True, n_good=int(use.sum()), n_clipped=int(good.sum() - use.sum()),
                         lo=lo, hi=hi, corr_lo=float(corr[lo]), resid_rms=rms, atlas_vs_overlap=float(agree)))
        rep.append(f'order {o:2d} {wave_static[o][0]:7.1f}-{wave_static[o][-1]:7.1f}: {use.sum():2d} chunks '
                   f'({good.sum() - use.sum()} clipped), valid px {lo}-{hi}, correction at px {lo} '
                   f'{corr[lo]:+5.2f} km/s, fit rms {rms:.2f}, |atlas-overlap| {agree:.2f}')
    return wave, valid_lo, valid_hi, rep, info


def compare(work_dir):
    """Offsets along the order before vs after the correction, using only
    each order's valid (kept) pixel range."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cand = np.load(os.path.join(work_dir, 'apf_wavelength_solution_candidate.npz'))
    info = pickle.load(open(os.path.join(work_dir, 'fit_info.pkl'), 'rb'))
    lo, hi = cand['valid_lo'], cand['valid_hi']
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    lines = []
    for ax, tag in zip(axs, ('static', 'corrected')):
        m = pickle.load(open(os.path.join(work_dir, f'offsets_{tag}.pkl'), 'rb'))
        x, a_rel, _, _ = relative_maps(m)
        med, _, _ = consensus(a_rel)
        for group, orders in (('calibrated', [i['order'] for i in info if i['fitted']]),
                              ('trimmed only', [i['order'] for i in info if not i['fitted']])):
            lines.append(f'[{tag}] {group} orders, kept pixels only: offset vs position (km/s)')
            for a, b in ((0, .025), (.025, .05), (.05, .1), (.1, .15), (.15, .2), (.2, .3), (.3, .7), (.7, .9), (.9, .95), (.95, 1)):
                vals = []
                for o in orders:
                    sel = (x >= a) & (x < b) & (x * N_PX >= lo[o]) & (x * N_PX <= hi[o])
                    vals += list(med[o][sel][np.isfinite(med[o][sel])])
                if len(vals) >= 5:
                    v = np.array(vals)
                    lines.append(f'    {a:4.2f}-{b:4.2f}: median {np.median(v):+5.2f}, 16-84% [{np.percentile(v, 16):+5.2f},'
                                 f' {np.percentile(v, 84):+5.2f}], n={len(v)}')
        for o in [i['order'] for i in info if i['fitted']]:
            sel = (x * N_PX >= lo[o]) & (x * N_PX <= hi[o])
            ax.plot(x[sel], med[o][sel], lw=0.6, alpha=0.6)
        ax.axhline(0, color='k', lw=0.5); ax.set_ylim(-12, 6)
        ax.set_title(f'{tag} solution: calibrated orders, kept pixels'); ax.set_xlabel('position along order')
    axs[0].set_ylabel('offset rel. to order middle (km/s)')
    plt.tight_layout(); plt.savefig(os.path.join(work_dir, 'before_after.png'), dpi=100)
    report = '\n'.join(lines)
    open(os.path.join(work_dir, 'compare_report.txt'), 'w').write(report)
    print(report)


def install(work_dir):
    cand = np.load(os.path.join(work_dir, 'apf_wavelength_solution_candidate.npz'))
    info = pickle.load(open(os.path.join(work_dir, 'fit_info.pkl'), 'rb'))
    calibrated = np.array([i['fitted'] for i in info])
    note = (NOTE_STATIC + ' -- wave_static is that solution. wave is it with the order-edge error '
            'removed (2026-10-12, derive_apf_wavelength_solution.py on the same 18 exposures: per-order '
            'cubic spline of the offset vs pixel, measured against the Kurucz solar flux atlas and '
            'adjacent-order overlaps, telluric bands masked, anchored at 0 in each order middle; '
            f'{int(calibrated.sum())} of {N_ORDERS} orders calibrated, the rest left as wave_static). '
            'valid_lo/valid_hi: inclusive pixel range of each order to keep -- the calibrated range, or '
            f'for uncalibrated orders a trim of {FALLBACK_TRIM[0]:.0%} (blue end; {FALLBACK_TRIM_RED[0]:.0%} '
            f'redward of {FALLBACK_RED_FROM_AA:.0f} A) and {FALLBACK_TRIM[1]:.1%} (red end).')
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'apf_wavelength_solution.npz')
    np.savez_compressed(out, wave=cand['wave'], wave_static=cand['wave_static'],
                        valid_lo=cand['valid_lo'], valid_hi=cand['valid_hi'], calibrated=calibrated,
                        note=np.array(note))
    print(f'wrote {out}')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'measure':
        measure(sys.argv[2], sys.argv[3])
    elif cmd == 'fit':
        work_dir = sys.argv[2]
        with np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'apf_wavelength_solution.npz')) as d:
            wave_static = d['wave_static'] if 'wave_static' in d else d['wave']
        wave, lo, hi, rep, info = fit(work_dir, wave_static)
        open(os.path.join(work_dir, 'fit_report.txt'), 'w').write('\n'.join(rep))
        np.savez_compressed(os.path.join(work_dir, 'apf_wavelength_solution_candidate.npz'),
                            wave=wave, wave_static=wave_static, valid_lo=lo, valid_hi=hi)
        with open(os.path.join(work_dir, 'fit_info.pkl'), 'wb') as fh:
            pickle.dump(info, fh)
        print('\n'.join(rep))
    elif cmd == 'validate':
        apf_dir, work_dir = sys.argv[2], sys.argv[3]
        with np.load(os.path.join(work_dir, 'apf_wavelength_solution_candidate.npz')) as d:
            measure(apf_dir, work_dir, wave_override=d['wave'], tag='corrected')
    elif cmd == 'compare':
        work_dir = sys.argv[2]
        compare(work_dir)
    elif cmd == 'install':
        work_dir = sys.argv[2]
        install(work_dir)
    elif cmd == 'epochs':
        with open(os.path.join(sys.argv[2], 'offsets_static.pkl'), 'rb') as fh:
            epoch_report(pickle.load(fh), sys.argv[2])
    else:
        raise SystemExit(f'unknown command {cmd}')
