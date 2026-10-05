"""
Universal spectrum file reader: given just a filename, inspects the file's
header/structure to figure out how to extract per-order wavelength
(Angstrom) and flux arrays -- no need for the caller to know or specify the
instrument.

Detection is a registry of (detector, extractor, name) entries, tried in
order. Adding support for a new instrument/pipeline means adding one new
entry, not redesigning the cascade -- this is deliberate, since we will
keep encountering spectrum files that don't match anything here yet.

Design preference throughout: choose RAW/uncalibrated flux (and
uncorrected wavelength) over any already-normalized, telluric-corrected, or
RV-corrected quantity when a format offers a choice, since Spectrum_Data
does its own continuum normalization and wavelength-shift correction --
we don't want to double up on that.

If no known format matches, read_spectrum() raises SpectrumFormatError with
a dump of the file's actual header/structure (so diagnosing a new format is
fast, not a mystery) and points at two ways to proceed without waiting on a
SPAE change: (1) extract (wavelength, flux) yourself and pass them via
Spectrum_Data(name, KECK_file=False, spectx=..., specty=...), the existing
raw-array path, unchanged and still fully supported; or (2) pass a
custom_reader= callable to read_spectrum()/Spectrum_Data() for a one-off
format without editing this file at all.
"""

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter


class SpectrumFormatError(Exception):
    pass


# ---------------------------------------------------------------------------
# FITS: GRACES / CFHT OPERA pipeline
# ---------------------------------------------------------------------------
# Self-documents its row layout via COL1..COLn header keywords (verified
# against a real reduced file, e.g. COL1='Order', COL5='Wavelength' (nm,
# uncorrected), COL9='RawFlux') -- not assumed from docs. Orders are split
# directly on the 'Order' row rather than by detecting wavelength resets.

_GRACES_WAVE_PREFERENCE = ['Wavelength', 'Tell', 'RVel']
_GRACES_FLUX_PREFERENCE = ['RawFlux', 'FcalFlux', 'NormalizedFlux']
_GRACES_ORDER_LABEL = 'Order'


def _graces_col_map(header):
    """{column label: data-row index} from COL1..COLn header keywords, or
    {} if this header doesn't have them."""
    col_map = {}
    n = 1
    while f'COL{n}' in header:
        col_map[str(header[f'COL{n}']).strip()] = n - 1
        n += 1
    return col_map


def _detect_graces(hdul):
    col_map = _graces_col_map(hdul[0].header)
    if not col_map:
        return False
    return (_GRACES_ORDER_LABEL in col_map
            and any(w in col_map for w in _GRACES_WAVE_PREFERENCE)
            and any(f in col_map for f in _GRACES_FLUX_PREFERENCE))


def _read_graces(hdul, **kwargs):
    header = hdul[0].header
    data = hdul[0].data
    col_map = _graces_col_map(header)

    order_row = data[col_map[_GRACES_ORDER_LABEL]]
    wave_label = next(w for w in _GRACES_WAVE_PREFERENCE if w in col_map)
    flux_label = next(f for f in _GRACES_FLUX_PREFERENCE if f in col_map)
    wave_row = data[col_map[wave_label]]
    flux_row = data[col_map[flux_label]]
    print(f"Detected GRACES/OPERA format -- wavelength column '{wave_label}' "
          f"(nm, converted to Angstrom), flux column '{flux_label}'")

    orders = np.unique(order_row)
    wavelength = np.empty(len(orders), dtype=object)
    flux = np.empty(len(orders), dtype=object)
    for i, o in enumerate(orders):
        mask = order_row == o
        wavelength[i] = wave_row[mask].astype(float) * 10.0  # nm -> Angstrom
        flux[i] = flux_row[mask].astype(float)
    return wavelength, flux, None


# ---------------------------------------------------------------------------
# FITS: GRACES / CFHT OPERA pipeline -- merged "intensity" (.i) format
# ---------------------------------------------------------------------------
# A different OPERA output variant from the per-order "spectrum" (.m) format
# above -- confirmed against a real file: no 'Order' column at all, every
# echelle order already merged into one continuous ~194000-point array. The
# file offers FOUR wavelength/intensity/errorbar column triplets (COL1-3,
# 4-6, 7-9, 10-12) -- Normalized/UnNormalized crossed with autowave-
# corrected/uncorrected -- distinguished only by each COLn's header COMMENT
# (e.g. "UnNormalized, no autowave correction"), since the VALUES are
# identical repeated 'Wavelength'/'Intensity'/'ErrorBar' labels. Per this
# module's "prefer raw/uncorrected" design principle, always uses the
# UnNormalized, no-autowave-correction triplet -- found by matching each
# triplet's comment text, not a hardcoded column position, in case a future
# file orders these differently.
#
# Since there's no Order column, orders are recovered by detecting
# wavelength BOUNDARIES -- the less-robust heuristic the .m reader above
# deliberately avoids in favor of its explicit Order column, used here only
# because this format gives us nothing better. A boundary is either a
# negative jump (adjacent orders OVERLAP, common at bluer wavelengths --
# confirmed: -2.7 to -3.7 A against a normal ~0.0026 A point spacing, three
# orders of magnitude of contrast) OR an abnormally large POSITIVE jump
# relative to the local point spacing (adjacent orders have a real
# wavelength GAP instead, confirmed to be how redder orders behave in a
# real file -- coverage stops overlapping past ~820 nm, so a negative-jump-
# only check silently merged the reddest ~5 orders into one, 39101-point,
# 2270-A-wide "order" at ~22x coarser effective sampling than everywhere
# else). Comparing each diff against a ROLLING local median (not a single
# global threshold) is needed either way, since real point spacing itself
# grows gradually and smoothly from blue to red across the whole spectrum.
# Confirmed against a real file: recovers exactly 35 orders, matching this
# instrument's known-good order count from the .m format's explicit Order
# column, with smoothly increasing order sizes throughout (no more merged
# segments).

_GRACES_INTENSITY_WANT_COMMENT_FRAGMENTS = ['unnormalized', 'no autowave']


def _graces_intensity_col_groups(header):
    """Return the (wave_idx, flux_idx, err_idx) row indices for the
    'UnNormalized, no autowave correction' triplet, found via each COLn's
    COMMENT text (not its value -- see module note above) -- or None if
    this header doesn't have the expected repeating-triplet structure."""
    n = 1
    labels = []
    while f'COL{n}' in header:
        labels.append((str(header[f'COL{n}']).strip(),
                        str(header.comments[f'COL{n}']).strip().lower()))
        n += 1
    if len(labels) < 3 or len(labels) % 3 != 0:
        return None
    for i in range(0, len(labels), 3):
        (wave_label, wave_comment), (flux_label, _), (err_label, _) = labels[i:i + 3]
        if (wave_label == 'Wavelength' and flux_label == 'Intensity' and err_label == 'ErrorBar'
                and all(frag in wave_comment for frag in _GRACES_INTENSITY_WANT_COMMENT_FRAGMENTS)):
            return i, i + 1, i + 2
    return None


def _detect_graces_intensity(hdul):
    header = hdul[0].header
    if _GRACES_ORDER_LABEL in _graces_col_map(header):
        return False  # the per-order '.m' format above already handles this
    return _graces_intensity_col_groups(header) is not None


def _read_graces_intensity(hdul, **kwargs):
    header = hdul[0].header
    data = hdul[0].data
    wave_idx, flux_idx, _ = _graces_intensity_col_groups(header)
    wave_row = data[wave_idx].astype(float) * 10.0  # nm -> Angstrom
    flux_row = data[flux_idx].astype(float)
    print("Detected GRACES/OPERA merged-intensity format -- using the "
          "UnNormalized, no-autowave-correction column triplet; splitting "
          "back into per-order chunks via wavelength-boundary detection (no "
          "Order column available in this format)")

    diffs = np.diff(wave_row)
    local_spacing = median_filter(diffs, size=101, mode='nearest')
    boundary = (diffs < 0) | (diffs > 5.0 * local_spacing)
    edges = np.concatenate(([0], np.where(boundary)[0] + 1, [len(wave_row)]))

    wavelength = np.empty(len(edges) - 1, dtype=object)
    flux = np.empty(len(edges) - 1, dtype=object)
    for i in range(len(edges) - 1):
        sl = slice(edges[i], edges[i + 1])
        wavelength[i] = wave_row[sl]
        flux[i] = flux_row[sl]
    return wavelength, flux, None


# ---------------------------------------------------------------------------
# FITS: binary table with named wave/flux columns -- e.g. this project's own
# KOA HIRES Level-1 consolidation output (Verification/Brewer2016_HIRES)
# ---------------------------------------------------------------------------

_TABLE_WAVE_NAMES = ['wave', 'wavelength']
_TABLE_FLUX_NAMES = ['flux', 'counts']


def _table_hdus(hdul):
    return [h for h in hdul[1:] if isinstance(h, (fits.BinTableHDU, fits.TableHDU))]


def _detect_binary_table(hdul):
    tables = _table_hdus(hdul)
    if not tables:
        return False
    names = [c.lower() for c in tables[0].columns.names]
    return (any(w in names for w in _TABLE_WAVE_NAMES)
            and any(f in names for f in _TABLE_FLUX_NAMES))


def _read_binary_table(hdul, **kwargs):
    tables = _table_hdus(hdul)
    names = [c.lower() for c in tables[0].columns.names]
    wave_name = tables[0].columns.names[next(i for i, n in enumerate(names) if n in _TABLE_WAVE_NAMES)]
    flux_name = tables[0].columns.names[next(i for i, n in enumerate(names) if n in _TABLE_FLUX_NAMES)]
    print(f"Detected FITS binary-table format -- wavelength column '{wave_name}', "
          f"flux column '{flux_name}' ({len(tables)} order/extension(s))")

    wavelength = np.empty(len(tables), dtype=object)
    flux = np.empty(len(tables), dtype=object)
    for i, t in enumerate(tables):
        wavelength[i] = np.asarray(t.data[wave_name], dtype=float)
        flux[i] = np.asarray(t.data[flux_name], dtype=float)
    return wavelength, flux, None


# ---------------------------------------------------------------------------
# FITS: classic Keck/MAKEE multi-order image
# ---------------------------------------------------------------------------
# Wavelength solution either per-order linear WCS (CRVL1_NN/CDLT1_NN header
# keys) or an IRAF-style multispec blob (WAT2_NNN keys), decoded by
# _wat_info(). Also the only format here with a usable gain estimate.

def _detect_keck_makee(hdul):
    header = hdul[0].header
    if header.get('NAXIS', 0) < 2:
        return False
    has_crvl = 'CRVL1_01' in header or 'CRVL1_1' in header
    has_wat = any(k.startswith('WAT2_') for k in header)
    return has_crvl or has_wat


def _wat_info(hdul):
    """Starting wavelength + spacing per order from IRAF-style WAT2_NNN
    multispec header keywords."""
    header = hdul[0].header
    spectrum_information = ''
    for key in header.keys():
        if 'WAT' in key:
            spectrum_information += header[key]

    min_wave = 3000
    starting_wavs = np.zeros(header['NAXIS2'])
    wave_spacing = np.zeros(header['NAXIS2'])
    count = 0
    for stuff in spectrum_information.split('spec'):
        if ' = ' in stuff:
            for i, item in enumerate(stuff.split(' ')):
                try:
                    floats = float(item)
                    spacing = 0.0
                    starting_wave = 0.0
                    if floats > min_wave:
                        starting_wave = floats
                        spacing = float(stuff.split(' ')[i + 1])
                        break
                except ValueError:
                    pass
            starting_wavs[count] = starting_wave
            wave_spacing[count] = spacing
            count += 1
    return starting_wavs, wave_spacing


def _keck_chip_gain(header, num_orders):
    """Per-chip CCD gain from Keck CCDGAIN/starting-wavelength header info
    (https://www2.keck.hawaii.edu/inst/hires/ccdgain.html)."""
    spec_info = header['WAT2_001'].split()
    num = None
    for item in spec_info:
        try:
            candidate = float(item)
            if candidate / 100 > 1.0:
                num = candidate
                break
        except ValueError:
            pass
    starting_wavelength = num
    low_high = header['CCDGAIN']
    if starting_wavelength < 5000.0:
        gain = 1.95 if low_high == 'low' else 0.78
    elif starting_wavelength < 6500.0:
        gain = 2.09 if low_high == 'low' else 0.84
    else:
        gain = 2.09 if low_high == 'low' else 0.89
    return np.ones(num_orders) * gain


def _read_keck_makee(hdul, **kwargs):
    print("Detected classic Keck/MAKEE multi-order format")
    header = hdul[0].header
    num_orders = header['NAXIS2']
    num_points = header['NAXIS1']

    wavelength = np.zeros((num_orders, num_points))

    try:
        gain = np.ones(num_orders) * header['CCDGN01']
    except KeyError:
        gain = _keck_chip_gain(header, num_orders)

    try:
        for i in range(num_orders):
            suffix = ('0' + str(i + 1)) if i + 1 < 10 else str(i + 1)
            start_wave_key = 'CRVL1_' + suffix
            spacing_wave_key = 'CDLT1_' + suffix
            wavelength[i][0] = header[start_wave_key]
            for j in range(1, num_points):
                wavelength[i][j] = wavelength[i][j - 1] + header[spacing_wave_key]
    except KeyError:
        starting_waves, wave_spacing = _wat_info(hdul)
        for i in range(num_orders):
            wavelength[i][0] = starting_waves[i]
            for j in range(1, num_points):
                wavelength[i][j] = wavelength[i][j - 1] + wave_spacing[i]

    flux = hdul[0].data
    return wavelength, flux, gain


# ---------------------------------------------------------------------------
# HDF5: MAROON-X
# ---------------------------------------------------------------------------
# spec_blue/spec_red are pandas DataFrames indexed by (Fiber, Order), with
# 'wavelengths' (nm) and 'optimal_extraction' (raw extracted counts, NOT
# continuum-normalized -- confirmed by inspecting a real file) columns.
#
# CAVEAT: which Fiber index is the science target is NOT reliably
# determinable from file metadata alone (checked: the header's FIBER1..5
# fields describe a different numbering than the Fiber index actually used
# in the DataFrame). Defaults to fiber 6, matching prior working usage
# (Verification/.../maroonx_xspect_ews.ipynb) -- pass fiber= to override if
# that's wrong for a given file.

_MAROONX_DEFAULT_FIBER = 6
_MAROONX_DEFAULT_FLUX_COLUMN = 'optimal_extraction'
_MAROONX_BANDS = ('spec_blue', 'spec_red')


def _is_hdf5(filename):
    try:
        with open(filename, 'rb') as f:
            return f.read(8) == b'\x89HDF\r\n\x1a\n'
    except OSError:
        return False


def _detect_maroonx(store):
    keys = set(store.keys())
    return {'/spec_blue', '/spec_red'}.issubset(keys)


def _read_maroonx(store, fiber=_MAROONX_DEFAULT_FIBER,
                   flux_column=_MAROONX_DEFAULT_FLUX_COLUMN, **kwargs):
    print(f"Detected MAROON-X format -- fiber={fiber} (ASSUMED science fiber, "
          f"not verifiable from file metadata -- pass fiber= to override if "
          f"wrong), flux column='{flux_column}'")
    wave_list = []
    flux_list = []
    for band in _MAROONX_BANDS:
        spec = store[band]
        wave_col = spec['wavelengths'].loc[fiber]
        flux_col = spec[flux_column].loc[fiber]
        for order in wave_col.index:
            wave_list.append(np.asarray(wave_col.loc[order], dtype=float) * 10.0)  # nm -> Angstrom
            flux_list.append(np.asarray(flux_col.loc[order], dtype=float))

    wavelength = np.empty(len(wave_list), dtype=object)
    flux = np.empty(len(flux_list), dtype=object)
    for i in range(len(wave_list)):
        wavelength[i] = wave_list[i]
        flux[i] = flux_list[i]
    return wavelength, flux, None


def get_maroonx_bands(filename, fiber=_MAROONX_DEFAULT_FIBER):
    """
    Return a list of 'blue'/'red' arm labels, one per order, in the same
    order _read_maroonx() (and therefore Spectrum_Data's wavelength/flux
    arrays) produces them -- all spec_blue orders first, then all
    spec_red orders.

    Needed for response correction specifically: MAROON-X's two arms
    physically overlap in wavelength near their dichroic split (confirmed
    against real files: spec_blue independently covers echelle orders
    91-124 and spec_red covers 67-94 -- both arms cover orders 91-94), so
    a science order in that overlap can look like a wavelength match for
    TWO different response chunks (one per arm) that describe different
    physical light paths. Wavelength overlap alone can't disambiguate
    them (their wavelength ranges are nearly identical); only knowing
    which arm the science data and the response chunk each came from can.
    See response_correction.apply_response_correction()'s band-matching.
    """
    import pandas as pd
    bands = []
    with pd.HDFStore(filename, 'r') as store:
        for band in _MAROONX_BANDS:
            spec = store[band]
            n = len(spec['wavelengths'].loc[fiber].index)
            bands.extend([band.replace('spec_', '')] * n)
    return bands


def load_maroonx_response(filename):
    """
    Load a MAROON-X PHOENIX-based instrument response/blaze correction
    file (e.g. MAROON-X_PHOENIX_RESPONSE_CORRECTORDERS.hd5) -- a SEPARATE
    calibration file, not embedded in individual science exposures
    (confirmed: the per-exposure 'blaze_blue'/'blaze_red' keys exist but
    are empty in real science files). Returns per-order (wave, response)
    arrays in the same list-of-arrays convention Spectrum_Data itself
    uses, ready to pass to Spectrum_Data.apply_response_correction(),
    plus a parallel list of 'blue'/'red' arm labels (see get_maroonx_bands()
    for why matching needs this, not just wavelength overlap).

    File structure (confirmed against a real file, not assumed): keys
    'wavelength_blue'/'wavelength_red' and 'response_blue'/'response_red',
    each a DataFrame with one column per echelle order. The two DataFrames'
    columns do NOT share label values -- wavelength_* uses plain positional
    labels (0, 1, 2, ...) while response_* uses the real echelle order
    numbers (e.g. 92, 93, ...) -- but they correspond POSITIONALLY
    (wavelength column 0 <-> response column 92, etc.), confirmed by
    checking that the wavelength ranges line up when paired that way.

    Returns
    -------
    wave_list, resp_list : as before
    band_list : 'blue'/'red' per chunk, parallel to wave_list/resp_list
    """
    import pandas as pd
    wave_list = []
    resp_list = []
    band_list = []
    with pd.HDFStore(filename, 'r') as store:
        for band in ('blue', 'red'):
            wave_df = store[f'wavelength_{band}']
            resp_df = store[f'response_{band}']
            for wave_col, resp_col in zip(wave_df.columns, resp_df.columns):
                wave_list.append(np.asarray(wave_df[wave_col], dtype=float) * 10.0)  # nm -> Angstrom
                resp_list.append(np.asarray(resp_df[resp_col], dtype=float))
                band_list.append(band)
    return wave_list, resp_list, band_list


# ---------------------------------------------------------------------------
# FITS: Gemini GHOST, DRAGONS orders-combined ("_dragons") 1D output
# ---------------------------------------------------------------------------
# Verified against real Gemini-S GHOST blue- and red-arm files (Standard
# resolution, DRAGONS combineOrders with scale='loglinear'): every echelle
# order is already merged into ONE continuous 1D spectrum per IFU, stored as
# SCI/VAR/DQ ImageHDU triplets (one triplet per IFU, EXTVER 1 and 2). The
# wavelength solution is a log-linear WCS (CTYPE1='AWAV-LOG', CUNIT1='nm'):
#     lambda = CRVAL1 * exp(CD1_1 * (p - CRPIX1) / CRVAL1)
# with p the 1-based FITS pixel -- NOT the linear CRVAL1 + CD1_1*(p-CRPIX1),
# which is off by ~10 nm at the ends of the blue arm. Evaluated via
# astropy.wcs rather than hand-coded.
#
# Which IFU holds the star: TARGETn in the primary header (2 = object, per
# every real file checked, where the object IFU was always ~30-60x brighter
# than the other; 0 = unused/sky). Falls back to the brightest SCI
# extension if no TARGETn==2 is present; ifu= overrides either.
#
# Departures from this module's "prefer raw" principle, forced by the
# format (DRAGONS offers no alternative in this product): wavelengths are
# ALREADY barycentric-corrected (barycentricCorrect step in the HISTORY
# table), and flux is ADU-scaled by order stacking rather than pure counts
# (VAR/flux ~0.33 in the bulk, up to 100x+ at the noisy blue end) --
# so sqrt(flux) is a poor error estimate and the reader returns sqrt(VAR)
# instead (see read_spectrum()'s return_err).
#
# Since there are no orders left to recover, the merged spectrum is cut
# into fixed-VELOCITY-width chunks on a GLOBAL ln(lambda) grid (edges at
# exp(k * chunk_kms / c), independent of each file's own start wavelength
# or log step, which differ slightly between reductions) -- so every
# exposure of the same arm is cut at identical wavelengths and
# combine_spectra() pairs matching chunks across epochs. Adjacent chunks
# overlap by 2 * pad_kms, so a line near one chunk's artificial boundary
# (where fit_als_continuum()'s edge_ignore_aa leaves the continuum
# extrapolated) sits comfortably inside its neighbour, and measure_all_
# ew()'s order-overlap handling keeps the measurement furthest from an edge.
#
# DQ-flagged pixels (flux exactly 0 in every real file checked): a leading/
# trailing run is trimmed off; isolated interior runs (1-4 px) have flux
# linearly interpolated across them and their error inflated to
# _GHOST_BAD_ERR_FACTOR x the chunk's median error, so they carry ~no
# weight without introducing NaN/inf into the continuum fit.

_C_KMS = 299792.458
_GHOST_BAD_ERR_FACTOR = 1e3

# Spectrum_Data.normalize_all() defaults for DRAGONS output. Its blaze is
# already removed, leaving only a smooth star x throughput trend, so a
# single STIFF fit over the whole merged arm replaces the flexible,
# adaptively-stiffened per-order fit tuned for raw echelle orders
# ('merged': fit once across all chunks, then split back -- fitting each
# artificial chunk separately left +/-10% continuum jumps wherever a chunk
# edge cut through a broad feature: Hgamma/G band, Hbeta, Ca H+K).
# lam=2e6 (~13 A smoothing length on GHOST's ~0.023 A/px grid) chosen over
# 2e5..2e8 on real data: 2e7 sat ~25% under the continuum peaks blueward of
# Ca H+K and drifted 4% high at the red end; 2e8 couldn't follow the real
# trend at all. calib_window_aa: see fit_als_continuum().
FORMAT_NORMALIZE_DEFAULTS = {
    'ghost': dict(lam=2e6, adaptive=False, calib_window_aa=50.0, merged=True),
}


def _ghost_sci_hdus(hdul):
    return [h for h in hdul[1:]
            if isinstance(h, fits.ImageHDU) and h.name == 'SCI' and h.data is not None]


def _detect_ghost(hdul):
    if str(hdul[0].header.get('INSTRUME', '')).strip().upper() != 'GHOST':
        return False
    sci = _ghost_sci_hdus(hdul)
    return bool(sci) and all(h.data.ndim == 1 for h in sci)


def _read_ghost(hdul, ifu='auto', chunk_kms=5000.0, pad_kms=200.0, **kwargs):
    """
    ifu : 'auto' (default; TARGETn==2, else brightest) or 1/2 to force a
        specific IFU's SCI/VAR/DQ extensions (matched on EXTVER).
    chunk_kms : chunk width in km/s (5000 km/s ~ 3100 px ~ 70 A at 4350 A,
        close to this package's usual ~3500-points-per-order target), or
        None to return the whole merged spectrum as a single "order".
    pad_kms : each chunk extends this far past its grid edges on both sides.
    """
    from astropy.wcs import WCS

    header = hdul[0].header
    sci = _ghost_sci_hdus(hdul)
    if ifu == 'auto':
        targets = [h.ver for h in sci if header.get(f'TARGET{h.ver}') == 2]
        if targets:
            ifu, why = targets[0], f'TARGET{targets[0]}=2'
        else:
            ifu = max(sci, key=lambda h: np.nanmedian(h.data)).ver
            why = 'brightest SCI extension (no TARGETn=2 in header)'
    else:
        why = 'requested'
    sci_hdu = next((h for h in sci if h.ver == ifu), None)
    if sci_hdu is None:
        raise ValueError(f"GHOST file has no SCI extension for IFU {ifu} "
                         f"(available: {[h.ver for h in sci]})")

    flux = np.asarray(sci_hdu.data, dtype=float)
    var = np.asarray(hdul['VAR', ifu].data, dtype=float)
    dq = np.asarray(hdul['DQ', ifu].data)

    wcs = WCS(sci_hdu.header)
    wave = np.asarray(wcs.pixel_to_world_values(np.arange(len(flux))), dtype=float)
    unit = str(sci_hdu.header.get('CUNIT1', 'nm')).strip().lower()
    # astropy.wcs normalizes spectral axes to SI (m); convert from whatever
    # it actually returned rather than trusting CUNIT1 alone.
    if np.nanmedian(wave) < 1e-3:
        wave = wave * 1e10                     # m -> Angstrom
    elif unit == 'nm':
        wave = wave * 10.0                     # nm -> Angstrom

    bad = (dq != 0) | ~np.isfinite(flux) | ~np.isfinite(var) | (var <= 0)
    good_idx = np.where(~bad)[0]
    if len(good_idx) == 0:
        raise ValueError(f"GHOST IFU {ifu}: every pixel is DQ-flagged/non-finite")
    keep = slice(good_idx[0], good_idx[-1] + 1)  # trim leading/trailing bad runs
    wave, flux, var, bad = wave[keep], flux[keep], var[keep], bad[keep]
    n_interior_bad = int(bad.sum())
    if n_interior_bad:
        flux[bad] = np.interp(wave[bad], wave[~bad], flux[~bad])
    err = np.sqrt(np.where(bad, np.nan, var))

    print(f"Detected Gemini GHOST DRAGONS 1D format ({header.get('CAMERA', sci_hdu.header.get('CAMERA', '?'))} arm) "
          f"-- IFU {ifu} ({why}), log-linear wavelength {wave[0]:.2f}-{wave[-1]:.2f} A "
          f"(already barycentric-corrected by DRAGONS), errors from VAR; "
          f"trimmed {keep.start} leading / {len(dq) - keep.stop} trailing bad px, "
          f"interpolated over {n_interior_bad} interior bad px")

    lnw = np.log(wave)
    if chunk_kms is None:
        sels = [np.ones(len(lnw), dtype=bool)]
    else:
        step = np.log1p(chunk_kms / _C_KMS)
        pad = np.log1p(pad_kms / _C_KMS)
        k0, k1 = int(np.floor(lnw[0] / step)), int(np.ceil(lnw[-1] / step))
        sels = [(lnw >= k * step - pad) & (lnw < (k + 1) * step + pad) for k in range(k0, k1)]
        sels = [sel for sel in sels if sel.any()]
        # fold a sliver at either end of the spectrum (under a quarter of a
        # full chunk) into its neighbour rather than dropping its pixels
        min_px = 0.25 * step / np.median(np.diff(lnw))
        if len(sels) > 1 and sels[0].sum() < min_px:
            sels = [sels[0] | sels[1]] + sels[2:]
        if len(sels) > 1 and sels[-1].sum() < min_px:
            sels = sels[:-2] + [sels[-2] | sels[-1]]
    wavelength, flux_out, err_out = [], [], []
    for sel in sels:
        e = err[sel]
        e = np.where(np.isfinite(e), e, _GHOST_BAD_ERR_FACTOR * np.nanmedian(e))
        wavelength.append(wave[sel])
        flux_out.append(flux[sel])
        err_out.append(e)

    def as_obj(lst):
        out = np.empty(len(lst), dtype=object)
        for i, a in enumerate(lst):
            out[i] = a
        return out
    return as_obj(wavelength), as_obj(flux_out), None, as_obj(err_out)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_FITS_FORMATS = [
    (_detect_graces, _read_graces, 'graces'),
    (_detect_graces_intensity, _read_graces_intensity, 'graces_intensity'),
    (_detect_ghost, _read_ghost, 'ghost'),
    (_detect_keck_makee, _read_keck_makee, 'keck_hires'),
    (_detect_binary_table, _read_binary_table, 'fits_table'),
]

_NAMED_READERS = {
    'graces': _read_graces,
    'graces_intensity': _read_graces_intensity,
    'ghost': _read_ghost,
    'keck_hires': _read_keck_makee,
    'fits_table': _read_binary_table,
    'maroonx': _read_maroonx,
}


def detect_format(filename, instrument='auto', custom_reader=None):
    """Name of the reader read_spectrum() would use for this file ('ghost',
    'keck_hires', ...), or None for a custom_reader or an unrecognized file."""
    if custom_reader is not None:
        return None
    if instrument != 'auto':
        return instrument
    try:
        with fits.open(filename) as hdul:
            for detect, _, name in _FITS_FORMATS:
                if detect(hdul):
                    return name
    except OSError:
        pass
    if _is_hdf5(filename):
        import pandas as pd
        with pd.HDFStore(filename, 'r') as store:
            if _detect_maroonx(store):
                return 'maroonx'
    return None


def read_spectrum(filename, instrument='auto', custom_reader=None, return_err=False, **kwargs):
    """
    Read a spectrum file, auto-detecting its format.

    Parameters
    ----------
    filename : str
    instrument : 'auto' (default) or one of 'graces', 'graces_intensity',
        'ghost', 'keck_hires', 'fits_table', 'maroonx' to skip detection and
        force a specific reader.
    custom_reader : callable(filename, **kwargs) -> (wavelength, flux),
        (wavelength, flux, gain) or (wavelength, flux, gain, err), optional
        -- for a one-off format without editing this file. wavelength/flux
        are object arrays of per-order ndarrays (Angstrom, counts).
    return_err : if True, also return per-order 1-sigma flux errors (the
        fourth element) -- only formats that ship a real error/variance
        array provide one (currently 'ghost'); None otherwise, in which case
        callers fall back to their own sqrt(flux) estimate. False (default)
        keeps the original 3-tuple return for existing callers.

    Returns
    -------
    wavelength, flux : object arrays of per-order ndarrays
    gain : ndarray or None
    err : object array of per-order ndarrays, or None -- only if return_err
    """
    result = _read_spectrum_raw(filename, instrument=instrument,
                                custom_reader=custom_reader, **kwargs)
    result = tuple(result) + (None,) * (4 - len(result))
    return result if return_err else result[:3]


def _read_spectrum_raw(filename, instrument='auto', custom_reader=None, **kwargs):
    """read_spectrum()'s dispatch; returns whatever tuple (2-4 elements) the
    matched reader produced."""
    if custom_reader is not None:
        return custom_reader(filename, **kwargs)

    if instrument != 'auto':
        if instrument not in _NAMED_READERS:
            raise ValueError(f"Unknown instrument '{instrument}'; "
                              f"choose from {list(_NAMED_READERS)} or 'auto'")
        if instrument == 'maroonx':
            import pandas as pd
            with pd.HDFStore(filename, 'r') as store:
                return _read_maroonx(store, **kwargs)
        with fits.open(filename) as hdul:
            return _NAMED_READERS[instrument](hdul, **kwargs)

    tried = []
    header_dump = None
    try:
        with fits.open(filename) as hdul:
            for detect, read, name in _FITS_FORMATS:
                tried.append(name)
                if detect(hdul):
                    return read(hdul, **kwargs)
            header_dump = repr(hdul[0].header)
    except OSError:
        pass

    if _is_hdf5(filename):
        import pandas as pd
        tried.append('maroonx')
        with pd.HDFStore(filename, 'r') as store:
            if _detect_maroonx(store):
                return _read_maroonx(store, **kwargs)

    msg = (f"Could not auto-detect the format of '{filename}'.\n"
           f"Tried: {tried}.\n")
    if header_dump is not None:
        msg += f"FITS primary header for reference:\n{header_dump}\n"
    msg += (
        "This spectrum uses a format not yet recognized. Two ways to "
        "proceed without waiting on a SPAE change:\n"
        "  1. Extract (wavelength, flux) arrays yourself and pass them via "
        "Spectrum_Data(name, KECK_file=False, spectx=..., specty=...).\n"
        "  2. Pass a custom_reader=<callable> to Spectrum_Data()/read_spectrum() "
        "for a one-off format.\n"
        "If this format will recur, add a new detector/extractor pair to "
        "SPAE/xspect_ew/readers.py (see the existing entries for the pattern) "
        "so future files are recognized automatically."
    )
    raise SpectrumFormatError(msg)
