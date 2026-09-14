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


def load_maroonx_response(filename):
    """
    Load a MAROON-X PHOENIX-based instrument response/blaze correction
    file (e.g. MAROON-X_PHOENIX_RESPONSE_CORRECTORDERS.hd5) -- a SEPARATE
    calibration file, not embedded in individual science exposures
    (confirmed: the per-exposure 'blaze_blue'/'blaze_red' keys exist but
    are empty in real science files). Returns per-order (wave, response)
    arrays in the same list-of-arrays convention Spectrum_Data itself
    uses, ready to pass to Spectrum_Data.apply_response_correction().

    File structure (confirmed against a real file, not assumed): keys
    'wavelength_blue'/'wavelength_red' and 'response_blue'/'response_red',
    each a DataFrame with one column per echelle order. The two DataFrames'
    columns do NOT share label values -- wavelength_* uses plain positional
    labels (0, 1, 2, ...) while response_* uses the real echelle order
    numbers (e.g. 92, 93, ...) -- but they correspond POSITIONALLY
    (wavelength column 0 <-> response column 92, etc.), confirmed by
    checking that the wavelength ranges line up when paired that way.
    """
    import pandas as pd
    wave_list = []
    resp_list = []
    with pd.HDFStore(filename, 'r') as store:
        for band in ('blue', 'red'):
            wave_df = store[f'wavelength_{band}']
            resp_df = store[f'response_{band}']
            for wave_col, resp_col in zip(wave_df.columns, resp_df.columns):
                wave_list.append(np.asarray(wave_df[wave_col], dtype=float) * 10.0)  # nm -> Angstrom
                resp_list.append(np.asarray(resp_df[resp_col], dtype=float))
    return wave_list, resp_list


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_FITS_FORMATS = [
    (_detect_graces, _read_graces, 'graces'),
    (_detect_keck_makee, _read_keck_makee, 'keck_hires'),
    (_detect_binary_table, _read_binary_table, 'fits_table'),
]

_NAMED_READERS = {
    'graces': _read_graces,
    'keck_hires': _read_keck_makee,
    'fits_table': _read_binary_table,
    'maroonx': _read_maroonx,
}


def read_spectrum(filename, instrument='auto', custom_reader=None, **kwargs):
    """
    Read a spectrum file, auto-detecting its format.

    Parameters
    ----------
    filename : str
    instrument : 'auto' (default) or one of 'graces', 'keck_hires',
        'fits_table', 'maroonx' to skip detection and force a specific reader.
    custom_reader : callable(filename, **kwargs) -> (wavelength, flux) or
        (wavelength, flux, gain), optional -- for a one-off format without
        editing this file. wavelength/flux are object arrays of per-order
        ndarrays (Angstrom, counts).

    Returns
    -------
    wavelength, flux : object arrays of per-order ndarrays
    gain : ndarray or None
    """
    if custom_reader is not None:
        result = custom_reader(filename, **kwargs)
        return result if len(result) == 3 else (result[0], result[1], None)

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
