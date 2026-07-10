"""
Matplotlib-based spectral plotting — replaces PGPLOT/X11 output from MOOG.

Functions
---------
plot_spectrum   — plot synthetic spectrum; optionally overlay observed
read_obs        — read observed spectrum (text or FITS)
"""
import numpy as np

# Column name candidates, checked in order (first match wins)
_WAVE_COLS = ['WAVELENGTH', 'WAVE', 'LAMBDA', 'LOGLAM', 'WAVE_VAC', 'WAVE_AIR',
              'WAV', 'WAVEL']
_FLUX_COLS = ['FLUX', 'NORMALIZED_FLUX', 'NORM_FLUX', 'FLUX_NORM', 'NFLUX',
              'SPEC', 'SPECTRUM', 'DATA']


def _to_angstrom(wave: np.ndarray, unit: str) -> np.ndarray:
    """Convert wavelength array to Å based on unit string."""
    u = unit.strip().lower()
    if u in ('nm', 'nanometer', 'nanometers'):
        return wave * 10.0
    if u in ('um', 'micron', 'microns', 'micrometer', 'micrometers'):
        return wave * 1e4
    return wave   # already Å (or unknown — leave as-is)


def _read_obs_fits(filename: str) -> tuple:
    """
    Read a 1-D spectrum from a FITS file.

    Tries, in order:
      1. Binary or ASCII table extension — looks for wavelength and flux
         columns by name (case-insensitive).  LOGLAM columns are treated as
         log10(Å).  Column units (nm, µm) are converted to Å.
      2. Image extension (primary or first IMAGE) — reconstructs the
         wavelength grid from the WCS keywords CRVAL1, CDELT1/CD1_1, CRPIX1.
         CTYPE1 containing 'LOG' triggers 10^wave conversion.
         CUNIT1 'nm' or 'µm' triggers unit conversion to Å.

    Returns
    -------
    wave : np.ndarray [Å]
    flux : np.ndarray
    """
    try:
        from astropy.io import fits
    except ImportError:
        raise ImportError(
            "astropy is required to read FITS spectra.  "
            "Install it with:  pip install astropy"
        )

    with fits.open(filename) as hdul:
        # ---- 1. Table extensions ----------------------------------------
        for hdu in hdul:
            if not hasattr(hdu, 'columns') or hdu.data is None:
                continue
            col_names = [c.name.upper() for c in hdu.columns]
            wave_name = next((c for c in _WAVE_COLS if c in col_names), None)
            flux_name = next((c for c in _FLUX_COLS if c in col_names), None)
            if wave_name is None or flux_name is None:
                continue

            wave = np.asarray(hdu.data[wave_name], dtype=np.float64).ravel()
            flux = np.asarray(hdu.data[flux_name], dtype=np.float64).ravel()

            # log10-wavelength (e.g. SDSS LOGLAM)
            if wave_name == 'LOGLAM':
                wave = 10.0 ** wave
            else:
                col_unit = hdu.columns[wave_name].unit or ''
                wave = _to_angstrom(wave, col_unit)
            return wave, flux

        # ---- 2. Image HDU with WCS --------------------------------------
        for hdu in hdul:
            if hdu.data is None:
                continue
            data = np.asarray(hdu.data, dtype=np.float64).ravel()
            if len(data) < 2:
                continue
            hdr = hdu.header
            if 'CRVAL1' not in hdr:
                continue
            crval1 = float(hdr['CRVAL1'])
            crpix1 = float(hdr.get('CRPIX1', 1))
            cdelt1 = float(hdr.get('CDELT1', hdr.get('CD1_1', 1.0)))
            pixels = np.arange(1, len(data) + 1, dtype=np.float64)
            wave   = crval1 + (pixels - crpix1) * cdelt1
            ctype1 = hdr.get('CTYPE1', '')
            if 'LOG' in ctype1.upper():
                wave = 10.0 ** wave
            cunit1 = hdr.get('CUNIT1', 'Angstrom')
            wave = _to_angstrom(wave, cunit1)
            return wave, data

    raise ValueError(
        f"Cannot read FITS spectrum from '{filename}': no recognised table "
        "columns (WAVELENGTH/FLUX variants) or image WCS (CRVAL1) found."
    )


def read_obs(filename: str) -> tuple:
    """
    Read an observed spectrum from a text or FITS file.

    Text format (default)
    ---------------------
    Whitespace-separated wavelength and flux columns (MONGO/MOOG style).
    Comment lines starting with '#' and a single non-numeric header line
    are skipped.  Wavelengths are sorted ascending on return.

    FITS format
    -----------
    Triggered automatically for files ending in .fits, .fit, .fits.gz, or
    .fit.gz.  Supports:
      - Binary/ASCII table with WAVELENGTH (or WAVE/LAMBDA/LOGLAM/…) and
        FLUX (or NORMALIZED_FLUX/NORM_FLUX/…) columns.
      - 1-D image extension with WCS keywords (CRVAL1, CDELT1, CRPIX1).
    Wavelengths in nm or µm are converted to Å automatically.

    Returns
    -------
    wave : np.ndarray [Å]
    flux : np.ndarray (normalised to continuum)
    """
    _FITS_EXTS = ('.fits', '.fit', '.fits.gz', '.fit.gz')
    if any(filename.lower().endswith(ext) for ext in _FITS_EXTS):
        wave, flux = _read_obs_fits(filename)
        order = np.argsort(wave)
        return wave[order], flux[order]

    # ---- plain text -------------------------------------------------------
    wave, flux = [], []
    with open(filename) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                w, f = float(parts[0]), float(parts[1])
                wave.append(w)
                flux.append(f)
            except ValueError:
                if i == 0:
                    continue   # skip title line
                raise

    wave = np.array(wave, dtype=np.float64)
    flux = np.array(flux, dtype=np.float64)
    if len(wave) > 1 and wave[1] < wave[0]:
        wave = wave[::-1]
        flux = flux[::-1]
    return wave, flux


def plot_spectrum(result: dict,
                 obs_wave=None, obs_flux=None,
                 title: str = '',
                 obs_label: str = 'Observed',
                 syn_label: str = 'Synthetic',
                 raw_label: str = 'Raw synthesis',
                 show_raw: bool = False,
                 xlim=None, ylim=None,
                 figsize=(12, 5),
                 syn_color: str = 'steelblue',
                 obs_color: str = 'black',
                 raw_color: str = 'lightblue',
                 ax=None):
    """
    Plot a synthetic spectrum from a synth() result dict.

    Parameters
    ----------
    result    : dict returned by synth() or synth_from_files()
    obs_wave  : array-like, observed wavelengths [Å]; None = no overlay
    obs_flux  : array-like, observed normalised flux; None = no overlay
    title     : plot title string
    obs_label : legend label for observed spectrum
    syn_label : legend label for smoothed synthetic spectrum
    raw_label : legend label for raw synthetic spectrum
    show_raw  : if True, also plot the unsmoothed depths
    xlim      : (wmin, wmax) in Å; None = auto
    ylim      : (fmin, fmax); None = auto
    figsize   : matplotlib figure size
    syn_color : colour for smoothed synthetic spectrum
    obs_color : colour for observed spectrum
    raw_color : colour for raw synthetic spectrum
    ax        : existing matplotlib Axes to draw into; None creates a new figure

    Returns
    -------
    fig, ax   : matplotlib Figure and Axes objects
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    wave = result['wave']
    flux_smooth = result['flux_smooth']
    flux_raw    = result.get('flux_raw', flux_smooth)

    if show_raw:
        ax.plot(wave, flux_raw, color=raw_color, lw=0.8,
                label=raw_label, zorder=1)

    ax.plot(wave, flux_smooth, color=syn_color, lw=1.2,
            label=syn_label, zorder=2)

    if obs_wave is not None and obs_flux is not None:
        ax.plot(obs_wave, obs_flux, color=obs_color, lw=0.8,
                label=obs_label, zorder=3)

    ax.set_xlabel('Wavelength (Å)', fontsize=11)
    ax.set_ylabel('Normalised flux', fontsize=11)
    if title:
        ax.set_title(title, fontsize=11)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(bottom=0.0)

    # Only show legend when there is more than one plotted element
    handles = ax.get_legend_handles_labels()
    if len(handles[0]) > 1:
        ax.legend(fontsize=10)

    fig.tight_layout()
    return fig, ax


def plot_abundance_grid(results: list,
                        labels=None,
                        obs_wave=None, obs_flux=None,
                        title: str = '',
                        figsize=(12, 5),
                        ax=None):
    """
    Overlay multiple synthetic spectra (abundance grid) on one plot.

    Parameters
    ----------
    results : list of synth() result dicts (one per abundance variation)
    labels  : list of strings; None → '(1)', '(2)', ...
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    colors = cm.viridis(np.linspace(0.15, 0.85, len(results)))

    for i, (res, col) in enumerate(zip(results, colors)):
        lbl = labels[i] if labels else f'({i + 1})'
        ax.plot(res['wave'], res['flux_smooth'], color=col, lw=1.2, label=lbl)

    if obs_wave is not None and obs_flux is not None:
        ax.plot(obs_wave, obs_flux, color='black', lw=0.8, label='Observed', zorder=99)

    ax.set_xlabel('Wavelength (Å)', fontsize=11)
    ax.set_ylabel('Normalised flux', fontsize=11)
    if title:
        ax.set_title(title, fontsize=11)
    ax.set_ylim(bottom=0.0)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig, ax
