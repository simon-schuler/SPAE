"""
Matplotlib-based spectral plotting — replaces PGPLOT/X11 output from MOOG.

Functions
---------
plot_spectrum   — plot synthetic spectrum; optionally overlay observed
read_obs        — read a two-column (wavelength, flux) observed spectrum
"""
import numpy as np


def read_obs(filename: str) -> tuple:
    """
    Read a MONGO-style two-column observed spectrum (whitespace-separated).

    First line is treated as a comment/title and skipped if it cannot be
    parsed as two floats.  Wavelengths are sorted ascending.

    Returns
    -------
    wave : np.ndarray [Å]
    flux : np.ndarray (normalised to continuum)
    """
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
