"""
XSpect-EW: user-guided equivalent-width measurement for stellar absorption
spectra. Originally developed by George Vejar (github.com/forgeousgeorge/XSpect),
incorporated into SPAE 2026-09-14.

Loads a spectrum -- auto-detecting its format (Keck/MAKEE, GRACES/OPERA,
MAROON-X, a FITS binary table with named wave/flux columns, or raw
wavelength/flux arrays; see readers.py) -- fits the continuum per echelle
order via Gaussian Process regression, optionally cross-correlates against
a reference spectrum to solve for wavelength shift/RV, measures each
line's equivalent width, and writes a MOOG-format linelist-with-EW file
directly usable by SPAE.moog.abfind.
"""

from .spectrum_data import Spectrum_Data
from .continuum import Continuum_scan
from .constants import ELEMENTS
from .plotting import plot_line_info, plot_comparison_res
from .combine import combine_files
from .io_utils import load_object
from .readers import read_spectrum, SpectrumFormatError
from .radial_velocity import measure_effective_rv, RV_REFERENCE_LINES
