"""
XSpect-EW: user-guided equivalent-width measurement for stellar absorption
spectra. Originally developed by George Vejar (github.com/forgeousgeorge/XSpect),
incorporated into SPAE 2026-09-14.

Loads a spectrum -- auto-detecting its format (Keck/MAKEE, GRACES/OPERA,
MAROON-X, a FITS binary table with named wave/flux columns, or raw
wavelength/flux arrays; see readers.py) -- fits the continuum per echelle
order via Asymmetric Least Squares smoothing, optionally cross-correlates
against a reference spectrum to solve for wavelength shift/RV, measures
each line's equivalent width, and writes a MOOG-format linelist-with-EW
file directly usable by SPAE.moog.abfind.
"""

from .spectrum_data import Spectrum_Data
from .continuum import fit_als_continuum
from .constants import ELEMENTS
from .plotting import plot_line_info, plot_comparison_res
from .combine import combine_files
from .io_utils import load_object
from .readers import read_spectrum, SpectrumFormatError, load_maroonx_response, get_maroonx_bands
from .radial_velocity import measure_effective_rv, measure_rv_from_linelist, RV_REFERENCE_LINES
from .response_correction import apply_response_correction
from .overlap_check import check_order_overlaps, print_overlap_report, flagged_overlap_ranges
from .line_identification import identify_line, identify_lines_in_spectrum
