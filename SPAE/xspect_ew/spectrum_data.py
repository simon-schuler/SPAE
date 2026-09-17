"""Spectrum_Data: the core driver class for loading a spectrum, fitting its
continuum, wave-shifting against a reference, and measuring line EWs."""

import copy
import glob
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simpson
from numpy.random import multivariate_normal

from .constants import ELEMENTS
from .continuum import fit_als_continuum
from .line_profile import get_line_window, gauss_model, gfit_simple, gauss_ew
from .gp_utils import SEKernel, Pred_GP
from .combine import make_line, parabolic_refine, measure_order_alignment
from .plotting import make_plots_folder
from .readers import read_spectrum
from .radial_velocity import measure_effective_rv, measure_rv_from_linelist, C_KMS
from .response_correction import apply_response_correction as _apply_response_correction
from .overlap_check import check_order_overlaps as _check_order_overlaps
from .overlap_check import flagged_overlap_ranges as _flagged_overlap_ranges
from .outlier_check import check_bad_orders as _check_bad_orders
from .outlier_check import detect_spikes as _detect_spikes
from .outlier_check import check_cross_exposure_spikes as _check_cross_exposure_spikes
from .line_identification import identify_lines_in_spectrum as _identify_lines_in_spectrum


class Spectrum_Data():
    def __init__(self, filename, KECK_file = True, spectx=False, specty=False, order_split = (False,3500),
                 instrument = 'auto', custom_reader = None, **reader_kwargs):
        """
            filenmae - used to name plots and files
            KECK_file - True (default) to read `filename` from disk, auto-detecting
                its format (Keck/MAKEE, GRACES/OPERA, MAROON-X, or a FITS binary
                table with named wave/flux columns -- see readers.py). False to
                supply already-extracted arrays directly via spectx/specty instead.
            spectx, specty - input spectrum if loading from arrays (KECK_file=False)
                if not using order_split:
                    input format: [[order1],[order2],...] orderN = [x1,x2,x3,...]
                if using order_split:
                    input format: [x1,x2,x3,...]
            instrument - 'auto' (default), or force a specific reader by name
                ('graces', 'keck_hires', 'fits_table', 'maroonx') -- see
                readers.read_spectrum(). Only used when KECK_file=True.
            custom_reader - callable(filename, **kwargs) -> (wavelength, flux)
                or (wavelength, flux, gain), for a one-off format not yet
                recognized by readers.py. Only used when KECK_file=True.
            reader_kwargs - passed through to the reader (e.g. fiber=... for
                MAROON-X to override the assumed science fiber).
        """
        self.filename = filename
        self.gain = None
        if KECK_file:
            self.wavelength, self.flux, self.gain = read_spectrum(
                filename, instrument=instrument, custom_reader=custom_reader, **reader_kwargs)
        else:
            #split input array into orders
            if order_split[0]:
                #Try targetting about 3500 points per order, too many points per order will slow code down
                order_count = int(len(spectx)/order_split[1])
                order_split_wave = np.array_split(spectx, order_count)
                order_split_flux = np.array_split(specty, order_count)
                self.wavelength = order_split_wave
                self.flux = order_split_flux
                del order_split_wave
                del order_split_flux
                print('If values need to be replaced, use replace_w() function')
            else:
                self.wavelength = spectx
                self.flux = specty
        self.normalized_flux = copy.deepcopy(self.flux) #deepcopy creates new object
        self.combined_flux = None
        self.shifted_wavelength = copy.deepcopy(self.wavelength)
        self.estimated_shift = np.zeros(len(self.wavelength))
        self.rv = None #km/s

        #continuum information
        self.continuum = np.full(len(self.wavelength), None)
        #print('cont array empty', self.continuum)
        self.pred_all = np.full(len(self.wavelength), None)
        self.pred_var_all = np.full(len(self.wavelength), None)
        self.obs_err = np.full(len(self.wavelength), None)
        for i in range(len(self.wavelength)):
            false_array = np.full(len(self.wavelength[i]), False)
            #print('false array created', false_array)
            self.continuum[i] = false_array
            self.pred_all[i] = np.full(len(self.wavelength[i]), 0)
            self.pred_var_all[i] = np.full(len(self.wavelength[i]), 0)
            # abs(), not a bare sqrt: real raw-electron counts (e.g. GRACES's
            # "UnNormalized" flux column) can be genuinely negative --
            # ordinary CCD read noise scattering a near-zero-signal pixel
            # below zero after bias/dark subtraction, not a data error. A
            # bare sqrt(negative) is NaN, and even a modest fraction of NaN
            # weights was confirmed to corrupt an ENTIRE order's AsLS
            # continuum fit, not just the affected points (real GRACES
            # order with 539/4055, ~13%, negative-flux points came back
            # fully NaN). abs() keeps the correct ORDER OF MAGNITUDE of
            # Poisson-like noise at that pixel either way.
            self.obs_err[i] = np.sqrt(np.abs(self.flux[i]))
        #print('empty continuum arrays created', self.continuum)
        #Old way of setting these variables (change back if above code causes problems)
        # self.continuum = np.full((len(self.wavelength),len(self.wavelength[0])), False)
        # self.pred_all = np.zeros((len(self.wavelength),len(self.wavelength[0])))
        # self.pred_var_all = np.zeros((len(self.wavelength),len(self.wavelength[0])))
        # self.obs_err = np.zeros((len(self.wavelength),len(self.wavelength[0])))

        #line information
        self.lines = None
        #line - extra parameters [0] - shift continuum, [1] - left boundary in Angstroms
        #[2] - right boundary in Angstroms, [3] - line center in Angstroms
        self.lines_exp = None
        #line - extra data [0] - element, [1] - excitation potential
        #[2] - gf, [3] - rad
        self.lines_exd = None
        #line - equivalent width
        self.lines_ew = None
        #line - equivalent width error
        self.lines_ew_err = None
        #line - equivalent width calculated by simpon's rule integration
        self.lines_ew_simp = None
        #line - equivalent width error from integraion
        self.lines_ew_simp_err = None
        #line - best fit parameters for gaussian fit
        self.lines_bf_params = None
        #line - X squared value for gaussian and data
        self.lines_gauss_Xsquare = None
        #line - X squared threshold value
        self.X_thresh = 0.003
        #line - the actually-identified line center (Angstrom, in
        #shifted_wavelength), set by measure_ew(). Comparing this against
        #the rest wavelength (self.lines) is how check_for_flags() catches
        #likely misidentification -- see its docstring.
        self.lines_found_position = None
        #line - max allowed |lines_found_position - lines| (Angstrom) before
        #check_for_flags() flags a line as a possible misidentification.
        #Note get_line_window()'s own search is bounded to +/-0.1 A around
        #the expected position by default, so this threshold only has teeth
        #below that. Calibrated empirically (see conversation/session notes,
        #the 78-line real-star test): apply_rv_shift()-corrected positions
        #had RMS ~37 mA and MAX ~63 mA offset from the physically-expected
        #position even when correctly identified (real per-line
        #astrophysical scatter, e.g. convective blueshift differences --
        #not error). A first attempt at 0.05 (50 mA) sat below that
        #legitimate max and produced false positives; 0.07 leaves margin
        #above the good method's observed worst case while remaining well
        #below where the unreliable clean_shift()-extrapolated positions
        #commonly landed (many 90-160 mA off in the same test).
        self.position_thresh = 0.07
        #line - X squared value above threshold or EW = 0
        self.lines_check_flag = None
        #line - human-readable reason(s) lines_check_flag was set, '' if not
        #flagged. Set by check_for_flags(); make_ew_doc() uses this to
        #explain why a line was routed to the flagged-lines file.
        self.lines_flag_reasons = None
        #wavelength ranges where two orders' overlap disagreed badly
        #enough to distrust either one there -- set by
        #flag_order_overlaps(), consulted by check_for_flags(). Empty
        #(no-op) until flag_order_overlaps() is called.
        self.overlap_flag_ranges = []
        #wavelength ranges of whole orders containing an extreme,
        #non-astrophysical raw-flux outlier (cosmic ray/detector defect)
        #-- set by flag_bad_orders(), consulted by check_for_flags().
        #Empty (no-op) until flag_bad_orders() is called.
        self.bad_order_ranges = []
        #wavelengths of individual pixels CORRECTED for an extreme,
        #non-astrophysical spike (cosmic ray/sky-emission contamination/
        #bad pixel) -- set by correct_spikes() and, when a second
        #exposure is available, combine_spectra()'s cross-exposure check.
        #Consulted by check_for_flags(): a line whose measurement window
        #contains a corrected pixel is excluded, since its measurement
        #reflects a corrected (not originally observed) value. Empty
        #(no-op) until correct_spikes()/combine_spectra() is called.
        self.corrected_pixel_wavelengths = []
        #used to switch between Adamow ew calculation and simpson's rule integration
        self.temp_line_ew = None
        self.temp_line_ew_err = None

    def apply_response_correction(self, response_wave, response, min_overlap_fraction=0.5,
                                   min_response_fraction=0.1, response_bands=None, science_bands=None):
        """
        Divide out an instrument response/blaze correction curve. General:
        works with a response curve from any source, matched to this
        spectrum's own orders by wavelength overlap -- see
        response_correction.py's module docstring for why. Run this
        BEFORE normalize()/normalize_all(); it corrects the raw counts
        (self.flux), and normalize()'s continuum fit will be far more
        robust on an already-flattened spectrum.

        For MAROON-X specifically: response_wave/response/response_bands
        can come from readers.load_maroonx_response('MAROON-X_PHOENIX_
        RESPONSE_...hd5') -- a separate calibration file, not embedded in
        individual science exposures (confirmed empty there). Pass
        science_bands=readers.get_maroonx_bands(sci_filename) too --
        MAROON-X's two arms physically overlap in wavelength near their
        dichroic split, so wavelength overlap alone can silently match a
        science order to the WRONG arm's response chunk there (confirmed:
        produces a spurious ~4x monotonic trend across the whole order,
        not just an edge artifact). Without both band arguments, matching
        falls back to wavelength-overlap-only, same as before.

        Returns
        -------
        corrected_orders : list of order indices that were actually
            corrected (others may have been skipped -- see
            response_correction.apply_response_correction()'s docstring,
            including its min_response_fraction edge-pixel guard).
        """
        return _apply_response_correction(self, response_wave, response,
                                           min_overlap_fraction=min_overlap_fraction,
                                           min_response_fraction=min_response_fraction,
                                           response_bands=response_bands, science_bands=science_bands)

    def check_order_overlaps(self, min_overlap_points=10):
        """
        Compare normalized_flux between every pair of orders whose
        wavelength ranges overlap -- see overlap_check.py's module
        docstring for why this is a useful, ground-truth-free
        consistency check (it's how the MAROON-X arm-mismatch bug was
        originally found, just automated across every overlapping pair
        instead of relying on noticing two orders share a line). Run
        this AFTER normalize_all() (and apply_rv_shift(), if used).

        Returns
        -------
        list of dicts, sorted worst-first -- see
        overlap_check.check_order_overlaps()'s docstring for the fields.
        """
        return _check_order_overlaps(self, min_overlap_points=min_overlap_points)

    def flag_order_overlaps(self, threshold_pct=2.0, min_overlap_points=10):
        """
        Run check_order_overlaps() and remember which wavelength ranges
        disagree badly enough to distrust (see
        overlap_check.flagged_overlap_ranges()'s docstring for the
        default threshold's rationale). check_for_flags() consults
        self.overlap_flag_ranges automatically -- call this once before
        check_for_flags() (or make_ew_doc(), which calls it for you) if
        you want that check included; otherwise it's a silent no-op,
        same as never calling it.

        A line is flagged if its REST wavelength falls inside a flagged
        range, regardless of which of the two disagreeing orders it's
        actually measured from -- deliberately conservative: this check
        doesn't try to decide WHICH of the two orders is at fault (often
        genuinely ambiguous), only that the region is in dispute.

        Returns
        -------
        list of dicts -- the flagged subset; see
        overlap_check.flagged_overlap_ranges()'s docstring for the
        fields. Also stored on self.overlap_flag_ranges.
        """
        results = _check_order_overlaps(self, min_overlap_points=min_overlap_points)
        self.overlap_flag_ranges = _flagged_overlap_ranges(results, threshold_pct=threshold_pct)
        return self.overlap_flag_ranges

    def flag_bad_orders(self, outlier_factor=20.0):
        """
        Detect whole orders whose RAW flux contains an extreme, non-
        astrophysical outlier (cosmic ray or a fixed detector/amplifier-
        boundary defect) -- see outlier_check.py's module docstring for
        why this matters: a single such point was confirmed to corrupt an
        ENTIRE order's AsLS continuum fit, not just the affected point.
        Works on raw flux, so it can (and should) be called BEFORE
        normalize_all() -- catches the problem at its source rather than
        only after it has already produced a bad fit.

        check_for_flags() consults self.bad_order_ranges automatically --
        call this once before check_for_flags() (or make_ew_doc(), which
        calls it for you) if you want this check included; otherwise it's
        a silent no-op, same as never calling it. Same deliberately-
        conservative philosophy as flag_order_overlaps(): flags the whole
        order's wavelength range rather than trying to salvage it (e.g.
        by masking just the bad point and re-fitting), since a defect
        this extreme is a data-quality problem this package can't fix,
        only report.

        Returns
        -------
        list of dicts -- the flagged subset; see
        outlier_check.check_bad_orders()'s docstring for the fields. Also
        stored on self.bad_order_ranges.
        """
        self.bad_order_ranges = _check_bad_orders(self, outlier_factor=outlier_factor)
        return self.bad_order_ranges

    def correct_spikes(self, window=5, factor=5.0, skip_orders=None):
        """
        Detect and CORRECT individual extreme, non-astrophysical spikes
        (cosmic rays, hot pixels, uncorrected sky-emission-line
        contamination) in this spectrum's own raw flux -- works on a
        single spectrum with no second exposure needed, unlike
        combine_spectra()'s more sensitive two-exposure comparison. See
        outlier_check.detect_spikes()'s docstring for why this uses a
        factor-based (not sigma-based) local test, and why it only
        catches EXTREME, unambiguous cases by design -- a moderate,
        possibly-real excursion is deliberately left alone here.

        Unlike flag_bad_orders() (which flags a WHOLE order and changes
        no data), this actually REPLACES each affected point's flux with
        its local median (and recomputes obs_err there to match), and
        records the wavelength in self.corrected_pixel_wavelengths --
        consulted by check_for_flags() to EXCLUDE any line whose
        measurement window contains a corrected pixel entirely, since
        its measurement would reflect a corrected, not originally
        observed, value. Run this BEFORE normalize_all(), same as
        flag_bad_orders().

        Parameters
        ----------
        window, factor : passed to outlier_check.detect_spikes().
        skip_orders : iterable of order indices to skip (e.g. orders
            already flagged by flag_bad_orders() -- correcting a handful
            of points in a whole-order defect isn't meaningful; call
            flag_bad_orders() first and pass its flagged order indices
            here if you're using both).

        Returns
        -------
        list of wavelengths corrected (also appended to
        self.corrected_pixel_wavelengths, not overwritten -- safe to
        call this alongside combine_spectra()'s own corrections).
        """
        skip_orders = set(skip_orders or [])
        corrected = []
        for i in range(len(self.flux)):
            if i in skip_orders:
                continue
            flux = np.asarray(self.flux[i], dtype=float)
            bad, local_median = _detect_spikes(flux, window=window, factor=factor)
            if not bad.any():
                continue
            flux[bad] = local_median[bad]
            self.flux[i] = flux
            self.obs_err[i] = np.sqrt(np.abs(flux))
            corrected.extend(self.wavelength[i][bad].tolist())
        self.corrected_pixel_wavelengths.extend(corrected)
        return corrected

    def normalize_all(self, lam = 2e3, p = 0.01, n_iter = 15, adaptive = True, **als_kwargs):
        #loop through orders
        for i in range(len(self.flux)):

            #fit continuum via Asymmetric Least Squares smoothing
            self.normalize(i, lam=lam, p=p, n_iter=n_iter, adaptive=adaptive, **als_kwargs)

            #Replace un-normalized points with value before it
            #This should only be replacing the last point in the
            #spectrum that is always missed by normalize
            # err_est = self.obs_err[i]/self.pred_all[i]
            # non_norm_points = np.where(self.normalized_flux[i] > np.average(self.normalized_flux[i][self.continuum[i]]+err_est[self.continuum[i]]*100))
            # #replace non norm points
            # self.normalized_flux[i][non_norm_points] = self.normalized_flux[i][non_norm_points[0]-1]

        return None

    def normalize(self, order, clip = [-999,-999], lam = 2e3, p = 0.01, n_iter = 15,
                  adaptive = True, **als_kwargs):
        """
        Fit and divide out the continuum for one order via Asymmetric
        Least Squares (AsLS) smoothing (fit_als_continuum(), continuum.py).

        Every previous approach here (Continuum_scan's local-window
        selection, a GP, a low-order polynomial, a piecewise spline) fit
        to a HARD selection of "continuum" points, decided by some
        threshold that wasn't locally aware -- and every one of them
        eventually broke on real data where that threshold left a
        stretch of an order with too few trusted points (a broad line, a
        real order edge, a real large-scale continuum slope, or just an
        ordinary-looking stretch that happened to sit a bit below the
        order's brightest region). AsLS has no hard mask at all: every
        point gets a soft, iteratively-updated weight (favoring points
        above the current fit as likely continuum, without ever fully
        discarding points below it), and a smoothness penalty (`lam`)
        constrains the whole curve continuously rather than through
        discrete segments -- so no stretch of an order can end up
        completely unconstrained the way a starved spline segment could.

        A single global `lam` still had a real limit -- not stiff enough
        to resist a real, densely-blended stretch (found on a real Keck
        order), but a `lam` stiff enough to resist that overshot a real
        GRACES order's genuine large-scale continuum decline. `adaptive`
        (default True) fixes this by stiffening the fit locally wherever
        the data itself shows an extended stretch lacking a true nearby
        continuum peak (whether from dense blending or one broad line),
        while leaving `lam` at its flexible base value everywhere else
        (including across a real large-scale slope, which doesn't trip
        this criterion) -- see fit_als_continuum()'s docstring for the
        detection method and its parameters (stiffen_factor,
        local_window, wide_window, severity_threshold, passed through
        via **als_kwargs).

        See continuum.py's module docstring for the full history of what
        this replaced and why each earlier attempt failed on real data.

        lam/p/n_iter/adaptive control the fit; see
        fit_als_continuum()'s docstring.
        """
        if clip[0] != -999 and clip[1] != -999:
            #clipped = True
            clipl = np.where(self.wavelength[order] <= clip[0])[0][-1]
            clipr = np.where(self.wavelength[order] >= clip[1])[0][0]
        else:
            #clipped = False
            clipl = 0
            clipr = len(self.flux[order])

        wave = self.wavelength[order][clipl:clipr]
        flux = self.flux[order][clipl:clipr]
        # self.obs_err carries the correct per-point error: sqrt(flux) as
        # set at __init__ time for an uncorrected order, or (if
        # apply_response_correction() has run) sqrt(raw_counts)/response --
        # NOT recomputed as sqrt(flux) here, which would silently drop the
        # response-division noise amplification (see
        # response_correction.py's module docstring for the real bug this
        # caused: an artificial continuum-fit gradient across an otherwise
        # clean, response-corrected MAROON-X order).
        err = self.obs_err[order][clipl:clipr]

        pred, pred_var = fit_als_continuum(wave, flux, err, lam=lam, p=p, n_iter=n_iter,
                                           adaptive=adaptive, **als_kwargs)

        self.continuum[order][clipl:clipr] = flux >= pred
        self.pred_all[order][clipl:clipr] = pred
        self.pred_var_all[order][clipl:clipr] = pred_var
        self.normalized_flux[order][clipl:clipr] = flux/pred
        return None

    def S_N(self, rows = 5, cols = 4, save_plot = False):
        #Siganl to noise is overestimated compared to MAKEE output
        plt.clf()
        f = plt.figure(figsize=(20,15))
        plt.suptitle(str(self.filename) + ' S/N', size = 20)
        for i in range(len(self.wavelength)):
            order = i
            i += 1
            ax = f.add_subplot(rows, cols,i)
            ax.plot(self.wavelength[order],np.sqrt(self.flux[order]*self.gain[order]))
            ax.set_title(str(i))
            ax.grid()
        plt.tight_layout()
        if save_plot:
            plt.savefig(str(self.filename)+'_SNR.pdf')
        plt.show()

    def load_normalized(self, name):
        pathnames = glob.glob(name+'*'+'.npy')
        for i in range(len(pathnames)):
            if '_flux' in pathnames[i]:
                self.normalized_flux = np.load(pathnames[i])
                print('flux loaded')
            elif '_wavelength' in pathnames[i]:
                self.wavelength = np.load(pathnames[i])
                print('wavelength loaded')
            elif '_cont' in pathnames[i]:
                self.continuum = np.load(pathnames[i])
                print('continuum loaded')
            elif '_obs_err' in pathnames[i]:
                self.obs_err = np.load(pathnames[i])
                print('errors loaded')
            elif '_pred' in pathnames[i]:
                self.pred_all = np.load(pathnames[i])
                print('all preds loaded')
            elif '_pred_var' in pathnames[i]:
                self.pred_var_all = np.load(pathnames[i])
                print('all pred vars loaded')

    def save_normalized(self, name):
        #arrays are separated by order
        np.save(name+'_flux',self.normalized_flux)
        np.save(name+'_wavelength',self.shifted_wavelength)
        np.save(name+'_cont',self.continuum)
        np.save(name+'_obs_err',self.obs_err)
        np.save(name+'_pred',self.pred_all)
        np.save(name+'_pred_var',self.pred_var_all)
        return None

    def wave_shift(self, order, shift):
        self.shifted_wavelength[order] = self.wavelength[order] + shift
        self.estimated_shift[order] = shift

    def combine_spectra(self, spectB, rv_A=None, rv_B=None, rv_shift=True,
                        min_overlap_fraction=0.5, align=True, align_search_radius=0.1,
                        align_warn_threshold=0.05, cross_exposure_check=True,
                        cross_exposure_significance=8.0, plot=False, plot_dir='.', verbose=False):
        """
        Co-add this spectrum with a second exposure of the SAME star
        (spectB), to reach a single, higher-S/N combined spectrum.
        normalize_all() (and, if desired, flag_bad_orders()) must already
        have been run on BOTH spectra -- apply_rv_shift() below needs
        normalized_flux to measure each spectrum's own RV.

        Redesigned from an earlier version (see DEVELOPMENT_LOG.md) that
        cross-correlated A directly against B via estimate_shift() (built
        for star-vs-solar-atlas matching, not two exposures of the SAME
        star) and combined flux via a fixed +/-5-INDEX nearest-neighbor
        window -- silently assuming A and B shared an identical pixel
        grid/length per order, fragile for two independently-reduced
        exposures. This version:

        1. Independently rest-frames BOTH spectra via apply_rv_shift()
           (each spectrum's own robust RV measurement) instead of cross-
           correlating them against each other -- reuses this package's
           already-validated, more robust RV machinery (named-line +
           linelist cross-check, no reference-coverage-gap failure mode)
           rather than estimate_shift()/clean_shift(). Set rv_shift=False
           to skip this (e.g. both spectra are already correctly rest-
           framed); pass rv_A/rv_B to apply a known RV directly instead
           of remeasuring it.
        2. Matches each of A's orders to whichever of B's orders has the
           MOST real wavelength-range overlap (requiring at least
           min_overlap_fraction of A's own range), not nearest-median --
           an order with no genuine match is left as A alone rather than
           silently paired with something unrelated.
        3. Properly interpolates B's flux/error onto A's rest-frame
           wavelength grid (scipy.interpolate.interp1d) over the real
           overlapping range only.
        4. Propagates error as a quadrature sum of each spectrum's own
           obs_err (interpolated), rather than re-deriving purely from
           sqrt(combined counts).
        5. Skips combining any order flagged by flag_bad_orders() in
           EITHER spectrum (self.bad_order_ranges/spectB.bad_order_ranges,
           silent no-op if never called, same convention as every other
           flag_*() method) -- combining a good order with a known-
           corrupted one would spread the corruption, not average it out;
           left as that spectrum's own order alone instead.
        6. Plotting is opt-in (plot=True) and non-blocking (saved as one
           PNG per combined order under plot_dir, not a blocking
           plt.show() per order).
        7. Per-order FINE alignment (align=True, default): independent
           RV correction alone leaves a real, measurable residual
           mismatch between A and B (confirmed on real data: 5-29 mA
           across different orders, from ordinary per-line RV
           measurement noise in each spectrum's own independent RV) --
           enough to measurably smear a naively-co-added spectrum. Before
           combining each order pair, directly cross-correlates A and B's
           own overlapping normalized flux (combine.measure_order_
           alignment()) and applies the measured residual to B's
           wavelength query before interpolating -- a much more precise,
           data-driven correction than trusting either spectrum's
           independent RV alone. Warns if the measured residual exceeds
           align_warn_threshold (default 0.05 A, well above every normal
           residual confirmed on real data).
        8. Cross-exposure spike detection (cross_exposure_check=True,
           default): at each aligned point, checks whether EITHER
           spectrum shows a significant EXCESS above continuum
           (outlier_check.check_cross_exposure_spikes()) -- catches real
           contamination (cosmic rays, uncorrected sky-emission lines --
           confirmed on real data: the classic 5577 A and 6300 A [OI]
           auroral lines) too modest for detect_spikes()'s single-
           spectrum-only test to safely catch alone, without confusing
           it for ordinary real line-core noise (tested and rejected: a
           plain two-spectrum DIFFERENCE test alone misidentifies dozens
           of real, ordinary line-core measurement differences as
           contamination -- see check_cross_exposure_spikes()'s
           docstring for why comparing each side's own excess above
           continuum, not the difference between them, is what actually
           works). A flagged point uses ONLY the LESS-contaminated
           side's raw flux/error (whichever has the lower excess
           significance, not necessarily zero -- a real sky line
           present in both exposures at different strength is still
           flagged at its most-contaminated point, not just its weaker
           edges) rather than the sum, and is recorded in
           self.corrected_pixel_wavelengths -- consulted by
           check_for_flags() to exclude any line measured there.

        Parameters
        ----------
        spectB : Spectrum_Data -- the second exposure of the same star.
        rv_A, rv_B : float, km/s, optional -- see point 1 above.
        rv_shift : bool -- see point 1 above.
        min_overlap_fraction : minimum fraction of A's order wavelength
            range a candidate B order must cover to be paired with it.
        align, align_search_radius, align_warn_threshold : see point 7.
        cross_exposure_check, cross_exposure_significance : see point 8;
            the latter passed to check_cross_exposure_spikes() as
            `significance`.
        plot, plot_dir : see point 6 above.
        verbose : passed through to apply_rv_shift(); also prints each
            order's match/skip decision, measured alignment residual, and
            any cross-exposure corrections here.

        Returns
        -------
        n_combined : int -- number of orders actually combined (out of
            len(self.flux); the rest are that order's A data unchanged).
            Call update_combined() once satisfied, to make this the
            spectrum's own self.flux/self.obs_err (self.combined_flux/
            self.combined_err hold the result until then).
        """
        if rv_shift:
            self.apply_rv_shift(rv=rv_A, verbose=verbose)
            spectB.apply_rv_shift(rv=rv_B, verbose=verbose)

        bad_A = {r['order'] for r in (self.bad_order_ranges or [])}
        bad_B = {r['order'] for r in (spectB.bad_order_ranges or [])}

        combined_flux = np.empty(len(self.flux), dtype=object)
        combined_err = np.empty(len(self.flux), dtype=object)
        n_combined = 0
        n_cross_corrected = 0

        if plot:
            os.makedirs(plot_dir, exist_ok=True)

        for i in range(len(self.shifted_wavelength)):
            wave_A, flux_A, err_A = self.shifted_wavelength[i], self.flux[i], self.obs_err[i]

            if i in bad_A:
                combined_flux[i], combined_err[i] = flux_A, err_A
                if verbose:
                    print(f"order {i}: flagged bad in A -- using A alone, not combined")
                continue

            a_lo, a_hi = wave_A.min(), wave_A.max()
            a_span = a_hi - a_lo
            best_j, best_overlap = None, 0.0
            for j in range(len(spectB.shifted_wavelength)):
                if j in bad_B:
                    continue
                wave_B = spectB.shifted_wavelength[j]
                overlap = max(0.0, min(a_hi, wave_B.max()) - max(a_lo, wave_B.min()))
                frac = overlap / a_span if a_span > 0 else 0.0
                if frac > best_overlap:
                    best_overlap, best_j = frac, j

            if best_j is None or best_overlap < min_overlap_fraction:
                combined_flux[i], combined_err[i] = flux_A, err_A
                if verbose:
                    print(f"order {i}: no usable B order overlap "
                          f"(best {best_overlap:.2f}) -- using A alone")
                continue

            wave_B, flux_B, err_B = (spectB.shifted_wavelength[best_j],
                                      spectB.flux[best_j], spectB.obs_err[best_j])

            residual = 0.0
            if align:
                measured = measure_order_alignment(
                    wave_A, self.normalized_flux[i], wave_B, spectB.normalized_flux[best_j],
                    search_radius=align_search_radius)
                if measured is not None:
                    residual = measured
                    if verbose:
                        print(f"order {i}: measured alignment residual = {residual*1000:.2f} mA")
                    if abs(residual) > align_warn_threshold:
                        print(f"WARNING: order {i}/{best_j} alignment residual "
                              f"({residual*1000:.1f} mA) exceeds align_warn_threshold "
                              f"({align_warn_threshold*1000:.1f} mA) -- larger than any "
                              f"residual confirmed normal on real data; inspect this order "
                              f"pair before trusting the combined result.")

            in_range = (wave_A >= wave_B.min() + abs(residual)) & (wave_A <= wave_B.max() - abs(residual))
            query = wave_A[in_range] + residual
            flux_B_interp = interp1d(wave_B, flux_B, kind='linear')(query)
            err_B_interp = interp1d(wave_B, err_B, kind='linear')(query)

            cflux, cerr = flux_A.copy(), err_A.copy()
            cflux[in_range] = flux_A[in_range] + flux_B_interp
            cerr[in_range] = np.sqrt(err_A[in_range]**2 + err_B_interp**2)

            if cross_exposure_check:
                cx_in_range, cx_bad, cx_a_higher = _check_cross_exposure_spikes(
                    wave_A, self.normalized_flux[i], self.pred_all[i], self.obs_err[i],
                    wave_B, spectB.normalized_flux[best_j], spectB.pred_all[best_j],
                    spectB.obs_err[best_j], residual=residual, significance=cross_exposure_significance)
                # cx_in_range uses the identical formula (same residual) as
                # this method's own in_range above, so they're the same
                # mask -- cx_bad/cx_a_higher already align positionally
                # with flux_B_interp/err_B_interp (both built from wave_A[in_range]).
                if cx_bad.any():
                    combo_idx = np.where(in_range)[0]
                    for pos in np.where(cx_bad)[0]:
                        idx = combo_idx[pos]
                        # use whichever side has the LOWER excess significance
                        # (the less-contaminated one), not the sum of both
                        if cx_a_higher[pos]:
                            cflux[idx], cerr[idx] = flux_B_interp[pos], err_B_interp[pos]
                        else:
                            cflux[idx], cerr[idx] = flux_A[idx], err_A[idx]
                        self.corrected_pixel_wavelengths.append(float(wave_A[idx]))
                        n_cross_corrected += 1
                    if verbose:
                        print(f"order {i}: {cx_bad.sum()} point(s) cross-exposure-corrected "
                              f"(significant excess above continuum in at least one spectrum)")

            combined_flux[i], combined_err[i] = cflux, cerr
            n_combined += 1
            if verbose:
                print(f"order {i}: combined with B order {best_j} (overlap={best_overlap:.2f})")

            if plot:
                fig, ax = plt.subplots(figsize=(10, 4))
                ax.plot(wave_A, flux_A, label='A', lw=0.7)
                ax.plot(wave_A[in_range], flux_B_interp, label='B (interp)', lw=0.7)
                ax.plot(wave_A, cflux, label='A+B', lw=0.7, color='k')
                ax.set_title(f'order {i} (A) <-> order {best_j} (B), overlap={best_overlap:.2f}')
                ax.legend()
                ax.grid()
                plt.tight_layout()
                fig.savefig(os.path.join(plot_dir, f'combine_order_{i:02d}.png'), dpi=100)
                plt.close(fig)

        self.combined_flux = combined_flux
        self.combined_err = combined_err
        print(f"combine_spectra(): combined {n_combined}/{len(self.flux)} orders "
              f"({len(bad_A)} flagged bad in A, {len(bad_B)} flagged bad in B, "
              f"{n_cross_corrected} point(s) cross-exposure-corrected)")
        print('Use self.update_combined() when you are happy with the combined flux/error '
              'to override self.flux/self.obs_err.')
        return n_combined

    def update_combined(self):
        self.flux = self.combined_flux
        self.obs_err = self.combined_err

    def estimate_shift(self, sun_spectra, shift_max = 5, shift_min = -5, shift_spacing = 100, verbose = False):
        """
        Per-order wavelength shift via grid-search cross-correlation against
        a reference spectrum (e.g. a solar atlas). Only reliable for orders
        that actually overlap the reference spectrum's own coverage -- see
        clean_shift()'s docstring for what happens otherwise.

        For preparing a spectrum for EW measurement specifically, prefer
        apply_rv_shift() instead: it needs no reference spectrum at all (so
        it has no coverage-gap failure mode), and was measured to track the
        true line position noticeably more precisely (RMS ~37 mA vs ~67 mA
        on a real test) where the two methods could be directly compared.
        estimate_shift() remains the right tool for reference-spectrum
        cross-correlation against a genuinely different spectrum (e.g. a
        solar atlas). combine_spectra() (co-adding two exposures of the
        SAME star) now uses apply_rv_shift() on each spectrum
        independently instead -- see its own docstring for why.
        """
        #setup num orders, place holder for chi min, shifts array
        orders = len(self.wavelength)
        chi = np.zeros(shift_spacing)
        shifts = np.linspace(shift_min, shift_max, shift_spacing)

        #begin looping through each order in star spectrum
        for q in range(orders):
            order = q
            order_found = False
            order_mean = self.shifted_wavelength[order].mean()

            #loop through each included solar spectrum to find a matching order
            for j in range(len(sun_spectra)):
                sun = sun_spectra[j]

                #Loop through sun orders to find similar order
                for i in range(len(sun.wavelength)):
                    #use average wavelength value in an order to look for
                    #matching solar order
                    if order_mean < sun.wavelength[i].max() and order_mean > sun.wavelength[i].min():
                        if verbose:
                            print('Order', q, 'matched to order in solar spectrum')
                        order_found = True
                        sun_range = sun.wavelength[i].max() - sun.wavelength[i].min()
                        sun_window = [sun.wavelength[i][0], sun.wavelength[i][-1]]

                        #loop through specified shift values to find best match
                        for k in range(len(shifts)):
                            shift = shifts[k]
                            self.wave_shift(order,shift)
                            order_range = self.shifted_wavelength[order].max() - self.shifted_wavelength[order].min()
                            order_window = [self.shifted_wavelength[order][0], self.shifted_wavelength[order][-1]]

                            #    |------------|     -->    |------------|
                            #|------------|         --> |xx|----------|
                            if order_window[0] > sun_window[0]:
                                sun_window[0] = sun.wavelength[i][np.where(sun.wavelength[i] > order_window[0])][0]

                            #|------------|         --> |------------|
                            #    |------------|     -->     |-------|xxxxx|
                            if order_window[1] < sun_window[1]:
                                sun_window[1] = sun.wavelength[i][np.where(sun.wavelength[i] < order_window[1])][-1]

                            compare_window_star = np.where((self.shifted_wavelength[order] >= order_window[0])&(self.shifted_wavelength[order] <= order_window[1]))
                            compare_window_sun = np.where((sun.wavelength[i] >= sun_window[0])&(sun.wavelength[i] <= sun_window[1]))
                            spect_interp = interp1d(self.shifted_wavelength[order][compare_window_star],self.normalized_flux[order][compare_window_star], kind = 'linear')
                            result_y = spect_interp(sun.wavelength[i][compare_window_sun])

                            diff = (result_y - sun.normalized_flux[i][compare_window_sun])**2
                            # mean, not sum: the overlap window (and hence the
                            # number of compared points) shrinks near the edges
                            # of the shift search range, which would otherwise
                            # bias the minimum toward whichever trial shift
                            # happens to produce the smallest overlap
                            chi[k] = np.mean(diff)

                            #chi[k] = chisquare(result_y, constraint_value)[0]
                            #chi[k] = chisquare(result_y, sun.normalized_flux[i][compare_window_sun])[0]

                        if verbose:
                            k_min = np.argmin(chi)
                            min_shift = parabolic_refine(shifts, chi, k_min)
                            print('The best shift is', min_shift)

                        break
                if order_found:
                    break
            if not order_found:
                #if order is not found in solar spectrum, value is set to -999
                #to be cleaned later
                print('order '+str(order)+' in star not found in solar specturm')
                print('missing values can be interpolated/extrapolated using clean_shift() method')
                self.estimated_shift[order] = -999.0
            else:
                try:
                    k_min = np.argmin(chi)
                    min_shift = parabolic_refine(shifts, chi, k_min)
                    self.estimated_shift[order] = min_shift
                    self.wave_shift(order,min_shift)
                except:
                    print('order ' + str(order) + ' experienced a problem')
                    print('missing values can be interpolated/extrapolated using clean_shift() method')

    def clean_shift(self):
        """
        NOTE (found by testing, see conversation/session notes): the
        interpolation/extrapolation this does for orders with no reference-
        spectrum match is only reliable *within* the wavelength range where
        real matches were actually found -- extrapolating beyond that range
        was measured to be substantially wrong (RMS ~67 mA, worst case
        ~119 mA off the true line position on a real cross-instrument test),
        enough to risk misidentifying a line during EW measurement. This
        method now warns when that happens. For EW-measurement prep, prefer
        apply_rv_shift() instead -- it doesn't depend on reference-spectrum
        coverage at all, so it doesn't have this failure mode.
        """
        #remove orders not found
        gd = np.where(self.estimated_shift != -999)
        bad_points = np.where(self.estimated_shift == -999)[0]

        #get mean wavelength values
        means = np.zeros(len(self.shifted_wavelength))
        for i in range(len(self.estimated_shift)):
            means[i] = self.shifted_wavelength[i].mean()

        #wavelength range actually covered by real reference-spectrum matches --
        #used below to tell interpolation (safer) apart from extrapolation
        #(risky, see docstring)
        gd_wave_min = means[gd].min()
        gd_wave_max = means[gd].max()

        #get standard deviation of good points
        stds = self.estimated_shift[gd].std()

        #fit to good points
        best_fit = np.polyfit(means[gd], self.estimated_shift[gd], 1)
        x_range = np.linspace(means[gd][0],means[gd][-1],len(self.shifted_wavelength))
        line = make_line(x_range, best_fit[0], best_fit[1])

        #residuals of good points
        residuals = self.estimated_shift - line

        #remove points > std/2
        good_points = np.where((residuals < stds/2)&(residuals > -10))
        bad_points = (np.append(bad_points, np.where(residuals > stds/2)[0]),)

        #fit again with best points
        best_fit = np.polyfit(means[good_points], residuals[good_points], 1)
        x_range = np.linspace(means[good_points][0],means[good_points][-1],len(residuals))

        #line to be used to interpolate and extrapolate missing or bad points
        line2 = make_line(x_range, best_fit[0], best_fit[1])

        #interpolate or extrapolate bad points using line 2 and replace values
        for i in range(len(self.estimated_shift[bad_points])):
            current_index = bad_points[0][i]
            wave = means[current_index]
            if wave < gd_wave_min or wave > gd_wave_max:
                print(f"WARNING: order {current_index} (mean wavelength {wave:.1f} A) has "
                      f"no reference-spectrum match and falls OUTSIDE the "
                      f"{gd_wave_min:.1f}-{gd_wave_max:.1f} A range where real matches were "
                      f"found -- its shift is an EXTRAPOLATION, which was measured to be "
                      f"unreliable (see clean_shift()'s docstring) and risks misidentifying "
                      f"lines in this order during EW measurement. Prefer apply_rv_shift() "
                      f"for this wavelength range instead.")
            self.estimated_shift[current_index] = make_line(wave, best_fit[0], best_fit[1]) + line[current_index]
            self.wave_shift(current_index, self.estimated_shift[current_index])

        #get radial velocity from slope of shifts
        best_fit, C = np.polyfit(means, self.estimated_shift*(-1), 1, cov=True)
        self.rv = (np.round(best_fit[0]*3e5,3), np.round(np.sqrt(np.diag(C))[1], 3))

    def apply_rv_shift(self, rv=None, lines=None, min_depth=0.02, verbose=False,
                       cross_check=True, linelist_search_radius=1.0,
                       linelist_min_significance=3.0, disagreement_kms=2.0):
        """
        RECOMMENDED default for preparing a spectrum for EW measurement.
        Shift every order by a single effective radial velocity, applied
        correctly as a multiplicative (1 + rv/c) wavelength scaling rather
        than a per-order independent Angstrom offset (see estimate_shift()).

        Unlike estimate_shift(), this needs no reference spectrum and no
        per-order overlap matching -- only a handful of radial_velocity.
        RV_REFERENCE_LINES need to fall somewhere in this spectrum's
        coverage -- so it has none of estimate_shift()/clean_shift()'s
        reference-coverage-gap failure mode (measured on a real test: RMS
        ~37 mA vs ~67 mA position error where the two methods could be
        directly compared, with clean_shift()'s extrapolated orders off by
        up to 119 mA -- enough to risk misidentifying a line). normalize()/
        normalize_all() must be run first (uses normalized_flux to locate
        line centers).

        RV_REFERENCE_LINES has two real weaknesses on its own: it can be
        entirely absent from a spectrum whose coverage happens to miss all
        of Ca II H&K/Balmer/Mg b/Na D, and even when present, mixing
        Balmer lines with metal lines in one average can be actively
        wrong, not just imprecise -- confirmed on real data, the two
        families disagreed by ~10 km/s (a real difference in line
        formation physics between H and metal lines, which naive sigma-
        clipping over only 3-4 lines has no way to separate from genuine
        measurement noise). If a science linelist is already loaded
        (self.lines, via load_lines()), this now also measures an
        independent RV from it (radial_velocity.measure_rv_from_linelist(),
        which reuses identify_lines()'s own detection-based centering --
        see its docstring) and uses it as follows: as the ONLY estimate if
        RV_REFERENCE_LINES found nothing usable at all; as a preferred
        replacement if the two estimates disagree by more than
        `disagreement_kms` (averaging over dozens of real linelist lines
        is more robust than 3-4 mixed-species reference lines); otherwise
        the (cheaper, already-computed) named-line RV is kept and the
        linelist estimate serves only as a passive cross-check. Call
        load_lines() before this if you want that cross-check available;
        it's a silent no-op (identical to the old behavior) if no linelist
        is loaded yet.

        Parameters
        ----------
        rv : float, km/s, optional -- apply this RV directly and skip line
            measurement (e.g. if the RV is already known from elsewhere).
        lines : {name: (rest_wavelength, window)}, optional -- defaults to
            radial_velocity.RV_REFERENCE_LINES.
        min_depth : float -- minimum line depth (in normalized flux) to
            trust a line's fitted center.
        verbose : bool -- print the measured RV(s), which lines were used,
            and the cross-check outcome.
        cross_check : bool -- if False, use RV_REFERENCE_LINES only, same
            as before this parameter existed.
        linelist_search_radius, linelist_min_significance : passed to
            measure_rv_from_linelist() as search_radius/min_significance.
        disagreement_kms : how far the two estimates must differ before
            the linelist-based one is preferred over the named-line one.

        Returns
        -------
        rv : float, km/s -- the RV actually applied.
        """
        if rv is None:
            measured_rv, rv_err, used = measure_effective_rv(self, lines=lines, min_depth=min_depth)

            linelist_rv, linelist_rv_err, linelist_n = None, None, 0
            if cross_check and self.lines is not None and len(self.lines) > 0:
                linelist_rv, linelist_rv_err, linelist_n = measure_rv_from_linelist(
                    self.lines, self.wavelength, self.normalized_flux, self.obs_err, self.pred_all,
                    search_radius=linelist_search_radius, min_significance=linelist_min_significance)

            if measured_rv is None and linelist_rv is None:
                raise ValueError(
                    "Could not measure an effective RV -- none of the reference "
                    "lines were found/usable in this spectrum's wavelength "
                    "coverage, and no usable linelist-based fallback was "
                    "available either (load_lines() first to enable that, or "
                    "pass rv= directly).")
            elif measured_rv is None:
                if verbose:
                    print(f"No usable RV_REFERENCE_LINES; falling back to linelist-based RV = "
                          f"{linelist_rv:.3f} +/- {linelist_rv_err:.3f} km/s from {linelist_n} line(s).")
                measured_rv, rv_err = linelist_rv, linelist_rv_err
            else:
                if verbose:
                    print(f"Effective RV = {measured_rv:.3f} +/- {rv_err:.3f} km/s, "
                          f"from {len(used)} line(s):")
                    for name, restw, order, v in used:
                        print(f"  {name} ({restw} A, order {order}): v={v:.3f} km/s")
                if linelist_rv is not None:
                    disagreement = abs(measured_rv - linelist_rv)
                    if disagreement > disagreement_kms:
                        if verbose:
                            print(f"WARNING: named-line RV ({measured_rv:.3f} km/s) and linelist RV "
                                  f"({linelist_rv:.3f} +/- {linelist_rv_err:.3f} km/s, {linelist_n} lines) "
                                  f"disagree by {disagreement:.3f} km/s (> {disagreement_kms}) -- "
                                  f"preferring the linelist RV as the more robust (larger-N) estimate.")
                        measured_rv, rv_err = linelist_rv, linelist_rv_err
                    elif verbose:
                        print(f"Linelist cross-check OK: {linelist_rv:.3f} +/- {linelist_rv_err:.3f} km/s "
                              f"from {linelist_n} line(s), within {disagreement_kms} km/s of the named-line RV.")

            self.rv = (round(measured_rv, 3), round(rv_err, 3))
            rv = measured_rv
        else:
            self.rv = (rv, 0.0)

        for order in range(len(self.wavelength)):
            # v = c*(observed-rest)/rest (measure_line_velocity()'s
            # convention: positive v = redshifted/receding), so
            # observed = rest*(1+v/c) is the FORWARD relation -- to
            # recover rest-frame wavelength from the observed spectrum
            # (the actual point of this method) needs the INVERSE,
            # rest = observed/(1+v/c) =~ observed*(1-v/c) for v << c.
            # This sign was wrong from when this method was introduced
            # (commit 8c3b32c): it used (1+v/c), which does not correct
            # the shift but DOUBLES it. That escaped detection because
            # the original validation compared two spectra DIFFERENTIALLY
            # (star RV minus reference-star RV), which partially cancels
            # a sign error applied consistently to both; confirmed wrong
            # directly on real data once compared against absolute
            # rest-wavelength positions from a real linelist (see
            # DEVELOPMENT_LOG.md): applying (1+v/c) took a -73.1 mA
            # pre-shift residual to -146.1 mA (doubled); (1-v/c) took it
            # to 0.0 mA.
            self.shifted_wavelength[order] = self.wavelength[order] * (1.0 - rv / C_KMS)
            # Angstrom-equivalent at the order's mean wavelength, kept for
            # reporting/consistency with estimate_shift()'s convention
            # (wave_shift(): shifted_wavelength = wavelength + shift) --
            # the actually-applied shift above is the correct multiplicative
            # one, not this per-order scalar approximation of it.
            self.estimated_shift[order] = self.wavelength[order].mean() * (-rv / C_KMS)

        return rv

    def load_lines(self, filename):
        self.lines = np.genfromtxt(filename, skip_header = 1, usecols = 0)
        elmnt = np.genfromtxt(filename, skip_header = 1, usecols = 1)
        ep = np.genfromtxt(filename, skip_header = 1, usecols = 2)
        gf = np.genfromtxt(filename, skip_header = 1, usecols = 3)
        rad = np.genfromtxt(filename, skip_header = 1, usecols = 4)
        self.lines_exd = np.zeros((len(self.lines),4))
        self.lines_exp = np.zeros((len(self.lines),4))
        self.lines_ew = np.zeros(len(self.lines))
        self.lines_ew_err = np.zeros(len(self.lines))
        self.lines_ew_simp = np.zeros(len(self.lines))
        self.lines_ew_simp_err = np.zeros(len(self.lines))
        self.lines_bf_params = np.array([None]*len(self.lines))
        self.lines_gauss_Xsquare = np.array([np.nan]*len(self.lines))
        self.lines_found_position = np.array([np.nan]*len(self.lines))
        self.lines_check_flag = np.array([False]*len(self.lines))
        self.lines_flag_reasons = np.array(['']*len(self.lines), dtype=object)
        #identify_lines() results -- a separate, prior step from EW
        #measurement (see line_identification.py's module docstring);
        #all default to "not run yet", distinct from lines_found_position
        #(set by measure_ew() during actual measurement)
        self.lines_id_detected = np.array([False]*len(self.lines))
        self.lines_id_position = np.array([np.nan]*len(self.lines))
        self.lines_id_order = np.array([None]*len(self.lines))
        self.lines_id_significance = np.array([0.0]*len(self.lines))
        self.lines_id_blended = np.array([False]*len(self.lines))
        #distinguishes "identify_lines() never called" from "called and
        #found nothing" -- lines_id_detected defaults to False either
        #way, so check_for_flags() needs this to avoid flagging every
        #single line as undetected when identify_lines() simply hasn't
        #run yet
        self.lines_id_run = False
        for i in range(len(self.lines)):
            self.lines_exd[i] = np.array([elmnt[i],ep[i],gf[i],rad[i]])

    def identify_lines(self, **kwargs):
        """
        Locate every loaded line (self.lines) in the spectrum, as a
        distinct step BEFORE any EW measurement is attempted -- see
        line_identification.py's module docstring for the detection-
        based approach and why it replaces "whatever the nearest local
        minimum happens to be" with an explicit real-detection test.

        Run this AFTER normalize_all() (and apply_rv_shift(), if used)
        and load_lines(). Populates, per line:
            lines_id_detected : bool -- a real absorption feature was
                found somewhere in the search window; False means no
                candidate cleared min_significance anywhere searched
                (too weak for this spectrum's S/N, or genuinely absent)
                -- distinct from a low-confidence detection, and NOT
                silently treated as "found at the rest wavelength" the
                way the older get_line_window()-based path would.
            lines_id_position : identified center (Angstrom), NaN if
                not detected.
            lines_id_order : which order the (best) detection came
                from, None if not detected in any candidate order.
            lines_id_significance : the detection's depth in units of
                local noise sigma (0 if not detected).
            lines_id_blended : a second, competitive candidate exists
                nearby -- the measured region isn't a clean, isolated
                feature even though a center was identified.

        **kwargs passed through to identify_line() (search_radius,
        min_significance, position_tolerance, etc.)

        Returns
        -------
        list of per-line result dicts (see
        line_identification.identify_lines_in_spectrum()'s docstring).
        """
        results = _identify_lines_in_spectrum(
            self.lines, self.shifted_wavelength, self.normalized_flux,
            self.obs_err, self.pred_all, **kwargs)
        for i, r in enumerate(results):
            self.lines_id_detected[i] = r['detected']
            self.lines_id_position[i] = r['center']
            self.lines_id_order[i] = r['order']
            self.lines_id_significance[i] = r['significance']
            self.lines_id_blended[i] = r['blended']
        self.lines_id_run = True
        return results

    # def switch_ew_values(self):
    #     if self.temp_line_ew == None:
    #         print("switching from Adamow calculation for EW to Simpson's rule integration")
    #         self.temp_line_ew = self.lines_ew
    #         self.temp_line_ew_err = self.lines_ew_err

    #         self.lines_ew = self.lines_ew_simp
    #         self.lines_ew_err = self.lines_ew_simp_err
    #     else:
    #         print("switching from Simpson's rule integration for EW to Adamow calculation for EW")
    #         self.lines_ew = self.temp_line_ew
    #         self.lines_ew_err = self.temp_line_ew_err

    #         self.temp_line_ew = self.lines_ew_simp
    #         self.temp_line_ew_err = self.lines_ew_simp_err

    def make_ew_doc(self, name, doc_title='STARNAME, PROJECT, YEAR; ', flagged_name=None):
        """
        Write the measured EWs as a MOOG-format linelist.

        For running unattended (little to no user interaction): lines that
        check_for_flags() marks as untrustworthy are NOT written to `name`
        -- a suspicious measurement should never silently end up in the
        science linelist. They're written instead to a separate file
        (flagged_name, default: `name` with "_flagged" inserted before the
        extension), each annotated with its flag reason(s), for later
        interactive/visual review -- not for automatic consumption.

        Calls check_for_flags() itself, so it reflects the current
        measurements even if you haven't called it explicitly.

        Returns
        -------
        removed_lines : ndarray -- rest wavelengths of lines never measured
            at all (EW == 0, e.g. excluded via measure_all_ew()'s
            exclude_lines). Distinct from flagged lines, which WERE
            measured but look untrustworthy.
        """
        if flagged_name is None:
            if '.' in name:
                base, ext = name.rsplit('.', 1)
                flagged_name = f'{base}_flagged.{ext}'
            else:
                flagged_name = f'{name}_flagged'

        self.check_for_flags()

        doc = open(name, 'w')
        doc.write(doc_title+'Extended Fe Linelist based on the SWP (2010) paper plus additions from Ivan\n')
        flagged_doc = open(flagged_name, 'w')
        flagged_doc.write(doc_title + 'FLAGGED lines -- excluded from the main linelist above; for '
                           'interactive/visual review, not automatic use. Trailing comment is the '
                           'flag reason(s) from check_for_flags().\n')
        removed_lines = []
        for i in range(len(self.lines)):
            if self.lines_ew[i] != 0.0:
                wave = "  "+str(self.lines[i])
                elmnt = str(self.lines_exd[i][0])
                if self.lines_exd[i][0] < 10:
                    elmnt = '0'+elmnt
                while len(elmnt) < 6:
                    elmnt = elmnt+'0'
                ep = str(self.lines_exd[i][1])
                gf = str(self.lines_exd[i][2])
                rad = str(self.lines_exd[i][3])
                ew = str(np.round(self.lines_ew[i],3))
                err = str(np.round(self.lines_ew_err[i],3))
                current_line = "{0:14s}{1:11s}{2:8s}{3:15s}{4:17s}{5:10s}{6:5s}\n".format(wave,elmnt,ep,gf,rad,ew,err)
                if self.lines_check_flag[i]:
                    flagged_doc.write(current_line.rstrip('\n') + '   # ' + self.lines_flag_reasons[i] + '\n')
                else:
                    doc.write(current_line)
            else:
                removed_lines.append(self.lines[i])
        doc.close()
        flagged_doc.close()
        return np.array(removed_lines)

    def measure_ew(self, i, order, plot = False, ex_params = [0,0,0,0], save_plot = False, window_size = 1.5):
        #extra parameters [0] - shift continuum
        #                 [1] - left boundary in Angstroms
        #                 [2] - right boundary in Angstroms
        #                 [3] - line center in Angstroms
        norm = 1.0
        wind, found_line, line_bound,dy = get_line_window(self.lines[i],self.shifted_wavelength[order],self.normalized_flux[order],ex_params[1],ex_params[2],ex_params[3], window_size)

        if [wind,found_line,line_bound,dy] == [1,1,0,0]:
            pass
            if order != 0:
                concat_wave = np.concatenate((self.shifted_wavelength[order-1],self.shifted_wavelength[order]), axis=None)
                concat_flux = np.concatenate((self.normalized_flux[order-1],self.normalized_flux[order]), axis=None)
                concat_err = np.concatenate((self.obs_err[order-1],self.obs_err[order]), axis=None)
                concat_pred = np.concatenate((self.pred_all[order-1],self.pred_all[order]), axis=None)
                wind, found_line, line_bound, dy = get_line_window(self.lines[i], concat_wave, concat_flux,ex_params[1],ex_params[2],ex_params[3], window_size)
                measure_x_array = concat_wave[wind]
                measure_y_array = concat_flux[wind]
                temp_err_array = concat_err[wind]
                temp_pred_array = concat_pred[wind]
        elif [wind,found_line,line_bound,dy] == [0,0,1,1]:
            pass
            if self.shifted_wavelength[order] != self.shifted_wavelength[-1]:
                concat_wave = np.concatenate((self.shifted_wavelength[order],self.shifted_wavelength[order+1]), axis=None)
                concat_flux = np.concatenate((self.normalized_flux[order],self.normalized_flux[order+1]), axis=None)
                concat_err = np.concatenate((self.obs_err[order],self.obs_err[order+1]), axis=None)
                concat_pred = np.concatenate((self.pred_all[order],self.pred_all[order+1]), axis=None)
                wind, found_line, line_bound, dy = get_line_window(self.lines[i], concat_wave, concat_flux,ex_params[1],ex_params[2],ex_params[3], window_size)
                measure_x_array = concat_wave[wind]
                measure_y_array = concat_flux[wind]
                temp_err_array = concat_err[wind]
                temp_pred_array = concat_pred[wind]
        else:
            measure_x_array = self.shifted_wavelength[order][wind]
            measure_y_array = self.normalized_flux[order][wind]
            temp_err_array = self.obs_err[order][wind]
            temp_pred_array = self.pred_all[order][wind]

        # record the actually-identified line center -- check_for_flags()
        # compares this against the rest wavelength to catch likely
        # misidentification (see its docstring)
        self.lines_found_position[i] = found_line

        other_than_line = np.where((measure_x_array <= line_bound[0])|(measure_x_array >= line_bound[1]))
        only_line = np.where((measure_x_array >= line_bound[0])|(measure_x_array <= line_bound[1]))
        flat_wing = measure_y_array.copy() + ex_params[0]
        flat_wing[other_than_line] = norm
        #highlight points within errors of continuum (or 1.0)
        upper_cont_bounds = measure_y_array+ ex_params[0] + 2*temp_err_array/temp_pred_array
        lower_cont_bounds = measure_y_array+ ex_params[0] - 2*temp_err_array/temp_pred_array
        points_within_norm = np.where((norm > lower_cont_bounds)&(norm < upper_cont_bounds))
        #GP fit
        xtest = np.linspace(measure_x_array[0], measure_x_array[-1], len(measure_x_array))
        m,C=Pred_GP(SEKernel,[1,100],measure_x_array,flat_wing,2*temp_err_array/temp_pred_array, xtest)
        try:
            samples = multivariate_normal(m,C,500)
        except:
            print('SVD did not converge, setting samples to 0')
            print('If line is close to an edge, try remeasuring line with a smaller window size')
            samples = np.array([0]*500)

        m_plot=m.copy()
        m = (-1)*(m-1)
        samp_ew = np.zeros(len(samples))
        a_values = np.zeros(len(samples))
        mu_values = np.zeros(len(samples))
        sig_values = np.zeros(len(samples))
        base_values = np.zeros(len(samples))
        simp_values = np.zeros(len(samples))
        plot_gaussian = False
        for j in range(len(samples)):
            #plt.plot(xtest,(-1)*(samples[j]-1), 'b--')
            bf, err, p0 = gfit_simple(xtest, (-1)*(samples[j]-1), found_line, 0.5,0)
            #print('best fit:', bf, err, p0)
            #if bf[0] > 0.0:
            if abs(gauss_ew(bf[0], bf[2]*2.355)) > 2 and abs(gauss_ew(bf[0], bf[2]*2.355)) < 200:
                samp_ew[j] = abs(gauss_ew(bf[0], bf[2]*2.355))
                a_values[j] = bf[0]
                mu_values[j] = bf[1]
                sig_values[j] = abs(bf[2])
                base_values[j] = bf[3]

                #simpson's rule integration
                # line_inpterp = interp1d(measure_x_array, flat_wing, kind='linear', bounds_error = False)
                # x = np.linspace(flat_wing[0], flat_wing[-1],100)
                # result_y = line_inpterp(x)

                #y =  gauss_model(x,bf[0],bf[1],bf[2],abs(bf[3]))-abs(bf[3])
                simp_values[j] = simpson((-1)*(samples[j]-1), xtest)*1000 #integrates each sample data
            else:
                samp_ew[j] = 0
                a_values[j] = 0
                mu_values[j] = 0
                sig_values[j] = 0
                base_values[j] = 0
                simp_values[j] = 0

        best_bf = np.array([a_values[np.where(a_values!=0)].mean(),mu_values[np.where(mu_values!=0)].mean(),sig_values[np.where(sig_values!=0)].mean(),base_values[np.where(base_values!=0)].mean()])
        fit_gauss = gauss_model(xtest,best_bf[0],best_bf[1],best_bf[2],best_bf[3])*(-1)+1
        #set values for line

        diff = (fit_gauss[only_line] - flat_wing[only_line])**2
        self.lines_gauss_Xsquare[i] = np.sum(diff)


        #self.lines_gauss_Xsquare[i] = chisquare(fit_gauss[only_line], flat_wing[only_line])[0]
        self.lines_bf_params[i] = best_bf
        if len(samp_ew[np.where(samp_ew==0)]) == len(samples):
            self.lines_ew[i] = 0
            self.lines_ew_err[i] = 0
            self.lines_ew_simp[i] = 0
            self.lines_ew_simp_err[i] = 0
        else:
            self.lines_ew[i] = samp_ew[np.where(samp_ew!=0)].mean()
            self.lines_ew_err[i] = samp_ew[np.where(samp_ew!=0)].std()
            self.lines_ew_simp[i] = simp_values[np.where(simp_values!=0)].mean()
            self.lines_ew_simp_err[i] = simp_values[np.where(simp_values!=0)].std()
        print('line to measure:', ELEMENTS[self.lines_exd[i][0]],self.lines[i], '- Line found:', found_line)
        print('EW:',np.round(self.lines_ew[i],2),u"±",np.round(self.lines_ew_err[i],2), 'simps-int:', np.round(self.lines_ew_simp[i],2),u"±", np.round(self.lines_ew_simp_err[i],2))

        #Plotting stuff
        if plot:
            fig = plt.figure(figsize=(12,5))
            fig.suptitle("Order: " + str(order) + " " + "(" + str(np.round(self.shifted_wavelength[order].min(),3)) + "-" + str(np.round(self.shifted_wavelength[order].max(),3)) + ")")
            fit_view = fig.add_subplot(121)
            fit_view.grid()
            fit_view.set_xlabel(r'$\rm Wavelength~(\AA)$', size = 14)
            fit_view.set_ylabel('Normalized Flux', size = 14)
            fit_view.errorbar(measure_x_array,measure_y_array + ex_params[0],
                 yerr=2*temp_err_array/temp_pred_array,capsize=0,fmt='.', color = 'k', label = 'cont', zorder = 2)
            fit_view.scatter(measure_x_array[points_within_norm],measure_y_array[points_within_norm] + ex_params[0], s = 10, c='#4daf4a', zorder = 3, alpha = 0.8)
            fit_view.fill_between(xtest,m_plot+2*np.sqrt(np.diag(C)),
                     m_plot-2*np.sqrt(np.diag(C)),color='#999999',alpha=0.5)
            fit_view.plot([self.lines[i],self.lines[i]],[norm,norm*0.95], '--', color = 'k', alpha = 0.75)
            fit_view.plot([found_line,found_line],[norm,norm*0.95], '-', color='k')
            fit_view.plot([line_bound[0],line_bound[0]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            fit_view.plot([line_bound[1],line_bound[1]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            fit_view.annotate(str(self.lines[i]), xy = [self.lines[i], norm*1.025])
            fit_view.plot(xtest, fit_gauss, '--', color = '#377eb8', lw= 2)
            fit_view.plot([xtest[0],xtest[-1]],[norm,norm], '--', color = '#4daf4a')

            data_view = fig.add_subplot(122)
            data_view.grid()
            data_view.set_xlabel(r'$\rm Wavelength~(\AA)$', size = 14)
            data_view.scatter(measure_x_array,measure_y_array+ ex_params[0], s = 5, c = 'k', zorder = 2)
            data_view.errorbar(measure_x_array,measure_y_array + ex_params[0],
                 yerr=2*temp_err_array/temp_pred_array,capsize=0,fmt='.', color = 'k', zorder = 3, alpha = 0.5)
            plt.tight_layout()


            # fig1, coarse_view = plt.subplots()
            # coarse_view.set_title("Order: " + str(order) + " " + "(" + str(np.round(self.shifted_wavelength[order].min(),3)) + "-" + str(np.round(self.shifted_wavelength[order].max(),3)) + ")")
            # coarse_view.grid()
            # coarse_view.set_xlabel(r'$\rm Wavelength~(\AA)$', size = 14)
            # coarse_view.set_ylabel('Normalized Flux', size = 14)
            #coarse_view.plot(xtest,m_plot, 'k--', alpha = 0.75)
            # coarse_view.errorbar(self.shifted_wavelength[order][wind],self.normalized_flux[order][wind] + ex_params[0],
            #      yerr=2*self.obs_err[order][wind]/self.pred_all[order][wind],capsize=0,fmt='.', color = 'k', label = 'cont', zorder = 2)
            # coarse_view.scatter(self.shifted_wavelength[order][wind][points_within_norm],self.normalized_flux[order][wind][points_within_norm] + ex_params[0], s = 10, c='#4daf4a', zorder = 3, alpha = 0.8)
            # coarse_view.fill_between(xtest,m_plot+2*np.sqrt(np.diag(C)),
            #          m_plot-2*np.sqrt(np.diag(C)),color='#999999',alpha=0.5)
            #coarse_view.plot(xtest,samples.T,alpha=0.1, color='#cccccc')
            # coarse_view.plot([self.lines[i],self.lines[i]],[norm,norm*0.95], '--', color = 'k', alpha = 0.75)
            # coarse_view.plot([found_line,found_line],[norm,norm*0.95], '-', color='k')
            # coarse_view.plot([line_bound[0],line_bound[0]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            # coarse_view.plot([line_bound[1],line_bound[1]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            # coarse_view.annotate(str(self.lines[i]), xy = [self.lines[i], norm*1.025])
            #coarse_view.plot(xtest,dy+norm, '--', color = '#e41a1c', lw = 2) #view gradient
            # if plot_gaussian:
            #     coarse_view.plot(xtest, fit_gauss, '--', color = '#377eb8', lw= 2)
            # coarse_view.plot([xtest[0],xtest[-1]],[norm,norm], '--', color = '#4daf4a')
            if save_plot:
                fig_title = ELEMENTS[self.lines_exd[i][0]] + '_' + str(self.lines[i]) + '_' + str(order) + '.pdf'
                plt.savefig('line_plots/'+fig_title)
            plt.show()

            print('#-----------------------#')

        #print extra parameter stuff
        if ex_params == [0,0,0,0]:
            pass
        else:
            self.lines_exp[i] = np.array(ex_params)
            print('extra params:',ex_params)

    def measure_all_ew(self, exclude_lines= [], plot_lines=[], ex_params = {}, window_size = 1.5, save_all = False):
        if save_all:
            make_plots_folder()

        for order in range(len(self.wavelength)):
            for i in range(len(self.lines)):
                if self.lines[i] in exclude_lines:
                    self.lines_ew[i] = 0.0
                    self.lines_ew_simp[i] = 0.0
                    self.lines_gauss_Xsquare[i] = np.nan
                    self.lines_bf_params[i] = None
                    self.lines_ew_err[i] = np.nan
                    self.lines_ew_simp_err[i] = np.nan
                    self.lines_exp[i] = [0,0,0,0]
                    #self.lines_check_flag[i] = False
                elif self.lines[i] >= self.shifted_wavelength[order][0] and self.lines[i] <= self.shifted_wavelength[order][-1]:
                    plot = False
                    exp = [0,0,0,0]
                    if self.lines[i] in plot_lines:
                        plot = True
                        if self.lines[i] in ex_params.keys():
                            exp = ex_params[self.lines[i]]
                    if save_all:
                        self.measure_ew(i,order, plot, exp, True, window_size)
                    else:
                        self.measure_ew(i,order, plot, exp, False, window_size)
        #self.lines_bf_params = np.array(self.lines_bf_params)

    def measure_line_ew(self,line,ex_params=[0,0,0,0], save_line = False, save_plot = False, window_size = 1.5):
        if save_plot:
            make_plots_folder()
        i = np.where(self.lines == line)[0][0]
        found = False
        for order in range(len(self.wavelength)):
            if self.lines[i] >= self.shifted_wavelength[order][0] and self.lines[i] <= self.shifted_wavelength[order][-1]:
                if not found:
                    self.lines_ew[i] = 0.0
                    self.lines_ew_simp[i] = 0.0
                    self.lines_gauss_Xsquare[i] = np.nan
                    self.lines_bf_params[i] = None
                    self.lines_ew_err[i] = np.nan
                    self.lines_ew_simp_err[i] = np.nan
                    #self.lines_check_flag[i] = False
                    self.measure_ew(i,order, True, ex_params, save_plot, window_size)
                    found = True
                    if save_line:
                        with open('line_'+str(line)+'.txt','w') as f:
                            for k in range(len(self.shifted_wavelength[order])):
                                f.write("{0:10f}\t{1:10f}\n".format(self.shifted_wavelength[order][k],self.normalized_flux[order][k]))
                            print('order', order, 'saved!')

    def check_for_flags(self):
        """
        Flag lines whose measurement looks untrustworthy: high EW error
        fraction, too shallow to trust, a poor Gaussian fit, or -- the
        check that matters most for avoiding a silently WRONG EW rather
        than just an imprecise one -- a found line center
        (lines_found_position, set by measure_ew()) that lands more than
        position_thresh away from the line's rest wavelength. Since
        get_line_window()'s search is itself bounded to +/-0.1 A by
        default, a found position near that boundary is a strong sign the
        code locked onto a neighboring feature/blend/noise dip rather than
        the intended line -- see the wavelength-shift robustness testing in
        conversation/session notes for how large this risk can be when a
        spectrum's wavelength correction is not well constrained. Also
        flags a line landing inside a disputed order-overlap range, if
        flag_order_overlaps() was called first (self.overlap_flag_ranges
        -- see overlap_check.py's module docstring) -- a continuum-
        placement problem this check catches even when nothing about the
        line's OWN fit looks wrong (e.g. a real order-edge droop found on
        GRACES this way, invisible to every check above). Also flags a
        line landing inside a whole order flagged by flag_bad_orders()
        (self.bad_order_ranges -- see outlier_check.py's module
        docstring), an extreme raw-flux outlier (cosmic ray/detector
        defect) confirmed to corrupt an entire order's continuum fit.
        Also EXCLUDES a line whose measurement window contains a pixel
        corrected by correct_spikes()/combine_spectra() (self.
        corrected_pixel_wavelengths) -- its measurement reflects a
        corrected, not originally observed, value, so it's excluded
        outright rather than just flagged on an uncorrected one. If
        identify_lines() was called first (self.lines_id_run -- see
        line_identification.py's module docstring), also flags: (a) a
        line identify_lines() never found a significant absorption
        feature for at all (lines_id_detected False) -- distinct from
        and prior to the position-offset check above, since a line
        get_line_window() "finds" is never allowed to silently be pure
        noise dressed up as a detection; (b) a line identify_lines()
        flagged as blended (a second, competitive candidate nearby).
        """
        for i in range(len(self.lines)):
            self.lines_check_flag[i] = False
            reasons = []
            #error measure check - above 10% is a problem
            if self.lines_ew_err[i]/self.lines_ew[i] >= .1:
                self.lines_check_flag[i] = True
                reasons.append(f'error>10% ({np.round(self.lines_ew_err[i],2)} mA)')
                print(self.lines[i], 'has more than a 10% error', np.round(self.lines_ew_err[i],2))
            #shallow line check - below 2 mA is a problem
            if self.lines_ew[i] < 2.0:
                self.lines_check_flag[i] = True
                reasons.append(f'too shallow ({np.round(self.lines_ew[i],2)} mA)')
                print(self.lines[i], 'might be too shallow', np.round(self.lines_ew[i],2))
            if np.isnan(self.lines_gauss_Xsquare[i]):
                self.lines_check_flag[i] = True
                reasons.append('no fit / line not found')
                print(self.lines[i], 'no fit, line may not be found in spectrum')
            #X square check - if fit above threshold (problem)
            elif self.lines_gauss_Xsquare[i] > self.X_thresh:
                self.lines_check_flag[i] = True
                reasons.append(f'bad fit (X^2={np.round(self.lines_gauss_Xsquare[i],4)})')
                print(self.lines[i], 'might have a bad fit', self.lines_gauss_Xsquare[i])
            #position check - found line far from its rest wavelength is a
            #likely misidentification, not just an imprecise measurement
            if not np.isnan(self.lines_found_position[i]):
                position_offset = abs(self.lines_found_position[i] - self.lines[i])
                if position_offset > self.position_thresh:
                    self.lines_check_flag[i] = True
                    reasons.append(f'position off by {np.round(position_offset*1000,1)} mA (possible misidentification)')
                    print(self.lines[i], 'found position is', np.round(position_offset*1000,1),
                          'mA from rest wavelength -- possible misidentification, inspect before trusting this EW')
            #order-overlap check - line sits in a wavelength range where
            #two orders' normalized flux disagreed badly enough to
            #distrust either one there (see flag_order_overlaps())
            for rng in (self.overlap_flag_ranges or []):
                if rng['wave_lo'] <= self.lines[i] <= rng['wave_hi']:
                    self.lines_check_flag[i] = True
                    reasons.append(f"order-overlap disagreement ({rng['median_pct']:+.1f}%, "
                                    f"orders {rng['order_i']}/{rng['order_j']})")
                    print(self.lines[i], 'sits in a disputed order-overlap range '
                          f"(orders {rng['order_i']}/{rng['order_j']}, {rng['median_pct']:+.1f}%)")
            #bad-order check - line sits in an order flagged for an
            #extreme raw-flux outlier (cosmic ray/detector defect) that
            #would corrupt the whole order's continuum fit (see
            #flag_bad_orders())
            for rng in (self.bad_order_ranges or []):
                if rng['wave_lo'] <= self.lines[i] <= rng['wave_hi']:
                    self.lines_check_flag[i] = True
                    reasons.append(f"order {rng['order']} has an extreme raw-flux outlier "
                                    f"({rng['outlier_factor']:.0f}x robust scale)")
                    print(self.lines[i], f"sits in order {rng['order']}, flagged for an extreme "
                          f"raw-flux outlier ({rng['outlier_factor']:.0f}x robust scale)")
            #corrected-pixel check - a spike (cosmic ray/sky-emission
            #contamination/bad pixel) was detected and CORRECTED within
            #this line's own measurement window (same +/-0.1 A default
            #search radius as get_line_window()) -- excluded outright
            #rather than just flagged, since the measurement there
            #reflects a corrected, not originally observed, value (see
            #correct_spikes()/combine_spectra())
            for cw in (self.corrected_pixel_wavelengths or []):
                if abs(cw - self.lines[i]) <= 0.1:
                    self.lines_check_flag[i] = True
                    reasons.append(f'corrected pixel within measurement window ({cw:.3f} A)')
                    print(self.lines[i], f'has a corrected pixel within its measurement window '
                          f'({cw:.3f} A) -- excluded, not just flagged on an uncorrected value')
            #identification checks - only meaningful once identify_lines()
            #has actually run (lines_id_detected defaults to False either
            #way, so this must be gated on lines_id_run to avoid flagging
            #every line as undetected when it simply hasn't run yet)
            if self.lines_id_run:
                if not self.lines_id_detected[i]:
                    self.lines_check_flag[i] = True
                    reasons.append('not detected during identification (no significant '
                                    'absorption feature found near rest wavelength)')
                    print(self.lines[i], 'was not detected during identify_lines() -- '
                          'too weak for this spectrum\'s S/N, or genuinely absent')
                elif self.lines_id_blended[i]:
                    self.lines_check_flag[i] = True
                    reasons.append(f'identified as blended (sig={self.lines_id_significance[i]:.1f}, '
                                    'a competing candidate sits nearby)')
                    print(self.lines[i], 'identified as blended -- a competing candidate '
                          'sits close enough to be comparably significant')
            self.lines_flag_reasons[i] = '; '.join(reasons)

    def check_spectra(self, norm=True, lines=False):
        orders = len(self.wavelength)

        if norm:
            for i in range(orders):
                plt.plot(self.shifted_wavelength[i],self.normalized_flux[i], linewidth = 0.5)
            if lines:
                for j in range(len(self.lines)):
                    plt.plot([self.lines[j],self.lines[j]],[0.95,1.0], 'k')
                    plt.annotate(str(self.lines[j]), xy=[self.lines[j],1.01])
            #plt.ylim([0.4,1.1])
            plt.show()
        else:
            for i in range(orders):
                fig = plt.figure(figsize=(8,4))
                ax = fig.add_subplot(111)
                ax.set_title('Order: '+str(i))
                ax.plot(self.wavelength[i],self.flux[i])
                #ax.set_xticks([])
                #ax.set_yticks([])
                #ax.set_ylabel(str(i+1))
                plt.tight_layout()
                plt.show()

    def save_object(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
