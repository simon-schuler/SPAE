"""Spectrum_Data: the core driver class for loading a spectrum, fitting its
continuum, wave-shifting against a reference, and measuring line EWs."""

import copy
import glob
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from .constants import ELEMENTS
from .continuum import fit_als_continuum
from .line_profile import (get_line_window, gauss_model, gfit_direct, gauss_ew, gauss_ew_err,
                            gauss_model_err, estimate_local_continuum,
                            estimate_local_continuum_sloped, EW_K)
from .combine import make_line, parabolic_refine
from .plotting import make_plots_folder
from .readers import read_spectrum
from .radial_velocity import measure_effective_rv, measure_rv_ccf, C_KMS
from .response_correction import apply_response_correction as _apply_response_correction
from .overlap_check import check_order_overlaps as _check_order_overlaps
from .overlap_check import flagged_overlap_ranges as _flagged_overlap_ranges
from .line_identification import identify_lines_in_spectrum as _identify_lines_in_spectrum
from .reference_atlas import (load_reference_atlas as _load_reference_atlas,
                               estimate_resolving_power, reference_continuum_mask)


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

        #reference-atlas cross-check (see reference_atlas.py) -- unset
        #until load_reference_atlas() is called; measure_ew() falls back
        #to its existing behavior while these are None
        self.ref_wave = None
        self.ref_flux = None
        self.ref_resolving_power = None

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
            # dtype=float (not the bare-int default from np.full(n, 0)) --
            # an int array here silently truncates any fractional
            # continuum value to 0, and wraps a NaN/Inf fit result (e.g.
            # from a singular AsLS solve, see fit_als_continuum()) to
            # int64's sentinel min (-9223372036854775808) instead of
            # propagating it as NaN. Both confirmed on real data: the
            # bundled suni.fits sample (order 9, a singular-matrix case)
            # and HD_102071's faintest blue orders (near-zero real
            # continuum, truncated to exact 0 -> normalized_flux=inf).
            self.pred_all[i] = np.full(len(self.wavelength[i]), 0, dtype=float)
            self.pred_var_all[i] = np.full(len(self.wavelength[i]), 0, dtype=float)
            #abs() guards against the occasional slightly-negative pixel
            #from background/bias subtraction (confirmed on the bundled
            #suni.fits sample: 4/4021 pixels in one order, ~-20 counts
            #against a ~150000 count median -- ordinary noise, not a real
            #data problem) -- sqrt() of a small negative value is NaN,
            #which silently breaks the downstream weighted continuum fit
            #(a singular-matrix crash, not a graceful failure). The
            #np.maximum(..., 1.0) floor guards the same fit against an
            #exact-zero-count (dead/masked) pixel instead making its
            #weight (1/err) infinite -- also confirmed on suni.fits: 146
            #pixels, non-contiguous, near one order's far edge.
            self.obs_err[i] = np.maximum(np.sqrt(np.abs(self.flux[i])), 1.0)
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
        #line - equivalent width, from the GLOBAL-continuum fit (bf_global
        #-- continuum pinned at norm, no local wing-based correction).
        #This is the REPORTED EW (what make_ew_doc()/check_for_flags()
        #use) -- see measure_ew()'s comment where bf_global is computed
        #for why the global fit, not the local-continuum-corrected one,
        #is the default.
        self.lines_ew = None
        #line - equivalent width error, GLOBAL-continuum fit
        self.lines_ew_err = None
        #line - equivalent width from the LOCAL-continuum-corrected fit
        #(best_bf/cont_offset -- estimate_local_continuum()'s flat,
        #photon-noise-weighted bias applied), set by measure_ew() whenever
        #fit_continuum=True (equal to lines_ew when fit_continuum=False,
        #since best_bf/bf_global are then the same fit). Kept as a
        #diagnostic/comparison value only -- NOT what's reported in
        #make_ew_doc()'s linelist or checked by check_for_flags().
        self.lines_ew_local = None
        #line - equivalent width error, LOCAL-continuum-corrected fit --
        #see lines_ew_local
        self.lines_ew_err_local = None
        #line - best fit parameters for gaussian fit
        self.lines_bf_params = None
        #line - sum-of-squared-residuals between the observed line core
        #and a Gaussian fit against the GLOBAL continuum (norm, no local
        #wing-based correction) -- see measure_ew()'s comment where this
        #is set for why the global fit, not the local-continuum-corrected
        #one, is what's checked here
        self.lines_gauss_Xsquare = None
        #line - X squared threshold value
        self.X_thresh = 0.003
        #line - the actually-identified line center (Angstrom, in
        #shifted_wavelength), set by measure_ew(). Comparing this against
        #the rest wavelength (self.lines) is how check_for_flags() catches
        #likely misidentification -- see its docstring.
        self.lines_found_position = None
        #line - per-side (blue/red wing) booleans from
        #estimate_local_continuum_sloped()'s internal near-vs-far trend
        #test, always set by measure_ew() regardless of fit_continuum (see
        #its docstring for why this diagnostic isn't gated by that flag).
        #check_for_flags() flags a line where BOTH are True -- see its
        #docstring and measure_ew()'s comment where these are set.
        self.lines_internal_trend_blue = None
        self.lines_internal_trend_red = None
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
        #line - True if a human reviewer explicitly decided to keep this
        #line despite check_for_flags() otherwise flagging it (see
        #pipeline.py's interactive review step). check_for_flags() always
        #recomputes lines_check_flag from scratch on every call, including
        #the one make_ew_doc() makes internally -- this is the one piece
        #of per-line state that check_for_flags() itself won't overwrite,
        #so a review decision actually survives to the final output
        #instead of being silently recomputed away.
        self.lines_human_keep = None
        #wavelength ranges where two orders' overlap disagreed badly
        #enough to distrust either one there -- set by
        #flag_order_overlaps(), consulted by check_for_flags(). Empty
        #(no-op) until flag_order_overlaps() is called.
        self.overlap_flag_ranges = []

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

        Each order's outer `edge_ignore_aa` Angstroms (default 2.0, also
        passed through via **als_kwargs) are excluded from influencing
        the fit -- real order edges are where detector/blaze artifacts
        concentrate (confirmed on a real Keck order), and this default
        avoids the fit getting dragged by one. A continuum value is still
        produced there (extrapolated from the trusted interior), so nothing
        downstream sees a gap. See fit_als_continuum()'s docstring for the
        real tradeoff this involves and when to lower it toward 0.

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
        # obs_err[order][clipl:clipr] is already `err` (read from self.obs_err
        # above, not recomputed) -- nothing to write back
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

    def combine_spectra(self, spectB, resolution = 1000, shift=True):
        print('Use self.update_combined() when you are happy with the combined flux to override self.flux')
        #find corresponding orders that match self in spectB
        med_A = [np.median(self.shifted_wavelength[i]) for i in range(len(self.shifted_wavelength))]
        med_B = [np.median(spectB.shifted_wavelength[i]) for i in range(len(spectB.shifted_wavelength))]
        b_order = []
        #find corresponding B order
        for k in range(len(med_A)):
            diff = abs(med_A[k] - med_B)
            loc = np.where(diff == diff.min())
            b_order.append(loc)

        combined_flux_orders = np.zeros_like(self.flux)
        if shift:
            #first shift A to match B with higher accuracy (higher resolution)
            #may have to include a try statement for errors
            self.estimate_shift([spectB], shift_spacing=resolution)
            self.clean_shift()

        #loop through orders
        for i in range(len(self.shifted_wavelength)):

            #combining flux values for each wavelength value
            combined_flux = np.zeros(len(self.shifted_wavelength[i]))
            combined_err = np.zeros(len(self.shifted_wavelength[i]))

            print('A order', i, 'B order', b_order[i][0][0])

            #loop through each shifted wavelength value
            for j in range(len(self.shifted_wavelength[i])):
                #difference between one shifted wavelength value and all B wavelength values
                #element closest to zero is location of closest wavelength values
                ed = 5
                if j < ed:
                    le = 0
                    re = j + ed
                elif j > len(self.shifted_wavelength[i])-ed:
                    le = j - ed
                    re = len(self.shifted_wavelength[i]) + 1
                else:
                    le = j - ed
                    re = j + ed
                near_point = spectB.wavelength[b_order[i][0][0]][le:re]
                diff_array = abs(self.shifted_wavelength[i][j] - near_point)
                loc = np.where(diff_array == diff_array.min())
                #add A flux with B flux at location where diff = 0
                combined_flux[j] = self.flux[i][j] + spectB.flux[b_order[i][0][0]][le:re][loc]
                #errors add in quadrature (independent measurements), using
                #each spectrum's own already-correct obs_err rather than
                #recomputing sqrt(combined_flux) -- which would be wrong
                #for the same reason it's wrong in normalize(): if either
                #spectrum has already been response-corrected, its flux is
                #no longer a raw Poisson count, so sqrt() of it is not its
                #true sigma (see apply_response_correction()'s docstring)
                combined_err[j] = np.sqrt(self.obs_err[i][j]**2
                                           + spectB.obs_err[b_order[i][0][0]][le:re][loc][0]**2)
            #print(combined_flux)
            #collect flux values for each order
            combined_flux_orders[i] = combined_flux
            self.obs_err[i] = combined_err

            plt.plot(self.shifted_wavelength[i], self.flux[i], label = 'A')
            plt.plot(spectB.wavelength[b_order[i][0][0]], spectB.flux[b_order[i][0][0]], label = 'B')
            plt.plot(self.shifted_wavelength[i], combined_flux, label = 'A+B')
            plt.xlim([np.median(self.shifted_wavelength[i]-(self.shifted_wavelength[i].max() - self.shifted_wavelength[i].min())/10),
                np.median(self.shifted_wavelength[i]+(self.shifted_wavelength[i].max() - self.shifted_wavelength[i].min())/10)])
            plt.grid()
            plt.legend()
            plt.show()
            print('#-----------------------#')
        #replace original flux for A with combined flux
        self.combined_flux = combined_flux_orders

    def update_combined(self):
        self.flux = self.combined_flux

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
        cross-correlation itself, e.g. combine_spectra()'s internal use.
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
                       cross_check=True, ccf_v_min=-250.0, ccf_v_max=250.0, ccf_v_step=0.5,
                       ccf_min_significance=5.0, disagreement_kms=2.0):
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

        RV_REFERENCE_LINES has real weaknesses on its own: it can be
        entirely absent from a spectrum whose coverage happens to miss all
        of Ca II H&K/Balmer/Mg b/Na D, mixing Balmer lines with metal
        lines in one average can be actively wrong (confirmed on real
        data, the two families disagreed by ~10 km/s -- a real difference
        in line formation physics, not just noise), and even a single
        named line's own per-line search window can contain more than one
        comparably strong absorption feature in a densely-blended region,
        letting an unweighted Gaussian fit lock onto the wrong one --
        confirmed on a real, genuinely high-velocity star (HD_10383, true
        RV ~+107 km/s): only 2 of 10 reference lines were usable at all,
        and those 2 disagreed with EACH OTHER by ~51 km/s.

        If a science linelist is already loaded (self.lines, via
        load_lines()), this also cross-checks against
        radial_velocity.measure_rv_ccf() -- a cross-correlation against
        the WHOLE linelist (typically dozens of lines) rather than a
        handful of named ones, searched over a wide, continuous velocity
        grid rather than a small per-line Angstrom window. This is used
        as the ONLY estimate if RV_REFERENCE_LINES found nothing usable,
        and PREFERRED over the named-line estimate whenever its
        significance clears `ccf_min_significance` (confirmed on real
        data: correctly recovers both a genuinely large shift the named-
        line method got badly wrong, HD_10383 at ~+107 km/s, and a small,
        already-well-known shift, the Sun at ~-2.6 km/s, from the exact
        same 78-line solar linelist with no per-target tuning -- see
        measure_rv_ccf()'s own docstring for the full validation and why
        it structurally can't hit either failure mode above). An earlier
        version of this cross-check (measure_rv_from_linelist(), still
        available but no longer used here) used a small per-line search
        radius for the same purpose; confirmed on the same HD_10383 case
        that this has its own hard ceiling on the shift it can ever find
        (set by the radius, translated through c/lambda) and silently
        locks onto an unrelated nearby feature instead of failing loudly
        once a real shift exceeds it. Call load_lines() before this if
        you want the cross-check available; it's a silent no-op
        (identical to the old behavior) if no linelist is loaded yet.

        Parameters
        ----------
        rv : float, km/s, optional -- apply this RV directly and skip line
            measurement (e.g. if the RV is already known from elsewhere,
            or confirmed independently -- see measure_rv_ccf()'s docstring
            for why HD_10383 needed this: a wide, unweighted, symmetric
            search window can flip a correct measurement to the wrong
            sign when a second strong feature falls within it, so a
            manual, corroborated check across >1 independent line is
            worth doing before trusting any automated value blindly).
        lines : {name: (rest_wavelength, window)}, optional -- defaults to
            radial_velocity.RV_REFERENCE_LINES.
        min_depth : float -- minimum line depth (in normalized flux) to
            trust a line's fitted center.
        verbose : bool -- print the measured RV(s), which lines were used,
            and the cross-check outcome.
        cross_check : bool -- if False, use RV_REFERENCE_LINES only, same
            as before this parameter existed.
        ccf_v_min, ccf_v_max, ccf_v_step : passed to measure_rv_ccf() as
            v_min/v_max/v_step -- the trial velocity grid searched.
        ccf_min_significance : minimum measure_rv_ccf() peak significance
            to trust/prefer it at all (default 5.0 -- both real validation
            cases scored 17-20, comfortable margin above this floor;
            below it, the CCF found no clean peak, e.g. too few of this
            linelist's lines actually present/detectable in this target).
        disagreement_kms : how far the named-line and CCF estimates must
            differ before printing a note about it (verbose only -- purely
            diagnostic now, doesn't affect which estimate is used; see
            ccf_min_significance for that).

        Returns
        -------
        rv : float, km/s -- the RV actually applied.
        """
        if rv is None:
            measured_rv, rv_err, used = measure_effective_rv(self, lines=lines, min_depth=min_depth)

            ccf_rv, ccf_significance = None, None
            if cross_check and self.lines is not None and len(self.lines) > 0:
                ccf_rv, ccf_significance, _, _ = measure_rv_ccf(
                    self.lines, self.wavelength, self.normalized_flux,
                    v_min=ccf_v_min, v_max=ccf_v_max, v_step=ccf_v_step)
            ccf_usable = ccf_rv is not None and ccf_significance >= ccf_min_significance

            if measured_rv is None and not ccf_usable:
                raise ValueError(
                    "Could not measure an effective RV -- none of the reference "
                    "lines were found/usable in this spectrum's wavelength "
                    "coverage, and the linelist-based cross-correlation (if a "
                    "linelist was even loaded) found no significant peak either "
                    "(load_lines() first to enable that, or pass rv= directly).")

            if verbose and measured_rv is not None:
                print(f"Named-line RV = {measured_rv:.3f} +/- {rv_err:.3f} km/s, "
                      f"from {len(used)} line(s):")
                for name, restw, order, v in used:
                    print(f"  {name} ({restw} A, order {order}): v={v:.3f} km/s")
            disagreement = (abs(measured_rv - ccf_rv)
                             if (measured_rv is not None and ccf_rv is not None) else None)
            if verbose and ccf_rv is not None:
                trust_note = "" if ccf_usable else f" (below ccf_min_significance={ccf_min_significance}, not used)"
                print(f"Linelist cross-correlation RV = {ccf_rv:.3f} km/s "
                      f"(significance={ccf_significance:.1f}){trust_note}")
            if verbose and disagreement is not None and disagreement > disagreement_kms:
                print(f"NOTE: named-line and cross-correlation RVs disagree by "
                      f"{disagreement:.3f} km/s (> {disagreement_kms}) -- likely "
                      f"the named-line fit locked onto a competing nearby feature "
                      f"(see measure_rv_ccf()'s docstring for a real confirmed case).")

            if ccf_usable:
                #preferred whenever it clears the significance floor, whether
                #or not it agrees with the named-line estimate -- confirmed
                #more robust in both directions on real data (see docstring)
                if verbose:
                    print("-> using the cross-correlation RV.")
                #measure_rv_ccf() doesn't (yet) produce a formal statistical
                #error the way the per-line Gaussian fits do. Reuse the
                #named-line estimate's error as a rough stand-in ONLY when
                #the two estimates actually agree -- a disagreement is
                #itself evidence the named-line fit (and, by extension, its
                #own error bar) isn't trustworthy here: confirmed on
                #HD_10383, whose named-line point estimate was ~53 km/s off,
                #so its formal 18.2 km/s error is not a meaningful precision
                #floor for the CCF value being adopted instead. Falls back
                #to a fixed placeholder otherwise, sized to both real
                #validation cases' actual method-to-method scatter (~1-2
                #km/s) rather than either 0 (false precision) or the
                #discredited named-line error (false confidence the wrong
                #way).
                trustworthy_error = (measured_rv is not None and disagreement is not None
                                      and disagreement <= disagreement_kms)
                measured_rv, rv_err = ccf_rv, (rv_err if trustworthy_error else 1.5)
            elif verbose:
                print("-> using the named-line RV.")

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
        self.lines_ew_local = np.zeros(len(self.lines))
        self.lines_ew_err_local = np.zeros(len(self.lines))
        self.lines_bf_params = np.array([None]*len(self.lines))
        self.lines_gauss_Xsquare = np.array([np.nan]*len(self.lines))
        self.lines_found_position = np.array([np.nan]*len(self.lines))
        self.lines_internal_trend_blue = np.array([False]*len(self.lines))
        self.lines_internal_trend_red = np.array([False]*len(self.lines))
        self.lines_check_flag = np.array([False]*len(self.lines))
        self.lines_flag_reasons = np.array(['']*len(self.lines), dtype=object)
        self.lines_human_keep = np.array([False]*len(self.lines))
        #line - which order's measurement is the one kept in lines_ew[i]
        #etc, and how far (Angstroms) that line sat from that order's
        #nearer edge -- set by measure_all_ew() when the same line falls
        #in more than one order (echelle order overlap); see its
        #docstring. None/NaN if the line hasn't been measured yet, or was
        #only ever seen in one order (no overlap to resolve).
        self.lines_order = np.array([None]*len(self.lines))
        self.lines_edge_distance = np.array([np.nan]*len(self.lines))
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
                d0 = '0.0'
                ew = str(np.round(self.lines_ew[i],3))
                err = str(np.round(self.lines_ew_err[i],3))
                # 8 fixed 10-char columns: wave1,atom1,e,gf,dampnum,d0,width,width_err
                # (the first 7 match MOOG's native 7e10.3 linelist format; the
                # trailing err column is a pymoog-only extension read by
                # moog/inlines.py -- real MOOG ignores anything past column 70)
                current_line = "{0:10s}{1:10s}{2:10s}{3:10s}{4:10s}{5:10s}{6:10s}{7:10s}\n".format(
                    wave, elmnt, ep, gf, rad, d0, ew, err)
                if self.lines_check_flag[i]:
                    flagged_doc.write(current_line.rstrip('\n') + '   # ' + self.lines_flag_reasons[i] + '\n')
                elif self.lines_human_keep[i] and self.lines_flag_reasons[i]:
                    #kept in the main linelist by an explicit human review
                    #decision, not because nothing looked wrong -- note why
                    #it would otherwise have been flagged, so this isn't
                    #silently indistinguishable from an uncontroversial line
                    doc.write(current_line.rstrip('\n') + '   # human-reviewed, kept despite: '
                              + self.lines_flag_reasons[i] + '\n')
                else:
                    doc.write(current_line)
            else:
                removed_lines.append(self.lines[i])
        doc.close()
        flagged_doc.close()
        return np.array(removed_lines)

    def load_reference_atlas(self, path, resolving_power=None):
        """Load an independent, high-S/N reference spectrum (e.g. the
        Kurucz solar flux atlas, fluxspliced.2005) that measure_ew() will
        cross-check a line's wing against -- see reference_atlas.py's
        module docstring for why: a real but shallow, gradual blend (not
        a sharp outlier) can bias estimate_local_continuum()'s median/MAD
        clip, and an independent high-S/N reference can catch that where
        this spectrum's own noise can't.

        Only useful when a suitable reference actually exists for this
        target (so far: the Sun) -- do not load a reference atlas of a
        different star. Nothing else changes if this is never called;
        measure_ew() falls back to its current behavior.

        Call load_lines() before this if resolving_power is left as None:
        the empirical estimate default-samples the already-loaded
        linelist (see reference_atlas.estimate_resolving_power()).

        resolving_power : R = lambda/FWHM for THIS spectrum (not the
            atlas). None (default) estimates it empirically from a robust
            LOW percentile of fitted width across the loaded linelist's
            weak-to-moderate, unsaturated lines (see
            reference_atlas.estimate_resolving_power()'s docstring for why
            a low percentile, not the median) -- pass an explicit value
            instead if you already know your instrument's R, since that's
            ordinarily a known setup property rather than something to
            infer from a handful of noisy fits.
        """
        self.ref_wave, self.ref_flux = _load_reference_atlas(path)
        if resolving_power is not None:
            self.ref_resolving_power = resolving_power
        else:
            self.ref_resolving_power = estimate_resolving_power(self)
            if self.ref_resolving_power is None:
                print('load_reference_atlas: could not empirically estimate a resolving '
                      'power (fewer than 5 usable lines) -- pass resolving_power '
                      'explicitly instead. Reference atlas NOT loaded.')
                self.ref_wave, self.ref_flux = None, None

    def measure_ew(self, i, order, plot = False, ex_params = [0,0,0,0], save_plot = False, window_size = 1.5, show_plot = True, fit_continuum = False, auto_widen = True, widen_window_size = 2.5, slope_sig_thresh = 3.0, slope_min_points = 5, slope_internal_sig_thresh = 3.0, plot_window_size = None, keep_result = True):
        #extra parameters [0] - shift continuum
        #                 [1] - left boundary in Angstroms
        #                 [2] - right boundary in Angstroms
        #                 [3] - line center in Angstroms
        #show_plot: set False to save/build the figure without blocking on
        #plt.show() -- used by measure_all_ew(save_all=True) so a QC plot
        #for every line doesn't pop up (and need closing) one at a time
        #fit_continuum: False (DEFAULT): assume the global normalization
        #already put this window's continuum exactly at norm -- the
        #reported EW (self.lines_ew) always comes from this GLOBAL-
        #continuum fit regardless of this flag (see bf_global below); this
        #flag only controls whether an ADDITIONAL, local-continuum-
        #corrected comparison fit is also computed. True: additionally
        #estimate a local flat continuum level (c0) from this line's own
        #wing data (see estimate_local_continuum()) and fit against that
        #instead -- stored as lines_ew_local/lines_ew_err_local, a
        #diagnostic/comparison value only, kept in the codebase but no
        #longer the default after real-star testing (HD_10383) showed the
        #GLOBAL-continuum fit was actually the more reliable default for
        #EW reporting (see check_for_flags()'s global-chi2 check, and
        #xspect-ew-continuum-flagging session notes) -- letting the local
        #continuum vary risked absorbing real blending into what should be
        #a continuum correction. The sloped diagnostic just below (which
        #feeds check_for_flags()'s "no flat continuum" flag) runs
        #UNCONDITIONALLY regardless of this setting -- it isn't gated by
        #fit_continuum any more, precisely so that flag keeps working with
        #the new global-only default.
        #slope_sig_thresh, slope_min_points: passed to
        #estimate_local_continuum_sloped() (see its docstring) -- always
        #computed (see above), reported/plotted alongside the reported EW,
        #but still does not itself drive the reported EW.
        #auto_widen: automatically retry once with widen_window_size if
        #this line's fit quality is poor at window_size -- see the retry
        #check below, right after ew/ew_err are finalized, for exactly
        #what "poor" means and why. Applies regardless of fit_continuum.
        #plot_window_size: PLOTTING ONLY -- show this much more
        #surrounding data/fit-curve extrapolation than the actual fit
        #window (window_size) used, for visual context (e.g. judging the
        #global continuum normalization against a wider stretch of the
        #order). None (default): plot exactly the fit window, unchanged
        #from before this parameter existed. Never affects the fit,
        #EW, or error -- those still only ever see window_size's data.
        #keep_result: True (default) writes this measurement into
        #self.lines_*[i] as usual. Set False to still run the full
        #measurement (fit, print, and -- if requested -- plot/save) but
        #WITHOUT overwriting self.lines_*[i] -- used by measure_all_ew()
        #when the same line falls in more than one order (echelle order
        #overlap): every order it appears in still gets measured and
        #logged/plotted for inspection, but only the measurement from the
        #order where the line sits FURTHEST from that order's edge is
        #kept as the one self.lines_ew[i]/check_for_flags()/make_ew_doc()
        #actually see -- an edge measurement is the one most exposed to
        #exactly the kind of order-boundary artifacts investigated this
        #session (see continuum.py's edge_ignore_aa). See measure_all_ew()
        #for how the winning order is chosen.
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
        if keep_result:
            self.lines_found_position[i] = found_line

        #in_line/other_than_line: boolean split of the fit window into the
        #line's own core (between its detected boundaries) and everything
        #else. other_than_line keeps the original inclusive (<=/>=)
        #formula; in_line is built as its exact complement (~) rather than
        #independently with its own inclusive bounds on both sides -- two
        #independently-inclusive formulas double-cover whichever real data
        #point happens to sit exactly ON a boundary (line_bound itself IS
        #an actual wavelength grid value, so this isn't just a theoretical
        #edge case), silently pulling that point into "the line" for
        #in_line's purposes while ALSO still getting pinned to continuum
        #by other_than_line -- confirmed to shift fit results measurably.
        #(only_line used to be built with `|` instead of `&`, which is
        #true for virtually every point regardless of line_bound -- fixed,
        #since it's what restricts the chi-square check below to the
        #line itself instead of the whole window.)
        other_than_line = (measure_x_array <= line_bound[0]) | (measure_x_array >= line_bound[1])
        in_line = ~other_than_line
        only_line = np.where(in_line)
        #highlight points within errors of continuum (or 1.0)
        upper_cont_bounds = measure_y_array+ ex_params[0] + 2*temp_err_array/temp_pred_array
        lower_cont_bounds = measure_y_array+ ex_params[0] - 2*temp_err_array/temp_pred_array
        points_within_norm = np.where((norm > lower_cont_bounds)&(norm < upper_cont_bounds))

        #Direct weighted Gaussian fit to the real data -- no GP smoothing,
        #no Monte Carlo resampling. Invert the continuum-normalized flux
        #into a positive-going bump first, since gauss_model/gauss_ew
        #expect a positive amplitude for an absorption line.
        xtest = measure_x_array
        full_y = measure_y_array + ex_params[0]  # real, unflattened data
        y_fit = norm - full_y
        y_err = 2*temp_err_array/temp_pred_array

        wing_idx = np.where(other_than_line)[0]
        c0, c0_err = norm, 0.  # "no correction": continuum assumed flat at norm

        #ref_keep: independent cross-check against a high-S/N reference
        #atlas (see load_reference_atlas()/reference_atlas.py): catches a
        #real but shallow, gradual blend that a median/MAD clip can't tell
        #from ordinary noise. Computed UNCONDITIONALLY (not just under
        #fit_continuum) -- the always-on sloped diagnostic below (which
        #sets lines_internal_trend_blue/red, what check_for_flags()'s "no
        #flat continuum" check reads) needs it regardless of whether
        #fit_continuum's own local-corrected EW is being computed.
        ref_keep = None
        if self.ref_wave is not None:
            ref_keep = reference_continuum_mask(
                xtest[wing_idx], y_err[wing_idx], self.ref_wave, self.ref_flux,
                self.ref_resolving_power)

        #line half-width / Gaussian sigma seed -- also needed
        #unconditionally now, by both the fit_continuum branch below and
        #the always-on sloped diagnostic
        line_hwidth = max((line_bound[1]-line_bound[0])/2.0, 0.01)
        sigma_guess = line_hwidth/1.5

        if fit_continuum:
            #Estimate the local continuum (a flat bias, no slope) from the
            #real wing data -- outside this line's own detected boundary --
            #robustly excluding points that look like a different, deeper
            #feature, rather than assuming this window's continuum is
            #already exactly at norm. This is a separate ESTIMATION step,
            #not a parameter fit jointly with the line: letting continuum
            #and amplitude trade off in one fit, seeded from only the
            #line's own handful of core points, was confirmed to overfit
            #badly on weaker lines (see estimate_local_continuum()'s
            #docstring) -- a local continuum should be set by the many
            #nearby continuum points, not the line's own few. No slope
            #term: a fitted slope was confirmed to let a one-sided
            #contaminating neighbor (its wing entering only one side of
            #the window) tilt the whole local continuum -- a flat bias
            #can't be tilted that way (see estimate_local_continuum()'s
            #docstring). estimate_local_continuum() works in real-flux
            #space (it clips LOW outliers, i.e. deeper-absorption
            #contamination) -- c0 here is the real continuum level, not
            #yet the inverted-space offset the rest of this function uses.
            c0, c0_err, _ = estimate_local_continuum(
                xtest[wing_idx], full_y[wing_idx], y_err[wing_idx], ref_keep=ref_keep)

        #convert to the inverted-space offset (norm - real continuum) that
        #y_fit/gauss_model operate in throughout the rest of this function
        cont_offset = norm - c0
        y_detrend = y_fit - cont_offset

        if fit_continuum:
            #fit the line against the now continuum-corrected data: the
            #line's own core, plus any wing points that don't still look
            #like a separate deeper feature after detrending -- giving the
            #fit real leverage on where the (now properly zeroed) baseline
            #sits, against a genuinely corrected local continuum
            fit_mask = in_line.copy()
            fit_mask[wing_idx] = np.abs(y_detrend[wing_idx]) < 5*np.median(y_err[wing_idx])
            if ref_keep is not None:
                #without this, a reference-atlas-flagged point can still be
                #included here even though estimate_local_continuum() above
                #excluded it from c0 -- confirmed on a real case (Fe I
                #5577.03): a still-declining neighbor-wing tail too close to
                #c0 to trip the |y_detrend|<5*median(y_err) cut leaked into
                #the Gaussian fit, pulling its free `baseline` parameter
                #away from zero and making the plotted fit curve's far-wing
                #level (c0 - baseline) visibly diverge from c0 itself.
                #Same fallback spirit as estimate_local_continuum(): don't
                #let this restriction alone starve the fit of points.
                restricted_fit_mask = fit_mask.copy()
                restricted_fit_mask[wing_idx] &= ref_keep
                if restricted_fit_mask.sum() >= 5:
                    fit_mask = restricted_fit_mask
            bf, pcov, p0 = gfit_direct(xtest[fit_mask], y_detrend[fit_mask], y_err[fit_mask],
                                        found_line, sigma_guess, 0.)
        else:
            #legacy: fit over the WHOLE window, with wing points pinned to
            #exactly 0 rather than excluded -- confirmed to matter, not
            #just be equivalent-but-wasteful: pinning gives the fit strong
            #baseline=0 leverage from dozens of points, and simply
            #excluding them instead (as the fit_continuum branch does,
            #appropriately, once they're genuinely detrended) measurably
            #hurt convergence/stability here where they're NOT detrended
            fit_y = y_detrend.copy()
            fit_y[other_than_line] = 0.
            bf, pcov, p0 = gfit_direct(xtest, fit_y, y_err, found_line, 0.5, 0.)
        fail_bf = np.array([0., found_line, 0., 0.])
        if bf is None:
            print('Gaussian fit did not converge')
            print('If line is close to an edge, try remeasuring line with a smaller window size')
            bf = fail_bf

        ew = abs(gauss_ew(bf[0], bf[2]*2.355))
        #propagate BOTH the Gaussian fit's own covariance AND the local
        #continuum estimate's uncertainty (c0_err, ~0 when fit_continuum is
        #False) -- a correction drawn from a poorly-sampled/noisy wing
        #shouldn't be reported as confidently as one from a clean wing, even
        #though it's applied identically to the central EW value either way
        ew_err = np.sqrt(gauss_ew_err(bf[0], bf[2], pcov)**2 + (EW_K*bf[2]*c0_err)**2)

        #sanity bounds -- below 2 mA the line is too shallow/undetected to
        #trust, above 200 mA the fit likely locked onto the wrong (blended
        #or saturated) feature
        if bf[0] == 0 or not (2 < ew < 200):
            bf = fail_bf
            ew = 0.
            ew_err = 0.
            pcov = None  # don't shade a fit band for a rejected/failed fit

        #Automated quality-triggered retry: an outright failed fit, or a
        #>10% relative EW error in the LOCAL-continuum fit (checked here
        #since bf_global isn't computed until after this point; uses the
        #same 10% threshold check_for_flags() applies downstream to the
        #reported GLOBAL-continuum EW, just against the other fit) usually
        #means this window's wing didn't leave enough clean, uncontaminated
        #points to trust -- confirmed on Fe I
        #5579.335: a strong neighbor ~0.6 A away left only 6 points to
        #constrain a 4-parameter Gaussian, EW error 25 mA on a 10 mA line,
        #and simply widening the window to 2.5 A (bringing in real clean
        #continuum further out, since the neighbor's own core was already
        #fully inside the 1.5 A window, not truncated at its edge) fixed
        #it: error dropped to 0.65 mA with an unchanged central value.
        #auto_widen=False on the recursive call below is what stops this
        #at a single retry rather than an unbounded escalation.
        quality_failed = (ew == 0) or (ew_err >= 0.1*ew)
        if auto_widen and quality_failed and window_size < widen_window_size:
            print(f'line {self.lines[i]}: EW {ew:.2f}+/-{ew_err:.2f} at window_size='
                  f'{window_size} looks unreliable -- retrying with window_size='
                  f'{widen_window_size}')
            self.measure_ew(i, order, plot, ex_params, save_plot, widen_window_size,
                             show_plot, fit_continuum, auto_widen=False,
                             plot_window_size=plot_window_size, keep_result=keep_result)
            return

        #DIAGNOSTIC: a parallel conditionally-sloped local continuum --
        #feeds only lines_ew_local (the local-continuum-corrected
        #comparison value), never self.lines_ew itself (the reported EW,
        #from the GLOBAL-continuum fit below). See
        #estimate_local_continuum_sloped()'s docstring: exists to let
        #slope_sig_thresh/slope_min_points be tuned against real lines.
        #Computed UNCONDITIONALLY (not gated by fit_continuum) -- unlike
        #the flat local-corrected EW above, this block's OTHER output,
        #lines_internal_trend_blue/red (see just below), is what
        #check_for_flags()'s "no flat continuum" check reads, and that
        #flag needs to keep working even in the fit_continuum=False
        #(now-default) global-only path.
        bf_slope, pcov_slope, ew_slope, ew_err_slope, cont_offset_slope = None, None, None, None, None
        used_slope, c1_slope, slope_diag = False, 0., {}
        #record whether EACH side independently shows a significant
        #internal near-vs-far trend (estimate_local_continuum_sloped()'s
        #own condition 2, computed regardless of used_slope) -- both sides
        #trending at once means neither wing settles into a flat,
        #trustworthy stretch anywhere in the window, unlike the one-sided
        #case that condition 2 was originally designed to catch (a
        #recovering contamination tail on just one side). check_for_flags()
        #uses BOTH being True as a "no reliable local continuum found"
        #flag -- see its docstring.
        c0_s, c1_slope, c0_err_s, _, used_slope, slope_diag = estimate_local_continuum_sloped(
            xtest[wing_idx], full_y[wing_idx], y_err[wing_idx], found_line,
            min_points=slope_min_points, ref_keep=ref_keep, sig_thresh=slope_sig_thresh,
            internal_sig_thresh=slope_internal_sig_thresh)
        cont_offset_slope = norm - (c0_s + c1_slope*(xtest-found_line))
        y_detrend_slope = y_fit - cont_offset_slope
        fit_mask_slope = in_line.copy()
        fit_mask_slope[wing_idx] = np.abs(y_detrend_slope[wing_idx]) < 5*np.median(y_err[wing_idx])
        if ref_keep is not None:
            restricted_slope = fit_mask_slope.copy()
            restricted_slope[wing_idx] &= ref_keep
            if restricted_slope.sum() >= 5:
                fit_mask_slope = restricted_slope
        bf_slope, pcov_slope, _ = gfit_direct(
            xtest[fit_mask_slope], y_detrend_slope[fit_mask_slope], y_err[fit_mask_slope],
            found_line, sigma_guess, 0.)
        if bf_slope is None:
            bf_slope = fail_bf
        ew_slope = abs(gauss_ew(bf_slope[0], bf_slope[2]*2.355))
        ew_err_slope = np.sqrt(gauss_ew_err(bf_slope[0], bf_slope[2], pcov_slope)**2
                                + (EW_K*bf_slope[2]*c0_err_s)**2)
        if bf_slope[0] == 0 or not (2 < ew_slope < 200):
            bf_slope = fail_bf
            ew_slope = 0.
            ew_err_slope = 0.
            pcov_slope = None
        sig_str = 'n/a' if slope_diag.get('significance') is None else f"{slope_diag['significance']:.2f}"
        isig_b = slope_diag.get('internal_sig_blue')
        isig_r = slope_diag.get('internal_sig_red')
        isig_b_str = 'n/a' if isig_b is None else f"{isig_b:.2f}"
        isig_r_str = 'n/a' if isig_r is None else f"{isig_r:.2f}"
        if keep_result:
            self.lines_internal_trend_blue[i] = bool(slope_diag.get('internal_trend_blue'))
            self.lines_internal_trend_red[i] = bool(slope_diag.get('internal_trend_red'))
        print(f'  [slope diagnostic] EW(flat)={ew:.2f}+/-{ew_err:.2f}  '
              f'EW(sloped)={ew_slope:.2f}+/-{ew_err_slope:.2f}  used_slope={used_slope}  '
              f'c1={c1_slope:.5f}  n_blue={slope_diag.get("n_blue")} n_red={slope_diag.get("n_red")}  '
              f'significance={sig_str}  internal_sig(blue/red)={isig_b_str}/{isig_r_str}')

        best_bf = bf

        #Comparison fit: assumes the GLOBAL continuum normalization is
        #already exact (continuum pinned at norm, no local wing-based
        #correction) -- i.e. the same one-window direct fit as the
        #fit_continuum=False code path above. When fit_continuum=False was
        #actually requested, `bf`/`ew` above ARE this fit already, so just
        #reuse them instead of refitting.
        if fit_continuum:
            global_line_hwidth = max((line_bound[1]-line_bound[0])/2.0, 0.01)
            global_sigma_guess = global_line_hwidth/1.5
            fit_y_global = y_fit.copy()
            fit_y_global[other_than_line] = 0.
            bf_global, pcov_global, _ = gfit_direct(xtest, fit_y_global, y_err, found_line,
                                                     global_sigma_guess, 0.)
            if bf_global is None:
                bf_global = fail_bf
            ew_global = abs(gauss_ew(bf_global[0], bf_global[2]*2.355))
            ew_err_global = gauss_ew_err(bf_global[0], bf_global[2], pcov_global)
            if bf_global[0] == 0 or not (2 < ew_global < 200):
                bf_global = fail_bf
                ew_global = 0.
                ew_err_global = 0.
                pcov_global = None
        else:
            bf_global, pcov_global, ew_global, ew_err_global = best_bf, pcov, ew, ew_err

        #lines_gauss_Xsquare (check_for_flags()'s fit-quality check) is
        #deliberately evaluated against the GLOBAL-continuum fit
        #(bf_global, offset 0), NOT the local-continuum-corrected one
        #(best_bf/cont_offset) -- the local correction can make an
        #otherwise-biased fit (e.g. amplitude pulled by nearby blending
        #leaking past the wing exclusion) look artificially clean by
        #construction, since it's fit on the same detrended data this
        #residual would be measured against. Checking against the raw,
        #uncorrected global assumption instead asks a more basic
        #question: does a single clean Gaussian, sitting on the
        #spectrum's plain normalization, actually match this line's
        #observed core shape at all.
        fit_gauss_global = norm - (gauss_model(xtest, *bf_global) + 0.)
        diff = (fit_gauss_global[only_line] - full_y[only_line])**2
        #DEFAULT REPORTED EW: the GLOBAL-continuum fit (bf_global/
        #ew_global), not the local-continuum-corrected one -- see the
        #comment above where bf_global is computed. lines_bf_params
        #follows the same choice so it stays consistent with lines_ew.
        #best_bf/ew/ew_err (the local-continuum-corrected fit) are kept
        #separately in lines_ew_local/lines_ew_err_local for comparison
        #only. All gated by keep_result (see its docstring above) -- when
        #this line was also measured in another, better-placed order
        #(measure_all_ew()'s order-overlap handling), this order's numbers
        #are still fit/printed/plotted below for inspection, just not
        #written into self.lines_*[i].
        if keep_result:
            self.lines_gauss_Xsquare[i] = np.sum(diff)
            self.lines_bf_params[i] = bf_global
            self.lines_ew[i] = ew_global
            self.lines_ew_err[i] = ew_err_global
            self.lines_ew_local[i] = ew
            self.lines_ew_err_local[i] = ew_err
        print('line to measure:', ELEMENTS[self.lines_exd[i][0]],self.lines[i], '- Line found:', found_line)
        print('EW:',np.round(ew_global,2),u"±",np.round(ew_err_global,2))

        #Plotting stuff
        if plot:
            #plot_window_size widens the DISPLAYED data/fit-curve range
            #beyond the actual fit window (xtest/measure_x_array above,
            #untouched) -- re-fetched fresh from this order's full arrays,
            #centered on found_line, purely for visual context
            if plot_window_size is not None and plot_window_size > window_size:
                plot_sel = np.abs(self.shifted_wavelength[order] - found_line) <= plot_window_size/2.0
                plot_x_array = self.shifted_wavelength[order][plot_sel]
                plot_y_array = self.normalized_flux[order][plot_sel]
                plot_err_array = self.obs_err[order][plot_sel]
                plot_pred_array = self.pred_all[order][plot_sel]
            else:
                plot_x_array = measure_x_array
                plot_y_array = measure_y_array
                plot_err_array = temp_err_array
                plot_pred_array = temp_pred_array
            plot_upper_cont_bounds = plot_y_array + ex_params[0] + 2*plot_err_array/plot_pred_array
            plot_lower_cont_bounds = plot_y_array + ex_params[0] - 2*plot_err_array/plot_pred_array
            plot_points_within_norm = np.where((norm > plot_lower_cont_bounds) & (norm < plot_upper_cont_bounds))

            #bf_global/pcov_global/ew_global/ew_err_global -- the GLOBAL-
            #continuum comparison fit -- are already computed above (also
            #now the basis for lines_gauss_Xsquare's fit-quality check, not
            #just this plot's left panel).

            fig = plt.figure(figsize=(12,5))
            title = ("Order: " + str(order) + " " + "(" + str(np.round(self.shifted_wavelength[order].min(),3))
                      + "-" + str(np.round(self.shifted_wavelength[order].max(),3)) + ")")
            #lines_check_flag/lines_flag_reasons only reflect the truth as
            #of the LAST check_for_flags() call -- for a plot generated
            #before that's been (re)run on this measurement, this is
            #whatever it was left at previously (default: unflagged/'')
            if self.lines_check_flag[i]:
                title += "\nFLAGGED: " + str(self.lines_flag_reasons[i])
                fig.suptitle(title, color='#e41a1c', fontsize=10)
            else:
                fig.suptitle(title)

            #plot the fit/band/local-continuum curves on a 5x denser
            #wavelength grid than the actual data -- xtest only has one
            #point per real pixel, which makes a narrow line's Gaussian
            #fit curve look faceted/low-resolution; this is purely
            #cosmetic (fitting and the chi-square check above still use
            #the real data grid, unchanged). Spans plot_x_array's (possibly
            #widened) range, not just xtest's -- so the fit curve/local-
            #continuum line visually extrapolate across the wider view too.
            xplot = np.linspace(plot_x_array[0], plot_x_array[-1], len(plot_x_array)*5)

            def _draw_window(ax):
                #shared data/window-markers drawing for both panels below --
                #only the overlaid fit curve differs between them
                ax.grid()
                ax.set_xlabel(r'$\rm Wavelength~(\AA)$', size = 14)
                ax.errorbar(plot_x_array,plot_y_array + ex_params[0],
                     yerr=2*plot_err_array/plot_pred_array,capsize=0,fmt='.', color = 'k', label = 'cont', zorder = 2)
                ax.scatter(plot_x_array[plot_points_within_norm],plot_y_array[plot_points_within_norm] + ex_params[0], s = 10, c='#4daf4a', zorder = 3, alpha = 0.8)
                ax.plot([self.lines[i],self.lines[i]],[norm,norm*0.95], '--', color = 'k', alpha = 0.75)
                ax.plot([found_line,found_line],[norm,norm*0.95], '-', color='k')
                ax.plot([line_bound[0],line_bound[0]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
                ax.plot([line_bound[1],line_bound[1]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
                ax.annotate(str(self.lines[i]), xy = [self.lines[i], norm*1.025])
                ax.plot([plot_x_array[0],plot_x_array[-1]],[norm,norm], '--', color = '#4daf4a', label = 'assumed continuum (norm)')

            #Left panel: fit assuming the global continuum normalization is
            #already exact (no local wing-based correction)
            fit_view = fig.add_subplot(121)
            _draw_window(fit_view)
            fit_view.set_ylabel('Normalized Flux', size = 14)
            fit_view.set_title(f'Global continuum (REPORTED) -- EW={ew_global:.2f}±{ew_err_global:.2f} mÅ', fontsize=10)
            fit_gauss_plot_global = norm - (gauss_model(xplot, *bf_global) + 0.)
            fit_view.plot(xplot, fit_gauss_plot_global, '--', color = '#377eb8', lw= 2, label = 'Gaussian fit')
            if pcov_global is not None:
                model_err_plot_global = gauss_model_err(xplot, bf_global, pcov_global)
                fit_view.fill_between(xplot, fit_gauss_plot_global-model_err_plot_global, fit_gauss_plot_global+model_err_plot_global,
                         color = '#377eb8', alpha = 0.25, zorder = 1, label = r'fit $\pm1\sigma$')
            fit_view.legend(loc='best', fontsize=8)

            #Right panel: fit against the per-line estimated LOCAL continuum
            #(see estimate_local_continuum()) -- this is the LOCAL-
            #continuum-corrected comparison fit (lines_ew_local), shown
            #only when fit_continuum=True was explicitly requested; it
            #never drives self.lines_ew itself (always the GLOBAL-
            #continuum fit, bf_global, regardless of fit_continuum -- see
            #where self.lines_ew is set below). When fit_continuum=False
            #(the default) no local estimate was made, so there's nothing
            #distinct to show here
            local_view = fig.add_subplot(122)
            _draw_window(local_view)
            if fit_continuum:
                local_view.set_title(f'Local continuum (diagnostic only) -- EW={ew:.2f}±{ew_err:.2f} mÅ', fontsize=10)
                fit_gauss_plot_local = norm - (gauss_model(xplot, *best_bf) + cont_offset)
                local_view.plot(xplot, fit_gauss_plot_local, '--', color = '#377eb8', lw= 2, label = 'Gaussian fit')
                if pcov is not None:
                    model_err_plot_local = gauss_model_err(xplot, best_bf, pcov)
                    local_view.fill_between(xplot, fit_gauss_plot_local-model_err_plot_local, fit_gauss_plot_local+model_err_plot_local,
                             color = '#377eb8', alpha = 0.25, zorder = 1, label = r'fit $\pm1\sigma$')
                #the estimated LOCAL continuum level (c0, flat/no slope), so
                #you can see directly how far the global normalization was
                #off here -- this is what fixes a fit biased by imperfect
                #normalization
                local_cont_plot = np.full_like(xplot, norm - cont_offset)
                local_view.plot(xplot, local_cont_plot, ':', color = '#ff7f00', lw = 2, label = 'estimated local continuum')
                local_view.legend(loc='best', fontsize=8)
            else:
                local_view.set_title('Local continuum not estimated (fit_continuum=False)', fontsize=10)
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
            if show_plot:
                plt.show()
            else:
                plt.close(fig)

            print('#-----------------------#')

        #print extra parameter stuff
        if ex_params == [0,0,0,0]:
            pass
        else:
            if keep_result:
                self.lines_exp[i] = np.array(ex_params)
            print('extra params:',ex_params)

    def measure_all_ew(self, exclude_lines= [], plot_lines=[], ex_params = {}, window_size = 1.5, save_all = False, fit_continuum = False, auto_widen = True, widen_window_size = 2.5, slope_sig_thresh = 3.0, slope_min_points = 5, slope_internal_sig_thresh = 3.0, plot_window_size = None):
        """
        Measure every loaded line's EW.

        save_all=True additionally saves a per-line fit-quality plot (data,
        best-fit Gaussian, +-1sigma shaded fit uncertainty) for EVERY line to
        line_plots/<element>_<wavelength>_<order>.pdf -- use this to review
        fit quality across a whole linelist. Those figures are built and
        saved without being shown interactively (so hundreds of lines don't
        block on hundreds of plot windows); list specific wavelengths in
        plot_lines as well if you also want those shown live as they're
        measured.

        fit_continuum=False (default) reports every EW from the GLOBAL-
        continuum fit (continuum normalization assumed exact) -- set True
        to ALSO compute a local-continuum-corrected comparison fit per
        line (lines_ew_local); still never changes the reported lines_ew
        itself. See measure_ew()'s docstring for the full rationale.

        auto_widen=True (default) automatically retries a line once at
        widen_window_size if it comes out of window_size looking
        unreliable -- see measure_ew()'s docstring for exactly what
        triggers a retry. Applies regardless of fit_continuum.

        slope_sig_thresh, slope_min_points: DIAGNOSTIC ONLY, passed through
        to measure_ew()'s parallel flat-vs-sloped local continuum -- see
        its docstring and estimate_local_continuum_sloped()'s. Does not
        change the reported EW.

        plot_window_size : PLOTTING ONLY, passed through to measure_ew()
            -- see its docstring. Never affects the fit, EW, or error.

        Order overlap: a line near two adjacent orders' shared boundary
        can fall inside BOTH orders' wavelength ranges, so this loop
        measures it once per order it appears in -- but only the
        measurement from the order where the line sits FURTHEST from
        that order's own edge (in Angstroms) is kept in self.lines_*[i]
        (via measure_ew()'s keep_result); every other order's measurement
        of the same line is still fit/printed/plotted (if save_all/
        plot_lines request it) but discarded, never overwriting the kept
        one -- and this holds regardless of which order this loop happens
        to reach first. lines_order[i]/lines_edge_distance[i] record which
        order won and by how much. An edge measurement is exactly the
        kind most exposed to real order-boundary artifacts (see
        continuum.py's edge_ignore_aa and its module docstring for a
        real confirmed case), so this is a real reliability choice, not
        just deduplication for its own sake.
        """
        if save_all:
            make_plots_folder()

        #running "best so far" edge-distance per line, for this call only
        #(not persisted -- self.lines_edge_distance below is the public,
        #persisted record of the WINNING order's distance)
        best_edge_dist = np.full(len(self.lines), -np.inf)

        for order in range(len(self.wavelength)):
            for i in range(len(self.lines)):
                if self.lines[i] in exclude_lines:
                    self.lines_ew[i] = 0.0
                    self.lines_gauss_Xsquare[i] = np.nan
                    self.lines_bf_params[i] = None
                    self.lines_ew_err[i] = np.nan
                    self.lines_exp[i] = [0,0,0,0]
                    #self.lines_check_flag[i] = False
                elif self.lines[i] >= self.shifted_wavelength[order][0] and self.lines[i] <= self.shifted_wavelength[order][-1]:
                    plot = False
                    exp = [0,0,0,0]
                    if self.lines[i] in plot_lines:
                        plot = True
                        if self.lines[i] in ex_params.keys():
                            exp = ex_params[self.lines[i]]
                    #edge_dist: how far (Angstroms) the line's REST
                    #wavelength sits from the nearer edge of this order --
                    #used (not the later-fitted center) so this decision
                    #is made once, up front, the same way regardless of
                    #which order the loop reaches first; a difference of
                    #up to ~0.1 A (typical RV-shift/line-position scatter)
                    #is negligible against orders that are tens of A wide
                    edge_dist = min(self.lines[i] - self.shifted_wavelength[order][0],
                                     self.shifted_wavelength[order][-1] - self.lines[i])
                    keep_result = edge_dist > best_edge_dist[i]
                    if keep_result:
                        best_edge_dist[i] = edge_dist
                        self.lines_order[i] = order
                        self.lines_edge_distance[i] = edge_dist
                    if save_all:
                        plot = True
                        self.measure_ew(i,order, plot, exp, True, window_size,
                                         show_plot=(self.lines[i] in plot_lines), fit_continuum=fit_continuum,
                                         auto_widen=auto_widen, widen_window_size=widen_window_size,
                                         slope_sig_thresh=slope_sig_thresh, slope_min_points=slope_min_points,
                                         slope_internal_sig_thresh=slope_internal_sig_thresh,
                                         plot_window_size=plot_window_size, keep_result=keep_result)
                    else:
                        self.measure_ew(i,order, plot, exp, False, window_size, fit_continuum=fit_continuum,
                                         auto_widen=auto_widen, widen_window_size=widen_window_size,
                                         slope_sig_thresh=slope_sig_thresh, slope_min_points=slope_min_points,
                                         slope_internal_sig_thresh=slope_internal_sig_thresh, keep_result=keep_result,
                                         plot_window_size=plot_window_size)
        #self.lines_bf_params = np.array(self.lines_bf_params)

    def measure_line_ew(self,line,ex_params=[0,0,0,0], save_line = False, save_plot = False, window_size = 1.5, fit_continuum = False, auto_widen = True, widen_window_size = 2.5, slope_sig_thresh = 3.0, slope_min_points = 5, slope_internal_sig_thresh = 3.0, plot_window_size = None):
        if save_plot:
            make_plots_folder()
        i = np.where(self.lines == line)[0][0]
        found = False
        for order in range(len(self.wavelength)):
            if self.lines[i] >= self.shifted_wavelength[order][0] and self.lines[i] <= self.shifted_wavelength[order][-1]:
                if not found:
                    self.lines_ew[i] = 0.0
                    self.lines_gauss_Xsquare[i] = np.nan
                    self.lines_bf_params[i] = None
                    self.lines_ew_err[i] = np.nan
                    #self.lines_check_flag[i] = False
                    self.measure_ew(i,order, True, ex_params, save_plot, window_size, fit_continuum=fit_continuum,
                                    auto_widen=auto_widen, widen_window_size=widen_window_size,
                                    slope_sig_thresh=slope_sig_thresh, slope_min_points=slope_min_points,
                                    slope_internal_sig_thresh=slope_internal_sig_thresh,
                                    plot_window_size=plot_window_size)
                    found = True
                    if save_line:
                        with open('line_'+str(line)+'.txt','w') as f:
                            for k in range(len(self.shifted_wavelength[order])):
                                f.write("{0:10f}\t{1:10f}\n".format(self.shifted_wavelength[order][k],self.normalized_flux[order][k]))
                            print('order', order, 'saved!')

    def check_for_flags(self):
        """
        Flag lines whose measurement looks untrustworthy: high EW error
        fraction, too shallow to trust, a poor Gaussian fit, both wings
        showing no flat continuum stretch at all (see the no-flat-
        continuum check below -- distinct from the EW-error check: this
        can trip even when the resulting EW error looks small, because a
        handful of mutually-consistent survivors after clipping can still
        give a deceptively tight formal error despite most of the wing
        having been pervasively contaminated), or -- the
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
        GRACES this way, invisible to every check above). If
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
            #no-flat-continuum check - BOTH wings independently show a
            #significant internal near-vs-far trend (see
            #estimate_local_continuum_sloped()'s condition 2 and the
            #comment where measure_ew() records these). A single side
            #tripping this already disables the sloped diagnostic, but
            #says nothing about whether the window has a trustworthy
            #continuum stretch anywhere at all -- a real, isolated
            #recovering tail on just one side still leaves the other side
            #flat and able to anchor a decent flat continuum. Both sides
            #tripping it at once means neither wing ever settles into a
            #flat, trustworthy stretch anywhere in the window -- confirmed
            #on a real case (HD_10383's Fe I 5054.643: min(internal_sig_
            #blue, internal_sig_red) of ~4-5 in both its order
            #measurements, clearly above the flagship clean case Fe I
            #5522.447's ~3.7, and well below confirmed-bad cases like
            #5587.574/5546.5/5234.625 at ~12-21) -- i.e. pervasive
            #weak-line blending throughout the window rather than one
            #identifiable contaminating neighbor. A window like that is
            #suspect for EITHER continuum assumption (local or global),
            #not just the local one -- this is a pure visibility flag: it
            #does not itself change lines_ew/lines_ew_err (the GLOBAL-
            #continuum fit, see measure_ew()).
            if self.lines_internal_trend_blue[i] and self.lines_internal_trend_red[i]:
                self.lines_check_flag[i] = True
                reasons.append('no flat continuum found in either wing (both sides show a '
                                'significant internal trend) -- likely pervasive weak-line blending')
                print(self.lines[i], 'both wings show a significant internal trend -- '
                      'no flat continuum stretch found anywhere in the window, local '
                      'continuum estimate may be unreliable')
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
            #a human reviewer's explicit "keep it anyway" (see pipeline.py)
            #overrides every check above -- reasons/lines_flag_reasons stay
            #populated either way, so make_ew_doc() can still note WHY this
            #line needed a human decision even though it ends up kept
            if self.lines_human_keep[i]:
                self.lines_check_flag[i] = False

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
