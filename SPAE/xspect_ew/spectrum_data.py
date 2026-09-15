"""Spectrum_Data: the core driver class for loading a spectrum, fitting its
continuum, wave-shifting against a reference, and measuring line EWs."""

import copy
import glob
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simpson

from .constants import ELEMENTS
from .continuum import Continuum_scan, iterative_continuum_select
from .line_profile import (get_line_window, gauss_model, gfit_direct, gauss_ew, gauss_ew_err,
                            gauss_model_err, estimate_local_continuum, EW_K)
from .combine import make_line, parabolic_refine
from .plotting import make_plots_folder
from .readers import read_spectrum
from .radial_velocity import measure_effective_rv, C_KMS
from .response_correction import apply_response_correction as _apply_response_correction


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
        #used to switch between Adamow ew calculation and simpson's rule integration
        self.temp_line_ew = None
        self.temp_line_ew_err = None

    def apply_response_correction(self, response_wave, response, min_overlap_fraction=0.5,
                                   min_response_fraction=0.1):
        """
        Divide out an instrument response/blaze correction curve. General:
        works with a response curve from any source, matched to this
        spectrum's own orders by wavelength overlap -- see
        response_correction.py's module docstring for why. Run this
        BEFORE normalize()/normalize_all(); it corrects the raw counts
        (self.flux), and normalize()'s continuum fit will be far more
        robust on an already-flattened spectrum.

        For MAROON-X specifically: response_wave/response can come from
        readers.load_maroonx_response('MAROON-X_PHOENIX_RESPONSE_...hd5')
        -- a separate calibration file, not embedded in individual science
        exposures (confirmed empty there).

        Returns
        -------
        corrected_orders : list of order indices that were actually
            corrected (others may have been skipped -- see
            response_correction.apply_response_correction()'s docstring,
            including its min_response_fraction edge-pixel guard).
        """
        return _apply_response_correction(self, response_wave, response,
                                           min_overlap_fraction=min_overlap_fraction,
                                           min_response_fraction=min_response_fraction)

    def normalize_all(self, window_width = 1.5, continuum_depth = 90, degree = 3,
                       n_iterations = 5, low_reject_sigma = 2.5, high_reject_sigma = 5.0):
        #loop through orders
        for i in range(len(self.flux)):

            #fit continuum with an iteratively sigma-clipped low-order polynomial
            self.normalize(i, window_width, continuum_depth, degree=degree,
                            n_iterations=n_iterations, low_reject_sigma=low_reject_sigma,
                            high_reject_sigma=high_reject_sigma)

            #Replace un-normalized points with value before it
            #This should only be replacing the last point in the
            #spectrum that is always missed by normalize
            # err_est = self.obs_err[i]/self.pred_all[i]
            # non_norm_points = np.where(self.normalized_flux[i] > np.average(self.normalized_flux[i][self.continuum[i]]+err_est[self.continuum[i]]*100))
            # #replace non norm points
            # self.normalized_flux[i][non_norm_points] = self.normalized_flux[i][non_norm_points[0]-1]

        return None

    def normalize(self, order, window_width = 1.5, continuum_depth = 90, clip = [-999,-999],
                  degree = 3, n_iterations = 5, low_reject_sigma = 2.5, high_reject_sigma = 5.0):
        """
        Fit and divide out the continuum for one order.

        Continuum_scan's local-window percentile selection (window_width,
        continuum_depth) is used only as an INITIAL guess; the final
        continuum-point selection and fit come from
        iterative_continuum_select() (continuum.py), which refines that
        guess via GLOBAL (order-wide), iteratively sigma-clipped low-order
        polynomial fit. This fixes a real failure mode of trusting the
        local-window selection as final: a window sitting entirely inside
        a moderately broad/strong line has no true-continuum points to
        select at all, and was measured (see conversation/session notes)
        to bias the fitted continuum low there by 20-70% -- silently
        making any such line (and anything sharing its local fit) look
        shallower than it really is, i.e. underestimating its EW. A
        first attempt at fixing this with a globally-refit Gaussian
        Process did NOT work (same bug, same magnitude) -- see
        continuum.py's module docstring for why a rigid low-order
        polynomial fit is structurally the right tool here instead.
        degree/n_iterations/low_reject_sigma/high_reject_sigma control
        the refinement; see iterative_continuum_select()'s docstring.
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
        # obs_err is maintained alongside self.flux -- set from raw Poisson
        # shot noise (sqrt(raw counts)) at construction, and correctly
        # rescaled by apply_response_correction() if that was applied.
        # Recomputing sqrt(flux) here instead would silently be wrong once
        # response correction has run, since flux is then raw/R and
        # sqrt(raw/R) != sigma(raw)/R -- see apply_response_correction()'s
        # docstring.
        err = self.obs_err[order][clipl:clipr]

        continuum_scan_obj = Continuum_scan(window_width, continuum_depth)
        continuum_scan_obj.load_data(wave, flux)
        continuum_scan_obj.scan()
        initial_select = continuum_scan_obj.get_selected()
        del continuum_scan_obj

        cont, pred, pred_var = iterative_continuum_select(
            wave, flux, err, initial_select, degree=degree,
            n_iterations=n_iterations, low_reject_sigma=low_reject_sigma,
            high_reject_sigma=high_reject_sigma)

        self.continuum[order][clipl:clipr] = cont
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

    def apply_rv_shift(self, rv=None, lines=None, min_depth=0.02, verbose=False):
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

        Parameters
        ----------
        rv : float, km/s, optional -- apply this RV directly and skip line
            measurement (e.g. if the RV is already known from elsewhere).
        lines : {name: (rest_wavelength, window)}, optional -- defaults to
            radial_velocity.RV_REFERENCE_LINES.
        min_depth : float -- minimum line depth (in normalized flux) to
            trust a line's fitted center.
        verbose : bool -- print the measured RV and which lines were used.

        Returns
        -------
        rv : float, km/s -- the RV actually applied.
        """
        if rv is None:
            measured_rv, rv_err, used = measure_effective_rv(self, lines=lines, min_depth=min_depth)
            if measured_rv is None:
                raise ValueError(
                    "Could not measure an effective RV -- none of the reference "
                    "lines were found/usable in this spectrum's wavelength "
                    "coverage. Pass rv= directly, or lines= with a custom set.")
            if verbose:
                print(f"Effective RV = {measured_rv:.3f} +/- {rv_err:.3f} km/s, "
                      f"from {len(used)} line(s):")
                for name, restw, order, v in used:
                    print(f"  {name} ({restw} A, order {order}): v={v:.3f} km/s")
            self.rv = (round(measured_rv, 3), round(rv_err, 3))
            rv = measured_rv
        else:
            self.rv = (rv, 0.0)

        for order in range(len(self.wavelength)):
            self.shifted_wavelength[order] = self.wavelength[order] * (1.0 + rv / C_KMS)
            # Angstrom-equivalent at the order's mean wavelength, kept for
            # reporting/consistency with estimate_shift()'s convention --
            # the actually-applied shift above is the correct multiplicative
            # one, not this per-order scalar approximation of it.
            self.estimated_shift[order] = self.wavelength[order].mean() * (rv / C_KMS)

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
        for i in range(len(self.lines)):
            self.lines_exd[i] = np.array([elmnt[i],ep[i],gf[i],rad[i]])

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

    def measure_ew(self, i, order, plot = False, ex_params = [0,0,0,0], save_plot = False, window_size = 1.5, show_plot = True, fit_continuum = True):
        #extra parameters [0] - shift continuum
        #                 [1] - left boundary in Angstroms
        #                 [2] - right boundary in Angstroms
        #                 [3] - line center in Angstroms
        #show_plot: set False to save/build the figure without blocking on
        #plt.show() -- used by measure_all_ew(save_all=True) so a QC plot
        #for every line doesn't pop up (and need closing) one at a time
        #fit_continuum: estimate a local linear continuum level+slope (c0,
        #c1) from this line's own wing data (see estimate_local_continuum)
        #rather than assuming the global normalization already put this
        #window's continuum exactly at norm -- fixes fits biased by
        #imperfect normalization. Set False to fall back to the old fixed-
        #continuum-at-norm behavior.
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
        #since it's what restricts the chi-square/Simpson checks below to
        #the line itself instead of the whole window.)
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
        c0, c1, c0_err = norm, 0., 0.  # "no correction": continuum assumed flat at norm
        if fit_continuum:
            #Estimate the local continuum (level+slope) from the real wing
            #data -- outside this line's own detected boundary -- robustly
            #excluding points that look like a different, deeper feature,
            #rather than assuming this window's continuum is already
            #exactly at norm. This is a separate ESTIMATION step, not a
            #parameter fit jointly with the line: letting continuum and
            #amplitude trade off in one fit, seeded from only the line's
            #own handful of core points, was confirmed to overfit badly on
            #weaker lines (see estimate_local_continuum()'s docstring) --
            #a local continuum should be set by the many nearby continuum
            #points, not the line's own few. estimate_local_continuum()
            #works in real-flux space (it clips LOW outliers, i.e.
            #deeper-absorption contamination) -- c0/c1 here are the real
            #continuum level/slope, not yet the inverted-space offset the
            #rest of this function uses.
            c0, c1, c0_err, _ = estimate_local_continuum(
                xtest[wing_idx], full_y[wing_idx], y_err[wing_idx], found_line)

        #convert to the inverted-space offset (norm - real continuum) that
        #y_fit/gauss_model operate in throughout the rest of this function
        cont_offset = norm - (c0 + c1*(xtest-found_line))
        y_detrend = y_fit - cont_offset

        if fit_continuum:
            #fit the line against the now continuum-corrected data: the
            #line's own core, plus any wing points that don't still look
            #like a separate deeper feature after detrending -- giving the
            #fit real leverage on where the (now properly zeroed) baseline
            #sits, against a genuinely corrected local continuum
            fit_mask = in_line.copy()
            fit_mask[wing_idx] = np.abs(y_detrend[wing_idx]) < 5*np.median(y_err[wing_idx])
            line_hwidth = max((line_bound[1]-line_bound[0])/2.0, 0.01)
            sigma_guess = line_hwidth/1.5
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
            fit_mask = in_line  # only used below to pick the Simpson integration domain
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

        best_bf = bf
        #predicted real flux = norm - (line dip + local continuum offset)
        fit_gauss = norm - (gauss_model(xtest, *best_bf) + cont_offset)
        #set values for line

        diff = (fit_gauss[only_line] - full_y[only_line])**2
        self.lines_gauss_Xsquare[i] = np.sum(diff)

        self.lines_bf_params[i] = best_bf
        self.lines_ew[i] = ew
        self.lines_ew_err[i] = ew_err
        if ew == 0:
            self.lines_ew_simp[i] = 0
            self.lines_ew_simp_err[i] = np.nan
        else:
            #Simpson's-rule integration of the continuum-corrected profile,
            #as a cross-check on the Gaussian EW above -- a point estimate,
            #not a resampled distribution, so it carries no error of its
            #own. Integrated over fit_mask (same points the fit itself
            #used), not just the line's own narrow core (only_line): that
            #core can be just a handful of points spanning well under the
            #Gaussian's full area for a line whose auto-detected boundary
            #undershoots its true width, which was confirmed to make this
            #cross-check read ~2x low on real lines in the bundled sample
            #even though the Gaussian fit itself was fine.
            self.lines_ew_simp[i] = simpson(y_detrend[fit_mask], xtest[fit_mask])*1000
            self.lines_ew_simp_err[i] = np.nan
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
            fit_view.plot([self.lines[i],self.lines[i]],[norm,norm*0.95], '--', color = 'k', alpha = 0.75)
            fit_view.plot([found_line,found_line],[norm,norm*0.95], '-', color='k')
            fit_view.plot([line_bound[0],line_bound[0]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            fit_view.plot([line_bound[1],line_bound[1]],[norm*1.025,norm*0.95], '--', color = '#e41a1c', alpha = 0.5)
            fit_view.annotate(str(self.lines[i]), xy = [self.lines[i], norm*1.025])
            fit_view.plot(xtest, fit_gauss, '--', color = '#377eb8', lw= 2, label = 'Gaussian fit')
            if pcov is not None:
                model_err = gauss_model_err(xtest, best_bf, pcov)
                fit_view.fill_between(xtest, fit_gauss-model_err, fit_gauss+model_err,
                         color = '#377eb8', alpha = 0.25, zorder = 1, label = r'fit $\pm1\sigma$')
            fit_view.plot([xtest[0],xtest[-1]],[norm,norm], '--', color = '#4daf4a', label = 'assumed continuum (norm)')
            if fit_continuum:
                #the estimated LOCAL continuum level (c0, c1), so you can
                #see directly how far the global normalization was off
                #here -- this is what fixes a fit biased by imperfect
                #normalization
                local_cont = norm - cont_offset
                fit_view.plot(xtest, local_cont, ':', color = '#ff7f00', lw = 2, label = 'estimated local continuum')
            fit_view.legend(loc='best', fontsize=8)

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
            if show_plot:
                plt.show()
            else:
                plt.close(fig)

            print('#-----------------------#')

        #print extra parameter stuff
        if ex_params == [0,0,0,0]:
            pass
        else:
            self.lines_exp[i] = np.array(ex_params)
            print('extra params:',ex_params)

    def measure_all_ew(self, exclude_lines= [], plot_lines=[], ex_params = {}, window_size = 1.5, save_all = False, fit_continuum = True):
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

        fit_continuum=True (default) corrects for imperfect global
        continuum normalization per-line -- see measure_ew()'s docstring.
        """
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
                        plot = True
                        self.measure_ew(i,order, plot, exp, True, window_size,
                                         show_plot=(self.lines[i] in plot_lines), fit_continuum=fit_continuum)
                    else:
                        self.measure_ew(i,order, plot, exp, False, window_size, fit_continuum=fit_continuum)
        #self.lines_bf_params = np.array(self.lines_bf_params)

    def measure_line_ew(self,line,ex_params=[0,0,0,0], save_line = False, save_plot = False, window_size = 1.5, fit_continuum = True):
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
                    self.measure_ew(i,order, True, ex_params, save_plot, window_size, fit_continuum=fit_continuum)
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
        spectrum's wavelength correction is not well constrained.
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
