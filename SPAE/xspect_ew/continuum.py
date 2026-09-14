"""Continuum point selection and fitting for Spectrum_Data.normalize().

Continuum_scan's local-window percentile selection provides a reasonable
INITIAL guess, but was found (see conversation/session notes -- a ground-
truth synthetic test) to badly underestimate the continuum near any
moderately broad/strong line: a local ~1.5 A window sitting mostly or
entirely inside a wider line has no true-continuum-level points to select
at all, so it picks the "least bad" points available -- measured 20-70%
continuum underestimate near a synthetic 3-A-wide line, vs <1.1% error in
clean regions or near narrow lines.

A first attempt fixed this by refining that initial guess via GLOBAL
(order-wide) iterative sigma-clipped rejection against a Gaussian Process
fit. That did NOT work (verified against the same synthetic test: error
near the broad line was unchanged, -26%/-19%, essentially identical to
the original bug) -- and the reason is structural, not a tuning mistake:
a GP's local flexibility (its fitted length scale) is exactly what lets
it bend down and track into a contaminated region, and sigma-clip
rejection from below can't self-correct once the very first fit is
already biased low across a whole contiguous stretch -- points in the
line's wings look "consistent" with the already-too-low prediction and
never get flagged as outliers. Fighting that by pinning the length scale
just replaces one arbitrary knob with another.

iterative_continuum_select() instead fits a LOW-ORDER polynomial
(Chebyshev basis, for numerical stability over an order's wavelength
range), the traditional approach used by IRAF's continuum task and most
EW pipelines (ARES, DAOSPEC, etc.). Its rigidity is the point: a
degree-3-ish polynomial fit to a whole order cannot bend down to follow
a single several-Angstrom-wide line, no matter how deep, without ruining
the fit everywhere else in the order -- a structural guarantee against
the GP's failure mode, not a masking heuristic. It also has no
hyperparameter-optimization cost, so unlike the GP version there's no
need to freeze anything between rejection iterations -- refitting every
iteration is cheap (closed-form weighted least squares)."""

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import Chebyshev


class Continuum_scan():
    '''Selects points at the continuum
    '''
    def __init__(self, distx, depth):
        #values currently viewed for selection
        self.select_window = None
        #standard deviation of selected window
        #self.current_sig = None
        #Size of selection box in x axis
        self.distx = distx
        self.points_in_window = None
        #Input spectra
        self.data = None
        #Points selected as part of the continuum
        self.select_points = None
        #Relates to how deeply to move selection box into data
        self.depth = depth
        return None

    def load_data(self,x,y):
        self.data = np.array([x,y])
        self.select_points = np.zeros(len(x))
        self.points_in_window = len(self.data[0][np.where(self.data[0] <= self.data[0][0]+self.distx)])
        return None

    def scan(self):
        split_order_into = int(np.ceil(len(self.data[0])/self.points_in_window))
        split_order_x = np.array_split(self.data[0], split_order_into)
        split_order_y = np.array_split(self.data[1], split_order_into)
        for i in range(len(split_order_y)):
            dex = np.where((self.data[0] >= split_order_x[i][0])&(self.data[0] <= split_order_x[i][-1]))
            percent = np.percentile(split_order_y[i], self.depth)
            self.select_points[dex] = (split_order_y[i] >= percent)
        return None

    def view_selected(self):
        fig = plt.figure(figsize=(15,5))
        ax = fig.add_subplot(111)
        ax.scatter(self.data[0], self.data[1], c = '#cccccc', alpha = 0.75, s = 5)
        bool_points = (self.select_points == 1)
        ax.scatter(self.data[0][bool_points],self.data[1][bool_points], c = 'g', s = 5)
        plt.show()
        return None

    def get_selected(self):
        bool_points = (self.select_points == 1)
        return bool_points


def fit_poly_continuum(wave, flux, err, select, degree=3):
    """
    Fit a weighted low-order Chebyshev polynomial to flux[select] and
    evaluate it at EVERY point in wave (not just the selected ones).

    The fit domain is pinned to [wave.min(), wave.max()] (the full
    order), not derived from wave[select] -- keeps the polynomial's
    internal x-rescaling (and therefore its behavior evaluated outside
    the selected points) consistent across rejection iterations as
    `select` changes.

    Returns
    -------
    pred : polynomial evaluated at every point in wave
    pred_var : formal variance of the fit residuals, broadcast to every
        point (a single order-wide number, not a per-point predictive
        variance like the old GP version -- nothing downstream currently
        consumes more than that; see continuum.py module docstring)
    """
    domain = [wave.min(), wave.max()]
    poly = Chebyshev.fit(wave[select], flux[select], degree, domain=domain, w=1.0 / err[select])
    pred = poly(wave)
    resid = flux[select] - poly(wave[select])
    pred_var = np.full_like(wave, np.var(resid))
    return pred, pred_var


def iterative_continuum_select(wave, flux, err, initial_select, degree=3,
                                n_iterations=5, low_reject_sigma=2.5,
                                high_reject_sigma=5.0, min_points=10):
    """
    Refine a continuum-point selection via GLOBAL (order-wide) iterative
    asymmetric sigma-clipped rejection against a low-order polynomial
    fit. See this module's docstring for why a rigid low-order fit,
    rather than a locally-flexible GP, is what actually fixes
    Continuum_scan's local-window bias near broad/strong lines.

    Each iteration: fit the polynomial to the current selection
    (fit_poly_continuum -- cheap, so refitting every iteration is fine,
    unlike the old GP version), then reselect as "continuum" every point
    within (-low_reject_sigma*err, +high_reject_sigma*err) of the fit.
    Low-side rejection is tight (absorption pulls flux down -- that's
    the failure mode of concern); high-side is loose by default, mainly
    to guard against a cosmic ray or emission spike rather than ordinary
    noise scatter above the continuum. Stops early if the selection
    stops changing.

    A min_points safety guard prevents runaway rejection (e.g. from a
    badly-chosen low_reject_sigma) from leaving too few points for the
    polynomial to fit at all -- if a proposed new selection would drop
    below min_points, the previous (working) selection is kept instead
    and iteration stops.

    Returns
    -------
    select : final boolean selection
    pred, pred_var : final polynomial evaluation at every point in wave
        (see fit_poly_continuum for what pred_var means here)
    """
    select = initial_select.copy()
    pred, pred_var = fit_poly_continuum(wave, flux, err, select, degree=degree)
    for _ in range(n_iterations):
        resid = flux - pred
        new_select = (resid > -low_reject_sigma * err) & (resid < high_reject_sigma * err)
        if new_select.sum() < min_points:
            break
        if np.array_equal(new_select, select):
            break
        select = new_select
        pred, pred_var = fit_poly_continuum(wave, flux, err, select, degree=degree)
    return select, pred, pred_var
