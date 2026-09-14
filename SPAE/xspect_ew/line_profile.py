"""Line-window finding and Gaussian profile fitting for EW measurement."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def get_line_window(line, wave, flux, left_bound, right_bound,
                     line_input, window_size=1.5):
    boundaries = [0, 0]
    #if no line is specified, auto fine best line guess
    if line_input == 0.0:
        #find line tip
        left_look = np.where((wave <= line)&(wave >= line - 0.1))
        right_look = np.where((wave >= line)&(wave <= line + 0.1))
        #find min
        mins = [flux[left_look].min(),flux[right_look].min()]
        best_line_guess = wave[np.where(flux == np.min(mins))][0]
    else:
        best_line_guess = line_input

    #get_window around line
    window = np.where((wave >= best_line_guess-window_size/2.0)&(wave <= best_line_guess+window_size/2.0))

    #calc derivative
    dy = np.gradient(flux[window])
    dy_std = dy.std()

    #if no left or right bound given set using std
    auto_bound_l = False
    auto_bound_r = False
    if left_bound == 0:
        dy_l = dy_std/2.0
        auto_bound_l = True
    if right_bound == 0:
        dy_r = dy_std/2.0
        auto_bound_r = True

    #if no line boundaries specified auto find boundaries
    if auto_bound_l:
        left_look = np.where(wave[window] <= best_line_guess - 0.05)
        dy1_left = np.where((dy[left_look] < dy_l)&(dy[left_look] > (-1)*dy_l))
        if len(wave[window][left_look][dy1_left]) ==0:
            print('line ',line,' very close to edge or dy selection value too small')
            print('will attempt to remeasure, if not possible, add line to exclude lines list in .measure_all_ew() function')
            plt.clf()
            plt.plot(wave[window],flux[window])
            plt.plot([line,line],[0.95,1.0], 'k')
            plt.annotate(str(line), xy=[line,1.01])
            plt.plot([best_line_guess,best_line_guess],[0.95,1.0], 'k--')
            plt.annotate(str(best_line_guess), xy=[best_line_guess,1.01])
            plt.show()
            return 1,1,0,0

        else:
            boundaries[0] = wave[window][left_look][dy1_left][-1]
    else:
        boundaries[0] = left_bound
    if auto_bound_r:
        right_look = np.where(wave[window] >= best_line_guess + 0.05)
        dy1_right = np.where((dy[right_look] < dy_r)&(dy[right_look] > (-1)*dy_r))
        if len(wave[window][right_look][dy1_right]) ==0:
            print('line ',line,' very close to edge or dy selection value too small')
            print('will attempt to remeasure, if not possible, add line to exclude lines list in .measure_all_ew() function')
            plt.clf()
            plt.plot(wave[window],flux[window])
            plt.plot([line,line],[0.95,1.0], 'k')
            plt.annotate(str(line), xy=[line,1.01])
            plt.plot([best_line_guess,best_line_guess],[0.95,1.0], 'k--')
            plt.annotate(str(best_line_guess), xy=[best_line_guess,1.01])
            plt.show()
            return 0,0,1,1

        else:
            boundaries[1] = wave[window][right_look][dy1_right][0]
    else:
        boundaries[1] = right_bound

    return window,best_line_guess, boundaries,dy


def gauss_model(x,A,mu,sigma, baseline):
    return A*np.exp(-(x-mu)**2/2/sigma**2) + baseline


def gfit(wav,flux,wav_cen, fwhm):
        sigma = fwhm/2.355
        #limit window of search center +- 2*fwhm to exclude other emission lines
        gwave = np.where((wav >= wav_cen-30)&(wav <= wav_cen+30))

        #find better center to account for small doppler shift within same window of search
        bet_cen = wav[np.where(flux == flux[gwave].max())[0][0]]

        #Initial value for guass max value guess from max of curve
        guess = flux[np.where(flux == flux[gwave].max())[0][0]]

        #Set parameters for gauss curve fit
        p0 = [guess,bet_cen,sigma, 0.]
        bf,cov = curve_fit(gauss_model,wav[gwave],flux[gwave],p0)

        #plt.plot(wav[gwave], flux[gwave], 'r')
        return bf, np.sqrt(np.diag(cov)), p0


def gfit_simple(x_array, y_array, mu, sigma, baseline):
    A = y_array.max()
    p0 = [A, mu, sigma, baseline]
    try:
        bf, cov = curve_fit(gauss_model, x_array, y_array, p0)
        return bf, np.sqrt(np.diag(cov)), p0
    except:
        bf, cov = [0,0,0,0],None
        return bf, cov, p0


def gauss_ew(a, fwhm):
    if a == 0 or fwhm == 0:
        return 0
    else:
        return 500.*a*np.sqrt(np.pi/np.log(2))*fwhm #From Adamow pyMOOG ew measure
