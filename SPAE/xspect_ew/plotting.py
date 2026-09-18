"""Diagnostic plots and the line_plots/ output folder helper."""

import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt


def plot_ew_fit(order, wave_min, wave_max, line_rest, found_line, line_bound,
                 measure_x_array, measure_y_array, temp_err_array, temp_pred_array,
                 points_within_norm, xtest, m_plot, C, fit_gauss, ex_params,
                 norm=1.0, axes=None):
    """
    Draw one line's EW-fit window -- the exact plot Spectrum_Data.measure_ew()
    shows when plot=True, factored out into its own function so it can also
    be drawn into EXISTING axes (interactive.EWWidget's live EW-measurement
    stage) instead of always creating a new figure. Pure drawing, no
    fitting -- every argument is already computed by measure_ew() itself
    (or, for the widget, whatever the current EW-measurement routine
    computes; this function only needs the same array shapes measure_ew()
    already produces, not any of its internals).

    Parameters
    ----------
    order : int, for the title.
    wave_min, wave_max : this order's wavelength range, for the title.
    line_rest : the linelist rest wavelength.
    found_line, line_bound : as returned by get_line_window().
    measure_x_array, measure_y_array, temp_err_array, temp_pred_array :
        the window's wavelength/normalized-flux/error/continuum arrays.
    points_within_norm : index array, points consistent with continuum.
    xtest, m_plot, C : the GP fit's test grid, mean, and covariance.
    fit_gauss : the Gaussian fit evaluated on xtest.
    ex_params : [continuum_shift, left_bound, right_bound, center] -- the
        same manual-adjustment parameters measure_ew() accepts.
    norm : the continuum level (always 1.0 for normalized flux).
    axes : (fit_view, data_view) existing Axes to draw into (cleared
        first), or None to create a new figure (measure_ew()'s original
        behavior).

    Returns
    -------
    fig, (fit_view, data_view)
    """
    title = f"Order: {order} ({wave_min:.3f}-{wave_max:.3f})"
    if axes is None:
        fig = plt.figure(figsize=(12, 5))
        fig.suptitle(title)
        fit_view = fig.add_subplot(121)
        data_view = fig.add_subplot(122)
    else:
        fit_view, data_view = axes
        fig = fit_view.get_figure()
        fit_view.clear()
        data_view.clear()
        fit_view.set_title(title, fontsize=10)

    fit_view.grid()
    fit_view.set_xlabel(r'$\rm Wavelength~(\AA)$', size=14)
    fit_view.set_ylabel('Normalized Flux', size=14)
    fit_view.errorbar(measure_x_array, measure_y_array + ex_params[0],
                       yerr=2 * temp_err_array / temp_pred_array, capsize=0, fmt='.',
                       color='k', label='cont', zorder=2)
    fit_view.scatter(measure_x_array[points_within_norm],
                      measure_y_array[points_within_norm] + ex_params[0],
                      s=10, c='#4daf4a', zorder=3, alpha=0.8)
    fit_view.fill_between(xtest, m_plot + 2 * np.sqrt(np.diag(C)),
                           m_plot - 2 * np.sqrt(np.diag(C)), color='#999999', alpha=0.5)
    fit_view.plot([line_rest, line_rest], [norm, norm * 0.95], '--', color='k', alpha=0.75)
    fit_view.plot([found_line, found_line], [norm, norm * 0.95], '-', color='k')
    fit_view.plot([line_bound[0], line_bound[0]], [norm * 1.025, norm * 0.95],
                  '--', color='#e41a1c', alpha=0.5)
    fit_view.plot([line_bound[1], line_bound[1]], [norm * 1.025, norm * 0.95],
                  '--', color='#e41a1c', alpha=0.5)
    fit_view.annotate(str(line_rest), xy=[line_rest, norm * 1.025])
    fit_view.plot(xtest, fit_gauss, '--', color='#377eb8', lw=2)
    fit_view.plot([xtest[0], xtest[-1]], [norm, norm], '--', color='#4daf4a')

    data_view.grid()
    data_view.set_xlabel(r'$\rm Wavelength~(\AA)$', size=14)
    data_view.scatter(measure_x_array, measure_y_array + ex_params[0], s=5, c='k', zorder=2)
    data_view.errorbar(measure_x_array, measure_y_array + ex_params[0],
                        yerr=2 * temp_err_array / temp_pred_array, capsize=0,
                        fmt='.', color='k', zorder=3, alpha=0.5)

    fig.tight_layout()
    return fig, (fit_view, data_view)


def plot_line_info(star, name, filt = None):
    fig = plt.figure(figsize = (10,5))
    info_plot = fig.add_subplot(111)
    if filt == None:
        info_plot.errorbar(star.lines, star.lines_ew, yerr = star.lines_ew_err, fmt='.', zorder = 2, ecolor='k', c = 'k')
        info_plot.scatter(star.lines, star.lines_ew, c= star.lines_gauss_Xsquare, cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        info_plot.scatter(star.lines[star.lines_check_flag], star.lines_ew[star.lines_check_flag], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)
    else:
        info_plot.errorbar(star.lines[filt], star.lines_ew[filt], yerr = star.lines_ew_err[filt], fmt='.', zorder = 2, ecolor='k', c = 'k')
        info_plot.scatter(star.lines[filt], star.lines_ew[filt], c= star.lines_gauss_Xsquare[filt], cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        info_plot.scatter(star.lines[filt][star.lines_check_flag[filt]], star.lines_ew[filt][star.lines_check_flag[filt]], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)

    info_plot.grid()
    info_plot.set_title(name, size = 20)
    info_plot.set_xlabel(r'$\rm Wavelength\ (nm)$', size = 15)
    info_plot.set_ylabel(r'$\rm Equivalent\ Width\ (mA)$', size = 15)
    #plt.savefig(name+'_line_ew_info.pdf')
    plt.show()


def plot_comparison_res(star,hand_measured, name,xy = [0,100], filt = None):
    fig = plt.figure(figsize = (10,5))
    #top plot
    info_plot = fig.add_subplot(211)
    if filt == None:
        info_plot.errorbar(hand_measured, star.lines_ew, yerr = star.lines_ew_err, fmt='.', zorder = 2,ecolor='k', c = 'k')
        info_plot.scatter(hand_measured, star.lines_ew, c= star.lines_gauss_Xsquare, cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        info_plot.scatter(hand_measured[star.lines_check_flag], star.lines_ew[star.lines_check_flag], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)
    else:
        info_plot.errorbar(hand_measured[filt], star.lines_ew[filt], yerr = star.lines_ew_err[filt], fmt='.', zorder = 2,ecolor='k', c = 'k')
        info_plot.scatter(hand_measured[filt], star.lines_ew[filt], c= star.lines_gauss_Xsquare[filt], cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        info_plot.scatter(hand_measured[filt][star.lines_check_flag[filt]], star.lines_ew[filt][star.lines_check_flag[filt]], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)

    info_plot.plot([xy[0],xy[1]],[xy[0],xy[1]], 'k--')
    info_plot.set_title(name, size = 20)
    info_plot.grid()
    info_plot.set_ylabel(r'$\rm Auto\ Measured\ (mA)$', size = 15)

    #residuals plot
    res_plot = fig.add_subplot(212, sharex=info_plot)
    if filt == None:
        star_res_values = star.lines_ew - hand_measured
        res_plot.errorbar(hand_measured, star_res_values, yerr = star.lines_ew_err, fmt='.', zorder = 2,ecolor='k', c = 'k')
        res_plot.scatter(hand_measured, star_res_values, c= star.lines_gauss_Xsquare, cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        res_plot.scatter(hand_measured[star.lines_check_flag], star_res_values[star.lines_check_flag], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)
    else:
        star_res_values = star.lines_ew - hand_measured
        res_plot.errorbar(hand_measured[filt], star_res_values[filt], yerr = star.lines_ew_err[filt], fmt='.', zorder = 2,ecolor='k', c = 'k')
        res_plot.scatter(hand_measured[filt], star_res_values[filt], c= star.lines_gauss_Xsquare[filt], cmap = plt.cm.Reds, edgecolors='k', zorder = 3, s = 30)
        res_plot.scatter(hand_measured[filt][star.lines_check_flag[filt]], star_res_values[star.lines_check_flag], c = 'w', edgecolors= 'r',  zorder = 0, s = 100)
    #plt.savefig(name+'_ew_comparison.pdf')
    res_plot.plot([xy[0],xy[1]],[0,0],'k--')
    res_plot.set_xlabel(r'$\rm Hand\ Measured\ (mA)$', size = 15)
    res_plot.grid()
    plt.tight_layout()
    plt.show()


def make_plots_folder():
    folder_name = 'line_plots'
    #check if folder exists
    filenames = glob.glob('*')
    if folder_name in filenames:
        pass
    else:
        #make folder if not
        cmd = 'mkdir '+folder_name
        subprocess.call(cmd, shell=True)
