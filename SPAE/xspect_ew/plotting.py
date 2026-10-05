"""Diagnostic plots and the line_plots/ output folder helper."""

import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from .line_profile import gauss_model, gauss_model_err


def plot_ew_fit(order, wave_min, wave_max, line_rest, found_line, line_bound,
                 plot_x_array, plot_y_array, plot_err_array, plot_pred_array,
                 plot_points_within_norm, ex_params, norm,
                 bf_global, pcov_global, ew_global, ew_err_global,
                 fit_continuum, best_bf, pcov, cont_offset, ew, ew_err,
                 flagged=False, flag_reasons='', axes=None):
    """
    Draw one line's EW-fit window: two panels -- "Global continuum
    (REPORTED)" (always shown, the fit that actually sets self.lines_ew)
    and "Local continuum (diagnostic only)" (shown only when
    fit_continuum=True was requested) -- factored out of
    Spectrum_Data.measure_ew() so it can also be drawn into EXISTING axes
    (interactive.EWWidget's live EW-measurement stage) instead of always
    creating a new figure. Pure drawing, no fitting -- every argument is
    already computed by measure_ew() itself.

    Parameters
    ----------
    order : int, for the title.
    wave_min, wave_max : this order's wavelength range, for the title.
    line_rest : the linelist rest wavelength.
    found_line, line_bound : as returned by get_line_window().
    plot_x_array, plot_y_array, plot_err_array, plot_pred_array :
        the (possibly plot_window_size-widened) wavelength/normalized-
        flux/error/continuum arrays to display.
    plot_points_within_norm : index array, points consistent with continuum.
    ex_params : [continuum_shift, left_bound, right_bound, center] -- the
        same manual-adjustment parameters measure_ew() accepts.
    norm : the continuum level (always 1.0 for normalized flux).
    bf_global, pcov_global, ew_global, ew_err_global : the GLOBAL-continuum
        fit (left panel) -- always shown, this is what self.lines_ew
        reports regardless of fit_continuum.
    fit_continuum : whether the LOCAL-continuum comparison fit (right
        panel) was actually computed -- if False, that panel is left
        as a placeholder.
    best_bf, pcov, cont_offset, ew, ew_err : the LOCAL-continuum-corrected
        comparison fit (right panel), only meaningful when fit_continuum.
    flagged, flag_reasons : self.lines_check_flag[i]/lines_flag_reasons[i]
        -- colors/labels the title when this line was flagged.
    axes : (fit_view, local_view) existing Axes to draw into (cleared
        first), or None to create a new figure (measure_ew()'s original
        plot=True behavior).

    Returns
    -------
    fig, (fit_view, local_view)
    """
    xplot = np.linspace(plot_x_array[0], plot_x_array[-1], len(plot_x_array) * 5)
    title = f"Order: {order} ({wave_min:.3f}-{wave_max:.3f})"
    if flagged:
        title += "\nFLAGGED: " + str(flag_reasons)

    if axes is None:
        fig = plt.figure(figsize=(12, 5))
        if flagged:
            fig.suptitle(title, color='#e41a1c', fontsize=10)
        else:
            fig.suptitle(title)
        fit_view = fig.add_subplot(121)
        local_view = fig.add_subplot(122)
    else:
        fit_view, local_view = axes
        fig = fit_view.get_figure()
        fit_view.clear()
        local_view.clear()

    def _draw_window(ax):
        #shared data/window-markers drawing for both panels below -- only
        #the overlaid fit curve differs between them
        ax.grid()
        ax.set_xlabel(r'$\rm Wavelength~(\AA)$', size=14)
        ax.errorbar(plot_x_array, plot_y_array + ex_params[0],
                    yerr=2 * plot_err_array / plot_pred_array, capsize=0, fmt='.',
                    color='k', label='cont', zorder=2)
        ax.scatter(plot_x_array[plot_points_within_norm],
                   plot_y_array[plot_points_within_norm] + ex_params[0],
                   s=10, c='#4daf4a', zorder=3, alpha=0.8)
        ax.plot([line_rest, line_rest], [norm, norm * 0.95], '--', color='k', alpha=0.75)
        ax.plot([found_line, found_line], [norm, norm * 0.95], '-', color='k')
        ax.plot([line_bound[0], line_bound[0]], [norm * 1.025, norm * 0.95],
                '--', color='#e41a1c', alpha=0.5)
        ax.plot([line_bound[1], line_bound[1]], [norm * 1.025, norm * 0.95],
                '--', color='#e41a1c', alpha=0.5)
        ax.annotate(str(line_rest), xy=[line_rest, norm * 1.025])
        ax.plot([plot_x_array[0], plot_x_array[-1]], [norm, norm], '--', color='#4daf4a',
                label='assumed continuum (norm)')

    #Left panel: fit assuming the global continuum normalization is
    #already exact (no local wing-based correction) -- always shown, this
    #is the REPORTED fit (self.lines_ew).
    _draw_window(fit_view)
    fit_view.set_ylabel('Normalized Flux', size=14)
    fit_title = f'Global continuum (REPORTED) -- EW={ew_global:.2f}±{ew_err_global:.2f} mÅ'
    if axes is not None:
        # no fig.suptitle available when embedded in a caller's own figure
        # (e.g. the widget) -- fold the order/flag title into this panel.
        fit_title = title.replace('\n', '  ') + '\n' + fit_title
    title_kwargs = {'fontsize': 9 if axes is not None else 10}
    if flagged and axes is not None:
        title_kwargs['color'] = '#e41a1c'
    fit_view.set_title(fit_title, **title_kwargs)
    fit_gauss_plot_global = norm - (gauss_model(xplot, *bf_global) + 0.)
    fit_view.plot(xplot, fit_gauss_plot_global, '--', color='#377eb8', lw=2, label='Gaussian fit')
    if pcov_global is not None:
        model_err_plot_global = gauss_model_err(xplot, bf_global, pcov_global)
        fit_view.fill_between(xplot, fit_gauss_plot_global - model_err_plot_global,
                               fit_gauss_plot_global + model_err_plot_global,
                               color='#377eb8', alpha=0.25, zorder=1, label=r'fit $\pm1\sigma$')
    fit_view.legend(loc='best', fontsize=8)

    #Right panel: fit against the per-line estimated LOCAL continuum --
    #diagnostic/comparison only (lines_ew_local), shown only when
    #fit_continuum=True was requested; never drives self.lines_ew itself.
    _draw_window(local_view)
    if fit_continuum:
        local_view.set_title(f'Local continuum (diagnostic only) -- EW={ew:.2f}±{ew_err:.2f} mÅ',
                              fontsize=9 if axes is not None else 10)
        fit_gauss_plot_local = norm - (gauss_model(xplot, *best_bf) + cont_offset)
        local_view.plot(xplot, fit_gauss_plot_local, '--', color='#377eb8', lw=2, label='Gaussian fit')
        if pcov is not None:
            model_err_plot_local = gauss_model_err(xplot, best_bf, pcov)
            local_view.fill_between(xplot, fit_gauss_plot_local - model_err_plot_local,
                                     fit_gauss_plot_local + model_err_plot_local,
                                     color='#377eb8', alpha=0.25, zorder=1, label=r'fit $\pm1\sigma$')
        #the estimated LOCAL continuum level (c0, flat/no slope), so you
        #can see directly how far the global normalization was off here
        local_cont_plot = np.full_like(xplot, norm - cont_offset)
        local_view.plot(xplot, local_cont_plot, ':', color='#ff7f00', lw=2,
                         label='estimated local continuum')
        local_view.legend(loc='best', fontsize=8)
    else:
        local_view.set_title('Local continuum not estimated (fit_continuum=False)',
                              fontsize=9 if axes is not None else 10)

    fig.tight_layout()
    return fig, (fit_view, local_view)


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
