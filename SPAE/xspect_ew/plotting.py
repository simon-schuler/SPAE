"""Diagnostic plots and the line_plots/ output folder helper."""

import glob
import subprocess
import matplotlib.pyplot as plt


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
