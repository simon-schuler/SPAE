"""Combining multiple Spectrum_Data objects, and helpers used by wavelength-
shift cleaning (Spectrum_Data.clean_shift())."""

import numpy as np


def make_line(x,m,b):
    return m*x+b


def combine_files(empty_obj,objects = []):
    final_wavelength = []
    final_flux = []
    final_norm_flux = []
    final_shifted_wavelength = []
    final_estimated_shift = []
    final_continuum = []
    final_obs_err = []
    final_pred_all = []
    final_pred_var_all = []
    final_gain = []

    for j in objects:

        for i in range(len(j.flux)):
            final_norm_flux.append(j.normalized_flux[i])
            final_shifted_wavelength.append(j.shifted_wavelength[i])
            final_wavelength.append(j.wavelength[i])
            final_flux.append(j.flux[i])
            final_estimated_shift.append(j.estimated_shift[i])
            final_continuum.append(j.continuum[i])
            final_obs_err.append(j.obs_err[i])
            final_pred_all.append(j.pred_all[i])
            final_pred_var_all.append(j.pred_var_all[i])
            final_gain.append(j.gain[i])

    empty_obj.wavelength = np.array(final_wavelength)
    empty_obj.flux = np.array(final_flux)
    empty_obj.shifted_wavelength = np.array(final_shifted_wavelength)
    empty_obj.normalized_flux = np.array(final_norm_flux)
    empty_obj.estimated_shift = np.array(final_estimated_shift)
    empty_obj.continuum = np.array(final_continuum)
    empty_obj.obs_err = np.array(final_obs_err)
    empty_obj.pred_all = np.array(final_pred_all)
    empty_obj.pred_var_all = np.array(final_pred_var_all)
    empty_obj.gain = np.array(final_gain)
    del final_wavelength
    del final_flux
    del final_norm_flux
    del final_shifted_wavelength
    del final_estimated_shift
    del final_continuum
    del final_obs_err
    del final_pred_all
    del final_pred_var_all
    del final_gain

    return empty_obj


def reduce_cc(x,y,lines,lines_removed,limit=0.12):
    #check correlation before going further
    cc = np.corrcoef(x,y)
    print('starting cc', cc[0,1])

    if abs(cc[0,1]) < limit:
        print('cc good enough')
        return lines,x,y,lines_removed

    check_ccs = np.zeros(len(x))

    #remove largest cc difference
    for i in range(len(x)):
        new_x = np.delete(x,i)
        new_y = np.delete(y,i)
        new_cc = np.corrcoef(new_x,new_y)
        check_ccs[i] = new_cc[0,1]

    #Calculate differences
    diffs = [abs(cc[0,1])- abs(j) for j in check_ccs]
    #which gives largest difference?
    biggest_diff = np.where(diffs == max(diffs))[0][0]
    #remove that one line
    lines_removed.append([lines[biggest_diff],x[biggest_diff],y[biggest_diff]])
    x = np.delete(x,biggest_diff)
    y = np.delete(y,biggest_diff)
    lines = np.delete(lines,biggest_diff)
    print('line removed: ', lines_removed)

    #recalculate cc
    cc = np.corrcoef(x,y)
    print('ending cc', cc[0,1])

    #Call function again to remove lines until 0.12 is passed
    lines,x,y,lines_removed = reduce_cc(x,y,lines,lines_removed)

    return lines,x,y,lines_removed
