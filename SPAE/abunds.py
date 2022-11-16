#Functions to derive absolute abundances with MOOGSILENT and abundances relative to the Sun ([x/H])

from scipy.stats import linregress
import numpy as np
import os

from .read_write import read_file
from . import atmos as atmos
from . import prior


#Function to derive abundances
def abunds_func(x, print_atmosphere=True, print_moog=False):
    teff, logg, feh, micro = x

    if not in_bounds(x):
        return -np.inf, -np.inf

    #will want to have path to MOOGSILENT to be a user input
    output = atmos.atmos(teff, logg, feh)
    if print_atmosphere:
        atmos.print_output(output, "star.mod", teff, logg, feh, micro)

    # Call moog
    os.system('/usr/local/moognov2019silent/MOOGSILENT') #helium
    el_found, abundances = read_file("moog_out.2")

    return el_found, abundances



#Function for deriving relative abundances
def rel_abunds(el_found,abundances,sun_el,sun_abs,el):
    rel_abunds = np.array([],dtype = abundances[0].dtype)

    for i in range(len(el_found)):
        if el_found[i] == el:
            break

    for j in range(len(sun_el)):
        if sun_el[j] == el:
            break

    star_abunds = abundances[i]
    sun_abunds = sun_abs[j]


    for k, line in enumerate(star_abunds):
        if line['wavelength'] - sun_abunds[k]['wavelength'] == 0:
            rel_abunds = np.append(rel_abunds, line)
            rel_abunds[-1]['abund'] -= sun_abunds[k]['abund']
        else:
            print('Stellar and solar linelists do not match at wavelength ' + str(sun_abunds[k]['wavelength']) + '!')


    return rel_abunds


#Function for deriving absolute abundances
def abs_abunds(el_found,abundances,el):

    for i in range(len(el_found)):
        if el_found[i] == el:
            break

    star_abunds = abundances[i]

    return star_abunds


def in_bounds(x):
    teff, logg, feh, micro = x

    if teff < 3500 or teff > 7000:
        return False
    if logg < 1 or logg > 5.0:
        return False
    if feh < -4.0 or feh >= 0.5:
        return False
    if micro < 0 or micro > 3:
        return False

    return True


def obj_func(x, n_elems, sun_el=None, sun_abs=None, include_prior=False):
    """Define the objective function."""
    teff, logg, feh, micro = x

    # Check if stellar parameters are within valid ranges
    if not in_bounds(x):
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))

    # Calculate the prior
    if include_prior:
        p = teff, logg
        ln_prior = prior.ln_prior(p)
    else:
        ln_prior = 0


    el_found, abundances = abunds_func(x)
    if len(abundances) < 2:
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))

    if sun_el is None or sun_abs is None:
        abunds_fe1 = abs_abunds(el_found, abundances, 'Fe I ')
        abunds_fe2 = abs_abunds(el_found, abundances, 'Fe II ')
        fe_mean = feh + 7.50
    else:
        abunds_fe1 = rel_abunds(el_found,abundances,sun_el,sun_abs,'Fe I ')
        abunds_fe2 = rel_abunds(el_found,abundances,sun_el,sun_abs,'Fe II ')
        fe_mean = feh

    fe_std = np.std(np.append(abunds_fe1['abund'], abunds_fe2['abund']))

    ep_slope, ep_intercept, ep_r, ep_p, ep_stderr = linregress(abunds_fe1['EP'], abunds_fe1['abund'])
    rew_slope, rew_intercept, rew_r, rew_p, rew_stderr = linregress(abunds_fe1['logRWin'], abunds_fe1['abund'])

    fe1_likely = np.sum(-(abunds_fe1['abund'] - fe_mean)**2 / (2*fe_std**2)) - np.log(fe_std) * len(abunds_fe1['abund'])
    fe2_likely = np.sum(-(abunds_fe2['abund'] - fe_mean)**2 / (2*fe_std**2)) - np.log(fe_std) * len(abunds_fe2['abund'])

    ln_likelihood = fe1_likely + fe2_likely

    ln_posterior = ln_prior + ln_likelihood

    # Calculate posterior probability
    params_obj = tuple((ln_posterior, ep_r, rew_r))

    # Cycle through other lines to get their mean an std
    for i, el in enumerate(el_found):
        if sun_el is None or sun_abs is None:
            abunds_element = abs_abunds(el_found,abundances,el)
        else:
            abunds_element = rel_abunds(el_found,abundances,sun_el,sun_abs,el)
        n_lines = len(abunds_element['abund'])
        if n_lines == 1:
            params_obj = params_obj + (np.mean(abunds_element['abund']), 0)
        else:
            params_obj = params_obj + (np.mean(abunds_element['abund']), np.std(abunds_element['abund'])/np.sqrt(n_lines - 1))

    return params_obj



#Function to derive line-by-line abundances for Sun; sun_linelist is path to Sun linelist
def sun_abs(sun_linelist, teff_sun=5777, logg_sun=4.44, feh_sun=0.00, micro_sun=1.38):
    param_file(sun_linelist)
    x_sun = teff_sun, logg_sun, feh_sun, micro_sun
    sun_el, sun_abs = abunds_func(x_sun)

    return sun_el, sun_abs
