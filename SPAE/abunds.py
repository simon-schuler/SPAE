#Functions to derive absolute abundances with MOOGSILENT and abundances relateive to the Sun ([x/H])

import numpy as np
import os

#Function to derive abundances
def abunds_func(x):
    teff, logg, feh, micro = x
    if teff < 3500 or teff > 7000:
        return -np.inf, -np.inf
    if logg < 1 or logg > 5.0:
        return -np.inf, -np.inf
    if feh < -4.0 or feh >= 0.5:
        return -np.inf, -np.inf
    if micro < 0 or micro > 3:
        return -np.inf, -np.inf
    
    model = '/usr/local/mspawn/mspawn -wstar.mod -t%.3f' % teff + ' -g%.3f' % logg
    
    if feh >= 0:
        model = model + ' -p%.3f' % feh
    else:
        model = model + ' -m%.3f' % abs(feh)
    
    model = model + ' -v%.2f' % micro

    #will want to have path to MOOGSILENT to be a user input
    os.system(model)
    os.system('/usr/local/moognov2019silent/MOOGSILENT') #helium
    el_found, abundances = read_file("moog_out.2")
      
    return el_found, abundances



#Function for deriving relative abundances
def rel_abs(el_found,abundances,sun_el,sun_abs,el):
    rel_abunds = np.array([],dtype = abundances[0].dtype)
    
    for i in range(len(el_found)):
        if el_found[i] == el:
            break
    
    for j in range(len(sun_el)):
        if sun_el[j] == el:
            break
            
    star_abunds = abundances[i]
    sun_abunds = sun_abs[j]
    

    for line in star_abunds:
        wave_dif = np.abs(line['wavelength'] - sun_abunds['wavelength'])
        idx = np.argmin(wave_dif)
        
        if wave_dif[idx] < 0.2:
            rel_abunds = np.append(rel_abunds, line)
            rel_abunds[-1]['abund'] -= sun_abunds['abund'][idx]
            
    return rel_abunds



#Function defining the objective function
def obj_func(x,n_elems,sun_el,sun_abs):
    teff, logg, feh, micro = x
    
    if teff < 3500 or teff > 7000:
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))
    if logg < 1 or logg > 5.0:
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))
    if feh < -4.0 or feh >= 0.5:
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))
    if micro < 0 or micro > 3:
        return (-np.inf,) + tuple(np.zeros(2*n_elems+2))


    el_found, abundances = abunds_func(x)
    
    abunds_rel_fe1 = rel_abs(el_found,abundances,sun_el,sun_abs,'Fe I ')
    abunds_rel_fe2 = rel_abs(el_found,abundances,sun_el,sun_abs,'Fe II ')
    
    ep_slope, ep_intercept, ep_r, ep_p, ep_stderr = linregress(abunds_rel_fe1['EP'], abunds_rel_fe1['abund'])
    rew_slope, rew_intercept, rew_r, rew_p, rew_stderr = linregress(abunds_rel_fe1['logRWin'], abunds_rel_fe1['abund'])

    fe_std = np.std(np.append(abunds_rel_fe1['abund'], abunds_rel_fe2['abund']))
    
    fe1_likely = np.sum(-(abunds_rel_fe1['abund'] - feh)**2 / (2*fe_std**2)) - np.log(fe_std) * len(abunds_rel_fe1['abund'])
    fe2_likely = np.sum(-(abunds_rel_fe2['abund'] - feh)**2 / (2*fe_std**2)) - np.log(fe_std) * len(abunds_rel_fe2['abund'])
    
    params_obj = tuple((fe1_likely + fe2_likely, ep_r, rew_r))
    
    for i, el in enumerate(el_found):
        abunds_relative = rel_abs(el_found,abundances,sun_el,sun_abs,el)
        n_lines = len(abunds_relative['abund'])
        if n_lines == 1:
            params_obj = params_obj + (np.mean(abunds_relative['abund']), 0)
        else:
            params_obj = params_obj + (np.mean(abunds_relative['abund']), np.std(abunds_relative['abund'])/np.sqrt(n_lines - 1))
        

    return params_obj
    


#Function to derive line-by-line abundances for Sun; sun_linelist is path to Sun linelist
def sun_abs(sun_linelist, teff_sun=5777, logg_sun=4.44, feh_sun=0.00, micro_sun=1.38):
    param_file(sun_linelist)
    x_sun = teff_sun, logg_sun, feh_sun, micro_sun
    sun_el, sun_abs = abunds_func(x_sun)

    return sun_el, sun_abs

