import os
import numpy as np
import emcee
import time
from multiprocessing import Pool
from .abunds import abunds_func, obj_func, in_bounds
from .moog.inlines import parse_linelist
from . import atmos as _atmos


def _worker_init():
    """Pre-load the Kurucz interpolator in each Pool worker at startup."""
    _atmos.load_data()


def run_spae(linelist, sun_el=None, sun_abs=None, x_0=(5777, 4.44, 0.01, 1.38),
             n_dim=4, n_walkers=40, n_steps=1000, include_prior=False,
             n_cores=None):
    """Run the SPAE MCMC sampler.

    Parameters
    ----------
    linelist     : str   — absolute path to the MOOG line list for the star
    sun_el       : list  — element labels from sun_abs() for differential mode
    sun_abs      : list  — solar abundances from sun_abs()
    x_0          : tuple — initial guess (teff, logg, feh, micro)
    n_walkers    : int   — number of emcee walkers (default 40)
    n_steps      : int   — MCMC steps per walker (default 1000)
    include_prior: bool  — include spectroscopic prior on Teff/logg
    n_cores      : int or None — CPU cores for parallel likelihood evaluation;
                   None uses all available cores, 1 disables multiprocessing
    """
    t_0 = time.time()

    if n_cores is None:
        n_cores = os.cpu_count()

    # Parse linelist once; the dict is passed through the entire MCMC run
    # so no file I/O occurs inside the hot loop
    if isinstance(linelist, str):
        linelist = parse_linelist(linelist)

    x_ball = np.array(x_0) * (1 + 1.0e-2 * (-0.5 + np.random.rand(n_walkers, n_dim)))

    el_found, abundances = abunds_func(x_0, linelist)
    blobs_dtype = [('ep_r', 'f8'), ('rew_r', 'f8'), ('ep_slope', 'f8'), ('rew_slope', 'f8')]
    for el in el_found:
        blobs_dtype.append((el, 'f8'))
        blobs_dtype.append((el + '_sigma_mean', 'f8'))

    sampler_kwargs = dict(
        blobs_dtype=blobs_dtype,
        args=(len(el_found), linelist, sun_el, sun_abs, include_prior),
    )

    if n_cores > 1:
        with Pool(n_cores, initializer=_worker_init) as pool:
            sampler = emcee.EnsembleSampler(
                n_walkers, n_dim, obj_func, pool=pool, **sampler_kwargs)
            sampler.run_mcmc(x_ball, n_steps, progress=True)
    else:
        sampler = emcee.EnsembleSampler(
            n_walkers, n_dim, obj_func, **sampler_kwargs)
        sampler.run_mcmc(x_ball, n_steps, progress=True)

    flat_blob = sampler.get_blobs().reshape((n_walkers * n_steps))

    log  = f'Runtime: {time.time() - t_0:.1f} s\n'
    log += f'Cores used: {n_cores}\n'
    log += f'Size of Full Array: {len(flat_blob)}\n'
    log += 'Acceptance Fractions: ' + np.array_str(sampler.acceptance_fraction) + '\n'

    return sampler, flat_blob, log


def run_one(teff, logg, feh, micro, linelist):
    """Derive line-by-line abundances for a single set of stellar parameters.

    Args:
        teff:     effective temperature (3500–7000 K)
        logg:     log surface gravity (1–5)
        feh:      metallicity [Fe/H] (−4 to 0.5)
        micro:    microturbulence (0–3 km/s)
        linelist: path to the MOOG line list file
    """
    params = (teff, logg, feh, micro)

    if not in_bounds(params):
        print('Your parameters are not in the bounds of the model grids.')
        print()
        print(run_one.__doc__)
        return

    el_found, abundances = abunds_func(params, linelist)
    return el_found, abundances
