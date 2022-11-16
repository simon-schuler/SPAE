import numpy as np
import emcee
import time
from .abunds import abunds_func, obj_func, in_bounds
from .read_write import param_file


def run_spae(sun_el=None,sun_abs=None,x_0=(5777,4.44,0.01,1.38),
		     n_dim=4, n_walkers=40, n_steps=1000, include_prior=False):

	t_0 = time.time()

	x_ball = x_0 * (1 + 1.0e-2 * (-0.5 + np.random.rand(n_walkers, n_dim)))

	el_found, abundances = abunds_func(x_0)
	blobs_dtype = [('ep_r','f8'), ('rew_r','f8')]

	for el in el_found:
	    blobs_dtype.append((el, 'f8'))
	    blobs_dtype.append((el+'_sigma_mean', 'f8'))


	sampler = emcee.EnsembleSampler(n_walkers, n_dim, obj_func, blobs_dtype=blobs_dtype, args=(len(el_found),sun_el,sun_abs,include_prior))
	results = sampler.run_mcmc(x_ball, n_steps, progress=True)

	flat_blob = sampler.get_blobs().reshape((n_walkers*n_steps))

	log = 'Runtime: ' + str(time.time() - t_0) + '\n'
	log = log + 'Size of Full Array: ' + str(len(flat_blob)) + '\n'
	log = log + 'Acceptance Fractions:' + np.array_str(sampler.acceptance_fraction) + '\n'

	return sampler, flat_blob, log

def run_one(teff, logg, feh, micro, linelist):
	"""Creates model atmosphere, MOOG parameter file, and runs abfind in MOOG.

	Creates model atmosphere, star.mod, by interpolating the Kurucz model atmosphere
	grids. Creates MOOG parameter file, batch.par. Outputs two files from MOOG, a
	summary of the abundance calculations and the derived abundances.

	Args:
		teff: effective temperature (3500-7000 K)
		logg: log of surface gravity (1-5)
		feh: metallicity relative to solar (-4 and 0.5 dex)
		micro: microturbulence (0-3 km s^-1)
		linelist: name of linelist file

	"""

	params = (teff, logg, feh, micro)

	if not in_bounds(params):
		print('Your parameters are not in the bounds of the model grids.')
		print()
		print(run_one.__doc__)
		return

	param_file(linelist)

	el_found, abundances = abunds_func(params)

	return
