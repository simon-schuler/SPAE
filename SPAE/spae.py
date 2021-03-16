import numpy as np
import emcee
import time
from .abunds import abunds_func, obj_func

def run_spae(sun_el=None,sun_abs=None,x_0=(5777,4.44,0.01,1.38), n_dim=4, n_walkers=40, n_steps=1000):

	t_0 = time.time()

	x_ball = x_0 * (1 + 1.0e-2 * (-0.5 + np.random.rand(n_walkers, n_dim)))

	el_found, abundances = abunds_func(x_0)
	blobs_dtype = [('ep_r','f8'), ('rew_r','f8')]

	for el in el_found:
	    blobs_dtype.append((el, 'f8'))
	    blobs_dtype.append((el+'_sigma_mean', 'f8'))


	sampler = emcee.EnsembleSampler(n_walkers, n_dim, obj_func, blobs_dtype=blobs_dtype, args=(len(el_found),sun_el,sun_abs))
	results = sampler.run_mcmc(x_ball, n_steps, progress=True)

	flat_blob = sampler.get_blobs().reshape((n_walkers*n_steps))

	log = 'Runtime: ' + str(time.time() - t_0) + '\n'
	log = log + 'Size of Full Array: ' + str(len(flat_blob)) + '\n'
	log = log + 'Acceptance Fractions:' + np.array_str(sampler.acceptance_fraction) + '\n'

	return sampler, flat_blob, log
