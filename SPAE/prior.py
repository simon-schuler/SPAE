import numpy as np
import kalepy as kale
from . import data

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

prior_kde = None


def load_prior(prior_kwargs={}):
    """Load the prior data and initialize the prior array."""
    global prior_kde

    # Load up prior data and create input array
    filename = "Brewer_data.npy"
    prior_data = np.load(pkg_resources.open_binary(data, filename))

    # Cycle through kwargs
    for key, value in prior_kwargs.items():
        if key == 'log_Rhk_min' and value is not None:
            prior_data = prior_data[prior_data['logRhk'] > value]
        if key == 'log_Rhk_max' and value is not None:
            prior_data = prior_data[prior_data['logRhk'] < value]

    # Create array of inputs for prior
    prior_inputs = np.array([prior_data['Teff'] / 1e3, prior_data['logg']])

    # Generate KDE from input data
    prior_kde = kale.KDE(prior_inputs, kernel=None, bandwidth=0.1)

    return


def ln_prior(p, prior_kwargs={}):
    """Calculate the log prior from the KDE representation."""
    global prior_kde

    # Check to make sure prior is loaded
    if prior_kde is None:
        load_prior(prior_kwargs)

    T_eff, logg = p

    prob = prior_kde.density([[T_eff/1e3], [logg]], probability=True)[1]
    if prob == 0:
        return -np.inf

    return np.log(prob)[0]
