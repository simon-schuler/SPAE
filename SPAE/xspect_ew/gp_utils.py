"""Gaussian Process regression used for per-line EW measurement.

Separate from the continuum-fitting GP in continuum.py/spectrum_data.py,
which uses the george package directly -- this is a from-scratch GP
implementation (kernel + prediction + negative log likelihood), taken from
an LSSTC DSFP notebook, used specifically inside measure_ew() to smooth the
flattened line profile before the Gaussian fit.
"""

import numpy as np
from scipy.spatial.distance import cdist
from numpy.linalg import inv
from numpy.linalg import slogdet


def SEKernel(par, x1, x2):
    A, Gamma = par
    D2 = cdist(x1.reshape(len(x1),1), x2.reshape(len(x2),1), metric = 'sqeuclidean')
    return A*np.exp(-Gamma*D2)


def Pred_GP(CovFunc, CovPar, xobs, yobs, eobs, xtest):
    # evaluate the covariance matrix for pairs of observed inputs
    K = CovFunc(CovPar, xobs, xobs)
    # add white noise
    K += np.identity(xobs.shape[0]) * eobs**2
    # evaluate the covariance matrix for pairs of test inputs
    Kss = CovFunc(CovPar, xtest, xtest)
    # evaluate the cross-term
    Ks = CovFunc(CovPar, xtest, xobs)
    # invert K
    Ki = inv(K)
    # evaluate the predictive mean
    m = np.dot(Ks, np.dot(Ki, yobs))
    # evaluate the covariance
    cov = Kss - np.dot(Ks, np.dot(Ki, Ks.T))
    return m, cov


def NLL_GP(p,CovFunc,x,y,e):
    # Evaluate the covariance matrix
    K = CovFunc(p,x,x)
    # Add the white noise term
    K += np.identity(x.shape[0]) * e**2
    # invert it
    Ki = inv(K)
    # evaluate each of the three terms in the NLL
    term1 = 0.5 * np.dot(y,np.dot(Ki,y))
    term2 = 0.5 * slogdet(K)[1]
    term3 = 0.5 * len(y) * np.log(2*np.pi)
    # return the total
    return term1 + term2 + term3
