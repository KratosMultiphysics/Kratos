from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
# import time

'''
This file contains all the functions to perform the Monte Carlo (MC) algorithm described in [PNL17]
References:
[PNL17] M. Pisaroni; F. Nobile; P. Leyland : A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
'''

'''
update mean and second moment values
M_{2,n} = sum_{i=1}^{n} (x_i - mean(x)_n)^2
M_{2,n} = M_{2,n-1} + (x_n - mean(x)_{n-1}) * (x_n - mean(x)_{n})
s_n^2 = M_{2,n} / (n-1)
'''
def update_onepass_M_VAR(sample, old_mean, old_M2, nsam):
    delta = np.subtract(sample, old_mean)
    if nsam == 1:
        new_mean = sample
        new_M2 = np.zeros(np.size(sample))
        new_M2 = np.asscalar(new_M2)
        new_sample_variance = np.zeros(np.size(sample))
        new_sample_variance = np.asscalar(new_sample_variance)
        '''do so to have a list of scalars, and not a list of arrays of one element'''
    else:
        new_mean = old_mean + np.divide(delta,nsam)
        new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
        new_sample_variance = compute_sample_variance_from_M2(new_M2,nsam)
    return new_mean, new_M2, new_sample_variance

def compute_sample_variance_from_M2(M2,nsam):
    sample_variance = np.divide(M2,np.subtract(nsam,1))
    return sample_variance