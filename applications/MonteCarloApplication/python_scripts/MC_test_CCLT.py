from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats

# Import Monte Carlo library
import mc as mc

def compute_random_number():
    size_rv = 1
    random_number = np.random.uniform(1,4,size_rv)
    return random_number**4 * np.exp(random_number) + 1.0/np.sqrt(random_number*54382)

if __name__ == '__main__':
    convergence = False
    nsam = 0
    old_mean = 0.0
    old_M2 = 0.0
    Xlist=[]
    true_mean = 0.5*(1+4)
    true_variance = 1.0/12.0 * (4-1)**2
    while convergence is not True:
        nsam = nsam +1
        X = compute_random_number()
        mean_X, M2_X, sample_variance_X = mc.update_onepass_M(X, old_mean, old_M2, nsam)
        rel_err = (1.96*np.sqrt(sample_variance_X))/(mean_X*np.sqrt(nsam))
        if nsam > 30 and rel_err < 3e-3:
            convergence = True
        elif nsam < 1000000:
            convergence = False
        else:
            raise Exception("Convergence not reached")
        old_mean = mean_X
        old_M2 = M2_X
        
        value = (mean_X)
        value = np.asscalar(value)
        Xlist.append(value)
    print("nsam = ",nsam)
    print("relative error = ", rel_err)
    print("mean = ",mean_X, "true mean = ",true_mean)
    print("sample variance = ", sample_variance_X, "true_variance (no sample) = ", true_variance)

    mu = true_mean
    mu = mean_X
    sigma = np.sqrt(true_variance/nsam)
    sigma = np.sqrt(sample_variance_X/nsam)
    
    ## compute bins
    aux_bins = 1/np.sqrt(nsam)
    aux_bins = aux_bins * sigma * 3.49
    num_bins = int(6*sigma/aux_bins)
    print(num_bins)

    # n, bins, patches = plt.hist(Xlist,density=True, bins=250, range=(2.49,2.50))
    n, bins, patches = plt.hist(Xlist,density=True, bins=num_bins)

    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    plt.plot(x,mlab.normpdf(x, mu, sigma))

    plt.show()
    