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
    return random_number

if __name__ == '__main__':
    convergence = False
    nsam = 0
    old_mean = 0.0
    old_M2 = 0.0
    Xlist=[]
    mean_list=[]
    variance_list=[]
    nsam_list=[]
    true_mean = 0.5*(1+4)
    true_variance = 1.0/12.0 * (4-1)**2
    while convergence is not True:
        nsam = nsam +1
        X = compute_random_number()
        mean_X, M2_X, sample_variance_X = mc.update_onepass_M(X, old_mean, old_M2, nsam)
        rel_err = (1.96*np.sqrt(sample_variance_X))/(mean_X*np.sqrt(nsam))
        if nsam > 5000 and rel_err < 1e-2:
            convergence = True
        elif nsam > 5000:
            convergence = True
        else:
            convergence = False
        old_mean = mean_X
        old_M2 = M2_X
        if nsam >= 2:
            value = mean_X
            value = np.asscalar(value)
            Xlist.append(value)
    print("nsam = ",nsam)
    print("relative error = ", rel_err)
    print("mean = ",mean_X)
    print("sample variance = ", sample_variance_X)


    n, bins, patches = plt.hist(Xlist,density=True, bins=50)

    y = mlab.normpdf( bins, true_mean, true_variance*nsam)
    l = plt.plot(bins, y, 'r--', linewidth=2)

    plt.show()
    