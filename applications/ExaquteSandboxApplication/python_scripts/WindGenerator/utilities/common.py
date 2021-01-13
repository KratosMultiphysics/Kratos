from math import *
import numpy as np
from tqdm import tqdm


#=====================================================================================================
#=====================================================================================================
#
#                                         KERNELS
#
#=====================================================================================================
#=====================================================================================================

from scipy.special import kv as Kv
import scipy.special

#######################################################################################################
#	Matérn kernel
#######################################################################################################

def Mv(nu, X):
    # R = np.abs(X)
    # c0 = 2**(nu-1) * gamma(nu)
    # return np.where(nu<1.e-4, np.where(r>0,0,1), np.where(r>1.e-3, (r**nu) * Kv(nu, r) / c0,1))
    if nu < 1.e-4:
        m = lambda x:  1 if x==0 else 0
    else:
        m0 = 1 / (2**(nu-1) * gamma(nu))
        # m = lambda x:  (x**nu) * Kv(nu, x) / m0 if x>1.e-3 else 1
        def m(x):
            if x > 1.e-4:
                y = m0 * (x**nu) * Kv(nu, x)
            else:
                y = 1
            # y = m0 * sqrt(pi) * 2**(-nu) * exp(-x) * scipy.special.hyperu(1/2-nu, 1-2*nu, 2*x)
            # res = abs(sqrt(pi) * 2**(-nu) * exp(-x) * scipy.special.hyperu(1/2-nu, 1-2*nu, 2*x) - (x**nu) * Kv(nu, x))
            if y>=0 and y<=1:
                return y
            else:
                # print()
                # print(m0, sqrt(pi) * 2**(-nu) * exp(-x) * scipy.special.hyperu(1/2-nu, 1-2*nu, 2*x), x, y, nu)
                # print()
                return 1
    try:
        return np.array( [m(x) for x in np.abs(X)] )
    except:
        return m(abs(X))
    
    
def Matern_kernel(x, nu=1, rho=1):
    r = (np.abs(x)).flatten()
    if nu < 170:
        kappa = sqrt(2*nu)/rho
        y = Mv(nu, kappa*r)
    else:
        y = np.exp(-0.5*(r/rho)**2)
    try:
        y = y.reshape(x.shape)
    except:
        pass
    return y


#######################################################################################################
#	Shifted Matern kernel
#######################################################################################################

def SM_kernel(x, a):
    nu, rho, cw = a[0], a[1], a[2:]
    n = len(cw)//2
    cw = cw.reshape([-1,2])
    c, w = cw[:,0], cw[:,1]
    cosx = np.cos(w*np.tile(x, (n, 1)).T)
    y = (1 + np.inner(c,cosx)) * Matern_kernel(x, nu, rho)
    return y


#######################################################################################################
#	Generalized Matern kernel
#######################################################################################################

def GM_kernel(x, nu, rho, a):
    absx  = (np.abs(x)).flatten()
    kappa = sqrt(2*nu)/rho
    r = kappa * absx
    P = 1
    for k in range(len(a)):
        P += a[k] * r**(k+2)
    g = Matern_kernel(x, nu, rho) * P
    return g
    # if np.amin(g)>-1 and np.amax(g)<=1:
    #     return g
    # else:
    #     print("g fails:",np.amin(g), np.amax(g), P[np.argmax(g)], np.argmax(g))
    #     raise Exception("Invalid value for a covariance kernel.")
    


#######################################################################################################
#	Exponential-polynomial kernel
#######################################################################################################

def EP_kernel(x, a):
    a0, an = a[0], a[1:]
    n = len(an) #an.size
    xn = (np.tile(x, (n, 1)).T) ** (1+np.arange(n))
    y = (1 + np.inner(an,xn)) * np.exp(-a0*x)
    return y







#=====================================================================================================
#=====================================================================================================
#
#                                      GAUSSIAN LEVEL-CUT
#
#=====================================================================================================
#=====================================================================================================

from scipy.special import erf, erfinv, erfc, owens_t

#######################################################################################################
#	Volume fraction to Tau (and vice-versa)
#######################################################################################################

def vf2tau(vf, sigma=1, strategy=0):
    if strategy:
        # |u|<tau
        return sqrt(2)*sigma*erfinv(1-vf)
    else:
        # u<tau
        return sqrt(2)*sigma*erfinv(1-2*vf)

def tau2vf(tau, sigma=1, strategy=0):
    if strategy:
        # |u|<tau
        return 1-erf(tau/sqrt(2)/sigma)
    else:
        # u<tau
        return 0.5*(1-erf(tau/sqrt(2)/sigma))

def Cov2S2(tau, g, strategy=0):
    vf = tau2vf(tau, strategy=strategy)
    x  = np.sqrt((1-g)/(1+g))
    if strategy:
        # |u|<tau
        x1 = np.where(np.abs(x)<1.e-6, 1.e-6, x)
        S2 = 2*vf - 4*owens_t(tau, x) - 4*owens_t(tau, 1/x1)             
    else:
        # u<tau
        S2 = vf - 2 * owens_t(tau, x)
    return S2


#######################################################################################################
#	Fourier Transform of Gaussian Noise
#######################################################################################################

def FourierOfGaussian(noise):
    a, b = noise, noise
    for j in range(noise.ndim):
        b = np.roll(np.flip(b, axis=j), 1, axis=j)
    noise_hat = 0.5*( (a + b) + 1j*(a-b) )
    # for j in range(noise.ndim):
    #     np.roll(noise_hat, noise.shape[0] // 2 , axis=j)    
    return noise_hat


#######################################################################################################
#	Inclusion geometry
#######################################################################################################

def compute_Sphericity(V, A):
    return pi**(1/3) * (6*V)**(2/3) / A



 
#=====================================================================================================
#=====================================================================================================
#
#                                       PROBABILITY
#
#=====================================================================================================
#=====================================================================================================

import scipy.fftpack as fft
import scipy.optimize

#######################################################################################################
#	Basic probability tools
#######################################################################################################

def Expectation(X):
    nvalues, nsamples = X.shape
    m = np.zeros((nvalues,1))
    for isample in range(nsamples):
        for i in range(nvalues):
            m[i] += X[i,isample]
    m = m/nsamples
    return m

def Variance(X, m=None):
    if m is None: m = Expectation(X)
    sigma2 = Expectation(X**2) - m**2
    return sigma2

def SpacialCovariance(X):
    nvalues, nsamples = X.shape
    XX = np.zeros((nvalues, nsamples))
    for isample in tqdm(range(nsamples)):
        XX[:,isample] = X[:,isample]*X[0,isample]
    C = Expectation(XX)
    return C


#######################################################################################################
#	Compute probability distribution (from data)
#######################################################################################################

def compute_ProbaDist(data, bins=None):
    if bins is None:
        bins = 'auto'
    elif bins is 'log':
        _, bins = np.histogram(np.log(np.array(data) + 1), bins='auto')
        bins = np.exp(bins)-1
    elif bins is 'exp':
        _, bins = np.histogram(np.exp(np.array(data))-1, bins='auto')
        bins = np.log(bins + 1)
    else:
        bins = bins
    p, bin_edges = np.histogram(data, bins=bins, density=True)
    x = 0.5*(bin_edges[:-1] + bin_edges[1:])    
    yield p
    yield x


#######################################################################################################
#	Fit a probability with LogNormal (or Normal)
#######################################################################################################

def fit_ProbaDist(x, p, type='LogNormal'):
    if type is 'LogNormal':
        density = lambda x, m, sigma: dens_LogNormal(x, m=m, sigma=sigma)
    elif type is 'Normal':
        density = lambda x, m, sigma: dens_Normal(x, m=m, sigma=sigma)
    def F(th):
        m, sigma = th
        delta = density(x, m, sigma) - p
        return delta
        # return 0.5*np.dot(delta, delta)
    th = [np.argmax(p), p.size/2]
    result = scipy.optimize.least_squares(F, th)
    # result = scipy.optimize.minimize(F, th)
    m, sigma = result.x
    p_fit = density(x, m, sigma)
    yield m
    yield sigma
    yield p_fit


###################################################################################################
#   Autocorrelation of an image
###################################################################################################

def autocorrelation(X):
    xp = (X - 0*X.mean())

    if X.ndim==2:
        n, m = X.shape
        normalize = np.tensordot(np.arange(n)[::-1]+1, np.arange(m)[::-1]+1, axes=0)
        xp = np.pad(xp, ((0, n), (0, m)), 'constant')
        # xp = np.pad(xp, ((0, n), (0, m)), 'symmetric')
    elif  X.ndim==3:
        n, m, l = X.shape
        normalize = np.tensordot(np.arange(n)[::-1]+1, np.arange(m)[::-1]+1, axes=0)
        normalize = np.tensordot(normalize, np.arange(l)[::-1]+1, axes=0)
        xp = np.pad(xp, ((0, n), (0, m), (0, l)), 'constant')

    f = fft.fftn(xp)
    f2 = np.absolute(f)**2
    y = fft.ifftn(f2)
    
    if X.ndim==2:
        return y.real[:n, :m] / normalize
    elif  X.ndim==3:
        return y.real[:n, :m, :l] / normalize

def slope_by_fft(C):
    d = C.ndim
    Slope = 0
    for i in range(d):
        ni = C.shape[i]
        w = fft.fftfreq(ni)
        idx =[None,]*d
        idx[i] = slice(None,None,None)
        idx = tuple(idx)
        dC_i = fft.ifftn(1j*w[idx]*fft.fftn(C, axes=(i,)), axes=(i,))
        Slope += dC_i.flatten()[0].real
    Slope /= C.ndim
    return Slope



#######################################################################################################
#	Probability densities
#######################################################################################################

def dens_Exponential(x, lmbda=1):
    assert lmbda >= 0
    # if x < 0:
    #     p = 0
    # else:
    p = lmbda*np.exp(-lmbda*x)
    return p

def dens_Normal(x, m=0, sigma=1):
    p = 1/(sqrt(2*pi)*sigma) * np.exp(-0.5*((x-m)/sigma)**2)
    return p

def dens_LogNormal(x, m=0, sigma=1):
    # if x < 1.e-6:
    #     p = 0
    # else:
    p = 1/(sqrt(2*pi)*sigma*x) * np.exp(-0.5*((np.log(x)-m)/sigma)**2)
    return p





#######################################################################################################
#	Estimate covariance using Monte-Carlo
#######################################################################################################

def MC_estimate_Covariance(RandomField, nsamples=100, nbins=None):
    ### Stationary covariance only !
    C = 0
    for isample in tqdm(range(nsamples)):
        field = RandomField.sample()
        C += autocorrelation(field)
    C = C[0,:]/nsamples
    return C if nbins is None else C[:nbins]


#=====================================================================================================