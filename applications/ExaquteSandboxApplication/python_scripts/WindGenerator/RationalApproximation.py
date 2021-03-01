
from math import *
import numpy as np
from scipy.optimize import fsolve, minimize, fminbound, root
from scipy.signal import residue
import scipy
# from scipy.special import kv as Kv
# from scipy.special import jv as Jv
# from scipy.special import iv as Iv
# from scipy.special import gamma, genlaguerre
import matplotlib.pyplot as plt
from time import time


############################################


def compute_RationalApproximation_AAA(alpha, MaxOrder=100, tol = 1.e-8, nPoints = 1000, verbose=False):

    # deal with divisions by zero
    np.seterr(divide='ignore', invalid='ignore')

    func = lambda x: x**(1-alpha)
    # func = lambda x: (1+x)**(-alpha) * x

    M = nPoints
    Z = np.linspace(1, 1000, M)
    # Z = (Z/Z[-1])**2 * Z[-1]
    Z = Z.reshape([M,1])
    M = Z.size
    F = func(Z)
    SF = scipy.sparse.spdiags(F.flatten(),0,M,M)
    J = np.arange(M)
    z = np.array([]).reshape([0,1])
    f = np.array([]).reshape([0,1])
    C = np.array([]).reshape([M,0])
    errvec = np.array([])
    R = np.mean(F)*np.ones_like(F)
    for m in range(MaxOrder):
        j = np.argmax(np.abs(F-R))
        z = np.vstack([z, Z[j]])
        f = np.vstack([f, F[j]])
        J = J[J!=j]
        C = np.hstack([ C, 1/(Z-Z[j]) ])
        Sf = np.diag(f.flatten())
        A = SF @ C - C @ Sf
        U,S,V = scipy.linalg.svd(A[J,:])
        w = V[m,:].reshape([-1,1])
        N = C @ (w*f)
        D = C @ w
        R[:] = F
        R[J] = N[J]/D[J]
        err = np.linalg.norm(F-R, ord=np.inf)
        errvec = np.append(errvec, err)
        if verbose: print(err)
        if err <= tol: break
    if verbose: print('degree',m)


    m = w.size
    B = np.eye(m+1)
    B[0,0] = 0
    E = np.block([ [ 0, w.reshape([1,m]) ], [ np.ones([m,1]), np.diag(z.flatten()) ] ])
    pol = scipy.linalg.eig(E,B, left=False, right=False)
    pol = pol[~np.isinf(pol)]
    E = np.block([ [ 0, (w*f).reshape([1,m]) ], [ np.ones([m,1]), np.diag(z.flatten()) ] ])
    zer = scipy.linalg.eig(E,B, left=False, right=False)
    zer = zer[~np.isinf(zer)]

    assert( np.all(np.isclose(np.imag(zer),0)) )
    assert( np.all(np.isclose(np.imag(pol),0)) )

    pc0 = np.sum(w*f)
    pd0 = np.sum(w)

    pc = np.poly(zer)
    pd = np.poly(pol)

    pd = np.append(pd, 0)

    c, d, k = residue(pc, pd)
    assert(all(k==0))

    d = -d
    c *= pc0/pd0

    ### Check RA error
    # x = np.linspace(1, 1000, 10000)
    # plt.plot(x, [ np.sum(c/(z+d)) - z**(-alpha) for z in x] )
    # plt.show()

    return c, d


def compute_RationalApproximation_AAA_new(alpha, beta=None, MaxOrder=100, tol = 1.e-12, nPoints = 1000, verbose=False):

    # deal with divisions by zero
    np.seterr(divide='ignore', invalid='ignore')

    if beta:
        func = lambda x: x**(1-alpha) * (x-1)**(-beta)
    else:
        func = lambda x: x**(1-alpha)

    M = nPoints
    Z = np.linspace(1+1.e-3, 1000, M).reshape([M,1])
    M = Z.size
    F = func(Z)
    SF = scipy.sparse.spdiags(F.flatten(),0,M,M)
    J = np.arange(M)
    z = np.array([]).reshape([0,1])
    f = np.array([]).reshape([0,1])
    C = np.array([]).reshape([M,0])
    errvec = np.array([])
    R = np.mean(F)*np.ones_like(F)
    for m in range(MaxOrder):
        j = np.argmax(np.abs(F-R))
        z = np.vstack([z, Z[j]])
        f = np.vstack([f, F[j]])
        J = J[J!=j]
        C = np.hstack([ C, 1/(Z-Z[j]) ])
        Sf = np.diag(f.flatten())
        A = SF @ C - C @ Sf
        U,S,V = scipy.linalg.svd(A[J,:])
        w = V[m,:].reshape([-1,1])
        N = C @ (w*f)
        D = C @ w
        R[:] = F
        R[J] = N[J]/D[J]
        err = np.linalg.norm(F-R, ord=np.inf)
        errvec = np.append(errvec, err)
        if verbose: print(err)
        if err <= tol: break
    if verbose: print('degree',m)


    m = w.size
    B = np.eye(m+1)
    B[0,0] = 0
    E = np.block([ [ 0, w.reshape([1,m]) ], [ np.ones([m,1]), np.diag(z.flatten()) ] ])
    pol = scipy.linalg.eig(E,B, left=False, right=False)
    pol = pol[~np.isinf(pol)]
    E = np.block([ [ 0, (w*f).reshape([1,m]) ], [ np.ones([m,1]), np.diag(z.flatten()) ] ])
    zer = scipy.linalg.eig(E,B, left=False, right=False)
    zer = zer[~np.isinf(zer)]

    assert( np.all(np.isclose(np.imag(zer),0)) )
    assert( np.all(np.isclose(np.imag(pol),0)) )

    pc0 = np.sum(w*f)
    pd0 = np.sum(w)

    pc = np.poly(zer)
    pd = np.poly(pol)

    pd = np.append(pd, 0)

    c, d, k = residue(pc, pd)
    assert(all(k==0))

    d = -d
    c *= pc0/pd0

    # x = np.linspace(100, 1000, 10000)
    # plt.plot(x, [ np.sum(c/(z+d)) - func(z)/z for z in x] )
    # plt.show()

    return c, d