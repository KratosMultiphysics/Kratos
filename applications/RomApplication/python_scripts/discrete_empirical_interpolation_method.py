# Discrete Empirical Interpolation Method DEIM
# Based on Chaturantabut & Sorensen, 2010, doi.org/10.1137/090766498

#Q-DEIM
#Based on Drmac & Gugercin, 2016, doi.org/10.1137/15M1019271


import numpy as np

try:
    from scipy.linalg import qr
    missing_scipy = False
except ImportError as e:
    missing_scipy = True



def DEIM(Basis):
    #find first point
    U = Basis[:,0].reshape(Basis.shape[0], -1)
    z = np.zeros(U.shape)
    P = z.copy()
    indexes = np.argmax( np.abs(Basis[:,0]) )
    P[indexes] = 1

    #find next points
    for i in range(1,Basis.shape[1]):
        c = np.linalg.solve(P.T @ U , P.T @ Basis[:,i] )
        residual = Basis[:,i] - U @ c
        U = np.c_[U,Basis[:,i]]
        index_i = np.argmax( np.abs(residual) )
        indexes = np.r_[indexes, index_i]
        P = np.c_[P, z]; P[index_i,i] = 1

    return indexes


def QDEIM(Basis):
    if missing_scipy==False:
        _,_,pivots = qr(Basis.T, pivoting=True)
        return pivots[:Basis.shape[1]]
    else:
        return DEIM(Basis)