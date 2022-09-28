import dislib as ds
import numpy as np
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import Type, COLLECTION_IN, Depth
from dislib.data.array import Array

from KratosMultiphysics.RomApplication.auxiliary_functions_workflow import load_blocks_array, load_blocks_rechunk


@constraint(computingUnits=2)
@task(Y_blocks={Type: COLLECTION_IN, Depth: 2}, returns=2)
def my_qr(Y_blocks):
    Y = np.block(Y_blocks)
    Q,R = np.linalg.qr(Y, mode='reduced')
    return Q,R


@constraint(computingUnits=1)
@task(B_blocks={Type: COLLECTION_IN, Depth: 2}, returns=2)
def my_svd(B_blocks):
    B = np.block(B_blocks)
    U_hat, s, _ = np.linalg.svd(B, full_matrices=False)
    return U_hat, s




def rsvd(A,desired_rank,oversampling=10,row_splits=10,column_splits=1):

#-----Dimensions--------
    k = desired_rank
    p = k + oversampling
    n = A.shape[0]
    m = A.shape[1]
    A_row_chunk_size = int( n / row_splits)
    A_column_chunk_size = int( m / column_splits)
    A = A.rechunk((A_row_chunk_size,A_column_chunk_size))
#-----Matrix Omega Initialization--------
    omega_column_chunk_size = p
    omega_row_chunk_size = A_column_chunk_size
    Omega = ds.random_array(shape=(m, p), block_size=(omega_row_chunk_size, omega_column_chunk_size) ) # Create a random projection matrix Omega of size mxp, for this test, p is of 110.
#----------------------------------

# STEP 1 (DISTRIBUTED): Sample the column space of A. Y results into an nxp matrix.
    Y = A @ Omega

# STEP 2 (SERIAL): Serial QR, to be done in distributed. Q results into a matrix of nxp
    Q,_ = my_qr(Y._blocks)
    Q=load_blocks_rechunk([Q], shape = (n, p), block_size= (n, p), new_block_size=(A_row_chunk_size, omega_column_chunk_size))

#"SERIAL STEPS FINISH"
    B = Q.T @ A
# STEP 3 (DISTRIBUTED): Project A into the orthonormal basis. B results into a matrix of pxm.

#""""SERIAL STEPS START"""
    U_hat, s = my_svd(B._blocks)
    U_hat = load_blocks_rechunk([U_hat], shape = (p, m), block_size = (p, m), new_block_size=(omega_column_chunk_size, omega_column_chunk_size))

#""""SERIAL STEPS FINISH"""
    U_hat = U_hat[:,:k] #U_hat results into a matrix of pxk.
    U = Q @ U_hat  #STEP 5 (DISTRIBUTED): Project the reduced basis into Q. U results into a matrix of nxk.
    U = U.collect() #STEP 6 (collecting the required basis):
    s = compss_wait_on(s)
    s = s[:k]
    return U, s

def truncated_svd(Matrix, epsilon=0):
    Matrix = compss_wait_on(Matrix)
    M,N=np.shape(Matrix)
    dimMATRIX = max(M,N)
    U, s, V = np.linalg.svd(Matrix, full_matrices=False) #U --> M xN, V --> N x N
    V = V.T
    tol = dimMATRIX*np.finfo(float).eps*max(s)/2
    R = np.sum(s > tol)  # Definition of numerical rank
    if epsilon == 0:
        K = R
    else:
        SingVsq = np.multiply(s,s)
        SingVsq.sort()
        normEf2 = np.sqrt(np.cumsum(SingVsq))
        epsilon = epsilon*normEf2[-1] #relative tolerance
        T = (sum(normEf2<epsilon))
        K = len(s)-T
    K = min(R,K)
    return U[:, :K], s[:K], V[:, :K]

