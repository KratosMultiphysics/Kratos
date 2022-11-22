import os
# os.environ['MKL_NUM_THREADS'] = '1'
# os.environ['OPENBLAS_NUM_THREADS']='1'
from mpi4py import MPI
import numpy as np
import numpy.linalg


def next_power_of_2(n):
    '''
    Find the next power of 2 of n
    '''
    p = 1
    if (n and not(n & (n - 1))):
        return n
    while (p < n):
        p <<= 1
    return p

# Ai and Bi are the local parts of A and B
# this function implements Ai.T@Bi and gathers it to root


def distributed_transpose_mult(Ai, Bi, comm, root):

    AitBi = Ai.T @ Bi  # multiplication of local blocks

    if(comm.rank == root):
        AtB = np.zeros(AitBi.shape, dtype=np.double)
    else:
        AtB = 0

    comm.Reduce([AitBi, MPI.DOUBLE], [AtB, MPI.DOUBLE], op=MPI.SUM, root=0)
    return AtB

# A -= Q Qt A

def OrthogonalProjector(Ai, Qi, comm):
    QtA = distributed_transpose_mult(Qi, Ai, comm, 0)
    QtA = comm.bcast(QtA, 0)
    Ai -= Qi @ QtA

# apply randomization

def Randomize(Ai, rand_size, comm):
    omega = 0
    np.random.seed(0)
    if(comm.Get_rank() == 0):    
        omega = np.random.random_sample((Ai.shape[1], rand_size))
    omega = comm.bcast(omega, 0)
    return Ai @ omega

def TotalNorm(Ai, comm):
    norm = np.array([0.0])
    local_norm2 = np.array([np.linalg.norm(Ai)**2],dtype=np.double)
    #comm.Allreduce([local_norm2, MPI.DOUBLE], [norm, MPI.DOUBLE], op=MPI.SUM)
    comm.Allreduce(local_norm2, norm, op=MPI.SUM)
    norm = np.sqrt(norm[0])
    return norm

def tsqr(Ai):
    '''
    Parallel QR factorization using Lapack
            Q(m,n) is the Q matrix
            R(n,n) is the R matrix
    '''

    # Recover rank and size
    MPI_COMM = MPI.COMM_WORLD      # Communications macro
    MPI_RANK = MPI_COMM.Get_rank()  # Who are you? who? who?
    MPI_SIZE = MPI_COMM.Get_size()  # Total number of processors used (workers)
    m, n = Ai.shape
    # Algorithm 1 from Demmel et al (2012)
    # 1: QR Factorization on Ai to obtain Q1i and Ri
    Q1i, R = numpy.linalg.qr(Ai)
    nextPower = next_power_of_2(MPI_SIZE)
    nlevels = int(np.log2(nextPower))
    QW = np.eye(n, dtype=np.double)
    C = np.zeros((2 * n, n), np.double)
    Q2l = np.zeros((2 * n * nlevels, n), np.double)
    blevel = 1
    for ilevel in range(nlevels):
        # Store R in the upper part of the C matrix
        C[:n, :] = R
        # Decide who sends and who recieves, use R as buffer
        prank = MPI_RANK ^ blevel
        if MPI_RANK & blevel:
            if prank < MPI_SIZE:
                MPI_COMM.send(R, prank, tag=0)
        else:
            if prank < MPI_SIZE:
                R = MPI_COMM.recv(source=prank, tag=0)
                # Store R in the lower part of the C matrix
                C[n:, :] = R
                # 2: QR from the C matrix, reuse C and R
                Q2i, R = np.linalg.qr(C)
                # Store Q2i from this level
                Q2l[2 * n * ilevel:2 * n * ilevel + 2 * n, :] = Q2i
        blevel <<= 1
    # At this point R is correct on processor 0
    # Broadcast R and its part of the Q matrix
    if MPI_SIZE > 1:    
        blevel = 1 << (nlevels - 1)
        mask = blevel - 1
    for ilevel in reversed(range(nlevels)):
        if MPI_RANK & mask == 0:
            # Obtain Q2i for this level - use C as buffer
            C = Q2l[2 * n * ilevel:2 * n * ilevel + 2 * n, :]
            # Multiply by QW either set to identity or allocated to a value
            # Store into Q2i
            Q2i = C @ QW  # matmul(C, QW)
            # Communications scheme
            prank = MPI_RANK ^ blevel
            if MPI_RANK & blevel:
                if prank < MPI_SIZE:
                    C = MPI_COMM.recv(source=prank, tag=0)
                    # Recover R from the upper part of C and QW from the lower
                    # part
                    R = C[:n, :]
                    QW = C[n:, :]
            else:
                if prank < MPI_SIZE:
                    # Set up C matrix for sending
                    # Store R in the upper part and Q2i on the lower part
                    # Store Q2i of this rank to QW
                    C[:n, :] = R
                    C[n:, :] = Q2i[n:, :]
                    QW = Q2i[:n, :]
                    MPI_COMM.send(C, prank, tag=0)
        blevel >>= 1
        mask >>= 1
    # Multiply Q1i and QW to obtain Qi
    Qi = Q1i @ QW  # matmul(Q1i, QW)
    return Qi



def randomized_orthogonalization(Ai, comm, mu=0, R=0):
    Mlocal, N = Ai.shape  # C has dimensions M by N
    
    #obtain M as the total size
    mlocal = np.array([Mlocal])
    mtotal = np.array([0])
    comm.Allreduce([mlocal, MPI.INT], [mtotal, MPI.INT], op=MPI.SUM)
    M = mtotal[0]
    
    c = nC = TotalNorm(Ai, comm)  # Norm of the initial residual

    if mu == 0:
        mu = max(M, N) * np.finfo(float).eps * nC / 2  # Machine precision parameter
    else:
        mu = nC*mu #it is relative to the initial norm

    dRmax = np.ceil(0.25 * (min(M, N)))
    dRmin = min(1, min(M, N))
    dRmin = max(dRmin, np.ceil(0.05 * min(M, N)))

    if(R==0):
        R = np.ceil(0.005 * (min(M, N)))  # Initial guess for the rank of C
    dR = int(R)
    TypeRankEstimate = 1  # Exponential
    i = 1  # iteration counter
    nC_old = c
    R_old = 0

    Q = None
    B = None

    while nC > mu:
        Qi = tsqr(Randomize(Ai, dR, comm))
        if Q is not None:
            for j in range(4):  # repeat application of projector to avoid numeric cancellation
                OrthogonalProjector(Qi, Q, comm)  # Qi=qr(Qi - Q @(Q.T @ Qi))
            Qi = tsqr(Qi)
        Bi = distributed_transpose_mult(Qi, Ai, comm, 0)
        Bi = comm.bcast(Bi, 0)  # broadcast Bi to all nodes
        
        Ai -= Qi @  Bi 
        if i == 1:
            Q = Qi
            if comm.Get_rank()==0:            
                B = Bi
        else:
            Q = np.c_[Q, Qi]  # da.concatenate([Q,Qi], axis=1)  #how to make this concatenation???
            if comm.Get_rank()==0:
                B = np.r_[B, Bi]  # da.concatenate([B,Bi], axis=0)
        nC = TotalNorm(Ai, comm)
        if(comm.Get_rank() == 0):
            print('iter = ', i, ' nC = ', nC, ' dR = ', dR, ' R = ', Q.shape[1])
        R_new = Q.shape[1]
        if TypeRankEstimate == 0:
            Rest = R_old + (R_new - R_old) / (nC - nC_old) * \
                (mu - nC_old)  # Estimated rank (linear)
        elif TypeRankEstimate == 1:
            Rest = R_old + (R_new - R_old) / (np.log(nC) - np.log(nC_old)) * (
                np.log(mu) - np.log(nC_old))  # Logarithmic Rank Estimation
        else:
            Rest = R_old + dRmin
        dR = np.ceil(Rest - R_new)
        dR = min(dR, dRmax)
        dR = max(dR, dRmin)
        Rest = R_new + dR
        dR = int(dR)
        if Rest >= N:
            print("FULL")
            Q = 'FULL'
            B = np.array([])
            # break
        i += 1
        R_old = R
        nC_old = nC
        R = Rest

    return Q, B



def svd_parallel(Q, B, comm, epsilon=0):

    if(comm.rank == 0):
        # u, s, _ = np.linalg.svd(B, full_matrices=False)
        M,N=np.shape(B)
        dimMATRIX = max(M,N)
        u, s, v = np.linalg.svd(B, full_matrices=False) #U --> M xN, V --> N x N
        v = v.T
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
        u = u[:, :K] 
        s = s[:K] 
        v = v[:, :K]
    else:
        u, s , v = [0,0,0]
    ub = comm.bcast(u, 0)
    s = comm.bcast(s, 0)
    U_final = Q @ ub

    return U_final,s,v

    

    
