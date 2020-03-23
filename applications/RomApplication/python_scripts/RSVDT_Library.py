### RANDOMIZED TRUNCATED SINGULAR VALUE DECOMPOSITION FUNCTIONS  ###

import scipy
from scipy import linalg
import numpy as np
import scipy.sparse.linalg as sp

###########################################################################################
def rsvdt( A, e0=0, mu=0, R=0, DATA={} ):    
    """
    RSVDT computes the truncated singular value decomposition (SVD) of matrix
    C = U*diag(S)*V' + E(e0,mu)
    
    INPUTS (mandatory)
    ------
    C --> M x N matrix whose factorization is to be calculated
    
    Optional inputs (with default values)
    ---------------

    1) e0 = 0
    2) mu = []
    3) R  = 0
    4) DATA.NITER = 10
    5) DATA.rho = 0.05
    6) DATA.COMPUTE_V_SVD = 0
    
    e0 --> Truncation threshold  (only the singular values and associated vectors that renders
                                  the frobenious norm of the residual E below this threshold are included in the factorization)

         If e0 = 0, we make e0 = mu (mu is defined below)
         
         mu --> Machine precision parameter
         If mu = [] , mu =  min(size(C))*eps(norm(C,'fro'))
         
         R is an estimation for an upper bound of rank(C) %
         if nargin == 1 or R=0,  then R = ceil(0.05*min(size(C))),
         
         DATA.NITER = Maximum number of iterations (rank revealing algorithm).
         
         DATA.COMPUTE_V_SVD --> Compute matrix of right singular vectors
         
         OUTPUT:
             ------
        U --> (truncated)  matrix of left singular vectors
        S ---> (truncated) vector of singular values
        V --> (truncated) matrix of right singular vectors
        e_svd --> Approximation error (frobenius norm)
        RankMatrix --> Rank matrix.
  ----------------------------------------------------------------
  Written in MATLAB by Joaquin  A. Hernandez, December 2016.
  UPC/CIMNE, Barcelona, Spain
  jhortega@cimne.upc.edu
  Python adaptation by Raul Bravo    
  jrbravo@cimne.upc.edu, June 2019    
  ------------------------------------------------------------------    
"""
    Rsup = min(np.shape(A))
    
    if 'NITER' not in DATA.keys():
        DATA['NITER'] = 10
    if 'rho' not in DATA.keys():
        DATA['rho'] = 0.05   
    if 'COMPUTE_V_SVD' not in DATA.keys():
        DATA['COMPUTE_V_SVD'] = 0
    if 'RELATIVE_SVD' not in DATA.keys():
        DATA['RELATIVE_SVD'] = False    
    if 'SVD' not in DATA.keys():
        DATA['SVD'] = {}
    if 'USE_ALWAYS_RANDOMIZATION' not in DATA.keys():
        DATA['USE_ALWAYS_RANDOMIZATION'] = 1
    
    
    if R >= Rsup and DATA['USE_ALWAYS_RANDOMIZATION']==0:
        Q='full'
    else:
        if R >= Rsup:
            #This means that the algorithm should take all possible modes.
            R=0
        Q, B, eORTH,a = rorth(A,mu,R,DATA)  #Randomized orthogonalization (machine precision parameter = mu)   
        if type(Q) != str:           
            Rsup = np.shape(Q)[1] 
    if len(Q)==0:
        U=S=V=np.array([])
        eSVD=0
    else:
        if type(Q) == str:
            #A appears to be full rank
            U,S,V,eSVD = svdt(A,e0,DATA)
        else:
            if DATA["RELATIVE_SVD"]==1 and e0>0:
                e0 = e0*a
            if e0 > eORTH:
                e=np.sqrt(e0**2 - eORTH**2)
            else:
                e=e0
            DATA["RELATIVE_SVD"] = 0
            DATA["SVD"]["MaxSize"] = max(np.shape(A))
            DATA["COMPUTE_V"]=DATA["COMPUTE_V_SVD"]
            U,S,V,eSVDb = svdt(B,e,DATA)
            U = Q @ U
            eSVD = np.sqrt(eORTH**2 + eSVDb**2)
            
    return U,S,V,eSVD

#--------------------------------------------------------------------------

def rorth(C, mu=0, R=0, DATA={}):
    """
    RORTH iteratively constructs an orthogonal basis matrix Q for the
    range of matrix C such that  norm(C-Q*B,'fro') <= mu, where B = Q^T*C,
    and Q^T*Q = I.
    
    Mandatory arguments:
        
        C --> M x N matrix
        
        Optional arguments
        
        R: estimation for  rank(C). Default value: R = ceil(DATA.rho_est*min(size(C)))
        mu:  Tolerance.  Default value: mu =   (max(size(C))*eps(nC))  where   nC = norm(C,'fro')
        DATA.NITER = Maximum number of iterations. Default = 10
        DATA.rho_est = Parameter defining the rank estimation R. Default value =
        0.05
        
        OUTPUT:
            ------
            Q --> Orthogonal basis matrix such that norm(C-Q*Q'*C,'fro') <= mu
            B = Q'*C
                                                        
        REMARK: If the estimated rank turns out to be greater than the number of
        columns of C, then randomization makes no sense (it is a full rank
        matrix) and the algorithm returns Q = 'FULL' and  B = [] ;

     This function is partially based on the randomized range finder algorithm
     proposed by Martisoon et al (2015) "A randomized blocked algorithm for efficiently computing
     rank-revealing factorizations of matrices".  The main novelty of our
     approach is the use of random matrices of varying size (the size of each random matrix is
     estimated based on linear rank predictions).

  ----------------------------------------------------------------
  Written in MATLAB by Joaquin  A. Hernandez, March 2017.
  UPC/CIMNE, Barcelona, Spain
  jhortega@cimne.upc.edu
  Python adaptation by Raul Bravo, June 2019    
  jrbravo@cimne.upc.edu   
  ------------------------------------------------------------------  
    
"""
    M,N=np.shape(C) # C has dimensions M by N
    c = nC = scipy.linalg.norm(C, 'fro')  # Norm of the initial residual
    if mu==0:
        mu = max(M,N)*np.finfo(float).eps*nC/2  # Machine presicion parameter
        
    if 'rho' not in DATA.keys():
        DATA['rho']=0.05        
    
    dRmax = np.ceil(0.25*(min(M,N)))
    dRmin = min(1,min(M,N))
    dRmin = max(dRmin, np.ceil(DATA['rho']*min(M,N) )  )
    
    if 'dRmax' not in DATA.keys():
        DATA['dRmax'] = dRmax
    if 'dRmin' not in DATA.keys():
        DATA['dRmin'] = dRmin   
    if 'R' not in DATA.keys():    
        DATA['R'] = np.ceil(0.005*(min(M,N))) #Initial guess for the rank of C
    if 'TypeRankEstimate' not in DATA.keys(): 
        DATA['TypeRankEstimate']= 1  #Exponential 
        
    #NITER=DATA["NITER"]
    if R==0:
        R=DATA['R']    
    dR = R
    i=1 #iteration counter
    nC_old = c  
    R_old = 0 
    Q = B = np.array([])
    DATA['COMPUTE_V']=False; DATA['RELATIVE_SVD']=False
    #No plotting
    while nC>mu:  
        Omega=np.random.RandomState(seed = 1234).normal(size=(N, int(dR))) # seeded random generator
        #Omega=np.random.RandomState().normal(size=(N, int(dR))) # Draw a N x dR random matrix 
        nOmega = np.sqrt(np.prod(np.shape(C)))
        factorRED = 10
        DATA['SVD']['MaxSize']=max(M,N)/factorRED
        Qi, _ = linalg.qr((C @ Omega)/nOmega, mode='economic') #Using QR to obtain the orthogonal basis
        #Qi,_,_,_ = svdt((C @ Omega)/nOmega,0,DATA) #Using trunctated svd to obtain the orthogonal basis (not suitable for the sparse implementation)
        
        if len(Qi)==0:
            break

        if len(Q) != 0:
            print('reorthogonalizing')
            DATA["SVD"]["MaxSize"]=max(np.shape(Qi))
            Qi, _ = linalg.qr(Qi - Q @(Q.T @ Qi), mode='economic')  #QR for re-orthogonalization
            #Qi,_,_,_= svdt(Qi - Q @(Q.T @ Qi),0,DATA)       #svdt for re-orthogonalization (not suitable for the sparse implementation)
        
        #Compute Residual        
        Bi = Qi.T @ C
        C = C - Qi @ Bi
    
        #Basis matrix is augmented with Qi
        if i == 1:
            Q = Qi
            B = Bi
        else:
            Q = np.c_[Q,Qi]
            B = np.r_[B,Bi]
        nC = scipy.linalg.norm(C, 'fro')  #Norm of the residual
        print('iter = ',i,' nC = ',nC,' dR = ',dR,' R = ', np.shape(Q)[1])  
        
        R_new=np.shape(Q)[1]
        
        if DATA["TypeRankEstimate"] ==0:
            Rest = R_old + (R_new-R_old)/(nC-nC_old)*(mu - nC_old)  #Estimated rank (linear)
        elif DATA["TypeRankEstimate"] ==1:
            Rest = R_old + (R_new -R_old)/(np.log(nC) - np.log(nC_old))*(np.log(mu) - np.log(nC_old))  #Logarithmic Rank Estimation
        else:
            Rest = R_old + DATA["dRmin"]    
    
        dR = np.ceil(Rest-R_new)
        dR = min(dR,DATA['dRmax'])
        dR = max(dR,DATA['dRmin']) 
        Rest = R_new+ dR 
        
        if Rest >= N:
            Q='FULL'; B=np.array([]); break
        
        i+=1
        R_old = R; nC_old = nC; R = Rest
        
        #if i >DATA[NITER]:
        #    break
    return Q, B, nC, c
#------------------------------------------------------------------------------------------



def svdt(B, epsilon=0, DATA={}):
    """ 
    SVDT(B,epsilon) computes a truncated SVD of B (with   tolerance epsilon)
    It employs the built-in matlab function "svd", but it only returns the
    singular values  truncated according to the   the specified tolerance

     S = SVDT(B,epsilon) gives the truncated vector of singular values

     [U,S] = SVDT(B,epsilon) furnishes the truncated matrix of left singular vectors and
     the associated singular values

     [U,S,V] = SVDT(B,epsilon)  returns the full factorization (truncated).

    By default,  norm(U*diag(S)*V'-B,'fro') <= epsilon
    If epsilon = 0, then the truncation criterion is the epsilon machine
    parameter
    If the third argument REL = 1, (by default REL = 0), then
    norm(U*diag(S)*V'-B,'fro') <= epsilon*norm(B,'fro')*
        
 
  ----------------------------------------------------------------
  Written in MATLAB by Joaquin  A. Hernandez, March 2017.
  UPC/CIMNE, Barcelona, Spain
  jhortega@cimne.upc.edu
  Python adaptation by Raul Bravo, June 2019    
  jrbravo@cimne.upc.edu  
  ------------------------------------------------------------------      
""" 
    M,N=np.shape(B)
    
    if 'COMPUTE_U' not in DATA.keys():
        DATA['COMPUTE_U'] = True
    if 'COMPUTE_V' not in DATA.keys():
        DATA['COMPUTE_V'] = False   
    if 'RELATIVE_SVD' not in DATA.keys():
        DATA['RELATIVE_SVD'] = 0
    if 'SVD' not in DATA.keys():
        DATA['SVD'] = {}
    if 'MaxSize' not in DATA['SVD'].keys():
        DATA['SVD']['MaxSize'] = max(M,N)
       
    CalcU=DATA['COMPUTE_U']
    CalcV=DATA['COMPUTE_V']
    dimMATRIX=DATA['SVD']['MaxSize']
    
    if DATA['TypeOfSVD'] == 0:   ## Using the regular scipy implementation (faster, but memory demanding)
        if M>=N:
            if CalcU==True and CalcV==True:
                U, s, V = linalg.svd(B, full_matrices=False) #U --> M xN, V --> N x N
                V = V.T
            elif CalcU==True and CalcV==False:
                U, s, _ = linalg.svd(B, full_matrices=False) #U --> M xN, V --> N x N           
            else:
                _, s, _ = linalg.svd(B, full_matrices=False)
        
        else: 
            # If N>M it proves more efficient to perform the SVD of B^T
            if CalcU==True and CalcV==True:
                V, s, U = linalg.svd(B.T, full_matrices=False) #U --> M x M, V --> Nx M 
                U = U.T
            elif CalcU==True and CalcV==False:
                _, s, U = linalg.svd(B.T, full_matrices=False) #U --> M x M, V --> Nx M            
                U = U.T
            else:
                _, s, _ = linalg.svd(B.T, full_matrices=False)

    elif DATA['TypeOfSVD'] == 1:   ## Using the sparse implementation (slower, but memory efficient. Based on ARPACK)
        r = min(B.shape)-1
        if r<1:
            r=2

        # Ontaining first approximation (Truncating potentially too many singular values)
        U, s, vt = sp.svds(B, k = min(B.shape)-1, tol=1e-12 )
      
        Truncate=0
        for i in range(len(s)):
            if s[i]==0:
                break
            Truncate+=1
        
        s= s[:Truncate]
        s= s[::-1]
        diagSigma = np.diag(s)
        U = U[:,:Truncate]
        U = U[:,::-1]
        vt = vt[:Truncate,:]
        vt = vt[::-1,:]

        B_reduced = U @ diagSigma @ vt
        B_star = B - B_reduced

        #normofstarmatrix = scipy.linalg.norm( B_star, 'fro') 
        #print(f'The norm of the matrix is: {normofstarmatrix }, it should be small')

        # Obtaining second approximation (to recover truncated singular values)
        if CalcU==True and CalcV==True: 
            U2, s2, vt2 = sp.svds(B_star, k = min(B_star.shape)-1 )
            Truncate=0
            for i in range(len(s2)):
                if s2[i]==0:
                    break
                Truncate+=1
            s2= s2[:Truncate]
            s2= s2[::-1]
            U2 = U2[:,:Truncate]
            U2 = U2[:,::-1]  
            vt2 = vt2[:Truncate,:]
            vt2 = vt2[::-1,:]            
            ## Concatenating the obtained basis and singular values
            s = np.r_[s,s2]
            U = np.c_[U,U2]  
            V = np.r_[vt,vt2]      
            V = V.T

        elif CalcU==True and CalcV==False:
            [U2, s2] = sp.svds(B_star, k = min(B_star.shape)-1 )[:2]
            Truncate=0
            for i in range(len(s2)):
                if s2[i]==0:
                    break
                Truncate+=1
            s2= s2[:Truncate]
            s2= s2[::-1]
            U2 = U2[:,:Truncate]
            U2 = U2[:,::-1]  
            ## Concatenating the obtained basis and singular values
            s = np.r_[s,s2]
            U = np.c_[U,U2]  

        else:
            [s2] = sp.svds(B_star, k = min(B_star.shape)-1 )[1]    
            Truncate=0
            for i in range(len(s2)):
                if s2[i]==0:
                    break
                Truncate+=1
            s2= s2[:Truncate]
            s2= s2[::-1]
            ## Concatenating the obtained singular values
            s = np.r_[s,s2]             

    tol = dimMATRIX*np.finfo(float).eps*max(s)/2    
    R = np.sum(s > tol)  # Definition of numerical rank
    eSVD = 0
    
    if epsilon == 0:
        K = R
    else:
        SingVsq = np.multiply(s,s)
        SingVsq.sort()
        normEf2 = np.sqrt(np.cumsum(SingVsq))
                
        if DATA['RELATIVE_SVD'] == 1:
            epsilon = epsilon*normEf2[-1]
            
        T = (sum(normEf2<epsilon))
        K = len(s)-T 
        
    K = min(R,K)
    
    if len(s)>K:
        eSVD = np.sqrt(sum( (s[K+1:len(s)-1])**2  ))


    if CalcU==True and CalcV==True:
        return U[:, :K], s[:K], V[:K, :], eSVD
    elif CalcU==True and CalcV==False:
        return U[:, :K], s[:K], np.nan , eSVD
    else:
        return np.nan, s[:K], np.nan , eSVD   
        
        
#----------------------------------------------------------------------------------------------


############################################################################################
# #Computes the orthonormal matrix whose range approximates the range of A. Halko et. al. 2009
# def randomized_range_finder(A, size, n_iter=5):
#     Q = np.random.normal(size=(A.shape[1], size))    
#     for i in range(n_iter):
#         Q, _ = linalg.lu(A @ Q, permute_l=True)
#         Q, _ = linalg.lu(A.T @ Q, permute_l=True)        
#     Q, _ = linalg.qr(A @ Q, mode='economic')
#     return Q

# def randomized_svd(M, n_components, n_oversamples=10, n_iter=4):    
#     n_random = n_components + n_oversamples    
#     Q = randomized_range_finder(M, n_random, n_iter)    
#     # project M to the (k + p) dimensional space using the basis vectors
#     B = Q.T @ M    
#     # compute the SVD on the thin matrix: (k + p) wide
#     Uhat, s, V = linalg.svd(B, full_matrices=False)
#     del B
#     U = Q @ Uhat    
#     return U[:, :n_components], s[:n_components], V[:n_components, :]
############################################################################################

























