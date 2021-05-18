import numpy as np

class RandomizedSingularValueDecomposition():
    """
    This class calculates the singular value decomposition of a matrix A (A = U@np.diag(S)@V.T + Error(truncation_tolerance)) using a randomized algorithm
    Reference: Halko et al 2009. "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions"
    """



    """
    Constructor setting up the parameters for calculation of the SVD
        COMPUTE_U: whether to return the matrix of left singular vectors U
        COMPUTE_V: whether to return the matrix of right singular vectors V
        RELATIVE_SVD: If true, the truncation_tolerance is multiplied by the norm of the original matrix
                    If false, the truncation_tolerance is taken as an absolute error tolerance
        USE_RANDOMIZATION: If false, the standard svd algorith of numpy is followed
    """
    def __init__(self, COMPUTE_U=True, COMPUTE_V=True, RELATIVE_SVD=True, USE_RANDOMIZATION=True):
        self.COMPUTE_U = COMPUTE_U
        self.COMPUTE_V = COMPUTE_V
        self.RELATIVE_SVD = RELATIVE_SVD
        self.USE_RANDOMIZATION = USE_RANDOMIZATION



    """
    Method for calculating the SVD
    input:  A: numpy array containing the matrix to decompose
            truncation_tolerance: this parameter is employed (as is or multiplied by the norm of A, depending on self.RELATIVE_SVD) to truncate the approximation
    output: U: numpy array containing the matrix of left singular vectors
            S: numpy array containing the matrix of singular values
            V: numpy array containing the matrix of right singular vectors
            eSVD : estimation of the error of the approximation
    """
    def Calculate(self, A, truncation_tolerance = 0):
        if self.USE_RANDOMIZATION == False:
            Q='full'
        else:
            Q, B, eORTH, a = self._RandomizedOrthogonalization(A)  #Randomized orthogonalization (machine precision parameter = mu)
        if len(Q)==0:
            U=S=V=np.array([])
            eSVD=0
        else:
            if isinstance(Q,str):
                #A appears to be full rank
                U,S,V,eSVD = self._SingularValueDecompostionTruncated( A, truncation_tolerance)
            else:
                if self.RELATIVE_SVD==1 and truncation_tolerance>0:
                    truncation_tolerance = truncation_tolerance*a
                if truncation_tolerance > eORTH:
                    e=np.sqrt(truncation_tolerance**2 - eORTH**2)
                else:
                    e=truncation_tolerance
                self.RELATIVE_SVD = 0
                self.SVD_MaxSize = max(np.shape(A))
                U,S,V,eSVDb = self._SingularValueDecompostionTruncated(B, e )
                U = Q @ U
                eSVD = np.sqrt(eORTH**2 + eSVDb**2)

        return U,S,V,eSVD



    """
    Method for obtaining an othonormal basis for the range of the matrix to decompose (stage A, from reference Halko et al 2009 )
    input:  C: numpy array containing the matrix from which the orthonornal basis will be obtained
            mu: machine precision parameter, if not specified it is estimated
            R: estimation for the rank of C, if not specified it is estimated
    output: Q: numpy array containing the orthonormal basis such that norm(C - Q@Q.T@C) <= mu
            B: numpy array Q.T@C
            nC: numpy array estimation of the orthogonalization error
            c : norm of C
    """
    def _RandomizedOrthogonalization(self, C , mu=0, R=0):

        M,N=np.shape(C) # C has dimensions M by N
        c = nC = np.linalg.norm(C, 'fro')  # Norm of the initial residual
        if mu==0:
            mu = max(M,N)*np.finfo(float).eps*nC/2  # Machine precision parameter

        dRmax = np.ceil(0.25*(min(M,N)))
        dRmin = min(1,min(M,N))
        dRmin = max(dRmin, np.ceil(0.05*min(M,N) )  )

        self.dRmax = dRmax
        self.dRmin = dRmin
        self.R = np.ceil(0.005*(min(M,N))) #Initial guess for the rank of C
        self.TypeRankEstimate = 1 #Exponential

        if R==0:
            R=self.R
        dR = R
        i=1 #iteration counter
        nC_old = c
        R_old = 0
        Q = B = np.array([])

        while nC>mu:
            #Omega=np.random.RandomState(seed = 1234).normal(size=(N, int(dR))) # seeded random generator
            Omega=np.random.RandomState().normal(size=(N, int(dR))) # Draw a N x dR random matrix
            nOmega = np.sqrt(np.prod(np.shape(C)))
            factorRED = 10
            self.SVD_MaxSize = max(M,N)/factorRED
            Qi, _ = np.linalg.qr((C @ Omega)/nOmega, mode='reduced') #Using QR to obtain the orthogonal basis
            #Qi,_,_,_ = self._SingularValueDecompostionTruncated((C @ Omega)/nOmega) #Using trunctated svd to obtain the orthogonal basis

            if len(Qi)==0:
                break

            if len(Q) != 0:
                print('reorthogonalizing')
                self.SVD_MaxSize = max(np.shape(Qi))
                Qi, _ = np.linalg.qr(Qi - Q @(Q.T @ Qi), mode='reduced')  #QR for re-orthogonalization
                #Qi,_,_,_= self._SingularValueDecompostionTruncated(Qi - Q @(Q.T @ Qi))       #svdt for re-orthogonalization

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
            nC = np.linalg.norm(C, 'fro')  #Norm of the residual
            print('iter = ',i,' nC = ',nC,' dR = ',dR,' R = ', np.shape(Q)[1])

            R_new=np.shape(Q)[1]

            if self.TypeRankEstimate == 0:
                Rest = R_old + (R_new-R_old)/(nC-nC_old)*(mu - nC_old)  #Estimated rank (linear)
            elif self.TypeRankEstimate == 1:
                Rest = R_old + (R_new -R_old)/(np.log(nC) - np.log(nC_old))*(np.log(mu) - np.log(nC_old))  #Logarithmic Rank Estimation
            else:
                Rest = R_old + self.dRmin

            dR = np.ceil(Rest-R_new)
            dR = min(dR,self.dRmax)
            dR = max(dR,self.dRmin)
            Rest = R_new+ dR

            if Rest >= N:
                Q='FULL'; B=np.array([]); break

            i+=1
            R_old = R; nC_old = nC; R = Rest

        return Q, B, nC, c


    """
    Method for calculating a truncated version of numpy's svd
    input:  B: numpy array containing a matrix to decompose
            epsilon: truncation tolerance for the svd
    output: U: numpy array containing the matrix of left singular vectors. If self.COMPUTE_U=False ==> U = np.nan
            S: numpy array containing the matrix of singular values
            V: numpy array containing the matrix of right singular vectors. If self.COMPUTE_V=False ==> V = np.nan
            eSVD : estimation of the error of the approximation
    """
    def _SingularValueDecompostionTruncated(self, B, epsilon = 0):

        M,N=np.shape(B)

        CalcU=self.COMPUTE_U
        CalcV=self.COMPUTE_V
        if hasattr(self, 'SVD_MaxSize'):
            dimMATRIX = self.SVD_MaxSize
        else:
            dimMATRIX = max(M,N)

        if M>=N:
            if CalcU==True and CalcV==True:
                U, s, V = np.linalg.svd(B, full_matrices=False) #U --> M xN, V --> N x N
                V = V.T
            elif CalcU==True and CalcV==False:
                U, s, _ = np.linalg.svd(B, full_matrices=False) #U --> M xN, V --> N x N
            else:
                _, s, _ = np.linalg.svd(B, full_matrices=False)

        else:
            # If N>M it proves more efficient to perform the SVD of B^T
            if CalcU==True and CalcV==True:
                V, s, U = np.linalg.svd(B.T, full_matrices=False) #U --> M x M, V --> Nx M
                U = U.T
            elif CalcU==True and CalcV==False:
                _, s, U = np.linalg.svd(B.T, full_matrices=False) #U --> M x M, V --> Nx M
                U = U.T
            else:
                _, s, _ = np.linalg.svd(B.T, full_matrices=False)

        tol = dimMATRIX*np.finfo(float).eps*max(s)/2
        R = np.sum(s > tol)  # Definition of numerical rank
        eSVD = 0

        if epsilon == 0:
            K = R
        else:
            SingVsq = np.multiply(s,s)
            SingVsq.sort()
            normEf2 = np.sqrt(np.cumsum(SingVsq))

            if self.RELATIVE_SVD == 1:
                epsilon = epsilon*normEf2[-1]

            T = (sum(normEf2<epsilon))
            K = len(s)-T

        K = min(R,K)

        if len(s)>K:
            eSVD = np.sqrt(sum( (s[K+1:len(s)-1])**2  ))


        if CalcU==True and CalcV==True:
            return U[:, :K], s[:K], V[:, :K], eSVD
        elif CalcU==True and CalcV==False:
            return U[:, :K], s[:K], np.nan , eSVD
        else:
            return np.nan, s[:K], np.nan , eSVD