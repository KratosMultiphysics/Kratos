import scipy
from scipy import linalg
import numpy as np

class RandomizedSingularValueDecomposition():

    def __init__(self, COMPUTE_U=True, COMPUTE_V=True, RELATIVE_SVD=True, USE_RANDOMIZATION=True):
        self.COMPUTE_U = COMPUTE_U
        self.COMPUTE_V = COMPUTE_V
        self.RELATIVE_SVD = RELATIVE_SVD
        self.USE_RANDOMIZATION = USE_RANDOMIZATION

    def Calculate(self, A, e0 = 0):
        Rsup = min(np.shape( A))
        if self.USE_RANDOMIZATION == False:
            Q='full'
        else:
            Q, B, eORTH, a = self._RandomizedOrthogonalization(A)  #Randomized orthogonalization (machine precision parameter = mu)
            if type(Q) != str:
                Rsup = np.shape(Q)[1]
        if len(Q)==0:
            U=S=V=np.array([])
            eSVD=0
        else:
            if type(Q) == str:
                #A appears to be full rank
                U,S,V,eSVD = self._SingularValueDecompostionTruncated( A, e0)
            else:
                if self.RELATIVE_SVD==1 and e0>0:
                    e0 = e0*a
                if e0 > eORTH:
                    e=np.sqrt(e0**2 - eORTH**2)
                else:
                    e=e0
                self.RELATIVE_SVD = 0
                self.SVD_MaxSize = max(np.shape(A))
                U,S,V,eSVDb = self._SingularValueDecompostionTruncated(B, e )
                U = Q @ U
                eSVD = np.sqrt(eORTH**2 + eSVDb**2)

        return U,S,V,eSVD

    def _RandomizedOrthogonalization(self, C , mu=0, R=0):

        M,N=np.shape(C) # C has dimensions M by N
        c = nC = scipy.linalg.norm(C, 'fro')  # Norm of the initial residual
        if mu==0:
            mu = max(M,N)*np.finfo(float).eps*nC/2  # Machine presicion parameter

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
            Qi, _ = linalg.qr((C @ Omega)/nOmega, mode='economic') #Using QR to obtain the orthogonal basis
            #Qi,_,_,_ = self._SingularValueDecompostionTruncated((C @ Omega)/nOmega) #Using trunctated svd to obtain the orthogonal basis (not suitable for the sparse implementation)

            if len(Qi)==0:
                break

            if len(Q) != 0:
                print('reorthogonalizing')
                self.SVD_MaxSize = max(np.shape(Qi))
                Qi, _ = linalg.qr(Qi - Q @(Q.T @ Qi), mode='economic')  #QR for re-orthogonalization
                #Qi,_,_,_= self._SingularValueDecompostionTruncated(Qi - Q @(Q.T @ Qi))       #svdt for re-orthogonalization (not suitable for the sparse implementation)

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
            return U[:, :K], s[:K], V[:K, :], eSVD
        elif CalcU==True and CalcV==False:
            return U[:, :K], s[:K], np.nan , eSVD
        else:
            return np.nan, s[:K], np.nan , eSVD




if __name__ == "__main__":

    # #load a test SnapshotMatrix

    # RSVDT_Object = RandomizedSingularValueDecomposition()
    # U,S,V,error = RSVDT_Object.Calculate(SnapshotMatrix, 1e-6)
    # APPROX_SnapshotMatrix = U@np.diag(S)@V

    # u,s,v = np.linalg.svd(SnapshotMatrix,full_matrices=False)
    # approx_SnapshotMatrix = u@np.diag(s)@v.T

    # print(f"\n\n\n Norm of the difference in implementations: { np.linalg.norm(APPROX_SnapshotMatrix-approx_SnapshotMatrix)}\n\n\n")
    # print(f"error accroding to RSVDT = {error}")
    pass

