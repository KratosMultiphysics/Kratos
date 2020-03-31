from KratosMultiphysics.RomApplication.element_selection_strategy import ElementSelectionStrategy
from KratosMultiphysics.RomApplication.RSVDT_Library import rsvdt

import numpy as np

try:
    from matplotlib import pyplot as plt
    missing_matplotlib = False
except ImportError as e:
    missing_matplotlib = True


class EmpiricalCubatureMethod(ElementSelectionStrategy):

    def __init__(self, ECM_tolerance = 1e-4, Filter_tolerance = 1e-16 ):
        super(EmpiricalCubatureMethod,self).__init__()
        self.ECM_tolerance = ECM_tolerance
        self.Filter_tolerance = Filter_tolerance
        self.Name = "EmpiricalCubature"

    def SetUp(self, ResidualSnapshots):
        super(EmpiricalCubatureMethod,self).SetUp()
        u , s  = self._ObtainBasis(ResidualSnapshots)
        self.W = np.ones(np.shape(u)[0])
        G = u[...,:] * np.ones(len(s))
        G = G.T
        G = np.vstack([ G , np.ones( np.shape(G)[1] )]  )
        b = G @ self.W
        bEXACT = b
        self.b = b
        self.G = G
        self.ExactNorm = np.linalg.norm(bEXACT)

    def Initialize(self):
        super(EmpiricalCubatureMethod,self).Initialize()
        self.Gnorm = np.sqrt(sum(np.multiply(self.G, self.G), 0))
        M = np.shape(self.G)[1]
        normB = np.linalg.norm(self.b)
        self.y = np.arange(0,M,1) # Set of candidate points (those whose associated column has low norm are removed)
        GnormNOONE = np.sqrt(sum(np.multiply(self.G[:-1,:], self.G[:-1,:]), 0))
        if self.Filter_tolerance > 0:
            TOL_REMOVE = self.Filter_tolerance * normB
            rmvpin = np.where(GnormNOONE[self.y] < TOL_REMOVE)
            self.y = np.delete(self.y,rmvpin)
        self.z = {}  # Set of intergration points
        self.mPOS = 0 # Number of nonzero weights
        self.r = self.b # residual vector
        self.m = len(self.b) # Default number of points
        self.nerror = np.linalg.norm(self.r)/normB
        self.nerrorACTUAL = self.nerror


    def Calculate(self):
        super(EmpiricalCubatureMethod,self).Calculate()

        k = 1 # number of iterations
        while self.nerrorACTUAL > self.ECM_tolerance and self.mPOS < self.m and len(self.y) != 0:

            #Step 1. Compute new point
            ObjFun = self.G[:,self.y].T @ self.r.T
            ObjFun = ObjFun.T / self.Gnorm[self.y]
            indSORT = np.argmax(ObjFun)
            i = self.y[indSORT]
            if k==1:
                alpha = np.linalg.lstsq(self.G[:, [i]], self.b)[0]
                H = 1/(self.G[:,i] @ self.G[:,i].T)
            else:
                H, alpha = self._UpdateWeightsInverse(self.G[:,self.z],H,self.G[:,i],alpha)

            #Step 3. Move i from set y to set z
            if k == 1:
                self.z = i
            else:
                self.z = np.r_[self.z,i]
            self.y = np.delete(self.y,indSORT)

            # Step 4. Find possible negative weights
            if any(alpha < 0):
                print("WARNING: NEGATIVE weight found")
                indexes_neg_weight = np.where(alpha <= 0.)[0]
                self.y = np.append(self.y, (self.z[indexes_neg_weight]).T)
                self.z = np.delete(self.z, indexes_neg_weight)
                H = self._MultiUpdateInverseHermitian(H, indexes_neg_weight)
                alpha = H @ (self.G[:, self.z].T @ self.b)
                alpha = alpha.reshape(len(alpha),1)

            #Step 6 Update the residual
            if len(alpha)==1:
                self.r = self.b - (self.G[:,self.z] * alpha)
            else:
                Aux = self.G[:,self.z] @ alpha
                self.r = np.squeeze(self.b - Aux.T)
            self.nerror = np.linalg.norm(self.r) / np.linalg.norm(self.b)  # Relative error (using r and b)
            self.nerrorACTUAL = self.nerror

            # STEP 7
            self.mPOS = np.size(self.z)
            print(f'k = {k}, m = {np.size(self.z)}, error n(res)/n(b) (%) = {self.nerror*100},  Actual error % = {self.nerrorACTUAL*100} ')

            if k == 1:
                ERROR_GLO = np.array([self.nerrorACTUAL])
                NPOINTS = np.array([np.size(self.z)])
            else:
                ERROR_GLO = np.c_[ ERROR_GLO , self.nerrorACTUAL]
                NPOINTS = np.c_[ NPOINTS , np.size(self.z)]

            k = k+1

        self.w = alpha.T * np.sqrt(self.W[self.z])

        print(f'Total number of iterations = {k}')

        if missing_matplotlib == False:
            plt.plot(NPOINTS[0], ERROR_GLO[0])
            plt.title('Element Selection Error Evolution')
            plt.xlabel('Number of elements')
            plt.ylabel('Error %')
            plt.show()

        return self.G, self.z, self.w

    def _UpdateWeightsInverse(self, A,Aast,a,xold):
        c = np.dot(A.T, a)
        d = np.dot(Aast, c).reshape(-1, 1)
        s = np.dot(a.T, a) - np.dot(c.T, d)
        aux1 = np.hstack([Aast + np.outer(d, d) / s, -d / s])
        if np.shape(-d.T / s)[1]==1:
            aux2 = np.squeeze(np.hstack([-d.T / s, 1 / s]))
        else:
            aux2 = np.hstack([np.squeeze(-d.T / s), 1 / s])
        Bast = np.vstack([aux1, aux2])
        v = np.dot(a.T, self.r) / s
        x = np.vstack([(xold - d * v), v])
        return Bast, x

    def _MultiUpdateInverseHermitian(self, invH, neg_indexes):
        neg_indexes = np.sort(neg_indexes)
        for i in range(np.size(neg_indexes)):
            neg_index = neg_indexes[i] - i
            invH = self._UpdateInverseHermitian(invH, neg_index)
        return invH

    def _UpdateInverseHermitian(self, invH, neg_index):
        if neg_index == np.shape(invH)[1]:
            aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
            invH_new = invH[:-1, :-1] - aux
        else:
            aux1 = np.hstack([invH[:, 0:neg_index], invH[:, neg_index + 1:], invH[:, neg_index].reshape(-1, 1)])
            aux2 = np.vstack([aux1[0:neg_index, :], aux1[neg_index + 1:, :], aux1[neg_index, :]])
            invH_new = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1, 0:-1]) / aux2[-1, -1]
        return invH_new

    def _ObtainBasis(self,ResidualSnapshots):
        ### Building the Snapshot matrix ####
        for i in range (len(ResidualSnapshots)):
            if i == 0:
                SnapshotMatrix = ResidualSnapshots[i]
            else:
                SnapshotMatrix = np.c_[SnapshotMatrix,ResidualSnapshots[i]]
        ### Taking the SVD ###  (randomized and truncated here)
        DATA = {}
        DATA['TypeOfSVD'] = 0
        u,s,_,_=rsvdt(SnapshotMatrix,0,0,0, DATA)
        return u, s


if __name__=='__main__':

    ECM_object = EmpiricalCubatureMethod()

