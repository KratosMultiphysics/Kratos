import numpy as np


try:
    from matplotlib import pyplot as plt
    missing_matplotlib = False
except ImportError as e:
    missing_matplotlib = True


class EmpiricalCubatureMethod():
    """
    This class selects a subset of elements and corresponding positive weights necessary for the construction of a hyper-reduced order model
    Reference: Hernandez 2020. "A multiscale method for periodic structures using domain decomposition and ECM-hyperreduction"
    """


    """
    Constructor setting up the parameters for the Element Selection Strategy
        ECM_tolerance: approximation tolerance for the element selection algorithm
        Filter_tolerance: parameter limiting the number of candidate points (elements) to those above this tolerance
        Plotting: whether to plot the error evolution of the element selection algorithm
    """
    def __init__(self, ECM_tolerance = 0, Filter_tolerance = 1e-16, Plotting = False):
        self.ECM_tolerance = ECM_tolerance
        self.Filter_tolerance = Filter_tolerance
        self.Name = "EmpiricalCubature"
        self.Plotting = Plotting



    """
    Method for setting up the element selection
    input:  - ResidualsBasis: numpy array containing a basis to the residuals projected
            - constrain_sum_of_weights: enable the user to constrain weights to be the sum of the number of entities.
            - constrain_conditions: enable the user to enforce weights to consider conditions (for specific boundary conditions).
    """
    def SetUp(self, ResidualsBasis, constrain_sum_of_weights=True, constrain_conditions = False, number_of_conditions = 0):

        self.W = np.ones(np.shape(ResidualsBasis)[0])
        self.G = ResidualsBasis.T
        self.add_constrain_count = None
        total_number_of_entities = np.shape(self.G)[1]
        elements_constraint = np.ones(total_number_of_entities)
        conditions_begin = total_number_of_entities - number_of_conditions
        elements_constraint[conditions_begin:] = 0

        if constrain_sum_of_weights and not constrain_conditions:
            """
            -This is necessary in case the sum of the columns of self.G equals the 0 vector,to avoid the trivial solution
            -It is enforcing that the sum of the weights equals the number of columns in self.G (total number of elements)
            """
            projection_of_constant_vector_elements = elements_constraint - self.G.T@( self.G @ elements_constraint)
            projection_of_constant_vector_elements/= np.linalg.norm(projection_of_constant_vector_elements)
            self.G = np.vstack([ self.G , projection_of_constant_vector_elements] )
            self.add_constrain_count = -1
        elif constrain_sum_of_weights and constrain_conditions:#Only for models which contains conditions
            projection_of_constant_vector_elements = elements_constraint - self.G.T@( self.G @ elements_constraint)
            projection_of_constant_vector_elements/= np.linalg.norm(projection_of_constant_vector_elements)
            self.G = np.vstack([ self.G , projection_of_constant_vector_elements] )
            # # # # # # # # #
            conditions_constraint = np.ones(total_number_of_entities)
            conditions_constraint[:conditions_begin] = 0
            projection_of_constant_vector_conditions = conditions_constraint - self.G.T@( self.G @ conditions_constraint)
            projection_of_constant_vector_conditions/= np.linalg.norm(projection_of_constant_vector_conditions)
            self.G = np.vstack([ self.G , projection_of_constant_vector_conditions ] )
            self.add_constrain_count = -2
        self.b = self.G @ self.W


    """
    Method performing calculations required before launching the Calculate method
    """
    def Initialize(self):
        self.Gnorm = np.sqrt(sum(np.multiply(self.G, self.G), 0))
        M = np.shape(self.G)[1]
        normB = np.linalg.norm(self.b)
        self.y = np.arange(0,M,1) # Set of candidate points (those whose associated column has low norm are removed)
        GnormNOONE = np.sqrt(sum(np.multiply(self.G[:self.add_constrain_count,:], self.G[:self.add_constrain_count,:]), 0))
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


    def Run(self):
        self.Initialize()
        self.Calculate()


    """
    Method launching the element selection algorithm to find a set of elements: self.z, and wiegths: self.w
    """
    def Calculate(self):

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

        self.w = alpha.T * np.sqrt(self.W[self.z]) #TODO FIXME cope with weights vectors different from 1

        print(f'Total number of iterations = {k}')

        if missing_matplotlib == False and self.Plotting == True:
            plt.plot(NPOINTS[0], ERROR_GLO[0])
            plt.title('Element Selection Error Evolution')
            plt.xlabel('Number of elements')
            plt.ylabel('Error %')
            plt.show()



    """
    Method for the quick update of weights (self.w), whenever a negative weight is found
    """
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



    """
    Method for the quick update of weights (self.w), whenever a negative weight is found
    """
    def _MultiUpdateInverseHermitian(self, invH, neg_indexes):
        neg_indexes = np.sort(neg_indexes)
        for i in range(np.size(neg_indexes)):
            neg_index = neg_indexes[i] - i
            invH = self._UpdateInverseHermitian(invH, neg_index)
        return invH



    """
    Method for the quick update of weights (self.w), whenever a negative weight is found
    """
    def _UpdateInverseHermitian(self, invH, neg_index):
        if neg_index == np.shape(invH)[1]:
            aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
            invH_new = invH[:-1, :-1] - aux
        else:
            aux1 = np.hstack([invH[:, 0:neg_index], invH[:, neg_index + 1:], invH[:, neg_index].reshape(-1, 1)])
            aux2 = np.vstack([aux1[0:neg_index, :], aux1[neg_index + 1:, :], aux1[neg_index, :]])
            invH_new = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1, 0:-1]) / aux2[-1, -1]
        return invH_new






