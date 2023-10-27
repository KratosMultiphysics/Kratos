import numpy as np
from KratosMultiphysics import Logger

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

    def __init__(
        self,
        ECM_tolerance = 0,
        Filter_tolerance = 0,
        Plotting = False,
        MaximumNumberUnsuccesfulIterations = 100
    ):
        """
        Constructor setting up the parameters for the Element Selection Strategy
            ECM_tolerance: approximation tolerance for the element selection algorithm
            Filter_tolerance: parameter limiting the number of candidate points (elements) to those above this tolerance
            Plotting: whether to plot the error evolution of the element selection algorithm
        """
        self.ECM_tolerance = ECM_tolerance
        self.Filter_tolerance = Filter_tolerance
        self.Name = "EmpiricalCubature"
        self.Plotting = Plotting
        self.MaximumNumberUnsuccesfulIterations = MaximumNumberUnsuccesfulIterations

    def SetUp(
        self,
        ResidualsBasis,
        InitialCandidatesSet = None,
        constrain_sum_of_weights=True,
        constrain_conditions = False,
        number_of_conditions = 0
    ):
        """
        Method for setting up the element selection
        input:  - ResidualsBasis: numpy array containing a basis to the residuals projected
                - constrain_sum_of_weights: enable the user to constrain weights to be the sum of the number of entities.
                - constrain_conditions: enable the user to enforce weights to consider conditions (for specific boundary conditions).
        """
        self.W = np.ones(np.shape(ResidualsBasis)[0])
        self.G = ResidualsBasis.T
        self.y = InitialCandidatesSet
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
        self.UnsuccesfulIterations = 0

    def Initialize(self):
        """
        Method performing calculations required before launching the Calculate method
        """
        self.GnormNOONE = np.linalg.norm(self.G[:self.add_constrain_count,:], axis = 0)
        M = np.shape(self.G)[1]
        normB = np.linalg.norm(self.b)

        if self.y is None:
            self.y = np.arange(0,M,1) # Set of candidate points (those whose associated column has low norm are removed)

            if self.Filter_tolerance > 0:
                TOL_REMOVE = self.Filter_tolerance * normB
                rmvpin = np.where(self.GnormNOONE[self.y] < TOL_REMOVE)
                #self.y_complement = self.y[rmvpin]
                self.y = np.delete(self.y,rmvpin)
        else:
            self.y_complement = np.arange(0, M, 1)  # Initialize complement with all points
            self.y_complement = np.delete(self.y_complement, self.y)  # Remove candidates from complement

            if self.Filter_tolerance > 0:
                TOL_REMOVE = self.Filter_tolerance * normB  # Compute removal tolerance

                # Filter out low-norm columns from complement
                rmvpin_complement = np.where(self.GnormNOONE[self.y_complement] < TOL_REMOVE)
                self.y_complement = np.delete(self.y_complement, rmvpin_complement)

                # Filter out low-norm columns from candidates
                rmvpin = np.where(self.GnormNOONE[self.y] < TOL_REMOVE)
                removed_count = np.size(rmvpin)
                self.y = np.delete(self.y, rmvpin)

                # Warning if some candidates were removed
                if removed_count > 0:
                    Logger.PrintWarning("EmpiricalCubatureMethod", f"Some of the candidates were removed ({removed_count} removed). To include all candidates (with 0 weights in the HROM model part) for visualization and projection, consider using 'include_elements_model_parts_list' and 'include_conditions_model_parts_list' in the 'hrom_settings'.")

                # Warning if all candidates were removed
                if np.size(self.y) == 0:
                    Logger.PrintWarning("EmpiricalCubatureMethod", "All candidates were removed because they have no contribution to the residual. To include them all (with 0 weights in the HROM model part) for visualization and projection, use 'include_elements_model_parts_list' and 'include_conditions_model_parts_list' in the 'hrom_settings'.")
                    self.y = self.y_complement  # Set candidates to complement

        self.z = {}  # Set of intergration points
        self.mPOS = 0 # Number of nonzero weights
        self.r = self.b.copy() # residual vector
        self.m = len(self.b) # Default number of points
        self.nerror = np.linalg.norm(self.r)/normB
        self.nerrorACTUAL = self.nerror

    def Run(self):
        self.Initialize()
        self.Calculate()

    def expand_candidates_with_complement(self):
        self.y = np.r_[self.y,self.y_complement]
        print('expanding set to include the complement...')
        ExpandedSetFlag = True
        return ExpandedSetFlag

    def Calculate(self):
        """
        Method launching the element selection algorithm to find a set of elements: self.z, and wiegths: self.w
        """
        MaximumLengthZ = 0
        ExpandedSetFlag = False
        k = 1 # number of iterations
        self.success = True
        while self.nerrorACTUAL > self.ECM_tolerance and self.mPOS < self.m and np.size(self.y) != 0:

            if  self.UnsuccesfulIterations >  self.MaximumNumberUnsuccesfulIterations and not ExpandedSetFlag and hasattr(self, 'y_complement'):
                ExpandedSetFlag = self.expand_candidates_with_complement()

            #Step 1. Compute new point
            if np.size(self.y)==1:
                #candidate set consists of a single element
                indSORT = 0
                i = int(self.y)
            else:
                ObjFun = self.G[:,self.y].T @ self.r.T
                ObjFun = ObjFun.T #/ self.GnormNOONE[self.y]
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

            #self.y = np.delete(self.y,indSORT)
            if np.size(self.y)==1:
                if hasattr(self, 'y_complement'):
                    self.expand_candidates_with_complement()
                    self.y = np.delete(self.y,indSORT)
                else:
                    self.success = False
                    break
            else:
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

            if np.size(self.z) > MaximumLengthZ :
                self.UnsuccesfulIterations = 0
            else:
                self.UnsuccesfulIterations += 1

            #Step 6 Update the residual
            if np.size(alpha)==1:
                self.r = self.b.reshape(-1,1) - (self.G[:,self.z] * alpha).reshape(-1,1)
                self.r = np.squeeze(self.r)
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

            MaximumLengthZ = max(MaximumLengthZ, np.size(self.z))
            k = k+1

            if k-MaximumLengthZ>1000 and ExpandedSetFlag:
                """
                this means using the initial candidate set, it was impossible to obtain a set of positive weights.
                Try again without constraints!!!
                TODO: incorporate this into greater workflow
                """
                self.success = False
                break

        self.w = alpha.T * np.sqrt(self.W[self.z]) #TODO FIXME cope with weights vectors different from 1

        print(f'Total number of iterations = {k}')

        if missing_matplotlib == False and self.Plotting == True:
            plt.plot(NPOINTS[0], ERROR_GLO[0])
            plt.title('Element Selection Error Evolution')
            plt.xlabel('Number of elements')
            plt.ylabel('Error %')
            plt.show()

    def _UpdateWeightsInverse(self, A,Aast,a,xold):
        """
        Method for the quick update of weights (self.w), whenever a negative weight is found
        """
        c = np.dot(A.T, a)
        d = np.dot(Aast, c).reshape(-1, 1)
        s = np.dot(a.T, a) - np.dot(c.T, d)
        aux1 = np.hstack([Aast + np.outer(d, d) / s, -d / s])
        if np.shape(-d.T / s)[1]==1:
            s = s.reshape(1,-1)
            aux2 = np.squeeze(np.hstack([-d.T / s, 1 / s]))
        else:
            aux2 = np.hstack([np.squeeze(-d.T / s), 1 / s])
        Bast = np.vstack([aux1, aux2])
        v = np.dot(a.T, self.r) / s
        x = np.vstack([(xold - d * v), v])
        return Bast, x

    def _MultiUpdateInverseHermitian(self, invH, neg_indexes):
        """
        Method for the quick update of weights (self.w), whenever a negative weight is found
        """
        neg_indexes = np.sort(neg_indexes)
        for i in range(np.size(neg_indexes)):
            neg_index = neg_indexes[i] - i
            invH = self._UpdateInverseHermitian(invH, neg_index)
        return invH

    def _UpdateInverseHermitian(self, invH, neg_index):
        """
        Method for the quick update of weights (self.w), whenever a negative weight is found
        """
        if neg_index == np.shape(invH)[1]:
            aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
            invH_new = invH[:-1, :-1] - aux
        else:
            aux1 = np.hstack([invH[:, 0:neg_index], invH[:, neg_index + 1:], invH[:, neg_index].reshape(-1, 1)])
            aux2 = np.vstack([aux1[0:neg_index, :], aux1[neg_index + 1:, :], aux1[neg_index, :]])
            invH_new = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1, 0:-1]) / aux2[-1, -1]
        return invH_new






