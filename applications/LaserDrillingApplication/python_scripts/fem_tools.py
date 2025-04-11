import numpy as np
from sympy import *

#TODO: Why is all this symbolic?
class SurfaceFEMProjector:
    def __init__(self, n_elements, R_far, sparse_option):
        n = n_elements
        self.n_elements = n # number of elements
        self.R_far = R_far
        self.sparse_option = sparse_option
        self.X = np.linspace(0.0, R_far, n+1)
        if not self.sparse_option:
            self.A = np.zeros((n+1,n+1))
        else:
            self.A_diag = np.zeros((n+1,1))
            self.A_diag_top = np.zeros((n,1))
            self.A_diag_bot = np.zeros((n,1))
        self.b = np.zeros((n+1,1))        
        self.solution = np.zeros((n+1,1))

    def NAsc(self, i, x):
        X = self.X
        return (x - X[i - 1]) / (X[i] - X[i - 1])

    def NDesc(self, i, x):
        X = self.X
        return (X[i + 1] - x) / (X[i + 1] - X[i])

    def N(self, i, x):
        NDesc = self.NDesc
        NAsc = self.NAsc        
        n = self.n_elements
        X = self.X

        if i == 0:
            if np.any(x >= X[i]) and np.any(x <= X[i + 1]):
                return NDesc(i, x)
            else:
                return 0.0
        if i == n:
            if np.any(x >= X[i - 1]) and np.any(x <= X[i]):
                return NAsc(i, x)
            else:
                return 0.0
        if np.any(x >= X[i - 1]) and np.any(x <= X[i]):
            return NAsc(i, x)
        elif np.any(x >= X[i]) and np.any(x <= X[i + 1]):
            return NDesc(i,x)
        else:
            return 0.0

    def EvaluateFEMFunction(self, nodal_values, x):
        if not self.sparse_option:
            result = 0.0
            for i, U in enumerate(nodal_values):
                result += U * self.N(i, x)
            return result
        else:
            X = self.X
            i_val = 0
            for i, _ in enumerate(X):
                if x >= X[i] and x <= X[i+1]:
                    i_val = i
                    break
            result = nodal_values[i_val] * self.N(i_val, x) + nodal_values[i_val+1] * self.N(i_val+1, x)
            return result

    def FillUpMassMatrix(self):
        X = self.X
        n = self.n_elements
        A = self.A
        NDesc = self.NDesc
        NAsc = self.NAsc

        x = symbols('x')
        for i in range(0, n+1):
            for j in range(0, n+1):
                if i == j:
                    if i == 0:
                        integrand = NDesc(i, x) * NDesc(j, x) * x
                        A[i, j] = integrate(integrand, (x, X[0], X[1]))
                    elif i == n:
                        integrand = NAsc(i, x) * NAsc(j, x) * x
                        A[i, j] = integrate(integrand, (x, X[n-1], X[n]))
                    else:
                        integrand1 = NAsc(i, x) * NAsc(i, x) * x
                        integrand2 = NDesc(i, x) * NDesc(i, x) * x
                        A[i, j] = integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
                elif j - i == 1:
                    integrand = NDesc(i, x) * NAsc(j, x) * x
                    A[i, j] = integrate(integrand, (x, X[i], X[j]))
        for i in range(0, n+1):
            for j in range(0, n+1):
                if i - j == 1:
                    A[i, j] = A[j, i]

        # We will always lump the matrix because it behaves much better numerically
        for i in range(0, n+1):
            sum_i = 0.0
            for j in range(0, n+1):
                sum_i += A[i, j]
            A[i, i] = sum_i
            for j in range(0, n+1):
                if i != j:
                    A[i, j] = 0.0
        A *= 2 * np.pi

    def FillUpSparseMassMatrix(self):
        X = self.X
        n = self.n_elements
        A_diag = self.A_diag
        A_diag_bot = self.A_diag_bot
        A_diag_top = self.A_diag_top

        A_diag[0] = (1.0 / (X[1]-X[0])**2) * (X[1]*X[1]*X[1]*X[1] * 0.083333333333333 - 0.5 * X[1]**2 * X[0]**2 - 0.25 * X[0]*X[0]*X[0]*X[0] + 0.66666666666666 * X[1] * X[0]*X[0]*X[0])
        A_diag_top[0] = (1.0 / (X[1]-X[0])**2) * (0.083333333333333*X[1]*X[1]*X[1]*X[1] - 0.1666666666666667 * X[1]*X[1]*X[1]*X[0] + 0.166666666666666*X[1]*X[0]*X[0]*X[0] - 0.083333333333333*X[0]*X[0]*X[0]*X[0])
        A_diag_bot[0] = A_diag_top[0]
        A_diag[0] = 2 * np.pi * (A_diag[0] + A_diag_bot[0])
        for i in range(1, n):
            integrand1 = (1.0 / (X[i]-X[i-1])**2) * (0.25 * X[i]*X[i]*X[i]*X[i] + 0.5 * X[i-1]**2 * X[i]**2 - 0.66666666666666 * X[i-1]*X[i]*X[i]*X[i] - X[i-1]*X[i-1]*X[i-1]*X[i-1] * 0.083333333333333)
            integrand2 = (1.0 / (X[i+1]-X[i])**2) * (X[i+1]*X[i+1]*X[i+1]*X[i+1] * 0.083333333333333 - 0.5 * X[i+1]**2 * X[i]**2 - 0.25 * X[i]*X[i]*X[i]*X[i] + 0.66666666666666 * X[i+1] * X[i]*X[i]*X[i])
            A_diag[i] = integrand1 + integrand2
            A_diag_top[i] = (1.0 / (X[i+1]-X[i])**2) * (0.083333333333333*X[i+1]*X[i+1]*X[i+1]*X[i+1] - 0.1666666666666667 * X[i+1]*X[i+1]*X[i+1]*X[i] + 0.166666666666666*X[i+1]*X[i]*X[i]*X[i] - 0.083333333333333*X[i]*X[i]*X[i]*X[i])
            A_diag_bot[i] = A_diag_top[i]
            A_diag[i] = 2 * np.pi * (A_diag[i] + A_diag_bot[i] + A_diag_top[i-1])
        A_diag[n] = (1.0 / (X[n]-X[n-1])**2) * (0.25 * X[n]*X[n]*X[n]*X[n] + 0.5 * X[n-1]**2 * X[n]**2 - 0.66666666666666 * X[n-1]*X[n]*X[n]*X[n] - X[n-1]*X[n-1]*X[n-1]*X[n-1] * 0.083333333333333)
        A_diag[n] = 2 * np.pi * (A_diag[n] + A_diag_top[n-1])

    def FillUpDeltasRHS(self, evap_element_centers, support_elements, evap_enthalpies):
        n = self.n_elements
        N = self.N
        self.b = np.zeros((n+1,1))
        b = self.b

        for i in range(n + 1):
            for k in support_elements[i]:
                b[i] += N(i, evap_element_centers[k]) * evap_enthalpies[k]

    def FillUpFEMRHS(self, j):
        X = self.X
        n = self.n_elements
        NDesc = self.NDesc
        NAsc = self.NAsc
        b = self.b

        x = symbols('x')
        for i in range(0, n+1):
            if i == j and i == 0:
                integrand = NDesc(i, x) * NDesc(j, x)
                b[i] = integrate(integrand, (x, X[0], X[1]))
            elif i == j and i == n:
                integrand = NAsc(i, x) * NAsc(j, x)
                b[i] = integrate(integrand, (x, X[n-1], X[n]))
            elif i == j:
                integrand1 = NAsc(i, x)**2
                integrand2 = NDesc(i, x)**2
                b[i] = integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
            elif j - i == 1:
                integrand = NDesc(i, x) * NAsc(j, x)
                b[i] = integrate(integrand, (x, X[i], X[j]))
            elif i - j == 1:
                integrand = NDesc(j, x) * NAsc(i, x)
                b[i] = integrate(integrand, (x, X[j], X[i]))                

    def CalculateEnergyOfFEMFunction(self, nodal_values):
        X = self.X
        n = self.n_elements
        #NDesc = self.NDesc
        #NAsc = self.NAsc
        total_energy = 0.0

        for i, u_value in enumerate(nodal_values):
            #x = symbols('x')
            U = u_value[0]

            if i == 0:
                #integrand = U * NDesc(i, x) * x
                #total_energy += integrate(integrand, (x, X[i], X[i+1]))
                total_energy += U * 1.0 / (X[i+1]-X[i]) * (0.1666666666666667 * X[i+1]*X[i+1]*X[i+1] - 0.5*X[i+1]*X[i]**2 + 0.3333333333333333333 * X[i]*X[i]*X[i])
            elif i < n:
                #integrand1 = U * NAsc(i, x) * x
                #integrand2 = U * NDesc(i, x) * x
                integrand1 = U * 1.0 / (X[i]-X[i-1]) * (0.3333333333333333333 * X[i]*X[i]*X[i] - 0.5*X[i-1]*X[i]**2 + 0.1666666666666667 * X[i-1]*X[i-1]*X[i-1])
                integrand2 = U * 1.0 / (X[i+1]-X[i]) * (0.1666666666666667 * X[i+1]*X[i+1]*X[i+1] - 0.5*X[i+1]*X[i]**2 + 0.3333333333333333333 * X[i]*X[i]*X[i])
                #total_energy += integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
                total_energy += integrand1 + integrand2
            else:
                #integrand = U * NAsc(i, x) * x
                #total_energy += integrate(integrand, (x, X[i-1], X[i]))
                total_energy += U * 1.0 / (X[i]-X[i-1]) * (0.3333333333333333333 * X[i]*X[i]*X[i] - 0.5*X[i-1]*X[i]**2 + 0.1666666666666667 * X[i-1]*X[i-1]*X[i-1])
                
        return float(2.0 * np.pi * total_energy)

    def InterpolateFunctionAndNormalize(self, f): #, normalization_value):
        # TODO: Check this function to see if it is normalizing the argument or not or what
        X = self.X
        n = self.n_elements

        nodal_values = np.zeros((n+1,1))
        for i, x in enumerate(X):
            nodal_values[[i]] = f(x)

        #obtained_energy = self.CalculateEnergyOfFEMFunction(nodal_values)
        #nodal_values *= normalization_value / obtained_energy
        return nodal_values     

    def Project(self):
        if not self.sparse_option:
            self.solution = np.linalg.solve(self.A, self.b)
        else:
            self.solution = self.b / self.A_diag
        return self.solution 

    def q(self, r):
        return 100 * np.exp(-r)

    def AssignDeltasToTestFunctionSupports(self, radii, support_elements):    
        X = self.X
        n = self.n_elements
        for j in range(n+1):
            for i_rad, r in enumerate(radii):
                if j == 0 and r <= X[1]:
                    support_elements[j].append(i_rad)
                elif j == n and r >= X[n-1]:
                    support_elements[j].append(i_rad)
                elif r >= X[j-1] and r <= X[j+1]:
                    support_elements[j].append(i_rad)
