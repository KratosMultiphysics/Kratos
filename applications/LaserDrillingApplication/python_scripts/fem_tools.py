import numpy as np
from sympy import *

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
        result = 0.0
        for i, U in enumerate(nodal_values):
            result += U * self.N(i, x)
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
        NDesc = self.NDesc
        NAsc = self.NAsc

        x = symbols('x')
        for i in range(0, n+1):
            if i == 0:
                integrand = NDesc(i, x)**2 * x
                A_diag[i] = integrate(integrand, (x, X[0], X[1]))
                integrand = NDesc(i, x) * NAsc(i+1, x) * x
                A_diag_top[i] = integrate(integrand, (x, X[i], X[i+1]))
                A_diag_bot[i] = A_diag_top[i]
                A_diag[i] = 2 * np.pi * (A_diag[i] + A_diag_bot[i])
            elif i < n:
                integrand1 = NAsc(i, x)**2 * x
                integrand2 = NDesc(i, x)**2 * x
                A_diag[i] = integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
                integrand = NDesc(i, x) * NAsc(i+1, x) * x
                A_diag_top[i] = integrate(integrand, (x, X[i], X[i+1]))
                A_diag_bot[i] = A_diag_top[i]
                A_diag[i] = 2 * np.pi * (A_diag[i] + A_diag_bot[i] + A_diag_top[i-1])
            else:
                integrand = NAsc(i, x)**2 * x
                A_diag[i] = integrate(integrand, (x, X[n-1], X[n]))
                A_diag[i] = 2 * np.pi * (A_diag[i] + A_diag_top[i-1])

    def FillUpDeltasRHS(self, evap_element_centers, support_elements, evap_enthalpies):
        n = self.n_elements
        N = self.N        
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
        NDesc = self.NDesc
        NAsc = self.NAsc

        total_energy = 0.0

        for i, u_value in enumerate(nodal_values):
            x = symbols('x')
            U = u_value[0]

            if i == 0:
                integrand = U * NDesc(i, x) * x
                total_energy += integrate(integrand, (x, X[i], X[i+1]))
            elif i == n:
                integrand = U * NAsc(i, x) * x
                total_energy += integrate(integrand, (x, X[i-1], X[i]))
            else:
                integrand1 = U * NAsc(i, x) * x
                integrand2 = U * NDesc(i, x) * x
                total_energy += integrate(integrand1, (x, X[i-1], X[i])) + integrate(integrand2, (x, X[i], X[i+1]))
        return float(2.0 * np.pi * total_energy)

    def InterpolateFunctionAndNormalize(self, f): #, normalization_value):
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
            '''for k, b in enumerate(self.b):
                self.solution[k] = b / self.A_diag[k]'''
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
