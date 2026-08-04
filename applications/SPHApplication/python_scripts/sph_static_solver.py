import numpy as np

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SPHApplication as SPH

from KratosMultiphysics.SPHApplication.sph_solver import SPHSolver

def CreateSolver(model, custom_settings):
    return StaticSPHSolver(model, custom_settings)

class StaticSPHSolver(SPHSolver):

    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticSPHSolver]:: ", "Construction finished")
    
    def AddVariables(self):
        super().AddVariables()

    def AddDofs(self):
        super().AddDofs()
        KratosMultiphysics.Logger.PrintInfo("::[StaticSPHSolver]::", "AddDofs called")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.ExposeSystemMatrix()
    
    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    
    ## DEBUG METHODS
    
    def ExposeSystemMatrix(self):
        A = self._GetSolutionStrategy().GetSystemMatrix() 
        vec = self._GetSolutionStrategy().GetSystemVector() 
        sol = self._GetSolutionStrategy().GetSolutionVector()

        n = A.Size1()
        m = A.Size2()
        B = np.zeros((n,m))
        vec1 = np.zeros(n)
        sol1 = np.zeros(n)

        for i in range(n):
            vec1[i] = vec[i]
            sol1[i] = sol[i]
            for j in range(m):
                B[i ,j] = A[i, j]
        
        self.system_matrix = A
        self.right_hand_side = vec
        self.solution = sol