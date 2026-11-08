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
    
    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()