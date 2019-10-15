from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return MPMStaticSolver(model, custom_settings)

class MPMStaticSolver(MPMSolver):
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MPMStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMStaticSolver]:: ", "Construction is finished.")

    def _CreateSolutionScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()