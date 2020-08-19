# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.ParticleMechanicsApplication.trilinos_mpm_solver import TrilinosMPMSolver

def CreateSolver(model, custom_settings):
    return TrilinosMPMStaticSolver(model, custom_settings)

class TrilinosMPMStaticSolver(TrilinosMPMSolver):
    """The trilinos mpm static solver.

    For more information see:
    mpm_solver.py
    trilinos_mpm_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(TrilinosMPMStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosMPMStaticSolver]:: ", "Construction finished")

    def _CreateSolutionScheme(self):
        return TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()