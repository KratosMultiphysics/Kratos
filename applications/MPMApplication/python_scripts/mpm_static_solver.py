
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM


# Import applications and dependencies
import KratosMultiphysics.MPMApplication as KratosMPM

# Importing the base class
from KratosMultiphysics.MPMApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return MPMStaticSolver(model, custom_settings)


class MPMStaticSolver(MPMSolver):
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MPMStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMStaticSolver]:: ", "Construction is finished.")

    def AddVariables(self):
        super(MPMStaticSolver, self).AddVariables()
        self._AddDynamicVariables(self.grid_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MPMStaticSolver]:: ", "Variables are all added.")

    def _CreateSolutionScheme(self):

        return KratosMPM.MPMResidualBasedIncrementalUpdateStaticScheme()
