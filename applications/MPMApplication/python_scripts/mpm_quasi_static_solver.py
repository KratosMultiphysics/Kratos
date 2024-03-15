
# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.MPMApplication.mpm_implicit_dynamic_solver import MPMImplicitDynamicSolver

def CreateSolver(model, custom_settings):
    return MPMQuasiStaticSolver(model, custom_settings)

class MPMQuasiStaticSolver(MPMImplicitDynamicSolver):

    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MPMQuasiStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMQuasiStaticSolver]:: ", "Construction is finished.")

    ### Protected functions ###

    def _IsDynamic(self):
        return False
