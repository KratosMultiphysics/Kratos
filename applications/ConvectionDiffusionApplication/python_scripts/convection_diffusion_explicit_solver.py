# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_base_solver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionExplicitSolver(model, custom_settings)

class ConvectionDiffusionExplicitSolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """
    The explicit class for convection-diffusion solvers.
    See convection_diffusion_base_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super(ConvectionDiffusionExplicitSolver, self).__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "explicit_parameters" : {
                "dynamic_tau": 1.0
            }
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    #### Private functions ####
    def _create_solution_scheme(self):
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.settings["explicit_parameters"]["dynamic_tau"].GetDouble()
