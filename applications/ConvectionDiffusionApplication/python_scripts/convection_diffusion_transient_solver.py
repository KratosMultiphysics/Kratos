from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_base_solver

def CreateSolver(model, custom_settings):
    return ConvectionDiffusionTransientSolver(model, custom_settings)

class ConvectionDiffusionTransientSolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """The transient class for convection-diffusion solvers.

    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """

    def __init__(self, model, custom_settings):
        # Construct the base solver and validate the remaining settings in the base class
        super(ConvectionDiffusionTransientSolver, self).__init__(model, custom_settings)

        # Overwrite the base solver minimum buffer size
        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "transient_parameters" : {
                "dynamic_tau": 1.0,
                "theta"    : 0.5
            }
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    #### Private functions ####
    def _create_solution_scheme(self):
        # Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME_INTEGRATION_THETA] = self.settings["transient_parameters"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.settings["transient_parameters"]["dynamic_tau"].GetDouble()
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme
