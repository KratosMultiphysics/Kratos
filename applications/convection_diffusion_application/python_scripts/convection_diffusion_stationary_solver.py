from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
import convection_diffusion_base_solver

def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionStationarySolver(main_model_part, custom_settings)

class ConvectionDiffusionStationarySolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """The stationary class for convection-diffusion solvers.

    Public member variables:
    stationary_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.stationary_settings = KratosMultiphysics.Parameters("""
        {
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.stationary_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ConvectionDiffusionStationarySolver, self).__init__(main_model_part, custom_settings)
        self.print_on_rank_zero("::[ConvectionDiffusionStationarySolver]:: ", "Construction finished")

    def GetMinimumBufferSize(self):
        if (self.settings["element_replace_settings"]["element_name"].GetString() == "EulerianConvDiff"):
            return 2
        else:
            return 1

    #### Private functions ####

    def _create_solution_scheme(self):
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[ConvectionDiffusionApplication.THETA] = 1.0
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = 0.0
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme
