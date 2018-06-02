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
    return ConvectionDiffusionTransientSolver(main_model_part, custom_settings)

class ConvectionDiffusionTransientSolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """The transient class for convection-diffusion solvers.

    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.transient_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "ConvectionDiffusionTransientSolver",
            "dynamic_tau": 1.0,
            "time_stepping" : {
                "theta"    : 1.0
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.transient_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(ConvectionDiffusionTransientSolver, self).__init__(main_model_part, custom_settings)
        self.print_on_rank_zero("::[ConvectionDiffusionTransientSolver]:: ", "Construction finished")

    #### Private functions ####

    def _create_solution_scheme(self):
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[ConvectionDiffusionApplication.THETA] = self.transient_settings["time_stepping"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.transient_settings["dynamic_tau"].GetDouble()
        mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return mechanical_scheme
