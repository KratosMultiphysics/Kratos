from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ConvectionDiffusionApplication")

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

# Import base class file
import convection_diffusion_base_solver

def CreateSolver(model, custom_settings):
    
    return ConvectionDiffusionTransientSolver(model, custom_settings)

class ConvectionDiffusionTransientSolver(convection_diffusion_base_solver.ConvectionDiffusionBaseSolver):
    """The transient class for convection-diffusion solvers.

    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.

    See convection_diffusion_base_solver.py for more information.
    """
     
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        


        self.transient_settings = KratosMultiphysics.Parameters("""
        {
            "transient_parameters" : {
                "dynamic_tau": 1.0,
                "theta"    : 0.5
            }
        }
        """)

              
        self.validate_and_transfer_matching_settings(custom_settings, self.transient_settings)

        # Validate the remaining settings in the base class.
        # Construct the base solver.
        #construct the linear solvers

        super(ConvectionDiffusionTransientSolver, self).__init__(model, custom_settings)
        
        self.print_on_rank_zero("::[ConvectionDiffusionTransientSolver]:: ", "Construction finished")

    def GetMinimumBufferSize(self):
        
        return 2

    #### Private functions ####

    def _create_solution_scheme(self):
        
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[ConvectionDiffusionApplication.THETA] = self.transient_settings["transient_parameters"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.transient_settings["transient_parameters"]["dynamic_tau"].GetDouble()
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme



