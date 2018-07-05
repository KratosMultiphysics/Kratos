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
     
    def __init__(self, thermal_model_part, custom_settings):
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


        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "laplacian_solver_newformat",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 10,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm" : false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                    "solver_type": "AMGCL",
                    "max_iteration": 100,
                    "tolerance": 1e-6
            }
        }""")
        # Validate the remaining settings in the base class.
        print(custom_settings)
        # Construct the base solver.
        #construct the linear solvers

        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(default_settings["linear_solver_settings"])
         
        #super(ConvectionDiffusionTransientSolver, self).__init__(main_model_part, custom_settings)
        
        self.print_on_rank_zero("::[ConvectionDiffusionTransientSolver]:: ", "Construction finished")

    def GetMinimumBufferSize(self):
        kjljkljklkjkljlkjljl
        return 2

    #### Private functions ####

    def _create_solution_scheme(self):
        oifgpipfgoipofgipgfoiphigfp 
        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.GetComputingModelPart().ProcessInfo[ConvectionDiffusionApplication.THETA] = self.transient_settings["transient_parameters"]["theta"].GetDouble()
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DYNAMIC_TAU] = self.transient_settings["transient_parameters"]["dynamic_tau"].GetDouble()
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return convection_diffusion_scheme



