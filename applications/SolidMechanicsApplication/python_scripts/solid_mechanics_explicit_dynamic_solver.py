from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics explicit dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):

        # Set defaults and validate custom settings.
        ##TODO : solving_strategy_settings must be time_integration_settings (GiD interface changes needed)
        explicit_solver_settings = KratosMultiphysics.Parameters("""
        {
           "solving_strategy_settings":{
               "time_step_prediction_level": 0,
               "max_delta_time": 1.0e-5,
               "fraction_delta_time": 0.9,
               "rayleigh_damping": false,
               "rayleigh_alpha": 0.0,
               "rayleigh_beta" : 0.0
           }

        }
        """)

        # Validate and transfer settings
        if( custom_settings.Has("solving_strategy_settings") ):
            self._validate_and_transfer_matching_settings(custom_settings["solving_strategy_settings"], explicit_solver_settings["solving_strategy_settings"])
        self.explicit_solver_settings = explicit_solver_settings["solving_strategy_settings"]

        # Validate the remaining settings in the base class.
        if not custom_settings["solving_strategy_settings"].Has("integration_method"): # Override defaults in the base class.
            custom_settings["solving_strategy_settings"].AddEmptyValue("integration_method")
            custom_settings["solving_strategy_settings"]["integration_method"].SetString("CentralDifferences")

        # Construct the base solver.
        super(ExplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Explicit_Scheme]:: "+self.time_integration_settings["integration_method"].GetString()+" Scheme Ready")


    def GetVariables(self):

        nodal_variables = super(ExplicitMechanicalSolver, self).GetVariables()

        time_integration = self.time_integration_settings["time_integration"].GetString()
        # Add specific variables for the explicit time integration scheme
        if(time_integration == "Explicit"):
            nodal_variables = nodal_variables + ['NODAL_MASS','MIDDLE_VELOCITY','FORCE_RESIDUAL']
            # Add specific variables for the explicit time integration scheme in rotations
            #if(self.time_integration_settings["integration_method"].GetString() == "ExplicitRotation"):
            #    nodal_variables = nodal_variables + ['INERTIA_DYADIC','MOMENT_RESIDUAL','POSITION_MOMENTUM','ROTATION_MOMENTUM', 'RESIDUAL_LYAPUNOV', 'TANGENT_LYAPUNOV']

        return nodal_variables

    #### Specific internal functions ####
    def _create_solution_scheme(self):

        integration_method   = self.time_integration_settings["integration_method"].GetString()

        if(integration_method == "CentralDifferences"):
            mechanical_scheme = KratosSolid.ExplicitCentralDifferencesScheme(self.explicit_solver_settings["max_delta_time"].GetDouble(),
                                                                             self.explicit_solver_settings["fraction_delta_time"].GetDouble(),
                                                                             self.explicit_solver_settings["time_step_prediction_level"].GetDouble(),
                                                                             self.explicit_solver_settings["rayleigh_damping"].GetBool())
        else:
            raise Exception("Unsupported integration_method: " + integration_method)

        return mechanical_scheme


    def _create_mechanical_solver(self):

        mechanical_solver = self._create_explicit_strategy()

        mechanical_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step

        return mechanical_solver


    def _create_explicit_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        return KratosSolid.ExplicitStrategy(self.model_part,
                                            mechanical_scheme,
                                            linear_solver,
                                            self.solving_strategy_settings["compute_reactions"].GetBool(),
                                            self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                            self.solving_strategy_settings["move_mesh_flag"].GetBool())
