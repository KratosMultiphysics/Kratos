from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_implicit_dynamic_solver as BaseSolver

def CreateSolver(custom_settings):
    return StaticMechanicalSolver(custom_settings)

class StaticMechanicalSolver(BaseSolver.ImplicitMechanicalSolver):
    """The solid mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    Public member variables:

    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, custom_settings):

        # Set defaults and validate custom settings.
        static_settings = KratosMultiphysics.Parameters("""
        {

        }
        """)

         # Validate and transfer settings
        self._validate_and_transfer_matching_settings(custom_settings, static_settings)
        time_integration_settings = custom_settings["time_integration_settings"]

        # Validate the remaining settings in the base class.
        if not time_integration_settings.Has("solution_type"):
            time_integration_settings.AddEmptyValue("solution_type")
            time_integration_settings["solution_type"].SetString("Static") # Override defaults in the base class.
        if not time_integration_settings.Has("integration_method"):
            time_integration_settings.AddEmptyValue("integration_method")
            time_integration_settings["integration_method"].SetString("Non-Linear") # Override defaults in the base class.

        # Construct the base solver.
        super(BaseSolver.ImplicitMechanicalSolver, self).__init__(custom_settings)

        print("::[Static_Scheme]:: "+self.time_integration_settings["integration_method"].GetString()+" Scheme Ready")

    #### Solver internal methods ####

    def _create_solution_scheme (self):

        self.linear_integration_method = None
        self.angular_integration_method = None
        
        integration_method   = self.time_integration_settings["integration_method"].GetString()

        if(integration_method == "Linear"):
            #mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            self.linear_integration_method = KratosSolid.StaticMethod()
            self.angular_integration_method = KratosSolid.StaticMethod() #shells
            mechanical_scheme = KratosSolid.DisplacementStaticScheme()
        elif(integration_method == "Non-Linear" ):
            #mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            self.linear_integration_method = KratosSolid.StaticMethod()
            self.angular_integration_method = KratosSolid.StaticMethod() #shells
            mechanical_scheme = KratosSolid.DisplacementStaticScheme()
        elif(integration_method == "RotationStatic"):
            #dynamic_factor = 0.0
            #damp_factor_m  = 0.0
            #mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, damp_factor_m)
            self.linear_integration_method = KratosSolid.StaticStepMethod()
            self.angular_integration_method = KratosSolid.StaticStepRotationMethod()
            mechanical_scheme = KratosSolid.DisplacementRotationStaticScheme()
        else:
            raise Exception("Unsupported integration_method: " + integration_method)

        # set integration parameters
        self._set_time_integration_methods()

        return mechanical_scheme

    def _create_mechanical_solver(self):
        if(self.solving_strategy_settings["line_search"].GetBool() == True):
            mechanical_solver = self._create_line_search_strategy()
        else:
            if(self.time_integration_settings["integration_method"].GetString() == "Linear"):
                mechanical_solver = self._create_linear_strategy()
            else:
                mechanical_solver = self._create_newton_raphson_strategy()

        return mechanical_solver


    def _create_linear_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        builder_and_solver = self._get_builder_and_solver()
        
        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.solving_strategy_settings["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool())

        return KratosSolid.LinearStrategy(self.model_part, mechanical_scheme, builder_and_solver, options)

        
        #return KratosMultiphysics.ResidualBasedLinearStrategy(self.model_part,
        #                                                      mechanical_scheme,
        #                                                      linear_solver,
        #                                                      builder_and_solver,
        #                                                      self.solving_strategy_settings["compute_reactions"].GetBool(),
        #                                                      self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
        #                                                      False,
        #                                                      self.solving_strategy_settings["move_mesh_flag"].GetBool())


