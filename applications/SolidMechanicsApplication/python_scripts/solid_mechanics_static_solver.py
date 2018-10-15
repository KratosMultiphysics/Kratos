from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_implicit_dynamic_solver as BaseSolver

def CreateSolver(custom_settings, Model):
    return StaticMonolithicSolver(Model, custom_settings)

class StaticMonolithicSolver(BaseSolver.ImplicitMonolithicSolver):
    """The solid mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    Public member variables:

    See solid_mechanics_monolithic_solver.py for more information.
    """
    def __init__(self, Model, custom_settings):

        # Set defaults and validate custom settings.
        static_settings = KratosMultiphysics.Parameters("""
        {

        }
        """)

        # Validate and transfer settings
        from json_settings_utility import JsonSettingsUtility
        JsonSettingsUtility.TransferMatchingSettingsToDestination(custom_settings, static_settings)
        time_integration_settings = custom_settings["time_integration_settings"]

        # Validate the remaining settings in the base class.
        if not time_integration_settings.Has("integration_method"):
            time_integration_settings.AddEmptyValue("integration_method")
            time_integration_settings["integration_method"].SetString("Static") # Override defaults in the base class.

        # Construct the base solver.
        # Calling base class of ImplicitMonolithicSolver it is ok.
        super(BaseSolver.ImplicitMonolithicSolver, self).__init__(Model, custom_settings)


    #### Solver internal methods ####

    def _set_scheme_process_info_parameters(self):
        pass

    def _create_mechanical_solver(self):
        if(self.settings["solving_strategy_settings"]["line_search"].GetBool() == True):
            mechanical_solver = self._create_line_search_strategy()
        else:
            if(self.settings["time_integration_settings"]["analysis_type"].GetString() == "Non-linear"):
                mechanical_solver = self._create_newton_raphson_strategy()
            else:
                mechanical_solver = self._create_linear_strategy()
        mechanical_solver.Set(KratosSolid.SolverLocalFlags.ADAPTIVE_SOLUTION,self.settings["solving_strategy_settings"]["adaptive_solution"].GetBool())
        return mechanical_solver


    def _create_linear_strategy(self):
        solution_scheme = self._get_solution_scheme()
        builder_and_solver = self._get_builder_and_solver()

        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.settings["solving_strategy_settings"]["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.settings["solving_strategy_settings"]["reform_dofs_at_each_step"].GetBool())

        return KratosSolid.LinearStrategy(self.model_part, solution_scheme, builder_and_solver, options)

    @classmethod
    def _class_prefix(self):
        header = "::[---Static_Solver---]::"
        return header
