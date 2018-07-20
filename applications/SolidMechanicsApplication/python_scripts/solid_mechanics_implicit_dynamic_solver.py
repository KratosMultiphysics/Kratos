from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_monolithic_solver as BaseSolver

def CreateSolver(custom_settings):
    return ImplicitMonolithicSolver(custom_settings)

class ImplicitMonolithicSolver(BaseSolver.MonolithicSolver):
    """The solid mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See solid_mechanics_monolithic_solver.py for more information.
    """
    def __init__(self, custom_settings):

        # Set defaults and validate custom settings.
        ##TODO : solving_strategy_settings must be time_integration_settings (GiD interface changes needed)
        implicit_solver_settings = KratosMultiphysics.Parameters("""
        {
            "solving_strategy_settings":{
                "bossak_factor" :-0.3,
                "dynamic_factor": 1.0,
                "lumped_mass_matrix" : true,
                "consistent_mass_matrix" : false,
                "rayleigh_damping": false,
                "rayleigh_alpha": 0.0,
                "rayleigh_beta" : 0.0
            }
        }
        """)

        # Validate and transfer settings
        if( custom_settings.Has("solving_strategy_settings") ):
            from json_settings_utility import JsonSettingsUtility
            JsonSettingsUtility.TransferMatchingSettingsToDestination(custom_settings["solving_strategy_settings"], implicit_solver_settings["solving_strategy_settings"])

        self.implicit_solver_settings = implicit_solver_settings["solving_strategy_settings"]

        # Construct the base solver.
        super(ImplicitMonolithicSolver, self).__init__(custom_settings)


    def GetVariables(self):

        nodal_variables = super(ImplicitMonolithicSolver, self).GetVariables()

        return nodal_variables

    #### Solver internal methods ####

    def _set_scheme_process_info_parameters(self):

        integration_method = self.settings["time_integration_settings"]["integration_method"].GetString()

        if( self.implicit_solver_settings["rayleigh_damping"].GetBool() == True ):
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = self.implicit_solver_settings["rayleigh_alpha"].GetDouble()
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = self.implicit_solver_settings["rayleigh_beta"].GetDouble()
        else:
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = 0.0
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = 0.0

        # compute dynamic tangent lhs and rhs
        self.process_info[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT] = False
        if( integration_method.find("Step") != -1 ):
            self.process_info[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT] = True

        # compute mass lumped matrix
        if( self.implicit_solver_settings["lumped_mass_matrix"].GetBool() == True ):
            self.process_info[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        else:
            self.process_info[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = False
            # compute consistent mass matrix
            if( self.implicit_solver_settings["consistent_mass_matrix"].GetBool() == True ):
                self.process_info[KratosSolid.COMPUTE_CONSISTENT_MASS_MATRIX] = True
            else:
                self.process_info[KratosSolid.COMPUTE_CONSISTENT_MASS_MATRIX] = False

        # set bossak factor
        if(integration_method.find("Bossak") != -1 or integration_method.find("Simo") != -1):
            bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
            self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;


    def _create_mechanical_solver(self):
        if(self.settings["solving_strategy_settings"]["line_search"].GetBool() == True):
            mechanical_solver = self._create_line_search_strategy()
        else:
            mechanical_solver = self._create_newton_raphson_strategy()
        return mechanical_solver


    def _create_line_search_strategy(self):
        solution_scheme = self._get_solution_scheme()
        #linear_solver = self._get_linear_solver()
        convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.settings["solving_strategy_settings"]["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.settings["solving_strategy_settings"]["reform_dofs_at_each_step"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.UPDATE_VARIABLES, self.settings["solving_strategy_settings"]["iterative_update"].GetBool())
        #options.Set(KratosSolid.SolverLocalFlags.MOVE_MESH, self.settings["solving_strategy_settings"]["move_mesh_flag"].GetBool())

        return KratosSolid.LineSearchStrategy(self.model_part, solution_scheme, builder_and_solver, convergence_criterion,
                                              options, self.settings["solving_strategy_settings"]["max_iteration"].GetInt())

    def _create_newton_raphson_strategy(self):
        solution_scheme = self._get_solution_scheme()
        #linear_solver = self._get_linear_solver()
        convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.settings["solving_strategy_settings"]["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.settings["solving_strategy_settings"]["reform_dofs_at_each_step"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.IMPLEX, self.settings["solving_strategy_settings"]["implex"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.UPDATE_VARIABLES, self.settings["solving_strategy_settings"]["iterative_update"].GetBool())
        #options.Set(KratosSolid.SolverLocalFlags.MOVE_MESH, self.settings["solving_strategy_settings"]["move_mesh_flag"].GetBool())

        return KratosSolid.NewtonRaphsonStrategy(self.model_part, solution_scheme, builder_and_solver, convergence_criterion,
                                                 options, self.settings["solving_strategy_settings"]["max_iteration"].GetInt())

    @classmethod
    def _class_prefix(self):
        header = "::[--Implicit_Solver--]::"
        return header
