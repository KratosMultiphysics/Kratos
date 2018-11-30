from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(model, custom_settings):
    return StaticMechanicalSolver(model, custom_settings)


class StaticMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis. It currently
    supports line search, linear, arc-length, form-finding and Newton-Raphson
    strategies.

    Public member variables:
    arc_length_settings -- settings for the arc length method.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        static_settings = KratosMultiphysics.Parameters("""
        {
            "arc_length_settings" : {
                "Ide": 5,
                "factor_delta_lmax": 1.00,
                "max_iteration": 20,
                "max_recursive": 50,
                "toler": 1.0E-10,
                "norm": 1.0E-7,
                "MaxLineSearchIterations": 20,
                "tolls": 0.000001,
                "amp": 1.618,
                "etmxa": 5,
                "etmna": 0.1
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, static_settings)
        self.arc_length_settings = static_settings["arc_length_settings"]
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(StaticMechanicalSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[StaticMechanicalSolver]:: ", "Construction finished")

    def Initialize(self):
        self.print_on_rank_zero("::[StaticMechanicalSolver]:: ", "Initializing ...")
        if self.settings["analysis_type"].GetString() == "arc_length":
            self.main_model_part.ProcessInfo[StructuralMechanicsApplication.LAMBDA] = 0.0
        super(StaticMechanicalSolver, self).Initialize() # The mechanical solver is created here.
        self.print_on_rank_zero("::[StaticMechanicalSolver]:: ", "Finished initialization.")

    def SolveSolutionStep(self):
        super(StaticMechanicalSolver, self).SolveSolutionStep()
        if self.settings["analysis_type"].GetString() == "arc_length":
            raise Exception('"arc_length is not available at the moment"')
            # it is not clear if this should be called after SolveSolutionStep or FinalizeSolutionStep
            # this has to be updated once the ArcLength-Method is working again
            lambda_value = self.main_model_part.ProcessInfo[StructuralMechanicsApplication.LAMBDA]
            if self.settings["echo_level"].GetInt() > 0:
                self.print_on_rank_zero("LAMBDA: ", lambda_value)
            self._update_arc_length_point_load(lambda_value)

    #### Private functions ####

    # If needed, these functions can be implemented with fct from "kratos/utilities/variable_utils"
    #def _IncreasePointLoad(forcing_nodes_list, Load):
    #    for node in forcing_nodes_list:
    #        node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD, 0, Load)

    #def _IncreaseDisplacement(forcing_nodes_list, disp):
    #    for node in forcing_nodes_list:
    #        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, disp)

    def _update_arc_length_point_load(self, lambda_value):
        force_x = 0.0
        force_y = 0.0
        force_z = 0.0
        for node in self.main_model_part.Nodes:
            new_load = node.GetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD, 0) * lambda_value;
            force_x += new_load[0]
            force_y += new_load[1]
            force_z += new_load[2]
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD, 0, new_load)

        if (self.settings["echo_level"].GetInt() > 0):
            self.print_on_rank_zero("Total Load Applied", "\nPOINT_LOAD_X: " + str(force_x) + "\nPOINT_LOAD_Y: " + str(force_y) + "\nPOINT_LOAD_Z: " + str(force_z))

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            mechanical_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                mechanical_solution_strategy = self._create_newton_raphson_strategy()
            else:
                mechanical_solution_strategy = self._create_line_search_strategy()
        elif analysis_type == "arc_length":
            raise Exception('"arc_length is not available at the moment"')
            mechanical_solution_strategy = self._create_arc_length_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\", \"arc_length\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_arc_length_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return StructuralMechanicsApplication.ResidualBasedArcLengthStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                linear_solver,
                                                                mechanical_convergence_criterion,
                                                                self.arc_length_settings["Ide"].GetInt(),
                                                                self.arc_length_settings["max_iteration"].GetInt(),
                                                                self.arc_length_settings["max_recursive"].GetInt(),
                                                                self.arc_length_settings["factor_delta_lmax"].GetDouble(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool())
