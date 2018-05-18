from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(main_model_part, custom_settings):
    return FormfindingMechanicalSolver(main_model_part, custom_settings)


class FormfindingMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics formfinding solver.

    This class creates the mechanical solver for formdinding.

    Public member variables:
    formfinding_settings -- settings for the formfinding solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.formfinding_settings = KratosMultiphysics.Parameters("""
        {
            "print_formfinding_iterations": false
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.formfinding_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(FormfindingMechanicalSolver, self).__init__(main_model_part, custom_settings)
        self.print_on_rank_zero("::[FormfindingMechanicalSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return StructuralMechanicsApplication.FormfindingUpdatedReferenceStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                linear_solver,
                                                                mechanical_convergence_criterion,
                                                                builder_and_solver,
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool(),
                                                                self.formfinding_settings["print_formfinding_iterations"].GetBool(),
                                                                self.settings["line_search"].GetBool())
