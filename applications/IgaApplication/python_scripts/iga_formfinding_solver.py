from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as IgaApplication

# Import base class file
from KratosMultiphysics.IgaApplication.iga_solver import IgaSolver


def CreateSolver(model, custom_settings):
    return IgaFormfindingSolver(model, custom_settings)


class IgaFormfindingSolver(IgaSolver):
    """The iga formfinding solver.

    This class creates the iga solver for formdinding.

    Public member variables:
    formfinding_settings -- settings for the formfinding solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        self.formfinding_settings = KratosMultiphysics.Parameters("""
        {
            "print_formfinding_iterations": false
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.formfinding_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(IgaFormfindingSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[IgaFormfindingSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_mechanical_solution_strategy(self):
        self.print_on_rank_zero("::[IgaFormfindingSolver]:: ", "_create_mechanical_solution_strategy")
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
