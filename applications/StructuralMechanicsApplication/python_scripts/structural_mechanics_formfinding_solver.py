from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(main_model_part, custom_settings):
    return FormfindingMechanicalSolver(main_model_part, custom_settings)

class FormfindingMechanicalSolver(MechanicalSolver):
    """The structural mechanics formfinding solver.

    This class creates the mechanical solver for formdinding.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super(FormfindingMechanicalSolver, self).__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[FormfindingMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "print_formfinding_iterations" : false
        }""")
        this_defaults.AddMissingParameters(super(FormfindingMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults

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
                                                                self.settings["print_formfinding_iterations"].GetBool(),
                                                                self.settings["line_search"].GetBool())
