# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

from KratosMultiphysics import eigen_solver_factory

def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)

class EigenSolver(MechanicalSolver):
    """The structural mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[EigenSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"         : "dynamic",
            "compute_modal_decomposition": false,
            "eigensolver_settings" : {
                "solver_type"           : "eigen_eigensystem",
                "max_iteration"         : 1000,
                "tolerance"             : 1e-6,
                "number_of_eigenvalues" : 5,
                "echo_level"            : 1
            },
            "eigensolver_diagonal_values" : { }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultSettings())
        return this_defaults

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme for the eigenvalue problem.

        The scheme determines the left- and right-hand side matrices in the
        generalized eigenvalue problem.
        """
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else: # here e.g. a stability scheme could be added
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)

        return solution_scheme

    def _create_linear_solver(self):
        """Create the eigensolver.

        This overrides the base class method and replaces the usual linear solver
        with an eigenvalue problem solver.
        """
        return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])

    def _create_mechanical_solution_strategy(self):
        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self.get_builder_and_solver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        solver_type = self.settings["eigensolver_settings"]["solver_type"].GetString()
        if solver_type == "eigen_eigensystem":
            mass_matrix_diagonal_value = 0.0
            stiffness_matrix_diagonal_value = 1.0
        elif solver_type == "feast":
            mass_matrix_diagonal_value = 1.0
            stiffness_matrix_diagonal_value = -1.0
        else:
            diag_values = self.settings["eigensolver_diagonal_values"]
            if not diag_values.Has("mass_matrix_diagonal_value") or not diag_values.Has("stiffness_matrix_diagonal_value"):
                err_msg  = 'For the used eigensolver "{}" no defaults for '.format(solver_type)
                err_msg += '"mass_matrix_diagonal_value" and "stiffness_matrix_diagonal_value" exist, '
                err_msg += 'please specify them under "eigensolver_diagonal_values"'
                raise Exception(err_msg)

            mass_matrix_diagonal_value = diag_values["mass_matrix_diagonal_value"].GetDouble()
            stiffness_matrix_diagonal_value = diag_values["stiffness_matrix_diagonal_value"].GetDouble()

        return StructuralMechanicsApplication.EigensolverStrategy(computing_model_part,
                                                                  eigen_scheme,
                                                                  builder_and_solver,
                                                                  mass_matrix_diagonal_value,
                                                                  stiffness_matrix_diagonal_value,
                                                                  self.settings["compute_modal_decomposition"].GetBool())
