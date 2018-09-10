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
    return EigenSolver(main_model_part, custom_settings)


class EigenSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.

    Public member variables:
    eigensolver_settings -- settings for the eigenvalue solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Validation of eigensolver_settings is done in the eigenvalue solvers
        self.eigensolver_settings = custom_settings["eigensolver_settings"]

        # Validate the remaining settings in the base class.
        structural_settings = custom_settings.Clone()
        structural_settings.RemoveValue("eigensolver_settings")

        self.structural_eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_type"   : "dynamic"
        }
        """)
        self.validate_and_transfer_matching_settings(structural_settings, self.structural_eigensolver_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(EigenSolver, self).__init__(main_model_part, structural_settings)
        self.print_on_rank_zero("::[EigenSolver]:: ", "Construction finished")

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme for the eigenvalue problem.

        The scheme determines the left- and right-hand side matrices in the
        generalized eigenvalue problem.
        """
        scheme_type = self.structural_eigensolver_settings["scheme_type"].GetString()
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
        import eigen_solver_factory
        return eigen_solver_factory.ConstructSolver(self.eigensolver_settings)

    def _create_mechanical_solution_strategy(self):
        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self.get_builder_and_solver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.EigensolverStrategy(computing_model_part,
                                                                  eigen_scheme,
                                                                  builder_and_solver)
