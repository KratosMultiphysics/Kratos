from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)


class EigenSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.
    It currently supports the Feast solver.

    Public member variables:
    eigensolver_settings -- settings for the eigenvalue solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            "eigensolver_settings" : {
                "solver_type": "FEAST",
                "print_feast_output": true,
                "perform_stochastic_estimate": true,
                "solve_eigenvalue_problem": true,
                "lambda_min": 0.0,
                "lambda_max": 1.0,
                "search_dimension": 10,
                "linear_solver_settings": {
                    "solver_type": "skyline_lu"
                }
            }
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, eigensolver_settings)
        self.eigensolver_settings = eigensolver_settings["eigensolver_settings"]
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("dynamic")
        
        # Construct the base solver.
        super(EigenSolver, self).__init__(main_model_part, custom_settings)
        print("::[EigenSolver]:: Construction finished")

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
        import eigen_solver_factory
        return eigen_solver_factory.ConstructSolver(self.eigensolver_settings)

    def _create_mechanical_solver(self):
        eigen_scheme = self.get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self.get_builder_and_solver() # The eigensolver is created here.
        computing_model_part = self.GetComputingModelPart()

        return StructuralMechanicsApplication.EigensolverStrategy(computing_model_part,
                                                                  eigen_scheme,
                                                                  builder_and_solver)
