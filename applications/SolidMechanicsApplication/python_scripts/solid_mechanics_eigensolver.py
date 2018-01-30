from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import solid_mechanics_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)


class EigenSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics eigen solver.

    This class creates the mechanical solvers for eigenvalue analysis.
    It currently supports the Feast solver.

    Public member variables:
    eigensolver_settings -- settings for the eigenvalue solvers.

    See solid_mechanics_solver.py for more information.
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
                "compute_modal_contribution": false,
                "lambda_min": 0.0,
                "lambda_max": 1.0,
                "search_dimension": 10,
                "linear_solver_settings": {
                    "solver_type": "complex_skyline_lu_solver"
                }
            }
        }
        """)
        self._validate_and_transfer_matching_settings(custom_settings, eigensolver_settings)
        self.eigensolver_settings = eigensolver_settings["eigensolver_settings"]
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Dynamic")

        self.compute_modal_contribution = self.eigensolver_settings["compute_modal_contribution"].GetBool()
        self.eigensolver_settings.RemoveValue("compute_modal_contribution")
        
        # Construct the base solver.
        super(EigenSolver, self).__init__(main_model_part, custom_settings)
        print("::[Eigen_Solver]:: Construction finished")

    #### Private functions ####

    def _create_solution_scheme(self):
        """Create the scheme for the eigenvalue problem.

        The scheme determines the left- and right-hand side matrices in the
        generalized eigenvalue problem. 
        """
        if self.settings["solution_type"].GetString() == "Dynamic":
            solution_scheme = KratosSolid.EigensolverDynamicScheme()
        else:
            raise Exception("Unsupported solution_type: " + self.settings["solution_type"])
        return solution_scheme

    def _create_linear_solver(self):
        """Create the eigensolver.
        
        This overrides the base class method and replaces the usual linear solver
        with an eigenvalue problem solver.
        """        
        if self.eigensolver_settings["solver_type"].GetString() == "FEAST":
            feast_system_solver_settings = self.eigensolver_settings["linear_solver_settings"]
            if feast_system_solver_settings["solver_type"].GetString() == "complex_skyline_lu_solver":
                # default built-in feast system solver
                linear_solver = ExternalSolversApplication.FEASTSolver(self.eigensolver_settings)
            elif feast_system_solver_settings["solver_type"].GetString() == "pastix":
                feast_system_solver = ExternalSolversApplication.PastixComplexSolver(feast_system_solver_settings)
                linear_solver = ExternalSolversApplication.FEASTSolver(self.eigensolver_settings, feast_system_solver)
            else:
                raise Exception("Unsupported FEAST system solver_type: " + feast_system_solver_settings["solver_type"].GetString())
        else:
            raise Exception("Unsupported eigensolver solver_type: " + self.eigensolver_settings["solver_type"].GetString())
        return linear_solver

    def _create_mechanical_solver(self):
        computing_model_part = self.GetComputingModelPart()
        eigen_scheme = self._get_solution_scheme() # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self._get_builder_and_solver() # The eigensolver is created here.

        return KratosSolid.EigensolverStrategy(computing_model_part,
                                               eigen_scheme,
                                               builder_and_solver,
                                               self.compute_modal_contribution)
