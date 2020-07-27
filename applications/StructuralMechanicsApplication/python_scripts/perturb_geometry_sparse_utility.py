
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics import eigen_solver_factory

import numpy as np

class PerturbGeometrySparseUtility(KratosMultiphysics.Process):
    """A process to perturb the initial geometry of a structure
    based on a sparse correlation matrix.
    """
    def __init__(self, mp, settings ):
        """Constructor of Process-Object

        Checks parameter settings and initializes the process.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "eigensolver_settings"  : {
                "solver_type"               : "eigen_eigensystem",
                "max_iteration"             : 1000,
                "tolerance"                 : 1e-6,
                "number_of_eigenvalues"     : 1,
		        "normalize_eigenvectors"    : false,
                "echo_level"                : 0
                },
            "perturbation_settings" : {
                "max_displacement"          : 1,
                "correlation_length"        : 100,
                "truncation_error"          : 1e-3,
                "echo_level"                : 0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        eigen_solver = eigen_solver_factory.ConstructSolver(settings["eigensolver_settings"])
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        perturbation_settings = settings["perturbation_settings"]
        # Initialize process
        self.process = StructuralMechanicsApplication.PerturbGeometrySparseUtility(mp, eigen_solver, perturbation_settings)
        # Generate perturbation matrix
        self.number_random_variables = self.process.CreateRandomFieldVectors()

    def PerturbGeometry(self, mp ):
        """ Apply perturbation matrix to geometry.
        Random field approach requires normal distributed random numbers (mean=0, sigma=1)
        """
        random_numbers = np.random.normal(0, 1, self.number_random_variables)
        self.process.ApplyRandomFieldVectorsToGeometry(mp, random_numbers)





