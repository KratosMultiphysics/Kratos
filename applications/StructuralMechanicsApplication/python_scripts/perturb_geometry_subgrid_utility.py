
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics import eigen_solver_factory

import random

class PerturbGeometrySubgridUtility():
    """An utility to perturb the initial geometry of a structure
    based on a reduced correlation matrix.
    """
    def __init__(self, mp, settings ):
        """Constructor of Utility-Object

        Checks parameter settings and initializes the utility.
        """
        default_settings = KratosMultiphysics.Parameters("""{
            "eigensolver_settings"  : {
                "solver_type"               : "dense_eigensolver",
                "ascending_order"           : false
                },
            "perturbation_settings" : {
                "min_distance_subgrid"      : 10,
                "max_displacement"          : 1,
                "correlation_length"        : 100,
                "truncation_error"          : 1e-3,
                "echo_level"                : 0
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        eigensolver_settings = settings["eigensolver_settings"]

        if( eigensolver_settings["ascending_order"].GetBool() ):
            warn_msg = 'Eigenvalues must be sorted in descending order. \n'
            warn_msg += '''Parameter: '"ascending_order" : true' specification is ignored.'''
            KratosMultiphysics.Logger.PrintWarning("PerturbGeometrySubgridUtility", warn_msg)
            eigensolver_settings["ascending_order"].SetBool(False)

        eigen_solver = eigen_solver_factory.ConstructSolver(eigensolver_settings)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        perturbation_settings = settings["perturbation_settings"]
        # Initialize utility
        self.utility = StructuralMechanicsApplication.PerturbGeometrySubgridUtility(mp, eigen_solver, perturbation_settings)
        # Generate perturbation matrix
        self.number_random_variables = self.utility.CreateRandomFieldVectors()


    def PerturbGeometry(self, mp ):
        """Apply perturbation matrix to geometry.
        Random field approach requires normal distributed random numbers (mean=0, sigma=1)
        """
        random_numbers = [random.gauss(0, 1) for i in range(self.number_random_variables)]
        self.utility.ApplyRandomFieldVectorsToGeometry(mp, random_numbers)
