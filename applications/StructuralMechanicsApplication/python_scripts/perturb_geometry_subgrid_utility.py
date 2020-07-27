
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics import eigen_solver_factory

import numpy as np

class PerturbGeometrySubgridUtility(KratosMultiphysics.Process):
    """A process to perturb the initial geometry of a structure
    based on a reduced correlation matrix.
    """
    def __init__(self, mp, settings ):
        """Constructor of Process-Object

        Checks parameter settings and initializes the process.
        """
        KratosMultiphysics.Process.__init__(self)

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
        # Initialize process
        self.process = StructuralMechanicsApplication.PerturbGeometrySubgridUtility(mp, eigen_solver, perturbation_settings)
        # Generate perturbation matrix
        self.number_random_variables = self.process.CreateRandomFieldVectors()


    def PerturbGeometry(self, mp ):
        """Apply perturbation matrix to geometry.
        Random field approach requires normal distributed random numbers (mean=0, sigma=1)
        """
        random_numbers = np.random.normal(0, 1, self.number_random_variables)
        self.process.ApplyRandomFieldVectorsToGeometry(mp, random_numbers)
