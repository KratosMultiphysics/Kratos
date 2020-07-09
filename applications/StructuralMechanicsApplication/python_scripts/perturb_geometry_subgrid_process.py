
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics import eigen_solver_factory

import numpy as np
# def Factory(settings, Model):
#     if(type(settings) != KratosMultiphysics.Parameters):
#         raise Exception("expected input shall be a Parameters object, encapsulating a json string")
#     return PerturbGeometrySparseProcessPython(Model, settings["Parameters"])

class PerturbGeometrySubgridProcess(KratosMultiphysics.Process):
    """A process to perturb the initial geometry of a structure based on a reduced correlation matrix."""

    def __init__(self, mp, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "eigensolver_settings"  : {
                "solver_type"               : "dense_eigensolver",
                "ascending_order"           : true
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

        #TODO: Add warning if sorted in wrong order
        eigen_solver = eigen_solver_factory.ConstructSolver(settings["eigensolver_settings"])
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        perturbation_settings = settings["perturbation_settings"]
        # Initialize process
        self.process = StructuralMechanicsApplication.PerturbGeometrySubgridProcess(mp, eigen_solver, perturbation_settings)
        # Generate perturbation matrix
        self.number_random_variables = self.process.CreateRandomFieldVectors()


    def PerturbGeometry(self, mp ):
        # Apply perturbation matrix to geometry
        # Random field approach requires normal distributed random numbers (mean=0, sigma=1)
        random_numbers = np.random.normal(0, 1, self.number_random_variables)
        self.process.ApplyRandomFieldVectorsToGeometry(mp, random_numbers)
