from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.TrilinosApplication as TrilinosApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.trilinos_structural_mechanics_solver import TrilinosMechanicalSolver

def CreateSolver(model, custom_settings):
    return TrilinosStaticMechanicalSolver(model, custom_settings)

class TrilinosStaticMechanicalSolver(TrilinosMechanicalSolver):
    """The trilinos structural mechanics static solver.

    For more information see:
    structural_mechanics_solver.py
    trilinos_structural_mechanics_solver.py
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(TrilinosStaticMechanicalSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[TrilinosStaticMechanicalSolver]:: ", "Construction finished")

    def _create_solution_scheme(self):
        return TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()
