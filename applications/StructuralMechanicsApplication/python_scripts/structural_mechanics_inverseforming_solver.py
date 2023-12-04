# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return InverseFormingMechanicalSolver(model, custom_settings)

class InverseFormingMechanicalSolver(MechanicalSolver):
    """The structural mechanics inverse forming solver.

    This class creates the mechanical solver for inverseforming.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[InverseFormingMechanicalSolver]:: ", "Construction finished")
    
    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()