import KratosMultiphysics
import KratosMultiphysics.TopologyOptimizationApplication as TopologyOptimizationApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):       
    return SIMPStaticMechanicalSolver(model, custom_settings)

    
class SIMPStaticMechanicalSolver(MechanicalSolver):

    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[SIMPStaticMechanicalSolver]:: ", "Construction finished")
        
    def _CreateScheme(self):
        return TopologyOptimizationApplication.ResidualBasedIncrementalUpdateStaticSIMPScheme()