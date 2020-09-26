# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_solver import ShallowWaterSolver
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(ShallowWaterSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.min_buffer_size = 2
        self.element_name = "LagrangianSWE"
        self.condition_name = "LineCondition"

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(SW.PROJECTED_SCALAR1)

    def Initialize(self):
        super().Initialize()
        self.bfecc = SW.BFECCConvectionUtility(self.GetComputingModelPart())

    def InitializeSolutionStep(self):
        self.bfecc.Convect(KM.MOMENTUM, KM.VELOCITY)
        self.bfecc.CopyVariableToPreviousTimeStep(KM.MOMENTUM)
        super().InitializeSolutionStep()
