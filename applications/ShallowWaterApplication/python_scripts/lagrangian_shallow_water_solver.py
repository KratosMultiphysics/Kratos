from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_solver import ShallowWaterSolver
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(LagrangianShallowWaterSolver, self).__init__(model, settings)
        self.min_buffer_size = 2
        self.element_name = "ConservedElement"
        self.condition_name = "LineCondition"
        self.thickness = 1e-2

    def AddVariables(self):
        super(LagrangianShallowWaterSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(SW.PROJECTED_SCALAR1)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

    def Initialize(self):
        super(LagrangianShallowWaterSolver, self).Initialize()
        self.bfecc = SW.BFECCConvectionUtility(self.GetComputingModelPart())

    def InitializeSolutionStep(self):
        self.bfecc.Convect(KM.MOMENTUM, KM.VELOCITY)
        self.bfecc.CopyVariableToPreviousTimeStep(KM.MOMENTUM)

        SW.ShallowWaterUtilities().IdentifyWetDomain(self.GetComputingModelPart(), KM.FLUID, self.thickness)
        SW.ShallowWaterUtilities().IdentifyWetDomain(self.GetComputingModelPart(), KM.ACTIVE, self.thickness)

        super(LagrangianShallowWaterSolver, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(LagrangianShallowWaterSolver, self).FinalizeSolutionStep()
        KM.VariableUtils().SetFlag(KM.ACTIVE, True, self.GetComputingModelPart().Elements)
        SW.ShallowWaterUtilities().ResetDryDomain(self.GetComputingModelPart(), self.thickness)
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.GetComputingModelPart())
        SW.ComputeVelocityProcess(self.GetComputingModelPart(), 1e-1).Execute()
        SW.ShallowWaterUtilities().ComputeAccelerations(self.GetComputingModelPart())
