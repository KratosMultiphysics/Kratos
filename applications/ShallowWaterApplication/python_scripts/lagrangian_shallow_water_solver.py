from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_solver import ShallowWaterSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(ShallowWaterSolver):
    def __init__(self, model, settings):
        super(LagrangianShallowWaterSolver, self).__init__(model, settings)
        self.element_name = "LagrangianSWE"

    def Initialize(self):
        super(LagrangianShallowWaterSolver, self).Initialize()
        self.bfecc = SW.BFECCConvectionUtility(self.GetComputingModelPart())

    def InitializeSolutionStep(self):
        self.bfecc.Convect(KM.MOMENTUM, KM.VELOCITY)
        self.bfecc.CopyVariableToPreviousTimeStep(KM.MOMENTUM)
        super(LagrangianShallowWaterSolver, self).InitializeSolutionStep()