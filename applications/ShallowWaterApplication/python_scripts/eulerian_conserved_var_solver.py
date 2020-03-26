from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return EulerianConservedVarSolver(model, custom_settings)

class EulerianConservedVarSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(EulerianConservedVarSolver, self).__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "ConservedElement"
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2

    def AddVariables(self):
        super(EulerianConservedVarSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo("::[EulerianConservedVarSolver]::", "Shallow water solver DOFs added correctly.")

    def FinalizeSolutionStep(self):
        super(EulerianConservedVarSolver, self).FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeVelocity(self.main_model_part)
