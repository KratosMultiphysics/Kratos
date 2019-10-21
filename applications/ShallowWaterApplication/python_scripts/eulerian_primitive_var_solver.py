from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return EulerianPrimitiveVarSolver(model, custom_settings)

class EulerianPrimitiveVarSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(EulerianPrimitiveVarSolver, self).__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "ReducedSWE"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.VELOCITY_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.VELOCITY_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo("::[EulerianPrimitiveVarSolver]::", "Shallow water solver DOFs added correctly.")

    def InitializeSolutionStep(self):
        KM.VariableUtils().CopyScalarVar(SW.HEIGHT, SW.PROJECTED_SCALAR1, self.main_model_part.Nodes)
        KM.VariableUtils().CopyVectorVar(KM.VELOCITY, SW.PROJECTED_VECTOR1, self.main_model_part.Nodes)
        super(EulerianPrimitiveVarSolver, self).InitializeSolutionStep()