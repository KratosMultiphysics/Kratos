from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return ShallowWaterSolver(model, custom_settings)

class ShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(ShallowWaterSolver, self).__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "SWE"
        self.condition_name = "Condition"
        self.min_buffer_size = 2
        self.advection_epsilon = self.settings["advection_epsilon"].GetDouble()

    def AddVariables(self):
        super(ShallowWaterSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.FREE_SURFACE_ELEVATION, self.main_model_part)

        KM.Logger.PrintInfo("::[ShallowWaterSolver]::", "Shallow water solver DOFs added correctly.")

    def FinalizeSolutionStep(self):
        super(ShallowWaterSolver, self).FinalizeSolutionStep()
        epsilon = max(self.advection_epsilon, self.main_model_part.ProcessInfo[SW.DRY_HEIGHT])
        SW.ShallowWaterUtilities().UpdatePrimitiveVariables(self.main_model_part, epsilon)
        SW.ShallowWaterUtilities().ComputeAccelerations(self.main_model_part)

    def InitializeSolutionStep(self):
        KM.VariableUtils().CopyVectorVar(KM.MOMENTUM, SW.PROJECTED_VECTOR1, self.main_model_part.Nodes)
        super(ShallowWaterSolver, self).InitializeSolutionStep()

    @classmethod
    def GetDefaultSettings(cls):
        default_settings = KM.Parameters("""
        {
            "advection_epsilon" : 1e-2
        }""")
        default_settings.AddMissingParameters(super(ShallowWaterSolver,cls).GetDefaultSettings())
        return default_settings