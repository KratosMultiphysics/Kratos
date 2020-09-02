from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver

def CreateSolver(model, custom_settings):
    return EulerianShallowWaterSolver(model, custom_settings)

class EulerianShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super(EulerianShallowWaterSolver, self).__init__(model, settings)

        # Set the element and condition names for the replace settings
        self.element_name = "ShallowWater"
        self.condition_name = "LineCondition"
        self.min_buffer_size = 2

    def AddVariables(self):
        super(EulerianShallowWaterSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(SW.ATMOSPHERIC_PRESSURE)

    def AddDofs(self):
        KM.VariableUtils().AddDof(KM.MOMENTUM_X, self.main_model_part)
        KM.VariableUtils().AddDof(KM.MOMENTUM_Y, self.main_model_part)
        KM.VariableUtils().AddDof(SW.HEIGHT, self.main_model_part)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shallow water solver DOFs added correctly.")

    def Initialize(self):
        super(EulerianShallowWaterSolver, self).Initialize()
        self.main_model_part.ProcessInfo[SW.LUMPED_MASS_FACTOR] = 0.0
        self.main_model_part.ProcessInfo[KM.STABILIZATION_FACTOR] = 0.1
        self.main_model_part.ProcessInfo[SW.SHOCK_STABILIZATION_FACTOR] = 0.01
        self.main_model_part.ProcessInfo[SW.GROUND_IRREGULARITY] = 0.01

    def FinalizeSolutionStep(self):
        super(EulerianShallowWaterSolver, self).FinalizeSolutionStep()
        SW.ShallowWaterUtilities().ComputeFreeSurfaceElevation(self.main_model_part)
        SW.ShallowWaterUtilities().ResetDryDomain(self.main_model_part, 0.005)
        SW.ComputeVelocityProcess(self.main_model_part, 0.01).Execute()
        self._CheckWaterLoss()

    def _InitializeWaterLoss(self):
        self.initial_water = KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.HEIGHT, self.main_model_part,0)
        self.initial_water /= self.main_model_part.NumberOfNodes()

    def _CheckWaterLoss(self):
        if not hasattr(self, 'initial_water'):
            self._InitializeWaterLoss()
        total_water = KM.VariableUtils().SumHistoricalNodeScalarVariable(SW.HEIGHT, self.main_model_part,0)
        total_water /= self.main_model_part.NumberOfNodes()
        water_loss = (total_water - self.initial_water) / self.initial_water
        if abs(water_loss) > 1e-3:
            msg = "Water loss : {} %"
            KM.Logger.PrintWarning(self.__class__.__name__, msg.format(water_loss*100))
