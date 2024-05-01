import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_adjoint_static_solver import StructuralMechanicsAdjointStaticSolver
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionDataLocation

class SensorSensitivityAdjointStaticSolver(StructuralMechanicsAdjointStaticSolver):
    def __init__(self, model: Kratos.Model, custom_settings: Kratos.Parameters):
        super().__init__(model, custom_settings)
        self.sensitivity_builder_initialized = False

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = Kratos.Parameters("""{
            "sensitivity_settings": {}
        }""")
        this_defaults.AddMissingParameters(super(StructuralMechanicsAdjointStaticSolver, cls).GetDefaultParameters())
        return this_defaults

    def Initialize(self) -> None:
        # create sensitivity builder
        self.sensitivity_builder = Kratos.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, None)

        # skip calling Initialize of StructuralMechanicsAdjointStaticSolver and call the its base class Initialize
        super(StructuralMechanicsAdjointStaticSolver, self).Initialize()

        # set the strategy to only build LHS once, and build RHS for all other solves
        self._GetSolutionStrategy().SetRebuildLevel(0)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished initialization.")

    def InitializeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def SetSensor(self, sensor: KratosSI.Sensors.Sensor) -> None:
        self.sensitivity_builder.SetResponseFunction(sensor)
        self._GetScheme().SetResponseFunction(sensor)
        if not self.sensitivity_builder_initialized:
            self.sensitivity_builder.Initialize()
            self.sensitivity_builder_initialized = True

    def GetSensitivtyVariables(self) -> 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]':
        sensitivity_settings = self.settings["sensitivity_settings"]
        return {
            ExpressionDataLocation.NodeHistorical: [Kratos.KratosGlobals.GetVariable(f"{var_name}_SENSITIVITY") for var_name in sensitivity_settings["nodal_solution_step_sensitivity_variables"].GetStringArray()],
            ExpressionDataLocation.Condition: [Kratos.KratosGlobals.GetVariable(f"{var_name}_SENSITIVITY") for var_name in sensitivity_settings["condition_data_value_sensitivity_variables"].GetStringArray()],
            ExpressionDataLocation.Element: [Kratos.KratosGlobals.GetVariable(f"{var_name}_SENSITIVITY") for var_name in sensitivity_settings["element_data_value_sensitivity_variables"].GetStringArray()]
        }

    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        return self.main_model_part.GetSubModelPart(self.settings["sensitivity_settings"]["sensitivity_model_part_name"].GetString())

    def SolveSolutionStep(self) -> bool:
        return super(StructuralMechanicsAdjointStaticSolver, self).SolveSolutionStep()

    def _CreateScheme(self):
        return Kratos.ResidualBasedAdjointStaticScheme(None)