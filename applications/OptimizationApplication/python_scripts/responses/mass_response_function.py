import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes

class MassResponseFunction(ResponseFunction):
    @classmethod
    def GetSensitivityFieldVariables(cls) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.SHAPE_SENSITIVITY, KratosOA.DENSITY_SENSITIVITY, KratosOA.THICKNESS_SENSITIVITY, KratosOA.CROSS_AREA_SENSITIVITY]

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_parts: 'list[Kratos.ModelPart]' = []
        for model_part_name in parameters["evaluated_model_part_names"].GetStringArray():
            self.model_parts.append(model[model_part_name])

        if len(self.model_parts) == 0:
            raise RuntimeError("No model parts were provided for MassResponseFunction.")

    def Check(self) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_parts)

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.MassResponseUtils.CalculateValue(self.model_parts)

    def CalculateSensitivity(self, sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]') -> None:
        KratosOA.ResponseUtils.MassResponseUtils.CalculateSensitivity(self.model_parts, sensitivity_model_part_variable_info)

