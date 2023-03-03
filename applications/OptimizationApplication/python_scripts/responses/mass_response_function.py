import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction

class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def Check(self) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_part)

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.MassResponseUtils.CalculateValue(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.CalculateSensitivity(sensitivity_model_part, sensitivity_variable)

