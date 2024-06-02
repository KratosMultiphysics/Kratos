import math, numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"CombinedResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"CombinedResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return CombinedResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class CombinedResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "list_of_responses": [
                {
                    "response_function_name": "",
                    "weight"                : 1.0
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.optimization_problem = optimization_problem

        self.list_of_response_settings: 'list[tuple[ResponseFunction, float]]' = []
        for response_params in parameters["list_of_responses"].values():
            response_params.ValidateAndAssignDefaults(default_settings["list_of_responses"].values()[0])

            response_function_name = response_params["response_function_name"].GetString()
            weight = response_params["weight"].GetDouble()

            self.list_of_response_settings.append((self.optimization_problem.GetResponse(response_function_name), weight))

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosDT.SENSOR_STATUS]

    def Initialize(self) -> None:
        component_data_view = ComponentDataView(self, self.optimization_problem)
        component_data_view.SetDataBuffer(1)
        self.buffer_data = component_data_view.GetBufferedData()

        for response, _ in self.list_of_response_settings:
            response.Initialize()

    def Check(self) -> None:
        for response, _ in self.list_of_response_settings:
            response.Check()


    def Finalize(self) -> None:
        for response, _ in self.list_of_response_settings:
            response.Finalize()

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        return self.list_of_response_settings[0][0].GetEvaluatedModelPart()

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return None

    def CalculateValue(self) -> float:
        value = 0.0
        for response, weight in self.list_of_response_settings:
            response_value = response.CalculateValue()
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"{response.GetName()} value = {response_value}")
            # self.buffer_data.SetValue(f"{response.GetName()}_value", response_value)
            value += weight * response_value
        return value

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        temp_collective_dict: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
        for k, v in physical_variable_collective_expressions.items():
            temp_collective_dict[k] = v.Clone()

        for response, weight in self.list_of_response_settings:
            response.CalculateGradient(temp_collective_dict)
            for k, v in temp_collective_dict.items():
                physical_variable_collective_expressions[k].GetContainerExpressions()[0].SetExpression(physical_variable_collective_expressions[k].GetContainerExpressions()[0].GetExpression() + v.GetContainerExpressions()[0].GetExpression() * weight)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = ]"