from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorCountResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorCountResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorCountResponse(parameters["name"].GetString(), model, parameters["settings"])


class SensorCountResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name": ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorCountResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        nodal_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_exp, KratosSI.SENSOR_STATUS, False)
        return Kratos.Expression.Utils.Sum(nodal_exp)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        nodal_exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(nodal_exp, 1.0)
        physical_variable_collective_expressions[KratosSI.SENSOR_STATUS].GetContainerExpressions()[0].SetExpression(nodal_exp.GetExpression())

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"