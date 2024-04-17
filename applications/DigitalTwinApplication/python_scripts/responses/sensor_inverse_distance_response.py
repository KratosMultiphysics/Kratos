from typing import Optional
from math import exp
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorInverseDistanceResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorInverseDistanceResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorInverseDistanceResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class SensorInverseDistanceResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for SensorInverseDistanceResponse. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosDT.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        # get the distance matrix
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = ComponentDataView("sensor", self.optimization_problem).GetUnBufferedData().GetValue("list_of_sensors")

        n = len(list_of_sensors)
        self.distance_matrix = Kratos.Matrix(n, n, 0)
        for i in range(n):
            for j in range(i + 1, n):
                distance = list_of_sensors[i].GetLocation() - list_of_sensors[j].GetLocation()
                distance = (distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2) ** (0.5)
                self.distance_matrix[i, j] = distance
                self.distance_matrix[j, i] = distance

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorInverseDistanceResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return None

    def CalculateValue(self) -> float:
        value = 0.0
        for i, node_i in enumerate(self.model_part.Nodes):
            for j, node_j in enumerate(self.model_part.Nodes):
                theta_ij = node_i.GetValue(KratosDT.SENSOR_STATUS) * node_j.GetValue(KratosDT.SENSOR_STATUS)
                value += exp(-theta_ij * self.distance_matrix[i, j])
        return value

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        gradient = np.zeros(self.distance_matrix.Size1(), dtype=np.float64)
        for k, node_k in enumerate(self.model_part.Nodes):
            theta_k = node_k.GetValue(KratosDT.SENSOR_STATUS)
            for i, node_i in enumerate(self.model_part.Nodes):
                theta_i = node_i.GetValue(KratosDT.SENSOR_STATUS)
                gradient[k] -= 2 * theta_i * self.distance_matrix[i, k] * exp (-theta_i * theta_k * self.distance_matrix[i, k])

        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)
        Kratos.Expression.CArrayExpressionIO.Read(physical_variable_collective_expressions[KratosDT.SENSOR_STATUS].GetContainerExpressions()[0], gradient)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"