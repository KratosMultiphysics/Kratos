from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorIsolationResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorIsolationResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorIsolationResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class SensorIsolationResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "isolation_distance": 0.0,
            "beta"              : 300.0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for SensorIsolationResponse. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.isolation_distance = parameters["isolation_distance"].GetDouble()
        self.beta = parameters["beta"].GetDouble()

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosDT.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        self.utils = KratosDT.SensorIsolationResponseUtils(self.model_part, int(1e+6), self.isolation_distance, self.beta)
        self.utils.Initialize()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorIsolationResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return None

    def CalculateValue(self) -> float:
        return self.utils.CalculateValue()

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            if physical_variable == KratosDT.SENSOR_STATUS:
                collective_expression.GetContainerExpressions()[0].SetExpression(self.utils.CalculateGradient().GetExpression())
            else:
                raise RuntimeError("Unsupported physical variable.")

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"