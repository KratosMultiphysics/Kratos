from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorLocalizationResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorLocalizationResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorLocalizationResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class SensorLocalizationResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "sensor_group_name"    : "",
            "sensor_mask_name"     : "",
            "echo_level"           : 0,
            "beta"                 : 4,
            "allowed_dissimilarity": 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.echo_level = parameters["echo_level"].GetInt()
        self.beta = parameters["beta"].GetDouble()
        self.allowed_dissimilarity = parameters["allowed_dissimilarity"].GetDouble()

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.sensor_group_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.mask_status_kd_tree: 'Optional[KratosSI.SensorMaskStatusKDTree]' = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for controller in GetMaskStatusControllers(sensor_group_data, self.sensor_mask_name):
            if isinstance(controller, KratosSI.SensorMaskStatusKDTree):
                self.mask_status_kd_tree = controller

        self.utils = KratosSI.Responses.SensorLocalizationResponseUtils(self.mask_status_kd_tree, self.beta, self.allowed_dissimilarity)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorCoverageResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return self.utils.CalculateValue()

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        physical_variable_collective_expressions[KratosSI.SENSOR_STATUS].GetContainerExpressions()[0].SetExpression(self.utils.CalculateGradient().GetExpression())

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"