from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import HasSensorStatusControlUpdater
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensorStatusControlUpdater
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import AddSensorStatusControlUpdater
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

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
            "sensor_mask_name"           : "",
            "echo_level"                 : 0,
            "beta"                       : 4,
            "allowed_dissimilarity"      : 1.0,
            "leaf_max_size"              : 10,
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.leaf_max_size = parameters["leaf_max_size"].GetInt()
        self.echo_level = parameters["echo_level"].GetInt()

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for SensorLocalizationResponse. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.optimization_problem = optimization_problem
        self.beta = parameters["beta"].GetDouble()
        self.allowed_dissimilarity = parameters["allowed_dissimilarity"].GetDouble()

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        if not HasSensorStatusControlUpdater(f"sensor_mask_status_{self.sensor_mask_name}", self.optimization_problem):
            list_of_sensors = GetSensors(self.optimization_problem)
            sensor_mask_status = KratosSI.SensorMaskStatus(self.model_part, [sensor.GetContainerExpression(self.sensor_mask_name) for sensor in list_of_sensors], self.echo_level)
            AddSensorStatusControlUpdater(f"sensor_mask_status_{self.sensor_mask_name}", sensor_mask_status, self.optimization_problem)
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created sensor mask status for \"{self.sensor_mask_name}\".")

        if not HasSensorStatusControlUpdater(f"sensor_mask_status_kd_tree_{self.sensor_mask_name}", self.optimization_problem):
            sensor_mask_status: KratosSI.SensorMaskStatus = GetSensorStatusControlUpdater(f"sensor_mask_status_{self.sensor_mask_name}", self.optimization_problem)
            sensor_mask_status_kd_tree = KratosSI.SensorMaskStatusKDTree(sensor_mask_status, self.leaf_max_size, self.echo_level)
            AddSensorStatusControlUpdater(f"sensor_mask_status_kd_tree_{self.sensor_mask_name}", sensor_mask_status_kd_tree, self.optimization_problem)
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created sensor mask status kd tree for \"{self.sensor_mask_name}\".")

        sensor_mask_status_kd_tree: KratosSI.SensorMaskStatusKDTree = GetSensorStatusControlUpdater(f"sensor_mask_status_kd_tree_{self.sensor_mask_name}", self.optimization_problem)
        self.utils = KratosSI.SensorLocalizationResponseUtils(sensor_mask_status_kd_tree, self.beta, self.allowed_dissimilarity)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call SensorLocalizationResponse::Initialize first.")
        return self.model_part

    def CalculateValue(self) -> float:
        return self.utils.CalculateValue()

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            if physical_variable == KratosSI.SENSOR_STATUS:
                collective_expression.GetContainerExpressions()[0].SetExpression(self.utils.CalculateGradient().GetExpression())
            else:
                raise RuntimeError("Unsupported physical variable.")

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"