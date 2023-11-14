from typing import Optional, Type
import csv

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_static_analysis import SensorSensitivityStaticAnalysis
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import SensorViewUnionType

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DamageDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DamageDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return DamageDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class DamageDetectionResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "adjoint_parameters": {},
            "test_analysis_list": [
                {
                    "primal_analysis_name": "Structure_static",
                    "sensor_measurement_data_file": "measurement_data.csv",
                    "sensor_computed_data_file": "computed_data.csv",
                    "weight": 1.0
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DamageDetectionResponse. [ response name = \"{self.GetName()}\"]")

        # reading test analsis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, str, str, float]]' = []
        for params in parameters["test_analysis_list"]:
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            sensor_measurement_data_file_name = params["sensor_measurement_data_file"].GetString()
            sensor_computed_data_file = params["sensor_computed_data_file"].GetString()
            weight = params["weight"].GetDouble()
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), sensor_measurement_data_file_name, sensor_computed_data_file, weight))

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _, _ in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None

        self.adjoint_analysis = SensorSensitivityStaticAnalysis(self.model, parameters["adjoint_parameters"])

        self.sensor_name_dict: 'dict[str, KratosDT.Sensors.Sensor]' = {}

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()
        self.list_of_sensors = self.adjoint_analysis.GetListOfSensors()
        for sensor in self.list_of_sensors:
            self.sensor_name_dict[sensor.GetName()] = sensor

    def Check(self) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_part)

    def Finalize(self) -> None:
        self.adjoint_analysis.Finalize()

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call DamageDetectionResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        if self.analysis_model_part is None:
            raise RuntimeError("Please call DamageDetectionResponse::Initialize first.")
        return self.analysis_model_part

    def CalculateValue(self) -> float:
        result = 0.0
        for exec_policy, sensor_measurement_data_file_name, sensor_computed_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # first run the primal analysis.
            exec_policy.Execute()
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Computed \"{exec_policy.GetName()}\".")

            # now open the data files generated from the primal analysis.
            csv_measurement_file = open(sensor_measurement_data_file_name, "r")
            csv_computed_file = open(sensor_computed_data_file_name, "r")

            csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
            csv_computed_stream = csv.reader(csv_computed_file, delimiter=",")

            measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)
            computed_name_index, computed_value_index = self.__GetHeaderIndices(csv_computed_stream)

            for measured_row, computed_row in zip(csv_measurement_stream, csv_computed_stream):
                measured_sensor_name = measured_row[measured_name_index].strip()
                computed_sensor_name = computed_row[computed_name_index].strip()
                if measured_sensor_name != computed_sensor_name:
                    raise RuntimeError(f"Mismatching sensor names found [ measured_sensor_name = {measured_sensor_name}, computed_sensor_name = {computed_sensor_name}]")

                sensor = self.__GetSensor(measured_sensor_name)
                measured_value = float(measured_row[measured_value_index])
                computed_value = float(computed_row[computed_value_index])
                result += sensor.GetWeight() * test_case_weight * ((measured_value - computed_value) ** 2) / 2.0

            csv_computed_file.close()
            csv_measurement_file.close()

        return result

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        self.adjoint_analysis.RunSolutionLoop()

        if len(physical_variable_collective_expressions.keys()) > 1:
            raise RuntimeError(f"Currently {self.__class__.__name__} only supports computing gradient w.r.t. one variable only.")

        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        for _, sensor_measurement_data_file_name, sensor_computed_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            csv_measurement_file = open(sensor_measurement_data_file_name, "r")
            csv_computed_file = open(sensor_computed_data_file_name, "r")

            csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
            csv_computed_stream = csv.reader(csv_computed_file, delimiter=",")

            measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)
            computed_name_index, computed_value_index = self.__GetHeaderIndices(csv_computed_stream)

            for physical_variable, _ in merged_model_part_map.items():
                list_of_container_expression = physical_variable_collective_expressions[physical_variable].GetContainerExpressions()
                if len(list_of_container_expression) > 1:
                    raise RuntimeError(f"Currently {self.__class__.__name__} only supports one model part.")

                cexp_gradient = list_of_container_expression[0]
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(cexp_gradient, physical_variable)
                sensor_view_type: 'Optional[Type[SensorViewUnionType]]' = None
                if isinstance(cexp_gradient, Kratos.Expression.NodalExpression):
                    sensor_view_type = KratosDT.Sensors.NodalSensorView
                elif isinstance(cexp_gradient, Kratos.Expression.ConditionExpression):
                    sensor_view_type = KratosDT.Sensors.ConditionSensorView
                elif isinstance(cexp_gradient, Kratos.Expression.ElementExpression):
                    sensor_view_type = KratosDT.Sensors.ElementSensorView
                else:
                    raise RuntimeError("Unsupported type.")

                for measured_row, computed_row in zip(csv_measurement_stream, csv_computed_stream):
                    measured_sensor_name = measured_row[measured_name_index].strip()
                    computed_sensor_name = computed_row[computed_name_index].strip()
                    if measured_sensor_name != computed_sensor_name:
                        raise RuntimeError(f"Mismatching sensor types found [ measured_sensor_name = {measured_sensor_name}, computed_sensor_name = {computed_sensor_name}]")

                    sensor = self.__GetSensor(measured_sensor_name)
                    measured_value = float(measured_row[measured_value_index])
                    computed_value = float(computed_row[computed_value_index])
                    sensor_view = sensor_view_type(sensor, physical_variable.Name() + "_SENSITIVITY")

                    cexp_gradient += sensor_view.GetContainerExpression() * (measured_value - computed_value) * test_case_weight * sensor.GetWeight()

                cexp_gradient.SetExpression(cexp_gradient.Flatten().GetExpression())

            csv_computed_file.close()
            csv_measurement_file.close()

    def __GetSensor(self, sensor_name: str) -> KratosDT.Sensors.Sensor:
        return self.sensor_name_dict[sensor_name]

    def __GetHeaderIndices(self, csv_stream: csv.reader) -> 'tuple[int, int]':
        headers = [s.strip() for s in next(csv_stream)]
        name_index = headers.index("name")
        value_index = headers.index("value")
        return name_index, value_index

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"