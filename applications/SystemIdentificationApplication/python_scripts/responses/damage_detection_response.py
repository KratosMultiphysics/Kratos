from typing import Optional
import csv
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_static_analysis import SystemIdentificationStaticAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

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
            "sensor_group_name"          : "",
            "p_coefficient"              : 1,
            "adjoint_parameters"         : {},
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "test_analysis_list": [
                {
                    "primal_analysis_name": "Structure_static",
                    "sensor_measurement_csv_file": "measurement_data.csv",
                    "weight": 1.0,
                    "variable_io_settings": {
                        "model_part_name"           : "Structure",
                        "nodal_hist_variable_names"    : [],
                        "nodal_non_hist_variable_names": [],
                        "condition_variable_names"     : [],
                        "element_variable_names"       : [],
                        "condition_property_names"     : [],
                        "element_property_names"       : []
                    }
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.sensor_group_name = parameters["sensor_group_name"].GetString()

        if self.sensor_group_name == "":
            raise RuntimeError(f"The sensor group name cannot be empty.")

        self.p_coefficient = parameters["p_coefficient"].GetDouble()

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DamageDetectionResponse. [ response name = \"{self.GetName()}\"]")

        # reading test analsis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, DataIO, str, float]]' = []
        for params in parameters["test_analysis_list"]:
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            sensor_measurement_data_file_name = params["sensor_measurement_csv_file"].GetString()
            weight = params["weight"].GetDouble()
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), sensor_measurement_data_file_name, weight))

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _ in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None

        self.adjoint_analysis = SystemIdentificationStaticAnalysis(self.model, parameters["adjoint_parameters"])

        self.sensor_name_dict: 'dict[str, KratosSI.Sensors.Sensor]' = {}
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()

        self.damage_response_function = KratosSI.Responses.MeasurementResidualResponseFunction(self.p_coefficient)

        if not self.optimization_problem.GetProblemDataContainer()["object"].HasValue(self.sensor_group_name):
            raise RuntimeError(f"The sensor group \"{self.sensor_group_name}\" not found. Followings are available: \n\t" + "\n\t".join(self.optimization_problem.GetProblemDataContainer()["object"].GetSubItems().keys()))

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        self.list_of_sensors = GetSensors(sensor_group_data)
        for sensor in self.list_of_sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, 0.0)
            self.damage_response_function.AddSensor(sensor)

        self.damage_response_function.Initialize()

        for sensor in self.list_of_sensors:
            self.sensor_name_dict[sensor.GetName()] = sensor

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.adjoint_analysis.Finalize()

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.analysis_model_part is None:
            raise RuntimeError("Please call DamageDetectionResponse::Initialize first.")
        return self.analysis_model_part

    def CalculateValue(self) -> float:
        result = 0.0
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # first run the primal analysis.
            exec_policy.Execute()

            self.__SetSensorMeasuredValue(sensor_measurement_data_file_name)

            result += test_case_weight * self.damage_response_function.CalculateValue(exec_policy.GetAnalysisModelPart())
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Computed \"{exec_policy.GetName()}\".")

        return result

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute sensitivities for each test scenario
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # read and replace the measurement data for each test scenario
            self.__SetSensorMeasuredValue(sensor_measurement_data_file_name)

            # run a single adjoint for each test scenario
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()
            sensitivities = self.adjoint_analysis.CalculateGradient(self.damage_response_function)

            for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                for container_expression in collective_expression.GetContainerExpressions():
                    sensitivity_variable = Kratos.KratosGlobals.GetVariable(f"{physical_variable.Name()}_SENSITIVITY")
                    container_expression.SetExpression((container_expression.GetExpression() - sensitivities[sensitivity_variable].GetExpression() * test_case_weight))
                    container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

    def __GetSensor(self, sensor_name: str) -> KratosSI.Sensors.Sensor:
        return self.sensor_name_dict[sensor_name]

    def __GetHeaderIndices(self, csv_stream: csv.reader) -> 'tuple[int, int]':
        headers = [s.strip() for s in next(csv_stream)]
        name_index = headers.index("name")
        value_index = headers.index("value")
        return name_index, value_index

    def __SetSensorMeasuredValue(self, sensor_measurement_data_file_name: str) -> None:
        with open(sensor_measurement_data_file_name, "r") as csv_measurement_file:
            csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
            measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)

            for measured_row in csv_measurement_stream:
                measured_sensor_name = measured_row[measured_name_index].strip()
                measured_value = float(measured_row[measured_value_index])
                self.__GetSensor(measured_sensor_name).GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, measured_value)


    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"