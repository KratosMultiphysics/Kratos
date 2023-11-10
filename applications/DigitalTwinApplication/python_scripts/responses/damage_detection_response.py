from typing import Optional, Type, Union

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
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import AddSensorSpecificationVariableData
from KratosMultiphysics.DigitalTwinApplication.sensor_specification_solvers.sensor_specification_static_analysis import SensorSpecificationStaticAnalysis
import pandas

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
                ],
            "list_of_specifications": [
                {
                    "name": "disp_x",
                    "id": 1,
                    "location":[0.0, 0.0, 0.0],
                    "variable_data": {}
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

        self.adjoint_analysis = SensorSpecificationStaticAnalysis(self.model, parameters["adjoint_parameters"])

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE, Kratos.DENSITY, Kratos.THICKNESS, KratosOA.CROSS_AREA]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()
        self.list_of_specifications = self.adjoint_analysis.GetListOfSpecifications()

    def Check(self) -> None:
        KratosOA.ResponseUtils.MassResponseUtils.Check(self.model_part)

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call DamageDetectionResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> None:
        return None

    def CalculateValue(self) -> float:
        result = 0.0
        for _, sensor_measurement_data_file_name, sensor_computed_data_file, test_case_weight in self.list_of_test_analysis_data:
            data = pandas.read_csv(sensor_measurement_data_file_name, delimiter=";")
            measured_names = data["name"]
            measured_ids = data["#"]
            measured_values = data["value"].to_numpy()

            data = pandas.read_csv(sensor_computed_data_file, delimiter=";")
            computed_values = data["value"]

            for sensor_name, sensor_id, measured_value, computed_value in zip(measured_names, measured_ids, measured_values, computed_values):
                sensor_weight = self.__GetSensorWeight(self.__GetSpecification(sensor_name, sensor_id))
                result += ((measured_value - computed_value) ** 2) / (2.0 * sensor_weight * test_case_weight)
        return result

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        self.adjoint_analysis.Run()

        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        for _, sensor_measurement_data_file_name, sensor_computed_data_file, test_case_weight in self.list_of_test_analysis_data:
            data = pandas.read_csv(sensor_measurement_data_file_name, delimiter=";")
            measured_names = data["name"]
            measured_ids = data["#"]
            measured_values = data["value"].to_numpy()

            data = pandas.read_csv(sensor_computed_data_file, delimiter=";")
            computed_values = data["value"]

            for physical_variable, merged_model_part in merged_model_part_map.items():
                list_of_container_expression = physical_variable_collective_expressions[physical_variable].GetContainerExpressions()
                if len(list_of_container_expression) > 0:
                    raise RuntimeError(f"Currently it only supports one model part.")

                cexp_gradient = list_of_container_expression[0]
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(cexp_gradient, physical_variable)
                specifcation_view_type: 'Union[Type[KratosDT.Sensors.NodalSensorSpecificationView], Type[KratosDT.Sensors.ConditionSensorSpecificationView], Type[KratosDT.Sensors.ElementSensorSpecificationView]]' = None
                if isinstance(cexp_gradient, Kratos.Expression.NodalExpression):
                    specifcation_view_type = KratosDT.Sensors.NodalSensorSpecificationView
                elif isinstance(cexp_gradient, Kratos.Expression.ConditionExpression):
                    specifcation_view_type = KratosDT.Sensors.ConditionSensorSpecificationView
                elif isinstance(cexp_gradient, Kratos.Expression.ElementExpression):
                    specifcation_view_type = KratosDT.Sensors.ElementSensorSpecificationView
                else:
                    raise RuntimeError("Unsupported type.")

                for sensor_name, sensor_id, measured_value, computed_value in zip(measured_names, measured_ids, measured_values, computed_values):
                    specification = self.__GetSpecification(sensor_name, sensor_id)
                    specicaition_view = specifcation_view_type(specification, physical_variable.Name() + "_SENSITIVITY")

                    sensor_weight = self.__GetSensorWeight(specification)
                    cexp_gradient += specicaition_view.GetContainerExpression() * (measured_value - computed_value) / (sensor_weight * test_case_weight)

                cexp_gradient.SetExpression(cexp_gradient.Flatten().GetExpression())


    def __GetSpecification(self, sensor_name: str, sensor_id: int) -> KratosDT.Sensors.SensorSpecification:
        found_specfication = False
        for specification in self.list_of_specifications:
            if specification.GetName() == sensor_name and specification.Id == sensor_id:
                found_specfication = True
                break
        if not found_specfication:
            raise RuntimeError(f"The sensor specification for {sensor_name} with {sensor_id} not found.")
        return specification

    def __GetSensorWeight(self, specification: KratosDT.Sensors.SensorSpecification) -> float:
        return specification.GetValue(KratosDT.SENSOR_WEIGHT)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"