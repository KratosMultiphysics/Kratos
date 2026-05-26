from typing import Optional
import csv
import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_rans_analysis import SystemIdentificationRANSAnalysis
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DamageDetectionRANSResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DamageDetectionRANSResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return DamageDetectionRANSResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class DamageDetectionRANSResponse(ResponseFunction):
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
                    "sensor_normalization_settings": {
                        "type": "none"
                    },
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
            raise RuntimeError(f"No model parts were provided for {self._GetResponsePrefix()}. [ response name = \"{self.GetName()}\"]")

        # reading test analsis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, DataIO, str, float]]' = []
        for params in parameters["test_analysis_list"].values():
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            sensor_measurement_data_file_name = params["sensor_measurement_csv_file"].GetString()
            sensor_normalization_settings = params["sensor_normalization_settings"]
            weight = params["weight"].GetDouble()
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), sensor_measurement_data_file_name, sensor_normalization_settings, weight))

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _, _ in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None

        self.adjoint_analysis = SystemIdentificationRANSAnalysis(self.model, parameters["adjoint_parameters"], optimization_problem, self.p_coefficient)
        self.model[self.model_part_operation.GetModelPartFullName()].AddNodalSolutionStepVariable(Kratos.VELOCITY)

        self.sensor_name_dict: 'dict[str, KratosSI.Sensors.Sensor]' = {}
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.VELOCITY_Z, 
                KratosRANS.TURBULENT_KINETIC_ENERGY, KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 
                KratosRANS.LOW_FIDELITY_INLET_1_SENSITIVITY, KratosRANS.LOW_FIDELITY_INLET_2_SENSITIVITY, KratosRANS.LOW_FIDELITY_INLET_3_SENSITIVITY]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()

        self.damage_response_function = KratosSI.Responses.MeasurementResidualResponseFunction(self.p_coefficient)

        if not self.optimization_problem.GetProblemDataContainer()["object"].HasValue(self.sensor_group_name):
            raise RuntimeError(f"The sensor group \"{self.sensor_group_name}\" not found. Followings are available: \n\t" + "\n\t".join(self.optimization_problem.GetProblemDataContainer()["object"].GetSubItems().keys()))

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        self.list_of_sensors = GetSensors(sensor_group_data)
        self.adjoint_analysis.SetSensors(self.sensor_group_name, self.list_of_sensors)
        for sensor in self.list_of_sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, 0.0)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
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
            raise RuntimeError(f"Please call {self._GetResponsePrefix()}::Initialize first.")
        return self.analysis_model_part

    def CalculateValue(self) -> float:
        result = 0.0
        for exec_policy, sensor_measurement_data_file_name, sensor_normalization_settings, test_case_weight in self.list_of_test_analysis_data:
            # first run the primal analysis.
            exec_policy.Execute()

            self.__SetSensorMeasuredValue(sensor_measurement_data_file_name)
            self.__SetSensorNormalizationFactor(sensor_normalization_settings)

            result += test_case_weight * self.damage_response_function.CalculateValue(exec_policy.GetAnalysisModelPart())
            Kratos.Logger.PrintInfo(self._GetResponsePrefix(), f"Computed \"{exec_policy.GetName()}\".")

        return result

    def CalculateGradient(self, physical_variable_gradient_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        # make everything zeros
        for physical_variable, cta in physical_variable_gradient_map.items():
            cta.data[:] = 0.0
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta, perform_store_data_recursively=False, copy=False).StoreData()

        # now compute sensitivities for each test scenario
        for exec_policy, sensor_measurement_data_file_name, sensor_normalization_settings, test_case_weight in self.list_of_test_analysis_data:
            # read and replace the measurement data for each test scenario
            self.__SetSensorMeasuredValue(sensor_measurement_data_file_name)
            self.__SetSensorNormalizationFactor(sensor_normalization_settings)

            self.adjoint_analysis.CalculateGradient(self.damage_response_function)

            for physical_variable, cta in physical_variable_gradient_map.items():
                # Temporary for low fidelity. Will provide Inlet Velocity gradients to control, who will aply chain rule during mapping.
                if physical_variable in [KratosRANS.LOW_FIDELITY_INLET_1_SENSITIVITY, KratosRANS.LOW_FIDELITY_INLET_2_SENSITIVITY, KratosRANS.LOW_FIDELITY_INLET_3_SENSITIVITY]:
                    physical_variable = Kratos.VELOCITY_X

                sensitivity_variable = Kratos.KratosGlobals.GetVariable(Kratos.SensitivityUtilities.GetSensitivityVariableName(physical_variable))
                for ta in cta.GetTensorAdaptors():
                    current_gradient = ta.Clone()
                    self.adjoint_analysis.GetGradient(sensitivity_variable, current_gradient)
                    ta.data[:] += current_gradient.data[:] * test_case_weight
                Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta, perform_collect_data_recursively=False, copy=False).CollectData()

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
                if measured_sensor_name in self.sensor_name_dict:
                    self.__GetSensor(measured_sensor_name).GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, measured_value)

    def __SetSensorNormalizationFactor(self, sensor_normalization_settings: Kratos.Parameters) -> None:
        sensor_normalization_type = sensor_normalization_settings["type"].GetString()
        if sensor_normalization_type == "none":
            normalization_factor = 1.0
            for sensor in self.list_of_sensors:
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        elif sensor_normalization_type == "maximum_measured_value":
            # Find groups and their maximums
            sensor_type_maximum_dict = {}
            for sensor in self.list_of_sensors:
                sensor_type = sensor.GetSensorParameters()["type"].GetString()
                if not sensor_type in sensor_type_maximum_dict:
                    sensor_type_maximum_dict[sensor_type] = abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE))
                else:
                    sensor_type_maximum_dict[sensor_type] = max(sensor_type_maximum_dict[sensor_type], abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)))
            # Set normalization factor
            for sensor in self.list_of_sensors:
                sensor_type = sensor.GetSensorParameters()["type"].GetString()
                normalization_factor = sensor_type_maximum_dict[sensor_type]
                if math.isclose(normalization_factor, 0.0, abs_tol=1e-16):
                    raise RuntimeError(f"The maximum measured value for sensor type \"{sensor_type}\" normalization is approximately zero ({normalization_factor}). Cannot normalize.")
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        elif sensor_normalization_type == "average_measured_value":
            normalization_factor = sum(abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)) for sensor in self.list_of_sensors) / len(self.list_of_sensors)
            if math.isclose(normalization_factor, 0.0, abs_tol=1e-16):
                raise RuntimeError(f"The average measured value for sensor normalization is approximately zero ({normalization_factor}). Cannot normalize.")
            for sensor in self.list_of_sensors:
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        elif sensor_normalization_type == "local_measured_value":
            for sensor in self.list_of_sensors:
                normalization_factor = abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE))
                if math.isclose(normalization_factor, 0.0, abs_tol=1e-16):
                    raise RuntimeError(f"The local measured value for sensor \"{sensor.GetName()}\" normalization is approximately zero ({normalization_factor}). Cannot normalize.")
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        elif sensor_normalization_type == "local_maximum_measured_value":
            if not sensor_normalization_settings.Has("epsilon"):
                raise RuntimeError(f"""The sensor normalization settings for type \"local_maximum_measured_value\"
                                   must have an \"epsilon\" value.""")
            epsilon = sensor_normalization_settings["epsilon"].GetDouble()
            max_measured_value = max(abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)) for sensor in self.list_of_sensors)
            for sensor in self.list_of_sensors:
                normalization_factor = max(abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)), epsilon * max_measured_value)
                if math.isclose(normalization_factor, 0.0, abs_tol=1e-16):
                    raise RuntimeError(f"The local maximum measured value for sensor \"{sensor.GetName()}\" normalization is approximately zero ({normalization_factor}). Cannot normalize.")
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        elif sensor_normalization_type == "maximum_velocity_local_dynamic_pressure_measured_value":
            maximum_velocity = 0.0
            for sensor in self.list_of_sensors:
                if sensor.GetSensorParameters()["type"].GetString() == "velocity_sensor":
                    maximum_velocity = max(maximum_velocity, abs(sensor.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE)))
            if math.isclose(maximum_velocity, 0.0, abs_tol=1e-16):
                raise RuntimeError(f"The maximum velocity measured value for sensor normalization is approximately zero ({maximum_velocity}). Cannot normalize.")

            # Set normalization factor
            for sensor in self.list_of_sensors:
                sensor_type = sensor.GetSensorParameters()["type"].GetString()
                if sensor_type == "velocity_sensor":
                    normalization_factor = maximum_velocity
                elif sensor_type == "pressure_sensor":
                    # Find corresponding velocity sensor for this pressure sensor
                    found_sensor = False
                    for sensor_1 in self.list_of_sensors:
                        if sensor_1.GetSensorParameters()["type"].GetString() == "velocity_sensor" and math.isclose(sensor_1.GetNode().X, sensor.GetNode().X) and math.isclose(sensor_1.GetNode().Y, sensor.GetNode().Y) and math.isclose(sensor_1.GetNode().Z, sensor.GetNode().Z):
                            velocity_at_pressure_sensor = abs(sensor_1.GetNode().GetValue(KratosSI.SENSOR_MEASURED_VALUE))
                            found_sensor = True
                            break
                    if not found_sensor:
                        raise RuntimeError(f"""Could not find a velocity sensor at the same location as pressure sensor \"{sensor.GetName()}\" for normalization.
                                           Make sure that there is a velocity sensor at the same location as each pressure sensor for this normalization type.""")
                    normalization_factor = 0.5 * 1.0* velocity_at_pressure_sensor * velocity_at_pressure_sensor

                if math.isclose(normalization_factor, 0.0, abs_tol=1e-16):
                    raise RuntimeError(f"The maximum measured value for sensor type \"{sensor_type}\" normalization is approximately zero ({normalization_factor}). Cannot normalize.")
                sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, normalization_factor)
            return
        else:
            raise RuntimeError(f"""Unsupported sensor normalization type \"{sensor_normalization_type}\".
                               Supported types are: none, maximum_measured_value, average_measured_value,
                               local_measured_value, local_maximum_measured_value.""")


    def _GetResponsePrefix(self) -> str:
        return "DamageDetectionRANSResponse"

    def __str__(self) -> str:
        return f"Response [type = {self._GetResponsePrefix()}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"