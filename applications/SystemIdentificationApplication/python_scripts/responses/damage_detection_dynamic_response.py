from typing import Optional
import csv
from pathlib import Path
import os, glob, math

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosDT
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_dynamic_analysis import SystemIdentificationDynamicAnalysis
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
                    "primal_analysis_name": "Structure_dynamic",
                    "weight": 1.0,
                    "sensor_measurement_hdf5_file": "../damaged_system/hdf5_output/<model_part_name>_T_<step>.h5",
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

        # reading test analysis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, str, float]]' = []
        for params in parameters["test_analysis_list"].values():
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            sensor_measurement_data_file_name = params["sensor_measurement_hdf5_file"].GetString()
            weight = params["weight"].GetDouble()
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), sensor_measurement_data_file_name, weight))
        

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _  in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None
    
        self.adjoint_analysis = SystemIdentificationDynamicAnalysis(self.model, parameters["adjoint_parameters"])

        self.num_of_timesteps = self.NumberOfTimeSteps()

        self.sensor_name_dict: 'dict[str, KratosDT.Sensors.Sensor]' = {}
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS]
    
    def NumberOfTimeSteps(self) -> int:
        return self.adjoint_analysis.NumberOfTimeSteps()

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()

        self.damage_response_function = KratosDT.Responses.MeasurementResidualDynamicResponseFunction(self.p_coefficient)

        if not self.optimization_problem.GetProblemDataContainer()["object"].HasValue(self.sensor_group_name):
            raise RuntimeError(f"The sensor group \"{self.sensor_group_name}\" not found. Followings are available: \n\t" + "\n\t".join(self.optimization_problem.GetProblemDataContainer()["object"].GetSubItems().keys()))

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        self.list_of_sensors = GetSensors(sensor_group_data)
        for sensor in self.list_of_sensors:
            sensor.GetNode().SetValue(KratosDT.SENSOR_MEASURED_VALUE, 0.0)
            sensor.GetNode().SetValue(KratosDT.SENSOR_COMPUTED_VALUE, 0.0)
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
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # first run the primal analysis.
            exec_policy.Execute()
            if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(exec_policy.GetAnalysisModelPart(), f"{exec_policy.GetName()}_COMPUTED_SENSOR_DATA_written"):
                raise RuntimeError(f"{exec_policy.GetName()} has not written down the computed sensor measurements")

            list_of_logged_statuses = KratosOA.OptAppModelPartUtils.GetModelPartStatusLog(exec_policy.GetAnalysisModelPart())
            search_str = f"{exec_policy.GetName()}_COMPUTED_SENSOR_DATA_Path:"
            for logged_status in list_of_logged_statuses:
                pos = logged_status.find(search_str)
                if pos != -1:
                    hdf5_path_data = logged_status[pos+len(search_str):]
                    hdf5_file_name = hdf5_path_data.split(":")[0]
                    hdf5_prefix = hdf5_path_data.split(":")[1]
                    break
            cost_function = []
            for step in range(1, self.num_of_timesteps+1):
                # Reading sensor measurement data
                params = Kratos.Parameters("""{
                    "file_access_mode": "read_only"
                }""")
                measurement_data_file_name_orig = sensor_measurement_data_file_name.replace("<model_part_name>", "SensorModelPart")
                measurement_data_file_name = measurement_data_file_name_orig.replace("<step>", f"{step}")
                print("sensor_measurement_data_file_name is ", measurement_data_file_name)
                params.AddString("file_name", measurement_data_file_name)
                with OpenHDF5File(params, self.model_part) as h5_file:
                    sensor_params = Kratos.Parameters("""{
                        "list_of_variables": ["SENSOR_MEASURED_VALUE"]
                    }""")
                    sensor_params.AddString("prefix", hdf5_prefix)
                    exp = KratosHDF5.HDF5NodalDataValueIO(sensor_params, h5_file)
                    exp.Read(self.model[self.sensor_group_name])

                # Write out data for the sensor measurement data to access later in gradient calculation
                KratosOA.OptAppModelPartUtils.LogModelPartStatus(exec_policy.GetAnalysisModelPart(), f"{exec_policy.GetName()}_MEASURED_SENSOR_DATA_written")
                KratosOA.OptAppModelPartUtils.LogModelPartStatus(exec_policy.GetAnalysisModelPart(), f"{exec_policy.GetName()}_MEASURED_SENSOR_DATA_Path:{measurement_data_file_name_orig}:{hdf5_prefix}")
                
                # Reading sensor computed data
                params = Kratos.Parameters("""{
                    "file_access_mode": "read_only"
                }""")
                params.AddString("file_name", hdf5_file_name.replace("<step>", f"{step}"))
                print("computed file name is ", hdf5_file_name.replace("<step>", f"{step}"))
                with OpenHDF5File(params, self.model_part) as h5_file:
                    sensor_params = Kratos.Parameters("""{
                        "list_of_variables": ["SENSOR_COMPUTED_VALUE"]
                    }""")
                    sensor_params.AddString("prefix", hdf5_prefix)
                    KratosHDF5.HDF5NodalDataValueIO(sensor_params, h5_file).Read(self.model[self.sensor_group_name])
                
                result += test_case_weight * self.damage_response_function.CalculateValue(exec_policy.GetAnalysisModelPart())
                cost_function.append(result)
                print("result ", result)
            Kratos.Logger.PrintInfo(self._GetResponsePrefix(), f"Computed \"{exec_policy.GetName()}\".")
            
        with open("cost_function.dat", "w") as fo:
            fo.write("cost_function")
            for i in range(len(cost_function)):
                fo.write(f"\n{cost_function[i]}")
                
        return result

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute sensitivities for each test scenario
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # Check if primal solution data was written
            if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(exec_policy.GetAnalysisModelPart(), f"{exec_policy.GetName()}_PRIMAL_DATA_written"):
                raise RuntimeError(f"{exec_policy.GetName()} has not written down the primal solution data")
            
            # run a single adjoint for each test scenario
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosDT.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()
            list_of_logged_statuses = KratosOA.OptAppModelPartUtils.GetModelPartStatusLog(exec_policy.GetAnalysisModelPart())
            for log_status in list_of_logged_statuses:
                KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.adjoint_analysis._GetSolver().GetComputingModelPart(), log_status)
            sensitivities = self.adjoint_analysis.CalculateGradient(self.damage_response_function)

            for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                sensitivity_variable = Kratos.KratosGlobals.GetVariable(Kratos.SensitivityUtilities.GetSensitivityVariableName(physical_variable))
                for container_expression in collective_expression.GetContainerExpressions():
                    container_expression.SetExpression((container_expression.GetExpression() - sensitivities[sensitivity_variable].GetExpression() * test_case_weight))
                    container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

                    with open(f"{Kratos.SensitivityUtilities.GetSensitivityVariableName(physical_variable)}_adjoint_sensitivites.dat", "w") as fo:
                        fo.write("element_id,gradient")
                        grad = container_expression.Evaluate()
                        for i in range(len(grad)):
                            fo.write(f"\n{i+1},{grad[i]}")

    def _GetResponsePrefix(self) -> str:
        return "DamageDetectionDynamicResponse"

    def __str__(self) -> str:
        return f"Response [type = {self._GetResponsePrefix()}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"