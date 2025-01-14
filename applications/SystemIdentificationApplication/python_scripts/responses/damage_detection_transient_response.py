from typing import Optional
import csv
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_transient_analysis import SystemIdentificationTransientAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"DamageDetectionTransientResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"DamageDetectionTransientResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return DamageDetectionTransientResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class DamageDetectionTransientResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "adjoint_parameters"         : {},
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "test_analysis_list": [
                {
                    "primal_analysis_name": "Structure_transient",
                    "primal_hdf5_data_file"      : "hdf5_data_<time>.h5",
                    "sensor_measurement_csv_file": "measurement_data_<step>.csv",
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

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for DamageDetectionTransientResponse. [ response name = \"{self.GetName()}\"]")

        # reading test analsis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, str, float]]' = []
        for params in parameters["test_analysis_list"]: # WARNING: deprecated, use "values" method
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            sensor_measurement_data_file_name = params["sensor_measurement_csv_file"].GetString()
            primal_hdf5_data_file = params["primal_hdf5_data_file"].GetString()
            weight = params["weight"].GetDouble()
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), sensor_measurement_data_file_name, primal_hdf5_data_file, weight))

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _, _ in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None

        self.adjoint_analysis = SystemIdentificationTransientAnalysis(self.model, parameters["adjoint_parameters"])
        self.adjoint_parameters = parameters["adjoint_parameters"]

        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()

        self.adjoint_analysis.Initialize()
        self.list_of_sensors = self.adjoint_analysis.GetListOfSensors()

        ComponentDataView(self, self.optimization_problem).GetUnBufferedData().SetValue("sensors", self.list_of_sensors)

    def Check(self) -> None:
        pass
    
    def Finalize(self) -> None:
        # self.adjoint_analysis.Finalize()
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.analysis_model_part is None:
            raise RuntimeError("Please call DamageDetectionTransientResponse::Initialize first.")
        return self.analysis_model_part

    def CalculateValue(self) -> float:
        value = 0.0
        for exec_policy, sensor_measurement_data_file_name, primal_hdf5_data_file, test_case_weight in self.list_of_test_analysis_data:
            # run the primal analysis first
            exec_policy.Execute()

            # get time stepping data from primal analysis
            start_time, end_time, time_step_settings = exec_policy.GetAnalysisTimeSteppingData()  
            if time_step_settings.Has("time_step_table"):
                raise Exception("DamageDetectionTransientResponse:: there is currently no variable time stepping possible!")
            time_step = time_step_settings["time_step"].GetDouble()

            # initialize time loop
            time = start_time
            step = 0
            while time < end_time:
                time += time_step
                step +=1

                # set measurement data
                self.__SetSensorMeasuredValue(sensor_measurement_data_file_name.replace("<step>", f"{step:06d}"))

                # set simulation data from primal case
                params = Kratos.Parameters("""{
                    "file_name" : "",
                    "file_access_mode" : "read_only"    
                }""")
                params["file_name"].SetString(primal_hdf5_data_file.replace("<time>", f"{time:0.4f}"))
                data_params = Kratos.Parameters("""{
                    "prefix": "/ResultsData",
                    "list_of_variables": ["ALL_VARIABLES_FROM_FILE"]
                }""")
                with OpenHDF5File(params, exec_policy.GetAnalysisModelPart()) as h5_file:
                    KratosHDF5.HDF5NodalSolutionStepDataIO(data_params, h5_file).ReadNodalResults(exec_policy.GetAnalysisModelPart(), 0)

                # calculate and add response value contribution
                value_contribution = self.adjoint_analysis.GetResponseFunction().CalculateValue(exec_policy.GetAnalysisModelPart()) * time_step           
                value += value_contribution * test_case_weight
            
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Computed \"{exec_policy.GetName()}\".")

        return value

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # clear container expressions
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)
        
        for exec_policy, sensor_measurement_data_file_name, _, test_case_weight in self.list_of_test_analysis_data:

            self.model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            # reset adjoint model part and adjoint analysis
            KratosOA.OptimizationUtils.ResetModelPartNodalSolutionStepData(self.model_part)
            self.model_part.ProcessInfo[Kratos.TIME] = self.adjoint_parameters["problem_data"]["start_time"].GetDouble()
            self.model_part.ProcessInfo[Kratos.STEP] = 0
            # self.adjoint_analysis = SystemIdentificationTransientAnalysis(self.model, self.adjoint_parameters.Clone())
            self.adjoint_analysis.Initialize()

            # set measurement data file name
            self.adjoint_analysis._SetSensorMeasurementDataFileName(sensor_measurement_data_file_name)

            # run a single adjoint for each test scenario
            self.adjoint_analysis.RunSolutionLoop()
            self.adjoint_analysis.Finalize()

            for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                sensitivity_variable = Kratos.KratosGlobals.GetVariable(f"{physical_variable.Name()}_SENSITIVITY")
                for container_expression in collective_expression.GetContainerExpressions():
                    sensitivities = self.adjoint_analysis.GetSensitivities(container_expression.GetModelPart())
                    container_expression.SetExpression((container_expression.GetExpression() + sensitivities[sensitivity_variable].GetExpression() * test_case_weight))
                    container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

        # for physical_variable, collective_expression in physical_variable_collective_expressions.items():
        #     sensitivity_variable = Kratos.KratosGlobals.GetVariable(f"{physical_variable.Name()}_SENSITIVITY")

        #     if sensitivity_variable ==KratosSA.YOUNG_MODULUS_SENSITIVITY:
        #         for container_expression in collective_expression.GetContainerExpressions():
        #             import h5py
        #             with h5py.File("sensor_data.h5", "r") as h5_output:
        #                 numpy_array = h5_output["/data"][:]
        #             sensitivities = Kratos.Expression.ElementExpression(self.model["all_nodes_elements_model_part"])
        #             Kratos.Expression.CArrayExpressionIO.Read(sensitivities, numpy_array)
        #             container_expression.SetExpression((container_expression.GetExpression() + sensitivities.GetExpression()))
        #             container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())
        #     else:
        #         raise RuntimeError("safsadda")

    def __GetSensor(self, sensor_name: str) -> KratosSI.Sensors.Sensor:
        return self.adjoint_analysis.GetSensorNameDictionary()[sensor_name]

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
                self.__GetSensor(measured_sensor_name).SetValue(KratosSI.SENSOR_MEASURED_VALUE, measured_value)

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"
