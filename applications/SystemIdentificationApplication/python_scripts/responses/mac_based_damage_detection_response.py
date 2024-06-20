from typing import Optional
import csv
from pathlib import Path
import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.system_identification_static_analysis import SystemIdentificationStaticAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionDataLocation

from KratosMultiphysics.SystemIdentificationApplication import MacResponseFunctionUtility

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"MacBasedDamageDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MacBasedDamageDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MacBasedDamageDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class MacBasedDamageDetectionResponse(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            
            "evaluated_model_part_names" : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "test_analysis_list": [
                {
                    "primal_analysis_name": "Structure_eigenvector",
                    "eigen_measurement_h5_file": "../damaged_system/measurement_eigenvector_data.h5",
                    "weight": 1.0,
                    "prefix" : "/EigenResults",
                    "variable_io_settings": {
                        "model_part_name"           : "Structure",
                        "nodal_hist_variable_names"    : [],
                        "nodal_non_hist_variable_names": [],
                        "condition_variable_names"     : [],
                        "element_variable_names"       : [],
                        "condition_property_names"     : [],
                        "element_property_names"       : []
                    },
                    "sensitivity_settings" :{
                        "gradient_mode"          : "semi_analytic",
                        "step_size"              : 1e-8,
                        "traced_eigenfrequencies"  : [1],
                        "weighting_method"       : "linear_scaling",
                        "weighting_factors"      : [1.0],
                        "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
                        "element_data_value_sensitivity_variables": [
                            "YOUNG_MODULUS"
                        ],
                        
                        "nodal_solution_step_sensitivity_variables" : [],
                        "condition_data_value_sensitivity_variables" : [],
                        "build_mode": "static"
                    }
                }
            ]        
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
             raise RuntimeError(f"No model parts were provided for MacBasedDamageDetectionResponse. [ response name = \"{self.GetName()}\"]")

        # reading test analysis list
        self.list_of_test_analysis_data: 'list[tuple[ExecutionPolicyDecorator, DataIO, str, float, str, Parameters]]' = []
        for params in parameters["test_analysis_list"]:
            params.ValidateAndAssignDefaults(default_settings["test_analysis_list"][0])
            primal_analysis_name = params["primal_analysis_name"].GetString()
            eigen_measurement_data_file_name = params["eigen_measurement_h5_file"].GetString()
            weight = params["weight"].GetDouble()
            prefix = params["prefix"].GetString()
            sensitivity_settings = params["sensitivity_settings"]
            self.list_of_test_analysis_data.append((optimization_problem.GetExecutionPolicy(primal_analysis_name), eigen_measurement_data_file_name, weight, prefix, sensitivity_settings))
        
        
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None
        self.sensitivity_model_part: Optional[Kratos.ModelPart] = None
        
        self.analysis_model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_test_analysis_{self.GetName()}", [exec.GetAnalysisModelPart().FullName() for exec, _, _, _, _ in self.list_of_test_analysis_data], False)
        self.analysis_model_part: Optional[Kratos.ModelPart] = None
        
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()
        self.evaluated_model_part = self.model_part_operation.GetModelPart()
        self.analysis_model_part = self.analysis_model_part_operation.GetModelPart()
        
        # sensitivity_model_part_name = self.parameters["test_analysis_list"]["sensitivity_settings"]["sensitivity_model_part_name"].GetString();
        # if (sensitivity_model_part_name != "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART"):
        #     self.sensitivity_model_part = self.model.GetModelPart(sensitivity_model_part_name)
        # else:
        self.sensitivity_model_part = self.model_part

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        #self.adjoint_analysis.Finalize()
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call MacBasedDamageDetectionResponse::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        if self.analysis_model_part is None:
            raise RuntimeError("Please call MacBasedDamageDetectionResponse::Initialize first.")
        return self.analysis_model_part
    
    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        if self.sensitivity_model_part is None:
            raise RuntimeError("Please call MacBasedDamageDetectionResponse::Initialize first.")
        return self.sensitivity_model_part

    def CalculateValue(self) -> float:
        result = 0.0
        for exec_policy, eigen_measurement_data_file_name, test_case_weight, prefix, sensitivity_setting in self.list_of_test_analysis_data:
            # first run the primal analysis.
            exec_policy.Execute()
            self.evaluated_model_part.SetValue(KratosSI.COMPUTED_EIGENVALUE_VECTOR, self.evaluated_model_part.ProcessInfo[KratosSI.EIGENVALUE_VECTOR])
            for node in self.evaluated_model_part.Nodes:
                node: Kratos.Node
                node.SetValue(KratosSI.COMPUTED_EIGENVECTOR_MATRIX, node.GetValue(KratosSI.EIGENVECTOR_MATRIX))
                #print(node.GetValue(KratosSI.COMPUTED_EIGENVECTOR_MATRIX))

            self.__SetMeasuredEigenvector(eigen_measurement_data_file_name, prefix)

            print(self.evaluated_model_part.GetValue(KratosSI.MEASURED_EIGENVALUE_VECTOR))
            print(self.evaluated_model_part.GetValue(KratosSI.COMPUTED_EIGENVALUE_VECTOR))
            
            #print(self.evaluated_model_part)
            #print(sensitivity_setting)
            utility = MacResponseFunctionUtility(self.evaluated_model_part, sensitivity_setting)

            result += test_case_weight * utility.CalculateValue()
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Computed \"{exec_policy.GetName()}\".")

        print("result is ",result)
        return result

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute sensitivities for each test scenario
        for exec_policy, eigen_measurement_data_file_name, test_case_weight, prefix, sensitivity_setting in self.list_of_test_analysis_data:
            # read and replace the measurement data for each test scenario
            self.__SetMeasuredEigenvector(eigen_measurement_data_file_name, prefix)

            # # run a single adjoint for each test scenario
            # self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            # self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()
            # self.adjoint_analysis.RunSolutionLoop()
            utility = MacResponseFunctionUtility(self.evaluated_model_part, sensitivity_setting)
            utility.CalculateGradient()
            sensitivities = self.__GetSensitivities(sensitivity_setting)
            #sensitivities = self.adjoint_analysis.GetSensitivities()

            for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                sensitivity_variable = Kratos.KratosGlobals.GetVariable(f"{physical_variable.Name()}_SENSITIVITY")
                for container_expression in collective_expression.GetContainerExpressions():
                    container_expression.SetExpression((container_expression.GetExpression() - sensitivities[sensitivity_variable].GetExpression() * test_case_weight))
                    container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

    # def __GetSensor(self, sensor_name: str) -> KratosSI.Sensors.Sensor:
    #     return self.sensor_name_dict[sensor_name]

    def __GetHeaderIndices(self, csv_stream: csv.reader) -> 'tuple[int, int]':
        headers = [s.strip() for s in next(csv_stream)]
        name_index = headers.index("name")
        value_index = headers.index("value")
        return name_index, value_index

    def __SetMeasuredEigenvector(self, eigen_measurement_data_file_name: str, prefix: str) -> None:
        hdf5_read_params = Kratos.Parameters("""{
            "file_name": "../damaged_system/measurement_eigenvector_data.h5",
            "file_access_mode": "read_only"
        }""")
        #hdf5_read_params["file_name"] = eigen_measurement_data_file_name
        
        with OpenHDF5File(hdf5_read_params, self.evaluated_model_part) as h5_file:
            
            #expio = KratosHDF5.HDF5NodalSolutionStepDataIO(prefix_settings, h5_file)
            #expio.ReadNodalResults(self.model.GetModelPart("Structure"), False)
            
            KratosHDF5.ReadDataValueContainer(h5_file, prefix, self.evaluated_model_part.ProcessInfo)
            self.evaluated_model_part.SetValue(KratosSI.MEASURED_EIGENVALUE_VECTOR, self.evaluated_model_part.ProcessInfo[KratosSI.MEASURED_EIGENVALUE_VECTOR])
            nodal_io_settings = Kratos.Parameters("""
                {
                    "list_of_variables": ["MEASURED_EIGENVECTOR_MATRIX"],
                    "prefix" : ""
                }
                """)
            nodal_io_settings["prefix"].SetString(prefix)
            nodal_data_value_io = KratosHDF5.HDF5NodalDataValueIO(nodal_io_settings, h5_file)
            nodal_data_value_io.ReadNodalResults(self.evaluated_model_part.Nodes, self.evaluated_model_part.GetCommunicator())
        

    def __GetSensitivtyVariables(self, sensitivity_setting) -> 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]':
        sensitivity_settings = sensitivity_setting
        return {
            ExpressionDataLocation.Element: [Kratos.KratosGlobals.GetVariable("YOUNG_MODULUS_SENSITIVITY")]
        }
    
    def __GetSensitivities(self, sensitivity_setting) -> 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]':
        sensitivity_model_part: Kratos.ModelPart = self.GetSensitivityModelPart()
        sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self.__GetSensitivtyVariables(sensitivity_setting)

        result: 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]' = {}
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                result[variable] = GetContainerExpression(sensitivity_model_part, data_location, variable)
        return result

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"