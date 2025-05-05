import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.SystemIdentificationApplication.responses.damage_detection_response import DamageDetectionResponse
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"MaskedDamageDetectionResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaskedDamageDetectionResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MaskedDamageDetectionResponse(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class MaskedDamageDetectionResponse(DamageDetectionResponse):
    @staticmethod
    def GetDefaultParameters() -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "sensor_group_name"          : "",
            "sensor_mask_name"           : "",
            "output_all_sensitivities"   : true,
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

    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name, model, parameters, optimization_problem)
        self.sensor_mask_name = parameters["sensor_mask_name"].GetString()
        self.output_all_sensitivities = parameters["output_all_sensitivities"].GetBool()

    def Initialize(self):
        super().Initialize()

        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        for controller in GetMaskStatusControllers(sensor_group_data, self.sensor_mask_name):
            if isinstance(controller, KratosSI.SensorMaskStatus):
                self.mask_status = controller

    def CalculateValue(self) -> float:
        self.damage_response_function.Clear()
        for sensor in self.list_of_sensors:
            self.damage_response_function.AddSensor(sensor)
        return super().CalculateValue()

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # make everything zeros
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for container_expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.LiteralExpressionIO.SetDataToZero(container_expression, physical_variable)

        # now compute sensitivities for each test scenario
        for exec_policy, sensor_measurement_data_file_name, test_case_weight in self.list_of_test_analysis_data:
            # read and replace the measurement data for each test scenario
            self._SetSensorMeasuredValue(sensor_measurement_data_file_name)

            # run a single adjoint for each test scenario
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = exec_policy.GetName()
            self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()

            # now compute the gradients of the each sensor error
            for sensor_index, sensor in enumerate(self.list_of_sensors):
                self.damage_response_function.Clear()
                self.damage_response_function.AddSensor(sensor)
                mask = sensor.GetContainerExpression(self.sensor_mask_name).Clone()
                self.adjoint_analysis._GetSolver().GetComputingModelPart().ProcessInfo[KratosSI.SENSOR_NAME] = f"{self.GetName()}({sensor.GetName()})"
                sensitivities = self.adjoint_analysis.CalculateGradient(self.damage_response_function)

                for physical_variable, collective_expression in physical_variable_collective_expressions.items():
                    for container_expression in collective_expression.GetContainerExpressions():
                        sensitivity_variable = Kratos.KratosGlobals.GetVariable(f"{physical_variable.Name()}_SENSITIVITY")

                        physical_var_sensitivities = sensitivities[sensitivity_variable] * mask

                        if self.output_all_sensitivities:
                            response_data = ComponentDataView(self, self.optimization_problem)
                            response_data.GetUnBufferedData().SetValue(f"d{self.GetName()}({sensor.GetName()})_d{physical_variable.Name()}", sensitivities[sensitivity_variable].Clone(), overwrite=True)
                            response_data.GetUnBufferedData().SetValue(f"dmasked_{self.GetName()}({sensor.GetName()})_d{physical_variable.Name()}", physical_var_sensitivities.Clone(), overwrite=True)

                        container_expression.SetExpression((container_expression.GetExpression() - physical_var_sensitivities.GetExpression() * test_case_weight))
                        container_expression.SetExpression(Kratos.Expression.Utils.Collapse(container_expression).GetExpression())

