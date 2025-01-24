import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorSensitivityAggregationProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorSensitivityAggregationProcess(model, parameters["settings"], optimization_problem)


class SensorSensitivityAggregationProcess(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"     : "",
            "expression_names"      : [],
            "output_expression_name": "<TEST_ANALYSIS_NAME>_<EXPRESSION_NAME>"
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.optimization_problem = optimization_problem
        self.model = model

        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.expression_names = parameters["expression_names"].GetStringArray()
        self.output_expression_name = parameters["output_expression_name"].GetString()

    def ExecuteFinalizeSolutionStep(self):
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        sensor_model_part = self.model[self.sensor_group_name]

        list_of_sensors = GetSensors(sensor_group_data)
        current_sensor_name: str = sensor_model_part.ProcessInfo[KratosSI.SENSOR_NAME]
        test_analysis_name: str = sensor_model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME]

        found_sensor = False
        for sensor in list_of_sensors:
            if current_sensor_name == sensor.GetName():
                found_sensor = True
                for expression_name in self.expression_names:
                    expression: ContainerExpressionTypes = sensor_group_data.GetUnBufferedData().GetValue(expression_name)
                    current_output_expression_name = self.output_expression_name.replace("<TEST_ANALYSIS_NAME>", test_analysis_name)
                    current_output_expression_name = current_output_expression_name.replace("<EXPRESSION_NAME>", expression_name)
                    sensor.AddContainerExpression(current_output_expression_name, expression)
                break

        if not found_sensor:
            raise RuntimeError(f"The sensor \"{current_sensor_name}\" not found.")


