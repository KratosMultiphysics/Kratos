import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorMaskComputationProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorMaskComputationProcess(parameters["settings"], optimization_problem)


class SensorMaskComputationProcess(Kratos.Process):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"    : "",
            "automatic_threshold"  : true,
            "sensitivity_threshold": -1,
            "expression_names"     : [
                {
                    "input" : "",
                    "output": ""
                }
            ]
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.sensor_group_name = parameters["sensor_group_name"].GetString()

        self.optimization_problem = optimization_problem

        list_of_expression_parameters: 'list[Kratos.Parameters]' = parameters["expression_names"].values()
        self.expressions_names_list: 'list[tuple[str, str]]' = []
        for expression_parameter in list_of_expression_parameters:
            expression_parameter.ValidateAndAssignDefaults(default_parameters["expression_names"].values()[0])
            self.expressions_names_list.append((expression_parameter["input"].GetString(), expression_parameter["output"].GetString()))

        self.automatic_threshold = parameters["automatic_threshold"].GetBool()
        self.sensitivity_threshold = parameters["sensitivity_threshold"].GetDouble()

    def ExecuteFinalizeSolutionStep(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for input_name, output_name in self.expressions_names_list:
            expression: ContainerExpressionTypes = sensor_group_data.GetUnBufferedData().GetValue(input_name)
            if self.automatic_threshold:
                mask_exp = KratosSI.MaskUtils.GetMask(expression)
            else:
                mask_exp = KratosSI.MaskUtils.GetMask(expression, self.sensitivity_threshold)
            sensor_group_data.GetUnBufferedData().SetValue(output_name, mask_exp.Clone(), overwrite=True)
