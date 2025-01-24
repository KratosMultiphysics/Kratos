import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"AbsoluteSensitivityComputationProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return AbsoluteSensitivityComputationProcess(parameters["settings"], optimization_problem)


class AbsoluteSensitivityComputationProcess(Kratos.Process):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name": "",
            "expression_names" : [
                {
                    "input" : "",
                    "output": ""
                }
            ]
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.sensor_group_name = parameters["sensor_group_name"].GetString()

        list_of_expression_parameters: 'list[Kratos.Parameters]' = parameters["expression_names"].values()
        self.expressions_names_list: 'list[tuple[str, str]]' = []
        for expression_parameter in list_of_expression_parameters:
            expression_parameter.ValidateAndAssignDefaults(default_parameters["expression_names"].values()[0])
            self.expressions_names_list.append((expression_parameter["input"].GetString(), expression_parameter["output"].GetString()))

        self.optimization_problem = optimization_problem

    def ExecuteFinalizeSolutionStep(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for input_name, output_name in self.expressions_names_list:
            expression: ContainerExpressionTypes = sensor_group_data.GetUnBufferedData().GetValue(input_name)
            abs_expression = Kratos.Expression.Utils.Abs(expression)
            sensor_group_data.GetUnBufferedData().SetValue(output_name, abs_expression.Clone(), overwrite=True)
