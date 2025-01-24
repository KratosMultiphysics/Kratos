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
            "expression_names" : []
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.sensor_group_name = parameters["sensor_group_name"].GetString()
        self.expression_names = parameters["expression_names"].GetStringArray()

        self.optimization_problem = optimization_problem

    def ExecuteFinalizeSolutionStep(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for expression_name in self.expression_names:
            expression: ContainerExpressionTypes = sensor_group_data.GetUnBufferedData().GetValue(expression_name)
            abs_expression = Kratos.Expression.Utils.Abs(expression)
            sensor_group_data.GetUnBufferedData().SetValue(f"{expression_name}_abs", abs_expression.Clone(), overwrite=True)
