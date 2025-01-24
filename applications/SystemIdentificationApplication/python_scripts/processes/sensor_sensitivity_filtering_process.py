import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorSensitivityFilteringProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorSensitivityFilteringProcess(model, parameters["settings"], optimization_problem)


class SensorSensitivityFilteringProcess(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""{
            "sensor_group_name"        : "",
            "filtering_model_part_name": "",
            "filtering_variable_name"  : "",
            "data_location"            : "",
            "filter_settings"          : {},
            "expression_names"         : [
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

        filtering_model_part_name = parameters["filtering_model_part_name"].GetString()
        filtering_variable = Kratos.KratosGlobals.GetVariable(parameters["filtering_variable_name"].GetString())

        data_location_name = parameters["data_location"].GetString()
        if data_location_name == "node_historical":
            data_location = Kratos.Globals.DataLocation.NodeHistorical
        elif data_location_name == "node_non_historical":
            data_location = Kratos.Globals.DataLocation.NodeNonHistorical
        elif data_location_name == "condition":
            data_location = Kratos.Globals.DataLocation.Condition
        elif data_location_name == "element":
            data_location = Kratos.Globals.DataLocation.Element
        else:
            raise RuntimeError(f"Unsupported data location = \"{data_location_name}\".")

        self.exp_filter = FilterFactory(model, filtering_model_part_name, filtering_variable, data_location, parameters["filter_settings"])

    def ExecuteBeforeSolutionLoop(self):
        self.exp_filter.SetComponentDataView(ComponentDataView(self.sensor_group_name, self.optimization_problem))
        self.exp_filter.Initialize()

    def Check(self):
        self.exp_filter.Check()

    def ExecuteFinalizeSolutionStep(self) -> None:
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)

        for input_name, output_name in self.expressions_names_list:
            expression: ContainerExpressionTypes = sensor_group_data.GetUnBufferedData().GetValue(input_name)
            abs_expression = self.exp_filter.ForwardFilterField(self.exp_filter.BackwardFilterIntegratedField(expression))
            sensor_group_data.GetUnBufferedData().SetValue(output_name, abs_expression.Clone(), overwrite=True)

    def ExecuteFinalize(self):
        self.exp_filter.Finalize()

