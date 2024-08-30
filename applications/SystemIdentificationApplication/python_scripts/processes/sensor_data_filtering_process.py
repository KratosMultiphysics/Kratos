import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.processes.execution_point_process import ExecutionPointProcess

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.Process:
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorDataFilteringProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorDataFilteringProcess(model, parameters["settings"], optimization_problem)

class SensorDataFilteringProcess(ExecutionPointProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "execution_point"              : "",
            "filtering_model_part_name"    : "",
            "filtering_expression_variable": "",
            "filtering_data_location"      : "",
            "filtering_settings"           : {},
            "filtering_expression_settings": [
                {
                    "input_name" : "",
                    "output_name": ""
                }
            ]
        }""")

        self.model = model
        self.optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(default_settings)

        super().__init__(parameters["execution_point"].GetString())

        self.filtering_expression_data: 'list[tuple[str, str]]' = []
        for filtering_expression_params in parameters["filtering_expression_settings"].values():
            filtering_expression_params.ValidateAndAssignDefaults(default_settings["filtering_expression_settings"].values()[0])
            self.filtering_expression_data.append((filtering_expression_params["input_name"].GetString(), filtering_expression_params["output_name"].GetString()))

        data_location_name = parameters["filtering_data_location"].GetString()
        if data_location_name == "node_historical":
            data_location = Kratos.Globals.DataLocation.NodeHistorical
        elif data_location_name == "node_non_historical":
            data_location = Kratos.Globals.DataLocation.NodeNonHistorical
        elif data_location_name == "condition":
            data_location = Kratos.Globals.DataLocation.Condition
        elif data_location_name == "element":
            data_location = Kratos.Globals.DataLocation.Element
        else:
            raise RuntimeError(f"Unsupported data location name = \"{data_location_name}\" provided. Followings are supported:\n\tnode_historical\n\tnode_non_historical\n\tcondition\n\telement")

        self.filtering_model_part_name = parameters["filtering_model_part_name"].GetString()
        self.variable = Kratos.KratosGlobals.GetVariable(parameters["filtering_expression_variable"].GetString())
        self.filter = FilterFactory(self.model, self.filtering_model_part_name, self.variable, data_location, parameters["filtering_settings"])

        self.filter_initialized = False

    def Execute(self) -> None:
        if not self.filter_initialized:
            self.filter_initialized = True
            self.filter.SetComponentDataView(ComponentDataView(f"{self.filtering_model_part_name}_{self.variable.Name()}", self.optimization_problem))
            self.filter.Initialize()

        list_of_sensors = GetSensors(self.optimization_problem)

        for sensor in list_of_sensors:
            for inp_exp_name, out_exp_name in self.filtering_expression_data:
                input_exp = sensor.GetContainerExpression(inp_exp_name)
                output_exp = self.filter.ForwardFilterField(self.filter.BackwardFilterIntegratedField(input_exp))
                sensor.AddContainerExpression(out_exp_name, output_exp)
