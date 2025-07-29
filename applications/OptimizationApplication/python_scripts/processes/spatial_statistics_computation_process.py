import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import SetComponentValueByFullName
import KratosMultiphysics.StatisticsApplication.spatial_utilities as SpatialUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    if not parameters.Has("settings"):
        raise RuntimeError(f"SpatialStatisticsComputationProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SpatialStatisticsComputationProcess(model, parameters["settings"], optimization_problem)

class SpatialStatisticsComputationProcess(Kratos.OutputProcess):
    """
    A process for computing statistical operations on variables within a Kratos model part
    and storing the results in an optimization problem's component data.

    Statistic Settings (each entry in "statistic_settings"):
        variable_name (str): Name of the variable to compute statistics on.
        method_name (str): Statistical method to apply (e.g., mean, max, min).
        norm_type (str): Type of norm to use (default: "magnitude").
        container (str): Data container type (e.g., "nodal_non_historical").
        output_value_name (str): Template for naming the output value.

    Refer documentation of the StatisticsApplication for further information
    on which methods available and how they can be used.

    StatisticsApplication documentation:
    https://kratosmultiphysics.github.io/Kratos/pages/Applications/Statistics_Application/General/Overview.html

    Usage:
        - Configure the process with appropriate settings for the desired statistics.
        - At each output step, call PrintOutput() to compute and store statistics.
    """

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name"      : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "output_component_name": "algorithm",
            "statistic_settings"   : []
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.__model_part = model[settings["model_part_name"].GetString()]
        self.__optimization_problem = optimization_problem
        self.__component_name = settings["output_component_name"].GetString()

        stat_default_settings = Kratos.Parameters("""{
            "variable_name"    : "",
            "method_name"      : "",
            "norm_type"        : "magnitude",
            "container"        : "nodal_non_historical",
            "output_value_name": "<METHOD_HEADER>_<VARIABLE_NAME>:historical"
        }""")

        self.__list_of_stat_operations = []
        for stat_settings in settings["statistic_settings"].values():
            stat_settings.ValidateAndAssignDefaults(stat_default_settings)
            stat_method_name = stat_settings["method_name"].GetString()
            norm_type = stat_settings["norm_type"].GetString()
            container = stat_settings["container"].GetString()
            variable_name = stat_settings["variable_name"].GetString()
            output_value_name = stat_settings["output_value_name"].GetString().replace("<VARIABLE_NAME>", variable_name)
            operation_type = getattr(SpatialUtilities, f"Spatial{Kratos.StringUtilities.ConvertSnakeCaseToCamelCase(stat_method_name)}Output")
            operation = operation_type(container, norm_type, variable_name)
            self.__list_of_stat_operations.append((operation, output_value_name))

    def IsOutputStep(self):
        return True

    def PrintOutput(self) -> None:
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        data_view = ComponentDataView(component, self.__optimization_problem)

        for operation, value_name in self.__list_of_stat_operations:
            operation.Evaluate(self.__model_part)
            data: 'list[float]' = operation.data
            for value, method_name in zip(data, operation.GetHeaders()):
                SetComponentValueByFullName(data_view, value_name.replace("<METHOD_HEADER>", method_name), value)
