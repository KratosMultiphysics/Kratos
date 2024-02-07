import matplotlib as mpl
mpl.use('Agg') # to have plotting possible without a Xserver running.
import matplotlib.pyplot as plt
from pathlib import Path

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetComponentHavingDataByFullName

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Kratos.OutputProcess:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemGraphOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemGraphOutputProcess(parameters["settings"], optimization_problem)

class GraphData:
    """A class to hold information about each graph.
    """
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        default_settings = Kratos.Parameters("""{
            "y_axis_label"   : "",
            "y_min"          : "0.0",
            "y_max"          : "<initial_max_value> * 1.5",
            "is_log"         : false,
            "legend_position": "upper right",
            "components"     : [""]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.y_axis_label = parameters["y_axis_label"].GetString()
        self.y_min = parameters["y_min"].GetString()
        self.y_max = parameters["y_max"].GetString()
        self.is_log = parameters["is_log"].GetBool()
        self.legend_position = parameters["legend_position"].GetString()
        self.component_paths: 'list[str]' = parameters["components"].GetStringArray()

        self.optimization_problem = optimization_problem
        self.data: 'list[list[float]]' = [[] for _ in self.component_paths]

        self.graph_names: 'list[str]' = []
        for component_path in self.component_paths:
            data = component_path.split(".")
            if data[0] in ["response_function", "control", "execution_policy"]:
                if len(data) < 3:
                    raise RuntimeError(f"\"{component_path}\" does not have the component value specified.")
                self.graph_names.append(".".join(data[1:]))
            else:
                if len(data) < 2:
                    raise RuntimeError(f"\"{component_path}\" does not have the component value specified.")
                self.graph_names.append(component_path)

    def GetYAxisLabel(self) -> str:
        return self.y_axis_label

    def GetLegendPosition(self) -> str:
        return self.legend_position

    def GetGraphNames(self) -> 'list[str]':
        return self.graph_names

    def IsLog(self) -> bool:
        return self.is_log

    def GetYLimits(self) -> 'tuple[float, float]':
        initial_max_value = -1e+100
        initial_min_value = 1e+100
        for data in self.data:
            initial_min_value = min(initial_min_value, data[0])
            initial_max_value = max(initial_max_value, data[0])

        y_min = Kratos.GenericFunctionUtility(self.y_min.replace("<initial_min_value>", str(initial_min_value))).CallFunction(0.0,0.0,0.0,0,0.0,0.0,0.0)
        y_max = Kratos.GenericFunctionUtility(self.y_max.replace("<initial_max_value>", str(initial_max_value))).CallFunction(0.0,0.0,0.0,0,0.0,0.0,0.0)
        return (y_min, y_max)

    def GetData(self) -> 'list[list[float]]':
        return self.data

    def Update(self) -> None:
        """Updates the internal data storage with the current values of the components.

        Raises:
            RuntimeError: If the component is not found.
        """
        for i, component_path in enumerate(self.component_paths):
            data = component_path.split(".")
            if data[0] in ["response_function", "control", "execution_policy"]:
                component = GetComponentHavingDataByFullName(f"{data[0]}.{data[1]}", self.optimization_problem)
                value_path = "/".join(data[2:])
            else:
                component = GetComponentHavingDataByFullName(data[0], self.optimization_problem)
                value_path = "/".join(data[1:])

            component_data_view = ComponentDataView(component, self.optimization_problem)
            component_buffered_data = component_data_view.GetBufferedData()
            self.data[i].append(component_buffered_data[value_path])

class OptimizationProblemGraphOutputProcess(Kratos.OutputProcess):
    """An output process which can be used to plot buffered data in the Optimization problem

    This process can be used to plot values from the buffered data in each component of the optimization problem.
    """
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "output_file_name": "SPECIFY_OUTPUT_FILE_NAME",
                "dpi"             : 300,
                "output_interval" : 1,
                "interval"        : [0, "End"],
                "max_iterations"  : 1e+3,
                "x_axis_label"    : "Iteration [-]",
                "graphs"          : []
            }
            """
        )

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        self.optimization_problem = optimization_problem
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.output_file_name = parameters["output_file_name"].GetString()
        self.x_axis_label = parameters["x_axis_label"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        output_file_name_data = Path(self.output_file_name).name.split(".")
        if len(output_file_name_data) == 1:
            self.output_file_name += ".png"

        self.max_iterations = parameters["max_iterations"].GetInt()
        self.output_interval = parameters["output_interval"].GetInt()
        self.dpi = parameters["dpi"].GetInt()

        self.list_of_graph_data: 'list[GraphData]' = []
        for graph_params in parameters["graphs"].values():
            self.list_of_graph_data.append(GraphData(graph_params, self.optimization_problem))

        self.iterations: 'list[int]' = []
        self.interval_utility = Kratos.IntervalUtility(parameters)

    def ExecuteFinalizeSolutionStep(self) -> None:
        list(map(lambda x: x.Update(), self.list_of_graph_data))
        self.iterations.append(self.optimization_problem.GetStep())

    def IsOutputStep(self) -> bool:
        current_step = self.optimization_problem.GetStep()
        if self.interval_utility.IsInInterval(current_step):
            begin_time = self.interval_utility.GetIntervalBegin()
            if ((current_step - begin_time) % self.output_interval) == 0:
                return True
        return False

    def PrintOutput(self) -> None:
        if self._IsWritingProcess():
            fig, axes = plt.subplots(len(self.list_of_graph_data), sharex=True)
            for i, graph_data in enumerate(self.list_of_graph_data):
                if len(self.list_of_graph_data) == 1:
                    ax = axes
                else:
                    ax = axes[i]
                for i, y_data in enumerate(graph_data.data):
                    if graph_data.IsLog():
                        ax.semilogy(self.iterations, y_data, label=graph_data.GetGraphNames()[i])
                    else:
                        ax.plot(self.iterations, y_data, label=graph_data.GetGraphNames()[i])

                    ax.set_xlim([min(self.iterations), self.max_iterations])
                    ax.set_ylim(graph_data.GetYLimits())
                    ax.legend(loc=graph_data.GetLegendPosition())
                    ax.set_ylabel(graph_data.GetYAxisLabel())
                    ax.grid(True)

            ax.set_xlabel(self.x_axis_label)
            output_path = Path(self.output_file_name.replace("<step>", str(self.optimization_problem.GetStep())))
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(str(output_path), dpi=self.dpi, bbox_inches="tight")
            plt.close()

    def _IsWritingProcess(self):
        if Kratos.IsDistributedRun():
            data_communicator: Kratos.DataCommunicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
            return data_communicator.Rank() == 0
        else:
            return True
