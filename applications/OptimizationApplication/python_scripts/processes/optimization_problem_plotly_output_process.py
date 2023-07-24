from datetime import datetime
from typing import Union, Any

import pandas as pd
import plotly.express as px

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetComponentHavingDataByFullName


def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemPlotlyOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemPlotlyOutputProcess(parameters["settings"], optimization_problem)


class OptimizationProblemPlotlyOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "list_of_output_components": ["all"],
                "display_plot": true,
                "write_html_output": false,
                "html_file_name": "iteration_summary",
                "plot_y_log_scale": false
            }
            """
        )

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        self.optimization_problem = optimization_problem
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        list_of_component_names = parameters["list_of_output_components"].GetStringArray()
        if len(list_of_component_names) == 1 and list_of_component_names[0] == "all":
            list_of_component_names = GetAllComponentFullNamesWithData(optimization_problem)

        for component_name in list_of_component_names:
            self.list_of_components.append(GetComponentHavingDataByFullName(component_name, optimization_problem))

        self.output_parameters = dict()
        self.output_parameters["write_html_output"] = parameters["write_html_output"].GetBool()
        self.output_parameters["display_plot"] = parameters["display_plot"].GetBool()
        self.output_parameters["html_file_name"] = parameters["html_file_name"].GetString()
        self.output_parameters["plot_y_log_scale"] = parameters["plot_y_log_scale"].GetBool()

        if self.output_parameters["html_file_name"].find(".html") == -1:
            self.output_parameters["html_file_name"] += ".html"

        self.output_dict = {}
        self.output_dataframe = pd.DataFrame()

    def IsOutputStep(self) -> bool:
        return True

    def PrintOutput(self) -> None:

        if self._IsWritingProcess():

            for component in self.list_of_components:
                componend_data_view = ComponentDataView(component, self.optimization_problem)
                buffered_dict = componend_data_view.GetBufferedData()
                values_map = buffered_dict.GetMap()
                for k, v in values_map.items():
                    if isinstance(v, (bool, int, float, str)):
                        if str(k) not in self.output_dict:
                            self.output_dict[str(k)] = []
                        self.output_dict[str(k)].append(v)

            self.output_dataframe = pd.DataFrame.from_dict(self.output_dict)

    def ExecuteFinalize(self):
        if self._IsWritingProcess():
            fig = px.line(self.output_dataframe, log_y=self.output_parameters["plot_y_log_scale"])

            if self.output_parameters["display_plot"]:
                fig.show()

            if self.output_parameters["write_html_output"]:
                fig.write_html(self.output_parameters["html_file_name"])

    def _IsWritingProcess(self):
        if Kratos.IsDistributedRun():
            data_communicator: Kratos.DataCommunicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
            return data_communicator.Rank() == 0
        else:
            return True
