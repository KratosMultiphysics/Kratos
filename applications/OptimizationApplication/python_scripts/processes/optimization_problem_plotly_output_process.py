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


# class Header:
#     def __init__(self, header_name: str, value: Any, format_info: dict):
#         header_name = header_name.strip()
#         header_length = len(header_name)

#         if isinstance(value, bool):
#             value_length = max(len(format_info[type(value)][0]), len(format_info[type(value)][1]))
#             value_format_post_fix = "s"
#             self.__value_converter = lambda x: format_info[type(value)][1] if x else format_info[type(value)][0]
#         elif isinstance(value, int):
#             value_length = len(("{:" + str(format_info[type(value)]) + "d}").format(value))
#             value_format_post_fix = "d"
#             self.__value_converter = lambda x: int(x)
#         elif isinstance(value, float):
#             value_length = len(("{:0." + str(format_info[type(value)]) + "e}").format(value))
#             value_format_post_fix = f".{format_info[type(value)]}e"
#             self.__value_converter = lambda x: float(x)
#         else:
#             value_length = format_info[str]
#             value_format_post_fix = "s"
#             self.__value_converter = lambda x: str(x)

#         if header_length > value_length:
#             self.__header_name = header_name
#             self.__value_format = "{:>" + str(header_length) + value_format_post_fix + "}"
#         else:
#             self.__header_name = ("{:>" + str(value_length) + "s}").format(header_name)
#             self.__value_format = "{:>" + str(value_length) + value_format_post_fix + "}"

#     def GetHeaderName(self) -> str:
#         return self.__header_name

#     def GetValueStr(self, value: Any) -> str:
#         return self.__value_format.format(self.__value_converter(value))


class OptimizationProblemPlotlyOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "list_of_output_components": ["all"]
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
            fig = px.line(self.output_dataframe)
            fig.show()

    def _IsWritingProcess(self):
        if Kratos.IsDistributedRun():
            data_communicator: Kratos.DataCommunicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
            return data_communicator.Rank() == 0
        else:
            return True
