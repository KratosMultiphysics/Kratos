from datetime import datetime
from typing import Union, Any

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetAllComponentFullNamesWithData
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetComponentHavingDataByFullName

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("settings"):
        raise RuntimeError(f"OptimizationProblemAsciiOutputProcess instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return OptimizationProblemAsciiOutputProcess(parameters["settings"], optimization_problem)

class Header:
    def __init__(self, header_name: str, value: Any, format_info: dict):
        header_name = header_name.strip()
        header_length = len(header_name)

        if isinstance(value, bool):
            value_length = max(len(format_info[type(value)][0]), len(format_info[type(value)][1]))
            value_format_post_fix = "s"
            self.__value_converter = lambda x: format_info[type(value)][1] if x else format_info[type(value)][0]
        elif isinstance(value, int):
            value_length = len(("{:" + str(format_info[type(value)]) + "d}").format(value))
            value_format_post_fix = "d"
            self.__value_converter = lambda x: int(x)
        elif isinstance(value, float):
            value_length = len(("{:0." + str(format_info[type(value)]) + "e}").format(value))
            value_format_post_fix = f".{format_info[type(value)]}e"
            self.__value_converter = lambda x: float(x)
        else:
            value_length = format_info[str]
            value_format_post_fix = "s"
            self.__value_converter = lambda x: str(x)

        if header_length > value_length:
            self.__header_name = header_name
            self.__value_format = "{:>" + str(header_length) + value_format_post_fix + "}"
        else:
            self.__header_name = ("{:>" + str(value_length) + "s}").format(header_name)
            self.__value_format = "{:>" + str(value_length) + value_format_post_fix + "}"

    def GetHeaderName(self) -> str:
        return self.__header_name

    def GetValueStr(self, value: Any) -> str:
        return self.__value_format.format(self.__value_converter(value))

class OptimizationProblemAsciiOutputProcess(Kratos.OutputProcess):
    def GetDefaultParameters(self):
        return Kratos.Parameters(
            """
            {
                "output_file_name"         : "SPECIFY_OUTPUT_FILE_NAME",
                "write_kratos_version"     : true,
                "write_time_stamp"         : true,
                "write_initial_values"     : true,
                "list_of_output_components": ["all"],
                "format_info": {
                    "int_length"     : 7,
                    "float_precision": 9,
                    "bool_values"    : ["no", "yes"],
                    "string_length"  : 10
                }
            }
            """
        )

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Kratos.OutputProcess.__init__(self)

        self.optimization_problem = optimization_problem
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.output_file_name = parameters["output_file_name"].GetString()
        if not self.output_file_name.endswith(".csv"):
            self.output_file_name += ".csv"

        self.write_kratos_version = parameters["write_kratos_version"].GetBool()
        self.write_time_stamp = parameters["write_time_stamp"].GetBool()
        self.write_initial_values = parameters["write_initial_values"].GetBool()

        self.format_info = {
            int  : parameters["format_info"]["int_length"].GetInt(),
            float: parameters["format_info"]["float_precision"].GetInt(),
            bool : parameters["format_info"]["bool_values"].GetStringArray(),
            str  : parameters["format_info"]["string_length"].GetInt(),
        }

        if len(self.format_info[bool]) != 2:
            raise RuntimeError("The \"bool_values\" should have only two strings corresponding to False and True values in the mentioned order.")

        self.list_of_components: 'list[Union[str, ResponseFunction, Control, ExecutionPolicy]]' = []
        self.list_of_headers: 'list[tuple[Any, dict[str, Header]]]' = []
        self.initialized_headers = False
        self.list_of_component_names = parameters["list_of_output_components"].GetStringArray()

    def IsOutputStep(self) -> bool:
        return True

    def PrintOutput(self) -> None:
        if not self.initialized_headers:
            if len(self.list_of_component_names) == 1 and self.list_of_component_names[0] == "all":
                self.list_of_component_names = GetAllComponentFullNamesWithData(self.optimization_problem)

            for component_name in self.list_of_component_names:
                self.list_of_components.append(GetComponentHavingDataByFullName(component_name, self.optimization_problem))

            # now get the buffered data headers
            self.list_of_headers = self._GetHeaders(lambda x: x.GetBufferedData() if x.HasDataBuffer() else BufferedDict())
            # write the ehader information
            self._WriteHeaders()
            self.initialized_headers = True

        if self._IsWritingProcess():
            # now write step data
            with open(self.output_file_name, "a") as file_output:
                # write the step
                file_output.write("{:>7d}".format(self.optimization_problem.GetStep()))

                # wrtie the values
                for component, header_info_dict in self.list_of_headers:
                    componend_data_view = ComponentDataView(component, self.optimization_problem)
                    buffered_dict = componend_data_view.GetBufferedData()
                    for k, header in header_info_dict.items():
                        file_output.write(", " + header.GetValueStr(buffered_dict[k]))

                file_output.write("\n")

    def ExecuteFinalize(self):
        if self._IsWritingProcess():
            with open(self.output_file_name, "a") as file_output:
                file_output.write("# End of file")

    def _IsWritingProcess(self):
        if Kratos.IsDistributedRun():
            data_communicator: Kratos.DataCommunicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
            return data_communicator.Rank() == 0
        else:
            return True

    def _WriteHeaders(self):
        if (self._IsWritingProcess()):
            kratos_version = "not_given"
            if (self.write_kratos_version):
                kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

            time_stamp = "not_specified"
            if (self.write_time_stamp):
                time_stamp = str(datetime.now())

            msg_header = ""
            msg_header = f"{msg_header}# Optimization problem ascii output\n"
            msg_header = f"{msg_header}# Kratos version: {kratos_version}\n"
            msg_header = f"{msg_header}# Timestamp     : {time_stamp}\n"
            msg_header = f"{msg_header}# -----------------------------------------------\n"

            if self.write_initial_values:
                msg_header = f"{msg_header}# --------------- Initial values ----------------\n"

                initial_headers = self._GetHeaders(lambda x: x.GetUnBufferedData())
                # now write the initial value container data
                for component, header_info_dict in initial_headers:
                    componend_data_view = ComponentDataView(component, self.optimization_problem)
                    buffered_dict = componend_data_view.GetUnBufferedData()
                    component_name = componend_data_view.GetComponentName()
                    # check if there are values to be written under the component name, if not skip the component.
                    if len(header_info_dict):
                        msg_header = f"{msg_header}# \t" + component_name + ":\n"
                        for k, header in header_info_dict.items():
                            component_name_header = header.GetHeaderName().strip()[len(component_name)+1:]
                            msg_header = f"{msg_header}# \t\t" + component_name_header + ": " + header.GetValueStr(buffered_dict[k]).strip() + "\n"

                msg_header = f"{msg_header}# ------------ End of initial values ------------\n"
                msg_header = f"{msg_header}# -----------------------------------------------\n"

            msg_header = f"{msg_header}# ------------ Start of step values -------------\n"
            msg_header = f"{msg_header}# Headers:\n"

            msg_header = f"{msg_header}#  STEP"

            for _, header_info_dict in self.list_of_headers:
                for header in header_info_dict.values():
                    msg_header = f"{msg_header}, " + header.GetHeaderName()

            msg_header = f"{msg_header}\n"

            # write the header
            with open(self.output_file_name, "w") as file_output:
                file_output.write(msg_header)

    def _GetHeaders(self, dict_getter_method) ->  'list[tuple[Any, dict[str, Header]]]':
        list_of_headers: 'list[tuple[Any, dict[str, Header]]]' = []
        for component in self.list_of_components:
            componend_data_view = ComponentDataView(component, self.optimization_problem)
            values_map = dict_getter_method(componend_data_view).GetMap()
            header_info_dict: 'dict[str, Header]' = {}
            component_name = componend_data_view.GetComponentName()
            for k, v in values_map.items():
                if isinstance(v, (bool, int, float, str)):
                    header_name = component_name + ":" + k[k.rfind("/") + 1:]
                    if header_name in [header.GetHeaderName().strip() for header in header_info_dict.values()]:
                        Kratos.Logger.PrintWarning(self.__class__.__name__, "Second value with same header name = \"" + header_name + "\" found.")
                    header_info_dict[k] = Header(header_name, v, self.format_info)
            if len(header_info_dict):
                list_of_headers.append([component, header_info_dict])
        return list_of_headers