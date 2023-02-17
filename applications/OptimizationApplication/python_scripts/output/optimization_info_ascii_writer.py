from datetime import datetime
from abc import ABC
from abc import abstractmethod
from typing import Union

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetOptimizationInfoAvailableKeysForType
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPatternMatchedOptimizationKeys
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPrefixWithoutDataName
from KratosMultiphysics.OptimizationApplication.output.output_utilities import GetPlaceHolderWithValues

class OptimizationAsciiInfo(ABC):
    @classmethod
    @abstractmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        pass

    @classmethod
    @abstractmethod
    def GetPattern(cls) -> str:
        pass

    def __init__(self, optimization_info: OptimizationInfo, settings: Kratos.Parameters):
        default_parameters = self.GetDefaultParameters()

        settings.ValidateAndAssignDefaults(default_parameters)
        self.header_prefix = settings["header_prefix"].GetString()

        self.optimization_info = optimization_info

        # now create the common prefix
        self.common_prefix = GetPrefixWithoutDataName(self.GetPattern(), settings)
        self.place_holder_values = GetPlaceHolderWithValues(self.GetPattern(), settings)

        if self.header_prefix == "":
            self.header_prefix = self.model_part_name + "_" + self.data_type_name

        self.data_list: 'list[dict[str, any]]' = []
        for data in settings["data_list"]:
            data.ValidateAndAssignDefaults(default_parameters["data_list"][0])
            self.data_list.append({
                "location" : f'{self.common_prefix}/{data["data_name"].GetString()}',
                "data_name": data["data_name"].GetString(),
                "format"   : data["format"].GetString(),
                "header"   : data["header_suffix"].GetString()
            })

    def Check(self):
        available_keys = GetOptimizationInfoAvailableKeysForType(self.optimization_info, Union[int, float, str, bool])

        # check for model part and response name availability
        if not self.optimization_info.HasValue(self.common_prefix):
            available_options = GetPatternMatchedOptimizationKeys(available_keys, self.GetPattern(), False)
            raise RuntimeError(f"The provided {self.place_holder_values} not found. Followings are available pairs:\n\t" + "\n\t".join(available_options))

        # now check for available data names
        for data in self.data_list:
            data_name = data["data_name"]
            location = data["location"]
            data_format = data["format"]

            if not self.optimization_info.HasValue(location):
                available_options = GetPatternMatchedOptimizationKeys(available_keys, self.common_prefix, True)
                raise RuntimeError(f"The provided \"data_name\" = \"{data_name}\" with {self.place_holder_values} not found. Followings are available data names:\n\t" + "\n\t".join(available_options))
            else:
                value = self.optimization_info.GetValue(location)
                value_format, value_format_exp = OptimizationAsciiInfo.GetFormat(value)
                if not data_format.endswith(value_format_exp):
                    raise RuntimeError(f"The provided \"data_name\" = \"{data_name}\" with {self.place_holder_values} is of type {value_format} which is not the requested type.")

    def GetFormattedData(self) -> 'list[str]':
        self.Check()
        return [OptimizationAsciiInfo.GetFormattedValue(self.optimization_info.GetValue(data["location"]), data["format"]) for data in self.data_list]

    def GetHeaders(self) -> 'list[str]':
        header_names = []
        for data in self.data_list:
            header = data["header"]
            if header == "":
                header_names.append(self.header_prefix + "_" + data["data_name"])
            else:
                header_names.append(self.header_prefix + "_" + header)

        return header_names

    def GetHeaderLengths(self) -> 'list[int]':
        lengths: 'list[int]' = []
        for header, data in zip(self.GetHeaders(), self.data_list):
            data_format: str = data["format"]
            data_type = data_format[-1]
            if data_type == "b":
                value = True
            elif data_type == "d":
                value = -1
            elif data_type == "e":
                value = -1.0
            elif data_type == "s":
                value = "      "
            else:
                raise RuntimeError(f"The header length cannot be determined for header {header}. [ Data format = {data_format}]")
            lengths.append(max(len(OptimizationAsciiInfo.GetFormattedValue(value, data_format)), len(header)))
        return lengths

    @staticmethod
    def GetFormattedValue(value: any, value_format: str) -> str:
        if isinstance(value, bool):
            if value:
                return "yes"
            else:
                return " no"
        else:
            return ("{:" + value_format + "}").format(value)

    @staticmethod
    def GetFormat(value: any) -> str:
        if isinstance(value, str):
            return "string", "s"
        elif isinstance(value, bool):
            return "bool", "b"
        elif isinstance(value, int):
            return "int", "d"
        elif isinstance(value, float):
            return "float", "e"
        else:
            raise RuntimeError(f"The format cannot be determined. [ Value = {str(value)}]")

class ResponseAsciiInfo(OptimizationAsciiInfo):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"           : "response_data",
            "model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME",
            "response_name"  : "PLEASE_PROVIDE_RESPONSE_FUNCTION_NAME",
            "header_prefix"  : "",
            "data_list"      : [
                {
                    "data_name"    : "PLEASE_PROVIDE_DATA_NAME",
                    "format"       : "PLEASE_SPECIFY_FORMAT",
                    "header_suffix": ""
                }
            ]
        }""")

    @classmethod
    def GetPattern(cls) -> str:
        return "problem_data/response_data/<model_part_name>/<response_name>"

class ControlAsciiInfo(OptimizationAsciiInfo):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"           : "control_data",
            "model_part_name": "PLEASE_PROVIDE_MODEL_PART_NAME",
            "control_name"   : "PLEASE_PROVIDE_CONTROL_NAME",
            "header_prefix"  : "",
            "data_list"      : [
                {
                    "data_name"    : "PLEASE_PROVIDE_DATA_NAME",
                    "format"       : "PLEASE_SPECIFY_FORMAT",
                    "header_suffix": ""
                }
            ]
        }""")

    @classmethod
    def GetPattern(cls) -> str:
        return "problem_data/control_data/<model_part_name>/<control_name>"

class AlgorithmAsciiInfo(OptimizationAsciiInfo):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"           : "algrithm_data",
            "header_prefix"  : "",
            "data_list"      : [
                {
                    "data_name"    : "PLEASE_PROVIDE_DATA_NAME",
                    "format"       : "PLEASE_SPECIFY_FORMAT",
                    "header_suffix": ""
                }
            ]
        }""")

    @classmethod
    def GetPattern(cls) -> str:
        return "problem_data/algorithm_data/info"

class OptimizationInfoAsciiWriter(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optmization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "module"                 : "KratosMultiphysics.OptimizationApplication.utilities",
            "type"                   : "OptimizationInfoAsciiWriter",
            "model_part_name"        : "",
            "interval"               : [0.0, "End"],
            "output_step_interval"   : 1,
            "write_kratos_version"   : false,
            "write_time_stamp"       : false,
            "output_value_precision" : 5,
            "output_value_length"    : 8,
            "output_file_settings"   : {
                "file_name"  : "<model_part_name>.dat",
                "output_path": "Optimization_Results",
                "write_buffer_size" : -1
            },
            "optimization_info_output_settings": []
        }""")

        self.optimization_info = optmization_info
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.write_kratos_version = parameters["write_kratos_version"].GetBool()
        self.write_time_stamp = parameters["write_time_stamp"].GetBool()
        self.output_file_settings = parameters["output_file_settings"]

        self.model_part = model[parameters["model_part_name"].GetString()]

        self.list_of_data_writers: 'list[OptimizationAsciiInfo]' = []
        for optimization_info_output_settings in parameters["optimization_info_output_settings"]:
            if optimization_info_output_settings.Has("type"):
                writer_type = optimization_info_output_settings["type"].GetString()
            else:
                raise RuntimeError(f"\"type\" is not found in \"optimization_info_output_settings\".")

            if writer_type == "response_data":
                self.list_of_data_writers.append(ResponseAsciiInfo(self.optimization_info, optimization_info_output_settings))
            elif writer_type == "control_data":
                self.list_of_data_writers.append(ControlAsciiInfo(self.optimization_info, optimization_info_output_settings))
            elif writer_type == "algorithm_data":
                self.list_of_data_writers.append(AlgorithmAsciiInfo(self.optimization_info, optimization_info_output_settings))
            else:
                raise RuntimeError(f"Unsupported \"{writer_type}\" found. Followings are available options:\n\tresponse_data\n\tcontrol_data\n\talgorithm_data")

        self.interval_utility = Kratos.IntervalUtility(parameters)
        self.output_interval = parameters["output_step_interval"].GetDouble()
        self.output_value_precision = parameters["output_value_precision"].GetInt()
        self.output_value_length = parameters["output_value_length"].GetInt()

        if self._IsWritingProcess():
            if self.output_file_settings.Has("file_name"):
                self.output_file_settings["file_name"].SetString(self.output_file_settings["file_name"].GetString().replace("<model_part_name>", self.model_part.FullName()))
            self._WriteHeaders()

    def IsOutputStep(self) -> bool:
        return True

    def PrintOutput(self):
        self._WriteData()

    def ExecuteFinalize(self):
        self.output_file.file.write("# End of file.")

    def _IsWritingProcess(self):
        return (self.model_part.GetCommunicator().MyPID() == 0)

    def _WriteHeaders(self):
        if (self._IsWritingProcess()):
            kratos_version = "not_given"
            if (self.write_kratos_version):
                kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

            time_stamp = "not_specified"
            if (self.write_time_stamp):
                time_stamp = str(datetime.now())

            msg_header = ""
            msg_header += "# Optimization info ascii output\n"
            msg_header += "# Kratos version               : " + kratos_version + "\n"
            msg_header += "# Timestamp                    : " + time_stamp + "\n"
            msg_header += "# Modelpart name               : " + self.model_part.FullName() + "\n"
            msg_header += "# ----------------------------------------------------------------------\n"
            msg_header += "# Headers:\n"

            msg_header += "# STEP"
            for writer in self.list_of_data_writers:
                for header_name, header_length in zip(writer.GetHeaders(), writer.GetHeaderLengths()):
                    msg_header += (", {:>" + str(header_length) + "s}").format(header_name)

            msg_header += "\n"

            self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, self.output_file_settings, msg_header)

            self.is_headers_written = True

    def _WriteData(self):
        if self._IsWritingProcess():
            msg_data = ""

            # writing the step
            msg_data += "{:>6d}".format(self.optimization_info["step"])
            for writer in self.list_of_data_writers:
                for data, header_length in zip(writer.GetFormattedData(), writer.GetHeaderLengths()):
                    msg_data += (", {:>" + str(header_length) + "s}").format(data)

            self.output_file.file.write(msg_data + "\n")


