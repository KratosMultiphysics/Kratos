from typing import Any

import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStat
from KratosMultiphysics.process_factory import KratosProcessFactory

from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.StatisticsApplication.spatial_statistics_operation import GetLengthAdjustedValue
from KratosMultiphysics.StatisticsApplication.spatial_statistics_operation import GetFormattedInt
from KratosMultiphysics.StatisticsApplication.spatial_statistics_operation import GetFormattedFloat
from KratosMultiphysics.StatisticsApplication.spatial_statistics_operation import SpatialStatisticsOperation
from KratosMultiphysics.StatisticsApplication.spatial_statistics_operation import GetSpatialStatisticsOperation

from datetime import datetime

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SpatialStatisticsProcess(model, settings["Parameters"])

def GetNorm(norm_name: str) -> Any:
    if norm_name == "none":
        return None
    elif norm_name == "l2":
        return KratosStat.Norms.L2()
    elif norm_name == "infinity":
        return KratosStat.Norms.Infinity()
    elif norm_name == "trace":
        return KratosStat.Norms.Trace()
    elif norm_name.startswith("pnorm_"):
        return KratosStat.Norms.P(float(norm_name[6:]))
    elif norm_name.startswith("lpqnorm_("):
        pos = norm_name.rfind(",")
        if pos != -1:
            p = float(norm_name[9:pos])
            q = float(norm_name[pos+1:-1])
        else:
            raise RuntimeError(f"Unsupported norm info provided for LPQ norm [ provided info = \"{norm_name}\", required info = \"lpqnorm_(p,q)\" ].")
        return KratosStat.Norms.LPQ(p, q)
    else:
        raise RuntimeError(f"Unsupported norm name requested [ requested norm name = \"{norm_name}\" ]. Followings are supported:\n\tnone\n\tl2\n\tinfinity\n\ttrace\n\tpnorm_p\n\tlpqnorm_(p,q)")

def GetContainerLocation(container_location_name: str) -> Kratos.Globals.DataLocation:
    locations_dict = {
        "nodal_historical"        : Kratos.Globals.DataLocation.NodeHistorical,
        "nodal_non_historical"    : Kratos.Globals.DataLocation.NodeNonHistorical,
        "condition_non_historical": Kratos.Globals.DataLocation.Condition,
        "element_non_historical"  : Kratos.Globals.DataLocation.Element
    }

    if container_location_name in locations_dict.keys():
        return locations_dict[container_location_name]
    else:
        raise RuntimeError(f"Unsupported container location name [ provided container location name = \"{container_location_name}\" ]. Followings are supported:\n\t" + "\n\t".join(locations_dict.keys()))


class SpatialStatisticsProcess(Kratos.OutputProcess):
    """A process to calculate spatial statistics on Kratos containers

    This process calculates spatial statistics for given variables in a given container.

    This process is compatible with OpenMP and MPI with restart

    Args:
        model (Kratos.Model): Model used in problem
        settings (Kratos.Parameters): Kratos parameter settings for process
    """
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        Kratos.OutputProcess.__init__(self)

        default_parameters = Kratos.Parameters("""
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"              : 0,
            "computation_processes"   : [],
            "computation_points"      : [
                "ExecuteInitialize",
                "ExecuteInitializeSolutionStep",
                "Execute",
                "ExecuteFinalizeSolutionStep"
            ],
            "input_variable_settings" : [
                {
                    "variable_names" : [],
                    "norm_type"      : "none",
                    "container"      : "nodal_historical"
                }
            ],
            "statistics_methods": [
                {
                    "method_name"    : "",
                    "method_settings": {}
                }
            ],
            "output_settings" : {
                "interval"               : [0.0, "End"],
                "output_control_variable": "STEP",
                "output_time_interval"   : 1,
                "write_kratos_version"   : false,
                "write_time_stamp"       : false,
                "output_value_precision" : 5,
                "output_file_settings"   : {
                    "file_name"  : "<model_part_name>.dat",
                    "output_path": "spatial_statistics_output",
                    "write_buffer_size" : -1
                }
            }
        }  """)

        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.variables_settings_list = self.settings["input_variable_settings"]

        self.echo_level = self.settings["echo_level"].GetInt()

        self.output_settings = self.settings["output_settings"]
        self.output_settings["output_file_settings"].ValidateAndAssignDefaults(default_parameters["output_settings"]["output_file_settings"])
        self.output_settings.RecursivelyValidateAndAssignDefaults(default_parameters["output_settings"])

        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        ## set the output controller settings
        self.interval_utility = Kratos.IntervalUtility(self.output_settings)
        output_control_variable_name = self.output_settings["output_control_variable"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(output_control_variable_name)):
            raise RuntimeError("Unknown output control variable. [ \"output_control_variable\" = \"" + output_control_variable_name + "\" ]")
        output_control_variable_type = Kratos.KratosGlobals.GetVariableType(output_control_variable_name)
        if (output_control_variable_type not in ["Integer", "Double"]):
            raise RuntimeError("Unsupported output control variable type for " + output_control_variable_name + " of " + output_control_variable_type + " type. Supported types are Integer and Double only")
        self.output_control_variable = Kratos.KratosGlobals.GetVariable(output_control_variable_name)
        self.output_interval = self.output_settings["output_time_interval"].GetDouble()
        self.output_value_precision = self.output_settings["output_value_precision"].GetInt()

        # get the list of variables
        list_of_variable_info = []
        for variable_settings in self.settings["input_variable_settings"].values():
            norm = GetNorm(variable_settings["norm_type"].GetString())
            data_location = GetContainerLocation(variable_settings["container"].GetString())
            for variable_name in variable_settings["variable_names"].GetStringArray():
                variable = Kratos.KratosGlobals.GetVariable(variable_name)
                list_of_variable_info.append([variable, data_location, norm])

        # set the input variavle x methods operations list
        self.list_of_operations: 'list[SpatialStatisticsOperation]' = []
        for method_settings in self.settings["statistics_methods"].values():
            method_settings.ValidateAndAssignDefaults(default_parameters["statistics_methods"].values()[0])
            for variable_info in list_of_variable_info:

                self.list_of_operations.append(GetSpatialStatisticsOperation(method_settings["method_name"].GetString(), self.model_part, variable_info[0], variable_info[1], variable_info[2], self.output_value_precision, method_settings["method_settings"].Clone()))

        ## set the execution process
        self.computation_points = self.settings["computation_points"].GetStringArray()
        if self.settings["computation_processes"].size() == 0 or len(self.computation_points) == 0:
            self.statistics_computation_processeses = None
        else:
            self.statistics_computation_processeses = KratosProcessFactory(self.model).ConstructListOfProcesses(self.settings["computation_processes"])
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Created statistics computation processes with following parameters: \n{:s}".format(self.settings["computation_processes"].PrettyPrintJsonString()))

        self.is_headers_written = False

        Kratos.Logger.PrintInfo("SpatialStatisticsProcess", "Initialized statistics process.")

    def Check(self):
        for statistics_operation in self.list_of_operations:
            statistics_operation.Check()

        if self.statistics_computation_processeses is not None:
            for process in self.statistics_computation_processeses:
                process.Check()

    def ExecuteInitialize(self):
        ## set output file path
        self.output_settings["output_file_settings"]["file_name"].SetString(self.output_settings["output_file_settings"]["file_name"].GetString().replace("<model_part_name>", self.model_part.FullName()))
        self.output_settings["output_file_settings"]["output_path"].SetString(self.output_settings["output_file_settings"]["output_path"].GetString().replace("<model_part_name>", self.model_part.FullName()))

        self.__write_headers()

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, "List of output operations:\n\t" + "\n\t".join([str(operation) for operation in self.list_of_operations]))

        self.previous_process_info_value = self.model_part.ProcessInfo[self.output_control_variable]

    def IsOutputStep(self) -> bool:
        n = self.model_part.ProcessInfo[self.output_control_variable]
        n_begin = self.interval_utility.GetIntervalBegin()
        n_intervals = int((n - n_begin) / self.output_interval)
        return abs(n_intervals*self.output_interval+n_begin - n) < 1e-10 and self.interval_utility.IsInInterval(n)

    def PrintOutput(self) -> None:
        # execute the statistics computation process
        if self.statistics_computation_processeses is not None:
            for process in self.statistics_computation_processeses:
                for computation_point in self.computation_points:
                    getattr(process, computation_point)()

        self.__write_values()

        if self.echo_level > 1:
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Calculated statistics on {:s}.".format(self.model_part.FullName()))

    def ExecuteFinalize(self):
        if (self.__is_writing_process()):
            self.output_file.file.write("# End of file.")
            self.output_file.file.close()

    def __write_headers(self):
        kratos_version = "not_given"
        if (self.output_settings["write_kratos_version"].GetBool()):
            kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

        time_stamp = "not_specified"
        if (self.output_settings["write_time_stamp"].GetBool()):
            time_stamp = str(datetime.now())

        operation_info_dict:'dict[str, dict[str, int]]' = {
            "variable": {},
            "method"  : {},
            "location": {},
            "norm"    : {}
        }
        for statistics_operation in self.list_of_operations:
            info = statistics_operation.GetVarianbleInfo()
            if info not in operation_info_dict["variable"].keys():
                operation_info_dict["variable"][info] = 0
            operation_info_dict["variable"][info] += 1

            info = statistics_operation.GetNormInfo()
            if info not in operation_info_dict["norm"].keys():
                operation_info_dict["norm"][info] = 0
            operation_info_dict["norm"][info] += 1

            info = statistics_operation.GetContainerInfo()
            if info not in operation_info_dict["location"].keys():
                operation_info_dict["location"][info] = 0
            operation_info_dict["location"][info] += 1

            info = statistics_operation.GetMethodInfo()
            if info not in operation_info_dict["method"].keys():
                operation_info_dict["method"][info] = 0
            operation_info_dict["method"][info] += 1

        header_details: 'dict[str, bool]' = {}
        for info_type, operation_info in operation_info_dict.items():
            header_details[info_type] = True
            for v in operation_info.values():
                if v == len(self.list_of_operations):
                    header_details[info_type] = False

        # now add the headers from all the methods in all ranks
        all_rank_headers = ""
        for statistics_operation in self.list_of_operations:
            all_rank_headers += f", {statistics_operation.GetHeadersString(header_details)}"

        if (self.__is_writing_process()):

            msg_header = ""
            msg_header += "# Spatial statistics process output\n"
            msg_header += "# Kratos version               : " + kratos_version + "\n"
            msg_header += "# Timestamp                    : " + time_stamp + "\n"
            msg_header += "# Modelpart name               : " + self.model_part.FullName() + "\n"
            msg_header += "# Output control variable name : " + self.output_control_variable.Name() + "\n"

            if not header_details["method"]:
                msg_header += "# Methon                       : " + list(operation_info_dict["method"].keys())[0] + "\n"
            if not header_details["variable"]:
                msg_header += "# Variable                     : " + list(operation_info_dict["variable"].keys())[0] + "\n"
            if not header_details["location"]:
                msg_header += "# Location                     : " + list(operation_info_dict["location"].keys())[0] + "\n"
            if not header_details["norm"]:
                msg_header += "# Norm                         : " + list(operation_info_dict["norm"].keys())[0] + "\n"

            msg_header += "# ----------------------------------------------------------------------\n"
            msg_header += "# Headers:\n"

            # add output control header
            if isinstance(self.output_control_variable, Kratos.IntegerVariable):
                v_str = GetFormattedInt(int(1+99), self.output_value_precision)
            elif isinstance(self.output_control_variable, Kratos.DoubleVariable):
                v_str = GetFormattedFloat(1+99, self.output_value_precision)
            if len(v_str) > len(self.output_control_variable.Name()):
                msg_header += ("#{:>" + str(len(v_str)) + "s}").format(self.output_control_variable.Name())
                self.output_control_value = lambda : GetLengthAdjustedValue(self.model_part.ProcessInfo[self.output_control_variable], len(v_str)+1, self.output_value_precision)
            else:
                msg_header += f"#{self.output_control_variable.Name()}"
                self.output_control_value = lambda : GetLengthAdjustedValue(self.model_part.ProcessInfo[self.output_control_variable], len(self.output_control_variable.Name()) + 1, self.output_value_precision)

            msg_header += all_rank_headers

            msg_header += "\n"

            self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, self.output_settings["output_file_settings"], msg_header)

            self.is_headers_written = True

    def __write_values(self):
        # now write the method values in all ranks
        all_rank_values = ""
        for statistics_operation in self.list_of_operations:
            all_rank_values += f", {statistics_operation.GetValueString()}"

        if (self.__is_writing_process()):
            # first write the output control value
            self.output_file.file.write(self.output_control_value())

            self.output_file.file.write(all_rank_values)

            self.output_file.file.write("\n")

    def __is_writing_process(self):
        return (self.model_part.GetCommunicator().MyPID() == 0)
