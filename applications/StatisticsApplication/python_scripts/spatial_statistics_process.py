import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

from KratosMultiphysics.StatisticsApplication.method_utilities import GetAvailableMethods
from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetFormattedValue
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

import KratosMultiphysics.StatisticsApplication.spatial_utilities as SpatialUtilities

from datetime import datetime


def Factory(settings, model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SpatialStatisticsProcess(model, settings["Parameters"])


class SpatialStatisticsProcess(Kratos.Process):
    """A process to calculate spatial statistics on Kratos containers

    This process calculates spatial statistics for given variables in a given container.

    This process is compatible with OpenMP and MPI with restart

    Args:
        model (Kratos.Model): Model used in problem
        settings (Kratos.Parameters): Kratos parameter settings for process
    """
    def __init__(self, model, settings):
        Kratos.Process.__init__(self)

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
                "output_value_length"    : 14,
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

        self.model_part_name = self.settings["model_part_name"].GetString()

        ## set the input variavle x methods operations list
        self.list_of_operations = []

        spatial_output_method_info = []
        for spatial_method_output in dir(SpatialUtilities):
            if (spatial_method_output.startswith("Spatial") and spatial_method_output.endswith("Output")):
                spatial_output_method_info.append([spatial_method_output.lower()[7:-6], spatial_method_output])

        spatial_method_output_names_ = [method_info[0] for method_info in spatial_output_method_info]
        spatial_method_output = [method_info[1] for method_info in spatial_output_method_info]

        for variable_settings in self.settings["input_variable_settings"].values():
            variable_settings.ValidateAndAssignDefaults(default_parameters["input_variable_settings"][0])

            container_type = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()

            container = GetItemContainer(container_type)
            norm_container = GetNormTypeContainer(container, norm_type)
            available_methods, _ = GetAvailableMethods(norm_container)

            for variable_name in variable_settings["variable_names"].GetStringArray():
                if (not Kratos.KratosGlobals.HasVariable(variable_name)):
                    raise RuntimeError("Variable not found. [ variable_name = " + variable_name + " ]")

                for method_settings in self.settings["statistics_methods"].values():
                    method_settings.ValidateAndAssignDefaults(default_parameters["statistics_methods"][0])

                    method_name = method_settings["method_name"].GetString()
                    if method_name == "":
                        raise RuntimeError("Found an empty method name. {:s}".format(method_settings))

                    if method_name not in spatial_method_output_names_:
                        msg = "Unknown method name [ \"method_name\" = \"" + method_name + "\" ]\n"
                        msg += "Allowed method names are:\n    "
                        msg += "\n    ".join(sorted(spatial_method_output_names_))
                        raise RuntimeError(msg)


                    if not method_name in available_methods:
                        Kratos.Logger.PrintInfo(self.__class__.__name__, "Skipping statistics output for {:s} method for {:s} in {:s} because {:s} norm type does not support it.".format(method_name, variable_name, container_type, norm_type))
                    else:
                        operation_type = getattr(SpatialUtilities, spatial_method_output[spatial_method_output_names_.index(method_name)])
                        self.list_of_operations.append(operation_type(container_type, norm_type, variable_name, method_settings["method_settings"]))

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
        self.output_value_length = self.output_settings["output_value_length"].GetInt()

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
        if (not self.model.HasModelPart(self.model_part_name)):
            raise RuntimeError(self.model_part_name + " not found.")

        for variable_settings in self.variables_settings_list.values():
            variable_names_list = variable_settings["variable_names"].GetStringArray()

            container = variable_settings["container"].GetString()
            if (container == "nodal_historical"):
                for variable_name in variable_names_list:
                    variable = Kratos.KratosGlobals.GetVariable(variable_name)
                    if (not self.__get_model_part().HasNodalSolutionStepVariable(variable)):
                        raise RuntimeError("Variable " + variable_name + " not found in nodal solution step data of " + self.__get_model_part().Name)

        if self.statistics_computation_processeses is not None:
            for process in self.statistics_computation_processeses:
                process.Check()

    def ExecuteInitialize(self):
        ## set output file path
        self.output_settings["output_file_settings"]["file_name"].SetString(self.output_settings["output_file_settings"]["file_name"].GetString().replace("<model_part_name>", self.__get_model_part().FullName()))
        self.output_settings["output_file_settings"]["output_path"].SetString(self.output_settings["output_file_settings"]["output_path"].GetString().replace("<model_part_name>", self.__get_model_part().FullName()))

        self.__write_headers()

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, "List of output operations:\n\t" + "\n\t".join([str(operation) for operation in self.list_of_operations]))

        self.previous_process_info_value = self.__get_model_part().ProcessInfo[self.output_control_variable]

    def ExecuteFinalizeSolutionStep(self):
        output_control_counter = self.__get_model_part().ProcessInfo[self.output_control_variable]

        if self.interval_utility.IsInInterval(output_control_counter):
            if (output_control_counter - self.previous_process_info_value) >= self.output_interval:

                # execute the statistics computation process
                if self.statistics_computation_processeses is not None:
                    for process in self.statistics_computation_processeses:
                        for computation_point in self.computation_points:
                            getattr(process, computation_point)()

                for operation in self.list_of_operations:
                    operation.Evaluate(self.model_part)

                self.__write_values()

                if self.echo_level > 1:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, "Calculated statistics on {:s}.".format(self.__get_model_part().FullName()))

                self.previous_process_info_value = output_control_counter

    def ExecuteFinalize(self):
        if (self.__is_writing_process()):
            self.output_file.file.write("# End of file.")
            self.output_file.file.close()

    def __write_headers(self):
        if (self.__is_writing_process()):
            kratos_version = "not_given"
            if (self.output_settings["write_kratos_version"].GetBool()):
                kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

            time_stamp = "not_specified"
            if (self.output_settings["write_time_stamp"].GetBool()):
                time_stamp = str(datetime.now())

            method_name = ""
            variable_info = ""
            for operation in self.list_of_operations:
                if operation.method_name not in method_name:
                    method_name += operation.method_name + ", "
                current_variable_info = "{:s}[{:s}, {:s}], ".format(operation.variable_name, operation.norm_type, operation.container_type)
                if current_variable_info not in variable_info:
                    variable_info += current_variable_info

            self.headers, self.header_lengths = self.__get_headers()

            msg_header = ""
            msg_header += "# Spatial statistics process output\n"
            msg_header += "# Kratos version               : " + kratos_version + "\n"
            msg_header += "# Timestamp                    : " + time_stamp + "\n"
            msg_header += "# Method name(s)               : " + method_name[:-2] + "\n"
            msg_header += "# Variable info(s)             : " + variable_info[:-2] + "\n"
            msg_header += "# Modelpart name               : " + self.model_part_name + "\n"
            msg_header += "# Output control variable name : " + self.output_control_variable.Name() + "\n"
            msg_header += "# ----------------------------------------------------------------------\n"
            msg_header += "# Headers:\n"
            msg_header += "#" + " ".join([("{:>" + str(header_length) + "s}").format(header) for header, header_length in zip(self.headers, self.header_lengths)]) + "\n"

            self.output_file = TimeBasedAsciiFileWriterUtility(self.__get_model_part(), self.output_settings["output_file_settings"], msg_header)

            self.is_headers_written = True

    def __get_headers(self):
        # check whether only one variable with one
        variable_info = ""
        for operation in self.list_of_operations:
            current_variable_info = "{:s}[{:s}, {:s}], ".format(operation.variable_name, operation.norm_type, operation.container_type)
            if current_variable_info not in variable_info:
                variable_info += current_variable_info

        headers = ["\"TIME\""]
        header_lengths = [max(6, self.output_value_length)]
        if self.output_control_variable != Kratos.TIME:
            headers.append("\"" + self.output_control_variable.Name() + "\"")
            header_lengths.append(max(len(self.output_control_variable.Name()) + 2, self.output_value_length))

        # found only one variable type info
        if variable_info.count(",") == 2:
            header_formatting_string = "\"<HEADER_NAME><COMPONENT>\""
        else:
            header_formatting_string = "\"<VARIABLE_NAME><COMPONENT> [ <NORM_TYPE> | <CONTAINER_TYPE> ] <HEADER_NAME>\""

        variable_components = {
            "Double": [""],
            "Array" : ["_X", "_Y", "_Z"]
        }

        for operation in self.list_of_operations:
            if operation.norm_type == "none":
                comps_list = variable_components[Kratos.KratosGlobals.GetVariableType(operation.variable.Name())]
            else:
                comps_list = variable_components["Double"]

            for current_header, current_header_length in zip(operation.GetHeaders(), operation.GetValueLengths(self.output_value_length)):
                for comp in comps_list:
                    current_header_name = str(header_formatting_string)

                    current_header_name = current_header_name.replace("<HEADER_NAME>", current_header)
                    current_header_name = current_header_name.replace("<VARIABLE_NAME>", operation.variable.Name())
                    current_header_name = current_header_name.replace("<COMPONENT>", comp)
                    current_header_name = current_header_name.replace("<NORM_TYPE>", operation.norm_type)
                    current_header_name = current_header_name.replace("<CONTAINER_TYPE>", operation.container_type)

                    headers.append(current_header_name)
                    header_lengths.append(max(len(current_header_name), current_header_length))

        return headers, header_lengths

    def __write_values(self):
        if (self.__is_writing_process()):
            values = [GetFormattedValue(self.__get_model_part().ProcessInfo[Kratos.TIME], self.output_value_length, self.output_value_precision)[0]]
            if self.output_control_variable != Kratos.TIME:
                values.append(GetFormattedValue(self.__get_model_part().ProcessInfo[self.output_control_variable], self.output_value_length, self.output_value_precision)[0])
            for operation in self.list_of_operations:
                values.extend(operation.GetValues(self.output_value_length, self.output_value_precision))
            self.output_file.file.write(" " + " ".join([("{:>" + str(header_length) + "s}").format(value) for value, header_length in zip(values, self.header_lengths)]) + "\n")

    def __get_model_part(self):
        if (not hasattr(self, "model_part")):
            self.model_part = self.model[self.model_part_name]

        return self.model_part

    def __is_writing_process(self):
        return (self.__get_model_part().GetCommunicator().MyPID() == 0)
