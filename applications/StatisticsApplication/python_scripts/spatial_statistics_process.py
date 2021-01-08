import KratosMultiphysics as Kratos

from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod

from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetMethodHeaders
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetMethodValues
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetVariableHeaders
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
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
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_settings" : [
                {
                    "method_name"    : "sum",
                    "norm_type"      : "none",
                    "container"      : "nodal_historical",
                    "variable_names" : [],
                    "method_settings": {}
                }
            ],
            "output_settings" : {
                "output_control_variable": "STEP",
                "output_time_interval"   : 1,
                "write_kratos_version"   : true,
                "write_time_stamp"       : true,
                "output_file_settings"   : {
                    "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
                    "folder_name": "spatial_statistics_output",
                    "write_buffer_size" : -1
                }
            }
        }  """)

        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.variables_settings_list = self.settings["input_variable_settings"]
        self.output_settings = self.settings["output_settings"]
        self.output_settings.RecursivelyValidateAndAssignDefaults(default_parameters["output_settings"])

        for variable_settings in self.variables_settings_list:
            variable_settings.ValidateAndAssignDefaults(default_parameters["input_variable_settings"][0])
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()

            item_container = GetItemContainer(container_name)
            item_norm_container = GetNormTypeContainer(item_container, norm_type)
            GetMethod(item_norm_container, method_name)

        self.model_part_name = self.settings["model_part_name"].GetString()

        output_control_variable_name = self.settings["output_settings"]["output_control_variable"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(output_control_variable_name)):
            raise Exception("Unknown output control variable. [ \"output_control_variable\" = \"" + output_control_variable_name + "\" ]")

        Kratos.Logger.PrintInfo("SpatialStatisticsProcess", "Initialized statistics process.")

    def Check(self):
        if (not self.model.HasModelPart(self.model_part_name)):
            raise Exception(self.model_part_name + " not found.")

        for variable_settings in self.variables_settings_list:
            variable_names_list = variable_settings["variable_names"].GetStringArray()
            for variable_name in variable_names_list:
                if (not Kratos.KratosGlobals.HasVariable(variable_name)):
                    raise Exception("Variable not found. [ variable_name = " + variable_name + " ]")

            container = variable_settings["container"].GetString()
            if (container == "nodal_historical"):
                for variable_name in variable_names_list:
                    variable = Kratos.KratosGlobals.GetVariable(variable_name)
                    if (not self.__get_model_part().HasNodalSolutionStepVariable(variable)):
                        raise Exception("Variable " + variable_name + " not found in nodal solution step data of " + self.__get_model_part().Name)

        process_info = self.__get_model_part().ProcessInfo
        output_settings = self.settings["output_settings"]
        output_control_variable_name = output_settings["output_control_variable"].GetString()
        output_control_variable = Kratos.KratosGlobals.GetVariable(output_control_variable_name)
        output_control_variable_type = Kratos.KratosGlobals.GetVariableType(output_control_variable_name)

        if (not Kratos.KratosGlobals.HasVariable(output_control_variable_name)):
            raise Exception("Unknown output control variable " + output_control_variable_name)

        if (output_control_variable_type not in ["Integer", "Double"]):
            raise Exception("Unsupported output control variable type for " + output_control_variable_name + " of " + output_control_variable_type + " type. Supported types are Integer and Double only")

        process_info_value = process_info[output_control_variable]
        self.output_control_counter = process_info_value
        self.previous_process_info_value = process_info_value

    def ExecuteInitialize(self):
        kratos_version = "not_given"
        if (self.output_settings["write_kratos_version"].GetBool()):
            kratos_version = str(Kratos.KratosGlobals.Kernel.Version())

        time_stamp = "not_specified"
        if (self.output_settings["write_time_stamp"].GetBool()):
            time_stamp = str(datetime.now())

        output_control_variable_name = self.output_settings["output_control_variable"].GetString()

        self.output_files = []
        for variable_settings in self.variables_settings_list:
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()

            msg_header = ""
            msg_header += "# Spatial statistics process output\n"
            msg_header += "# Kratos version               : " + kratos_version + "\n"
            msg_header += "# Timestamp                    : " + time_stamp + "\n"
            msg_header += "# Method Name                  : " + method_name + "\n"
            msg_header += "# Norm type                    : " + norm_type + "\n"
            msg_header += "# Container type               : " + container_name + "\n"
            msg_header += "# Modelpart name               : " + self.model_part_name + "\n"
            msg_header += "# Output control variable name : " + output_control_variable_name + "\n"
            msg_header += "# ----------------------------------------------------------------------\n"
            msg_header += "# Headers:\n"

            output_file_settings = self.output_settings["output_file_settings"]
            output_file_name_syntax = output_file_settings["file_name"].GetString()
            output_file_name = output_file_name_syntax.replace("<model_part_name>", self.model_part_name)
            output_file_name = output_file_name.replace("<container>", container_name)
            output_file_name = output_file_name.replace("<norm_type>", norm_type)
            output_file_name = output_file_name.replace("<method_name>", method_name)

            current_output_file_settings = Kratos.Parameters("""{}""")
            current_output_file_settings.AddEmptyValue("file_name")
            current_output_file_settings["file_name"].SetString(output_file_name)
            current_output_file_settings.AddEmptyValue("folder_name")
            current_output_file_settings["folder_name"].SetString(output_file_settings["folder_name"].GetString())
            # restarting is not supported if STEP is used as the control variable
            if (self.__is_writing_process()):
                self.output_files.append(TimeBasedAsciiFileWriterUtility(self.__get_model_part(), current_output_file_settings, msg_header))
            else:
                self.output_files.append("dummy")

    def ExecuteFinalizeSolutionStep(self):
        output_settings = self.settings["output_settings"]
        output_control_variable_name = output_settings["output_control_variable"].GetString()
        output_control_variable = Kratos.KratosGlobals.GetVariable(output_control_variable_name)
        output_control_variable_type = Kratos.KratosGlobals.GetVariableType(output_control_variable_name)

        if (output_control_variable_type == "Integer"):
            current_output_control_variable_value = output_settings["output_time_interval"].GetInt()
        elif (output_control_variable_type == "Double"):
            current_output_control_variable_value = output_settings["output_time_interval"].GetDouble()

        self.process_info_value = self.__get_model_part().ProcessInfo[output_control_variable]

        self.output_control_counter += (self.process_info_value - self.previous_process_info_value)

        if (self.output_control_counter >= current_output_control_variable_value):
            self.CalculateOutput()
            self.output_control_counter = 0
            if (self.__is_writing_process()):
                for output_file in self.output_files:
                    output_file.file.write("\n")
                    output_file.file.flush()

        self.previous_process_info_value = self.process_info_value

    def ExecuteFinalize(self):
        if (self.__is_writing_process()):
            for output_file in self.output_files:
                output_file.file.write("# End of file\n")
                output_file.file.close()

    def CalculateOutput(self):
        for variable_settings, output_file in zip(self.variables_settings_list, self.output_files):
            self.is_output_control_variable_value_written = False
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()

            item_container = GetItemContainer(container_name)
            item_norm_container = GetNormTypeContainer(item_container, norm_type)
            method = GetMethod(item_norm_container, method_name)

            variable_list = []
            variable_names_list = variable_settings["variable_names"].GetStringArray()
            for variable_name in variable_names_list:
                variable_list.append(Kratos.KratosGlobals.GetVariable(variable_name))

            if (norm_type == "none"):
                for index, variable in enumerate(variable_list):
                    output = method(self.__get_model_part(), variable)
                    method_headers = GetMethodHeaders(method_name, variable_settings["method_settings"])
                    self.__write_output(output,  variable_names_list[index], variable_names_list, norm_type, method_name, method_headers, output_file)
            else:
                for index, variable in enumerate(variable_list):
                    output = method(self.__get_model_part(), variable, norm_type, variable_settings["method_settings"])
                    method_headers = GetMethodHeaders(method_name, variable_settings["method_settings"])
                    self.__write_output(output, variable_names_list[index], variable_names_list, norm_type, method_name, method_headers, output_file)


    def __write_output(self, output, variable_name, variable_names_list, norm_type, method_name, method_headers, output_file):
        if (self.__is_writing_process()):
            self.__write_headers(norm_type, variable_names_list, method_headers, output_file)
            if (not self.is_output_control_variable_value_written):
                output_file.file.write(str(self.process_info_value))
                self.is_output_control_variable_value_written = True
            output_file.file.write(" ")
            output_file.file.write(GetMethodValues(method_name, norm_type, variable_name, output))


    def __write_headers(self, norm_type, variable_names_list, method_headers, output_file):
        if (self.__is_writing_process()):
            if (not hasattr(output_file, "is_variable_headers_written")):
                output_file.is_variable_headers_written = False

            if (not output_file.is_variable_headers_written):
                msg = "# OutputControlVariableValue"
                for variable_name in variable_names_list:
                    variable_sub_list = GetVariableHeaders(norm_type, variable_name)
                    for header in method_headers:
                        for variable_sub_name in variable_sub_list:
                            msg += " " + variable_sub_name + header
                msg += "\n"
                output_file.file.write(msg)
                output_file.is_variable_headers_written = True

    def __get_model_part(self):
        if (not hasattr(self, "model_part")):
            self.model_part = self.model[self.model_part_name]

        return self.model_part

    def __is_writing_process(self):
        return (self.__get_model_part().GetCommunicator().MyPID() == 0)
