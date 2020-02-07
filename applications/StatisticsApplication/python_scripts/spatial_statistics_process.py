import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as Statistics
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetMethod
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetMethodHeaders
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility


def Factory(settings, model):
    if (type(model) != Kratos.Model):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (type(settings) != Kratos.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SpatialStatisticsProcess(model, settings["Parameters"])


class SpatialStatisticsProcess(Kratos.Process):
    def __init__(self, model, settings):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_settings" : [
                {
                    "method_name"   : "sum",
                    "norm_type"     : "none",
                    "container"     : "nodal_historical",
                    "variable_names": []
                }
            ],
            "output_settings" : {
                "output_control_variable": "STEP",
                "output_time_interval"   : 1,
                "output_to_screen"       : true,
                "output_file_settings" : {}
            }
        }  """)

        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.variables_settings_list = self.settings["input_variable_settings"]

        for variable_settings in self.variables_settings_list:
            variable_settings.ValidateAndAssignDefaults(default_parameters["input_variable_settings"][0])
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()

            item_container = GetItemContainer(container_name)
            item_norm_container = GetNormTypeContainer(item_container, norm_type)
            method = GetMethod(item_norm_container, method_name)

        self.model_part_name = self.settings["model_part_name"].GetString()

        output_control_variable_name = self.settings["output_settings"]["output_control_variable"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(output_control_variable_name)):
            raise Exception("Unknown output control variable. [ \"output_control_variable\" = \"" + output_control_variable_name + "\" ]")

        Kratos.Logger.PrintInfo("SpatialStatisticsProcess", "Initialized statistics process.")

    def Check(self):
        if (not self.model.HasModelPart(self.model_part_name)):
            raise Exception(self.model_part_name + " not found.")

        self.model_part = self.model[self.model_part_name]

        for variable_settings in self.variables_settings_list:
            variable_names_list = variable_settings["variable_names"].GetStringArray()
            for variable_name in variable_names_list:
                if (not Kratos.KratosGlobals.HasVariable(variable_name)):
                    raise Exception("Variable not found. [ variable_name = " + variable_name + " ]")

            container = variable_settings["container"].GetString()
            if (container == "nodal_historical"):
                for variable_name in variable_names_list:
                    variable = Kratos.KratosGlobals.GetVariable(variable_name)
                    if (not self.model_part.HasNodalSolutionStepVariable(variable)):
                        raise Exception("Variable " + variable_name + " not found in nodal solution step data of " + self.model_part.Name)

        process_info = self.model_part.ProcessInfo
        output_settings = self.settings["output_settings"]
        output_control_variable_name = output_settings["output_control_variable"].GetString()
        output_control_variable = Kratos.KratosGlobals.GetVariable(output_control_variable_name)
        output_control_variable_type = Kratos.KratosGlobals.GetVariableType(output_control_variable_name)

        if (not Kratos.KratosGlobals.HasVariable(output_control_variable_name)):
            raise Exception("Unknown output control variable " + output_control_variable_name)

        if (output_control_variable_type not in ["Integer", "Double"]):
            raise Exception("Unsupported output control variable type for " + output_control_variable_name + " of " + output_control_variable_type + " type. Supported types are Integer and Double only")

        if (not process_info.Has(output_control_variable)):
            raise Exception("Output control variable " + output_control_variable_name + " not found in process info of " + self.model_part_name)

        process_info_value = self.model_part.ProcessInfo[output_control_variable]
        self.output_control_counter = process_info_value
        self.previous_process_info_value = process_info_value

    def ExecuteFinalizeSolutionStep(self):
        output_settings = self.settings["output_settings"]
        output_control_variable_name = output_settings["output_control_variable"].GetString()
        output_control_variable = Kratos.KratosGlobals.GetVariable(output_control_variable_name)
        output_control_variable_type = Kratos.KratosGlobals.GetVariableType(output_control_variable_name)

        if (output_control_variable_type == "Integer"):
            current_output_control_variable_value = output_settings["output_time_interval"].GetInt()
        elif (output_control_variable_type == "Double"):
            current_output_control_variable_value = output_settings["output_time_interval"].GetDouble()

        process_info_value = self.model_part.ProcessInfo[output_control_variable]

        self.output_control_counter += (process_info_value - self.previous_process_info_value)

        if (self.output_control_counter >= current_output_control_variable_value):
            self.CalculateOutput()
            self.output_control_counter = 0

        self.previous_process_info_value = process_info_value

    def CalculateOutput(self):
        for variable_settings in self.variables_settings_list:
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()
            method_headers = GetMethodHeaders(method_name)

            item_container = GetItemContainer(container_name)
            item_norm_container = GetNormTypeContainer(item_container, norm_type)
            method = GetMethod(item_norm_container, method_name)

            variable_list = []
            max_variable_length = 0
            for variable_name in variable_settings["variable_names"].GetStringArray():
                max_variable_length = max(max_variable_length, len(variable_name))
                variable_list.append(Kratos.KratosGlobals.GetVariable(variable_name))

            max_variable_length = self.__write_headers(method_headers, container_name, method_name, norm_type, max_variable_length + 3)
            if (norm_type == "none"):
                for variable in variable_list:
                    output = method(self.model_part, variable)
                    self.__write_values(variable, output, max_variable_length, method_headers)
            else:
                for variable in variable_list:
                    output = method(norm_type, self.model_part, variable)
                    self.__write_values(variable, output, max_variable_length, method_headers)


    def __write_headers(self, method_headers, container_name, method_name, norm_type, max_variable_length):
        if (self.settings["output_settings"]["output_to_screen"].GetBool()):
            Kratos.Logger.PrintInfo("SpatialStatisticsProcess", "Spatial statistical results for " + self.model_part_name + "'s " + container_name + " container under method " + method_name + " using " + norm_type + " norm:")
            msg = "Variable Name"
            max_variable_length = max(len(msg) + 3, max_variable_length)
            msg = msg.rjust(max_variable_length)
            for header_name in method_headers:
                msg += "  " + header_name.ljust(30)
            Kratos.Logger.PrintInfo("SpatialStatisticsProcess", msg)
            return max_variable_length

    def __write_values(self, variable, output, max_variable_length, method_headers):
        if (self.settings["output_settings"]["output_to_screen"].GetBool()):
            msg = ""
            if len(method_headers) == 1:
                msg += variable.Name().rjust(max_variable_length)
                msg += str(output)
            else:
                msg += variable.Name().rjust(max_variable_length)
                for value in output:
                    msg += "  " + str(value).ljust(30)
            Kratos.Logger.PrintInfo("SpatialStatisticsProcess", msg)