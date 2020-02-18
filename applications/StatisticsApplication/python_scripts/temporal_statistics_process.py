import KratosMultiphysics as Kratos

from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod
from KratosMultiphysics.StatisticsApplication.temporal_utilities import GetItemContainer

def Factory(settings, model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return TemporalStatisticsProcess(model, settings["Parameters"])


class TemporalStatisticsProcess(Kratos.Process):
    def __init__(self, model, settings):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_settings" : [
                {
                    "method_name"     : "sum",
                    "norm_type"       : "none",
                    "container"       : "nodal_historical_non_historical",
                    "echo_level"      : 0,
                    "method_settings" : {}
                }
            ],
            "echo_level" : 0,
            "statistics_start_point_control_variable_name" : "TIME",
            "statistics_start_point_control_value"         : 0.0
        }  """)

        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        self.echo_level = self.settings["echo_level"].GetInt()
        self.variables_settings_list = self.settings["input_variable_settings"]

        for variable_settings in self.variables_settings_list:
            variable_settings.ValidateAndAssignDefaults(default_parameters["input_variable_settings"][0])
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()

            item_container = GetItemContainer(container_name)
            _ = GetNormTypeContainer(item_container, norm_type)

        self.model_part_name = self.settings["model_part_name"].GetString()

        statistics_control_variable_name = self.settings["statistics_start_point_control_variable_name"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(statistics_control_variable_name)):
            raise Exception("Unknown statistics control variable. [ \"statistics_control_variable_name\" = \"" + statistics_control_variable_name + "\" ]")

        self.statistics_control_variable = Kratos.KratosGlobals.GetVariable(statistics_control_variable_name)
        statistics_control_variable_type = Kratos.KratosGlobals.GetVariableType(statistics_control_variable_name)
        if (statistics_control_variable_type == "Integer"):
            self.statistics_control_value = self.settings["statistics_start_point_control_value"].GetInt()
        elif (statistics_control_variable_type == "Double"):
            self.statistics_control_value = self.settings["statistics_start_point_control_value"].GetDouble()
        else:
            raise Exception("Unsupported statistics control variable type. Only supports Integer, and Double.")

        if (self.echo_level > 0):
            Kratos.Logger.PrintInfo("TemporalStatisticsProcess", "Initialized statistics process.")

    def Check(self):
        if (not self.model.HasModelPart(self.model_part_name)):
            raise Exception(self.model_part_name + " not found.")

    def ExecuteInitialize(self):
        self.method_list = []
        for variable_settings in self.variables_settings_list:
            container_name = variable_settings["container"].GetString()
            norm_type = variable_settings["norm_type"].GetString()
            method_name = variable_settings["method_name"].GetString()
            echo_level = variable_settings["echo_level"].GetInt()

            item_container = GetItemContainer(container_name)
            method = GetMethod(item_container, method_name)
            method_objects = method(self.__get_model_part(), norm_type, echo_level, variable_settings["method_settings"])

            self.method_list.extend(method_objects)

        for method in self.method_list:
            method.InitializeStatisticsMethod()
        self.previous_value = 0.0

    def ExecuteFinalizeSolutionStep(self):
        current_value = self.model_part.ProcessInfo[self.statistics_control_variable]
        if (current_value >= self.statistics_control_value):
            delta_time = current_value - self.previous_value
            for method in self.method_list:
                method.CalculateStatistics(delta_time)
        self.previous_value = current_value

    def __get_model_part(self):
        if (not hasattr(self, "model_part")):
            self.model_part = self.model[self.model_part_name]

        return self.model_part