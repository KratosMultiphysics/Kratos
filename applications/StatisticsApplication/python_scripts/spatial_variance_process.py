import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as Statistics
from KratosMultiphysics.StatisticsApplication.spatial_utilities import CheckDefaultVariableSettings
from KratosMultiphysics.StatisticsApplication.spatial_utilities import CheckVariables
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetMethodsContainer
from KratosMultiphysics.StatisticsApplication.spatial_utilities import SpatialOutput
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
    return SpatialVarianceProcess(model, settings["Parameters"])


class SpatialVarianceProcess(Kratos.Process):
    def __init__(self, model, settings):
        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters("""
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_settings" : {},
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

        CheckDefaultVariableSettings(self.settings["input_settings"])
        self.model_part_name = self.settings["model_part_name"].GetString()

    def Check(self):
        if (not self.model.HasModelPart(self.model_part_name)):
            raise Exception(self.model_part_name + " not found.")

        self.model_part = self.model[self.model_part_name]
        CheckVariables(self.settings["input_settings"], self.model_part)
        self.method = GetMethodsContainer(
            self.settings["input_settings"]).Variance

    def ExecuteFinalizeSolutionStep(self):
        output = SpatialOutput(["mean", "variance"])
        for variable_name in self.settings["input_settings"][
                "variable_names"].GetStringArray():
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            variable_mean, variable_variance = self.method(
                self.model_part, variable)
            output.PrintOutput(
                [variable_name, str(variable_mean), str(variable_variance)])
