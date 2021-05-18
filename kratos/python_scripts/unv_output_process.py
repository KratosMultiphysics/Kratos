import KratosMultiphysics
from KratosMultiphysics import kratos_utilities
from  KratosMultiphysics.deprecation_management import DeprecationManager
import os

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return UnvOutputProcess(model, settings["Parameters"])


class UnvOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # IMPORTANT: when "output_control_type" is "time",
        # then paraview will not be able to group them
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_interval"                    : 1.0,
            "output_control_type"                : "step",
            "nodal_results"                      : []
        }""")

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        self.settings = settings

        # Warning: we may be changing the parameters object here:
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)

        self.settings.ValidateAndAssignDefaults(default_parameters)

        # if self.settings["save_output_files_in_folder"].GetBool():
        #     if self.model_part.GetCommunicator().MyPID() == 0:
        #         folder_name = self.settings["folder_name"].GetString()
        #         if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
        #             import KratosMultiphysics.kratos_utilities as kratos_utils
        #             kratos_utils.DeleteDirectoryIfExisting(folder_name)
        #         if not os.path.isdir(folder_name):
        #             os.mkdir(folder_name)
        #     self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        # self.unv_io = KratosMultiphysics.UnvOutput(self.model_part, self.settings)
        self.unv_io = KratosMultiphysics.UnvOutput(self.model_part, "nxout")
        self.unv_io.InitializeMesh()
        self.unv_io.WriteMesh()

        self.output_interval = self.settings["output_interval"].GetDouble()
        self.output_control = self.settings["output_control_type"].GetString()
        self.next_output = 0.0
        self.step_count = 0

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):

        context_string = type(self).__name__

        old_name = 'output_frequency'
        new_name = 'output_interval'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def ExecuteInitialize(self):
        if self.output_control == "time":
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        self.nodal_variables = kratos_utilities.GenerateVariableListFromInput(self.settings["nodal_results"])

    def ExecuteInitializeSolutionStep(self):
        self.step_count += 1

    def PrintOutput(self):
        for variable in self.nodal_variables:
            self.unv_io.PrintOutput(variable, self.next_output)

        # Schedule next output
        time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while GetPrettyTime(self.next_output) <= time:
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_interval

    def IsOutputStep(self):
        if self.output_control == "time":
            time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= GetPrettyTime(self.next_output))
        else:
            return ( self.step_count >= self.next_output )


def GetPrettyTime(time):
    pretty_time = "{0:.12g}".format(time)
    pretty_time = float(pretty_time)
    return pretty_time
