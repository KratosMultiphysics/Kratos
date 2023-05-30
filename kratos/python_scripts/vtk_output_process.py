import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")

    return VtkOutputProcess(model, settings["Parameters"])

class VtkOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters):
        super().__init__()

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        # Warning: we may be changing the parameters object here:
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)
        if settings.Has("write_properties_id"):
            KratosMultiphysics.Logger.PrintWarning("VtkOutputProcess", "The setting `write_properties_id` is deprecated, use `write_ids` instead!")
            if not settings.Has("write_ids"):
                settings.AddEmptyValue("write_ids").SetBool(settings["write_properties_id"].GetBool())
            settings.RemoveValue("write_properties_id")

        # default settings can be found in "vtk_output.cpp"
        self.vtk_io = KratosMultiphysics.VtkOutput(self.model_part, settings) # this also validates the settings

        if settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                output_path = settings["output_path"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(output_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        output_control_type = settings["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_variable = KratosMultiphysics.TIME
            self.output_control_utility = KratosMultiphysics.DoubleFixedIntervalRecurringEventUtility(self.model_part.ProcessInfo[self.output_control_variable], settings["output_interval"].GetDouble())
        elif output_control_type == "step":
            self.output_control_variable = KratosMultiphysics.STEP
            self.output_control_utility = KratosMultiphysics.IntegerFixedIntervalRecurringEventUtility(self.model_part.ProcessInfo[self.output_control_variable], settings["output_interval"].GetInt())
        else:
            raise RuntimeError(f"Unsupported output control type = \"{output_control_type}\" requested. Supported control types are:\n\ttime\n\tstep")

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        old_name = 'output_frequency'
        new_name = 'output_interval'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        old_name = 'write_properties_id'
        new_name = 'write_ids'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        old_name = 'folder_name'
        new_name = 'output_path'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def PrintOutput(self):
        self.vtk_io.PrintOutput()
        self.output_control_utility.ScheduleNextEvent(self.model_part.ProcessInfo[self.output_control_variable])

    def IsOutputStep(self):
        return self.output_control_utility.IsEventExpected(self.model_part.ProcessInfo[self.output_control_variable])

