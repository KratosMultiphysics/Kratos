import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcess(model, settings["Parameters"])


class VtkOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.OutputProcess.__init__(self)

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

        self.__controller = KratosMultiphysics.OutputController(model, settings)

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings: KratosMultiphysics.Parameters) -> None:
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

    def Check(self) -> int:
        return self.__controller.Check()

    def PrintOutput(self) -> None:
        self.vtk_io.PrintOutput()
        self.__controller.Update()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()
