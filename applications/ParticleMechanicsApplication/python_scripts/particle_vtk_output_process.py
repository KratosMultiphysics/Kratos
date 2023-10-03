import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ParticleVtkOutputProcess(model, settings["Parameters"])

class ParticleVtkOutputProcess(KratosMultiphysics.OutputProcess):

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)

        # Validate settings using default parameters defined in ParticleVtkOutput process
        default_settings = KratosParticle.ParticleVtkOutput.GetDefaultParameters()
        settings.ValidateAndAssignDefaults(default_settings)

        # Default settings can be found in "custom_io/particle_vtk_output.cpp"
        self.vtk_io = KratosParticle.ParticleVtkOutput(self.model_part, settings)

        if settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                output_path = settings["output_path"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(output_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        self.__controller = KratosMultiphysics.OutputController(model, settings)

        # Print background grid model part
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            grid_settings = KratosMultiphysics.Parameters()
            grid_settings.AddString("model_part_name","Background_Grid")
            grid_settings.AddString("file_format",settings["file_format"].GetString())
            grid_settings.AddDouble("output_precision",settings["output_precision"].GetDouble())
            grid_settings.AddBool("output_sub_model_parts",settings["output_sub_model_parts"].GetBool())
            grid_settings.AddString("output_path",settings["output_path"].GetString())
            grid_settings.AddBool("save_output_files_in_folder",settings["save_output_files_in_folder"].GetBool())
            background_grid = model["Background_Grid"]
            KratosMultiphysics.VtkOutput(background_grid, grid_settings).PrintOutput()

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings: KratosMultiphysics.Parameters) -> None:
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

        if settings.Has("gauss_point_results"):
            KratosMultiphysics.Logger.PrintWarning("VtkOutputProcess", "Setting `gauss_point_results` is deprecated, use `gauss_point_variables_in_elements` instead!")
            if not settings.Has("gauss_point_variables_in_elements"):
                settings.AddEmptyValue("gauss_point_variables_in_elements").SetStringArray(settings["gauss_point_results"].GetStringArray())
            settings.RemoveValue("gauss_point_results")

    def ExecuteBeforeSolutionLoop(self) -> None:
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.vtk_io.PrintOutput()

    def Check(self) -> int:
        return self.__controller.Check()

    def PrintOutput(self) -> None:
        self.vtk_io.PrintOutput()
        self.__controller.Update()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()
