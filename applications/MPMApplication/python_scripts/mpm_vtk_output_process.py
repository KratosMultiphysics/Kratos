import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMVtkOutputProcess(model, settings["Parameters"])

class MPMVtkOutputProcess(KratosMultiphysics.OutputProcess):

    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()

        self.model = model
        self.settings = settings

        self.full_model_part_name = self.settings["model_part_name"].GetString()
        if self.full_model_part_name.startswith('Background_Grid'):
            self.full_model_part_name = self.full_model_part_name.replace('Background_Grid','MPM_Material')
        main_model_part_name = self.full_model_part_name.split(".",1)[0]
        self.main_model_part = self.model[main_model_part_name]

        # Change name deprecated settings
        self.TranslateLegacyVariablesAccordingToCurrentStandard(self.settings)
        # Validate settings using default parameters defined in MPMVtkOutput process
        # Default settings can be found in "custom_io/mpm_vtk_output.cpp"
        default_settings = KratosMPM.MPMVtkOutput.GetDefaultParameters()
        self.settings.ValidateAndAssignDefaults(default_settings)

        if self.settings["save_output_files_in_folder"].GetBool():
            if self.main_model_part.GetCommunicator().MyPID() == 0:
                output_path = self.settings["output_path"].GetString()
                if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(output_path)
            self.main_model_part.GetCommunicator().GetDataCommunicator().Barrier()

        # Print background grid model part
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            grid_settings = KratosMultiphysics.Parameters()
            grid_settings.AddString("model_part_name","Background_Grid")
            grid_settings.AddString("file_format",self.settings["file_format"].GetString())
            grid_settings.AddDouble("output_precision",self.settings["output_precision"].GetDouble())
            grid_settings.AddBool("output_sub_model_parts",self.settings["output_sub_model_parts"].GetBool())
            grid_settings.AddString("output_path",self.settings["output_path"].GetString())
            grid_settings.AddBool("save_output_files_in_folder",self.settings["save_output_files_in_folder"].GetBool())
            background_grid = self.model["Background_Grid"]
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

        old_name = 'gauss_point_results'
        new_name = 'gauss_point_variables_in_elements'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def ExecuteBeforeSolutionLoop(self) -> None:
        vtk_model_part = self.model[self.full_model_part_name]
        self.vtk_io = KratosMPM.MPMVtkOutput(vtk_model_part, self.settings)

        self.__controller = KratosMultiphysics.OutputController(self.model, self.settings)

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.vtk_io.PrintOutput()

    def PrintOutput(self) -> None:
        self.vtk_io.PrintOutput()
        self.__controller.Update()

    def IsOutputStep(self) -> bool:
        return self.__controller.Evaluate()
