import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ParticleVtkOutputProcess(model, settings["Parameters"])

class ParticleVtkOutputProcess(KratosMultiphysics.OutputProcess):

    def __init__(self, model, settings):
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

        self.output_interval = settings["output_interval"].GetDouble()
        self.output_control = settings["output_control_type"].GetString()
        self.next_output = 0.0

        self.__ScheduleNextOutput()

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
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):
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

    def ExecuteBeforeSolutionLoop(self):
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.vtk_io.PrintOutput()

    def PrintOutput(self):
        self.vtk_io.PrintOutput()
        self.__ScheduleNextOutput()

    def IsOutputStep(self):
        if self.output_control == "time":
            return self.__GetTime() >= self.next_output
        else:
            return self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output

    def __ScheduleNextOutput(self):
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while self.next_output <= self.__GetTime():
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
                    self.next_output += self.output_interval

    def __GetTime(self):
        # Remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))
