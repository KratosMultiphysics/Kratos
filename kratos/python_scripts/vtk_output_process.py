import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
import os

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcess(model, settings["Parameters"])


class VtkOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

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
                folder_name = settings["folder_name"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(folder_name)
                if not os.path.isdir(folder_name):
                    os.mkdir(folder_name)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        self.output_interval = settings["output_interval"].GetDouble()
        self.output_control = settings["output_control_type"].GetString()
        self.next_output = 0.0

        self.__ScheduleNextOutput() # required here esp for restart

    @staticmethod
    def HasDeprecatedVariable(param, old_variable_name, new_variable_name):

        if param.Has(old_variable_name):
            if not param.Has(new_variable_name):
                KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m(DEPRECATED INPUT PARAMETERS)\x1b[0m',
                    'The input variable \'' + old_variable_name + '\' is deprecated; use \'' + new_variable_name + '\' instead.'
                    )
                return True
            else:
                raise NameError('Conflicting input variable names: Both the deprecated variable \''
                                + old_variable_name + '\' and its current standard replacement \''
                                + new_variable_name + '\' were found. Please, remove \'' + old_variable_name + '\'.')
        return False

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, param):

        old_name = 'output_frequency'
        new_name = 'output_interval'

        if VtkOutputProcess.HasDeprecatedVariable(param, old_name, new_name):
            param.AddEmptyValue(new_name)
            if param.Has('output_control_type'):
                control_type = param['output_control_type'].GetString()
            else:
                control_type = self.defaults['output_control_type'].GetString()

            if control_type == 'step':
                param[new_name].SetInt(param[old_name].GetInt())
            else:
                param[new_name].SetDouble(param[old_name].GetDouble())

            param.RemoveValue(old_name)

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
        # remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))
