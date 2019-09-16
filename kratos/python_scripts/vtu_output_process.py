import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
import os

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtuOutputProcess(model, settings["Parameters"])


class VtuOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        # default settings can be found in "vtu_output.cpp"
        self.vtu_output = KratosMultiphysics.VtuOutput(self.model_part, settings) # this also validates the settings

        if settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                folder_name = settings["folder_name"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(folder_name)
                if not os.path.isdir(folder_name):
                    os.mkdir(folder_name)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        self.output_frequency = settings["output_frequency"].GetDouble()
        self.output_control = settings["output_control_type"].GetString()
        self.next_output = 0.0

        self.__ScheduleNextOutput() # required here esp for restart

    def PrintOutput(self):
        self.vtu_output.PrintOutput()

        self.__ScheduleNextOutput()

    def IsOutputStep(self):
        if self.output_control == "time":
            return self.__GetTime() >= self.next_output
        else:
            return self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output

    def __ScheduleNextOutput(self):
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while self.next_output <= self.__GetTime():
                    self.next_output += self.output_frequency
            else:
                while self.next_output <= self.model_part.ProcessInfo[KratosMultiphysics.STEP]:
                    self.next_output += self.output_frequency

    def __GetTime(self):
        # remove rounding errors that mess with the comparison
        # e.g. 1.99999999999999999 => 2.0
        return float("{0:.12g}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME]))
