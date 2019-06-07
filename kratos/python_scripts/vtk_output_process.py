import KratosMultiphysics
import os

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcess(model, settings["Parameters"])


class VtkOutputProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # IMPORTANT: when "output_control_type" is "time",
        # then paraview will not be able to group them
        default_parameters = KratosMultiphysics.Parameters("""{
            "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "file_format"                        : "ascii",
            "output_precision"                   : 7,
            "output_control_type"                : "step",
            "output_frequency"                   : 1.0,
            "output_sub_model_parts"             : false,
            "folder_name"                        : "VTK_Output",
            "custom_name_prefix"                 : "",
            "save_output_files_in_folder"        : true,
            "nodal_solution_step_data_variables" : [],
            "nodal_data_value_variables"         : [],
            "nodal_flags"                        : [],
            "element_data_value_variables"       : [],
            "element_flags"                      : [],
            "condition_data_value_variables"     : [],
            "condition_flags"                    : [],
            "gauss_point_variables"              : []
        }""")

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        if self.settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                folder_name = self.settings["folder_name"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    import KratosMultiphysics.kratos_utilities as kratos_utils
                    kratos_utils.DeleteDirectoryIfExisting(folder_name)
                if not os.path.isdir(folder_name):
                    os.mkdir(folder_name)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        self.vtk_io = KratosMultiphysics.VtkOutput(self.model_part, self.settings)

        self.output_frequency = self.settings["output_frequency"].GetDouble()
        self.output_control = self.settings["output_control_type"].GetString()
        self.next_output = 0.0

        self.__ScheduleNextOutput() # required here esp for restart

    def PrintOutput(self):
        self.vtk_io.PrintOutput()

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
