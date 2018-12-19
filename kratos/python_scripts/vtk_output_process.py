import KratosMultiphysics
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcess(Model, settings)


## All the processes python should be derived from "Process"
class VtkOutputProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "file_format"                        : "ASCII",
            "output_control_type"                : "step",
            "output_frequency"                   : 1.0,
            "output_sub_model_parts"             : true,
            "folder_name"                        : "",
            "save_output_files_in_folder"        : true,
            "nodal_solution_step_data_variables" : [],
            "nodal_data_value_variables"         : [],
            "element_data_value_variables"       : []
        }
        """)

        self.model = Model
        model_part_name = settings["model_part_name"].GetString()
        self.model_part = self.model[model_part_name]

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        folder_name = self.settings["folder_name"].GetString()

        if(self.settings["folder_name"].GetString() == ""):
            self.settings["folder_name"].SetString("VTK_Output_test")

        folder_name = self.settings["folder_name"].GetString()

        if self.settings["save_output_files_in_folder"].GetBool() :
            if(os.path.isdir(folder_name)):
                if(not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]):
                    import shutil
                    shutil.rmtree(folder_name)
            os.mkdir(folder_name)

        self.cpp_process = KratosMultiphysics.VtkOutput(self.model_part, self.settings)

        self.output_frequency = self.settings["output_frequency"].GetDouble()
        #
        self.output_control = self.settings["output_control_type"].GetString()
        self.next_output = 0.0
        self.step_count = 0

    def ExecuteInitializeSolutionStep(self):
        self.step_count += 1

    def ExecuteInitialize(self):
        if self.output_control == "time":
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

    def PrintOutput(self):
        if(self.IsOutputStep()):
            self.cpp_process.PrintOutput()

            # Schedule next output
            time = self.__GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
                if self.output_control == "time":
                    while self.__GetPrettyTime(self.next_output) <= time:
                        self.next_output += self.output_frequency
                else:
                    while self.next_output <= self.step_count:
                        self.next_output += self.output_frequency

    def IsOutputStep(self):
        if self.output_control == "time":
            time = self.__GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= self.__GetPrettyTime(self.next_output))
        else:
            return ( self.step_count >= self.next_output )


    def __GetPrettyTime(self,time):
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time
