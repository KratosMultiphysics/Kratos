import KratosMultiphysics
import os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcessPython(Model, settings)


## All the processes python should be derived from "Process"
class VtkOutputProcessPython(KratosMultiphysics.Process):
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
        model_part = self.model[model_part_name]

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)
        folder_name = self.settings["folder_name"].GetString()

        if(self.settings["folder_name"].GetString() == ""):
            self.settings["folder_name"].SetString("VTK_Output_test")

        folder_name = self.settings["folder_name"].GetString()

        if self.settings["save_output_files_in_folder"].GetBool() :
            if(os.path.isdir(folder_name)):
                if(not model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]):
                    import shutil
                    shutil.rmtree(folder_name)
            os.mkdir(folder_name)

        self.cpp_process = KratosMultiphysics.VtkOutput(model_part, self.settings)

    def PrintOutput(self):
        self.cpp_process.PrintOutput()
