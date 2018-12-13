import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkOutputProcessPython(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class VtkOutputProcessPython(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                    : "PLEASE_SPECIFY_MOEL_PART_NAME",
            "file_name"                          : "",
            "output_control_type"                : "step",
            "output_frequency"                   : 1.0,
            "save_output_files_in_folder"        : true,
            "nodal_solution_step_data_variables" : [],
            "nodal_data_value_variables"         : [],
            "element_data_value_variables"       : []
        }
        """)

        self.model = Model

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        if settings["save_output_files_in_folder"].GetBool() :
            os.mkdir("VTK_Output")

        model_part_name = settings["model_part_name"].GetString()
        model_part = self.model[model_part_name]

        self.cpp_process = KratosMultiphysics.VtkOutputProcess(model_part, self.settings)

    def ExecuteInitialize(self):
        self.cpp_process.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        self.cpp_process.ExecuteFinalizeSolutionStep()
