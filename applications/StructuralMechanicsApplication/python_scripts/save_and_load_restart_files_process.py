from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
# Import utilities
from KratosMultiphysics.restart_utility import RestartUtility

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SaveAndLoadRestartFilesProcess(Model, settings["Parameters"])


class SaveAndLoadRestartFilesProcess(KratosMultiphysics.Process):
    """This process compares saves restart files
    It works both in OpenMP and MPI
    see the "default_settings" for available options
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""{
            "input_filename"               : "test_restart_file",
            "help"                         : "This process is used in order to save the problem databse with the serializer the current problem",
            "model_part_name"              : "SPECIFY_MODEL_PART_NAME",
            "echo_level"                   : 0,
            "serializer_trace"             : "no_trace",
            "restart_save_frequency"       : 1.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : false,
            "load_restart_files_from_folder" : false,
            "save_restart_files"           : false,
            "load_restart_files"           : false,
            "restart_load_file_label"        : ""
        }""")

        ## Overwrites the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_settings)
        params.RemoveValue("help")

        self.model_part = model[params["model_part_name"].GetString()]
        self.save_model_part = params["save_restart_files"].GetBool()
        self.load_model_part = params["load_restart_files"].GetBool()
        self.model_part_name = params["model_part_name"].GetString()

        params.RemoveValue("model_part_name")
        params.RemoveValue("save_restart_files")
        params.RemoveValue("load_restart_files")

        self.parameters = params

    def ExecuteInitializeSolutionStep(self):
        if self.load_model_part == True:
            loaded_model = KratosMultiphysics.Model()
            self.loaded_model_part = loaded_model.CreateModelPart(self.model_part_name)
            self.parameters["restart_load_file_label"].SetString( str(round(self.model_part.ProcessInfo[KratosMultiphysics.TIME], 3) ))
            rest_utility_load = RestartUtility(self.loaded_model_part, self.parameters)
            rest_utility_load.LoadRestart()
            StructuralMechanicsApplication.ReplaceElementsWithSerializedElementsProcess(self.model_part, self.loaded_model_part).Execute()

    def ExecuteFinalizeSolutionStep(self):
        if self.save_model_part == True:
            self.restart_utility = RestartUtility(self.model_part, self.parameters)
            self.restart_utility.SaveRestart()
