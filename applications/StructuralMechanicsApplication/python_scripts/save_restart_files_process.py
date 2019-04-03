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
    return SaveRestartFilesProcess(Model, settings["Parameters"])


class SaveRestartFilesProcess(KratosMultiphysics.Process):
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
            "restart_save_frequency"       : 0.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : true
        }""")

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_settings)
        params.RemoveValue("help")

        model_part = model[params["model_part_name"].GetString()]

        #params.AddValue("input_filename", params["model_part_name"])
        params.RemoveValue("model_part_name")

        self.restart_utility = RestartUtility(model_part, params)

        self.restart_utility.SaveRestart()

    def ExecuteFinalizeSolutionStep(self):
        self.restart_utility.SaveRestart()
