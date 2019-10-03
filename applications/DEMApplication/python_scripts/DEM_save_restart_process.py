from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics.save_restart_process import SaveRestartProcess
from KratosMultiphysics.DEMApplication import DEM_restart_utility
import json

def Factory(settings, Model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return DEMSaveRestartProcess(Model, settings["Parameters"])

class DEMSaveRestartProcess(SaveRestartProcess):
    """This process compares saves restart files
    It works both in OpenMP and MPI
    see the "default_settings" for available options
    """
    def __init__(self, model, params):
        Kratos.Process.__init__(self)
        ## Settings string in json format
        default_settings = Kratos.Parameters("""{
            "help"                         : "This process is used in order to save the problem databse with the serializer the current problem",
            "model_part_names"              : ["SPECIFY_MODEL_PART_NAMES"],
            "echo_level"                   : 0,
            "serializer_trace"             : "no_trace",
            "restart_save_frequency"       : 0.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : true,
            "remove_restart_folder_after_reading" : false
        }""")

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_settings)
        self.remove_files = params["remove_restart_folder_after_reading"].GetBool()
        params.RemoveValue("help")
        params.RemoveValue("remove_restart_folder_after_reading")
        params_dict = json.loads(params.PrettyPrintJsonString())
        params_dict["input_filenames"] = [mp.GetString() for mp in params["model_part_names"]]
        params = Kratos.Parameters(json.dumps(params_dict))

        params.RemoveValue("model_part_names")

        self.restart_utility = DEM_restart_utility.DEMRestartUtility(model, params)

        # already create the folder now to avoid problems on slow file-systems
        self.restart_utility.CreateOutputFolder()
