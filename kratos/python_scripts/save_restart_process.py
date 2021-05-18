# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SaveRestartProcess(Model, settings["Parameters"])


class SaveRestartProcess(KratosMultiphysics.Process):
    """This process compares saves restart files
    It works both in OpenMP and MPI
    see the "default_settings" for available options
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""{
            "model_part_name"              : "SPECIFY_MODEL_PART_NAME",
            "echo_level"                   : 0,
            "serializer_trace"             : "no_trace",
            "restart_save_frequency"       : 0.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : true,
            "output_path"                  : "",
            "max_files_to_keep"            : -1
        }""")

        ## Overwrite the default settings with user-provided parameters
        if params.Has("io_foldername"):
            params.AddValue("output_path",params["io_foldername"])
            params.RemoveValue("io_foldername")
            KratosMultiphysics.Logger.PrintWarning('SaveRestartProcess', '"io_foldername" key is deprecated. Use "output_path" instead.')

        params.ValidateAndAssignDefaults(default_settings)

        model_part = model[params["model_part_name"].GetString()]

        if model_part.IsDistributed(): # mpi-execution
            from KratosMultiphysics.mpi.distributed_restart_utility import DistributedRestartUtility as RestartUtility
        else:
            from KratosMultiphysics.restart_utility import RestartUtility

        if params["output_path"].GetString() == '':
            output_path = params["model_part_name"].GetString() + "__restart_files"
            info_msg  = 'No entry found for "output_path"\n'
            info_msg += 'Using the default "' + output_path + '"'
            KratosMultiphysics.Logger.PrintInfo("SaveRestartProcess", info_msg)

        params.AddValue("input_filename", params["model_part_name"])
        params.RemoveValue("model_part_name")

        params.AddValue("input_output_path",params["output_path"])
        params.RemoveValue("output_path")

        self.restart_utility = RestartUtility(model_part, params)

        # already create the folder now to avoid problems on slow file-systems
        self.restart_utility.CreateOutputFolder()

    def IsOutputStep(self):
        return self.restart_utility.IsRestartOutputStep()

    def PrintOutput(self):
        self.restart_utility.SaveRestart()
