from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SaveRestartProcess(Model, settings["Parameters"])

class SaveRestartProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        """This process compares files that are written during a simulation
        against reference files.
        Please see the "ExecuteFinalize" functions for details about the
        available file-formats
        """
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"              : "",
            "parallel_type"                : "OpenMP",
            "echo_level"                   : 0,
            "serializer_trace"             : "no_trace",
            "restart_save_frequency"       : 0.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : true
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)
        model_part_name = params["model_part_name"].GetString()

        parallel_type = params["parallel_type"].GetString()

        if parallel_type == "OpenMP":
            from restart_utility import RestartUtility as Restart
        elif parallel_type == "MPI":
            import KratosMultiphysics.TrilinosApplication
            from trilinos_restart_utility import TrilinosRestartUtility as Restart
        else:
            err

        params.RemoveValue("model_part_name")
        params.RemoveValue("parallel_type")

        self.restart_utility = Restart(model[model_part_name], params)

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        self.restart_utility.SaveRestart()

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass