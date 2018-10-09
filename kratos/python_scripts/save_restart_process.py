from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
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
            "help"                         : "This process is used in order to save/load the problem databse with the serializer the current problem",
            "model_part_name"              : "",
            "echo_level"                   : 0,
            "serializer_trace"             : "no_trace",
            "restart_save_frequency"       : 0.0,
            "restart_control_type"         : "time",
            "save_restart_files_in_folder" : true
        }""")

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_settings)
        params.RemoveValue("help")
        self.params = params
        self.model = model

        if self.params["model_part_name"].GetString() == "":
            raise Exception('No "model_part_name" was specified!')

    def ExecuteInitialize(self):
        model_part = self.model[self.params["model_part_name"].GetString()]

        is_mpi_execution = (model_part.GetCommunicator().TotalProcesses() > 1)

        if is_mpi_execution:
            import KratosMultiphysics.TrilinosApplication
            from trilinos_restart_utility import TrilinosRestartUtility as Restart
        else:
            from restart_utility import RestartUtility as Restart

        self.params.AddValue("input_filename", self.params["model_part_name"])
        self.params.RemoveValue("model_part_name")

        self.restart_utility = Restart(model_part, self.params)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def IsOutputStep(self):
        return self.restart_utility.IsRestartOutputStep()

    def PrintOutput(self):
        self.restart_utility.SaveRestart()

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
