from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics.restart_utility import RestartUtility
from KratosMultiphysics import kratos_utilities
# Other imports
import os

class DEMRestartUtility(RestartUtility):
    """
    This class collects the common functionalities needed for
    saving / loading restart files.

    It can either be integrated into python-solvers or used directly
    in the main-script
    """
    def __init__(self, model, settings, restart_save_location='', restart_load_location=''):
        default_settings = Kratos.Parameters("""
        {
            "input_filenames"                 : [],
            "echo_level"                     : 0,
            "serializer_trace"               : "no_trace",
            "restart_load_file_label"        : "",
            "load_restart_files_from_folder" : true,
            "restart_save_frequency"         : 0.0,
            "restart_control_type"           : "time",
            "save_restart_files_in_folder"   : true,
            "set_mpi_communicator"           : true
        }
        """)

        settings.ValidateAndAssignDefaults(default_settings)
        self.file_names = []
        self.restart_utilities = dict()
        self.model_parts = dict()
        for name in settings["input_filenames"]:
            self.file_names.append(name.GetString())
            settings_copy = settings.Clone()
            settings_copy.AddValue("input_filename", name)
            settings_copy.RemoveValue("input_filenames")
            model_part = model.GetModelPart(name.GetString())
            self.restart_utilities[name.GetString()] = RestartUtility(model_part, settings_copy)
            self.model_parts[name.GetString()] = model_part
        self.restart_save_location = restart_save_location
        super(DEMRestartUtility, self).__init__(list(self.model_parts.values())[0], settings_copy)
        # self.restart_load_location = restart_load_location

    def SaveRestart(self):
        for name in self.file_names:
            restart_utility = self.restart_utilities[name]

            if restart_utility.IsRestartOutputStep():
                restart_utility.raw_path, restart_utility.raw_file_name = os.path.split(name)
                #self.raw_path = os.path.join(os.getcwd(), self.raw_path)
                restart_utility.raw_path = os.path.join(os.getcwd(), self.raw_path)
                restart_utility.SaveRestart()

    def LoadRestart(self,  restart_file_name=""):
        for name in self.file_names:
            restart_utility = self.restart_utilities[name]
            restart_utility.raw_path, restart_utility.raw_file_name = os.path.split(name)
            restart_utility.raw_path = os.path.join(os.getcwd(), self.raw_path)

            restart_utility.LoadRestart()

            kratos_utilities.DeleteDirectoryIfExisting(restart_utility._RestartUtility__GetFolderPathLoad())
