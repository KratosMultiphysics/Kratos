# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os
from pathlib import Path

class RestartUtility:
    """
    This class collects the common functionalities needed for
    saving / loading restart files.

    It can either be integrated into python-solvers or used directly
    in the main-script
    """
    def __init__(self, model_part, settings):
        __serializer_flags = {
            "no_trace":           KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE,     # binary
            "trace_error":        KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR,  # ascii
            "trace_all":          KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL     # ascii
        }

        if settings.Has("io_foldername"):
            settings.AddValue("input_output_path",settings["io_foldername"])
            settings.RemoveValue("io_foldername")
            KratosMultiphysics.Logger.PrintWarning('RestartUtility', '"io_foldername" key is deprecated. Use "input_output_path" instead.')

        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())

        self.model_part = model_part
        self.model_part_name = model_part.Name

        # the path is splitted in case it already contains a path (needed if files are moved to a folder)
        self.raw_path, self.raw_file_name = os.path.split(settings["input_filename"].GetString())
        self.raw_path = os.path.join(os.getcwd(), self.raw_path)

        if settings["input_output_path"].GetString() == '':
            self.input_output_path = self.raw_file_name + "__restart_files"
            info_msg  = 'No entry found for "input_output_path"\n'
            info_msg += 'Using the default "' + self.input_output_path + '"'
            KratosMultiphysics.Logger.PrintInfo("RestartUtility", info_msg)

        else:
            self.input_output_path = settings["input_output_path"].GetString()
            info_msg  = 'Found entry found for "input_output_path"\n'
            info_msg += 'Using the user-defined value "' + self.input_output_path + '"'
            KratosMultiphysics.Logger.PrintInfo("RestartUtility", info_msg)

        serializer_trace = settings["serializer_trace"].GetString()
        if not serializer_trace in __serializer_flags.keys():
            err_msg =  'The requested serializer_trace "' + serializer_trace + '" is not available!\n'
            err_msg += 'Available options are: "no_trace", "trace_error", "trace_all"'
            raise Exception(err_msg)
        self.serializer_flag = __serializer_flags[serializer_trace]

        self.echo_level = settings["echo_level"].GetInt()

        # load settings
        self.input_file_label = settings["restart_load_file_label"].GetString()
        self.load_restart_files_from_folder = settings["load_restart_files_from_folder"].GetBool()

        # save settings
        self.restart_save_frequency = settings["restart_save_frequency"].GetDouble()

        restart_control_type = settings["restart_control_type"].GetString()
        if restart_control_type == "time":
            self.restart_control_type_is_time = True
        elif restart_control_type == "step":
            self.restart_control_type_is_time = False
            self.restart_save_frequency = int(self.restart_save_frequency) # STEP is an integer
        else:
            err_msg =  'The requested restart_control_type "' + restart_control_type + '" is not available!\n'
            err_msg += 'Available options are: "time", "step"'
            raise Exception(err_msg)

        self.next_output = self.restart_save_frequency # Schedule the first output to avoid printing in first step

        self.save_restart_files_in_folder   = settings["save_restart_files_in_folder"].GetBool()
        self.max_files_to_keep              = settings["max_files_to_keep"].GetInt()
        if (self.max_files_to_keep < -1) or (self.max_files_to_keep == 0):
            err_msg  = 'Specifier for \'max_files_to_keep\' with value ' + str(self.max_files_to_keep) +' invalid\n'
            err_msg += 'Use -1 or any non-negative values'
            raise Exception(err_msg)


    #### Public functions ####

    def LoadRestart(self,  restart_file_name=""):
        """
        This function loads a restart file into a ModelPart
        """
        if restart_file_name == "": # Using the default restart file name
            # Get file name
            restart_path = self.__GetFileNameLoad()

        else: # Using a custom restart file name
            if restart_file_name.endswith('.rest'):
                restart_file_name = restart_file_name[:-5] # removing ".rest" from the file name
            restart_path = restart_file_name

        # Check path
        if not os.path.exists(restart_path+".rest"):
            raise FileNotFoundError("Restart file not found: " + restart_path + ".rest")
        KratosMultiphysics.Logger.PrintInfo("Restart Utility", "Loading restart file:", restart_path + ".rest")

        # Load the ModelPart
        serializer = KratosMultiphysics.FileSerializer(restart_path, self.serializer_flag)
        serializer.Set(self._GetSerializerFlags())
        serializer.Load(self.model_part_name, self.model_part)

        self._ExecuteAfterLoad()

        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        KratosMultiphysics.Logger.PrintInfo("Restart Utility", "Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        Use "IsRestartOutputStep" to check if a restart file should be written in this time-step
        """
        if self.restart_control_type_is_time:
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            control_label = self.__GetPrettyTime(time)
        else:
            control_label = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        self.CreateOutputFolder()

        if not os.path.isdir(self.__GetFolderPathSave()):
            err_msg  = 'The directory for saving the restart-files of modelpart "'
            err_msg += self.model_part_name + '" does not exist!\n'
            if self.save_restart_files_in_folder:
                err_msg += 'Something went wrong with the creation of the folder "'
                err_msg += self.__GetFolderPathSave()+ '"!'
            raise Exception(err_msg)

        file_name = self.__GetFileNameSave(control_label)

        # Save the ModelPart
        serializer = KratosMultiphysics.FileSerializer(file_name, self.serializer_flag)
        serializer.Set(self._GetSerializerFlags())
        serializer.Save(self.model_part.Name, self.model_part)
        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("Restart Utility", "Saved restart file", file_name + ".rest")

        # Schedule next output
        if self.restart_save_frequency > 0.0: # Note: if == 0, we'll just always print
            while self.next_output <= control_label:
                self.next_output += self.restart_save_frequency

        # Cleanup
        self._ClearObsoleteRestartFiles()

    def IsRestartOutputStep(self):
        """
        This function checks and returns whether a restart file should be written in this time-step
        """
        if self.restart_control_type_is_time:
            return (self.model_part.ProcessInfo[KratosMultiphysics.TIME] > self.next_output)
        else:
            return (self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output)

    def CreateOutputFolder(self):
        if self.save_restart_files_in_folder:
            KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(self.__GetFolderPathSave())

    #### Protected functions ####

    def _GetRestartFiles(self):
        """Return a dictionary of stepID - restart_file_list dictionary that stores sets of restart files for each step."""
        restart_path    = Path(self.__GetFolderPathSave())
        restart_files   = {}
        if restart_path.is_dir():
            file_name_data_collector = KratosMultiphysics.FileNameDataCollector(self.model_part, os.path.join(self.__GetFolderPathSave(), self._GetFileNamePattern()), {})

            for file_name_data in file_name_data_collector.GetFileNameDataList():
                # Get step id
                if self.restart_control_type_is_time:
                    step_id = file_name_data.GetTime()
                else:
                    step_id = file_name_data.GetStep()

                self._UpdateRestartFilesMap(restart_files, step_id, file_name_data)

            # barrier is necessary to avoid having some ranks deleting files while other ranks still detect them in the same directory
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

        return restart_files

    def _UpdateRestartFilesMap(self, restart_files_map, step_id, file_name_data):
        restart_files_map[step_id] = file_name_data.GetFileName()

    def _ClearObsoleteRestartFiles(self):
        """Delete restart files that are no longer needed."""
        if self.max_files_to_keep > -1:
            self.restart_files = self._GetRestartFiles()

            number_of_obsolete_files = len(self.restart_files) - self.max_files_to_keep
            restart_file_keys = sorted(self.restart_files)

            for i in range(number_of_obsolete_files):
                i_key = restart_file_keys[i]
                file_path = os.path.join(self.__GetFolderPathSave(), self.restart_files[i_key])
                KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(file_path)

    def _GetFileLabelLoad(self):
        return self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(file_label)

    def _ExecuteAfterLoad(self):
        """This function creates the communicators in MPI/trilinos"""
        pass

    def _GetFileNamePattern(self):
        """Return the pattern of flags in the file name for FileNameDataCollector."""
        file_name_pattern = "<model_part_name>"
        if self.restart_control_type_is_time:
            file_name_pattern += "_<time>"
        else:
            file_name_pattern += "_<step>"
        file_name_pattern += ".rest"
        return file_name_pattern

    def _GetSerializerFlags(self):
        return KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION

    @classmethod
    def _GetDefaultParameters(cls):
        # max_files_to_keep    : max number of restart files to keep
        #                       - negative value keeps all restart files (default)
        return KratosMultiphysics.Parameters("""{
            "input_filename"                 : "",
            "input_output_path"              : "",
            "echo_level"                     : 0,
            "serializer_trace"               : "no_trace",
            "restart_load_file_label"        : "",
            "load_restart_files_from_folder" : true,
            "restart_save_frequency"         : 0.0,
            "restart_control_type"           : "time",
            "save_restart_files_in_folder"   : true,
            "max_files_to_keep"              : -1
        }""")

    #### Private functions ####

    def __GetFolderPathLoad(self):
        if self.load_restart_files_from_folder:
            return os.path.join(self.raw_path, self.input_output_path)
        else:
            return self.raw_path

    def __GetFolderPathSave(self):
        if self.save_restart_files_in_folder:
            return os.path.join(self.raw_path, self.input_output_path)
        else:
            return self.raw_path

    def __GetFileNameLoad(self):
        restart_file_name = self.raw_file_name + "_" + self._GetFileLabelLoad()

        return os.path.join(self.__GetFolderPathLoad(), restart_file_name)

    def __GetFileNameSave(self, file_label):
        restart_file_name = self.raw_file_name + '_' + self._GetFileLabelSave(file_label)

        return os.path.join(self.__GetFolderPathSave(), restart_file_name)

    def __GetPrettyTime(self, time):
        """This functions reduces the digits of a number to a relevant precision
        """
        pretty_time = "{0:.12g}".format(time)
        pretty_time = float(pretty_time)
        return pretty_time
