from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class RestartUtility(object):
    """
    This class collects the common functionalities needed for
    saving / loading restart files.

    It can either be integrated into python-solvers or used directly
    in the main-script
    """
    def __init__(self, model_part, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "input_filename"                 : "",
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

        __serializer_flags = {
            "no_trace":           KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE,     # binary
            "trace_error":        KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR,  # ascii
            "trace_all":          KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL     # ascii
        }

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model_part
        self.model_part_name = model_part.Name

        # the path is splitted in case it already contains a path (neeeded if files are moved to a folder)
        self.raw_path, self.raw_file_name = os.path.split(settings["input_filename"].GetString())
        self.raw_path = os.path.join(os.getcwd(), self.raw_path)
        self.folder_name = self.raw_file_name + "__restart_files"

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

        self.save_restart_files_in_folder = settings["save_restart_files_in_folder"].GetBool()

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
        if (os.path.exists(restart_path+".rest") == False):
            raise Exception("Restart file not found: " + restart_path + ".rest")
        self._PrintOnRankZero("::[Restart Utility]::", "Loading restart file:", restart_path + ".rest")

        # Load the ModelPart
        serializer = KratosMultiphysics.Serializer(restart_path, self.serializer_flag)
        serializer.Load(self.model_part_name, self.model_part)

        self._ExecuteAfterLoad()

        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        self._PrintOnRankZero("::[Restart Utility]::", "Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        Use "IsRestartOutputStep" to check if a restart file should be written in this time-step
        """
        if self.save_restart_files_in_folder:
            folder_path = self.__GetFolderPathSave()
            if not os.path.isdir(folder_path) and self.model_part.GetCommunicator().MyPID() == 0:
                os.makedirs(folder_path)
            self.model_part.GetCommunicator().Barrier()

        if self.restart_control_type_is_time:
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            control_label = self.__GetPrettyTime(time)
        else:
            control_label = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        file_name = self.__GetFileNameSave(control_label)

        # Save the ModelPart
        serializer = KratosMultiphysics.Serializer(file_name, self.serializer_flag)
        serializer.Save(self.model_part.Name, self.model_part)
        if self.echo_level > 0:
            self._PrintOnRankZero("::[Restart Utility]::", "Saved restart file", file_name + ".rest")

        # Schedule next output
        if self.restart_save_frequency > 0.0: # Note: if == 0, we'll just always print
            while self.next_output <= control_label:
                self.next_output += self.restart_save_frequency

    def IsRestartOutputStep(self):
        """
        This function checks and returns whether a restart file should be written in this time-step
        """
        if self.restart_control_type_is_time:
            return (self.model_part.ProcessInfo[KratosMultiphysics.TIME] > self.next_output)
        else:
            return (self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output)

    #### Protected functions ####

    def _GetFileLabelLoad(self):
        return self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(file_label)

    def _ExecuteAfterLoad(self):
        """This function creates the communicators in MPI/trilinos"""
        pass

    def _PrintOnRankZero(self, *args):
        # This function will be overridden in the trilinos-version
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    #### Private functions ####

    def __GetFolderPathLoad(self):
        if self.load_restart_files_from_folder:
            return os.path.join(self.raw_path, self.folder_name)
        else:
            return self.raw_path

    def __GetFolderPathSave(self):
        if self.save_restart_files_in_folder:
            return os.path.join(self.raw_path, self.folder_name)
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
