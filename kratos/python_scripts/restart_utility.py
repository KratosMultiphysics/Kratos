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
        # number_of_restart_files    : max number of restart files to keep
        #                             - negative value keeps all restart files (default)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "input_filename"                 : "",
            "io_foldername"                  : "",
            "echo_level"                     : 0,
            "serializer_trace"               : "no_trace",
            "restart_load_file_label"        : "",
            "load_restart_files_from_folder" : true,
            "restart_save_frequency"         : 0.0,
            "restart_control_type"           : "time",
            "save_restart_files_in_folder"   : true,
            "set_mpi_communicator"           : true,
            "number_of_restart_files"        : -1
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

        if settings["io_foldername"].GetString() == '':
            self.io_foldername = self.raw_file_name + "__restart_files"
            info_msg  = 'No entry found for "io_foldername"\n'
            info_msg += 'Using the default "' + self.io_foldername + '"'
            KratosMultiphysics.Logger.PrintInfo("RestartUtility", info_msg)

        else:
            self.io_foldername = settings["io_foldername"].GetString()
            info_msg  = 'Found entry found for "io_foldername"\n'
            info_msg += 'Using the user-defined value "' + self.io_foldername + '"'
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
        self.number_of_restart_files        = settings["number_of_restart_files"].GetInt()

        # Check if restarted, and initialize the list of restart is true
        self.restart_files = {}
        if self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.restart_files = self.GetRestartFiles()


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
        serializer.Save(self.model_part.Name, self.model_part)
        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("Restart Utility", "Saved restart file", file_name + ".rest")

        # Schedule next output
        if self.restart_save_frequency > 0.0: # Note: if == 0, we'll just always print
            while self.next_output <= control_label:
                self.next_output += self.restart_save_frequency

        # Add current file to stored dictionary
        label = self._GetFileLabelSave(control_label)
        if label in self.restart_files:
            self.restart_files[label].append( file_name + ".rest" )
        else:
            self.restart_files[label] = [ file_name + ".rest" ]

        # Cleanup
        self.ClearObsoleteRestartFiles()

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
            folder_path = self.__GetFolderPathSave()
            if not os.path.isdir(folder_path) and self.model_part.GetCommunicator().MyPID() == 0:
                os.makedirs(folder_path)
            self.model_part.GetCommunicator().GetDataCommunicator().Barrier()

    def GetRestartFiles(self):
        """
        Return a dictionary of stepID - restart_file_list dictionary that stores sets of restart
        files for each step. 
        """
        restart_files = {}
        if os.path.isdir(self.__GetFolderPathSave()):
            number_of_restart_files = 0

            with os.scandir(path=self.__GetFolderPathSave()) as contents:
                for entry in contents:
                    if self.__IsRestartFile(entry):

                        label = self.__ExtractFileLabel(entry.name) # Get thread ID and step ID
                        
                        if label[1] in restart_files:               # Check if this step has entries already
                            restart_files[label[1]].append( entry.name )
                        else:
                            restart_files[label[1]] = [ entry.name ]

                        number_of_restart_files += 1

            # Throw a warning if the number of restart files is too large
            if ( number_of_restart_files > 1000 ):
                message =   "Detected " + str(len(restart_files)) + " restart files. "
                message +=  "Consider restricting the number of restart files to keep."
                print( message )

        return restart_files

    def ClearObsoleteRestartFiles(self):
        """
        Collect all restart files from the current restart directory, sort them by date(time) modified
        and delete the oldest ones such that only number_of_restart_files remain.
        Note: a secondary sorting is performed based on the labels of the file names to resolve equalities.
        """
        if self.number_of_restart_files > -1:               # <-- number of restart files is limited
            number_of_obsolete_files = len(self.restart_files) - self.number_of_restart_files
            for _ in range(number_of_obsolete_files):

                # Get oldest restart file set
                oldest_step_id, oldest_file_set = min( self.restart_files.items(), key=lambda item: float(item[0]) )

                # Try to delete every file in the set
                for file_name in oldest_file_set:
                    file_path = os.path.join( self.__GetFolderPathSave(), file_name )
                    try:
                        if os.path.isfile( file_path ):
                            os.remove( file_path )
                    except OSError:
                        message =   'Failed to delete restart file "'
                        message +=  file_path
                        message += '"'
                        print( message )                # <-- TODO: decide whether to throw a non-blocking exception or display a warning
                        #raise Exception(message)       #
                
                # Update stored dictionary
                del self.restart_files[oldest_step_id]


    #### Protected functions ####

    def _GetFileLabelLoad(self):
        return self.input_file_label

    def _GetFileLabelSave(self, file_label):
        return str(file_label)

    def _ExecuteAfterLoad(self):
        """This function creates the communicators in MPI/trilinos"""
        pass

    #### Private functions ####

    def __GetFolderPathLoad(self):
        if self.load_restart_files_from_folder:
            return os.path.join(self.raw_path, self.io_foldername)
        else:
            return self.raw_path

    def __GetFolderPathSave(self):
        if self.save_restart_files_in_folder:
            return os.path.join(self.raw_path, self.io_foldername)
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

    def __IsRestartFile(self, dir_entry):
        """
        Check whether the input os.DirEntry object refers to a file
        and has an appropriate file name for a restart file.
        """
        if not isinstance( dir_entry, os.DirEntry ):
            message =   "Expecting an "
            message +=  type(os.DirEntry).__name__
            message +=  ", got "
            message +=  type(dir_entry).__name__
            raise ValueError( message )

        if dir_entry.is_file(follow_symlinks=False):
            if dir_entry.name.endswith('.rest') and self.raw_file_name in os.path.basename(dir_entry.name):
                # additional checks might have to be performed here if multiple simulations
                # save their restart files in the same directory.
                return True
        return False

    def __ExtractFileLabel(self, file_name):
        """
        Return a list of labels attached to a file name (list entries separated by '_').
        Expecting one of two possible file name formats:
            1) serial execution     : <file_base_name>_<step_id>.rest
            2) parallel execution   : <file_base_name>_<thread_id>_<step_id>.rest
        The returned list always has 2 entries, but the first component will be '' for
        serial execution.
        """
        label_begin = file_name.find(self.raw_file_name + '_') + len(self.raw_file_name + '_')
        label_end   = file_name.find('.rest')

        labels      = file_name[label_begin:label_end].split("_")
        if len(labels) < 2:
            labels = [""] + labels

        return labels