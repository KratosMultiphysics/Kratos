from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class RestartUtility(object):
    def __init__(self, model_part, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "input_filename"                : "",
            "load_restart"                  : false,
            "restart_load_file_label"       : "",
            "save_restart"                  : false,
            "restart_save_frequency"        : 0.0,
            "restart_control_type"          : "time",
            "move_restart_files_to_folder"  : true,
            "serializer_trace"              : "no_trace",
            "echo_level"                    : 0
        }
        """)

        __serializer_flags = {
            "no_trace":           KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE,     # binary
            "trace_error":        KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR,  # ascii
            "trace_all":          KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL     # ascii
        }

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model_part

        # load settings
        self.input_filename = settings["input_filename"].GetString()
        self.input_file_label = settings["restart_load_file_label"].GetString()

        # save settings
        self.restart_save_frequency = settings["restart_save_frequency"].GetDouble()
        serializer_trace = settings["serializer_trace"].GetString()

        restart_control_type = settings["restart_control_type"].GetString()
        if restart_control_type == "time":
            self.restart_control_type_is_time = True
        elif restart_control_type == "step":
            self.restart_control_type_is_time = False
        else:
            err_msg =  'The requested restart_control_type "' + restart_control_type + '" is not available!\n'
            err_msg += 'Available options are: "time", "step"'
            raise Exception(err_msg)

        serializer_trace = settings["serializer_trace"].GetString()
        if not serializer_trace in __serializer_flags.keys():
            err_msg =  'The requested serializer_trace "' + serializer_trace + '" is not available!\n'
            err_msg += 'Available options are: "no_trace", "trace_error", "trace_all"'
            raise Exception(err_msg)
        self.serializer_flag = __serializer_flags[serializer_trace]

        self.next_output = 0.0

        self.echo_level = settings["echo_level"].GetDouble()

        self.move_restart_files_to_folder = settings["move_restart_files_to_folder"].GetBool()

    def LoadRestart(self):
        # Get file name
        restart_path = self._GetFileNameLoad()
        # Check path
        if (os.path.exists(restart_path+".rest") == False):
            raise Exception("Restart file not found: " + restart_path + ".rest")
        self._PrintOnRankZero("::[Restart Utility]::", "Loading restart file:", restart_path + ".rest")

        # Load the ModelPart
        serializer = KratosMultiphysics.Serializer(restart_path, self.serializer_flag)
        serializer.Load(self.model_part.Name, self.model_part)

        self._ExecuteAfterLaod()

        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        self._PrintOnRankZero("::[Restart Utility]::", "Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        """
        if self.__IsRestartOutputStep():
            if self.restart_control_type_is_time:
                control_label = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
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

    def _GetFileNameLoad(self):
        problem_path = os.getcwd()
        return os.path.join(problem_path, self.input_filename + "_" + self._GetFileLoadLabel())

    def _GetFileLoadLabel(self):
        return self.input_file_label

    def _GetFileSaveLabel(self, file_label):
        return str(file_label)

    def _ExecuteAfterLaod(self):
        """This function creates the communicators in MPI/trilinos"""
        pass

    def _PrintOnRankZero(self, *args):
        # This function will be overridden in the trilinos-version
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    def _GetOutputFolderName(self):
        if self.move_restart_files_to_folder:
            folder_name = self.input_filename + "__restart_files"
            if not os.path.isdir(folder_name):
                os.makedirs(folder_name)
            return folder_name
        else:
            return ""

    def __GetFileNameSave(self, file_label):
        return os.path.join(self._GetOutputFolderName(), self.input_filename + '_' + self._GetFileSaveLabel(file_label))

    def __IsRestartOutputStep(self):
        if self.restart_control_type_is_time:
            return (self.model_part.ProcessInfo[KratosMultiphysics.TIME] > self.next_output)
        else:
            return (self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.next_output)