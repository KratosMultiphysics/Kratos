from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class RestartUtility(object):
    def __init__(self, model_part, settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "input_filename"          : "",
            "load_restart"            : false,
            "restart_load_file_label" : 0.0,
            "save_restart"            : false,
            "restart_save_frequency"  : 0.0,
            "serializer_trace"        : "no_trace"
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
        self.input_file_label = settings["restart_load_file_label"].GetDouble()

        # save settings
        self.output_frequency = settings["restart_save_frequency"].GetDouble()
        serializer_trace = settings["serializer_trace"].GetString()

        serializer_trace = settings["serializer_trace"].GetString()
        if not serializer_trace in __serializer_flags.keys():
            err_msg =  'The requested load serializer_trace "' + serializer_trace + '" is not available!\n'
            err_msg += 'Available options are: "no_trace", "trace_error", "trace_all"'
            raise Exception(err_msg)
        self.serializer_flag = __serializer_flags[serializer_trace]

        self.next_output = 0.0

    def LoadRestart(self):
        # Get file name
        restart_path = self._GetFileNameLoad()
        # Check path
        if (os.path.exists(restart_path+".rest") == False):
            raise Exception("Restart file not found: " + restart_path + ".rest")
        self._PrintOnRankZero("Restart Utility", "Loading Restart file: ", restart_path + ".rest")

        # Load the ModelPart
        serializer = KratosMultiphysics.Serializer(restart_path, self.serializer_flag)
        serializer.Load(self.model_part.Name, self.model_part)

        self._ExecuteAfterLaod()

        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        self._PrintOnRankZero("Restart Utility", "Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        """
        if self.__IsRestartOutputStep():
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            file_name = self._GetFileNameSave(time)

            # Save the ModelPart
            serializer = KratosMultiphysics.Serializer(file_name, self.serializer_flag)
            serializer.Save(self.model_part.Name, self.model_part)

            # Schedule next output
            if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
                while self.next_output <= time:
                    self.next_output += self.output_frequency

    def _GetFileNameLoad(self):
        problem_path = os.getcwd()
        return os.path.join(problem_path, self.input_filename + "_" + str(self.input_file_label))

    def _GetFileNameSave(self, time):
        return self.input_filename + '_' + str(time)

    def _ExecuteAfterLaod(self):
        """This function creates the communicators in MPI/trilinos"""
        pass

    def _PrintOnRankZero(self, *args):
        # This function will be overridden in the trilinos-version
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    def __IsRestartOutputStep(self):
        return (self.model_part.ProcessInfo[KratosMultiphysics.TIME] > self.next_output)