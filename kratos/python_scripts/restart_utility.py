from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class RestartUtility(object):
    def __init__(self, model_part, load_settings, save_settings):
        self.model_part = model_part
        self.input_filename = load_settings["input_filename"].GetString()
        self.input_file_label = load_settings["input_file_label"].GetDouble()
        self.output_frequency = save_settings["restart_time_frequency"].GetDouble()
        self.next_output = 0.0

    def LoadRestart(self):
        # Get file name
        restart_path = self._get_load_file_name()
        # Check path
        if (os.path.exists(restart_path+".rest") == False):
            raise Exception("Restart file not found: " + restart_path + ".rest")        
        print("    Loading Restart file: ", restart_path + ".rest")

        # Load the ModelPart
        serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE      # binary
        # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
        # serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii
        serializer = KratosMultiphysics.Serializer(restart_path, serializer_flag)
        serializer.Load(self.model_part.Name, self.model_part)

        self._execute_after_load()

        # Set restart-info in ProcessInfo
        print(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        err
        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        print("    Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        """
        if self.is_restart_output_step():
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            file_name = self._get_save_file_name(time)
            serializer = KratosMultiphysics.Serializer(file_name)
            serializer.Save(self.model_part.Name, self.model_part)

            # Schedule next output
            if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
                while self.next_output <= time:
                    self.next_output += self.output_frequency

    def is_restart_output_step(self):
        return (self.model_part.ProcessInfo[KratosMultiphysics.TIME] > self.next_output)

    def _get_load_file_name(self):
        problem_path = os.getcwd()
        return os.path.join(problem_path, self.input_filename + "_" + str(self.input_file_label))

    def _get_save_file_name(self, time):
        return self.input_filename + '_' + str(time)

    def _execute_after_load(self):
        """This function creates the communicators in MPI/trilinos"""        
        pass
