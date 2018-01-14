from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import os

class RestartUtility(object):
    def __init__(self, model_part, load_settings, save_settings):
        self.model_part = model_part
        self.input_filename = load_settings["input_filename"].GetString()
        self.input_file_label = load_settings["input_file_label"].GetString()
        self.restart_frequency = save_settings["restart_time_frequency"].GetDouble()

    def LoadRestart(self):
        # Get file name
        problem_path = os.getcwd()
        restart_path = os.path.join(problem_path, self.input_filename + "_" + self.input_file_label)
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

        # Set restart-info in ProcessInfo
        self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        load_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + 1
        self.model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = load_step

        print("    Finished loading model part from restart file.")

    def SaveRestart(self):
        """
        This function saves the restart file. It should be called at the end of a time-step.
        """
        if self.is_restart_output_step():
            serializer = Serializer(self.input_filename + '_' + self.input_file_label)
            serializer.Save(self.model_part.Name, self.model_part)

            self.out = 0.0 # TODO check if this indexing is working as expected!
        else:
            self.out += self.model_part.ProcessInfo[DELTA_TIME]


    def is_restart_output_step(self):
        return (self.out >= self.restart_frequency)




    

    # def is_restart_output_step(self):
    ### THIS IS WRONG, JUST COPY-PASTE !!!
    #     if (self.settings["model_import_settings"]["input_type"].GetString() == "rest"):
    #         return True
    #     else:
    #         return False

    '''
    def __init__(self, problem_name, model_part, settings):
        Process.__init__(self)
        self.problem_name = problem_name
        self.model_part = model_part
        self.load_restart = settings["LoadRestart"].GetBool()
        self.restart_step = settings["RestartStep"].GetString()
        self.restart_frequency = settings["RestartFrequency"].GetDouble()

    def IsRestart(self):
        return self.load_restart

    def LoadRestart(self):
        serializer = Serializer(self.problem_name + '_' + str(KratosMPI.mpi.rank) + '_' + self.restart_step)
        serializer.Load('ModelPart', self.model_part)
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, Parameters('''{"model_import_settings":[]}'''))
        TrilinosModelPartImporter.CreateCommunicators() # parallel fill communicator

    def ExecuteInitializeSolutionStep(self):
        self.out += self.model_part.ProcessInfo[DELTA_TIME]

    def ExecuteBeforeSolutionLoop(self):
        self.out = 0.0

    def IsOutputStep(self):
        return (self.out >= self.restart_frequency)

    def WriteRestartFile(self):
        time = self.model_part.ProcessInfo[TIME]
        serializer = Serializer(self.problem_name + '_' + str(KratosMPI.mpi.rank) + '_' + str(time))
        serializer.Save('ModelPart', self.model_part)
        self.out = 0.0
    '''