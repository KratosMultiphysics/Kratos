from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Other imports
import restart_utility
import os

class TrilinosRestartUtility(RestartUtility):
    def __init__(self, model_part, load_settings, save_settings):
        # Construct the base class
        super(TrilinosRestartUtility, self).__init__(model_part, load_settings, save_settings)

    def _get_load_file_name(self):
        problem_path = os.getcwd()
        return os.path.join(problem_path, self.input_filename + '_' + str(KratosMPI.mpi.rank) + "_" + self.input_file_label)

    def _get_save_file_name(self, time):
        return self.input_filename + '_' + str(KratosMPI.mpi.rank) + '_' + str(time)

    def _execute_after_load(self):
        import trilinos_import_model_part_utility
        trilinos_model_part_importer = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, Parameters('''{"model_import_settings":[]}'''))
        trilinos_model_part_importer.CreateCommunicators() # parallel fill communicator




    

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