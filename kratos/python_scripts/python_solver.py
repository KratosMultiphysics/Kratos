from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

# Other imports
import os

class PythonSolver(object):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        settings -- The solver settings used
        """
        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(settings) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.settings = settings

        self.echo_level = self.settings["echo_level"].GetInt()

    def AddVariables(self):
        """This function add the Variables needed by this PythonSolver to the the ModelPart
        It has to be called BEFORE the ModelPart is read!
        """
        pass

    def AddDofs(self):
        """This function add the Dofs needed by this PythonSolver to the the ModelPart
        It has to be called AFTER the ModelPart is read!
        """
        pass

    def ImportModelPart(self):
        """This function reads the ModelPart
        """
        raise Exception("This function has to be implemented in the derived class")

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        pass

    def GetMinimumBufferSize(self):
        """This function returns the minimum buffer size needed for this PythonSolver
        """
        raise Exception('Please implement "GetMinimumBufferSize" in your derived solver')

    def ExportModelPart(self):
        """This function exports the ModelPart to and mdpa-file
        """
        raise Exception("This function has to be implemented in the derived class")

    def AdvanceInTime(self, current_time):
        """This function advances the PythonSolver in time
        Usage: It is designed to be called once per solution step, before performing the solution
        """
        pass

    def Initialize(self):
        """This function initializes the PythonSolver
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        pass

    def Finalize(self):
        """This function finalizes the PythonSolver
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        pass

    def Predict(self):
        """This function performs all the required operations that should be executed
        (for each step) ONCE, AFTER initializing the solution step.
        """
        pass

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        pass

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        pass

    def SolveSolutionStep(self):
        """This function solves the current step.
        It can be called multiple times within one solution step
        """
        pass

    def Check(self):
        """This function checks the PythonSolver. It usually calls the "Check" function of a solving strategy
        """
        pass

    def Clear(self):
        """This function clears the PythonSolver
        """
        pass

    def Solve(self):
        warning_msg  = 'Using "Solve" is deprecated and will be removed in the future!\n'
        warning_msg += 'Use the separate calls to "Initialize", "InitializeSolutionStep", "Predict", '
        warning_msg += '"SolveSolutionStep" and "FinalizeSolutionStep"'
        KratosMultiphysics.Logger.PrintWarning("::[PythonSolver]::", warning_msg)
        self.Initialize()
        self.Predict()
        self.InitializeSolutionStep()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def GetComputingModelPart(self):
        raise Exception("This function has to be implemented in the derived class")

    def _ImportModelPart(self, model_part, model_part_import_settings):
        """This function imports the ModelPart
        """
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part.")
        problem_path = os.getcwd()
        input_filename = model_part_import_settings["input_filename"].GetString()
        input_type = model_part_import_settings["input_type"].GetString()

        if (input_type == "mdpa"):
            # Import model part from mdpa file.
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(model_part)
            if (model_part_import_settings.Has("reorder") and model_part_import_settings["reorder"].GetBool()):
                tmp = KratosMultiphysics.Parameters("{}")
                KratosMultiphysics.ReorderAndOptimizeModelPartProcess(model_part, tmp).Execute()
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished reading model part from mdpa file.")
        elif (input_type == "rest"):
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Loading model part from restart file.")
            from restart_utility import RestartUtility
            RestartUtility(model_part, self._GetRestartSettings(model_part_import_settings)).LoadRestart()
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished loading model part from restart file.")
        else:
            raise Exception("Other model part input options are not yet implemented.")
        KratosMultiphysics.Logger.PrintInfo("ModelPart", model_part)
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]:: ", "Finished reading model part.")

    def _GetRestartSettings(self, model_part_import_settings):
        restart_settings = model_part_import_settings.Clone()
        restart_settings.RemoveValue("input_type")
        if not restart_settings.Has("restart_load_file_label"):
            raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')
        if model_part_import_settings.Has("echo_level"):
            restart_settings.AddValue("echo_level", model_part_import_settings["echo_level"])

        return restart_settings