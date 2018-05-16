from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

# Other imports
import os

class PythonSolver(object):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model_part, solver_settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(solver_settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model_part -- The ModelPart to be used
        solver_settings -- The solver settings used
        """
        if (type(model_part) != KratosMultiphysics.ModelPart):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(solver_settings) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.main_model_part = model_part
        self.solver_settings = solver_settings

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

    def ReadModelPart(self):
        """This function reads the ModelPart
        """
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part.")
        problem_path = os.getcwd()
        input_filename = self.solver_settings["model_import_settings"]["input_filename"].GetString()
        input_type = self.solver_settings["model_import_settings"]["input_type"].GetString()

        if (input_type == "mdpa"):
            # Import model part from mdpa file.
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished reading model part from mdpa file.")
        elif (input_type == "rest"):
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Loading model part from restart file.")
            from restart_utility import RestartUtility
            RestartUtility(self.main_model_part, self._GetRestartSettings()).LoadRestart()
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished loading model part from restart file.")
        else:
            raise Exception("Other model part input options are not yet implemented.")
        KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]:: ", "Finished reading model part.")

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        pass

    def GetMinimumBufferSize(self):
        """This function returns the minimum buffer size needed for this PythonSolver
        """
        raise Exception('Please implement "GetMinimumBufferSize" in your derived solver')

    def ImportModelPart(self):
        warning_msg  = 'Using "ImportModelPart" is deprecated and will be removed in the future!\n'
        warning_msg  = 'Use "ReadModelPart" + "PrepareModelPart" instead'
        KratosMultiphysics.Logger.PrintWarning("::[PythonSolver]::", warning_msg)
        self.ReadModelPart()
        self.PrepareModelPart()

    def ExportModelPart(self):
        """This function exports the ModelPart to and mdpa-file
        """
        name_out_file = self.solver_settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def AdvanceInTime(self):
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
        """This function checks the PythonSolver
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

    def Clear(self):
        """This function clears the PythonSolver
        """
        pass

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.solver_settings["computing_model_part_name"].GetString())


    def _GetRestartSettings(self):
        restart_settings = self.settings["model_import_settings"].Clone()
        restart_settings.RemoveValue("input_type")
        if not restart_settings.Has("restart_load_file_label"):
            raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')
        if self.solver_settings.Has("echo_level"):
            restart_settings.AddValue("echo_level", self.solver_settings["echo_level"])

        return restart_settings