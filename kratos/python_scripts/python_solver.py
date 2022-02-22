# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.restart_utility import RestartUtility

# Other imports
import os

class PythonSolver:
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super().__init__(settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        settings -- The solver settings used
        """
        if not isinstance(model, KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if not isinstance(settings, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.settings = settings

        self.ValidateSettings()

        self.echo_level = self.settings["echo_level"].GetInt()

    @classmethod
    def GetDefaultParameters(cls):
        """This function returns the default-settings used by this class
        """
        return KratosMultiphysics.Parameters("""{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "use_input_model_part"
            }
        }""")

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        default_settings = self.GetDefaultParameters()
        self.settings.ValidateAndAssignDefaults(default_settings)

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

    def GetDofsList(self):
        """This function returns a list containing the DOFs required by this PythonSolver
        In needs to be reimplemented in the derived solvers according to the physics to be solved
        Ideally, DOFs list must be retrieved from the elements and conditions GetSpecification method
        (see GetDofsListFromSpecifications in the SpecificationUtilities)
        """
        raise Exception("This function has to be implemented in the derived class")

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
        Returns whether the problem is converged
        """
        return True

    def Check(self):
        """This function checks the PythonSolver. It usually calls the "Check" function of a solving strategy
        """
        pass

    def Clear(self):
        """This function clears the PythonSolver
        """
        pass

    def GetComputingModelPart(self):
        raise Exception("This function has to be implemented in the derived class")

    def _ImportModelPart(self, model_part, model_part_import_settings):
        """This function imports the ModelPart
        """
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part.")
        input_type = model_part_import_settings["input_type"].GetString()

        if input_type == "mdpa": # NOTE: Add more types in the future
            KratosMultiphysics.SingleImportModelPart.Import(model_part, model_part_import_settings, input_type)
        elif input_type == "rest":
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Loading model part from restart file.")
            RestartUtility(model_part, self._GetRestartSettings(model_part_import_settings)).LoadRestart()
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished loading model part from restart file.")

        elif input_type == "use_input_model_part":
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Using already imported model part - no reading necessary.")
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
