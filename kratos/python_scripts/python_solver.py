from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

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
        pass

    def AddDofs(self):
        pass

    def GetMinimumBufferSize(self):
        pass

    def ReadModelPart(self):
        # TODO replace the functions in the solvers in Fluid and Structure
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part.")
        problem_path = os.getcwd()
        input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
        if self.is_restarted():
            self.get_restart_utility().LoadRestart()
        elif(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]::", "Finished reading model part from mdpa file.")
        else:
            raise Exception("Other model part input options are not yet implemented.")
        KratosMultiphysics.Logger.PrintInfo("ModelPart", self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[PythonSolver]:: ", "Finished reading model part.")

    def PrepareModelPartForSolver(self):
        pass

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def AdvanceInTime(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def Initialize(self):
        pass

    def Predict(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def SolveSolutionStep(self):
        pass

    def Check(self):
        pass

    def Clear(self):
        pass

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())
