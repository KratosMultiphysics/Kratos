from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class BaseKratosSolver(object):
    """The base class for the python solvers in the applications
    This class knows abt the solution of the physical problem
    """ 
    def __init__(self, main_model_part, project_parameters):
        """The constructor of the Solver.
        It obtains the project parameters
        This function is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(project_parameters)
        """
        if (type(main_model_part) != KratosMultiphysics.ModelPart):
            raise Exception("Input is expected to be provided as a Kratos ModelPart object")

        if (type(custom_settings) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.main_model_part = main_model_part
        self.project_parameters = project_parameters

        self._ValidateParameters()

    #### Public functions to run the Analysis ####
    def AddVariables(self):
        pass
    def AddDofs(self):
        pass
    def GetMinimumBufferSize(self):
        pass
    def PrepareModelPartForSolver(self):
        pass
    def _ValidateParameters(self):
        pass
    def ComputeDeltaTime(self):
        pass
    def Inizialize(self):
        pass
    def GetComputingModelPart(self):
        pass
    def SetEchoLevel(self):
        pass
    def Clear(self):
        pass
    def Check(self):
        pass

    def Solve(self):
        """This function solves one step
        It can be called several times during one time-step
        This is equivalent to calling "solving_strategy.Solve()" (without "Initialize")
        This function is NOT intended to be overridden in deriving classes!
        """
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be done
        (for each step) before solving the solution step.
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def Predict(self):
        """This function predicts the solution
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def SolveSolutionStep(self):
        """This function solves the current step
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def FinalizeSolutionStep(self):
        """This function Performs all the required operations that should be done
        (for each step) after solving the solution step.
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")