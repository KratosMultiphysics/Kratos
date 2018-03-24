from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class BaseKratosAnalysis(object):
    """The base class for the analysis classes in the applications

    Deriving classes have to implement the _GetSolvingStrategy method
    """
    def __init__(self, project_parameters):
        """The constructor of the Analysis-Object.
        It obtains the project parameters used for the analysis
        This function is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)
        """
        if (type(project_parameters) == str): # a file name is provided
            with open(project_parameters,'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
        elif (type(project_parameters) == KratosMultiphysics.Parameters): # a Parameters object is provided
            self.ProjectParameters = project_parameters
        else:
            raise Exception("Input is expected to be provided as a Kratos Parameters object or a file name")

    #### Public functions to run the Analysis ####
    def Run(self):
        """This function executes the entire analysis
        It is NOT intended to be overridden in deriving classes!
        """
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def RunMainTemporalLoop(self):
        """This function executes the temporal loop of the analysis
        It is NOT intended to be overridden in deriving classes!
        """
        while self.time < self.end_time:
            self.InitializeTimeStep()
            self.SolveTimeStep()
            self.FinalizeTimeStep()

    def Initialize(self):
        """This function initializes the analysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def Finalize(self):
        """This function finalizes the analysis
        Usage: It is designed to be called ONCE, AFTER the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def InitializeTimeStep(self):
        """This function initializes the time-step
        Usage: It is designed to be called once at the beginning of EACH time-step
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def SolveStep(self):
        """This function solves one step
        It can be called several times during one time-step
        This is equivalent to calling "solving_strategy.Solve()" (without "Initialize")
        This function is NOT intended to be overridden in deriving classes!
        """
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def FinalizeTimeStep(self):
        """This function finalizes the time-step
        Usage: It is designed to be called once at the end of EACH time-step
        This function IS intended to be overridden in deriving classes!
        """
        pass

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