from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class BaseKratosAnalysisStage(object):
    """The base class for the analysis classes in the applications
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        This function is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The ProjectParameters used
        """
        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(project_parameters) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object or a file name")

        self.model = model
        self.project_parameters = project_parameters

    #### Public functions to run the Analysis ####
    def Run(self):
        """This function executes the entire analysis
        It is NOT intended to be overridden in deriving classes! TODO change
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def RunSolutionLoop(self):
        """This function executes the temporal loop of the analysis
        It is NOT intended to be overridden in deriving classes!
        """
        while self.time < self.end_time:
            self.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self.Predict()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputStep()

    def Initialize(self):
        """This function initializes the analysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        At the end of this function the StageAnalysis is ready for the time-loop
        It should be called AFTER the ModelPart used for this AnalysisStage is used
        """
        pass

    def Finalize(self):
        """This function finalizes the analysis
        Usage: It is designed to be called ONCE, AFTER the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def AdvanceInTime(self, new_time):
        """This function prepares the database for the new time-step
        It is designed to be called once at the beginning of a time-step
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

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

    def Check(self):
        """This function Performs all the required operations that should be done
        (for each step) after solving the solution step.
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")
