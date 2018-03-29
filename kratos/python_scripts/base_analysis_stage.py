from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class BaseKratosAnalysisStage(object):
    """The base class for the AnalysisStage-classes in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        It is intended to be called from the constructor
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
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.time < self.end_time:
            self.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self.Predict()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived classes")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived classes")

    def AdvanceInTime(self, new_time):
        """This function prepares the database for the new time-step
        It is designed to be called once at the beginning of a time-step
        It typically calles the CloneTimeStep function of a ModelPart
        """
        pass

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        pass

    def Predict(self):
        """This function predicts the solution and should be executed
        (for each step) BEFORE solving the solution step.
        """
        pass

    def SolveSolutionStep(self):
        """This function solves the current step
        It can be called multiple times during one time-step
        """
        raise NotImplementedError("This function has to be implemented by derived classes")

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        pass

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        pass

    def Check(self):
        """This function checks the AnalysisStage
        """
        pass
