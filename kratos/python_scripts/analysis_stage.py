from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class AnalysisStage(object):
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



        ##HERE WE SHOULD CONSTRUCT A SOLVER - stages should contain at least one solver
        ##self._GetSolver() = ... HERE WE CONSTRUCT THE SOLVER

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
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()


    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        self._GetSolver().ImportModelPart()
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()


    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()


    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        for process in self._GetListOfProcesses():
            process.ExecuteBeforeOutputStep()

        #here the output should be done when needed

        for process in self._GetListOfProcesses():
            process.ExecuteAfterOutputStep()


    def Check(self):
        """This function checks the AnalysisStage
        """
        pass

    def ModifyInitialProperties(self):
        """this is the place to eventually modify material properties in the stage """
        pass

    def ModifyInitialGeometry(self):
        """this is the place to eventually modify geometry (for example moving nodes) in the stage """
        pass

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        #other operations as needed

    def ChangeMaterialProperties(self):
        """this function is where the user could change material parameters as a part of the solution step """
        pass

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    def _GetListOfProcesses(self):
        if not hasattr(self, '_list_of_processes'):
            self._list_of_processes = self._CreateListOfProcesses()
        return self._list_of_processes

    def _CreateListOfProcesses(self):
        """Create the list of Processes
        """
        raise Exception("Creation of the processes must be implemented in the derived class.")
