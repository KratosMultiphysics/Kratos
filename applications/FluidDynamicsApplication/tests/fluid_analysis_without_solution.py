from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidAnalysisWithoutSolution(FluidDynamicsAnalysis):
    """A fluid dynamics analysis variant that skips calls to solver.Predict() and solver.SolveSolutionStep().
       To be used in testing, it allows skipping the actual solution of the problem.
    """

    def __init__(self,model,parameters):
        super(FluidAnalysisWithoutSolution,self).__init__(model,parameters)

    def RunSolutionLoop(self):
        """ Run the solution loop, but do not actually solve anything.
        """

        while self.time < self.end_time:
            self.time = self._solver.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            #self.solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        # self._GetSolver().Predict()
        self._GetSolver().InitializeSolutionStep()