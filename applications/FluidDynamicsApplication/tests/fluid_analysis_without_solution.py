from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class FluidAnalysisWithoutSolution(FluidDynamicsAnalysis):
    """A fluid dynamics analyis variant that skips calls to solver.Predict() and solver.SolveSolutionStep().
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
            #self.solver.Predict()
            #self.solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()