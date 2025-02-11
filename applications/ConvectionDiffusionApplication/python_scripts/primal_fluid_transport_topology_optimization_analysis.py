# Importing the Kratos Library
import KratosMultiphysics


# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.ConvectionDiffusionApplication import fluid_transport_topology_optimization_solver

class PrimalFluidTransportTopologyOptimizationAnalysis(ConvectionDiffusionAnalysis):
    def __init__(self,model,parameters):
        super().__init__( model,parameters )
        self._DefineTransportCustomization()

    def _CreateSolver(self):
        return fluid_transport_topology_optimization_solver.CreateSolver(self.model, self.project_parameters)
    
    def _DefineTransportCustomization(self, decay=0.0):
        # self._GetSolver()._GetTransportSolver()._DefineNodalProperties(decay)
        dumb = 0

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        self._SetTopologyOptimizationStage(problem_stage=1)
        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()
            # self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def _SetTopologyOptimizationStage(self, problem_stage):
        self.topology_optimization_stage = problem_stage
        self._GetSolver()._SetTopologyOptimizationStage(self.topology_optimization_stage)
    
