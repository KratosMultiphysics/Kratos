from KratosMultiphysics.analysis_stage_with_solver import AnalysisStageWithSolver
from KratosMultiphysics.CompressiblePotentialFlowApplication import python_solvers_wrapper_compressible_potential

class PotentialFlowAnalysis(AnalysisStageWithSolver):
    '''Main script for potential flow simulations.'''

    def _CreateSolver(self):
        return python_solvers_wrapper_compressible_potential.CreateSolver(self.model, self.project_parameters)
