from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.CompressiblePotentialFlowApplication.python_solvers_wrapper_compressible_potential import CreateSolverByParameters

class PotentialFlowAnalysis(AnalysisStage):
    '''Main script for potential flow simulations.'''

    def _CreateSolver(self):
        return CreateSolverByParameters(self.model, self.project_parameters["solver_settings"], self.parallel_type)