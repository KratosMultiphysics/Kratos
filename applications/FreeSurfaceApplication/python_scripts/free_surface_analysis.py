# Kratos imports
import KratosMultiphysics
from KratosMultiphysics.FreeSurfaceApplication.edge_based_level_set_solver import EdgeBasedLevelSetSolver
from KratosMultiphysics.analysis_stage import AnalysisStage


class FreeSurfaceAnalysis(AnalysisStage):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        """"""
        # Base class constructor defines:
        # - self.model
        # - self.project_parameters
        # - self.echo_level
        # - self.parallel_type
        AnalysisStage.__init__(self, model, parameters)


    def _CreateSolver(self) -> EdgeBasedLevelSetSolver:
        return EdgeBasedLevelSetSolver(self.model, self.project_parameters["solver_settings"].Clone())


    def _GetSimulationName(self) -> str:
        return "Free Surface Analysis"