import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import RetrieveObject

from KratosMultiphysics.OptimizationApplication.mesh_controllers.mesh_controller import MeshController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import CreateResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.optimization_solver import OptimizationSolver

class OptimizationAnalysis(AnalysisStage):
    def _CreateSolver(self):
        return OptimizationSolver(self.model, self.project_parameters["optimization_settings"])

    def _GetSimulationName(self):
        return "OptimizationAnalysis"

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetOptimizationInfo()["step"])

    def KeepAdvancingSolutionLoop(self):
        return self._GetSolver().GetOptimizationInfo()["step"] < self.end_time and not self._GetSolver().IsConverged()