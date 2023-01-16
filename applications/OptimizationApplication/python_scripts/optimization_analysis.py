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
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.optimization_info = OptimizationInfo()
        super().__init__(model, project_parameters)
        self.optimization_info.SetBufferSize(self._GetSolver().GetMinimumBufferSize())

        default_settings = Kratos.Parameters("""{
            "problem_data" : {},
            "meshes"       : [],
            "analyses"     : [],
            "responses"    : [],
            "controls"     : [],
            "optimizations": {}
        }""")
        project_parameters.AddMissingParameters(default_settings)
        self.__list_of_meshers = []

        self._CreateMeshes()
        self._CreateAnalyses()

    def Initialize(self):
        # first we initialize mesh controlelrs
        for mesh_controller in self.__list_of_meshers:
            mesh_controller.Initialize()

        # we execute secondly the execution policy wrappers because they involve
        # solving problems using analsis (or not).
        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "Initialize")

        # now we create other types
        self._CreateResponses()
        self._CreateControls()

        # now we initialize rest
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "Initialize")
        self.__ExecuteRoutinesMethod("ControlWrapper", "Initialize")

        # then we initialize the optimization solver
        super().Initialize()

    def InitializeSolutionStep(self):
        for mesh_controller in self.__list_of_meshers:
            mesh_controller.InitializeSolutionStep()

        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "InitializeSolutionStep")

        super().InitializeSolutionStep()

        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "InitializeSolutionStep")
        self.__ExecuteRoutinesMethod("ControlWrapper", "InitializeSolutionStep")

    def FinalizeSolutionStep(self):
        for mesh_controller in self.__list_of_meshers:
            mesh_controller.FinalizeSolutionStep()

        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "FinalizeSolutionStep")
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "FinalizeSolutionStep")
        self.__ExecuteRoutinesMethod("ControlWrapper", "FinalizeSolutionStep")
        super().FinalizeSolutionStep()

    def Finalize(self):
        for mesh_controller in self.__list_of_meshers:
            mesh_controller.Finalize()

        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "Finalize")

        super().Finalize()

        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "Finalize")
        self.__ExecuteRoutinesMethod("ControlWrapper", "Finalize")

    def _CreateMeshes(self):
        for mesher_settings in self.project_parameters["meshes"]:
            self.__list_of_meshers.append(RetrieveObject(self.model, mesher_settings, self.optimization_info, MeshController))

    def _CreateAnalyses(self):
        for analyses_settings in self.project_parameters["analyses"]:
            self.optimization_info.AddRoutine(ExecutionPolicyWrapper(self.model, analyses_settings))

    def _CreateResponses(self):
        for response_settings in self.project_parameters["responses"]:
            self.optimization_info.AddRoutine(CreateResponseFunctionWrapper(self.model, response_settings, self.optimization_info))

    def _CreateControls(self):
        for control_settings in self.project_parameters["controls"]:
            self.optimization_info.AddRoutine(ControlWrapper(self.model, control_settings, self.optimization_info))

    def _CreateSolver(self):
        return OptimizationSolver(self.model, self.project_parameters["optimizations"], self.optimization_info)

    def __ExecuteRoutinesMethod(self, routine_class_type_name: str, execution_method: str):
        if self.optimization_info.HasRoutineType(routine_class_type_name):
            for routine in self.optimization_info.GetRoutines(routine_class_type_name):
                getattr(routine, execution_method)()
        else:
            Kratos.Logger.PrintWarning(self.__class__.__name__, f"No routine type \"{routine_class_type_name}\" is found in optimization info. Hence, skipping running \"{execution_method}\" methods on those types.")

    def _GetSimulationName(self):
        return "OptimizationAnalysis"

    def PrintAnalysisStageProgressInformation(self):
        Kratos.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self.optimization_info["step"])

    def KeepAdvancingSolutionLoop(self):
        return self.optimization_info["step"] < self.end_time and not self._GetSolver().IsConverged()