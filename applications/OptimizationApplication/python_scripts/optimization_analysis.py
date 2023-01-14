import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import RetrieveObject

from KratosMultiphysics.OptimizationApplication.mesh_controllers.mesh_controller import MeshController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import CreateResponseFunctionWrapper

class OptimizationAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)

        default_settings = Kratos.Parameters("""{
            "problem_data": {},
            "meshes"      : [],
            "analyses"    : [],
            "responses"   : [],
            "controls"    : [],
            "solvers"     : []
        }""")
        project_parameters.ValidateAndAssignDefaults(default_settings)
        self.__list_of_meshers = []

        self.optimization_info = OptimizationInfo()

        self._CreateMeshes()
        self._CreateAnalyses()

    def Initialize(self):
        super().Initialize()

        self._CreateResponses()
        self._CreateControls()

        # now set the buffers of the model parts
        root_model_part_buffer_sizes = {}
        for execution_policy_wrapper in self.optimization_info.GetExecutionPolicyWrappers():
            execution_policy_wrapper: ExecutionPolicyWrapper = execution_policy_wrapper
            analysis = execution_policy_wrapper.GetExecutionPolicy().GetAnalysis()
            if analysis is not None:
                analsis_solver = analysis._GetSolver()
                analysis_model_part: Kratos.ModelPart = analsis_solver.GetComputingModelPart()
                root_model_part_name = analysis_model_part.GetRootModelPart().FullName()

                if root_model_part_name not in root_model_part_buffer_sizes.keys():
                    root_model_part_buffer_sizes[root_model_part_name] = -1

                root_model_part_buffer_sizes[root_model_part_name] = max(root_model_part_buffer_sizes[root_model_part_name], analsis_solver.GetMinimumBufferSize())

        for mesh_controller in self.__list_of_meshers:
            mesh_model_part_name = mesh_controller.GetModelPart().FullName()
            if mesh_model_part_name in root_model_part_buffer_sizes.keys():
                mesh_controller.GetModelPart().SetBufferSize(root_model_part_buffer_sizes[mesh_model_part_name])
            mesh_controller.Initialize()

        self.optimization_info.Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.optimization_info.InitializeSolutionStep()

        for mesh_controller in self.__list_of_meshers:
            mesh_controller.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        for mesh_controller in self.__list_of_meshers:
            mesh_controller.FinalizeSolutionStep()

        self.optimization_info.FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

        for mesh_controller in self.__list_of_meshers:
            mesh_controller.Finalize()

        self.optimization_info.Finalize()

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
        pass

    def _CreateSolver(self):
        pass
        # return python_solvers_wrapper_fluid.CreateSolver(self.model, self.project_parameters)