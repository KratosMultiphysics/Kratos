import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.python_solver import PythonSolver

from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import RetrieveObject

from KratosMultiphysics.OptimizationApplication.mesh_controllers.mesh_controller import MeshController
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import CreateResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper

class OptimizationSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "model_part_name": "optimization_model_part",
            "meshes"         : [],
            "analyses"       : [],
            "responses"      : [],
            "controls"       : [],
            "algorithms"     : [],
            "echo_level"     : 0
        }""")

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        super().__init__(model, settings)

        # creates the optimization info data holder
        self.optimization_info = OptimizationInfo()

        default_solver_settings = Kratos.Parameters("""{
            "module"  : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"    : "PLEASE_PROVIDE_SOLVER_CLASS_NAME",
            "settings": {}
        }""")

        self.__list_of_meshers = []

        self._CreateMeshes()
        self._CreateAnalyses()

        self.__list_of_algorithms = []
        for algorithm_settings in settings["algorithms"]:
            algorithm_settings.ValidateAndAssignDefaults(default_solver_settings)
            self.__list_of_algorithms.append(RetrieveObject(self.model, algorithm_settings, self.optimization_info, Algorithm))

        if len(self.__list_of_algorithms) == 0:
            raise RuntimeError("No optimization solvers provided.")

        # set the optimization info buffer size
        self.optimization_info.SetBufferSize(self.GetMinimumBufferSize())

    def GetOptimizationInfo(self):
        return self.optimization_info

    def GetMinimumBufferSize(self):
        buffer_size = -1e+30
        for algorithm in self.__list_of_algorithms:
            if algorithm.GetMinimumBufferSize() > buffer_size:
                buffer_size = algorithm.GetMinimumBufferSize()
        return buffer_size

    def AddVariables(self):
        self.__ExecuteMethod(self.__list_of_algorithms, "AddVariables")

    def AddDofs(self):
        self.__ExecuteMethod(self.__list_of_algorithms, "AddDofs")

    def ImportModelPart(self):
        self.__ExecuteMethod(self.__list_of_meshers, "ImportModelPart")

        # now we create other types because, responses and
        # controls can be used on submodel parts which are
        # only available after importing the whole model part
        self._CreateResponses()
        self._CreateControls()

    def Initialize(self):
        # set the current step to 0
        self.optimization_info["step"] = 0

        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "Initialize")
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "Initialize")
        self.__ExecuteRoutinesMethod("ControlWrapper", "Initialize")
        self.__ExecuteMethod(self.__list_of_algorithms, "Initialize")

    def InitializeSolutionStep(self):
        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "InitializeSolutionStep")
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "InitializeSolutionStep")
        self.__ExecuteRoutinesMethod("ControlWrapper", "InitializeSolutionStep")
        self.__ExecuteMethod(self.__list_of_algorithms, "InitializeSolutionStep")

    def SolveSolutionStep(self):
        self.__ExecuteMethod(self.__list_of_algorithms, "SolveSolutionStep")
        return False

    def IsConverged(self):
        is_converged = True

        for algorithm in self.__list_of_algorithms:
            is_converged = is_converged and algorithm.IsConverged()

        return is_converged

    def FinalizeSolutionStep(self):
        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "FinalizeSolutionStep")
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "FinalizeSolutionStep")
        self.__ExecuteRoutinesMethod("ControlWrapper", "FinalizeSolutionStep")
        self.__ExecuteMethod(self.__list_of_algorithms, "FinalizeSolutionStep")

    def Finalize(self):
        self.__ExecuteRoutinesMethod("ExecutionPolicyWrapper", "Finalize")
        self.__ExecuteRoutinesMethod("ResponseFunctionWrapper", "Finalize")
        self.__ExecuteRoutinesMethod("ControlWrapper", "Finalize")
        self.__ExecuteMethod(self.__list_of_algorithms, "Finalize")

    def GetComputingModelPart(self):
        model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(model_part_name):
            self.model.CreateModelPart(model_part_name)

        return self.model[model_part_name]

    def AdvanceInTime(self, _):
        self.optimization_info.AdvanceSolutionStep()
        self.optimization_info["step"] = self.optimization_info.GetSolutionStepData(1)["step"] + 1
        return self.optimization_info["step"]

    def __ExecuteMethod(self, items, method_name: str):
        for itr in items:
            getattr(itr, method_name)()

    def __ExecuteRoutinesMethod(self, routine_class_type_name: str, execution_method: str):
        if self.optimization_info.HasOptimizationRoutineType(routine_class_type_name):
            for routine in self.optimization_info.GetOptimizationRoutines(routine_class_type_name):
                getattr(routine, execution_method)()
        else:
            Kratos.Logger.PrintWarning(self.__class__.__name__, f"No routine type \"{routine_class_type_name}\" is found in optimization info. Hence, skipping running \"{execution_method}\" methods on those types.")

    def _CreateMeshes(self):
        for mesher_settings in self.settings["meshes"]:
            self.__list_of_meshers.append(RetrieveObject(self.model, mesher_settings, self.optimization_info, MeshController))

    def _CreateAnalyses(self):
        for analyses_settings in self.settings["analyses"]:
            self.optimization_info.AddOptimizationRoutine(ExecutionPolicyWrapper(self.model, analyses_settings))

    def _CreateResponses(self):
        for response_settings in self.settings["responses"]:
            self.optimization_info.AddOptimizationRoutine(CreateResponseFunctionWrapper(self.model, response_settings, self.optimization_info))

    def _CreateControls(self):
        for control_settings in self.settings["controls"]:
            self.optimization_info.AddOptimizationRoutine(ControlWrapper(self.model, control_settings, self.optimization_info))