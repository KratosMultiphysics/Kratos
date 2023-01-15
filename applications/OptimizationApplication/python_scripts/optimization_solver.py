import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.python_solver import PythonSolver

from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import RetrieveObject
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm

class OptimizationSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "echo_level"     : 0,
            "model_part_name": "optimization_data_holder_model_part",
            "algorithms"     : []
        }""")

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, settings)
        self.optimization_info = optimization_info

        default_solver_settings = Kratos.Parameters("""{
            "module"  : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"    : "PLEASE_PROVIDE_SOLVER_CLASS_NAME",
            "settings": {}
        }""")

        self.__list_of_algorithms = []
        for algorithm_settings in settings["algorithms"]:
            algorithm_settings.ValidateAndAssignDefaults(default_solver_settings)
            self.__list_of_algorithms.append(RetrieveObject(self.model, algorithm_settings, self.optimization_info, Algorithm))

        if len(self.__list_of_algorithms) == 0:
            raise RuntimeError("No optimization solvers provided.")

        # create the global data holding model part
        self.model_part = self.model.CreateModelPart(self.settings["model_part_name"].GetString())

    def GetMinimumBufferSize(self):
        buffer_size = -1e+30
        for algorithm in self.__list_of_algorithms:
            if algorithm.GetMinimumBufferSize() > buffer_size:
                buffer_size = algorithm.GetMinimumBufferSize()
        return buffer_size

    def AddVariables(self):
        super().AddVariables()
        self.__ExecuteMethod("AddVariables")

    def AddDofs(self):
        super().AddDofs()
        self.__ExecuteMethod("AddDofs")

    def ImportModelPart(self):
        pass

    def Initialize(self):
        self.optimization_info["step"] = 0
        super().Initialize()
        self.__ExecuteMethod("Initialize")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.__ExecuteMethod("InitializeSolutionStep")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.__ExecuteMethod("FinalizeSolutionStep")

    def Finalize(self):
        super().Finalize()
        self.__ExecuteMethod("Finalize")

    def GetComputingModelPart(self):
        return self.model_part

    def AdvanceInTime(self, _):
        self.optimization_info.AdvanceSolutionStep()
        self.optimization_info["step"] = self.optimization_info["step", 1] + 1
        return self.optimization_info["step"]

    def __ExecuteMethod(self, method_name: str):
        for algorithm in self.__list_of_algorithms:
            getattr(algorithm, method_name)()