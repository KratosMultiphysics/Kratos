from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class Algorithm(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info

    def AddVariables(self):
        pass

    def AddDofs(self):
        pass

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    @abstractmethod
    def GetMinimumBufferSize(self) -> int:
        pass

    @abstractmethod
    def SolveSolutionStep(self):
        pass

    @abstractmethod
    def IsConverged(self) -> bool:
        pass