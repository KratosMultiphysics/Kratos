from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos

class OptimizationRoutine(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        self.model = model
        self.parameters = parameters

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    @abstractmethod
    def GetName(self) -> str:
        pass