from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos

class ModelPartController(ABC):
    def __init__(self):
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
    def ImportModelPart(self):
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass


