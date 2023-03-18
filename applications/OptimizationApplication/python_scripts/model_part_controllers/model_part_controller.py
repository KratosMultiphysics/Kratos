from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class ModelPartController(ABC):
    def Initialize(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    @abstractmethod
    def ImportModelPart(self):
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass


