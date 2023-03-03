from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class ResponseFunction(ABC):
    def Initialize(self) -> None:
        pass

    def InitializeSolutionStep(self) -> None:
        pass

    def FinalizeSolutionStep(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    @abstractmethod
    def Check(self) -> None:
        pass

    @abstractmethod
    def CalculateValue(self) -> float:
        pass

    @abstractmethod
    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart) -> None:
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass
