from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class ExecutionPolicy(ABC):
    def __init__(self, execution_policy_name: str) -> None:
        self.__name = execution_policy_name

    def GetName(self) -> str:
        return self.__name

    @abstractmethod
    def Initialize(self) -> None:
        pass

    @abstractmethod
    def Check(self) -> None:
        pass

    @abstractmethod
    def Finalize(self) -> None:
        pass

    @abstractmethod
    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        pass

    @abstractmethod
    def Execute(self) -> None:
        pass