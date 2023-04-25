from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunction(ABC):
    @classmethod
    @abstractmethod

    def Initialize(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    @abstractmethod
    def ComputeResponseValue(self, MC) -> float:
        pass

    @abstractmethod
    def ComputeResponseGradients(self, gradients) -> None:
        pass