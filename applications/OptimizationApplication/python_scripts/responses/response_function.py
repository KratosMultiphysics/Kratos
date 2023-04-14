from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunction(ABC):
    @classmethod
    @abstractmethod
    def GetSensitivityFieldVariables(cls) -> 'list[SupportedSensitivityFieldVariableTypes]':
        pass

    def Initialize(self) -> None:
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
    def CalculateSensitivity(self, sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]') -> None:
        pass