from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedControlVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityVariableTypes

class ResponseFunction(ABC):
    def Initialize(self) -> None:
        pass

    def InitializeSolutionStep(self) -> None:
        pass

    def FinalizeSolutionStep(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetRequiredSensitivityVariablesListForControlVariable(self, control_variable: SupportedControlVariableTypes) -> 'list[SupportedSensitivityVariableTypes]':
        return [Kratos.KratosGlobals.GetVariable(f"{control_variable.Name()}_SENSITIVITY")]

    @abstractmethod
    def Check(self) -> None:
        pass

    @abstractmethod
    def CalculateValue(self) -> float:
        pass

    @abstractmethod
    def CalculateSensitivity(self, sensitivity_variable: SupportedSensitivityVariableTypes, sensitivity_model_part: Kratos.ModelPart) -> None:
        pass