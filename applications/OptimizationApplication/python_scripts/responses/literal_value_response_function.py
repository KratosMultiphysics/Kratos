import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class LiteralValueResponseFunction(ResponseFunction):
    def __init__(self, value: float):
        super().__init__(str(value))

        self.value = value

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return []

    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        raise RuntimeError(f"The literal value response function does not have an influencing model part.")

    def CalculateValue(self) -> float:
        return self.value

    def CalculateGradient(self, _: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        raise RuntimeError(f"The literal value response function does not depend on any variable, hence no gradients.")

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}]"