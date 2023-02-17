import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ResponseFunctionImplementor

class Control(Kratos.Process):
    def __init__(self):
        super().__init__()
        self.__name = None

    def SetName(self, name: str):
        self.__name = name

    def GetName(self) -> str:
        return self.__name

    def CalculateSensitivity(self, response_function: ResponseFunctionImplementor, output_sensitivities: KratosOA.CollectiveVariableDataHolder):
        raise NotImplementedError("Calling base class Control::CalculateSensitivity method. Please implement it in the derrived class.")

    def UpdateControl(self, update: KratosOA.CollectiveVariableDataHolder):
        raise NotImplementedError("Calling base class Control::UpdateControl method. Please implement it in the derrived class.")

