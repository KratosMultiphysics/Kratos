from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunction(ABC):
    @abstractmethod
    def ComputeResponseValue(des_variables, info )-> float:
        pass
    
    @abstractmethod
    def ComputeResponseGradients(des_variables, info)-> 'List[float]':
        pass