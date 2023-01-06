# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos

class ResponesFunction(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        self.model = model
        self.parameters = parameters

    @abstractmethod
    def Check(self, optimization_info: dict):
        pass

    @abstractmethod
    def Initialize(self, optimization_info: dict):
        pass

    @abstractmethod
    def InitializeIteration(self, optimization_info: dict):
        pass

    @abstractmethod
    def FinalizeIteration(self, optimization_info: dict):
        pass

    @abstractmethod
    def Finalize(self, optimization_info: dict):
        pass

    @abstractmethod
    def CalculateValue(self, optimization_info: dict) -> float:
        pass

    @abstractmethod
    def CalculateSensitivity(self, optimization_info: dict) -> dict:
        pass

    @abstractmethod
    def GetResponseFunctionName(self):
        pass