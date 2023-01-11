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

class ExecutionPolicy(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        self.model = model
        self.parameters = parameters

    def Initialize(self, optimization_info: dict):
        pass

    def InitializeIteration(self, optimization_info: dict):
        pass

    def FinalizeIteration(self, optimization_info: dict):
        pass

    def Finalize(self, optimization_info: dict):
        pass

    @abstractmethod
    def Execute(self, optimization_info: dict):
        pass

    @abstractmethod
    def GetAnalysis(self, optimization_info: dict):
        pass