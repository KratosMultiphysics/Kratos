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

    def Initialize(self):
        pass

    def InitializeIteration(self):
        pass

    def FinalizeIteration(self):
        pass

    def Finalize(self):
        pass

    @abstractmethod
    def Execute(self):
        pass

    @abstractmethod
    def GetAnalysis(self):
        pass