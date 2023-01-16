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
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class MeshController:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info

    @abstractmethod
    def ImportModelPart(self):
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass


