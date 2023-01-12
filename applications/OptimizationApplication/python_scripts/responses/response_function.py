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
from enum import Enum

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo

class ContainerEnum(Enum):
    NODES = 1
    ELEMENTS = 2
    CONDITIONS = 3
    ELEMENT_PROPERTIES = 4
    CONDITION_PROPERTIES = 5

def GetSensitivityContainer(model_part: Kratos.ModelPart, container_type: ContainerEnum):
    match (container_type):
        case ContainerEnum.NODES:
            return model_part.Nodes
        case ContainerEnum.ELEMENTS | ContainerEnum.ELEMENT_PROPERTIES:
            return model_part.Elements
        case ContainerEnum.CONDITIONS | ContainerEnum.CONDITION_PROPERTIES:
            return model_part.Conditions
        case _:
            raise RuntimeError("Unsupported container type requested.")

class ResponseFunction(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    @abstractmethod
    def Check(self):
        pass

    @abstractmethod
    def CalculateValue(self) -> float:
        pass

    @abstractmethod
    def CalculateSensitivity(self, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum) -> None:
        pass

