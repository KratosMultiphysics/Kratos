# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

from enum import Enum
import KratosMultiphysics as Kratos

class ContainerEnum(Enum):
    NODES = 1
    ELEMENTS = 2
    CONDITIONS = 3
    ELEMENT_PROPERTIES = 4
    CONDITION_PROPERTIES = 5

class ResponseSentivity:
    def __init__(self, variable, sensitivity_container_type: ContainerEnum, model_part: Kratos.ModelPart):
        self.__sensitivity_container_type = sensitivity_container_type
        self.__variable = variable
        self.__model__part = model_part

    def GetSensitivityContainerType(self):
        return  self.__sensitivity_container_type

    def GetSensitivityVariable(self):
        return self.__variable

    def GetValues(self):
        return Kratos.Matrix()



