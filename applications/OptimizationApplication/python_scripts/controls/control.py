# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class Control(ABC):
    def __init__(self, name: str, type: str, model: Kratos.Model, settings: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.name = name
        self.type = type
        self.model = model
        self.settings = settings
        self.optimization_info = optimization_info

    def Initialize(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def MapFirstDerivative(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]') -> None:
        pass

    @abstractmethod
    def Update(self,update_collective_expression: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        pass


