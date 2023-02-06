from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine

class ModelPartController(OptimizationRoutine, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

    @abstractmethod
    def ImportModelPart(self):
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass


