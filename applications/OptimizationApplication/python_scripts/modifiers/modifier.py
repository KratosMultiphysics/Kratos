from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum

class Modifier(ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info

    @abstractmethod
    def ModifySensitivities(self, sensitivities: Kratos.Vector, model_part: Kratos.ModelPart, container_type: ContainerEnum):
        pass

    @abstractmethod
    def ModifyControlUpdates(self, update_change: Kratos.Vector, model_part: Kratos.ModelPart, container_type: ContainerEnum):
        pass