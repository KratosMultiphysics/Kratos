from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion

class TransformationTechnique(OptimizationRoutine, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

    @abstractmethod
    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        pass

    @abstractmethod
    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        pass