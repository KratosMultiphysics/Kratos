from abc import ABC
from abc import abstractmethod

from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion

class TransformationTechnique(ABC):
    def __init__(self):
        pass

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    @abstractmethod
    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        pass

    @abstractmethod
    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        pass