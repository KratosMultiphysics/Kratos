from abc import ABC
from abc import abstractmethod

from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class Modifier(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def ModifyControl(self, controls: ContainerData):
        pass

    @abstractmethod
    def ModifySensitivities(self, sensitivities: ContainerData):
        pass

    @abstractmethod
    def ModifyControlUpdates(self, control_updates: ContainerData):
        pass