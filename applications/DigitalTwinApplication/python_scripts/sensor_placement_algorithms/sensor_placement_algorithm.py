import abc
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

class SensorPlacementAlgorithm(abc.ABC):
    @classmethod
    @abc.abstractmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        pass

    @abc.abstractmethod
    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        pass
