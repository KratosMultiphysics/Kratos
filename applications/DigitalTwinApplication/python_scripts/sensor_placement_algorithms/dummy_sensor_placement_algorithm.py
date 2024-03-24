import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm

class DummySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        pass

    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{}""")

    def Execute(self, _: 'list[KratosDT.Sensors.Sensor]') -> None:
        pass
