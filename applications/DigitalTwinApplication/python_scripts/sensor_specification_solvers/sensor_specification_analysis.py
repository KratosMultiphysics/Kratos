import abc

import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.analysis_stage import AnalysisStage

class SensorSpecificationAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def GetListOfSpecifications(self) -> 'list[KratosDT.Sensors.SensorSpecification]':
        pass