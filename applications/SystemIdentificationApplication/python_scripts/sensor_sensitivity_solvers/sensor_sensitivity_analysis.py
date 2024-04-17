import abc

import KratosMultiphysics.SystemIdentificationApplication as KratosDT
from KratosMultiphysics.analysis_stage import AnalysisStage

class SensorSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def GetListOfSensors(self) -> 'list[KratosDT.Sensors.Sensor]':
        pass