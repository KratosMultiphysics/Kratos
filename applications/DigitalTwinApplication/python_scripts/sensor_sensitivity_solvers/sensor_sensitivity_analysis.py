import abc

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.analysis_stage import AnalysisStage

class SensorSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def GetListOfSensors(self) -> 'list[KratosDT.Sensors.Sensor]':
        pass

    @abc.abstractmethod
    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        pass