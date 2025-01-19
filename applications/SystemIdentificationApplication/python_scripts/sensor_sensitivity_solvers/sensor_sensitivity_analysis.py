import abc

import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.analysis_stage import AnalysisStage

class SensorSensitivityAnalysis(AnalysisStage, abc.ABC):
    @abc.abstractmethod
    def GetListOfSensors(self) -> 'list[KratosSI.Sensors.Sensor]':
        """Returns the list of sensors used in the sensor sensitivity analysis.

        Returns:
            list[KratosSI.Sensors.Sensor]: List of sensors used in the sensor sensitivity analysis.
        """
        pass