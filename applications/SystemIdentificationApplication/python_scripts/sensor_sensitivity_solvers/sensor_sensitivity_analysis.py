import abc

import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class SensorSensitivityAnalysis(AnalysisStage, abc.ABC):
    def GetProcessesOrder(self) -> 'list[str]':
        return ["auxiliary_processes", "output_processes"]

    @abc.abstractmethod
    def CalculateGradient(self, sensor: KratosSI.Sensors.Sensor) -> 'dict[str, ContainerExpressionTypes]':
        pass
