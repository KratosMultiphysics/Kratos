import importlib

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_analysis import SensorSensitivityAnalysis
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.cosine_similarity_sensor_placement_algorithm import CosineSimilaritySensorPlacementAlgorithm

class SensorPlacementAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"      : {},
            "model_parts"       : [],
            "test_analyses"     : [],
            "sensor_analyses"   : [],
            "algorithm_settings": {}
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__list_of_model_part_controllers: 'list[MdpaModelPartController]' = []

        self._CreateModelPartControllers()

        self.__list_of_test_analyses: 'list[AnalysisStage]' = self._CreateAnalyses(self.project_parameters["test_analyses"])
        self.__list_of_sensor_analyses: 'list[SensorSensitivityAnalysis]' = self._CreateAnalyses(self.project_parameters["sensor_analyses"])

        algorithm_type = self.project_parameters["algorithm_settings"]["type"].GetString()
        if algorithm_type == "cosine_similarity_sensor_placement":
            self.algorithm = CosineSimilaritySensorPlacementAlgorithm(self.model, self.project_parameters["algorithm_settings"])
        else:
            raise RuntimeError(f"Unsupported algorithm type = \"{algorithm_type}\" requested.")

    def Initialize(self):
        list(map(MdpaModelPartController.ImportModelPart, self.__list_of_model_part_controllers))
        list(map(MdpaModelPartController.Initialize, self.__list_of_model_part_controllers))

    def Run(self):
        self.Initialize()

        list(map(AnalysisStage.Run, self.__list_of_test_analyses))
        list(map(SensorSensitivityAnalysis.Run, self.__list_of_sensor_analyses))

        sensors_list: 'list[KratosDT.Sensors.SensorSpecification]' = []
        for sensor_analysis in self.__list_of_sensor_analyses:
            sensors_list.extend(sensor_analysis.GetListOfSensors())

        self.algorithm.Execute(sensors_list)

    def _CreateModelPartControllers(self):
        for model_part_controller_settings in self.project_parameters["model_parts"]:
            self.__list_of_model_part_controllers.append(MdpaModelPartController(self.model, model_part_controller_settings))

    def _CreateAnalyses(self, parameters: Kratos.Parameters) -> 'list[SensorSensitivityAnalysis]':
        list_of_analyses: 'list[SensorSensitivityAnalysis]' = []
        for analyses_settings in parameters:

            analysis_stage_module_name = analyses_settings["analysis_stage"].GetString()
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

            analysis_stage_module = importlib.import_module(analysis_stage_module_name)

            list_of_analyses.append(getattr(analysis_stage_module, analysis_stage_class_name)(self.model, analyses_settings))

        return list_of_analyses


