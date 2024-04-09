import importlib
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_analysis import SensorSensitivityAnalysis
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_io import OpenSensorFile
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class SensorPlacementAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"      : {},
            "model_parts"       : [],
            "test_analyses"     : [],
            "sensor_analyses"   : [],
            "vtu_output_path"   : "",
            "filter_settings"   : {
                "YOUNG_MODULUS_SENSITIVITY": {}
            },
            "algorithm_settings": {},
            "restart_settings"  : {
                "sensor_sensitivity_state": "none",
                "file_name"               : "sensor_sensitivities.h5"
            }
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.project_parameters["restart_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["restart_settings"])
        self.optimization_problem = OptimizationProblem(0)

        self.__list_of_model_part_controllers: 'list[MdpaModelPartController]' = []

        self._CreateModelPartControllers()

        self.__list_of_test_analyses: 'list[AnalysisStage]' = self._CreateAnalyses(self.project_parameters["test_analyses"])
        self.__list_of_sensor_analyses: 'list[SensorSensitivityAnalysis]' = self._CreateAnalyses(self.project_parameters["sensor_analyses"])

        algorithm_type = self.project_parameters["algorithm_settings"]["type"].GetString()
        algorithm_module = importlib.import_module(f"KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.{algorithm_type}")
        self.algorithm: SensorPlacementAlgorithm = getattr(algorithm_module, Kratos.StringUtilities.ConvertSnakeCaseToCamelCase(algorithm_type))(self.model, self.project_parameters["algorithm_settings"])

    def Initialize(self):
        list(map(MdpaModelPartController.ImportModelPart, self.__list_of_model_part_controllers))
        list(map(MdpaModelPartController.Initialize, self.__list_of_model_part_controllers))

    def Run(self):
        self.Initialize()

        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []

        restart_settings = self.project_parameters["restart_settings"]
        sensitivity_state = restart_settings["sensor_sensitivity_state"].GetString()

        self.filters: 'dict[tuple[str, str], Filter]' = {}

        if sensitivity_state in ["save", "none"]:
            list(map(AnalysisStage.Run, self.__list_of_test_analyses))
            list(map(SensorSensitivityAnalysis.Run, self.__list_of_sensor_analyses))
            for sensor_analysis in self.__list_of_sensor_analyses:
                list_of_sensors = sensor_analysis.GetListOfSensors()

                for sensor in list_of_sensors:
                    for var_name, _ in sensor.GetNodalExpressionsMap().items():
                        self.__PostProcessSensorViewData(KratosDT.Sensors.NodalSensorView(sensor, var_name), Kratos.Globals.DataLocation.NodeHistorical)
                    for var_name, _ in sensor.GetConditionExpressionsMap().items():
                        self.__PostProcessSensorViewData(KratosDT.Sensors.ConditionSensorView(sensor, var_name), Kratos.Globals.DataLocation.Condition)
                    for var_name, _ in sensor.GetElementExpressionsMap().items():
                        self.__PostProcessSensorViewData(KratosDT.Sensors.ElementSensorView(sensor, var_name), Kratos.Globals.DataLocation.Element)

                sensitivity_mp = sensor_analysis.GetSensitivityModelPart()
                if sensitivity_state == "save":
                    with OpenSensorFile(sensitivity_mp, self.project_parameters["restart_settings"]["file_name"].GetString(), f"/{sensitivity_mp.FullName()}", "a") as sensor_io:
                        for sensor in list_of_sensors:
                            sensor_io.Write(sensor)

                vtu_output_path = self.project_parameters["vtu_output_path"].GetString()
                if vtu_output_path != "":
                    Path(vtu_output_path).mkdir(exist_ok=True, parents=True)
                    vtu_output = Kratos.VtuOutput(sensitivity_mp)
                    for sensor in list_of_sensors:
                        vtu_output.ClearCellContainerExpressions()
                        vtu_output.ClearNodalContainerExpressions()
                        for k, v in sensor.GetNodalExpressionsMap().items():
                            vtu_output.AddContainerExpression(k, v)
                        for k, v in sensor.GetConditionExpressionsMap().items():
                            vtu_output.AddContainerExpression(k, v)
                        for k, v in sensor.GetElementExpressionsMap().items():
                            vtu_output.AddContainerExpression(k, v)
                        vtu_output.PrintOutput(f"{vtu_output_path}/{sensor.GetName()}")

                sensors_list.extend(list_of_sensors)
        elif sensitivity_state == "load":
            sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
            for sensor_analysis in self.__list_of_sensor_analyses:
                sensor_analysis.Initialize()
                sensitivity_mp: Kratos.ModelPart = sensor_analysis.GetSensitivityModelPart()
                current_sensors_list = sensor_analysis.GetListOfSensors()

                with OpenSensorFile(sensitivity_mp, self.project_parameters["restart_settings"]["file_name"].GetString(), f"/{sensitivity_mp.FullName()}", "r") as sensor_io:
                    for sensor in current_sensors_list:
                        sensor_io.Read(sensor)
                        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Read {sensor.GetName()} data.")
                sensors_list.extend(current_sensors_list)
        else:
            raise RuntimeError(f"Unsupported sensor_sensitivity_state = {sensitivity_state}.")

        KratosDT.SensorUtils.AssignSensorIds(sensors_list)
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

    def __PostProcessSensorViewData(self, sensor_view: SensorViewUnionType, data_location: Kratos.Globals.DataLocation) -> None:
        exp = sensor_view.GetContainerExpression()
        var_name = sensor_view.GetExpressionName()
        model_part_name = exp.GetModelPart().FullName()

        if not (var_name, model_part_name) in self.filters.keys():
            parameters = self.project_parameters["filter_settings"][var_name]
            new_filter = FilterFactory(self.model, exp.GetModelPart().FullName(), Kratos.KratosGlobals.GetVariable(var_name), data_location, parameters)
            new_filter.SetComponentDataView(ComponentDataView(var_name, self.optimization_problem))
            new_filter.Initialize()
            new_filter.Check()
            new_filter.Update()
            self.filters[(var_name, model_part_name)] = new_filter

        current_filter = self.filters[(var_name, model_part_name)]
        abs_exp = Kratos.Expression.Utils.Abs(exp)
        sensor_view.AddAuxiliaryExpression("abs", abs_exp)
        filtered_exp = current_filter.ForwardFilterField(current_filter.BackwardFilterIntegratedField(abs_exp))
        sensor_view.AddAuxiliaryExpression("filtered", filtered_exp)
        l2_norm = Kratos.Expression.Utils.NormL2(filtered_exp)
        if l2_norm > 0.0:
            normalized_filtered_exp = filtered_exp / l2_norm
        else:
            normalized_filtered_exp = filtered_exp.Clone()
        mask_exp = KratosDT.MaskUtils.GetMask(normalized_filtered_exp)
        sensor_view.AddAuxiliaryExpression("mask", mask_exp)
