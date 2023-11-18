import importlib
import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.OptimizationApplication.model_part_controllers.mdpa_model_part_controller import MdpaModelPartController
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_analysis import SensorSensitivityAnalysis
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.cosine_similarity_sensor_placement_algorithm import CosineSimilaritySensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionDataLocation
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
class SensorPlacementAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"      : {},
            "model_parts"       : [],
            "test_analyses"     : [],
            "sensor_analyses"   : [],
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

        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []

        restart_settings = self.project_parameters["restart_settings"]
        sensitivity_state = restart_settings["sensor_sensitivity_state"].GetString()

        if sensitivity_state in ["save", "none"]:
            list(map(AnalysisStage.Run, self.__list_of_test_analyses))
            list(map(SensorSensitivityAnalysis.Run, self.__list_of_sensor_analyses))
            for sensor_analysis in self.__list_of_sensor_analyses:
                sensors_list.extend(sensor_analysis.GetListOfSensors())

            # if save is prescribed then save it
            if sensitivity_state == "save":
                self.__Save(restart_settings["file_name"].GetString(), sensors_list)
        elif sensitivity_state == "load":
            sensors_list = self.__Load(restart_settings["file_name"].GetString())
        else:
            raise RuntimeError(f"Unsupported sensor_sensitivity_state = {sensitivity_state}.")

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

    def __Load(self, file_name: str) -> 'list[KratosDT.Sensors.Sensor]':
        h5_file_params = Kratos.Parameters("""{
            "file_name": "",
            "file_access_mode": "read_only"
        }""")
        h5_file_params["file_name"].SetString(file_name)

        sensors_list: 'list[KratosDT.Sensors.Sensor]' = []
        for sensor_analysis in self.__list_of_sensor_analyses:
            sensor_analysis.Initialize()
            current_sensors_list = sensor_analysis.GetListOfSensors()
            sensitivity_model_part: Kratos.ModelPart = sensor_analysis._GetSolver().GetSensitivityModelPart()
            sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = sensor_analysis._GetSolver().GetSensitivtyVariables()
            with OpenHDF5File(h5_file_params, sensitivity_model_part) as h5_file:
                for sensor in current_sensors_list:
                    h5_path = f"/{sensitivity_model_part.FullName()}/{sensor.GetName()}"
                    sensor.SetSensorValue(h5_file.ReadAttribute(h5_path)["value"].GetDouble())

                    exp_params = Kratos.Parameters("""{
                        "prefix": ""
                    }""")
                    for data_location, vars in sensitivity_variables.items():
                        if data_location in [ExpressionDataLocation.NodeHistorical, ExpressionDataLocation.NodeNonHistorical]:
                            exp_params["prefix"].SetString(f"{h5_path}/nodal/")
                            expio = KratosHDF5.ExpressionIO(exp_params, h5_file)
                            exp = Kratos.Expression.NodalExpression(sensitivity_model_part)
                            for var in vars:
                                expio.Read(var.Name(), exp)
                                sensor.AddNodalExpression(var.Name(), exp.Clone())
                        elif data_location in [ExpressionDataLocation.Condition, ExpressionDataLocation.ConditionProperties]:
                            exp_params["prefix"].SetString(f"{h5_path}/condition/")
                            expio = KratosHDF5.ExpressionIO(exp_params, h5_file)
                            exp = Kratos.Expression.ConditionExpression(sensitivity_model_part)
                            for var in vars:
                                expio.Read(var.Name(), exp)
                                sensor.AddConditionExpression(var.Name(), exp.Clone())
                        elif data_location in [ExpressionDataLocation.Element, ExpressionDataLocation.ElementProperties]:
                            exp_params["prefix"].SetString(f"{h5_path}/element/")
                            expio = KratosHDF5.ExpressionIO(exp_params, h5_file)
                            exp = Kratos.Expression.ElementExpression(sensitivity_model_part)
                            for var in vars:
                                expio.Read(var.Name(), exp)
                                sensor.AddElementExpression(var.Name(), exp.Clone())
            sensors_list.extend(current_sensors_list)
            sensor_analysis.Finalize()

        return sensors_list

    def __Save(self, file_name: str, sensors_list: 'list[KratosDT.Sensors.Sensor]'):
        h5_file_params = Kratos.Parameters("""{
            "file_name": "",
            "file_access_mode": "read_write"
        }""")
        h5_file_params["file_name"].SetString(file_name)

        # first get list of model parts
        set_of_model_parts: 'set[Kratos.ModelPart]' = set()
        for sensor in sensors_list:
            for exp in sensor.GetNodalExpressionsMap().values():
                set_of_model_parts.add(exp.GetModelPart())
            for exp in sensor.GetConditionExpressionsMap().values():
                set_of_model_parts.add(exp.GetModelPart())
            for exp in sensor.GetElementExpressionsMap().values():
                set_of_model_parts.add(exp.GetModelPart())

        for model_part in set_of_model_parts:
            with OpenHDF5File(h5_file_params, model_part) as h5_file:
                exp_params = Kratos.Parameters("""{
                    "prefix": ""
                }""")
                for sensor in sensors_list:
                    h5_path = f"/{model_part.FullName()}/{sensor.GetName()}"

                    # create the path recursively.
                    __temp = h5_path.split("/")
                    for i, _ in enumerate(__temp[1:]):
                        __temp_full_name = "/".join(__temp[0:i+2])
                        if not h5_file.HasPath(__temp_full_name):
                            h5_file.CreateGroup(__temp_full_name)

                    # write sensor value
                    sensor_value_parameters = Kratos.Parameters("""{
                        "value" : 0.0
                    }""")
                    sensor_value_parameters["value"].SetDouble(sensor.GetSensorValue())
                    h5_file.WriteAttribute(h5_path, sensor_value_parameters)

                    # write nodal sensitivities
                    exp_params["prefix"].SetString(f"{h5_path}/nodal/")
                    exp_io = KratosHDF5.ExpressionIO(exp_params, h5_file)
                    for name, exp in sensor.GetNodalExpressionsMap().items():
                        if exp.GetModelPart() == model_part:
                            exp_io.Write(name, exp)

                    # write condition sensitivities
                    exp_params["prefix"].SetString(f"{h5_path}/condition/")
                    exp_io = KratosHDF5.ExpressionIO(exp_params, h5_file)
                    for name, exp in sensor.GetConditionExpressionsMap().items():
                        if exp.GetModelPart() == model_part:
                            exp_io.Write(name, exp)

                    # write element sensitivities
                    exp_params["prefix"].SetString(f"{h5_path}/element/")
                    exp_io = KratosHDF5.ExpressionIO(exp_params, h5_file)
                    for name, exp in sensor.GetElementExpressionsMap().items():
                        if exp.GetModelPart() == model_part:
                            exp_io.Write(name, exp)

