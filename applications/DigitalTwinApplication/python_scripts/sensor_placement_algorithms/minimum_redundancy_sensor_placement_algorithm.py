import typing
import abc
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from pathlib import Path
import pyomo.environ as pyomo
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances, GetEuclideanDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetFilter, GetDistance, GetBestCoverageSensorView
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionFilterUnionType, ExpressionUnionType

ClusterUnionType = typing.Union[
                            KratosDT.ClusterUtils.NodalSensorViewCluster,
                            KratosDT.ClusterUtils.ConditionSensorViewCluster,
                            KratosDT.ClusterUtils.ElementSensorViewCluster
                        ]
class MinimumRedundancySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "required_minimum_redundancy"      : 4,
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.required_minimum_redundancy = self.parameters["required_minimum_redundancy"].GetInt()

    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        self.list_of_sensors = list_of_sensors

        if len(self.list_of_sensors) == 0:
            raise RuntimeError("No specifications are given.")

        self.SmoothenSensitivityFields()

        first_specification = self.list_of_sensors[0]

        for k in first_specification.GetNodalExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.NodalSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.NodalSensorViewCluster
            self.ComputeSensorPlacement(f"nodal/{k}", [KratosDT.Sensors.NodalSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetConditionExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.ConditionSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.ConditionSensorViewCluster
            self.ComputeSensorPlacement(f"condition/{k}", [KratosDT.Sensors.ConditionSensorView(sensor, k) for sensor in list_of_sensors])

        for k in first_specification.GetElementExpressionsMap().keys():
            self.domain_sensor_view_cluster_type = KratosDT.ClusterUtils.ElementSensorViewClusterData
            self.cluster_type = KratosDT.ClusterUtils.ElementSensorViewCluster
            self.ComputeSensorPlacement(f"element/{k}", [KratosDT.Sensors.ElementSensorView(sensor, k) for sensor in list_of_sensors])

    def ComputeSensorPlacement(self, name: str, list_of_sensor_views: 'list[SensorViewUnionType]'):
        normalized_sensor_views = GetNormalizedSensorViews(list_of_sensor_views)
        distances = squareform(GetEuclideanDistances(normalized_sensor_views))
        best_coverage_sensor = GetBestCoverageSensorView(normalized_sensor_views)
        best_coverage_index = normalized_sensor_views.index(best_coverage_sensor)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(normalized_sensor_views[0].GetContainerExpression().GetModelPart())

        threshold = 0.5
        self.iteration = 0
        while (True):
            self.iteration += 1
            resulting_sensor_views = self.ComputeRedundancyForThreshold(threshold, normalized_sensor_views, distances, best_coverage_index,vtu_output)
            if len(resulting_sensor_views) == 0:
                threshold /= 2.0
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Reducing threshold to {threshold}")
            else:
                break

        # prepare for output
        if self.is_vtu_output or self.is_csv_output:
            output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
            output_path.mkdir(parents=True, exist_ok=True)

        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in resulting_sensor_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    def ComputeRedundancyForThreshold(self, threshold: float, normalized_sensor_views: 'list[SensorViewUnionType]', distances: np.ndarray, initial_sensor_index: int, vtu_output: Kratos.VtuOutput) -> 'list[SensorViewUnionType]':
        list_of_threshold_values: 'list[np.ndarray]' = []
        for sensor_view in normalized_sensor_views:
            list_of_threshold_values.append(np.array(np.where(sensor_view.GetContainerExpression().Evaluate() < threshold, 0, 1), dtype=np.int32))
            cexp = Kratos.Expression.ElementExpression(vtu_output.GetModelPart())
            Kratos.Expression.CArrayExpressionIO.Read(cexp, list_of_threshold_values[-1])
            vtu_output.AddContainerExpression(f"{sensor_view.GetSensor().GetName()}", cexp)
        vtu_output.PrintOutput(f"sensor_placement/sensors/threshold_{self.iteration:05d}")

        list_of_found_sensor_view_indices: 'list[int]' = [initial_sensor_index]
        entity_redundancy = np.array(list_of_threshold_values[initial_sensor_index])
        minimum_redundancy = np.min(entity_redundancy)
        remaining_indices = list(range(len(normalized_sensor_views)))
        del remaining_indices[remaining_indices.index(initial_sensor_index)]
        local_iteration = 0
        while (minimum_redundancy < self.required_minimum_redundancy):
            local_iteration += 1
            # find the indices corresponds to current minimum redundancy
            minimum_redundant_mask = np.array(np.where(entity_redundancy == minimum_redundancy, 1, 0), dtype=np.int32)

            vtu_output.ClearCellContainerExpressions()
            cexp = Kratos.Expression.ElementExpression(vtu_output.GetModelPart())
            Kratos.Expression.CArrayExpressionIO.Read(cexp, minimum_redundant_mask)
            vtu_output.AddContainerExpression("minimum_redundant_mask", cexp.Clone())

            sensor_indices_coverage_dict: 'dict[int, list[int]]' = {}
            for i in remaining_indices:
                threshold_value = list_of_threshold_values[i]
                number_of_activations = np.dot(threshold_value, minimum_redundant_mask)
                if number_of_activations not in sensor_indices_coverage_dict.keys():
                    sensor_indices_coverage_dict[number_of_activations] = []
                sensor_indices_coverage_dict[number_of_activations].append(i)

            max_number_of_activations = max(sensor_indices_coverage_dict.keys())
            if max_number_of_activations == 0:
                return []

            # now sort the max_number_of_activations group based max min cosine distance
            max_distance = 0
            max_min_potential_sensor_index = -1
            # print(sensor_indices_coverage_dict[max_number_of_activations])
            for potential_sensor_index in sensor_indices_coverage_dict[max_number_of_activations]:
                minimum_distance = 1e+9
                for found_sensor_index in list_of_found_sensor_view_indices:
                    current_distance = distances[found_sensor_index, potential_sensor_index]
                    if current_distance < minimum_distance:
                        minimum_distance = current_distance
                if minimum_distance > max_distance:
                    max_distance = minimum_distance
                    max_min_potential_sensor_index = potential_sensor_index

            # now add the potential sensor index
            list_of_found_sensor_view_indices.append(max_min_potential_sensor_index)
            del remaining_indices[remaining_indices.index(max_min_potential_sensor_index)]

            Kratos.Logger.PrintInfo("", f"--- Found sensor {normalized_sensor_views[max_min_potential_sensor_index].GetSensor().GetName()}.")

            # update the minimum redundancy
            entity_redundancy += list_of_threshold_values[max_min_potential_sensor_index]
            minimum_redundancy = min(entity_redundancy)

            cexp = Kratos.Expression.ElementExpression(vtu_output.GetModelPart())
            Kratos.Expression.CArrayExpressionIO.Read(cexp, entity_redundancy)
            vtu_output.AddContainerExpression("entity_redundancy", cexp.Clone())
            for i in list_of_found_sensor_view_indices:
                Kratos.Expression.CArrayExpressionIO.Read(cexp, list_of_threshold_values[i])
                vtu_output.AddContainerExpression(f"{normalized_sensor_views[i].GetSensor().GetName()}", cexp.Clone())
            vtu_output.PrintOutput(f"sensor_placement/sensors/output_{self.iteration:05d}_{local_iteration:05d}")

        return [normalized_sensor_views[i] for i in list_of_found_sensor_view_indices]

    def SmoothenSensitivityFields(self) -> None:
        mp_nodal, nodal_filter = self.__GetFilter(Kratos.Globals.DataLocation.NodeHistorical)
        mp_condition, condition_filter = self.__GetFilter(Kratos.Globals.DataLocation.Condition)
        mp_element, element_filter = self.__GetFilter(Kratos.Globals.DataLocation.Element)

        output_path = Path(self.parameters["output_folder"].GetString())

        model_part = mp_nodal
        if model_part is None:
            model_part = mp_condition
        if model_part is None:
            model_part = mp_element

        if mp_nodal not in [None, model_part] or mp_condition not in [None, model_part] or mp_element not in [None, model_part]:
            raise RuntimeError(f"Found mismatching model parts.")

        if self.is_vtu_output:
            vtu_output_path = output_path / "sensitivities"
            vtu_output_path.mkdir(parents=True, exist_ok=True)
            vtu_output = Kratos.VtuOutput(model_part)

        for sensor in self.list_of_sensors:
            if self.is_vtu_output:
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()
            for k, v in sensor.GetNodalExpressionsMap().items():
                filtered_field = nodal_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in sensor.GetConditionExpressionsMap().items():
                filtered_field = condition_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            for k, v in sensor.GetElementExpressionsMap().items():
                filtered_field = element_filter.FilterIntegratedField(v.Abs())
                if self.is_vtu_output:
                    vtu_output.AddContainerExpression(f"{k}_raw", v.Clone())
                    vtu_output.AddContainerExpression(f"{k}_filtered", filtered_field.Clone())

                v.SetExpression(filtered_field.GetExpression())

            if self.is_vtu_output:
                vtu_output.PrintOutput(str(vtu_output_path / f"{sensor.GetName()}"))

    def __GetFilter(self, filter_field_type: Kratos.Globals.DataLocation) -> 'tuple[Kratos.ModelPart, ExpressionFilterUnionType]':
        field_filter: ExpressionFilterUnionType = None
        model_part: Kratos.ModelPart = None
        for sensor in self.list_of_sensors:
            if filter_field_type in [Kratos.Globals.DataLocation.NodeHistorical, Kratos.Globals.DataLocation.NodeNonHistorical]:
                expressions_dict = sensor.GetNodalExpressionsMap()
            elif filter_field_type == Kratos.Globals.DataLocation.Condition:
                expressions_dict = sensor.GetConditionExpressionsMap()
            elif filter_field_type == Kratos.Globals.DataLocation.Element:
                expressions_dict = sensor.GetElementExpressionsMap()

            for v in expressions_dict.values():
                if field_filter is None:
                    model_part = v.GetModelPart()
                    field_filter = GetFilter(model_part, filter_field_type, self.parameters["filtering"])
                    break

            if field_filter is not None:
                break

        return model_part, field_filter





