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
class RedundancyGapSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "number_of_redundancy_groups"      : 15,
            "number_of_sensors_per_group"      : 2,
            "group_search_method"              : "most_distanced",
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.number_of_redundancy_groups = self.parameters["number_of_redundancy_groups"].GetInt()
        self.number_of_sensors_per_group = self.parameters["number_of_sensors_per_group"].GetInt()
        self.group_search_method = self.parameters["group_search_method"].GetString()

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
        if self.group_search_method == "most_distanced":
            distances = squareform(GetEuclideanDistances(normalized_sensor_views))
        elif self.group_search_method == "most_orthogonal":
            distances = squareform(GetCosineDistances(normalized_sensor_views))
        else:
            raise RuntimeError(f"Unuspported group_search_method = {self.group_search_method}.")
        best_coverage_sensor = GetBestCoverageSensorView(normalized_sensor_views)
        best_coverage_sensor_dist = best_coverage_sensor.GetContainerExpression().Evaluate()
        best_coverage_index = normalized_sensor_views.index(best_coverage_sensor)

        list_of_sensor_distributes: 'list[np.ndarray]' = [sensor_view.GetContainerExpression().Evaluate() for sensor_view in normalized_sensor_views]

        # prepare for output
        if self.is_vtu_output or self.is_csv_output:
            output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
            output_path.mkdir(parents=True, exist_ok=True)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(normalized_sensor_views[0].GetContainerExpression().GetModelPart())

        redundancy_gap = 1.0 / self.number_of_redundancy_groups
        list_of_available_sensor_indices = list(range(len(normalized_sensor_views)))
        list_of_sensor_view_indices: 'list[int]' = [best_coverage_index]
        del list_of_available_sensor_indices[best_coverage_index]
        for iteration, i in enumerate(reversed(range(self.number_of_redundancy_groups))):
            max_redundancy = (i+1) * redundancy_gap
            min_redundancy = i * redundancy_gap

            # get all the sensors between the redundancy gap
            list_of_sensor_views_in_redundancy_gap: 'list[int]' = []
            for sensor_dist_index in list_of_available_sensor_indices:
                redundancy = list_of_sensor_distributes[sensor_dist_index].dot(best_coverage_sensor_dist)
                if redundancy >= min_redundancy and redundancy < max_redundancy:
                    list_of_sensor_views_in_redundancy_gap.append(sensor_dist_index)

            # skip the rest of the code if not sensor_view found in the redundancy gap.
            if len(list_of_sensor_views_in_redundancy_gap) == 0:
                continue

            # now remove potential sensors in the gap from the available sensors
            # becaues they will not be potentital sensors for rest of the redundancy
            # gaps.
            for i in sorted(list_of_sensor_views_in_redundancy_gap, reverse=True):
                del list_of_available_sensor_indices[list_of_available_sensor_indices.index(i)]

            for _ in range(self.number_of_sensors_per_group):
                max_distance = 0.0
                for potential_sensor_i in list_of_sensor_views_in_redundancy_gap:
                    minimum_distance = 1e+9
                    for found_sensor_j in list_of_sensor_view_indices:
                        current_distance = distances[potential_sensor_i, found_sensor_j]
                        if minimum_distance > current_distance:
                            minimum_distance = current_distance
                    if minimum_distance > max_distance:
                        max_distance = minimum_distance
                        most_distanced_index = potential_sensor_i

                sensor_view = normalized_sensor_views[most_distanced_index]
                Kratos.Logger.PrintInfo("", f"---- Found {sensor_view.GetSensor().GetName()} for gap [{min_redundancy:0.3f}, {max_redundancy:0.3f}).")
                list_of_sensor_view_indices.append(most_distanced_index)
                del list_of_sensor_views_in_redundancy_gap[list_of_sensor_views_in_redundancy_gap.index(most_distanced_index)]

            PrintSensorListToCSV(output_path / f"overal_sensors_{iteration:05d}.csv", [normalized_sensor_views[sensor_view_index].GetSensor() for sensor_view_index in list_of_sensor_view_indices], ["type", "name", "location", "value"])
            PrintSensorListToCSV(output_path / f"sensors_in_redundancy_group_{iteration:05d}.csv", [normalized_sensor_views[sensor_view_index].GetSensor() for sensor_view_index in list_of_sensor_views_in_redundancy_gap], ["type", "name", "location", "value"])

        list_of_best_sensors = [normalized_sensor_views[sensor_view_index].GetSensor() for sensor_view_index in list_of_sensor_view_indices]
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_best_sensors)} sensors.")
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

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





