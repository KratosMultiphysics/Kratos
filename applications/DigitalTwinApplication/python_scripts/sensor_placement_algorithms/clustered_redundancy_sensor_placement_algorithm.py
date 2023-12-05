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
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetFilter, GetDistance, GetBestCoverageSensorView, GetSimilarSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionFilterUnionType, ExpressionUnionType

ClusterUnionType = typing.Union[
                            KratosDT.ClusterUtils.NodalSensorViewCluster,
                            KratosDT.ClusterUtils.ConditionSensorViewCluster,
                            KratosDT.ClusterUtils.ElementSensorViewCluster
                        ]
class ClusteredRedundancySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "number_of_base_sensors"           : 6,
            "redundancy_threshold"             : 0.7,
            "base_sensor_redundancy"           : 5,
            "redundancy_search_method"         : "most_distanced",
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.number_of_base_sensors = self.parameters["number_of_base_sensors"].GetInt()
        self.base_sensor_redundancy = self.parameters["base_sensor_redundancy"].GetInt()
        self.redundancy_threshold = self.parameters["redundancy_threshold"].GetDouble()
        self.redundancy_search_method = self.parameters["redundancy_search_method"].GetString()

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
        # prepare for output
        if self.is_vtu_output or self.is_csv_output:
            output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
            output_path.mkdir(parents=True, exist_ok=True)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(list_of_sensor_views[0].GetContainerExpression().GetModelPart())

        normalized_sensor_views = GetNormalizedSensorViews(list_of_sensor_views)

        if self.is_vtu_output:
            vtu_output.ClearCellContainerExpressions()
            vtu_output.ClearNodalContainerExpressions()

        # identify number of covered entities in each sensor view
        sensor_view_coverage_indices: 'list[np.ndarray]' = []
        for sensor_view in normalized_sensor_views:
            number_of_entities_covered, coverage_mask = KratosDT.SensorUtils.GetEntityCoverageMask(sensor_view)
            coverage = number_of_entities_covered / len(coverage_mask.GetContainer())
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Coverage of \"{sensor_view.GetSensor().GetName()}\" is {coverage * 100:6.3f}%")
            sensor_view.GetSensor().SetValue(KratosDT.SENSOR_COVERAGE, coverage)

            if self.is_vtu_output:
                vtu_output.AddContainerExpression(f"{sensor_view.GetSensor().GetName()}", coverage_mask)

        if self.is_csv_output:
            PrintSensorListToCSV(output_path / f"sensor_coverage.csv", [sensor_view.GetSensor() for sensor_view in normalized_sensor_views], ["type", "name", "location", "value", "SENSOR_COVERAGE"])

        if self.is_vtu_output:
            vtu_output.PrintOutput(str(output_path / "sensor_coverage"))

        raise RuntimeError(123)

        cosine_distances = GetCosineDistances(normalized_sensor_views)
        if self.redundancy_search_method == "most_distanced":
            distances = squareform(GetEuclideanDistances(normalized_sensor_views))
        elif self.redundancy_search_method == "most_orthogonal":
            distances = squareform(cosine_distances)
        else:
            raise RuntimeError(f"Unuspported redundancy_search_method = {self.redundancy_search_method}.")

        # now find the base sensors
        sensor_view_coverage_map = np.zeros((number_of_entities))
        list_of_selected_sensor_views: 'list[SensorViewUnionType]' = []
        for base_sensor_index in range(self.number_of_base_sensors):
            # find the best sensor which increases the coverage
            # by finding which maximizes the min_value
            min_value: float = min(sensor_view_coverage_map)
            current_coverage = (sensor_view_coverage_map > min_value).sum()
            best_coverage_index = -1
            for i, coverage_indices in enumerate(sensor_view_coverage_indices):
                current_coverage_map = np.array(sensor_view_coverage_map)
                current_coverage_map[coverage_indices] += 1
                new_coverage = (current_coverage_map > min_value).sum()
                if new_coverage > current_coverage:
                    current_coverage = new_coverage
                    best_coverage_index = i

            if best_coverage_index != -1:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found \"{normalized_sensor_views[best_coverage_index].GetSensor().GetName()}\" as a base sensor.")
                sensor_view_coverage_map[sensor_view_coverage_indices[best_coverage_index]] += 1
                # now find the redundant sensors for this sensor
                list_of_sensors_with_redundancy = self.GetRedundantSensors(self.base_sensor_redundancy, distances, normalized_sensor_views[best_coverage_index], normalized_sensor_views)
                list_of_selected_sensor_views.extend(list_of_sensors_with_redundancy)
                if self.is_csv_output:
                    PrintSensorListToCSV(output_path / f"base_sensor_{base_sensor_index:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_sensors_with_redundancy], ["type", "name", "location", "value"])
            else:
                break

        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list(set(list_of_selected_sensor_views))]
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_best_sensors)} sensors.")
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value", "SENSOR_COVERAGE"])

    def GetRedundantSensors(self, required_redundancy:int, distances: np.ndarray, base_sensor: SensorViewUnionType, normalized_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
        base_coverage_sensor_dist = base_sensor.GetContainerExpression().Evaluate()
        base_coverage_sensor_index = normalized_sensor_views.index(base_sensor)

        list_of_sensor_distributes: 'list[np.ndarray]' = [sensor_view.GetContainerExpression().Evaluate() for sensor_view in normalized_sensor_views]

        list_of_available_sensor_indices = list(range(len(normalized_sensor_views)))
        list_of_sensor_view_indices: 'list[int]' = [base_coverage_sensor_index]
        del list_of_available_sensor_indices[base_coverage_sensor_index]

        # get all the sensors between the redundancy gap
        list_of_sensor_views_in_redundancy_gap: 'list[int]' = []
        for sensor_dist_index in list_of_available_sensor_indices:
            redundancy = list_of_sensor_distributes[sensor_dist_index].dot(base_coverage_sensor_dist)
            if redundancy >= self.redundancy_threshold:
                list_of_sensor_views_in_redundancy_gap.append(sensor_dist_index)

        # skip the rest of the code if not sensor_view found in the redundancy gap.
        if len(list_of_sensor_views_in_redundancy_gap) == 0:
            return [base_sensor]

        for _ in range(required_redundancy):
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
            Kratos.Logger.PrintInfo("", f"---- Found {sensor_view.GetSensor().GetName()}.")
            list_of_sensor_view_indices.append(most_distanced_index)
            del list_of_sensor_views_in_redundancy_gap[list_of_sensor_views_in_redundancy_gap.index(most_distanced_index)]

        return [normalized_sensor_views[i] for i in list_of_sensor_view_indices]

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





