from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSmoothenedAbsoluteSensitivityFieldSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_utils import AddSensorViewMasks
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_utils import UpdateClusters
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_utils import OutputSensorViewExpressions
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_utils import PrintSensorViewsListToCSV


class RobustDistancedSensorSelection:
    def __init__(self, , params: Kratos.Parameters) -> None:
        default_params = Kratos.Parameters("""{
            "sensor_selection_method" : "robust_distanced",
            "distancing_coefficient"  : 0.5,
            "robustness_coefficient"  : 0.5,
            "sensor_group_id_varaible": "SENSOR_GROUP_ID",
            "robustness_settings"    : {
                "search_radius"           : 10.0,
                "max_number_of_neighbours": 1000
            }
        }""")
        params.RecursivelyValidateAndAssignDefaults(default_params)

        robustness_settings = params["robustness_settings"]
        self.search_radius = robustness_settings["search_radius"].GetDouble()
        self.max_number_of_neighbours = robustness_settings["max_number_of_neighbours"].GetInt()
        self.coeff_robustness = params["robustness_coefficient"].GetDouble()
        self.coeff_distancing = params["distancing_coefficient"].GetDouble()
        self.sensor_group_id_var = Kratos.KratosGlobals.GetVariable(params["sensor_group_id_varaible"].GetString())

        self.sensor_groups: 'dict[int, tuple[list[KratosDT.Sensors.Sensor], list[ExpressionUnionType]]]' = {}
        self.sensor_distance_matrices: 'dict[int, KratosDT.SensorDistanceMatrix]' = {}
        self.list_of_selected_sensor_views: 'list[SensorViewUnionType]' = []

    def Initialize(self, list_of_sensor_views: 'list[SensorViewUnionType]', list_of_selected_sensor_views: 'list[SensorViewUnionType]') -> None:
        self.list_of_selected_sensor_views = list_of_selected_sensor_views

        # first identify the different sensor groups / types
        for sensor_view in list_of_sensor_views:
            group_id = sensor_view.GetSensor().GetValue(self.sensor_group_id_var)
            if group_id not in self.sensor_groups.keys():
                self.sensor_groups[group_id] = ([], [])
            self.sensor_groups[group_id][0].append(sensor_view.GetSensor())
            self.sensor_groups[group_id][1].append(sensor_view.GetAuxiliaryExpression("filtered"))

        for group_id, (sensors_list, exp_list) in self.sensor_groups.items():
            # now compute robustness
            KratosDT.SensorUtils.AssignSensorNeighbours(sensors_list, self.search_radius, self.max_number_of_neighbours, KratosDT.SENSOR_NEIGHBOUR_IDS, KratosDT.SENSOR_ID)
            KratosDT.SensorUtils.ComputeSensorRobustness(sensors_list, exp_list, KratosDT.SENSOR_NEIGHBOUR_IDS, KratosDT.SENSOR_ID, KratosDT.SENSOR_LOCATION_ROBUSTNESS)

            Kratos.Logger.PrintInfo(self.__class__.__name__, "Sensor location based robustness values computed.")

            # now compute distance matrix
            self.sensor_distance_matrices[group_id] = KratosDT.SensorDistanceMatrix(sensors_list)
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Sensor distances computed.")

    def Sort(self, potential_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
        potential_sensors = [sensor_view.GetSensor() for sensor_view in potential_sensor_views]
        selected_sensors = [sensor_view.GetSensor() for sensor_view in self.list_of_selected_sensor_views]

        # compute the distances of potential and current sensor views
        distances = Kratos.Matrix(len(potential_sensor_views), len(self.list_of_selected_sensor_views), -1.0)
        for _, sensor_distance_matrix in self.sensor_distance_matrices.items():
            current_distance_matrix = sensor_distance_matrix.GetDistance(potential_sensors, selected_sensors)
            for i in range(distances.Size1()):
                for j in range(distances.Size2()):
                    distances[i, j] = max(distances[i, j], current_distance_matrix[i, j])

        # now compute the minimum distances for each potential sensor
        potential_sensor_minimum_distances: 'list[float]' = []
        for i in range(distances.Size1()):
            min_distance = 1e+100
            for j in range(distances.Size2()):
                if min_distance > distances[i, j] and distances[i, j] >= 0.0:
                    min_distance = distances[i, j]

            if min_distance == 1e+100:
                potential_sensor_minimum_distances.append(-1)
            else:
                potential_sensor_minimum_distances.append(min_distance)

        # now normalize the minimum distances
        max_minimum_distance = max(potential_sensor_minimum_distances)
        potential_sensor_minimum_distances = [d / max_minimum_distance for d in  potential_sensor_minimum_distances]

        # now get the robustness values
        potential_sensor_robustness = [sensor_view.GetSensor().GetValue(KratosDT.SENSOR_LOCATION_ROBUSTNESS) for sensor_view in potential_sensor_views]
        max_robustness = max(potential_sensor_robustness)
        potential_sensor_robustness = [d / max_robustness for d in  potential_sensor_robustness]

        # now compute the ranking
        ranking: 'list[tuple[SensorViewUnionType, float]]' = []
        for i, sensor_view in enumerate(potential_sensor_views):
            min_distance = potential_sensor_minimum_distances[i]
            robustness = potential_sensor_robustness[i]

            if min_distance > 0.0:
                ranking.append((sensor_view, self.coeff_distancing * min_distance + self.coeff_robustness * robustness))
            else:
                ranking.append((sensor_view, self.coeff_distancing + self.coeff_robustness * robustness))

        sorted_ranks = sorted(ranking, key=lambda x: x[1])
        return [sensor_view for sensor_view, _ in sorted_ranks]

class MaskBasedClusteringSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                                   : "mask_based_clustering_sensor_placement_algorithm",
            "output_folder"                          : "output",
            "required_maximum_overall_coverage_ratio": 0.95,
            "required_maximum_cluster_coverage_ratio": 0.01,
            "top_percentile"                         : 10,
            "maximum_number_of_iterations"           : 100,
            "filtering"                              : {},
            "sensor_selection_settings"              : {},
            "masking": {
                "masking_type": "automatic"
            },
            "vtu_output_settings": {
                "sensor_data"   : true,
                "iteration_data": true
            },
            "csv_output_settings": {
                "sensor_data"   : true,
                "iteration_data": true
            }
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.required_maximum_overall_coverage_ratio = self.parameters["required_maximum_overall_coverage_ratio"].GetDouble()
        self.required_maximum_cluster_coverage_ratio = self.parameters["required_maximum_cluster_coverage_ratio"].GetDouble()
        self.top_percentile = self.parameters["top_percentile"].GetInt()
        self.maximum_number_of_iterations = self.parameters["maximum_number_of_iterations"].GetInt()

        vtu_output_settings = self.parameters["vtu_output_settings"]
        vtu_output_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["vtu_output_settings"])
        self.vtu_output_sensor_data = vtu_output_settings["sensor_data"].GetBool()
        self.vtu_output_iteration_data = vtu_output_settings["iteration_data"].GetBool()

        csv_output_settings = self.parameters["csv_output_settings"]
        csv_output_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["csv_output_settings"])
        self.csv_output_sensor_data = csv_output_settings["sensor_data"].GetBool()
        self.csv_output_iteration_data = csv_output_settings["iteration_data"].GetBool()

        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.cluster_details: 'dict[int, tuple[list[int], ExpressionUnionType]]' = {}

        sensor_selection_method = self.parameters["sensor_selection_settings"]["sensor_selection_method"].GetString()
        if sensor_selection_method == "robust_distanced":
            self.sensor_selection_method = RobustDistancedSensorSelection(self.parameters["sensor_selection_settings"])
        else:
            raise RuntimeError("Unsupported sensor selection method")

    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        self.list_of_sensors = list_of_sensors

        KratosDT.SensorUtils.AssignConsecutiveSensorIds(self.list_of_sensors, KratosDT.SENSOR_ID)

        if len(self.list_of_sensors) == 0:
            raise RuntimeError("No sensors are given.")

        first_sensor = self.list_of_sensors[0]

        for k in first_sensor.GetNodalExpressionsMap().keys():
            list_of_sensor_views = GetSmoothenedAbsoluteSensitivityFieldSensorViews(self.list_of_sensors, k, Kratos.Globals.DataLocation.NodeHistorical, self.parameters["filtering"])
            self.ComputeSensorPlacement(f"nodal/{k}", list_of_sensor_views)

        for k in first_sensor.GetConditionExpressionsMap().keys():
            list_of_sensor_views = GetSmoothenedAbsoluteSensitivityFieldSensorViews(self.list_of_sensors, k, Kratos.Globals.DataLocation.Condition, self.parameters["filtering"])
            self.ComputeSensorPlacement(f"condition/{k}", list_of_sensor_views)

        for k in first_sensor.GetElementExpressionsMap().keys():
            list_of_sensor_views = GetSmoothenedAbsoluteSensitivityFieldSensorViews(self.list_of_sensors, k, Kratos.Globals.DataLocation.Element, self.parameters["filtering"])
            self.ComputeSensorPlacement(f"element/{k}", list_of_sensor_views)

    def ComputeSensorPlacement(self, name: str, list_of_sensor_views: 'list[SensorViewUnionType]'):
        dummy_exp = list_of_sensor_views[0].GetContainerExpression()
        self.cluster_details = {}

        if self.vtu_output_sensor_data or self.vtu_output_iteration_data:
            vtu_output = Kratos.VtuOutput(dummy_exp.GetModelPart())

        # add the masks lists for the sensor views
        AddSensorViewMasks(self.parameters["masking"], list_of_sensor_views)

        list_of_selected_sensor_views: 'list[SensorViewUnionType]' = []
        self.sensor_selection_method.Initialize(list_of_sensor_views, list_of_selected_sensor_views)

        # get the total domain size
        domain_size_exp = dummy_exp.Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(domain_size_exp)
        total_domain_size = KratosDT.SensorUtils.Sum(domain_size_exp)

        if self.vtu_output_sensor_data or self.csv_output_sensor_data:
            sensor_data_output_path = Path(self.parameters["output_folder"].GetString()) / f"{name}/sensor_data"
            sensor_data_output_path.mkdir(parents=True, exist_ok=True)

            if self.vtu_output_sensor_data:
                for sensor_view in list_of_sensor_views:
                    OutputSensorViewExpressions(sensor_data_output_path / f"{sensor_view.GetSensor().GetName()}", vtu_output, sensor_view)

            if self.csv_output_sensor_data:
                for sensor_view in list_of_sensor_views:
                    sensor_view.GetSensor().SetValue(KratosDT.SENSOR_COVERAGE, KratosDT.SensorUtils.Sum(domain_size_exp.Scale(sensor_view.GetAuxiliaryExpression("mask"))))
                PrintSensorViewsListToCSV(sensor_data_output_path / "sensor_mask_coverage.csv", list_of_sensor_views, ["type", "name", "location", "value", "SENSOR_COVERAGE"])

        # create iteration folder structure
        if self.vtu_output_iteration_data or self.csv_output_iteration_data:
            iteration_data_output_path = Path(self.parameters["output_folder"].GetString()) / f"{name}/iteration_data"
            iteration_data_output_path.mkdir(parents=True, exist_ok=True)

        # stage 1: Find the sensors which covers up to the required_maximum_overall_coverage_ratio
        iteration = 1
        max_coverage_ratio = 0.0
        sensor_redundancy = (dummy_exp * 0.0).Flatten()
        while (max_coverage_ratio < self.required_maximum_overall_coverage_ratio and iteration <= self.maximum_number_of_iterations):
            potential_sensor_views: 'list[tuple[SensorViewUnionType, float]]' = []
            for sensor_view in list_of_sensor_views:
                new_coverage_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, KratosDT.MaskUtils.Union(sensor_redundancy, sensor_view.GetAuxiliaryExpression("mask")))) / total_domain_size
                if new_coverage_ratio > max_coverage_ratio:
                    potential_sensor_views.append((sensor_view, new_coverage_ratio))

            potential_sensor_views = self.sensor_selection_method.GetSorted(potential_sensor_views, list_of_selected_sensor_views)
            list_of_selected_sensor_views.append(max(sensor_view))

            sensor_redundancy = (sensor_redundancy + sensor_view.GetAuxiliaryExpression("mask")).Flatten()
            cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks([sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_selected_sensor_views])
            clustering_exp = UpdateClusters(self.cluster_details, list_of_selected_sensor_views, cluster_data)

            if self.vtu_output_iteration_data:
                OutputSensorViewExpressions(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}", vtu_output, sensor_view, [("sensor_redundancy", sensor_redundancy), ("cluster", clustering_exp)])
            if self.csv_output_iteration_data:
                PrintSensorViewsListToCSV(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}.csv", list_of_selected_sensor_views, ["type", "name", "location", "value"])

            max_coverage_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, sensor_redundancy)) / total_domain_size
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Found sensor {sensor_view.GetSensor().GetName()} increasing coverage to {max_coverage_ratio * 100.0:6.3f}%.")

            del list_of_sensor_views[list_of_sensor_views.index(sensor_view)]
            iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished searching sensors to maximize coverage.")

        # stage 2: Find sensors which reduces the max cluster size
        list_of_cluster_ids_not_divisible: 'list[int]' = []
        while (iteration <= self.maximum_number_of_iterations):
            max_cluster_ids: 'list[int]' = []
            for cluster_id, (_, mask) in self.cluster_details.items():
                if cluster_id not in list_of_cluster_ids_not_divisible:
                    cluster_size_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) / total_domain_size
                    if cluster_size_ratio > self.required_maximum_cluster_coverage_ratio:
                        max_cluster_ids.append(cluster_id)

            if len(max_cluster_ids) == 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Finished dividing clusters.")
                break

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Found {len(max_cluster_ids)} clusters above threshold. Trying to divide them...")

            potential_masks = [sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views]
            potential_sensor_ranking: 'dict[int, int]' = {}
            for max_cluster_id in max_cluster_ids:
                potential_indices = KratosDT.MaskUtils.GetMasksDividingReferenceMask(self.cluster_details[max_cluster_id][1], potential_masks)
                for potential_index in potential_indices:
                    if potential_index not in potential_sensor_ranking.keys():
                        potential_sensor_ranking[potential_index] = 0
                    potential_sensor_ranking[potential_index] += 1

            potential_sensor_views: 'list[tuple[SensorViewUnionType, float]]' = [(list_of_sensor_views[k], v) for k, v in potential_sensor_ranking.items()]

            if len(potential_sensor_views) == 0:
                list_of_cluster_ids_not_divisible.extend(max_cluster_ids)
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tFound {len(max_cluster_ids)} cluster not divisible.")
            else:
                sensor_view = self.sensor_selection_method.Select(self.__GetTopPercentile(potential_sensor_views), list_of_selected_sensor_views)
                list_of_selected_sensor_views.append(sensor_view)

                sensor_redundancy = (sensor_redundancy + sensor_view.GetAuxiliaryExpression("mask")).Flatten()
                cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks([sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_selected_sensor_views])
                clustering_exp = UpdateClusters(self.cluster_details, list_of_selected_sensor_views, cluster_data)

                if self.vtu_output_iteration_data:
                    OutputSensorViewExpressions(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}", vtu_output, sensor_view, [("sensor_redundancy", sensor_redundancy), ("cluster", clustering_exp)])
                if self.csv_output_iteration_data:
                    PrintSensorViewsListToCSV(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}.csv", list_of_selected_sensor_views, ["type", "name", "location", "value"])

                Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tFound {sensor_view.GetSensor().GetName()} dividing cluster {max_cluster_id} resulting with {len(self.cluster_details.keys())} clusters.")
                iteration += 1

        # stage 3: remove some unnecessary sensors
        # first compute the current coverage and clustering.
        overall_coverage_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, sensor_redundancy)) / total_domain_size
        max_cluster_size_ratio = max([KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) for _, (_, mask) in self.cluster_details.items()]) / total_domain_size

        while (iteration <= self.maximum_number_of_iterations):
            # find sensors which does not affect current clustering or the overall coverage

            if len(list_of_selected_sensor_views) == 1:
                # at least one sensor is required.
                break

            potential_sensor_views: 'list[SensorViewUnionType]' = []
            for test_sensor_view in reversed(list_of_selected_sensor_views):
                masks_list: 'list[ExpressionUnionType]' = []
                for sensor_view in list_of_selected_sensor_views:
                    if sensor_view != test_sensor_view:
                        masks_list.append(sensor_view.GetAuxiliaryExpression("mask"))

                current_redundancy = sensor_redundancy * 0.0
                for mask in masks_list:
                    current_redundancy += mask

                current_coverage_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, current_redundancy)) / total_domain_size
                if current_coverage_ratio >= overall_coverage_ratio:
                    # the coverage is not reduced. Hence check for the clustering
                    cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks(masks_list)
                    current_max_cluster_size_ratio = max([KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) for _, mask in cluster_data]) / total_domain_size

                    if current_max_cluster_size_ratio <= max_cluster_size_ratio:
                        # the clustering is also ok. Then this sensor can be removed safely without degrading
                        # current sensor positioning
                        potential_sensor_views.append(test_sensor_view)

            if len(potential_sensor_views) == 0:
                break

            select_sensor_view = min(potential_sensor_views, key=lambda x: x.GetSensor().GetValue(KratosDT.SENSOR_LOCATION_ROBUSTNESS))
            del list_of_selected_sensor_views[list_of_selected_sensor_views.index(select_sensor_view)]

            sensor_redundancy = (sensor_redundancy - select_sensor_view.GetAuxiliaryExpression("mask")).Flatten()
            cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks([sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_selected_sensor_views])
            clustering_exp = UpdateClusters(self.cluster_details, list_of_selected_sensor_views, cluster_data)

            if self.vtu_output_iteration_data:
                OutputSensorViewExpressions(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}", vtu_output, select_sensor_view, [("sensor_redundancy", sensor_redundancy), ("cluster", clustering_exp)])
            if self.csv_output_iteration_data:
                PrintSensorViewsListToCSV(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}.csv", list_of_selected_sensor_views, ["type", "name", "location", "value"])

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Removed {select_sensor_view.GetSensor().GetName()}.")
            iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Summary:")

        total_coverage = domain_size_exp * 0.0
        for cluster_id, (_, mask) in self.cluster_details.items():
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"\t\tCluster id: {cluster_id} with coverage ratio = {KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) * 100.0 / total_domain_size:6.3f}%")
            total_coverage += mask

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"\t\tTotal coverage ratio = {KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, total_coverage)) * 100.0 / total_domain_size:6.3f}%")
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"\t\tFound {len(list_of_selected_sensor_views)} sensors.")

        output_path = Path(self.parameters["output_folder"].GetString())
        output_path.mkdir(parents=True, exist_ok=True)
        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_selected_sensor_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    def __GetTopPercentile(self, potential_sensor_view_values: 'list[tuple[SensorViewUnionType, float]]') -> 'list[SensorViewUnionType]':
        potential_sensor_view_values = sorted(potential_sensor_view_values, reverse=True, key=lambda x: x[1])
        upper_bound = min(len(potential_sensor_view_values), self.top_percentile)
        result: 'list[SensorViewUnionType]' = []
        for i, (sensor_view, value) in enumerate(potential_sensor_view_values):
            if i < upper_bound:
                latest_value = value
                result.append(sensor_view)
            else:
                if value == latest_value:
                    result.append(sensor_view)
                else:
                    break
        return result