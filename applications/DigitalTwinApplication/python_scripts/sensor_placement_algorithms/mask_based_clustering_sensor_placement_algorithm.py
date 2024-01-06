from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSmoothenedAbsoluteSensitivityFieldSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType

class RobustDistancedSensorSelection:
    def __init__(self, params: Kratos.Parameters) -> None:
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

    def Initialize(self, list_of_sensor_views: 'list[SensorViewUnionType]') -> None:
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

    def Select(self, potential_sensor_views: 'list[SensorViewUnionType]', selected_sensor_views: 'list[SensorViewUnionType]') -> SensorViewUnionType:
        potential_sensors = [sensor_view.GetSensor() for sensor_view in potential_sensor_views]
        selected_sensors = [sensor_view.GetSensor() for sensor_view in selected_sensor_views]

        # compute the distances of potential and current sensor views
        distances = Kratos.Matrix(len(potential_sensor_views), len(selected_sensor_views), -1.0)
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

        sensor_view, _  = max(ranking, key=lambda x: x[1])
        return sensor_view

class MaskBasedClusteringSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                                   : "mask_based_clustering_sensor_placement_algorithm",
            "output_folder"                          : "output",
            "required_maximum_overall_coverage_ratio": 0.9,
            "required_maximum_cluster_coverage_ratio": 0.1,
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
        self.__AddSensorViewMasks(list_of_sensor_views)

        self.sensor_selection_method.Initialize(list_of_sensor_views)

        # get the total domain size
        domain_size_exp = dummy_exp.Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(domain_size_exp)
        total_domain_size = KratosDT.SensorUtils.Sum(domain_size_exp)

        if self.vtu_output_sensor_data or self.csv_output_sensor_data:
            sensor_data_output_path = Path(self.parameters["output_folder"].GetString()) / f"{name}/sensor_data"
            sensor_data_output_path.mkdir(parents=True, exist_ok=True)

            if self.vtu_output_sensor_data:
                for sensor_view in list_of_sensor_views:
                    self.__OutputSensorViewExpressions(sensor_data_output_path / f"{sensor_view.GetSensor().GetName()}", vtu_output, sensor_view)

            if self.csv_output_sensor_data:
                for sensor_view in list_of_sensor_views:
                    sensor_view.GetSensor().SetValue(KratosDT.SENSOR_COVERAGE, KratosDT.SensorUtils.Sum(domain_size_exp.Scale(sensor_view.GetAuxiliaryExpression("mask"))))
                self.__PrintSensorViewsListToCSV(sensor_data_output_path / "sensor_mask_coverage.csv", list_of_sensor_views, ["type", "name", "location", "value", "SENSOR_COVERAGE"])

        # create iteration folder structure
        if self.vtu_output_iteration_data or self.csv_output_iteration_data:
            iteration_data_output_path = Path(self.parameters["output_folder"].GetString()) / f"{name}/iteration_data"
            iteration_data_output_path.mkdir(parents=True, exist_ok=True)

        # stage 1: Find the sensors which covers up to the required_maximum_overall_coverage_ratio
        iteration = 1
        max_coverage_ratio = 0.0
        sensor_redundancy = (dummy_exp * 0.0).Flatten()
        list_of_selected_sensor_views: 'list[SensorViewUnionType]' = []
        while (max_coverage_ratio < self.required_maximum_overall_coverage_ratio and iteration <= self.maximum_number_of_iterations):
            sensor_view: SensorViewUnionType = max(list_of_sensor_views, key=lambda x: KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, sensor_redundancy + x.GetAuxiliaryExpression("mask"))))
            list_of_selected_sensor_views.append(sensor_view)

            sensor_redundancy = (sensor_redundancy + sensor_view.GetAuxiliaryExpression("mask")).Flatten()
            self.__OutputIteration(iteration, iteration_data_output_path, vtu_output, sensor_redundancy, list_of_selected_sensor_views)

            max_coverage_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, sensor_redundancy)) / total_domain_size
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Found sensor {sensor_view.GetSensor().GetName()} increasing coverage to {max_coverage_ratio * 100.0:6.3f}%.")

            del list_of_sensor_views[list_of_sensor_views.index(sensor_view)]
            iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished searching sensors to maximize coverage.")

        # stage 2: Find sensors which reduces the max cluster size
        list_of_cluster_ids_not_divisible: 'list[int]' = []
        while (iteration <= self.maximum_number_of_iterations):
            max_cluster_id = -1
            max_cluster_size_ratio = 0.0
            for cluster_id, (_, mask) in self.cluster_details.items():
                if cluster_id not in list_of_cluster_ids_not_divisible:
                    cluster_size_ratio = KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) / total_domain_size
                    if cluster_size_ratio > max_cluster_size_ratio:
                        max_cluster_size_ratio = cluster_size_ratio
                        max_cluster_id = cluster_id

            if max_cluster_size_ratio < self.required_maximum_cluster_coverage_ratio:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Finished dividing clusters.")
                break

            if max_cluster_id == -1:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Clusters are not divisible further.")
                break

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Cluser {max_cluster_id} size ratio is above threshold ({max_cluster_size_ratio * 100.0:6.3f}% > {self.required_maximum_cluster_coverage_ratio * 100.0:6.3f}%). Trying to divide it...")

            potential_masks = [sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views]
            potential_indices = KratosDT.MaskUtils.GetMasksDividingReferenceMask(self.cluster_details[max_cluster_id][1], potential_masks)

            if len(potential_indices) == 0:
                list_of_cluster_ids_not_divisible.append(max_cluster_id)
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tCluster {max_cluster_id} is not divisible.")
            else:
                sensor_view = self.sensor_selection_method.Select([list_of_sensor_views[i] for i in potential_indices], list_of_selected_sensor_views)
                list_of_selected_sensor_views.append(sensor_view)

                sensor_redundancy = (sensor_redundancy + sensor_view.GetAuxiliaryExpression("mask")).Flatten()
                self.__OutputIteration(iteration, iteration_data_output_path, vtu_output, sensor_redundancy, list_of_selected_sensor_views)

                Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tFound {sensor_view.GetSensor().GetName()} dividing cluster {max_cluster_id} resulting with {len(self.cluster_details.keys())} clusters.")
                iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_selected_sensor_views)} sensors.")

        output_path = Path(self.parameters["output_folder"].GetString())
        output_path.mkdir(parents=True, exist_ok=True)
        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_selected_sensor_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    def __AddSensorViewMasks(self, list_of_sensor_views: 'list[SensorViewUnionType]') -> None:
        masking_type = self.parameters["masking"]["masking_type"].GetString()
        if masking_type == "automatic":
            for sensor_view in list_of_sensor_views:
                sensor_view.AddAuxiliaryExpression("mask", KratosDT.MaskUtils.GetMask(sensor_view.GetAuxiliaryExpression("filtered")))
        else:
            raise RuntimeError(f"Unsupported masking_type={masking_type} requested.")

    def __UpdateClusters(self, list_of_sensor_views: 'list[SensorViewUnionType]') -> ExpressionUnionType:
        list_of_masks = [sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views]
        cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks(list_of_masks)
        new_cluster_sensor_ids_list: 'list[list[int]]' = [[list_of_sensor_views[i].GetSensor().GetValue(KratosDT.SENSOR_ID) for i in indices_list] for indices_list, _ in cluster_data]

        list_of_cluster_ids_to_delete: 'list[int]' = list(self.cluster_details.keys())

        # first add ll new clusters
        for i, current_cluster_sensor_ids in enumerate(new_cluster_sensor_ids_list):
            old_cluster_id = -1
            for cluster_id, (cluster_sensor_ids, _) in self.cluster_details.items():
                if cluster_sensor_ids == current_cluster_sensor_ids:
                    old_cluster_id = cluster_id
                    del list_of_cluster_ids_to_delete[list_of_cluster_ids_to_delete.index(cluster_id)]
                    break

            if old_cluster_id == -1:
                self.cluster_details[max(self.cluster_details.keys(), default=0) + 1] = (current_cluster_sensor_ids, cluster_data[i][1])
            else:
                self.cluster_details[old_cluster_id][1].SetExpression(cluster_data[i][1].GetExpression())

        # now remove clusters which are not anymore valid
        for k in list_of_cluster_ids_to_delete:
            del self.cluster_details[k]

        clustering_exp = list_of_masks[0] * 0.0
        for cluster_id, (_, cluster_mask_exp) in self.cluster_details.items():
            clustering_exp += cluster_mask_exp * cluster_id
        return clustering_exp

    def __OutputIteration(self, iteration: int, output_path: Path, vtu_output: Kratos.VtuOutput, sensor_redundancy: ExpressionUnionType, list_of_selected_sensor_views: 'list[SensorViewUnionType]') -> None:
        last_sensor_view = list_of_selected_sensor_views[-1]
        clustering_exp = self.__UpdateClusters(list_of_selected_sensor_views)

        if self.vtu_output_iteration_data:
            self.__OutputSensorViewExpressions(
                        output_path / f"sensor_placement_iteration_{iteration:05d}",
                        vtu_output, last_sensor_view,
                        [
                            ("sensor_redundancy", sensor_redundancy),
                            ("cluster", clustering_exp)
                        ])

        if self.csv_output_iteration_data:
            self.__PrintSensorViewsListToCSV(output_path / f"sensor_placement_iteration_{iteration:05d}.csv", list_of_selected_sensor_views, ["type", "name", "location", "value"])

    @staticmethod
    def __OutputSensorViewExpressions(output_path: Path, vtu_output: Kratos.VtuOutput, sensor_view: SensorViewUnionType, additional_expressions: 'list[tuple[str, ExpressionUnionType]]' = []) -> None:
        vtu_output.ClearCellContainerExpressions()
        vtu_output.ClearNodalContainerExpressions()

        # add the main expression
        vtu_output.AddContainerExpression(sensor_view.GetExpressionName(), sensor_view.GetContainerExpression())

        # now add the auxiliary expression
        for auxiliary_name in sensor_view.GetAuxiliarySuffixes():
            vtu_output.AddContainerExpression(f"{sensor_view.GetExpressionName()}_{auxiliary_name}", sensor_view.GetAuxiliaryExpression(auxiliary_name))

        # now add the additional expressions
        for exp_name, exp in additional_expressions:
            vtu_output.AddContainerExpression(exp_name, exp)

        vtu_output.PrintOutput(str(output_path))

    @staticmethod
    def __PrintSensorViewsListToCSV(output_file_name: Path, list_of_sensor_views: 'list[SensorViewUnionType]', list_of_sensor_properties: 'list[str]') -> None:
        PrintSensorListToCSV(output_file_name, [sensor_view.GetSensor() for sensor_view in list_of_sensor_views], list_of_sensor_properties)

