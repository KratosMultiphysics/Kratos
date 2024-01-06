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


class LeastRobustSensorSelection:
    def __init__(self, params: Kratos.Parameters) -> None:
        default_params = Kratos.Parameters("""{
            "sensor_selection_method" : "least_robust",
            "sensor_group_id_varaible": "SENSOR_GROUP_ID",
            "search_radius"           : 10.0,
            "max_number_of_neighbours": 1000
        }""")
        params.RecursivelyValidateAndAssignDefaults(default_params)

        self.search_radius = params["search_radius"].GetDouble()
        self.max_number_of_neighbours = params["max_number_of_neighbours"].GetInt()
        self.sensor_group_id_var = Kratos.KratosGlobals.GetVariable(params["sensor_group_id_varaible"].GetString())

        self.sensor_groups: 'dict[int, tuple[list[KratosDT.Sensors.Sensor], list[ExpressionUnionType]]]' = {}

    def Initialize(self, list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
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

        return sorted(list_of_sensor_views, key=lambda x: x.GetSensor().GetValue(KratosDT.SENSOR_LOCATION_ROBUSTNESS))

    def Select(self, potential_sensor_views: 'list[SensorViewUnionType]', selected_sensor_views: 'list[SensorViewUnionType]') -> SensorViewUnionType:
        # now get the least robustness sensor view
        return min(potential_sensor_views, key=lambda x: x.GetSensor().GetValue(KratosDT.SENSOR_LOCATION_ROBUSTNESS))

class ClusterReductionSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                                   : "cluster_reduction_sensor_placement_algorithm",
            "output_folder"                          : "output",
            "required_maximum_cluster_coverage_ratio": 0.1,
            "maximum_number_of_iterations"           : 100,
            "filtering"                              : {},
            "sensor_removal_settings"                : {},
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

        sensor_removal_method = self.parameters["sensor_removal_settings"]["sensor_selection_method"].GetString()
        if sensor_removal_method == "least_robust":
            self.sensor_selection_method = LeastRobustSensorSelection(self.parameters["sensor_removal_settings"])
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

        list_of_sensor_views = self.sensor_selection_method.Initialize(list_of_sensor_views)

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

        # first compute the clustering with all the sensors
        cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks([sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views])
        _ = UpdateClusters(self.cluster_details, list_of_sensor_views, cluster_data)

        # now we have potential list of sensors which can be removed. Try removing one by one and computing
        # clusters sizes
        iteration = 1
        cluster_domain_size_ratios = self.__ComputeClusterDomainSizeRatios([mask for _, (_, mask) in self.cluster_details.items()], total_domain_size, domain_size_exp)
        while (True):

            found_sensor_to_remove = False
            for index, test_sensor_view in enumerate(list_of_sensor_views):
                is_domain_size_ratios_degraded = False

                # first check if this sensor covers some uniue areas. If so skip from removal
                for _, (sensor_ids, _) in self.cluster_details.items():
                    if len(sensor_ids) == 1 and sensor_ids[0] == test_sensor_view.GetSensor().GetValue(KratosDT.SENSOR_ID):
                        is_domain_size_ratios_degraded = True
                        break

                # now check whether removal of this sensor degrades the clustering
                if not is_domain_size_ratios_degraded:
                    masks: 'list[ExpressionUnionType]' = []
                    for sensor_view in list_of_sensor_views:
                        if sensor_view != test_sensor_view:
                            masks.append(sensor_view.GetAuxiliaryExpression("mask"))

                    # now cluster again
                    cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks(masks)

                    # now calculate cluster domain size ratios
                    current_cluster_domain_size_ratios = self.__ComputeClusterDomainSizeRatios([mask for _, mask in cluster_data], total_domain_size, domain_size_exp)

                    # now check whether clusters have become larger
                    if len(current_cluster_domain_size_ratios) <= len(cluster_domain_size_ratios):
                        # new clustering has found either similar number of clusters
                        # or fewer clusters which are above the threshold.

                        # now checking whether the cluster domain sizes has degraded.
                        for i, p_v in current_cluster_domain_size_ratios:
                            if cluster_domain_size_ratios[i] < p_v:
                                is_domain_size_ratios_degraded = True
                                break

                if not is_domain_size_ratios_degraded:
                    # found a sensor which can be removed.
                    sensor_to_remove = test_sensor_view
                    found_sensor_to_remove = True
                    break

            if not found_sensor_to_remove:
                break

            # now we found list of sensors which can be safely removed.
            index = list_of_sensor_views.index(sensor_to_remove)
            del list_of_sensor_views[index]

            cluster_masks = [sensor_view.GetAuxiliaryExpression("mask") for sensor_view in list_of_sensor_views]
            cluster_data: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.MaskUtils.ClusterMasks(cluster_masks)
            clustering_exp = UpdateClusters(self.cluster_details, list_of_sensor_views, cluster_data)
            cluster_domain_size_ratios = self.__ComputeClusterDomainSizeRatios(cluster_masks, total_domain_size, domain_size_exp)

            if self.vtu_output_iteration_data:
                sensor_redundancy = list_of_sensor_views[0].GetAuxiliaryExpression("mask")
                for sensor_view in list_of_sensor_views[1:]:
                    sensor_redundancy += sensor_view.GetAuxiliaryExpression("mask")
                OutputSensorViewExpressions(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}", vtu_output, sensor_to_remove, [("sensor_redundancy", sensor_redundancy), ("cluster", clustering_exp)])
            if self.csv_output_iteration_data:
                PrintSensorViewsListToCSV(iteration_data_output_path / f"sensor_placement_iteration_{iteration:05d}.csv", list_of_sensor_views, ["type", "name", "location", "value"])

            cluster_sizes = [KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) / total_domain_size for _, (_, mask) in self.cluster_details.items()]
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration: {iteration} - Removed {sensor_to_remove.GetSensor().GetName()} resulting with {len(self.cluster_details.keys())} clusters having domain size ratios between [ {min(cluster_sizes) * 100.0:6.3f}%, {max(cluster_sizes)*100.0:6.3f}% ].")
            iteration += 1

        output_path = Path(self.parameters["output_folder"].GetString())
        output_path.mkdir(parents=True, exist_ok=True)
        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_sensor_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    def __ComputeClusterDomainSizeRatios(self, cluster_masks: 'list[ExpressionUnionType]', total_domain_size: float, domain_size_exp: ExpressionUnionType) -> 'list[float]':
        cluster_domain_size_ratios = [KratosDT.SensorUtils.Sum(KratosDT.MaskUtils.Scale(domain_size_exp, mask)) / total_domain_size for mask in cluster_masks]
        cluster_domain_size_ratios = sorted(cluster_domain_size_ratios, reverse=True)

        for i, v in enumerate(cluster_domain_size_ratios):
            if v <= self.required_maximum_cluster_coverage_ratio:
                return cluster_domain_size_ratios[:i]

        return cluster_domain_size_ratios

