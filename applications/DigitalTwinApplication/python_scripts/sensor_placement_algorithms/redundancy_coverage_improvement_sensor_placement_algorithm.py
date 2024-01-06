import typing
from pathlib import Path
import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetFilter, GetSensorCoverageMasks, GetMostDistancedMax, GetMostDistancedMin
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionFilterUnionType, ExpressionUnionType

ClusterUnionType = typing.Union[
                            KratosDT.ClusterUtils.NodalSensorViewCluster,
                            KratosDT.ClusterUtils.ConditionSensorViewCluster,
                            KratosDT.ClusterUtils.ElementSensorViewCluster
                        ]
class RedundancyCoverageImprovementSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement",
            "required_maximum_overall_coverage": 0.9,
            "required_maximum_cluster_coverage": 0.1,
            "relaxation_max_coverage"          : 1e-9,
            "relaxation_min_coverage"          : 0.05,
            "maximum_number_of_iterations"     : 100,
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.required_maximum_overall_coverage = self.parameters["required_maximum_overall_coverage"].GetDouble()
        self.required_maximum_cluster_coverage = self.parameters["required_maximum_cluster_coverage"].GetDouble()
        self.relaxation_max_coverage = self.parameters["relaxation_max_coverage"].GetDouble()
        self.relaxation_min_coverage = self.parameters["relaxation_min_coverage"].GetDouble()
        self.maximum_number_of_iterations = self.parameters["maximum_number_of_iterations"].GetInt()

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

        normalized_sensor_views = GetNormalizedSensorViews(list_of_sensor_views)
        dummy_exp = normalized_sensor_views[0].GetContainerExpression()
        domain_size_exp = dummy_exp.Clone()
        Kratos.Expression.EntityDomainSizeExpressionIO.Read(domain_size_exp)
        coverage_masks = GetSensorCoverageMasks(normalized_sensor_views, domain_size_exp)
        total_domain_size = KratosDT.SensorUtils.Sum(domain_size_exp)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(normalized_sensor_views[0].GetContainerExpression().GetModelPart())
            for coverage_mask, sensor_view in zip(coverage_masks, normalized_sensor_views):
                vtu_output.AddContainerExpression(sensor_view.GetSensor().GetName(), coverage_mask)
            vtu_output.PrintOutput(str(output_path / "sensor_coverage"))

        if self.is_csv_output:
            PrintSensorListToCSV(output_path / "sensor_coverage.csv", [sensor_view.GetSensor() for sensor_view in normalized_sensor_views], ["type", "name", "location", "value", "SENSOR_COVERAGE"])

        overall_coverage_mask = (dummy_exp * 0.0).Flatten()
        max_coverage = 0.0
        iteration = 1
        list_of_selected_views: 'list[SensorViewUnionType]' = []
        list_of_selected_coverage_masks: 'list[ExpressionUnionType]' = []
        clusters: 'dict[int, ExpressionUnionType]' = {}
        while (max_coverage < self.required_maximum_overall_coverage and iteration <= self.maximum_number_of_iterations):
            # find the sensor mask which increases coverage to the max
            list_of_coverages: 'list[float]' = []
            for i, coverage_mask in enumerate(coverage_masks):
                # computing the additional area which is covered by this coverage mask
                modified_overall_coverage_mask = overall_coverage_mask + coverage_mask
                KratosDT.ControlUtils.ClipContainerExpression(modified_overall_coverage_mask, 0, 1)
                coverage = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(modified_overall_coverage_mask)) / total_domain_size
                list_of_coverages.append(coverage)

            best_coverage_index = GetMostDistancedMax(self.relaxation_max_coverage, list_of_coverages, normalized_sensor_views, list_of_selected_views)
            list_of_selected_views.append(normalized_sensor_views[best_coverage_index])
            list_of_selected_coverage_masks.append(coverage_masks[best_coverage_index])
            overall_coverage_mask += coverage_masks[best_coverage_index]
            overall_coverage_mask = overall_coverage_mask.Flatten()
            max_coverage = list_of_coverages[best_coverage_index]

            if self.is_vtu_output:
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()
                self.UpdateClusters(clusters, KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(list_of_selected_coverage_masks))
                cluster_mask = dummy_exp * 0.0
                for cluster_id, mask_expression in clusters.items():
                    cluster_mask += mask_expression * cluster_id
                vtu_output.AddContainerExpression(f"clusters_overall", cluster_mask)
                vtu_output.AddContainerExpression("sensor_coverage", coverage_masks[best_coverage_index].Clone())
                vtu_output.AddContainerExpression("redundancy", overall_coverage_mask.Clone())
                vtu_output.PrintOutput(str(output_path / f"redundancy_iteration_{iteration:05d}"))

            if self.is_csv_output:
                PrintSensorListToCSV(output_path / f"redundancy_iteration_{iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_selected_views], ["type", "name", "location", "value"])

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found sensor {normalized_sensor_views[best_coverage_index].GetSensor().GetName()} increasing coverage to {max_coverage * 100.0:6.3f}%.")

            # now remove them from the remaining searches
            del normalized_sensor_views[best_coverage_index]
            del coverage_masks[best_coverage_index]

            iteration += 1

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_selected_views)} sensors covering {max_coverage * 100.0:6.3f}% domain.")

        iteration = 0
        clusters: 'dict[int, ExpressionUnionType]' = {}
        self.UpdateClusters(clusters, KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(list_of_selected_coverage_masks))
        maxmum_allowed_cluster_size = self.required_maximum_cluster_coverage * total_domain_size
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Allowed maximum cluster coverage area is {maxmum_allowed_cluster_size}")
        cluster_domain_sizes = {cluster_id: KratosDT.SensorUtils.Sum(domain_size_exp.Scale(exp)) for cluster_id, exp in clusters.items()}
        cluster_id__max_domain_size = max(cluster_domain_sizes.items(), key=lambda x: x[1])
        while (iteration <= self.maximum_number_of_iterations and cluster_id__max_domain_size[1] > maxmum_allowed_cluster_size):
            # now check for all the remaining potential sensors and store cluster coverages for each potential cluster
            list_of_values: 'list[float]' = []
            for i, coverage_mask in enumerate(coverage_masks):
                # find the potential clustering with new coverage mask
                list_of_selected_coverage_masks.append(coverage_mask)
                potential_cluster_mask_indices_mask_expressions_list: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(list_of_selected_coverage_masks)
                del list_of_selected_coverage_masks[-1]

                if max([KratosDT.SensorUtils.Sum(domain_size_exp.Scale(mask)) for _, mask in  potential_cluster_mask_indices_mask_expressions_list]) < cluster_id__max_domain_size[1]:
                    list_of_values.append(0.0)
                else:
                    list_of_values.append(1e+100)

            if min(list_of_values) == 1e+100:
                # we have not found any sensor which reduces the max cluster size. hence terminating.
                break

            minimum_cov_sensor_index = GetMostDistancedMin(self.relaxation_min_coverage, list_of_values, normalized_sensor_views, list_of_selected_views)
            list_of_selected_views.append(normalized_sensor_views[minimum_cov_sensor_index])
            list_of_selected_coverage_masks.append(coverage_masks[minimum_cov_sensor_index])
            overall_coverage_mask += coverage_masks[minimum_cov_sensor_index]
            overall_coverage_mask = overall_coverage_mask.Flatten()

            self.UpdateClusters(clusters, KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(list_of_selected_coverage_masks))
            cluster_domain_sizes = {cluster_id: KratosDT.SensorUtils.Sum(domain_size_exp.Scale(exp)) for cluster_id, exp in clusters.items()}
            cluster_id__max_domain_size = max(cluster_domain_sizes.items(), key=lambda x: x[1])

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found a sensor {normalized_sensor_views[minimum_cov_sensor_index].GetSensor().GetName()} resulting with {cluster_id__max_domain_size[1]:0.3f} max cluster domain size with {len(clusters.keys())} clusters.")

            del normalized_sensor_views[minimum_cov_sensor_index]
            del coverage_masks[minimum_cov_sensor_index]

            # now vtu output of the current clustering
            if self.is_vtu_output:
                vtu_output.ClearCellContainerExpressions()
                vtu_output.ClearNodalContainerExpressions()
                cluster_mask = dummy_exp * 0.0
                for cluster_id, mask_expression in clusters.items():
                    cluster_mask += mask_expression * cluster_id
                vtu_output.AddContainerExpression(f"clusters_overall", cluster_mask)
                vtu_output.AddContainerExpression("redundancy", overall_coverage_mask.Clone())
                vtu_output.AddContainerExpression("sensor_coverage", list_of_selected_coverage_masks[-1].Clone())
                vtu_output.PrintOutput(str(output_path / f"clustering_iteration_{iteration:05d}"))

            # now csv output of the current clustering
            if self.is_csv_output:
                PrintSensorListToCSV(output_path / f"clustering_iteration_{iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_selected_views], ["type", "name", "location", "value"])
                PrintSensorListToJson(output_path / f"clustering_iteration_{iteration:05d}.json", [sensor_view.GetSensor() for sensor_view in list_of_selected_views])

            # compute area of each cluster and find the clusters which needs to be divided further
            iteration += 1

        # now the phase 3. Try to remove sensors and check whether we can keep the same overall coverage and
        # largest cluster size
        overall_coverage_mask_exp = overall_coverage_mask.Clone()
        KratosDT.ControlUtils.ClipContainerExpression(overall_coverage_mask_exp, 0, 1.0)
        overall_coverage = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(overall_coverage_mask_exp))
        max_cluster_size = max([KratosDT.SensorUtils.Sum(domain_size_exp.Scale(mask)) for _, mask in clusters.items()])

        iteration = 0
        while True:
            iteration += 1
            # first find all possible sensors to remove
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration {iteration}: Try to remove sensors...")
            list_of_sensor_indices_to_remove: 'list[int]' = []
            index = 0
            while (index < len(list_of_selected_views)):
                temp_selected_views = list(list_of_selected_views)
                temp_selected_masks = list(list_of_selected_coverage_masks)

                # remove the index
                del temp_selected_views[index]
                del temp_selected_masks[index]

                #now check the clustering
                potential_cluster_mask_indices_mask_expressions_list: 'list[tuple[list[int], ExpressionUnionType]]' = KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(temp_selected_masks)

                # now check if current clusters are within the bounds
                current_max_cluster_size = max([KratosDT.SensorUtils.Sum(domain_size_exp.Scale(mask)) for _, mask in potential_cluster_mask_indices_mask_expressions_list])
                current_coverage_exp = dummy_exp * 0.0
                for _, mask in potential_cluster_mask_indices_mask_expressions_list:
                    current_coverage_exp += mask
                KratosDT.ControlUtils.ClipContainerExpression(current_coverage_exp, 0, 1)
                current_coverage = KratosDT.SensorUtils.Sum(domain_size_exp.Scale(current_coverage_exp))

                if current_coverage >= overall_coverage and current_max_cluster_size <= max_cluster_size:
                    list_of_sensor_indices_to_remove.append(index)
                index += 1

            # if there are no more sensors to be removed, exit the loop
            if len(list_of_sensor_indices_to_remove) == 0:
                break

            # now remove the most insignificant sensor from the potential sensors to remove.
            # which will be the most closest one to existing sensors
            min_distance = 1e+9
            min_distance_index = -1
            for i in list_of_sensor_indices_to_remove:
                max_distance = 0.0
                for sensor_view in list_of_selected_views:
                    current_distance = sensor_view.GetSensor().GetLocation() - list_of_selected_views[i].GetSensor().GetLocation()
                    max_distance = max(max_distance, current_distance[0]**2 + current_distance[1]**2 + current_distance[2]**2)
                if min_distance > max_distance:
                    min_distance = max_distance
                    min_distance_index = i

            if min_distance_index != -1:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"--- Removing {list_of_selected_views[min_distance_index].GetSensor().GetName()}.")
                removal_mask = list_of_selected_coverage_masks[min_distance_index]
                overall_coverage_mask -= list_of_selected_coverage_masks[min_distance_index]
                del list_of_selected_views[min_distance_index]
                del list_of_selected_coverage_masks[min_distance_index]

                self.UpdateClusters(clusters, KratosDT.SensorUtils.ClusterBasedOnCoverageMasks(list_of_selected_coverage_masks))

                if self.is_vtu_output:
                    vtu_output.ClearCellContainerExpressions()
                    vtu_output.ClearNodalContainerExpressions()
                    cluster_mask = dummy_exp * 0.0
                    for cluster_id, mask_expression in clusters.items():
                        cluster_mask += mask_expression * cluster_id
                    vtu_output.AddContainerExpression(f"clusters_overall", cluster_mask)
                    vtu_output.AddContainerExpression("redundancy", overall_coverage_mask.Clone())
                    vtu_output.AddContainerExpression("sensor_coverage", removal_mask.Clone())
                    vtu_output.PrintOutput(str(output_path / f"removal_iteration_{iteration:05d}"))

                # now csv output of the current clustering
                if self.is_csv_output:
                    PrintSensorListToCSV(output_path / f"removal_iteration_{iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_selected_views], ["type", "name", "location", "value"])
                    PrintSensorListToJson(output_path / f"removal_iteration_{iteration:05d}.json", [sensor_view.GetSensor() for sensor_view in list_of_selected_views])

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_selected_views)} sensors.")
        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_selected_views]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value"])

    @staticmethod
    def UpdateClusters(current_clusters: 'dict[int, ExpressionUnionType]', cluster_mask_indices_mask_expressions_list: 'list[tuple[list[int], ExpressionUnionType]]') -> None:
        found_cluster_ids: 'list[int]' = []
        for _, new_mask in cluster_mask_indices_mask_expressions_list:
            # first try to find the mask indices in the existing
            found_cluster = False
            for cluster_id, cluster_mask in current_clusters.items():
                if KratosOA.ExpressionUtils.NormInf(cluster_mask - new_mask) == 0:
                    found_cluster = True
                    found_cluster_ids.append(cluster_id)
                    break

            if not found_cluster:
                if len(current_clusters.keys()) > 0:
                    new_cluster_id = max(current_clusters.keys()) + 1
                else:
                    new_cluster_id = 1
                current_clusters[new_cluster_id] = new_mask
                found_cluster_ids.append(new_cluster_id)

        # check and remove all unnecessary clusters in current_clusters
        list_of_clusters_to_remove: 'list[int]' = [cluster_id for cluster_id in current_clusters.keys() if cluster_id not in found_cluster_ids]
        for cluster_id in list_of_clusters_to_remove:
            del current_clusters[cluster_id]

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





