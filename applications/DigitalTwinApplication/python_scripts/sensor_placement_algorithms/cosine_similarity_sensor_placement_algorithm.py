import typing
import abc
import scipy.cluster.hierarchy as sch
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetFilter, GetDistance
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionFilterUnionType, ExpressionUnionType

ClusterUnionType = typing.Union[
                            KratosDT.ClusterUtils.NodalSensorViewCluster,
                            KratosDT.ClusterUtils.ConditionSensorViewCluster,
                            KratosDT.ClusterUtils.ElementSensorViewCluster
                        ]

class ClusterSensorIdentification(abc.ABC):
    @abc.abstractmethod
    def GetClusterSensors(self, list_of_clusters: 'list[ClusterUnionType]') -> 'list[SensorViewUnionType]':
        pass

class MostOrthogonalSensorIdentification(ClusterSensorIdentification):
    def GetClusterSensors(self, list_of_clusters: list[ClusterUnionType]) -> list[SensorViewUnionType]:
        most_covered_sensitivities: ExpressionUnionType = list_of_clusters[0].GetSensorViews()[0].GetContainerExpression().Clone()
        most_covered_sensitivities = most_covered_sensitivities * 0.0 + 1.0
        most_covered_sensitivities /= KratosOA.ExpressionUtils.NormL2(most_covered_sensitivities)
        most_covered_sensitivities = most_covered_sensitivities.Flatten()

        # first sort list of clusters with the least cosine distance to the most covered
        list_of_value_cluster_pair: 'list[tuple[float, ClusterUnionType]]' = []
        best_coverage_sensor: 'typing.Optional[SensorViewUnionType]' = None
        best_coverage_sensor_inner_product = 0.0
        for cluster in list_of_clusters:
            max_inner_product = 0.0
            current_sensor_max_inner_product: 'typing.Optional[SensorViewUnionType]' = None
            for sensor_view in cluster.GetSensorViews():
                current_inner_product = KratosOA.ExpressionUtils.InnerProduct(most_covered_sensitivities, sensor_view.GetContainerExpression())
                if current_inner_product > max_inner_product:
                    max_inner_product = current_inner_product
                    current_sensor_max_inner_product = sensor_view
            if max_inner_product > best_coverage_sensor_inner_product:
                best_coverage_sensor_inner_product = max_inner_product
                best_coverage_sensor = current_sensor_max_inner_product
            list_of_value_cluster_pair.append((max_inner_product, cluster))

        list_of_sorted_clusters = sorted(list_of_value_cluster_pair, key=lambda x: x[0], reverse=True)

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found best coverage sensor \"{best_coverage_sensor.GetSensor().GetName()}\" with coverage of {best_coverage_sensor_inner_product} from cluster {list_of_sorted_clusters[0][1].Id}.")
        best_sensor_cluster = list_of_clusters[0].Clone()
        best_sensor_cluster.Clear()
        best_sensor_cluster.SetSensorViews([best_coverage_sensor])

        dummy_cluster = best_sensor_cluster.Clone()
        for _, cluster in list_of_sorted_clusters[1:]:
            max_distance = 0.0
            potential_sensor = None
            for sensor_view in cluster.GetSensorViews():
                dummy_cluster.Clear()
                dummy_sensor_views = best_sensor_cluster.GetSensorViews()
                dummy_sensor_views.append(sensor_view)
                dummy_cluster.SetSensorViews(dummy_sensor_views)
                dummy_sensor_views = dummy_cluster.GetSensorViews()
                dummy_distances = dummy_cluster.GetDistances("cosine_distance")
                i = dummy_sensor_views.index(sensor_view)
                n = len(dummy_sensor_views)
                min_sensor_index = min([j for j in range(n) if i != j], key=lambda j: GetDistance(i, j, dummy_distances))
                min_distance = GetDistance(i, min_sensor_index, dummy_distances)
                if min_distance > max_distance:
                    max_distance = min_distance
                    potential_sensor = sensor_view

            current_sensor_views = best_sensor_cluster.GetSensorViews()
            current_sensor_views.append(potential_sensor)
            best_sensor_cluster.SetSensorViews(current_sensor_views)
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found sensor \"{potential_sensor.GetSensor().GetName()}\" from cluster {cluster.Id}.")

        return best_sensor_cluster.GetSensorViews()


class CosineSimilaritySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "clustering_method"                : "average",
            "cluster_representation_method"    : "average_vector",
            "sensor_coverage_percentage"       : 5.0,
            "max_clustering_iterations"        : 100,
            "best_sensor_identification_method": "highest_minimum",
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()
        self.best_sensor_identification_method = self.parameters["best_sensor_identification_method"].GetString()
        self.cluster_representation_method = self.parameters["cluster_representation_method"].GetString()

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

        domain_sensor_view_cluster_data = self.domain_sensor_view_cluster_type(normalized_sensor_views)
        domain_sensor_view_cluster_data.AddDistances("cosine_distance", GetCosineDistances(normalized_sensor_views))

        model_part = domain_sensor_view_cluster_data.GetModelPart()
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()

        # prepare for output
        if self.is_vtu_output or self.is_csv_output:
            output_path = Path(self.parameters["output_folder"].GetString()) / f"sensors/{name}"
            output_path.mkdir(parents=True, exist_ok=True)

        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(model_part)
            if isinstance(domain_sensor_view_cluster_data, KratosDT.ClusterUtils.NodalSensorViewClusterData):
                vtu_output.AddNonHistoricalVariable(KratosDT.SENSOR_CLUSTER_ID, vtu_output.NODES)
            elif isinstance(domain_sensor_view_cluster_data, KratosDT.ClusterUtils.ConditionSensorViewClusterData):
                vtu_output.AddNonHistoricalVariable(KratosDT.SENSOR_CLUSTER_ID, vtu_output.CONDITIONS)
            elif isinstance(domain_sensor_view_cluster_data, KratosDT.ClusterUtils.ElementSensorViewClusterData):
                vtu_output.AddNonHistoricalVariable(KratosDT.SENSOR_CLUSTER_ID, vtu_output.ELEMENTS)

        dummy_exp = normalized_sensor_views[0].GetContainerExpression()

        # create the first cluster
        first_cluster = self.cluster_type(1, domain_sensor_view_cluster_data)
        first_cluster.SetSensorViews(normalized_sensor_views)
        first_cluster.SetEntities(dummy_exp.GetContainer())
        maximum_allowed_domain_size_per_cluster = KratosDT.SensorUtils.GetDomainSize(first_cluster.GetEntities(), data_communicator) * self.parameters["sensor_coverage_percentage"].GetDouble() / 100.0

        list_of_clusters: 'list[ClusterUnionType]' = [first_cluster]
        list_of_clusters_to_divide: 'list[int]' = [0]
        list_of_cluster_ids_not_divisible: list[int] = []
        for clustering_iteration in range(1, self.parameters["max_clustering_iterations"].GetInt() + 1):
            # no more clusters to divide, then exit the while loop
            if len(list_of_clusters_to_divide) == 0:
                break

            Kratos.Logger.PrintInfo(self.__class__.__name__,f"Clustering iteration {clustering_iteration}")

            # now try to divide the clusters
            for cluster_index in list_of_clusters_to_divide:
                cluster = list_of_clusters[cluster_index]

                cosine_distances = cluster.GetDistances("cosine_distance")
                new_clusters = self.DivideCluster(cluster, [cluster.Id, len(list_of_clusters) + 1], cosine_distances)

                if len(new_clusters) == 2:
                    list_of_clusters[cluster_index] = new_clusters[0]
                    list_of_clusters.append(new_clusters[1])
                else:
                    list_of_cluster_ids_not_divisible.append(cluster.Id)

            # now compute the average for each cluster
            cluster_expressions_list: 'list[ExpressionUnionType]' = []
            for cluster in list_of_clusters:
                current_exp = dummy_exp.Clone() * 0.0
                for sensor_view in cluster.GetSensorViews():
                    current_exp += sensor_view.GetContainerExpression()
                if self.cluster_representation_method == "unit_vector":
                    cluster_expressions_list.append((current_exp / KratosOA.ExpressionUtils.NormL2(current_exp)).Flatten())
                elif self.cluster_representation_method == "average_vector":
                    cluster_expressions_list.append((current_exp / len(cluster.GetSensorViews())).Flatten())
                else:
                    raise RuntimeError(f"Unsupported cluster representation method.")

            KratosDT.SensorUtils.AssignEntitiesToClusters(list_of_clusters, cluster_expressions_list)

            # now check whether the cluster is too large
            list_of_clusters_to_divide.clear()
            for i, cluster in enumerate(list_of_clusters):
                if cluster.Id not in list_of_cluster_ids_not_divisible and KratosDT.SensorUtils.GetDomainSize(cluster.GetEntities(), data_communicator) > maximum_allowed_domain_size_per_cluster:
                    list_of_clusters_to_divide.append(i)

            if self.is_vtu_output:
                # now set cluster id in every entity
                for cluster in list_of_clusters:
                    Kratos.VariableUtils().SetNonHistoricalVariable(KratosDT.SENSOR_CLUSTER_ID, cluster.Id, cluster.GetEntities())
                vtu_output.PrintOutput(str(output_path / f"{model_part.FullName()}_{clustering_iteration:05d}"))
            if self.is_csv_output:
                for cluster in list_of_clusters:
                    for sensor_view in cluster.GetSensorViews():
                        sensor_view.GetSensor().SetValue(KratosDT.SENSOR_CLUSTER_ID, cluster.Id)
                PrintSensorListToCSV(output_path / f"sensor_cluster_iteration_{clustering_iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in normalized_sensor_views], ["type", "name", "location", "value", "SENSOR_CLUSTER_ID"])

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_clusters)} clusters.")

        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in MostOrthogonalSensorIdentification().GetClusterSensors(list_of_clusters)]
        PrintSensorListToJson(output_path / "best_sensor_data.json", list_of_best_sensors)
        PrintSensorListToCSV(output_path / "best_sensor_data.csv", list_of_best_sensors, ["type", "name", "location", "value", "SENSOR_CLUSTER_ID"])

    def DivideCluster(self, cluster_to_divide: ClusterUnionType, new_cluster_ids: 'list[int]', cosine_distances: 'list[float]') -> 'list[ClusterUnionType]':
        number_of_divisions = len(new_cluster_ids)
        list_of_sensor_views = cluster_to_divide.GetSensorViews()
        number_of_sensor_views = len(list_of_sensor_views)

        if number_of_sensor_views == 1:
            # not possible to divide. Hence returning the original cluster
            return [cluster_to_divide]

        Kratos.Logger.PrintInfo("", f"\tDividing the cluster {cluster_to_divide.Id} having {number_of_sensor_views} sensors to {number_of_divisions} clusters...")

        sensor_cosine_linkage = sch.linkage(cosine_distances, method=self.parameters["clustering_method"].GetString())
        sensor_cosine_clusters = sch.fcluster(sensor_cosine_linkage, len(new_cluster_ids), 'maxclust')

        cluster_id_map: 'dict[int, list[SensorViewUnionType]]' = {}
        for i, cluster_id in enumerate(sensor_cosine_clusters):
            if cluster_id not in cluster_id_map.keys():
                cluster_id_map[cluster_id] = []
            cluster_id_map[cluster_id].append(list_of_sensor_views[i])

        list_of_clusters: 'list[ClusterUnionType]' = []
        for cluster_id, sensor_views in cluster_id_map.items():
            new_cluster = self.cluster_type(new_cluster_ids[cluster_id-1], cluster_to_divide.GetDataContainer())
            new_cluster.SetSensorViews(sensor_views)
            list_of_clusters.append(new_cluster)
        return list_of_clusters

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





