import typing
import abc
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from pathlib import Path
import pyomo.environ as pyomo

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetNormalizedSensorViews
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
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
class CosineSimilaritySensorPlacementAlgorithm(SensorPlacementAlgorithm):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                             : "cosine_similarity_sensor_placement",
            "output_to_vtu"                    : true,
            "output_to_csv"                    : true,
            "output_folder"                    : "sensor_placement/cluster_average",
            "clustering_method"                : "average",
            "sensor_coverage_percentage"       : 5.0,
            "max_clustering_iterations"        : 100,
            "filtering"                        : {}
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        self.model = model
        self.parameters = parameters
        self.list_of_sensors:'list[KratosDT.Sensors.Sensor]' = []
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.is_vtu_output = self.parameters["output_to_vtu"].GetBool()
        self.is_csv_output = self.parameters["output_to_csv"].GetBool()

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
        compressed_distances = GetCosineDistances(normalized_sensor_views)
        domain_sensor_view_cluster_data.AddDistances("cosine_distance", compressed_distances)
        distances = 1.0 - squareform(compressed_distances)

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
        list_of_cluster_representative_sensors: 'list[SensorViewUnionType]' = [GetBestCoverageSensorView(normalized_sensor_views)]
        list_of_clusters_to_divide: 'list[int]' = [0]
        list_of_cluster_ids_not_divisible: list[int] = []

        for clustering_iteration in range(1, self.parameters["max_clustering_iterations"].GetInt() + 1):
            # first assign entities to clusters
            KratosDT.SensorUtils.AssignEntitiesToClustersBasedOnOptimalSensor(list_of_clusters, [sensor_view.GetContainerExpression() for sensor_view in list_of_cluster_representative_sensors])

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
                PrintSensorListToCSV(output_path / f"cluster_representative_sensor_iteration_{clustering_iteration:05d}.csv", [sensor_view.GetSensor() for sensor_view in list_of_cluster_representative_sensors], ["type", "name", "location", "value", "SENSOR_CLUSTER_ID"])

            # now check whether the cluster is too large
            list_of_clusters_to_divide.clear()
            for i, cluster in enumerate(list_of_clusters):
                if cluster.Id not in list_of_cluster_ids_not_divisible and KratosDT.SensorUtils.GetDomainSize(cluster.GetEntities(), data_communicator) > maximum_allowed_domain_size_per_cluster:
                    list_of_clusters_to_divide.append(i)

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

            # create the mixed integer non-linear programming model
            model = pyomo.ConcreteModel()
            model.N = pyomo.Set(initialize=range(len(normalized_sensor_views)))
            model.w = pyomo.Var()
            model.x_list = pyomo.Var(model.N, domain=pyomo.Binary)
            model.objective = pyomo.Objective(expr=model.w, sense=pyomo.minimize)

            # add the equality constraints
            model.constraints_list = pyomo.ConstraintList()
            for cluster in list_of_clusters:
                model.constraints_list.add(sum([model.x_list[i] for i in cluster.GetSensorViewIndices()]) == 1)

            # add the inequality constraints
            for i, cluster_i in enumerate(list_of_clusters):
                for cluster_j in list_of_clusters[i+1:]:
                    z_ineuqality_constraint = 0.0
                    for i_x_index in cluster_i.GetSensorViewIndices():
                        for j_x_index in cluster_j.GetSensorViewIndices():
                            z_ineuqality_constraint += distances[i_x_index, j_x_index] * model.x_list[i_x_index] * model.x_list[j_x_index]
                    model.constraints_list.add(expr=z_ineuqality_constraint <= model.w)

            pyomo.SolverFactory("mindtpy").solve(model, tee=False)

            list_of_cluster_representative_sensors.clear()
            for cluster in list_of_clusters:
                for cluster_index, sensor_view in zip(cluster.GetSensorViewIndices(), cluster.GetSensorViews()):
                    if pyomo.value(model.x_list[cluster_index]) == 1:
                        list_of_cluster_representative_sensors.append(sensor_view)

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Found {len(list_of_clusters)} clusters.")

        list_of_best_sensors = [sensor_view.GetSensor() for sensor_view in list_of_cluster_representative_sensors]
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





