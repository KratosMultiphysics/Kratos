import importlib
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.sensor_placement_algorithm import SensorPlacementAlgorithm
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetMostCoveringSensorView
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import Cluster
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import GetClusters
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import GetReassignedClusters
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorListToJson
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.clustering_method import ClusteringMethod

class ClusteringSensorPlacementAlgorithm(SensorPlacementAlgorithm):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        default_parameters = Kratos.Parameters("""{
            "type"                      : "clustering_sensor_placement_algorithm",
            "output_path"               : "sensor_data",
            "vtu_output"                : false,
            "clustering_method_settings": {},
            "convergence_criteria"      : {
                "max_cluster_size_ratio" : 0.1,
                "max_number_of_iteration": 100
            }
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        parameters["convergence_criteria"].ValidateAndAssignDefaults(default_parameters["convergence_criteria"])

        self.output_path = parameters["output_path"].GetString()
        self.is_vtu_output = parameters["vtu_output"].GetBool()
        self.max_cluster_size_ratio = parameters["convergence_criteria"]["max_cluster_size_ratio"].GetDouble()
        self.max_number_of_iteration = parameters["convergence_criteria"]["max_number_of_iteration"].GetInt()

        clustering_type = parameters["clustering_method_settings"]["type"].GetString()
        clustering_module = importlib.import_module(f"KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.{clustering_type}")
        self.clustering_method: ClusteringMethod = getattr(clustering_module, Kratos.StringUtilities.ConvertSnakeCaseToCamelCase(clustering_type))(model, parameters["clustering_method_settings"])

    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{}""")

    def Execute(self, list_of_sensors: 'list[KratosDT.Sensors.Sensor]') -> None:
        if len(list_of_sensors) == 0:
            raise RuntimeError("No sensors are given.")

        first_sensor = list_of_sensors[0]

        for k in first_sensor.GetNodalExpressionsMap().keys():
            if k.endswith("_mask"):
                list_of_sensor_views = [KratosDT.Sensors.NodalSensorView(sensor, k) for sensor in list_of_sensors]
                self.__AnalyseSensorViews(k[:-5], list_of_sensor_views)

        for k in first_sensor.GetConditionExpressionsMap().keys():
            if k.endswith("_mask"):
                list_of_sensor_views = [KratosDT.Sensors.ConditionSensorView(sensor, k) for sensor in list_of_sensors]
                self.__AnalyseSensorViews(k[:-5], list_of_sensor_views)

        for k in first_sensor.GetElementExpressionsMap().keys():
            if k.endswith("_mask"):
                list_of_sensor_views = [KratosDT.Sensors.ElementSensorView(sensor, k) for sensor in list_of_sensors]
                self.__AnalyseSensorViews(k[:-5], list_of_sensor_views)

    def __AnalyseSensorViews(self, prefix: str, list_of_sensor_masks: 'list[SensorViewUnionType]') -> None:
        best_sensor_views: 'list[SensorViewUnionType]' = []

        base_path = Path(self.output_path) / prefix
        base_path.mkdir(exist_ok=True, parents=True)

        # first find the most coverage sensor as a starting point
        most_covering_sensor_view = GetMostCoveringSensorView(list_of_sensor_masks)
        best_sensor_views.append(most_covering_sensor_view)
        del list_of_sensor_masks[list_of_sensor_masks.index(most_covering_sensor_view)]

        # create the cluster data mask
        cluster_data_mask = most_covering_sensor_view.GetContainerExpression().Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(cluster_data_mask, 1.0)
        clusters = [Cluster(1, [], cluster_data_mask.Clone())]

        # get the initial clusters
        clusters = GetReassignedClusters(clusters, GetClusters([most_covering_sensor_view]))

        if isinstance(cluster_data_mask, Kratos.Expression.ConditionExpression) or isinstance(cluster_data_mask, Kratos.Expression.ElementExpression):
            domain_exp = cluster_data_mask.Clone()
            Kratos.Expression.DomainSizeExpressionIO.Read(domain_exp)
        else:
            raise RuntimeError("Unsupported mask type.")

        total_domain_size = Kratos.Expression.Utils.Sum(domain_exp)
        max_allowed_cluster_size = total_domain_size * self.max_cluster_size_ratio

        cluster_domain_sizes = [(cluster, cluster.GetDomainSize(domain_exp)) for cluster in clusters]
        cluster_domain_sizes = sorted(cluster_domain_sizes, key=lambda x: x[1], reverse=True)
        list_of_large_clusters = [cluster for (cluster, cluster_size) in cluster_domain_sizes if cluster_size > max_allowed_cluster_size]

        iteration = 1
        if self.is_vtu_output:
            vtu_output = Kratos.VtuOutput(most_covering_sensor_view.GetContainerExpression().GetModelPart())
            self.__OutputClustersToVtu(vtu_output, str(base_path / f"iteration_{iteration}"), clusters, most_covering_sensor_view)

        list_of_sensors: list[KratosDT.Sensors.Sensor] = [most_covering_sensor_view.GetSensor()]
        PrintSensorListToCSV(base_path / f"iteration_{iteration}.csv", list_of_sensors, ["name", "location", "SENSOR_ID"])

        while iteration <= self.max_number_of_iteration and len(list_of_large_clusters) > 0:
            iteration += 1
            large_cluster_mask = list_of_large_clusters[0].GetMask().Clone()
            for large_cluster in list_of_large_clusters:
                large_cluster_mask = KratosDT.MaskUtils.Union(large_cluster_mask, large_cluster.GetMask())

            list_of_potential_sensor_mask_indices: 'list[int]' = KratosDT.MaskUtils.GetMasksDividingReferenceMask(large_cluster_mask, [sensor_view.GetContainerExpression() for sensor_view in list_of_sensor_masks])
            list_of_potential_sensor_masks = [list_of_sensor_masks[index] for index in list_of_potential_sensor_mask_indices]

            if len(list_of_potential_sensor_masks) == 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Cannot further subdivide large clusters.")
                break

            new_sensor_view, clusters = self.clustering_method.Execute(clusters, list_of_large_clusters, list_of_potential_sensor_masks)

            cluster_domain_sizes = [(cluster, cluster.GetDomainSize(domain_exp)) for cluster in clusters]
            cluster_domain_sizes = sorted(cluster_domain_sizes, key=lambda x: x[1], reverse=True)
            list_of_large_clusters = [cluster for (cluster, cluster_size) in cluster_domain_sizes if cluster_size > max_allowed_cluster_size]

            if self.is_vtu_output:
                self.__OutputClustersToVtu(vtu_output, str(base_path / f"iteration_{iteration}"), clusters, new_sensor_view)

            list_of_sensors.append(new_sensor_view.GetSensor())
            PrintSensorListToCSV(base_path / f"iteration_{iteration}.csv", list_of_sensors, ["name", "location", "SENSOR_ID"])

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Iteration {iteration}: Added sensor {new_sensor_view.GetSensor().GetName()} resulting with {len(clusters)} clusters")
            if len(list_of_large_clusters) > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Clusters requiring further sub-division:")
                for large_cluster in list_of_large_clusters:
                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tCluster id = {large_cluster.GetId()}, domain size ratio = {large_cluster.GetDomainSize(domain_exp) / total_domain_size:0.4e}")

        PrintSensorListToJson(base_path / "best_sensor_data.json", list_of_sensors)

    def __OutputClustersToVtu(self, vtu_output: Kratos.VtuOutput, prefix: str, clusters: 'list[Cluster]', new_sensor_view: SensorViewUnionType) -> None:
        if len(clusters) == 0:
            raise RuntimeError("Np clusters to output")

        cluster_coverage_exp = clusters[0].GetMask() * 0.0
        for cluster in clusters:
            cluster_coverage_exp += cluster.GetMask() * cluster.GetId()

        vtu_output.ClearCellContainerExpressions()
        vtu_output.ClearNodalContainerExpressions()
        vtu_output.AddContainerExpression("clustering", cluster_coverage_exp.Clone())
        sensor_view_exp_name = new_sensor_view.GetExpressionName()
        base_name = sensor_view_exp_name[:-5]
        sensor_view = type(new_sensor_view)(new_sensor_view.GetSensor(), base_name)
        for auxiliary_suffix in sensor_view.GetAuxiliarySuffixes():
            vtu_output.AddContainerExpression(auxiliary_suffix, sensor_view.GetAuxiliaryExpression(auxiliary_suffix))
        vtu_output.PrintOutput(prefix)

