import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.sensor_placement_algorithms.clustering_method import ClusteringMethod
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import Cluster
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import GetClusters
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import GetReassignedClusters

class LargestSplittingClusteringMethod(ClusteringMethod):
    def __init__(self, _: Kratos.Model, parameters: Kratos.Parameters) -> None:
        defaults = Kratos.Parameters("""{
            "type"      : "largest_splitting_clustering_method",
            "echo_level": 0
        }""")
        parameters.ValidateAndAssignDefaults(defaults)
        self.echo_level = parameters["echo_level"].GetInt()

    def Execute(self, current_clusters: list[Cluster], large_clusters: 'list[Cluster]', potential_sensor_views: 'list[SensorViewUnionType]') -> 'tuple[SensorViewUnionType, list[Cluster]]':
        sensor_view_splits: 'dict[SensorViewUnionType, int]' = {}
        for potential_sensor_view in potential_sensor_views:
            sensor_view_splits[potential_sensor_view] = 0

        potential_masks = [sensor_view.GetContainerExpression() for sensor_view in potential_sensor_views]
        indices: list[int] = KratosDT.MaskUtils.GetMasksDividingReferenceMask(large_clusters[0].GetMask(), potential_masks)
        list_of_most_splitting_sensors = [potential_sensor_views[index] for index in indices]

        # now find the most distance with the current sensor views.
        origin_sensors: list[KratosDT.Sensors.Sensor] = []
        origin_sensor_ids: 'list[int]' = []
        origin_sensor_views: 'list[SensorViewUnionType]' = []
        for cluster in current_clusters:
            for sensor_view in cluster.GetSensorViews():
                sensor_view_id = sensor_view.GetSensor().GetValue(KratosDT.SENSOR_ID)
                if sensor_view_id not in origin_sensor_ids:
                    origin_sensors.append(sensor_view.GetSensor())
                    origin_sensor_views.append(sensor_view)
                    origin_sensor_ids.append(sensor_view_id)

        test_sensors = [sensor_view.GetSensor() for sensor_view in list_of_most_splitting_sensors]
        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"\tFound following potential sensor:")
            for test_sensor in test_sensors:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"\t\t{test_sensor.GetName()}")

        index = KratosDT.SensorUtils.GetMostDistanced(origin_sensors, test_sensors)
        sensor_view = list_of_most_splitting_sensors[index]

        origin_sensor_views.append(sensor_view)
        return sensor_view, GetReassignedClusters(current_clusters, GetClusters(origin_sensor_views))

