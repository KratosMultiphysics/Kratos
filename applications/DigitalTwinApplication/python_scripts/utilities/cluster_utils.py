import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class Cluster:
    def __init__(self, cluster_id: int, sensor_views: 'list[SensorViewUnionType]', mask: ContainerExpressionTypes) -> None:
        self.__cluster_id = cluster_id
        self.__sensor_views = sorted(sensor_views, key=lambda x: x.GetSensor().GetValue(KratosDT.SENSOR_ID))
        self.__mask = mask

    def GetId(self) -> int:
        return self.__cluster_id

    def GetSensorViews(self) -> 'list[SensorViewUnionType]':
        return list(self.__sensor_views)

    def GetMask(self) -> ContainerExpressionTypes:
        return self.__mask.Clone()

    def GetDomainSize(self, domain_size_exp: ContainerExpressionTypes) -> float:
        return Kratos.Expression.Utils.Sum(domain_size_exp * self.__mask)

    def GetSensorIds(self) -> 'list[int]':
        return [sensor_view.GetSensor().GetValue(KratosDT.SENSOR_ID) for sensor_view in self.__sensor_views]

    def __eq__(self, other) -> bool:
        if isinstance(other, Cluster):
            return self.GetSensorIds() == other.GetSensorIds()
        else:
            return False

def GetClusters(list_of_sensor_masks: 'list[SensorViewUnionType]') -> 'list[Cluster]':
    clusters_with_indices: 'list[tuple[list[int], ContainerExpressionTypes]]' = KratosDT.MaskUtils.ClusterMasks([sensor_mask.GetContainerExpression() for sensor_mask in list_of_sensor_masks])

    clusters: 'list[Cluster]' = []
    for i, (indices, mask) in enumerate(clusters_with_indices):
        sensor_views = [list_of_sensor_masks[index] for index in indices]
        clusters.append(Cluster(i, sensor_views, mask))

    return clusters

def GetReassignedClusters(origin_clusters: 'list[Cluster]', new_clusters: 'list[Cluster]') -> 'list[Cluster]':
    new_cluster_id = max(origin_clusters, key=lambda x: x.GetId()).GetId() + 1

    clusters: 'list[Cluster]' = []
    for new_cluster in new_clusters:
        if new_cluster in origin_clusters:
            index = origin_clusters.index(new_cluster)
            clusters.append(Cluster(origin_clusters[index].GetId(), new_cluster.GetSensorViews(), new_cluster.GetMask()))
        else:
            clusters.append(Cluster(new_cluster_id, new_cluster.GetSensorViews(), new_cluster.GetMask()))
            new_cluster_id += 1

    return clusters

def GetClusterSensorViews(list_of_clusters: 'list[Cluster]') -> 'list[SensorViewUnionType]':
    added_sensor_view_ids: 'list[int]' = []
    sensor_views: 'list[SensorViewUnionType]' = []
    for cluster in list_of_clusters:
        for sensor_view in cluster.GetSensorViews():
            sensor_id = sensor_view.GetSensor().GetValue(KratosDT.SENSOR_ID)
            if sensor_id not in added_sensor_view_ids:
                sensor_views.append(sensor_view)
                added_sensor_view_ids.append(sensor_id)
    return sensor_views

