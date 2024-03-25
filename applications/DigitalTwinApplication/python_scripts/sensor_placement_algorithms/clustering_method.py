import abc
from KratosMultiphysics.DigitalTwinApplication.utilities.cluster_utils import Cluster
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType

class ClusteringMethod(abc.ABC):
    @abc.abstractmethod
    def Execute(self, current_clusters: 'list[Cluster]', large_clusters: 'list[Cluster]', potential_sensor_views: 'list[SensorViewUnionType]') -> 'tuple[SensorViewUnionType, list[Cluster]]':
        pass